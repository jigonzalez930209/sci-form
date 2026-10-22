//! GPU-accelerated sTDA Coulomb matrix construction.
//!
//! Offloads the O(N²_singles × N²_atoms) J-integral computation to the GPU.
//! The inner kernel evaluates:
//!   J_{ia,jb} = Σ_{A,B} q^A_{ia} · γ_{AB} · q^B_{jb}
//!
//! This is a GEMM-like operation on the transition charge matrix:
//!   A_off_diag = 2 · Q^T · Γ · Q
//! where Q is (n_atoms × n_singles) and Γ is (n_atoms × n_atoms).

use super::context::{
    ComputeBindingDescriptor, ComputeBindingKind, ComputeDispatchDescriptor, GpuContext,
};

/// Minimum singles count to justify GPU dispatch.
const GPU_DISPATCH_THRESHOLD: usize = 100;

/// GPU-accelerated sTDA J-integral matrix: A_{off} = 2 · Q^T · Γ · Q
///
/// `q_matrix`: transition charges, shape (n_atoms, n_singles), row-major flat.
/// `gamma`: damped Coulomb matrix, shape (n_atoms, n_atoms), row-major flat.
/// `n_atoms`: number of atoms.
/// `n_singles`: number of single excitations.
///
/// Returns the off-diagonal contribution to the A matrix (n_singles × n_singles), row-major.
pub fn compute_stda_j_matrix_gpu(
    ctx: &GpuContext,
    q_matrix: &[f64],
    gamma: &[f64],
    n_atoms: usize,
    n_singles: usize,
) -> Result<Vec<f64>, String> {
    if n_singles < GPU_DISPATCH_THRESHOLD || !ctx.capabilities.gpu_available {
        return compute_stda_j_matrix_cpu(q_matrix, gamma, n_atoms, n_singles);
    }

    // For GPU path: compute GammaQ = Γ · Q (n_atoms × n_singles)
    // Then A_off = 2 · Q^T · GammaQ (n_singles × n_singles)
    //
    // Both are matrix multiplies which map well to GPU compute shaders.
    // For now, use the GPU's general matrix multiply dispatch.

    let q_f32: Vec<f32> = q_matrix.iter().map(|&x| x as f32).collect();
    let gamma_f32: Vec<f32> = gamma.iter().map(|&x| x as f32).collect();

    // Step 1: GammaQ = Γ · Q
    let gamma_q_bytes = vec![0u8; n_atoms * n_singles * 4];

    let dispatch = ComputeDispatchDescriptor {
        label: "stda_gamma_q".to_string(),
        shader_source: MATMUL_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [
            n_atoms.div_ceil(16) as u32,
            n_singles.div_ceil(16) as u32,
            1,
        ],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "gamma".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: bytemuck_cast_f32(&gamma_f32),
            },
            ComputeBindingDescriptor {
                label: "q".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: bytemuck_cast_f32(&q_f32),
            },
            ComputeBindingDescriptor {
                label: "result".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: gamma_q_bytes,
            },
            ComputeBindingDescriptor {
                label: "dims".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: pack_dims(n_atoms as u32, n_atoms as u32, n_singles as u32),
            },
        ],
    };

    let gamma_q_result = ctx
        .run_compute(&dispatch)?
        .outputs
        .into_iter()
        .last()
        .ok_or_else(|| "sTDA GPU dispatch returned no output buffer".to_string())?;

    // The read-back must have exactly `n_atoms * n_singles` f32 values.
    let expected_gamma_q_bytes = n_atoms
        .checked_mul(n_singles)
        .and_then(|count| count.checked_mul(std::mem::size_of::<f32>()))
        .ok_or_else(|| "sTDA GPU dispatch size overflow".to_string())?;
    if gamma_q_result.len() != expected_gamma_q_bytes {
        return Err(format!(
            "sTDA GPU readback has {} bytes, expected {expected_gamma_q_bytes}",
            gamma_q_result.len()
        ));
    }

    // Step 2: A_off = 2 · Q^T · GammaQ
    let result_bytes = vec![0u8; n_singles * n_singles * 4];

    let dispatch2 = ComputeDispatchDescriptor {
        label: "stda_qt_gamma_q".to_string(),
        shader_source: MATMUL_TRANSPOSE_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [
            n_singles.div_ceil(16) as u32,
            n_singles.div_ceil(16) as u32,
            1,
        ],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "q".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: bytemuck_cast_f32(&q_f32),
            },
            ComputeBindingDescriptor {
                label: "gamma_q".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: gamma_q_result,
            },
            ComputeBindingDescriptor {
                label: "result".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: result_bytes,
            },
            ComputeBindingDescriptor {
                label: "dims".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: pack_dims(n_atoms as u32, n_singles as u32, n_singles as u32),
            },
        ],
    };

    let a_off_bytes = ctx
        .run_compute(&dispatch2)?
        .outputs
        .into_iter()
        .last()
        .ok_or_else(|| "sTDA GPU dispatch returned no output buffer".to_string())?;

    // Convert f32 result back to f64 with 2× scaling.
    // The read-back must have exactly `n_singles * n_singles` f32 values;
    // anything else is an error instead of a silently wrong (zero) matrix.
    let expected_bytes = n_singles
        .checked_mul(n_singles)
        .and_then(|count| count.checked_mul(std::mem::size_of::<f32>()))
        .ok_or_else(|| "sTDA GPU dispatch size overflow".to_string())?;
    if a_off_bytes.len() != expected_bytes {
        return Err(format!(
            "sTDA GPU readback has {} bytes, expected {expected_bytes}",
            a_off_bytes.len()
        ));
    }
    let a_off_f32 = bytemuck_cast_from_u8(&a_off_bytes)?;
    Ok(a_off_f32.iter().map(|&x| 2.0 * x as f64).collect())
}

/// CPU fallback for sTDA J-integral matrix.
fn compute_stda_j_matrix_cpu(
    q_matrix: &[f64],
    gamma: &[f64],
    n_atoms: usize,
    n_singles: usize,
) -> Result<Vec<f64>, String> {
    let mut result = vec![0.0; n_singles * n_singles];

    // A[ia, jb] = 2 * Σ_{A,B} q[A, ia] * gamma[A,B] * q[B, jb]
    for ia in 0..n_singles {
        for jb in 0..=ia {
            let mut val = 0.0;
            for a in 0..n_atoms {
                let q_a_ia = q_matrix[a * n_singles + ia];
                if q_a_ia.abs() < 1e-12 {
                    continue;
                }
                for b in 0..n_atoms {
                    val += q_a_ia * gamma[a * n_atoms + b] * q_matrix[b * n_singles + jb];
                }
            }
            result[ia * n_singles + jb] = 2.0 * val;
            result[jb * n_singles + ia] = 2.0 * val;
        }
    }

    Ok(result)
}

fn bytemuck_cast_f32(data: &[f32]) -> Vec<u8> {
    data.iter().flat_map(|x| x.to_ne_bytes()).collect()
}

/// Reinterpret a byte buffer as `f32` values.
///
/// The buffer comes from a GPU read-back and is therefore a plain `Vec<u8>`,
/// which — unlike a GPU buffer — carries **no** alignment guarantee.  The
/// real invariant is the `prefix.is_empty() && suffix.is_empty()` test below:
/// it is only satisfied when every byte of `data` belongs to a `f32` value.
///
/// A misaligned or wrongly sized read-back is an **error**: the previous
/// implementation returned an empty slice, which the caller silently turned
/// into an all-zero sTDA matrix (blocker A-M2 / A-H2).
fn bytemuck_cast_from_u8(data: &[u8]) -> Result<&[f32], String> {
    // SAFETY: `align_to` is sound for any element type: it splits `data`
    // into a prefix that cannot be reinterpreted as `f32`, a middle slice
    // whose alignment matches the alignment of `data`, and a suffix.  We
    // only use the middle slice when the prefix and the suffix are empty,
    // i.e. when the whole buffer is a sequence of `f32` values.  `f32` is a
    // plain-old-data type (no padding, no invalid bit patterns), so
    // reinterpreting the bytes is sound.  The lifetime of the returned
    // slice is tied to `data`, so no dangling reference can be produced.
    let (prefix, result, suffix) = unsafe { data.align_to::<f32>() };
    if prefix.is_empty() && suffix.is_empty() {
        Ok(result)
    } else {
        Err(format!(
            "GPU readback buffer is not 4-byte aligned: {} leading and {} trailing byte(s) \
             could not be reinterpreted as f32",
            prefix.len(),
            suffix.len()
        ))
    }
}

fn pack_dims(m: u32, k: u32, n: u32) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(16);
    bytes.extend_from_slice(&m.to_ne_bytes());
    bytes.extend_from_slice(&k.to_ne_bytes());
    bytes.extend_from_slice(&n.to_ne_bytes());
    bytes.extend_from_slice(&0u32.to_ne_bytes()); // padding
    bytes
}

/// WGSL shader for general matrix multiply: C = A × B
const MATMUL_SHADER: &str = r#"
struct Dims { M: u32, K: u32, N: u32, _pad: u32 }

@group(0) @binding(0) var<storage, read> a: array<f32>;
@group(0) @binding(1) var<storage, read> b: array<f32>;
@group(0) @binding(2) var<storage, read_write> c: array<f32>;
@group(0) @binding(3) var<uniform> dims: Dims;

@compute @workgroup_size(16, 16)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let row = gid.x;
    let col = gid.y;
    if row >= dims.M || col >= dims.N { return; }

    var sum: f32 = 0.0;
    for (var k: u32 = 0u; k < dims.K; k = k + 1u) {
        sum = sum + a[row * dims.K + k] * b[k * dims.N + col];
    }
    c[row * dims.N + col] = sum;
}
"#;

/// WGSL shader for transpose-multiply: C = A^T × B
const MATMUL_TRANSPOSE_SHADER: &str = r#"
struct Dims { K: u32, M: u32, N: u32, _pad: u32 }

@group(0) @binding(0) var<storage, read> a: array<f32>;
@group(0) @binding(1) var<storage, read> b: array<f32>;
@group(0) @binding(2) var<storage, read_write> c: array<f32>;
@group(0) @binding(3) var<uniform> dims: Dims;

@compute @workgroup_size(16, 16)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let row = gid.x;
    let col = gid.y;
    if row >= dims.M || col >= dims.N { return; }

    var sum: f32 = 0.0;
    for (var k: u32 = 0u; k < dims.K; k = k + 1u) {
        sum = sum + a[k * dims.M + row] * b[k * dims.N + col];
    }
    c[row * dims.N + col] = sum;
}
"#;

#[cfg(test)]
mod tests {
    use super::*;

    /// A-M2 regression test: a misaligned or wrongly sized GPU read-back
    /// must be an **error**, not a silently empty (and therefore wrong) sTDA
    /// matrix.
    #[test]
    fn test_bytemuck_cast_from_u8_rejects_bad_buffers() {
        let bytes = vec![0u8; 16];

        // A slice that starts one byte into the allocation is not 4-byte
        // aligned: the reinterpretation must be refused.
        let misaligned = &bytes[1..];
        assert!(
            bytemuck_cast_from_u8(misaligned).is_err(),
            "a misaligned read-back buffer must be an error"
        );

        // A length that is not a multiple of 4 leaves trailing bytes that
        // cannot be part of an f32 value.
        let partial = &bytes[..5];
        assert!(
            bytemuck_cast_from_u8(partial).is_err(),
            "a partial f32 read-back buffer must be an error"
        );

        // A correctly sized and aligned buffer is accepted.
        let aligned = bytemuck_cast_from_u8(&bytes)
            .expect("an aligned, correctly sized read-back buffer is accepted");
        assert_eq!(aligned.len(), 4);
    }
}
