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
        .unwrap_or_default();

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
        .unwrap_or_default();

    // Convert f32 result back to f64 with 2× scaling
    let a_off_f32: &[f32] = bytemuck_cast_from_u8(&a_off_bytes);
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

fn bytemuck_cast_from_u8(data: &[u8]) -> &[f32] {
    // Safety: data alignment is guaranteed by GPU buffer alignment
    let (prefix, result, suffix) = unsafe { data.align_to::<f32>() };
    if prefix.is_empty() && suffix.is_empty() {
        result
    } else {
        &[]
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
