//! GPU-accelerated SCC-DFTB gamma-matrix construction.
//!
//! Builds the (undamped, inverse-distance) second-order perturbation matrix
//! **Γ** used by the self-consistent-charge DFTB SCC cycle.
//!
//! # Contract with callers
//!
//! The CPU implementation always works; this module is the optional
//! accelerator. The kernel refuses systems with fewer than two atoms.
//!
//! # Units
//!
//! `positions_bohr` are in **Bohr**, so the off-diagonal entries `1/r_ij` are
//! in **Bohr⁻¹**; the diagonal entries are the per-atom chemical hardness
//! `eta` (in the unit the caller supplies — the SCC DFTB codes of this crate
//! work in **eV**).

use nalgebra::DMatrix;

use super::context::{
    bytes_to_f64_vec_from_f32, ceil_div_u32, f32_slice_to_bytes, pack_uniform_values,
    pack_vec3_positions_f32, ComputeBindingDescriptor, ComputeBindingKind,
    ComputeDispatchDescriptor, GpuContext, UniformValue,
};

const GPU_DISPATCH_THRESHOLD: usize = 2;

/// Builds the inverse-distance gamma matrix on the GPU.
///
/// `eta` holds the per-atom chemical hardness and `positions_bohr` the atom
/// coordinates in **Bohr**. The returned matrix has `1/r_ij` (in
/// **Bohr⁻¹**) off the diagonal and `eta_i` on it.
///
/// # Errors
///
/// - `"Matrix too small for GPU dispatch"` when fewer than two atoms are
///   supplied.
/// - `"eta/position length mismatch"` when the two slices differ in length.
/// - `"No output from gamma kernel"` when the dispatch returned no read-back
///   payload.
/// - `"Gamma output size mismatch: expected {n}, got {m}"` when the kernel
///   produced the wrong number of bytes.
/// - Any error from [`GpuContext::run_compute`].
pub fn build_gamma_gpu(
    ctx: &GpuContext,
    eta: &[f64],
    positions_bohr: &[[f64; 3]],
) -> Result<DMatrix<f64>, String> {
    let n = eta.len();
    if n < GPU_DISPATCH_THRESHOLD {
        return Err("Matrix too small for GPU dispatch".to_string());
    }
    if positions_bohr.len() != n {
        return Err("eta/position length mismatch".to_string());
    }

    let eta_f32: Vec<f32> = eta.iter().map(|value| *value as f32).collect();
    let params = pack_uniform_values(&[
        UniformValue::U32(n as u32),
        UniformValue::U32(0),
        UniformValue::U32(0),
        UniformValue::U32(0),
    ]);

    let output_seed = vec![0.0f32; n * n];
    let wg_size = 16u32;
    let wg_x = ceil_div_u32(n, wg_size);
    let wg_y = wg_x;

    let descriptor = ComputeDispatchDescriptor {
        label: "gamma matrix".to_string(),
        shader_source: GAMMA_MATRIX_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [wg_x, wg_y, 1],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "eta".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&eta_f32),
            },
            ComputeBindingDescriptor {
                label: "positions".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: pack_vec3_positions_f32(positions_bohr),
            },
            ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: params,
            },
            ComputeBindingDescriptor {
                label: "output".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: f32_slice_to_bytes(&output_seed),
            },
        ],
    };

    let mut result = ctx.run_compute(&descriptor)?;
    let bytes = result.outputs.pop().ok_or("No output from gamma kernel")?;
    if bytes.len() != n * n * 4 {
        return Err(format!(
            "Gamma output size mismatch: expected {}, got {}",
            n * n * 4,
            bytes.len()
        ));
    }

    let values = bytes_to_f64_vec_from_f32(&bytes);
    Ok(DMatrix::from_row_slice(n, n, &values))
}

/// WGSL source of the inverse-distance gamma-matrix kernel.
///
/// One invocation per matrix element (`@workgroup_size(16, 16, 1)`); the
/// kernel writes `1/r_ij` (in **Bohr⁻¹**) for `i ≠ j` and `eta_i` on the
/// diagonal. The uniform block carries `n_atoms` and three padding words.
pub const GAMMA_MATRIX_SHADER: &str = r#"
struct AtomPos {
    x: f32, y: f32, z: f32, _pad: f32,
};

struct Params {
    n_atoms: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
};

@group(0) @binding(0) var<storage, read> eta: array<f32>;
@group(0) @binding(1) var<storage, read> positions: array<AtomPos>;
@group(0) @binding(2) var<uniform> params: Params;
@group(0) @binding(3) var<storage, read_write> output: array<f32>;

@compute @workgroup_size(16, 16, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let a = gid.x;
    let b = gid.y;
    let n = params.n_atoms;

    if (a >= n || b >= n) {
        return;
    }

    if (a == b) {
        output[a * n + b] = eta[a];
        return;
    }

    let pa = positions[a];
    let pb = positions[b];
    let dx = pa.x - pb.x;
    let dy = pa.y - pb.y;
    let dz = pa.z - pb.z;
    let r = sqrt(dx * dx + dy * dy + dz * dz);

    output[a * n + b] = 1.0 / r;
}
"#;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gamma_gpu_threshold() {
        let ctx = GpuContext::cpu_fallback();
        let result = build_gamma_gpu(&ctx, &[0.5], &[[0.0, 0.0, 0.0]]);
        assert!(result.is_err());
    }
}
