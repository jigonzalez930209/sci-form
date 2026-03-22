//! GPU-accelerated PM3 NDDO Fock matrix build.
//!
//! Offloads the two-center Coulomb G-matrix construction to GPU when
//! the `experimental-gpu` feature is enabled and the basis is large enough.

#[cfg(feature = "experimental-gpu")]
use crate::gpu::context::{
    bytes_to_f64_vec_from_f32, ceil_div_u32, f32_slice_to_bytes, pack_uniform_values,
    ComputeBindingDescriptor, ComputeBindingKind, ComputeDispatchDescriptor, GpuContext,
    UniformValue,
};

/// Minimum basis size to justify GPU dispatch for PM3.
#[cfg(feature = "experimental-gpu")]
const GPU_DISPATCH_THRESHOLD: usize = 16;

/// Build the PM3 NDDO G-matrix on GPU.
///
/// Computes the two-center electron-electron Coulomb contribution to the
/// Fock matrix diagonal: G_ii += Σ_{B≠A_i} P_BB · γ_{A_i,B}
///
/// Falls back to error if basis is too small or GPU unavailable.
#[cfg(feature = "experimental-gpu")]
pub fn build_pm3_g_matrix_gpu(
    ctx: &GpuContext,
    density_diag: &[f64],
    atom_of_basis: &[u32],
    gamma_ab: &[f64],
    n_basis: usize,
    n_atoms: usize,
) -> Result<Vec<f64>, String> {
    if n_basis < GPU_DISPATCH_THRESHOLD {
        return Err("Basis too small for GPU dispatch".to_string());
    }

    let dens_f32: Vec<f32> = density_diag.iter().map(|v| *v as f32).collect();
    let atom_f32: Vec<f32> = atom_of_basis.iter().map(|v| *v as f32).collect();
    let gamma_f32: Vec<f32> = gamma_ab.iter().map(|v| *v as f32).collect();

    let params = pack_uniform_values(&[
        UniformValue::U32(n_basis as u32),
        UniformValue::U32(n_atoms as u32),
        UniformValue::U32(0),
        UniformValue::U32(0),
    ]);

    let output_seed = vec![0.0f32; n_basis];
    let wg_size = 64u32;
    let wg_x = ceil_div_u32(n_basis, wg_size);

    let descriptor = ComputeDispatchDescriptor {
        label: "pm3 two-center coulomb".to_string(),
        shader_source: PM3_G_MATRIX_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [wg_x, 1, 1],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "density_diag".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&dens_f32),
            },
            ComputeBindingDescriptor {
                label: "atom_of_basis".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&atom_f32),
            },
            ComputeBindingDescriptor {
                label: "gamma".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&gamma_f32),
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
    let bytes = result
        .outputs
        .pop()
        .ok_or("No output from PM3 G-matrix kernel")?;

    Ok(bytes_to_f64_vec_from_f32(&bytes))
}

/// WGSL shader for PM3 two-center Coulomb contribution.
///
/// Each thread computes G_diag[i] = Σ_{B≠atom(i)} P_B · γ_{atom(i),B}
/// where P_B = Σ_{k ∈ B} density_diag[k].
#[cfg(feature = "experimental-gpu")]
const PM3_G_MATRIX_SHADER: &str = r#"
struct Params {
    n_basis: u32,
    n_atoms: u32,
    _pad0: u32,
    _pad1: u32,
};

@group(0) @binding(0) var<storage, read> density_diag: array<f32>;
@group(0) @binding(1) var<storage, read> atom_of_basis: array<f32>;
@group(0) @binding(2) var<storage, read> gamma: array<f32>;
@group(0) @binding(3) var<uniform> params: Params;
@group(0) @binding(4) var<storage, read_write> output: array<f32>;

@compute @workgroup_size(64, 1, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= params.n_basis) { return; }

    let atom_i = u32(atom_of_basis[i]);
    var g_val: f32 = 0.0;

    // Sum density on atom B, multiply by gamma(atom_i, B)
    for (var b: u32 = 0u; b < params.n_atoms; b = b + 1u) {
        if (b == atom_i) { continue; }
        var p_b: f32 = 0.0;
        for (var k: u32 = 0u; k < params.n_basis; k = k + 1u) {
            if (u32(atom_of_basis[k]) == b) {
                p_b += density_diag[k];
            }
        }
        g_val += p_b * gamma[atom_i * params.n_atoms + b];
    }

    output[i] = g_val;
}
"#;

#[cfg(test)]
mod tests {
    #[test]
    fn test_pm3_gpu_module_compiles() {
        // Verify the module compiles with or without gpu feature
        assert!(true);
    }
}
