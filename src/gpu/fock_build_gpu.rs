//! GPU-accelerated Fock matrix construction.
//!
//! Computes the two-electron part G(P) of the Fock matrix on the GPU:
//!   G(μν) = Σ_{λσ} P_{λσ} [(μν|λσ) - 0.5·(μλ|νσ)]
//!   F = H_core + G(P)
//!
//! Each GPU thread handles one (μ,ν) matrix element, reading the full
//! density matrix and ERI tensor to compute the Coulomb-exchange contribution.

use super::context::{
    bytes_to_f64_vec_from_f32, ceil_div_u32, f32_slice_to_bytes, pack_uniform_values,
    ComputeBindingDescriptor, ComputeBindingKind, ComputeDispatchDescriptor, GpuContext,
    UniformValue,
};

/// Minimum matrix dimension to justify GPU dispatch.
const GPU_DISPATCH_THRESHOLD: usize = 4;

/// Build the Fock matrix on the GPU.
///
/// Inputs:
/// - `h_core`: Core Hamiltonian (N×N), row-major flattened f64
/// - `density`: Density matrix P (N×N), row-major flattened f64
/// - `eris`: Two-electron integrals (N⁴), flattened f64
/// - `n_basis`: Number of basis functions
///
/// Returns the Fock matrix F = H + G(P) as flattened f64 (N×N).
pub fn build_fock_gpu(
    ctx: &GpuContext,
    h_core: &[f64],
    density: &[f64],
    eris: &[f64],
    n_basis: usize,
) -> Result<Vec<f64>, String> {
    if n_basis < GPU_DISPATCH_THRESHOLD {
        return Err("Basis too small for GPU dispatch".to_string());
    }

    let n2 = n_basis * n_basis;
    if h_core.len() != n2 || density.len() != n2 {
        return Err("Matrix dimension mismatch".to_string());
    }

    // Pack matrices as f32
    let h_core_f32: Vec<f32> = h_core.iter().map(|v| *v as f32).collect();
    let density_f32: Vec<f32> = density.iter().map(|v| *v as f32).collect();
    let eris_f32: Vec<f32> = eris.iter().map(|v| *v as f32).collect();

    // Params: [n_basis: u32, pad: u32, pad: u32, pad: u32] = 16 bytes
    let params = pack_uniform_values(&[
        UniformValue::U32(n_basis as u32),
        UniformValue::U32(0),
        UniformValue::U32(0),
        UniformValue::U32(0),
    ]);

    // Output: Fock matrix (N×N f32)
    let output_seed = vec![0.0f32; n2];

    let wg_size = 16u32;
    let wg_x = ceil_div_u32(n_basis, wg_size);
    let wg_y = wg_x;

    let descriptor = ComputeDispatchDescriptor {
        label: "fock matrix build".to_string(),
        shader_source: FOCK_BUILD_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [wg_x, wg_y, 1],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "h_core".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&h_core_f32),
            },
            ComputeBindingDescriptor {
                label: "density".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&density_f32),
            },
            ComputeBindingDescriptor {
                label: "eris".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&eris_f32),
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
        .ok_or("No output from Fock build kernel")?;

    if bytes.len() != n2 * 4 {
        return Err(format!(
            "Fock output size mismatch: expected {}, got {}",
            n2 * 4,
            bytes.len()
        ));
    }

    let fock = bytes_to_f64_vec_from_f32(&bytes);

    Ok(fock)
}

/// WGSL compute shader for Fock matrix construction.
///
/// F(μ,ν) = H_core(μ,ν) + Σ_{λσ} P(λ,σ) [(μν|λσ) - 0.5·(μλ|νσ)]
///
/// Workgroup size: (16, 16, 1) = 256 threads.
/// Each thread computes one F(μ,ν) element.
pub const FOCK_BUILD_SHADER: &str = r#"
struct Params {
    n_basis: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
};

@group(0) @binding(0) var<storage, read> h_core: array<f32>;
@group(0) @binding(1) var<storage, read> density: array<f32>;
@group(0) @binding(2) var<storage, read> eris: array<f32>;
@group(0) @binding(3) var<uniform> params: Params;
@group(0) @binding(4) var<storage, read_write> fock: array<f32>;

@compute @workgroup_size(16, 16, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let mu = gid.x;
    let nu = gid.y;
    let n = params.n_basis;

    if (mu >= n || nu >= n) {
        return;
    }

    let n2 = n * n;
    var g_mn: f32 = 0.0;

    // G(μ,ν) = Σ_{λσ} P(λ,σ) · [(μν|λσ) - 0.5·(μλ|νσ)]
    for (var lam: u32 = 0u; lam < n; lam = lam + 1u) {
        for (var sig: u32 = 0u; sig < n; sig = sig + 1u) {
            let p_ls = density[lam * n + sig];

            // Coulomb: (μν|λσ)
            let j_idx = mu * n * n2 + nu * n2 + lam * n + sig;
            let j_val = eris[j_idx];

            // Exchange: (μλ|νσ)
            let k_idx = mu * n * n2 + lam * n2 + nu * n + sig;
            let k_val = eris[k_idx];

            g_mn += p_ls * (j_val - 0.5 * k_val);
        }
    }

    fock[mu * n + nu] = h_core[mu * n + nu] + g_mn;
}
"#;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_f32_slice_to_bytes() {
        let data = vec![1.0f32, 2.0, 3.0];
        let bytes = f32_slice_to_bytes(&data);
        assert_eq!(bytes.len(), 12);
    }

    #[test]
    fn test_gpu_threshold() {
        let ctx = GpuContext::cpu_fallback();
        let n = 2;
        let h = vec![0.0f64; n * n];
        let d = vec![0.0f64; n * n];
        let e = vec![0.0f64; n * n * n * n];
        let result = build_fock_gpu(&ctx, &h, &d, &e, n);
        assert!(result.is_err());
    }

    #[test]
    fn test_dimension_mismatch() {
        let ctx = GpuContext::cpu_fallback();
        // n=5 passes threshold but matrices have wrong sizes
        let result = build_fock_gpu(&ctx, &[0.0; 25], &[0.0; 16], &[0.0; 625], 5);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("mismatch"));
    }
}
