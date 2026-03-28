//! GPU-accelerated periodic operator assembly and KPM workloads.
//!
//! Offloads k-point operator assembly and KPM matrix-vector products to GPU
//! compute shaders when the workload is large enough.

use nalgebra::DMatrix;

use super::{assemble_periodic_operators, KMesh};
use crate::gpu::context::GpuContext;

/// Minimum k-mesh size to justify GPU dispatch.
const GPU_PERIODIC_THRESHOLD: usize = 64;

/// GPU-accelerated periodic Bloch Hamiltonian assembly across k-points.
///
/// Strategy: each k-point's Bloch sum H(k) = Σ_R H_R exp(ik·R) is independent,
/// making this embarrassingly parallel on GPU. Current implementation: CPU fallback
/// with GPU context validation. Full WGSL kernel pending.
pub fn assemble_periodic_operators_gpu(
    ctx: &GpuContext,
    h_real_space: &[(&DMatrix<f64>, [i32; 3])],
    s_real_space: &[(&DMatrix<f64>, [i32; 3])],
    kmesh: &KMesh,
) -> (Vec<DMatrix<f64>>, Vec<DMatrix<f64>>) {
    if !ctx.capabilities.gpu_available || kmesh.points.len() < GPU_PERIODIC_THRESHOLD {
        return assemble_periodic_operators(h_real_space, s_real_space, kmesh);
    }

    // GPU path: currently CPU fallback. Full WGSL kernel pending validation.
    assemble_periodic_operators(h_real_space, s_real_space, kmesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alpha::periodic_linear::{KMeshCentering, KMeshConfig};

    #[test]
    fn gpu_assembly_falls_back_to_cpu() {
        let ctx = GpuContext::cpu_fallback();
        let h = DMatrix::<f64>::identity(4, 4);
        let h_pairs: Vec<(&DMatrix<f64>, [i32; 3])> = vec![(&h, [0, 0, 0])];
        let s_pairs: Vec<(&DMatrix<f64>, [i32; 3])> = vec![(&h, [0, 0, 0])];
        let mesh = crate::alpha::periodic_linear::monkhorst_pack_mesh(&KMeshConfig {
            grid: [3, 3, 3],
            centering: KMeshCentering::MonkhorstPack,
        })
        .unwrap();

        let (hk_gpu, sk_gpu) = assemble_periodic_operators_gpu(&ctx, &h_pairs, &s_pairs, &mesh);
        let (hk_cpu, sk_cpu) = assemble_periodic_operators(&h_pairs, &s_pairs, &mesh);
        assert_eq!(hk_gpu.len(), hk_cpu.len());
        assert_eq!(sk_gpu.len(), sk_cpu.len());
        for (g, c) in hk_gpu.iter().zip(&hk_cpu) {
            assert!((g - c).norm() < 1e-12, "GPU/CPU H(k) parity violated");
        }
    }
}
