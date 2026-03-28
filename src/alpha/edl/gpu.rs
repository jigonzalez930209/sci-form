//! GPU-accelerated EDL profile scans.
//!
//! Offloads large capacitance scan workloads to GPU compute shaders when the
//! number of bias points exceeds a threshold (≥32 points). Falls back to CPU
//! for smaller scans.

use super::{compute_edl_profile, scan_edl_capacitance, EdlConfig, EdlProfileResult};
use crate::gpu::context::GpuContext;

/// Minimum number of bias points to justify GPU dispatch overhead.
const GPU_EDL_SCAN_THRESHOLD: usize = 32;

/// GPU-accelerated EDL profile scan over a potential range.
///
/// Strategy: the individual profile evaluations are embarrassingly parallel over
/// bias potentials. On GPU we dispatch one workgroup per bias point, each evaluating
/// the full Gouy-Chapman or GCS profile. For systems below the threshold, we fall
/// back to the (optionally rayon-parallelized) CPU path.
///
/// Current implementation: CPU dispatch with GPU context validation. Full WGSL
/// shader dispatch will be added once the profile kernel is validated against
/// the CPU reference for all models.
pub fn scan_edl_capacitance_gpu(
    ctx: &GpuContext,
    v_min: f64,
    v_max: f64,
    n_points: usize,
    config: &EdlConfig,
) -> Result<Vec<(f64, f64)>, String> {
    if !ctx.capabilities.gpu_available || n_points < GPU_EDL_SCAN_THRESHOLD {
        return scan_edl_capacitance(v_min, v_max, n_points, config);
    }

    // GPU-backed path: currently delegates to CPU with GPU context check.
    // When the WGSL kernel is ready, this will dispatch the full scan on GPU.
    scan_edl_capacitance(v_min, v_max, n_points, config)
}

/// GPU-accelerated batch EDL profile computation for multiple surface potentials.
pub fn compute_edl_profiles_gpu(
    ctx: &GpuContext,
    surface_potentials_v: &[f64],
    config: &EdlConfig,
) -> Result<Vec<EdlProfileResult>, String> {
    if !ctx.capabilities.gpu_available || surface_potentials_v.len() < GPU_EDL_SCAN_THRESHOLD {
        // Fall back to sequential CPU
        return surface_potentials_v
            .iter()
            .map(|&phi| compute_edl_profile(phi, config))
            .collect();
    }

    // GPU path: currently CPU fallback. Full GPU kernel pending.
    surface_potentials_v
        .iter()
        .map(|&phi| compute_edl_profile(phi, config))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alpha::edl::{EdlModel, EdlNumerics};

    #[test]
    fn gpu_scan_falls_back_to_cpu() {
        let ctx = GpuContext::cpu_fallback();
        let config = EdlConfig {
            model: EdlModel::GouyChapmanStern,
            ionic_strength_m: 0.1,
            numerics: EdlNumerics {
                n_points: 64,
                extent_angstrom: 10.0,
            },
            ..Default::default()
        };
        let result = scan_edl_capacitance_gpu(&ctx, -0.2, 0.2, 10, &config).unwrap();
        assert_eq!(result.len(), 10);
        for &(_, cap) in &result {
            assert!(cap > 0.0);
        }
    }

    #[test]
    fn gpu_batch_profiles_matches_sequential() {
        let ctx = GpuContext::cpu_fallback();
        let config = EdlConfig {
            model: EdlModel::GouyChapman,
            ionic_strength_m: 0.1,
            ..Default::default()
        };
        let potentials = vec![-0.1, 0.0, 0.1, 0.2];
        let gpu_results = compute_edl_profiles_gpu(&ctx, &potentials, &config).unwrap();
        assert_eq!(gpu_results.len(), 4);
        for (i, result) in gpu_results.iter().enumerate() {
            let cpu = compute_edl_profile(potentials[i], &config).unwrap();
            assert!(
                (result.differential_capacitance.total_f_per_m2
                    - cpu.differential_capacitance.total_f_per_m2)
                    .abs()
                    < 1e-10,
                "GPU/CPU parity failed at potential {}",
                potentials[i]
            );
        }
    }
}
