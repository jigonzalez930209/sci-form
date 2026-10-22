//! GPU-accelerated HTST temperature sweeps.
//!
//! Offloads large batched temperature sweeps to GPU compute when the number
//! of temperature points exceeds a threshold (≥256). Falls back to CPU
//! (optionally rayon-parallelized) for smaller sweeps.
//!
//! NOTE: this is currently a **CPU fallback only** — the WGSL kernel it
//! advertises does not exist yet, so the "GPU" path executes exactly the same
//! CPU code and the parity test compares CPU against CPU. The audit
//! (reports/auditoria-total.md, P2-32) flagged this; it is kept as a stable
//! call site for the kernel to be added.
//!
//! # Stability tier
//!
//! **Alpha** — gated behind `alpha-kinetics` (plus `experimental-gpu` for the
//! GPU context); the validation gate is `tests/alpha_kinetics_validation.rs`
//! (F87).
//!
//! # Units
//!
//! Temperatures in **K**, pressure in **bar**, activation barriers in
//! **kcal/mol** and the returned rates in **s⁻¹**.

use super::{evaluate_htst_temperature_sweep, ElementaryRateResult, ElementaryStep};
use crate::gpu::context::GpuContext;

/// Minimum sweep size to justify GPU dispatch.
const GPU_HTST_SWEEP_THRESHOLD: usize = 256;

/// GPU-accelerated HTST rate evaluation across a temperature grid.
///
/// Strategy: each temperature point is independent — the exponential
/// Arrhenius/Eyring evaluation is trivially parallel on GPU. Current
/// implementation: CPU fallback. Full WGSL kernel pending.
///
/// # Errors
///
/// Propagates the errors of [`evaluate_htst_temperature_sweep`]: a
/// non-positive or non-finite temperature, a negative or non-finite
/// activation free energy, or a non-finite reaction free energy.
pub fn evaluate_htst_sweep_gpu(
    ctx: &GpuContext,
    step: &ElementaryStep,
    temperatures_k: &[f64],
    pressure_bar: f64,
) -> Result<Vec<ElementaryRateResult>, String> {
    if !ctx.capabilities.gpu_available || temperatures_k.len() < GPU_HTST_SWEEP_THRESHOLD {
        return evaluate_htst_temperature_sweep(step, temperatures_k, pressure_bar);
    }

    // GPU path: currently CPU fallback. Full WGSL kernel pending.
    evaluate_htst_temperature_sweep(step, temperatures_k, pressure_bar)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gpu_sweep_falls_back_to_cpu() {
        let ctx = GpuContext::cpu_fallback();
        let step = ElementaryStep {
            step_id: "test".into(),
            activation_free_energy_ev: 0.5,
            reaction_free_energy_ev: -0.2,
            prefactor_s_inv: None,
        };
        let temps: Vec<f64> = (300..=800).step_by(50).map(|t| t as f64).collect();
        let gpu_result = evaluate_htst_sweep_gpu(&ctx, &step, &temps, 1.0).unwrap();
        let cpu_result = evaluate_htst_temperature_sweep(&step, &temps, 1.0).unwrap();
        assert_eq!(gpu_result.len(), cpu_result.len());
        for (g, c) in gpu_result.iter().zip(&cpu_result) {
            assert!(
                (g.forward_rate_s_inv - c.forward_rate_s_inv).abs() < 1e-10,
                "GPU/CPU parity: {} vs {}",
                g.forward_rate_s_inv,
                c.forward_rate_s_inv
            );
        }
    }
}
