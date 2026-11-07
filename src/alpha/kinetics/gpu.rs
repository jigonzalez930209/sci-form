//! GPU-accelerated HTST temperature sweeps.
//!
//! Offloads large batched temperature sweeps to GPU compute when the number
//! of temperature points exceeds a threshold (≥256). Falls back to CPU
//! (optionally rayon-parallelized) for smaller sweeps.
//!
//! The WGSL kernel evaluates the Eyring/Arrhenius expressions for every
//! temperature in parallel; the CPU path remains the reference and is used
//! below 256 points or without an adapter.
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

use super::{
    evaluate_htst_temperature_sweep, ElementaryRateResult, ElementaryStep, ThermodynamicState,
    BOLTZMANN_EV_PER_K, PLANCK_EV_S,
};
use crate::gpu::context::{
    bytes_to_f32_vec, ceil_div_u32, f32_slice_to_bytes, pack_uniform_values,
    ComputeBindingDescriptor, ComputeBindingKind, ComputeDispatchDescriptor, GpuContext,
    UniformValue,
};

/// Minimum sweep size to justify GPU dispatch.
const GPU_HTST_SWEEP_THRESHOLD: usize = 256;

/// GPU-accelerated HTST rate evaluation across a temperature grid.
///
/// The kernel evaluates the same Eyring/Arrhenius expressions as
/// [`super::evaluate_htst_rate`] for every temperature in parallel: forward
/// rate, equilibrium constant and (on the host) reverse rate. Validation runs
/// before the dispatch, so an invalid temperature or free energy returns the
/// same error as the CPU path. Systems below 256 points use the CPU path.
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

    // Validate exactly like the CPU path before touching the device.
    // Pressure is part of the thermodynamic state but not of the ideal-gas
    // Eyring expression used here, so it only has to be finite.
    if !pressure_bar.is_finite() {
        return Err("pressure must be finite".into());
    }
    if !step.activation_free_energy_ev.is_finite() || step.activation_free_energy_ev < 0.0 {
        return Err("activation free energy must be finite and non-negative".into());
    }
    if !step.reaction_free_energy_ev.is_finite() {
        return Err("reaction free energy must be finite".into());
    }
    for &t in temperatures_k {
        if !t.is_finite() || t <= 0.0 {
            return Err("temperature must be positive and finite".into());
        }
    }

    let n = temperatures_k.len();
    let params = pack_uniform_values(&[
        UniformValue::U32(n as u32),
        UniformValue::U32(step.prefactor_s_inv.is_some() as u32),
        UniformValue::F32(step.activation_free_energy_ev as f32),
        UniformValue::F32(step.reaction_free_energy_ev as f32),
        UniformValue::F32(step.prefactor_s_inv.unwrap_or(0.0) as f32),
        UniformValue::F32(BOLTZMANN_EV_PER_K as f32),
        UniformValue::F32(PLANCK_EV_S as f32),
        UniformValue::F32(pressure_bar as f32),
    ]);

    let descriptor = ComputeDispatchDescriptor {
        label: "htst temperature sweep".to_string(),
        shader_source: HTST_SWEEP_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [ceil_div_u32(n, 64), 1, 1],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "temperatures".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(
                    &temperatures_k.iter().map(|&t| t as f32).collect::<Vec<_>>(),
                ),
            },
            ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: params,
            },
            ComputeBindingDescriptor {
                label: "forward_rate".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: f32_slice_to_bytes(&vec![0.0f32; n]),
            },
            ComputeBindingDescriptor {
                label: "equilibrium".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: f32_slice_to_bytes(&vec![0.0f32; n]),
            },
        ],
    };

    let outputs = ctx.run_compute(&descriptor)?.outputs;
    if outputs.len() != 2 {
        return Err(format!(
            "HTST kernel returned {} read-back buffers, expected 2",
            outputs.len()
        ));
    }
    let forward = bytes_to_f32_vec(&outputs[0]);
    let equilibrium = bytes_to_f32_vec(&outputs[1]);
    if forward.len() != n || equilibrium.len() != n {
        return Err(format!(
            "HTST output size mismatch: forward {}, equilibrium {}, expected {n}",
            forward.len(),
            equilibrium.len()
        ));
    }

    Ok((0..n)
        .map(|i| {
            let fwd = forward[i] as f64;
            let eq = equilibrium[i] as f64;
            ElementaryRateResult {
                step_id: step.step_id.clone(),
                forward_rate_s_inv: fwd,
                reverse_rate_s_inv: if eq > 0.0 { fwd / eq } else { 0.0 },
                equilibrium_constant: eq,
                state: ThermodynamicState {
                    temperature_k: temperatures_k[i],
                    pressure_bar,
                },
            }
        })
        .collect())
}

/// WGSL source of the HTST temperature-sweep kernel.
///
/// One invocation per temperature; the prefactor is either the physical
/// `k_B T / h` (Eyring) or the per-step override carried in the uniform block.
pub const HTST_SWEEP_SHADER: &str = r#"
struct Params {
    n: u32,
    use_prefactor: u32,
    activation_ev: f32,
    reaction_ev: f32,
    prefactor: f32,
    kb_ev_k: f32,
    planck_ev_s: f32,
    pressure_bar: f32,
};

@group(0) @binding(0) var<storage, read> temperatures: array<f32>;
@group(0) @binding(1) var<uniform> params: Params;
@group(0) @binding(2) var<storage, read_write> forward_rate: array<f32>;
@group(0) @binding(3) var<storage, read_write> equilibrium: array<f32>;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= params.n) { return; }
    let t = temperatures[i];
    let kt = params.kb_ev_k * t;
    let prefactor = select(t * params.kb_ev_k / params.planck_ev_s, params.prefactor, params.use_prefactor != 0u);
    forward_rate[i] = prefactor * exp(-params.activation_ev / kt);
    equilibrium[i] = exp(-params.reaction_ev / kt);
}
"#;

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
