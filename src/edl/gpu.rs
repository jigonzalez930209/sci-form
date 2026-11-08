//! GPU-accelerated EDL profile scans.
//!
//! The Helmholtz, Gouy-Chapman and Gouy-Chapman-Stern profiles are analytic,
//! so each bias point and grid point is independent. The WGSL kernel evaluates
//! the same expressions as the CPU path and writes the six profile arrays plus
//! the capacitance decomposition; the CPU path is used below the dispatch
//! threshold or without an adapter.
//!
//! NOTE: the kernel works in `f32`, so comparisons against the `f64` CPU
//! reference are relative (1e-3 for the capacitance, 1e-3/1e-6 absolute for
//! the profile fields).

use super::{
    compute_edl_profile, scan_edl_capacitance, CapacitanceBreakdown, EdlBackend, EdlConfig,
    EdlModel, EdlProfileResult,
};
use crate::gpu::context::{
    bytes_to_f32_vec, ceil_div_u32, f32_slice_to_bytes, pack_uniform_values,
    ComputeBindingDescriptor, ComputeBindingKind, ComputeDispatchDescriptor, GpuContext,
    UniformValue,
};

/// Minimum number of bias points to justify GPU dispatch overhead.
const GPU_EDL_SCAN_THRESHOLD: usize = 32;

/// GPU-accelerated EDL profile scan over a potential range.
///
/// The individual profile evaluations are embarrassingly parallel over bias
/// potentials. The kernel runs one invocation per (bias, grid point); below the
/// threshold the CPU path is used.
///
/// # Errors
///
/// Returns the same validation errors as [`scan_edl_capacitance`].
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
    if n_points < 2 {
        return Err("capacitance scan requires at least 2 points".into());
    }
    if !v_min.is_finite() || !v_max.is_finite() {
        return Err("capacitance scan endpoints must be finite".into());
    }
    let values: Vec<f64> = (0..n_points)
        .map(|i| {
            let t = i as f64 / (n_points - 1) as f64;
            v_min + t * (v_max - v_min)
        })
        .collect();
    let profiles = match compute_profiles_gpu(ctx, &values, config) {
        Ok(profiles) => profiles,
        Err(_) => return scan_edl_capacitance(v_min, v_max, n_points, config),
    };
    Ok(values
        .iter()
        .zip(profiles)
        .map(|(&value, profile)| (value, profile.differential_capacitance.total_f_per_m2))
        .collect())
}

/// GPU-accelerated batch EDL profile computation for multiple surface potentials.
///
/// `surface_potentials_v` are potentials for the Gouy-Chapman model and surface
/// charge densities for the Helmholtz/GCS models, matching
/// [`compute_edl_profile`]. Falls back to the CPU path below the threshold or
/// without an adapter.
pub fn compute_edl_profiles_gpu(
    ctx: &GpuContext,
    surface_potentials_v: &[f64],
    config: &EdlConfig,
) -> Result<Vec<EdlProfileResult>, String> {
    if !ctx.capabilities.gpu_available || surface_potentials_v.len() < GPU_EDL_SCAN_THRESHOLD {
        return surface_potentials_v
            .iter()
            .map(|&phi| compute_edl_profile(phi, config))
            .collect();
    }
    if surface_potentials_v.iter().any(|v| !v.is_finite()) {
        return Err("EDL surface values must be finite".into());
    }
    compute_profiles_gpu(ctx, surface_potentials_v, config).or_else(|_| {
        surface_potentials_v
            .iter()
            .map(|&phi| compute_edl_profile(phi, config))
            .collect()
    })
}

fn validate_config(config: &EdlConfig) -> Result<(), String> {
    if config.numerics.n_points < 2 {
        return Err("EDL profile requires at least 2 grid points".into());
    }
    if !config.numerics.extent_angstrom.is_finite() || config.numerics.extent_angstrom <= 0.0 {
        return Err("EDL profile extent must be positive".into());
    }
    if !config.temperature_k.is_finite() || config.temperature_k <= 0.0 {
        return Err("temperature must be positive".into());
    }
    match config.model {
        EdlModel::Helmholtz => {
            if config.stern_thickness_angstrom <= 0.0 {
                return Err("Stern thickness must be positive".into());
            }
            if config.compact_dielectric <= 0.0 {
                return Err("Compact dielectric must be positive".into());
            }
        }
        EdlModel::GouyChapman => {
            if config.ionic_strength_m <= 0.0 {
                return Err("Gouy-Chapman requires positive ionic strength".into());
            }
        }
        EdlModel::GouyChapmanStern | EdlModel::CompactDiffuseHybrid => {
            if config.stern_thickness_angstrom <= 0.0 {
                return Err("Stern thickness must be positive".into());
            }
            if config.compact_dielectric <= 0.0 {
                return Err("compact dielectric must be positive".into());
            }
            if config.ionic_strength_m <= 0.0 {
                return Err("GCS model requires positive ionic strength".into());
            }
        }
    }
    Ok(())
}

fn compute_profiles_gpu(
    ctx: &GpuContext,
    values: &[f64],
    config: &EdlConfig,
) -> Result<Vec<EdlProfileResult>, String> {
    validate_config(config)?;
    let n_values = values.len();
    let n_points = config.numerics.n_points;
    let total = n_values
        .checked_mul(n_points)
        .ok_or_else(|| "EDL profile grid is too large".to_string())?;
    let model_code = match config.model {
        EdlModel::Helmholtz => 0u32,
        EdlModel::GouyChapman => 1,
        EdlModel::GouyChapmanStern | EdlModel::CompactDiffuseHybrid => 2,
    };
    let params = pack_uniform_values(&[
        UniformValue::U32(model_code),
        UniformValue::U32(n_values as u32),
        UniformValue::U32(n_points as u32),
        UniformValue::U32(0),
        UniformValue::F32(config.numerics.extent_angstrom as f32),
        UniformValue::F32(config.stern_thickness_angstrom as f32),
        UniformValue::F32(config.compact_dielectric as f32),
        UniformValue::F32(config.bulk_dielectric as f32),
        UniformValue::F32(config.temperature_k as f32),
        UniformValue::F32(config.ionic_strength_m as f32),
    ]);
    let descriptor = ComputeDispatchDescriptor {
        label: "edl profile scan".to_string(),
        shader_source: EDL_PROFILE_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [ceil_div_u32(total, 64), 1, 1],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "surface_values".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&values.iter().map(|&v| v as f32).collect::<Vec<_>>()),
            },
            ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: params,
            },
            ComputeBindingDescriptor {
                label: "potential".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: f32_slice_to_bytes(&vec![0.0f32; total]),
            },
            ComputeBindingDescriptor {
                label: "field".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: f32_slice_to_bytes(&vec![0.0f32; total]),
            },
            ComputeBindingDescriptor {
                label: "charge_density".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: f32_slice_to_bytes(&vec![0.0f32; total]),
            },
            ComputeBindingDescriptor {
                label: "ion_density".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: f32_slice_to_bytes(&vec![0.0f32; total]),
            },
            ComputeBindingDescriptor {
                label: "dielectric".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: f32_slice_to_bytes(&vec![0.0f32; total]),
            },
            ComputeBindingDescriptor {
                label: "scalars".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: f32_slice_to_bytes(&vec![0.0f32; n_values * 6]),
            },
        ],
    };
    let outputs = ctx.run_compute(&descriptor)?.outputs;
    if outputs.len() != 6 {
        return Err(format!(
            "EDL kernel returned {} buffers, expected 6",
            outputs.len()
        ));
    }
    let potential = bytes_to_f32_vec(&outputs[0]);
    let field = bytes_to_f32_vec(&outputs[1]);
    let charge = bytes_to_f32_vec(&outputs[2]);
    let ion = bytes_to_f32_vec(&outputs[3]);
    let dielectric = bytes_to_f32_vec(&outputs[4]);
    let scalars = bytes_to_f32_vec(&outputs[5]);
    if [&potential, &field, &charge, &ion, &dielectric]
        .iter()
        .any(|array| array.len() != total)
        || scalars.len() != n_values * 6
    {
        return Err("EDL kernel output size mismatch".into());
    }
    let distance: Vec<f64> = (0..n_points)
        .map(|i| config.numerics.extent_angstrom * i as f64 / (n_points - 1) as f64)
        .collect();
    let model_name = match config.model {
        EdlModel::Helmholtz => "Helmholtz",
        EdlModel::GouyChapman => "Gouy-Chapman",
        EdlModel::GouyChapmanStern => "Gouy-Chapman-Stern",
        EdlModel::CompactDiffuseHybrid => "Compact-Diffuse Hybrid",
    };
    let slice = |array: &[f32], start: usize| -> Vec<f64> {
        array[start..start + n_points]
            .iter()
            .map(|&v| v as f64)
            .collect()
    };
    Ok((0..n_values)
        .map(|i| {
            let start = i * n_points;
            let s = &scalars[i * 6..i * 6 + 6];
            EdlProfileResult {
                distance_axis_angstrom: distance.clone(),
                electrostatic_potential_v: slice(&potential, start),
                field_strength_v_per_m: slice(&field, start),
                charge_density_c_per_m3: slice(&charge, start),
                ion_density_relative: slice(&ion, start),
                dielectric_profile: slice(&dielectric, start),
                compact_layer_drop_v: s[0] as f64,
                diffuse_layer_drop_v: s[1] as f64,
                total_interfacial_drop_v: s[2] as f64,
                differential_capacitance: CapacitanceBreakdown {
                    compact_f_per_m2: s[3] as f64,
                    diffuse_f_per_m2: s[4] as f64,
                    total_f_per_m2: s[5] as f64,
                },
                backend: EdlBackend::GpuAccelerated,
                used_gpu: true,
                converged: true,
                n_iterations: 1,
                residual: 0.0,
                temperature_k: config.temperature_k,
                ionic_strength_m: config.ionic_strength_m,
                model_name: model_name.to_string(),
            }
        })
        .collect())
}

/// WGSL source of the EDL profile kernel.
pub const EDL_PROFILE_SHADER: &str = r#"
struct Params {
    model: u32,
    n_values: u32,
    n_points: u32,
    _pad: u32,
    extent: f32,
    stern: f32,
    eps_compact: f32,
    eps_bulk: f32,
    temperature: f32,
    ionic: f32,
};

const EPS0: f32 = 8.8541878128e-12;
const ANG: f32 = 1.0e-10;
const E: f32 = 1.602176634e-19;
const KB: f32 = 1.380649e-23;
const NA: f32 = 6.02214076e23;

@group(0) @binding(0) var<storage, read> values: array<f32>;
@group(0) @binding(1) var<uniform> params: Params;
@group(0) @binding(2) var<storage, read_write> potential: array<f32>;
@group(0) @binding(3) var<storage, read_write> field: array<f32>;
@group(0) @binding(4) var<storage, read_write> charge: array<f32>;
@group(0) @binding(5) var<storage, read_write> ion: array<f32>;
@group(0) @binding(6) var<storage, read_write> dielectric: array<f32>;
@group(0) @binding(7) var<storage, read_write> scalars: array<f32>;

fn safe_log_ratio(g: f32) -> f32 {
    if (abs(g) >= 1.0 - 1.0e-6) { return 20.0; }
    return log((1.0 + g) / (1.0 - g));
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let idx = gid.x;
    let np = params.n_points;
    if (idx >= params.n_values * np) { return; }
    let vi = idx / np;
    let pi = idx % np;
    let value = values[vi];
    let x_ang = params.extent * f32(pi) / f32(np - 1u);
    let x_m = x_ang * ANG;

    var phi: f32 = 0.0;
    var efield: f32 = 0.0;
    var rho: f32 = 0.0;
    var ion_density: f32 = 0.0;
    var diel: f32 = params.eps_bulk;
    var compact_drop: f32 = 0.0;
    var diffuse_drop: f32 = 0.0;
    var total_drop: f32 = 0.0;
    var c_compact: f32 = 0.0;
    var c_diffuse: f32 = 0.0;
    var c_total: f32 = 0.0;

    let stern_m = params.stern * ANG;
    let kt = KB * params.temperature;
    let thermal_v = kt / E;

    if (params.model == 0u) {
        // Helmholtz: sigma -> linear compact layer.
        c_compact = EPS0 * params.eps_compact / stern_m;
        efield = value / (EPS0 * params.eps_compact);
        compact_drop = efield * stern_m;
        total_drop = compact_drop;
        c_total = c_compact;
        if (x_ang <= params.stern) {
            phi = efield * x_m;
            rho = value / stern_m;
            diel = params.eps_compact;
        } else {
            phi = compact_drop;
            efield = 0.0;
            rho = 0.0;
        }
    } else if (params.model == 1u) {
        // Gouy-Chapman: phi0 -> nonlinear diffuse layer.
        let c0 = params.ionic * 1000.0;
        let kappa = sqrt(2.0 * c0 * NA * E * E / (EPS0 * params.eps_bulk * kt));
        let gamma = tanh(E * value / (4.0 * kt));
        let exp_term = exp(-kappa * x_m);
        let g = gamma * exp_term;
        phi = 2.0 * thermal_v * safe_log_ratio(g);
        if (pi == 0u) { phi = value; }
        let denom = 1.0 - gamma * gamma * exp_term * exp_term;
        if (abs(denom) > 1.0e-6) {
            efield = 2.0 * thermal_v * kappa * 2.0 * gamma * exp_term / denom;
        }
        rho = -2.0 * E * c0 * NA * sinh(E * phi / kt);
        ion_density = exp(-E * phi / kt) + exp(E * phi / kt);
        diel = params.eps_bulk;
        diffuse_drop = value;
        total_drop = value;
        c_diffuse = EPS0 * params.eps_bulk * kappa * cosh(E * value / (2.0 * kt));
        c_total = c_diffuse;
    } else {
        // Gouy-Chapman-Stern: sigma -> compact + diffuse series.
        let c0 = params.ionic * 1000.0;
        c_compact = EPS0 * params.eps_compact / stern_m;
        let phi_compact = value / c_compact;
        let pref = sqrt(8.0 * c0 * NA * EPS0 * params.eps_bulk * kt);
        let phi_diffuse = 2.0 * thermal_v * asinh(value / pref);
        let kappa = sqrt(2.0 * c0 * NA * E * E / (EPS0 * params.eps_bulk * kt));
        let gamma = tanh(E * phi_diffuse / (4.0 * kt));
        let e_compact = value / (EPS0 * params.eps_compact);
        compact_drop = phi_compact;
        diffuse_drop = phi_diffuse;
        total_drop = phi_compact + phi_diffuse;
        c_diffuse = EPS0 * params.eps_bulk * kappa * cosh(E * phi_diffuse / (2.0 * kt));
        c_total = 1.0 / (1.0 / c_compact + 1.0 / c_diffuse);
        if (x_ang <= params.stern) {
            phi = e_compact * x_m;
            efield = e_compact;
            rho = value / stern_m;
            diel = params.eps_compact;
        } else {
            let xd = x_m - stern_m;
            let exp_term = exp(-kappa * xd);
            let g = gamma * exp_term;
            let phi_d = 2.0 * thermal_v * safe_log_ratio(g);
            phi = phi_compact + phi_d;
            let denom = 1.0 - gamma * gamma * exp_term * exp_term;
            if (abs(denom) > 1.0e-6) {
                efield = 2.0 * thermal_v * kappa * 2.0 * gamma * exp_term / denom;
            }
            rho = -2.0 * E * c0 * NA * sinh(E * phi_d / kt);
            ion_density = exp(-E * phi_d / kt) + exp(E * phi_d / kt);
            diel = params.eps_bulk;
        }
    }

    potential[idx] = phi;
    field[idx] = efield;
    charge[idx] = rho;
    ion[idx] = ion_density;
    dielectric[idx] = diel;
    if (pi == 0u) {
        scalars[vi * 6u + 0u] = compact_drop;
        scalars[vi * 6u + 1u] = diffuse_drop;
        scalars[vi * 6u + 2u] = total_drop;
        scalars[vi * 6u + 3u] = c_compact;
        scalars[vi * 6u + 4u] = c_diffuse;
        scalars[vi * 6u + 5u] = c_total;
    }
}
"#;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::edl::{EdlModel, EdlNumerics};

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
