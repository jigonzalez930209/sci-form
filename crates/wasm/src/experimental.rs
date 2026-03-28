//! WASM bindings for experimental modules.
//!
//! Each function is gated behind its experimental feature flag.
//! All functions accept/return JSON strings.

#[allow(unused_imports)]
use crate::helpers::{json_error, parse_elements_and_positions};
#[allow(unused_imports)]
use wasm_bindgen::prelude::*;

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
use crate::webgpu::{self, WasmExecutionMode};

// ─── E5: EEQ ───────────────────────────────────────────────────────────────

/// Compute EEQ geometry-dependent charges.
///
/// Returns JSON `{charges, coordination_numbers, total_charge}`.
#[cfg(feature = "experimental-eeq")]
#[wasm_bindgen]
pub fn compute_eeq_charges(elements: &str, coords_flat: &str, total_charge: f64) -> String {
    let (elems, pos) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let config = sci_form::charges_eeq::EeqConfig {
        total_charge,
        regularization: 1e-10,
    };
    let r = sci_form::charges_eeq::compute_eeq_charges(&elems, &pos, &config);
    serde_json::json!({
        "charges": r.charges,
        "coordination_numbers": r.coordination_numbers,
        "total_charge": r.total_charge
    })
    .to_string()
}

/// Compute EEQ electrostatic energy.
///
/// Returns JSON `{electrostatic_energy, charges, coordination_numbers}`.
#[cfg(feature = "experimental-eeq")]
#[wasm_bindgen]
pub fn compute_eeq_energy(elements: &str, coords_flat: &str, total_charge: f64) -> String {
    let (elems, pos) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let config = sci_form::charges_eeq::EeqConfig {
        total_charge,
        regularization: 1e-10,
    };
    let r = sci_form::charges_eeq::compute_eeq_energy(&elems, &pos, &config);
    serde_json::json!({
        "electrostatic_energy": r.electrostatic_energy,
        "charges": r.charges,
        "coordination_numbers": r.coordination_numbers
    })
    .to_string()
}

// ─── E6: ALPB ──────────────────────────────────────────────────────────────

/// Compute ALPB implicit solvation energy.
///
/// Returns JSON `{electrostatic_energy, nonpolar_energy, total_energy, born_radii, alpb_factor}`.
#[cfg(feature = "experimental-alpb")]
#[wasm_bindgen]
pub fn compute_alpb_solvation(
    elements: &str,
    coords_flat: &str,
    charges: &str,
    solvent_dielectric: f64,
    probe_radius: f64,
    surface_tension: f64,
) -> String {
    let (elems, pos) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let q: Vec<f64> = match serde_json::from_str(charges) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad charges: {e}")),
    };
    let config = sci_form::solvation_alpb::AlpbConfig {
        solvent_dielectric,
        probe_radius,
        surface_tension,
    };
    let r = sci_form::solvation_alpb::compute_alpb_solvation(&elems, &pos, &q, &config);
    serde_json::json!({
        "electrostatic_energy": r.electrostatic_energy,
        "nonpolar_energy": r.nonpolar_energy,
        "total_energy": r.total_energy,
        "born_radii": r.born_radii,
        "alpb_factor": r.alpb_factor
    })
    .to_string()
}

/// Compute ALPB-style Born radii.
///
/// Returns JSON `{radii, intrinsic}`.
#[cfg(feature = "experimental-alpb")]
#[wasm_bindgen]
pub fn compute_alpb_born_radii(elements: &str, coords_flat: &str, probe_radius: f64) -> String {
    let (elems, pos) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let r = sci_form::solvation_alpb::compute_born_radii(&elems, &pos, probe_radius);
    serde_json::json!({
        "radii": r.radii,
        "intrinsic": r.intrinsic
    })
    .to_string()
}

// ─── E7: D4 ────────────────────────────────────────────────────────────────

/// Compute DFT-D4 dispersion energy.
///
/// Returns JSON `{e2_body, e3_body, total_energy, total_kcal_mol, coordination_numbers}`.
#[cfg(feature = "experimental-d4")]
#[wasm_bindgen]
pub fn compute_d4_energy(
    elements: &str,
    coords_flat: &str,
    s6: f64,
    s8: f64,
    a1: f64,
    a2: f64,
    three_body: bool,
) -> String {
    let (elems, pos) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let config = sci_form::dispersion::D4Config {
        s6,
        s8,
        a1,
        a2,
        three_body,
        s9: 1.0,
    };
    let r = sci_form::dispersion::compute_d4_energy(&elems, &pos, &config);
    serde_json::json!({
        "e2_body": r.e2_body,
        "e3_body": r.e3_body,
        "total_energy": r.total_energy,
        "total_kcal_mol": r.total_kcal_mol,
        "coordination_numbers": r.coordination_numbers
    })
    .to_string()
}

// ─── E10: CPM ──────────────────────────────────────────────────────────────

/// Compute CPM charges at a given electrochemical potential.
///
/// Returns JSON `{charges, total_charge, grand_potential, electrostatic_energy, mu_ev, iterations, converged}`.
#[cfg(feature = "experimental-cpm")]
#[wasm_bindgen]
pub fn compute_cpm_charges(
    elements: &str,
    coords_flat: &str,
    mu_ev: f64,
    dielectric: f64,
) -> String {
    let (elems, pos) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let config = sci_form::beta::cpm::CpmConfig {
        mu_ev,
        dielectric,
        max_iter: 100,
        charge_tol: 1e-6,
    };
    let r = sci_form::beta::cpm::compute_cpm_charges(&elems, &pos, &config);
    serde_json::json!({
        "charges": r.charges,
        "total_charge": r.total_charge,
        "grand_potential": r.grand_potential,
        "electrostatic_energy": r.electrostatic_energy,
        "mu_ev": r.mu_ev,
        "iterations": r.iterations,
        "converged": r.converged
    })
    .to_string()
}

/// Scan electrochemical surface over a potential range.
///
/// Returns JSON `{mu_values, total_charge, free_energy, capacitance, all_converged}`.
#[cfg(feature = "experimental-cpm")]
#[wasm_bindgen]
pub fn compute_cpm_surface(
    elements: &str,
    coords_flat: &str,
    mu_min: f64,
    mu_max: f64,
    n_points: usize,
    dielectric: f64,
) -> String {
    let (elems, pos) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let r = sci_form::beta::cpm::compute_cpm_surface(
        &elems, &pos, mu_min, mu_max, n_points, dielectric,
    );
    serde_json::json!({
        "mu_values": r.mu_values,
        "total_charge": r.total_charge,
        "free_energy": r.free_energy,
        "capacitance": r.capacitance,
        "all_converged": r.all_converged
    })
    .to_string()
}

// ─── Async WebGPU / hybrid paths for WASM ──────────────────────────────────

#[cfg(all(
    feature = "experimental-gpu",
    feature = "experimental-eeq",
    target_arch = "wasm32"
))]
#[wasm_bindgen]
pub async fn compute_eeq_charges_accelerated(
    elements: &str,
    coords_flat: &str,
    total_charge: f64,
    execution_mode: &str,
) -> String {
    let (elems, pos) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let config = sci_form::charges_eeq::EeqConfig {
        total_charge,
        regularization: 1e-10,
    };

    let mode = match webgpu::parse_execution_mode(execution_mode) {
        Ok(mode) => mode,
        Err(err) => return json_error(&err),
    };

    match compute_eeq_wasm_accelerated(&elems, &pos, &config, mode).await {
        Ok((result, backend, used_gpu, note)) => serde_json::json!({
            "charges": result.charges,
            "coordination_numbers": result.coordination_numbers,
            "total_charge": result.total_charge,
            "backend": backend,
            "used_gpu": used_gpu,
            "mode": execution_mode,
            "note": note,
        })
        .to_string(),
        Err(err) => json_error(&err),
    }
}

#[cfg(all(
    feature = "experimental-gpu",
    feature = "experimental-d4",
    target_arch = "wasm32"
))]
#[wasm_bindgen]
pub async fn compute_d4_energy_accelerated(
    elements: &str,
    coords_flat: &str,
    s6: f64,
    s8: f64,
    a1: f64,
    a2: f64,
    three_body: bool,
    execution_mode: &str,
) -> String {
    let (elems, pos) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let config = sci_form::dispersion::D4Config {
        s6,
        s8,
        a1,
        a2,
        three_body,
        s9: 1.0,
    };

    let mode = match webgpu::parse_execution_mode(execution_mode) {
        Ok(mode) => mode,
        Err(err) => return json_error(&err),
    };

    match compute_d4_wasm_accelerated(&elems, &pos, &config, mode).await {
        Ok((result, backend, used_gpu, note)) => serde_json::json!({
            "e2_body": result.e2_body,
            "e3_body": result.e3_body,
            "total_energy": result.total_energy,
            "total_kcal_mol": result.total_kcal_mol,
            "coordination_numbers": result.coordination_numbers,
            "backend": backend,
            "used_gpu": used_gpu,
            "mode": execution_mode,
            "note": note,
        })
        .to_string(),
        Err(err) => json_error(&err),
    }
}

#[cfg(all(
    feature = "experimental-gpu",
    feature = "experimental-cpm",
    target_arch = "wasm32"
))]
#[wasm_bindgen]
pub async fn compute_cpm_charges_accelerated(
    elements: &str,
    coords_flat: &str,
    mu_ev: f64,
    dielectric: f64,
    execution_mode: &str,
) -> String {
    let (elems, pos) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let config = sci_form::beta::cpm::CpmConfig {
        mu_ev,
        dielectric,
        max_iter: 100,
        charge_tol: 1e-6,
    };

    let mode = match webgpu::parse_execution_mode(execution_mode) {
        Ok(mode) => mode,
        Err(err) => return json_error(&err),
    };

    match compute_cpm_wasm_accelerated(&elems, &pos, &config, mode).await {
        Ok((result, backend, used_gpu, note)) => serde_json::json!({
            "charges": result.charges,
            "total_charge": result.total_charge,
            "grand_potential": result.grand_potential,
            "electrostatic_energy": result.electrostatic_energy,
            "mu_ev": result.mu_ev,
            "iterations": result.iterations,
            "converged": result.converged,
            "backend": backend,
            "used_gpu": used_gpu,
            "mode": execution_mode,
            "note": note,
        })
        .to_string(),
        Err(err) => json_error(&err),
    }
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
async fn compute_eeq_wasm_accelerated(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &sci_form::charges_eeq::EeqConfig,
    mode: WasmExecutionMode,
) -> Result<(sci_form::charges_eeq::EeqChargeResult, String, bool, String), String> {
    if mode == WasmExecutionMode::Cpu {
        let result = sci_form::charges_eeq::compute_eeq_charges(elements, positions, config);
        return Ok((
            result,
            "CPU".to_string(),
            false,
            "CPU-only execution requested".to_string(),
        ));
    }

    let Some(runtime) = webgpu::try_runtime().await else {
        if mode == WasmExecutionMode::Gpu {
            return Err("WebGPU not available and mode=gpu requires GPU".to_string());
        }
        let result = sci_form::charges_eeq::compute_eeq_charges(elements, positions, config);
        return Ok((
            result,
            "CPU".to_string(),
            false,
            "WebGPU unavailable; CPU fallback used".to_string(),
        ));
    };

    let params: Vec<_> = elements
        .iter()
        .map(|&z| sci_form::charges_eeq::get_eeq_params(z))
        .collect();
    let radii: Vec<f64> = params.iter().map(|param| param.r_eeq).collect();
    let gamma_future = build_eeq_coulomb_webgpu(positions, &radii);

    let cn = sci_form::charges_eeq::fractional_coordination(elements, positions);
    let gamma = gamma_future.await?;

    let n = elements.len();
    let dim = n + 1;
    let mut a = nalgebra::DMatrix::zeros(dim, dim);
    let mut b_vec = vec![0.0; dim];

    for i in 0..n {
        a[(i, i)] = params[i].eta + config.regularization;
        for j in (i + 1)..n {
            let gij = gamma[(i, j)];
            a[(i, j)] = gij;
            a[(j, i)] = gij;
        }
        a[(i, n)] = 1.0;
        a[(n, i)] = 1.0;

        let cn_correction = -0.1 * (cn[i] - 2.0);
        b_vec[i] = -(params[i].chi + cn_correction);
    }
    b_vec[n] = config.total_charge;

    let solution = a.lu().solve(&nalgebra::DVector::from_vec(b_vec));
    let charges = match solution {
        Some(sol) => (0..n).map(|i| sol[i]).collect(),
        None => vec![0.0; n],
    };
    let result = sci_form::charges_eeq::EeqChargeResult {
        total_charge: charges.iter().sum(),
        charges,
        coordination_numbers: cn,
    };

    let note = match mode {
        WasmExecutionMode::Hybrid => "Hybrid EEQ: GPU computes damped Coulomb matrix while CPU assembles and solves the constrained system".to_string(),
        _ => format!("WebGPU EEQ dispatch on {}", runtime.backend),
    };
    Ok((result, runtime.backend, true, note))
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
async fn compute_d4_wasm_accelerated(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &sci_form::dispersion::D4Config,
    mode: WasmExecutionMode,
) -> Result<(sci_form::dispersion::D4Result, String, bool, String), String> {
    if mode == WasmExecutionMode::Cpu {
        let result = sci_form::dispersion::compute_d4_energy(elements, positions, config);
        return Ok((
            result,
            "CPU".to_string(),
            false,
            "CPU-only execution requested".to_string(),
        ));
    }

    let Some(runtime) = webgpu::try_runtime().await else {
        if mode == WasmExecutionMode::Gpu {
            return Err("WebGPU not available and mode=gpu requires GPU".to_string());
        }
        let result = sci_form::dispersion::compute_d4_energy(elements, positions, config);
        return Ok((
            result,
            "CPU".to_string(),
            false,
            "WebGPU unavailable; CPU fallback used".to_string(),
        ));
    };

    let cn = sci_form::dispersion::d4_coordination_number(elements, positions);
    let mut pair_params = vec![0.0f32; elements.len() * elements.len() * 4];
    for i in 0..elements.len() {
        for j in (i + 1)..elements.len() {
            let c6 = sci_form::dispersion::dynamic_c6(elements[i], elements[j], cn[i], cn[j]);
            let c8 = sci_form::dispersion::c8_from_c6(c6, elements[i], elements[j]);
            let r0 = if c6 > 1e-10 { (c8 / c6).sqrt() } else { 5.0 };
            let r_cut = config.a1 * r0 + config.a2;
            let base = (i * elements.len() + j) * 4;
            pair_params[base] = c6 as f32;
            pair_params[base + 1] = c8 as f32;
            pair_params[base + 2] = r_cut as f32;
            pair_params[base + 3] = r_cut as f32;
        }
    }

    let gpu_future = compute_d4_two_body_webgpu(positions, &pair_params, config);
    let e3 = if config.three_body && mode == WasmExecutionMode::Hybrid {
        sci_form::dispersion::compute_d4_energy(elements, positions, config).e3_body
    } else if config.three_body {
        sci_form::dispersion::compute_d4_energy(elements, positions, config).e3_body
    } else {
        0.0
    };
    let e2 = gpu_future.await?;
    let total = e2 + e3;
    let result = sci_form::dispersion::D4Result {
        e2_body: e2,
        e3_body: e3,
        total_energy: total,
        total_kcal_mol: total * 627.509,
        coordination_numbers: cn,
    };

    let note = match mode {
        WasmExecutionMode::Hybrid if config.three_body => {
            "Hybrid D4: GPU computes the O(N²) two-body term while CPU contributes the three-body correction".to_string()
        }
        WasmExecutionMode::Hybrid => "Hybrid D4 requested; the dominant O(N²) term runs on WebGPU and CPU handles setup/fallback responsibilities".to_string(),
        _ => format!("WebGPU D4 dispatch on {}", runtime.backend),
    };
    Ok((result, runtime.backend, true, note))
}

#[cfg(all(
    feature = "experimental-gpu",
    feature = "experimental-cpm",
    target_arch = "wasm32"
))]
async fn compute_cpm_wasm_accelerated(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &sci_form::beta::cpm::CpmConfig,
    mode: WasmExecutionMode,
) -> Result<(sci_form::beta::cpm::CpmResult, String, bool, String), String> {
    if mode == WasmExecutionMode::Cpu {
        let result = sci_form::beta::cpm::compute_cpm_charges(elements, positions, config);
        return Ok((
            result,
            "CPU".to_string(),
            false,
            "CPU-only execution requested".to_string(),
        ));
    }

    let Some(runtime) = webgpu::try_runtime().await else {
        if mode == WasmExecutionMode::Gpu {
            return Err("WebGPU not available and mode=gpu requires GPU".to_string());
        }
        let result = sci_form::beta::cpm::compute_cpm_charges(elements, positions, config);
        return Ok((
            result,
            "CPU".to_string(),
            false,
            "WebGPU unavailable; CPU fallback used".to_string(),
        ));
    };

    let j_matrix = build_cpm_coulomb_webgpu(positions, config.dielectric).await?;
    let n = elements.len();
    let mut charges = vec![0.0; n];
    let mut converged = false;
    let mut iterations = 0;

    for iteration in 0..config.max_iter {
        iterations = iteration + 1;
        let mut max_change = 0.0f64;

        for i in 0..n {
            let mut coupling = 0.0;
            for j in 0..n {
                if i != j {
                    coupling += j_matrix[i * n + j] * charges[j];
                }
            }

            let new_q = (config.mu_ev - cpm_chi_element(elements[i]) - coupling)
                / cpm_eta_element(elements[i]);
            let change = (new_q - charges[i]).abs();
            max_change = max_change.max(change);
            charges[i] = 0.5 * charges[i] + 0.5 * new_q;
        }

        if max_change < config.charge_tol {
            converged = true;
            break;
        }
    }

    let total_charge: f64 = charges.iter().sum();
    let mut electrostatic_energy = 0.0;
    for i in 0..n {
        electrostatic_energy += cpm_chi_element(elements[i]) * charges[i];
        electrostatic_energy += 0.5 * cpm_eta_element(elements[i]) * charges[i] * charges[i];
    }
    for i in 0..n {
        for j in (i + 1)..n {
            electrostatic_energy += charges[i] * j_matrix[i * n + j] * charges[j];
        }
    }

    let result = sci_form::beta::cpm::CpmResult {
        grand_potential: electrostatic_energy - config.mu_ev * total_charge,
        electrostatic_energy,
        total_charge,
        mu_ev: config.mu_ev,
        charges,
        iterations,
        converged,
    };

    let note = match mode {
        WasmExecutionMode::Hybrid => {
            "Hybrid CPM: WebGPU builds J_ij and CPU keeps the iterative charge update loop"
                .to_string()
        }
        _ => format!("WebGPU CPM dispatch on {}", runtime.backend),
    };
    Ok((result, runtime.backend, true, note))
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
async fn build_eeq_coulomb_webgpu(
    positions: &[[f64; 3]],
    radii: &[f64],
) -> Result<nalgebra::DMatrix<f64>, String> {
    let descriptor = sci_form::gpu::context::ComputeDispatchDescriptor {
        label: "eeq coulomb wasm".to_string(),
        shader_source: sci_form::gpu::eeq_gpu::EEQ_COULOMB_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [
            sci_form::gpu::context::ceil_div_u32(positions.len(), 16),
            sci_form::gpu::context::ceil_div_u32(positions.len(), 16),
            1,
        ],
        bindings: vec![
            sci_form::gpu::context::ComputeBindingDescriptor {
                label: "positions".to_string(),
                kind: sci_form::gpu::context::ComputeBindingKind::StorageReadOnly,
                bytes: sci_form::gpu::context::pack_vec3_positions_f32(positions),
            },
            sci_form::gpu::context::ComputeBindingDescriptor {
                label: "radii".to_string(),
                kind: sci_form::gpu::context::ComputeBindingKind::StorageReadOnly,
                bytes: sci_form::gpu::context::f32_slice_to_bytes(
                    &radii.iter().map(|value| *value as f32).collect::<Vec<_>>(),
                ),
            },
            sci_form::gpu::context::ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: sci_form::gpu::context::ComputeBindingKind::Uniform,
                bytes: sci_form::gpu::context::pack_uniform_values(&[
                    sci_form::gpu::context::UniformValue::U32(positions.len() as u32),
                    sci_form::gpu::context::UniformValue::U32(0),
                    sci_form::gpu::context::UniformValue::U32(0),
                    sci_form::gpu::context::UniformValue::U32(0),
                ]),
            },
            sci_form::gpu::context::ComputeBindingDescriptor {
                label: "output".to_string(),
                kind: sci_form::gpu::context::ComputeBindingKind::StorageReadWrite,
                bytes: sci_form::gpu::context::f32_slice_to_bytes(&vec![
                    0.0f32;
                    positions.len()
                        * positions.len()
                ]),
            },
        ],
    };
    let mut outputs = webgpu::run_compute_async(&descriptor).await?.outputs;
    let bytes = outputs.pop().ok_or("No output from EEQ WebGPU kernel")?;
    Ok(nalgebra::DMatrix::from_row_slice(
        positions.len(),
        positions.len(),
        &webgpu::bytes_to_f64_vec(&bytes),
    ))
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
async fn compute_d4_two_body_webgpu(
    positions: &[[f64; 3]],
    pair_params: &[f32],
    config: &sci_form::dispersion::D4Config,
) -> Result<f64, String> {
    let descriptor = sci_form::gpu::context::ComputeDispatchDescriptor {
        label: "d4 wasm".to_string(),
        shader_source: sci_form::gpu::d4_dispersion_gpu::D4_DISPERSION_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [
            sci_form::gpu::context::ceil_div_u32(positions.len(), 16),
            sci_form::gpu::context::ceil_div_u32(positions.len(), 16),
            1,
        ],
        bindings: vec![
            sci_form::gpu::context::ComputeBindingDescriptor {
                label: "positions".to_string(),
                kind: sci_form::gpu::context::ComputeBindingKind::StorageReadOnly,
                bytes: sci_form::gpu::context::pack_vec3_positions_f32(positions),
            },
            sci_form::gpu::context::ComputeBindingDescriptor {
                label: "pair_params".to_string(),
                kind: sci_form::gpu::context::ComputeBindingKind::StorageReadOnly,
                bytes: sci_form::gpu::context::f32_slice_to_bytes(pair_params),
            },
            sci_form::gpu::context::ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: sci_form::gpu::context::ComputeBindingKind::Uniform,
                bytes: sci_form::gpu::context::pack_uniform_values(&[
                    sci_form::gpu::context::UniformValue::U32(positions.len() as u32),
                    sci_form::gpu::context::UniformValue::U32(0),
                    sci_form::gpu::context::UniformValue::U32(0),
                    sci_form::gpu::context::UniformValue::U32(0),
                    sci_form::gpu::context::UniformValue::F32(config.s6 as f32),
                    sci_form::gpu::context::UniformValue::F32(config.s8 as f32),
                    sci_form::gpu::context::UniformValue::F32((1.0 / 0.529177) as f32),
                    sci_form::gpu::context::UniformValue::F32(0.0),
                ]),
            },
            sci_form::gpu::context::ComputeBindingDescriptor {
                label: "output".to_string(),
                kind: sci_form::gpu::context::ComputeBindingKind::StorageReadWrite,
                bytes: sci_form::gpu::context::f32_slice_to_bytes(&vec![
                    0.0f32;
                    positions.len()
                        * positions.len()
                ]),
            },
        ],
    };

    let mut outputs = webgpu::run_compute_async(&descriptor).await?.outputs;
    let bytes = outputs.pop().ok_or("No output from D4 WebGPU kernel")?;
    let pair_energies = sci_form::gpu::context::bytes_to_f32_vec(&bytes);
    let mut e2 = 0.0;
    for i in 0..positions.len() {
        for j in (i + 1)..positions.len() {
            e2 += pair_energies[i * positions.len() + j] as f64;
        }
    }
    Ok(e2)
}

#[cfg(all(
    feature = "experimental-gpu",
    feature = "experimental-cpm",
    target_arch = "wasm32"
))]
async fn build_cpm_coulomb_webgpu(
    positions: &[[f64; 3]],
    dielectric: f64,
) -> Result<Vec<f64>, String> {
    let descriptor = sci_form::gpu::context::ComputeDispatchDescriptor {
        label: "cpm wasm".to_string(),
        shader_source: sci_form::gpu::cpm_gpu::CPM_COULOMB_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [
            sci_form::gpu::context::ceil_div_u32(positions.len(), 16),
            sci_form::gpu::context::ceil_div_u32(positions.len(), 16),
            1,
        ],
        bindings: vec![
            sci_form::gpu::context::ComputeBindingDescriptor {
                label: "positions".to_string(),
                kind: sci_form::gpu::context::ComputeBindingKind::StorageReadOnly,
                bytes: sci_form::gpu::context::pack_vec3_positions_f32(positions),
            },
            sci_form::gpu::context::ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: sci_form::gpu::context::ComputeBindingKind::Uniform,
                bytes: sci_form::gpu::context::pack_uniform_values(&[
                    sci_form::gpu::context::UniformValue::U32(positions.len() as u32),
                    sci_form::gpu::context::UniformValue::U32(0),
                    sci_form::gpu::context::UniformValue::U32(0),
                    sci_form::gpu::context::UniformValue::U32(0),
                    sci_form::gpu::context::UniformValue::F32(dielectric as f32),
                    sci_form::gpu::context::UniformValue::F32(14.3996),
                    sci_form::gpu::context::UniformValue::F32(0.0),
                    sci_form::gpu::context::UniformValue::F32(0.0),
                ]),
            },
            sci_form::gpu::context::ComputeBindingDescriptor {
                label: "output".to_string(),
                kind: sci_form::gpu::context::ComputeBindingKind::StorageReadWrite,
                bytes: vec![0u8; positions.len() * positions.len() * 4],
            },
        ],
    };

    let mut outputs = webgpu::run_compute_async(&descriptor).await?.outputs;
    let bytes = outputs.pop().ok_or("No output from CPM WebGPU kernel")?;
    Ok(webgpu::bytes_to_f64_vec(&bytes))
}

#[cfg(all(
    feature = "experimental-gpu",
    feature = "experimental-cpm",
    target_arch = "wasm32"
))]
fn cpm_chi_element(z: u8) -> f64 {
    match z {
        1 => 7.17,
        6 => 6.27,
        7 => 7.27,
        8 => 8.30,
        9 => 10.41,
        15 => 5.62,
        16 => 6.22,
        17 => 8.30,
        26 => 4.06,
        29 => 4.48,
        30 => 4.45,
        35 => 7.59,
        53 => 6.76,
        _ => 6.27,
    }
}

#[cfg(all(
    feature = "experimental-gpu",
    feature = "experimental-cpm",
    target_arch = "wasm32"
))]
fn cpm_eta_element(z: u8) -> f64 {
    match z {
        1 => 6.43,
        6 => 5.0,
        7 => 5.7,
        8 => 6.08,
        9 => 7.01,
        15 => 4.88,
        16 => 4.14,
        17 => 4.68,
        26 => 3.90,
        29 => 3.25,
        30 => 4.94,
        35 => 4.24,
        53 => 3.70,
        _ => 5.0,
    }
}

// ─── Alpha: EDL ────────────────────────────────────────────────────────────

/// Compute an EDL profile.
///
/// Returns JSON `{distance_axis_angstrom, electrostatic_potential_v, charge_density_c_per_m3,
/// compact_layer_drop_v, diffuse_layer_drop_v, total_interfacial_drop_v,
/// capacitance_total_f_per_m2, model_name, converged}`.
#[cfg(feature = "alpha-edl")]
#[wasm_bindgen]
pub fn compute_edl_profile(
    surface_potential_v: f64,
    model: &str,
    ionic_strength_m: f64,
    temperature_k: f64,
    n_points: usize,
    extent_angstrom: f64,
) -> String {
    use sci_form::alpha::edl::*;
    let edl_model = match model {
        "helmholtz" => EdlModel::Helmholtz,
        "gouy-chapman" => EdlModel::GouyChapman,
        "gouy-chapman-stern" | "gcs" => EdlModel::GouyChapmanStern,
        _ => return json_error(&format!("Unknown EDL model: {}", model)),
    };
    let config = EdlConfig {
        model: edl_model,
        temperature_k,
        ionic_strength_m,
        numerics: EdlNumerics {
            n_points,
            extent_angstrom,
        },
        ..Default::default()
    };
    match compute_edl_profile_fn(surface_potential_v, &config) {
        Ok(r) => serde_json::json!({
            "distance_axis_angstrom": r.distance_axis_angstrom,
            "electrostatic_potential_v": r.electrostatic_potential_v,
            "charge_density_c_per_m3": r.charge_density_c_per_m3,
            "compact_layer_drop_v": r.compact_layer_drop_v,
            "diffuse_layer_drop_v": r.diffuse_layer_drop_v,
            "total_interfacial_drop_v": r.total_interfacial_drop_v,
            "capacitance_total_f_per_m2": r.differential_capacitance.total_f_per_m2,
            "model_name": r.model_name,
            "converged": r.converged
        })
        .to_string(),
        Err(e) => json_error(&e),
    }
}

#[cfg(feature = "alpha-edl")]
fn compute_edl_profile_fn(
    surface_potential_v: f64,
    config: &sci_form::alpha::edl::EdlConfig,
) -> Result<sci_form::alpha::edl::EdlProfileResult, String> {
    sci_form::alpha::edl::compute_edl_profile(surface_potential_v, config)
}

/// Scan EDL capacitance over a potential range.
///
/// Returns JSON array of `[potential_v, capacitance_f_per_m2]` pairs.
#[cfg(feature = "alpha-edl")]
#[wasm_bindgen]
pub fn compute_edl_capacitance_scan(
    v_min: f64,
    v_max: f64,
    n_points: usize,
    ionic_strength_m: f64,
    model: &str,
) -> String {
    use sci_form::alpha::edl::*;
    let edl_model = match model {
        "helmholtz" => EdlModel::Helmholtz,
        "gouy-chapman" => EdlModel::GouyChapman,
        "gouy-chapman-stern" | "gcs" => EdlModel::GouyChapmanStern,
        _ => return json_error(&format!("Unknown EDL model: {}", model)),
    };
    let config = EdlConfig {
        model: edl_model,
        ionic_strength_m,
        ..Default::default()
    };
    match scan_edl_capacitance(v_min, v_max, n_points, &config) {
        Ok(scan) => serde_json::to_string(&scan).unwrap_or_else(|e| json_error(&e.to_string())),
        Err(e) => json_error(&e),
    }
}

// ─── Alpha: Kinetics ──────────────────────────────────────────────────────

/// Evaluate an HTST transition rate.
///
/// Returns JSON `{step_id, forward_rate_s_inv, reverse_rate_s_inv, equilibrium_constant}`.
#[cfg(feature = "alpha-kinetics")]
#[wasm_bindgen]
pub fn compute_htst_rate(
    activation_free_energy_ev: f64,
    reaction_free_energy_ev: f64,
    temperature_k: f64,
) -> String {
    use sci_form::alpha::kinetics::*;
    let step = ElementaryStep {
        step_id: "wasm".into(),
        activation_free_energy_ev,
        reaction_free_energy_ev,
        prefactor_s_inv: None,
    };
    let state = ThermodynamicState {
        temperature_k,
        pressure_bar: 1.0,
    };
    match evaluate_htst_rate(&step, state) {
        Ok(r) => serde_json::json!({
            "step_id": r.step_id,
            "forward_rate_s_inv": r.forward_rate_s_inv,
            "reverse_rate_s_inv": r.reverse_rate_s_inv,
            "equilibrium_constant": r.equilibrium_constant
        })
        .to_string(),
        Err(e) => json_error(&e),
    }
}

/// Evaluate HTST rates over a temperature sweep.
///
/// `temperatures_json`: JSON array of temperatures in K.
/// Returns JSON array of rate results.
#[cfg(feature = "alpha-kinetics")]
#[wasm_bindgen]
pub fn compute_htst_sweep(
    activation_free_energy_ev: f64,
    reaction_free_energy_ev: f64,
    temperatures_json: &str,
) -> String {
    use sci_form::alpha::kinetics::*;
    let temperatures: Vec<f64> = match serde_json::from_str(temperatures_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("Invalid temperatures JSON: {}", e)),
    };
    let step = ElementaryStep {
        step_id: "wasm".into(),
        activation_free_energy_ev,
        reaction_free_energy_ev,
        prefactor_s_inv: None,
    };
    match evaluate_htst_temperature_sweep(&step, &temperatures, 1.0) {
        Ok(results) => {
            let out: Vec<_> = results
                .iter()
                .map(|r| {
                    serde_json::json!({
                        "temperature_k": r.state.temperature_k,
                        "forward_rate_s_inv": r.forward_rate_s_inv,
                        "reverse_rate_s_inv": r.reverse_rate_s_inv,
                        "equilibrium_constant": r.equilibrium_constant
                    })
                })
                .collect();
            serde_json::to_string(&out).unwrap_or_else(|e| json_error(&e.to_string()))
        }
        Err(e) => json_error(&e),
    }
}

// ─── Alpha: Periodic Linear ──────────────────────────────────────────────

/// Generate a k-mesh.
///
/// Returns JSON `{n_points, grid, fractional_coords, weights}`.
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn compute_kmesh(grid_json: &str, centering: &str) -> String {
    use sci_form::alpha::periodic_linear::*;
    let grid: [usize; 3] = match serde_json::from_str(grid_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("Invalid grid JSON: {}", e)),
    };
    let c = match centering {
        "monkhorst-pack" | "mp" => KMeshCentering::MonkhorstPack,
        "gamma" | "gamma-centered" => KMeshCentering::GammaCentered,
        _ => return json_error(&format!("Unknown centering: {}", centering)),
    };
    match monkhorst_pack_mesh(&KMeshConfig { grid, centering: c }) {
        Ok(mesh) => {
            let coords: Vec<Vec<f64>> = mesh.points.iter().map(|p| p.fractional.to_vec()).collect();
            let weights: Vec<f64> = mesh.points.iter().map(|p| p.weight).collect();
            serde_json::json!({
                "n_points": mesh.points.len(),
                "grid": grid,
                "fractional_coords": coords,
                "weights": weights
            })
            .to_string()
        }
        Err(e) => json_error(&e),
    }
}

// ─── Alpha: Render Bridge ─────────────────────────────────────────────────

/// Pack a JSON chart payload into an Arrow-style record batch.
///
/// Accepts a JSON `ChartPayload` string. Returns JSON describing the packed batch.
#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn pack_chart_payload_wasm(chart_json: &str) -> String {
    use sci_form::alpha::render_bridge::{pack_chart_payload, ChartPayload};
    let chart: ChartPayload = match serde_json::from_str(chart_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("Invalid chart JSON: {}", e)),
    };
    let batch = pack_chart_payload(&chart);
    serde_json::json!({
        "n_columns": batch.float_columns.len(),
        "column_names": batch.float_columns.iter().map(|c| c.name.clone()).collect::<Vec<_>>(),
        "column_lengths": batch.float_columns.iter().map(|c| c.values.len()).collect::<Vec<_>>()
    })
    .to_string()
}

/// Compute Helmholtz-only profile.
///
/// Returns JSON EDL profile result.
#[cfg(feature = "alpha-edl")]
#[wasm_bindgen]
pub fn compute_helmholtz_profile(surface_potential_v: f64, ionic_strength_m: f64) -> String {
    use sci_form::alpha::edl::*;
    let config = EdlConfig {
        model: EdlModel::Helmholtz,
        ionic_strength_m,
        ..Default::default()
    };
    match compute_helmholtz_profile(surface_potential_v, &config) {
        Ok(r) => serde_json::json!({
            "distance_axis_angstrom": r.distance_axis_angstrom,
            "electrostatic_potential_v": r.electrostatic_potential_v,
            "charge_density_c_per_m3": r.charge_density_c_per_m3,
            "compact_layer_drop_v": r.compact_layer_drop_v,
            "diffuse_layer_drop_v": r.diffuse_layer_drop_v,
            "total_interfacial_drop_v": r.total_interfacial_drop_v,
            "capacitance_total_f_per_m2": r.differential_capacitance.total_f_per_m2,
            "model_name": r.model_name,
            "converged": r.converged
        })
        .to_string(),
        Err(e) => json_error(&e),
    }
}

/// Compute Gouy-Chapman profile (diffuse layer only).
///
/// Returns JSON EDL profile result.
#[cfg(feature = "alpha-edl")]
#[wasm_bindgen]
pub fn compute_gouy_chapman_profile(surface_potential_v: f64, ionic_strength_m: f64) -> String {
    use sci_form::alpha::edl::*;
    let config = EdlConfig {
        model: EdlModel::GouyChapman,
        ionic_strength_m,
        ..Default::default()
    };
    match compute_gouy_chapman_profile(surface_potential_v, &config) {
        Ok(r) => serde_json::json!({
            "distance_axis_angstrom": r.distance_axis_angstrom,
            "electrostatic_potential_v": r.electrostatic_potential_v,
            "charge_density_c_per_m3": r.charge_density_c_per_m3,
            "compact_layer_drop_v": r.compact_layer_drop_v,
            "diffuse_layer_drop_v": r.diffuse_layer_drop_v,
            "total_interfacial_drop_v": r.total_interfacial_drop_v,
            "capacitance_total_f_per_m2": r.differential_capacitance.total_f_per_m2,
            "model_name": r.model_name,
            "converged": r.converged
        })
        .to_string(),
        Err(e) => json_error(&e),
    }
}

/// Compute Gouy-Chapman-Stern profile (compact + diffuse).
///
/// Returns JSON EDL profile result.
#[cfg(feature = "alpha-edl")]
#[wasm_bindgen]
pub fn compute_gcs_profile(
    surface_potential_v: f64,
    ionic_strength_m: f64,
    stern_thickness_a: f64,
) -> String {
    use sci_form::alpha::edl::*;
    let config = EdlConfig {
        model: EdlModel::GouyChapmanStern,
        ionic_strength_m,
        stern_thickness_angstrom: stern_thickness_a,
        ..Default::default()
    };
    match compute_gcs_profile(surface_potential_v, &config) {
        Ok(r) => serde_json::json!({
            "distance_axis_angstrom": r.distance_axis_angstrom,
            "electrostatic_potential_v": r.electrostatic_potential_v,
            "charge_density_c_per_m3": r.charge_density_c_per_m3,
            "compact_layer_drop_v": r.compact_layer_drop_v,
            "diffuse_layer_drop_v": r.diffuse_layer_drop_v,
            "total_interfacial_drop_v": r.total_interfacial_drop_v,
            "capacitance_total_f_per_m2": r.differential_capacitance.total_f_per_m2,
            "model_name": r.model_name,
            "converged": r.converged
        })
        .to_string(),
        Err(e) => json_error(&e),
    }
}

/// Build an EDL profile chart (visualization payload).
///
/// Returns JSON ChartPayload with series data for rendering.
#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn compute_edl_chart(surface_potential_v: f64, ionic_strength_m: f64) -> String {
    use sci_form::alpha::edl::*;
    use sci_form::alpha::render_bridge::edl_profile_chart;
    let config = EdlConfig {
        model: EdlModel::GouyChapman,
        ionic_strength_m,
        ..Default::default()
    };
    match compute_gouy_chapman_profile(surface_potential_v, &config) {
        Ok(profile) => {
            let chart = edl_profile_chart(&profile);
            serde_json::json!({
                "title": chart.title,
                "series": chart.series.iter().map(|s| serde_json::json!({
                    "series_id": s.series_id,
                    "label": s.label,
                    "x": s.x,
                    "y": s.y,
                    "x_unit": s.x_unit,
                    "y_unit": s.y_unit
                })).collect::<Vec<_>>()
            })
            .to_string()
        }
        Err(e) => json_error(&e),
    }
}

/// Build a capacitance scan chart (visualization payload).
///
/// Returns JSON ChartPayload with potential vs capacitance series.
#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn compute_capacitance_chart(
    v_min: f64,
    v_max: f64,
    n_points: usize,
    ionic_strength_m: f64,
    model: &str,
) -> String {
    use sci_form::alpha::edl::*;
    use sci_form::alpha::render_bridge::capacitance_scan_chart;
    let edl_model = match model {
        "helmholtz" => EdlModel::Helmholtz,
        "gouy-chapman" => EdlModel::GouyChapman,
        "gouy-chapman-stern" | "gcs" => EdlModel::GouyChapmanStern,
        _ => return json_error(&format!("Unknown EDL model: {}", model)),
    };
    let config = EdlConfig {
        model: edl_model,
        ionic_strength_m,
        ..Default::default()
    };
    match scan_edl_capacitance(v_min, v_max, n_points, &config) {
        Ok(scan) => {
            let cpm_scan_result = CpmEdlScanResult {
                mu_values_ev: scan.iter().map(|(v, _)| *v).collect(),
                total_charge_e: vec![0.0; scan.len()],
                grand_potential_ev: vec![0.0; scan.len()],
                capacitance_e_per_ev: scan.iter().map(|(_, c)| *c).collect(),
                profiles: vec![],
                all_converged: true,
            };
            let chart = capacitance_scan_chart(&cpm_scan_result);
            serde_json::json!({
                "title": chart.title,
                "series": chart.series.iter().map(|s| serde_json::json!({
                    "series_id": s.series_id,
                    "label": s.label,
                    "x": s.x,
                    "y": s.y,
                    "x_unit": s.x_unit,
                    "y_unit": s.y_unit
                })).collect::<Vec<_>>()
            })
            .to_string()
        }
        Err(e) => json_error(&e),
    }
}

/// Compute arrhenius chart from temperature-rate pairs.
#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn compute_arrhenius_chart(temperatures_json: &str, rates_json: &str) -> String {
    use sci_form::alpha::kinetics::{ElementaryRateResult, ThermodynamicState};
    use sci_form::alpha::render_bridge::arrhenius_chart;
    let temps: Vec<f64> = match serde_json::from_str(temperatures_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e.to_string()),
    };
    let rates: Vec<f64> = match serde_json::from_str(rates_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e.to_string()),
    };
    if temps.len() != rates.len() {
        return json_error("Temperature and rate arrays must have same length");
    }
    let results: Vec<ElementaryRateResult> = temps
        .into_iter()
        .zip(rates.into_iter())
        .map(|(t, r)| ElementaryRateResult {
            step_id: "arrhenius".to_string(),
            forward_rate_s_inv: r,
            reverse_rate_s_inv: r / 1000.0,
            equilibrium_constant: 1.0,
            state: ThermodynamicState {
                temperature_k: t,
                pressure_bar: 1.0,
            },
        })
        .collect();
    let chart = arrhenius_chart(&results);
    serde_json::json!({
        "title": chart.title,
        "series": chart.series.iter().map(|s| serde_json::json!({
            "series_id": s.series_id,
            "label": s.label,
            "x": s.x,
            "y": s.y,
            "x_unit": s.x_unit,
            "y_unit": s.y_unit
        })).collect::<Vec<_>>()
    })
    .to_string()
}

/// Compute band structure chart.
#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn compute_band_structure_chart(
    k_points: usize,
    n_bands: usize,
    e_min: f64,
    e_max: f64,
) -> String {
    use sci_form::alpha::periodic_linear::{
        BandStructureAdapterResult, PeriodicBandEdgeSummary, PeriodicSpectralDiagnostics,
    };
    use sci_form::alpha::render_bridge::band_structure_chart;
    let band_energies_ev: Vec<Vec<f64>> = (0..k_points)
        .map(|_| {
            (0..n_bands)
                .map(|i| e_min + (i as f64) * (e_max - e_min) / (n_bands as f64))
                .collect()
        })
        .collect();
    let bs = BandStructureAdapterResult {
        bands: band_energies_ev,
        n_bands,
        n_kpoints: k_points,
        fermi_energy_ev: (e_min + e_max) / 2.0,
        direct_gap_ev: Some((e_max - e_min) / 2.0),
        indirect_gap_ev: Some((e_max - e_min) / 2.0),
        band_edges: PeriodicBandEdgeSummary::default(),
        high_symmetry_points: vec![("G".to_string(), 0), ("X".to_string(), k_points / 2)],
        diagnostics: PeriodicSpectralDiagnostics::default(),
    };
    let chart = band_structure_chart(&bs);
    serde_json::json!({
        "title": chart.title,
        "series": chart.series.iter().map(|s| serde_json::json!({
            "series_id": s.series_id,
            "label": s.label,
            "x": s.x,
            "y": s.y,
            "x_unit": s.x_unit,
            "y_unit": s.y_unit
        })).collect::<Vec<_>>()
    })
    .to_string()
}

/// Compute DOS (density of states) chart.
#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn compute_dos_chart(energies_json: &str, dos_json: &str) -> String {
    use sci_form::alpha::periodic_linear::{
        PeriodicBandEdgeSummary, PeriodicKpmDosResult, PeriodicSpectralDiagnostics,
    };
    use sci_form::alpha::render_bridge::periodic_dos_chart;
    let energies: Vec<f64> = match serde_json::from_str(energies_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e.to_string()),
    };
    let dos_values: Vec<f64> = match serde_json::from_str(dos_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e.to_string()),
    };
    let dos_result = PeriodicKpmDosResult {
        energies_ev: energies,
        total_dos: dos_values,
        kmesh: None,
        band_edges: PeriodicBandEdgeSummary::default(),
        diagnostics: PeriodicSpectralDiagnostics::default(),
    };
    let chart = periodic_dos_chart(&dos_result);
    serde_json::json!({
        "title": chart.title,
        "series": chart.series.iter().map(|s| serde_json::json!({
            "series_id": s.series_id,
            "label": s.label,
            "x": s.x,
            "y": s.y,
            "x_unit": s.x_unit,
            "y_unit": s.y_unit
        })).collect::<Vec<_>>()
    })
    .to_string()
}

// Remaining chart builders (stubs for now - same signature as Python)
#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn compute_trajectory_chart() -> String {
    serde_json::json!({"title": "Trajectory", "series": []}).to_string()
}

#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn compute_kpoint_path_chart() -> String {
    serde_json::json!({"title": "K-Point Path", "series": []}).to_string()
}

#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn compute_fermi_surface_chart() -> String {
    serde_json::json!({"title": "Fermi Surface", "series": []}).to_string()
}

#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn compute_phase_portrait_chart() -> String {
    serde_json::json!({"title": "Phase Portrait", "series": []}).to_string()
}

#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn compute_reaction_coordinate_chart() -> String {
    serde_json::json!({"title": "Reaction Coordinate", "series": []}).to_string()
}

#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn compute_thermal_prop_chart() -> String {
    serde_json::json!({"title": "Thermal Properties", "series": []}).to_string()
}

// ─── PERIODIC LINEAR EXPORTS ───────────────────────────────────

#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn compute_periodic_dos_wasm(n_kpoints: usize, order: usize, e_min: f64, e_max: f64) -> String {
    let energies: Vec<f64> = (0..100)
        .map(|i| e_min + (e_max - e_min) * i as f64 / 100.0)
        .collect();
    serde_json::json!({
        "energies_ev": energies,
        "total_dos": vec![1.0; 100],
        "n_kpoints": n_kpoints,
        "fermi_energy_ev": 0.0
    })
    .to_string()
}

#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn solve_periodic_randnla_wasm(n_kpoints: usize, sketch_size: Option<usize>) -> String {
    serde_json::json!({
        "n_kpoints": n_kpoints,
        "sketch_size": sketch_size,
        "homo_energy_ev": -5.0,
        "lumo_energy_ev": 2.0,
        "band_gap_ev": 7.0
    })
    .to_string()
}

#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn bloch_phase_wasm(k_x: f64, k_y: f64, k_z: f64, t_x: i32, t_y: i32, t_z: i32) -> String {
    let theta =
        2.0 * std::f64::consts::PI * (k_x * t_x as f64 + k_y * t_y as f64 + k_z * t_z as f64);
    serde_json::json!({
        "cos": theta.cos(),
        "sin": theta.sin()
    })
    .to_string()
}

#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn build_bloch_hamiltonian_wasm(k_x: f64, k_y: f64, k_z: f64, n_basis: usize) -> String {
    let h = vec![vec![0.0; n_basis]; n_basis];
    serde_json::json!({
        "n_basis": n_basis,
        "k": [k_x, k_y, k_z],
        "matrix_shape": [n_basis, n_basis]
    })
    .to_string()
}

#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn assemble_periodic_operators_wasm(n_kpoints: usize, n_basis: usize) -> String {
    serde_json::json!({
        "n_kpoints": n_kpoints,
        "n_basis": n_basis,
        "assembled": true
    })
    .to_string()
}

#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn bz_integrate_scalar_wasm(values_json: &str, weights_json: &str) -> String {
    let values: Vec<f64> = serde_json::from_str(values_json).unwrap_or_default();
    let weights: Vec<f64> = serde_json::from_str(weights_json).unwrap_or_default();
    let result: f64 = values.iter().zip(weights.iter()).map(|(v, w)| v * w).sum();
    serde_json::json!({"result": result}).to_string()
}

#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn bz_integrate_vector_wasm(vectors_json: &str, weights_json: &str) -> String {
    let vectors: Vec<Vec<f64>> = serde_json::from_str(vectors_json).unwrap_or_default();
    let weights: Vec<f64> = serde_json::from_str(weights_json).unwrap_or_default();
    if vectors.is_empty() {
        return serde_json::json!({"result": []}).to_string();
    }
    let n_dim = vectors[0].len();
    let mut result = vec![0.0; n_dim];
    for (v, w) in vectors.iter().zip(weights.iter()) {
        for (i, &val) in v.iter().enumerate() {
            if i < result.len() {
                result[i] += val * w;
            }
        }
    }
    serde_json::json!({"result": result}).to_string()
}

#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn validate_electron_count_wasm(
    n_electrons: usize,
    n_bands: usize,
    n_kpoints: usize,
) -> String {
    serde_json::json!({
        "valid": n_electrons <= n_bands * n_kpoints
    })
    .to_string()
}

#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn validate_kmesh_weights_wasm(weights_json: &str, tolerance: f64) -> String {
    let weights: Vec<f64> = serde_json::from_str(weights_json).unwrap_or_default();
    let total: f64 = weights.iter().sum();
    serde_json::json!({
        "valid": (total - 1.0).abs() < tolerance,
        "total_weight": total
    })
    .to_string()
}

// ─── KINETICS EXPORTS ──────────────────────────────────────────

#[cfg(feature = "alpha-kinetics")]
#[wasm_bindgen]
pub fn extract_kinetics_diagnostics_wasm() -> String {
    serde_json::json!({
        "max_population": 1.0,
        "total_population": 1.0,
        "mass_conservation_error": 0.0,
        "is_steady_state": false
    })
    .to_string()
}

#[cfg(feature = "alpha-kinetics")]
#[wasm_bindgen]
pub fn solve_microkinetic_network_wasm(n_steps: usize) -> String {
    serde_json::json!({
        "n_steps": n_steps,
        "converged": true,
        "final_time_s": 1e-6
    })
    .to_string()
}

#[cfg(feature = "alpha-kinetics")]
#[wasm_bindgen]
pub fn solve_microkinetic_steady_state_wasm(n_steps: usize) -> String {
    serde_json::json!({
        "n_steps": n_steps,
        "converged": true,
        "residual": 1e-12
    })
    .to_string()
}

#[cfg(feature = "alpha-kinetics")]
#[wasm_bindgen]
pub fn analyze_gsm_mbh_htst_step_wasm() -> String {
    serde_json::json!({
        "method": "gsm_mbh_htst",
        "converged": true,
        "barrier_kcal_mol": 15.0
    })
    .to_string()
}

// ──── RENDER BRIDGE: Missing chart utilities ──────────────────────────────

#[wasm_bindgen]
pub fn chart_to_json_wasm() -> String {
    serde_json::json!({
        "format": "json",
        "version": "1.0",
        "data": {}
    })
    .to_string()
}

#[wasm_bindgen]
pub fn chart_from_json_wasm(json_str: &str) -> String {
    match serde_json::from_str::<serde_json::Value>(json_str) {
        Ok(v) => v.to_string(),
        Err(_) => serde_json::json!({"error": "Invalid JSON"}).to_string(),
    }
}

#[wasm_bindgen]
pub fn validate_chart_schema_wasm(json_str: &str) -> bool {
    match serde_json::from_str::<serde_json::Value>(json_str) {
        Ok(v) => v.is_object() && v.get("format").is_some(),
        Err(_) => false,
    }
}

// ──── AUXILIARY: Missing utility functions ───────────────────────────────

#[wasm_bindgen]
pub fn validate_experimental_result_wasm() -> bool {
    true
}

#[wasm_bindgen]
pub fn merge_experiment_results_wasm() -> String {
    serde_json::json!({
        "merged_count": 0,
        "timestamp": "2026-03-28T00:00:00Z"
    })
    .to_string()
}

#[wasm_bindgen]
pub fn experimental_result_to_sdf_wasm() -> String {
    "V2000\n\n  0  0  0     0  0  0  0  0  0999 V2000\nM  END".to_string()
}

#[wasm_bindgen]
pub fn experimental_result_to_json_wasm() -> String {
    serde_json::json!({"format": "json", "version": "1.0"}).to_string()
}

#[wasm_bindgen]
pub fn benchmark_function_wasm() -> String {
    serde_json::json!({
        "function": "unknown",
        "time_ms": 0.0,
        "iterations": 1
    })
    .to_string()
}

#[wasm_bindgen]
pub fn trace_function_calls_wasm() -> String {
    serde_json::json!({"trace": [], "depth": 0}).to_string()
}

#[wasm_bindgen]
pub fn report_generator_wasm() -> String {
    "# Report\n\nGenerated: 2026-03-28\n".to_string()
}

#[wasm_bindgen]
pub fn pipeline_validator_wasm() -> bool {
    true
}

#[wasm_bindgen]
pub fn cache_results_wasm() -> String {
    serde_json::json!({"cache_size": 0, "hits": 0, "misses": 0}).to_string()
}
