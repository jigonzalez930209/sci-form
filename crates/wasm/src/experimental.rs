//! WASM bindings for experimental modules.
//!
//! Each function is gated behind its experimental feature flag.
//! All functions accept/return JSON strings.

#[allow(unused_imports)]
use crate::helpers::{bad_json, invalid_value, json_error, parse_elements_and_positions};
#[allow(unused_imports)]
use wasm_bindgen::prelude::*;

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
use crate::webgpu::{self, WasmExecutionMode};

// ─── E5: EEQ ───────────────────────────────────────────────────────────────

/// Compute EEQ geometry-dependent charges.
///
/// Returns JSON `{charges, coordination_numbers, total_charge}`.
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
        Err(e) => return json_error(&bad_json("'charges'", charges, &e)),
    };
    let config = sci_form::solvation_alpb::AlpbConfig {
        solvent_dielectric,
        probe_radius,
        surface_tension,
    };
    let r = match sci_form::solvation_alpb::compute_alpb_solvation(&elems, &pos, &q, &config) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
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
    let r = match sci_form::solvation_alpb::compute_born_radii(&elems, &pos, probe_radius) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
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
        mu_per_atom: None,
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

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
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
        _ => {
            return json_error(&invalid_value(
                "'model'",
                model,
                &["helmholtz", "gouy-chapman", "gouy-chapman-stern", "gcs"],
            ))
        }
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
        _ => {
            return json_error(&invalid_value(
                "'model'",
                model,
                &["helmholtz", "gouy-chapman", "gouy-chapman-stern", "gcs"],
            ))
        }
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
        Err(e) => return json_error(&bad_json("'temperatures'", temperatures_json, &e)),
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
        Err(e) => return json_error(&bad_json("'grid'", grid_json, &e)),
    };
    let c = match centering {
        "monkhorst-pack" | "mp" => KMeshCentering::MonkhorstPack,
        "gamma" | "gamma-centered" => KMeshCentering::GammaCentered,
        _ => {
            return json_error(&invalid_value(
                "'centering'",
                centering,
                &["monkhorst-pack", "mp", "gamma", "gamma-centered"],
            ))
        }
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

// ─── Shared helpers ───────────────────────────────────────────────────────

/// Owned real-space operator blocks: matrix plus lattice translation.
#[cfg(feature = "alpha-periodic-linear")]
type RealSpaceBlocks = Vec<(nalgebra::DMatrix<f64>, [i32; 3])>;

/// Borrowed real-space operator blocks in the shape expected by the core API.
#[cfg(feature = "alpha-periodic-linear")]
type RealSpaceBlockRefs<'a> = Vec<(&'a nalgebra::DMatrix<f64>, [i32; 3])>;

/// Parse a JSON array of real-space operator blocks.
///
/// A block is `{"matrix": [[...], ...], "translation": [i, j, k]}` where the
/// matrix is square and the translation is expressed in lattice-vector units.
#[cfg(feature = "alpha-periodic-linear")]
fn parse_blocks(json_str: &str) -> Result<RealSpaceBlocks, String> {
    let value: serde_json::Value = serde_json::from_str(json_str)
        .map_err(|e| format!("invalid real-space block JSON: {e}"))?;
    let array = value
        .as_array()
        .ok_or_else(|| "real-space blocks must be a JSON array".to_string())?;
    let mut blocks = Vec::with_capacity(array.len());
    for item in array {
        let matrix = item
            .get("matrix")
            .and_then(|m| m.as_array())
            .ok_or_else(|| "each real-space block needs a 'matrix' array".to_string())?;
        let translation = item
            .get("translation")
            .and_then(|t| t.as_array())
            .ok_or_else(|| "each real-space block needs a 'translation' array".to_string())?;
        if translation.len() != 3 {
            return Err("lattice translations need exactly three components".into());
        }
        let mut trans = [0i32; 3];
        for (slot, value) in trans.iter_mut().zip(translation.iter()) {
            *slot = value
                .as_i64()
                .ok_or_else(|| "lattice translations must be integers".to_string())?
                as i32;
        }
        let n = matrix.len();
        if n == 0 {
            return Err("real-space matrices must be non-empty".into());
        }
        let mut flat = Vec::with_capacity(n * n);
        for row in matrix {
            let row = row
                .as_array()
                .ok_or_else(|| "real-space matrix rows must be arrays".to_string())?;
            if row.len() != n {
                return Err("real-space matrices must be square".into());
            }
            for entry in row {
                flat.push(
                    entry
                        .as_f64()
                        .ok_or_else(|| "real-space matrix entries must be numbers".to_string())?,
                );
            }
        }
        blocks.push((nalgebra::DMatrix::from_row_slice(n, n, &flat), trans));
    }
    if blocks.is_empty() {
        return Err("at least one real-space block is required".into());
    }
    Ok(blocks)
}

/// Borrow parsed blocks in the shape expected by the core k-mesh API.
#[cfg(feature = "alpha-periodic-linear")]
fn as_core_refs(blocks: &RealSpaceBlocks) -> RealSpaceBlockRefs<'_> {
    blocks
        .iter()
        .map(|(matrix, trans)| (matrix, *trans))
        .collect()
}

/// Parse a k-mesh centering label.
#[cfg(feature = "alpha-periodic-linear")]
fn parse_centering(centering: &str) -> Result<sci_form::kmesh::KMeshCentering, String> {
    use sci_form::kmesh::KMeshCentering;
    match centering {
        "monkhorst-pack" | "mp" => Ok(KMeshCentering::MonkhorstPack),
        "gamma" | "gamma-centered" => Ok(KMeshCentering::GammaCentered),
        _ => Err(invalid_value(
            "'centering'",
            centering,
            &["monkhorst-pack", "mp", "gamma", "gamma-centered"],
        )),
    }
}

/// Row-major representation of a square matrix.
#[cfg(feature = "alpha-periodic-linear")]
fn matrix_rows(matrix: &nalgebra::DMatrix<f64>) -> Vec<Vec<f64>> {
    let n = matrix.nrows();
    (0..n)
        .map(|i| (0..n).map(|j| matrix[(i, j)]).collect())
        .collect()
}

/// Ascending eigenvalues of a symmetric matrix.
#[cfg(feature = "alpha-periodic-linear")]
fn symmetric_eigenvalues(matrix: &nalgebra::DMatrix<f64>) -> Vec<f64> {
    let eigen = nalgebra::SymmetricEigen::new(matrix.clone());
    let mut values: Vec<f64> = eigen.eigenvalues.iter().copied().collect();
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    values
}

/// Convert a flat coordinate array into Cartesian positions.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
fn positions_from_flat(coords: &[f64]) -> Vec<[f64; 3]> {
    coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect()
}

/// Parse the nine lattice-vector components into a 3×3 lattice.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
fn parse_lattice(values: &[f64]) -> Result<[[f64; 3]; 3], String> {
    if values.len() != 9 {
        return Err("lattice needs exactly nine components".into());
    }
    let mut lattice = [[0.0_f64; 3]; 3];
    for (row, chunk) in values.chunks_exact(3).enumerate() {
        lattice[row] = [chunk[0], chunk[1], chunk[2]];
    }
    Ok(lattice)
}

/// Interpolate a high-symmetry k-path and return the cumulative path distance
/// together with the sampled fractional k-points.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
fn interpolate_path(
    path: &[([f64; 3], String)],
    n_points_per_segment: usize,
) -> (Vec<f64>, Vec<[f64; 3]>) {
    let n = n_points_per_segment.max(1);
    let mut distance = Vec::new();
    let mut points = Vec::new();
    let mut accumulated = 0.0_f64;
    let mut previous: Option<[f64; 3]> = None;
    for index in 0..path.len().saturating_sub(1) {
        let start = path[index].0;
        let end = path[index + 1].0;
        for step in 0..=n {
            if index > 0 && step == 0 {
                continue;
            }
            let t = step as f64 / n as f64;
            let mut k = [0.0_f64; 3];
            for axis in 0..3 {
                k[axis] = start[axis] + t * (end[axis] - start[axis]);
            }
            if let Some(prev) = previous {
                accumulated += ((k[0] - prev[0]).powi(2)
                    + (k[1] - prev[1]).powi(2)
                    + (k[2] - prev[2]).powi(2))
                .sqrt();
            }
            previous = Some(k);
            distance.push(accumulated);
            points.push(k);
        }
    }
    (distance, points)
}

/// Compute the Fermi-surface contour of the selected band in the `kz = 0` plane.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
fn fermi_contour(
    refs: &RealSpaceBlockRefs<'_>,
    fermi_energy_ev: Option<f64>,
    n_kx: usize,
    n_ky: usize,
) -> (Vec<f64>, Vec<f64>, f64) {
    let n_kx = n_kx.max(2);
    let n_ky = n_ky.max(2);
    let kx_values: Vec<f64> = (0..n_kx).map(|i| -0.5 + i as f64 / n_kx as f64).collect();
    let ky_values: Vec<f64> = (0..n_ky).map(|j| -0.5 + j as f64 / n_ky as f64).collect();

    let mut samples: Vec<Vec<Vec<f64>>> = Vec::with_capacity(n_ky);
    let mut all_energies: Vec<f64> = Vec::new();
    for &ky in &ky_values {
        let mut row = Vec::with_capacity(n_kx);
        for &kx in &kx_values {
            let hamiltonian = sci_form::kmesh::build_bloch_hamiltonian(refs, [kx, ky, 0.0]);
            let eigenvalues = symmetric_eigenvalues(&hamiltonian);
            all_energies.extend(eigenvalues.iter().copied());
            row.push(eigenvalues);
        }
        samples.push(row);
    }

    let fermi = match fermi_energy_ev {
        Some(value) => value,
        None => {
            let mut sorted = all_energies.clone();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            if sorted.is_empty() {
                0.0
            } else {
                sorted[sorted.len() / 2]
            }
        }
    };

    let band = samples
        .first()
        .and_then(|row| row.first())
        .map(|energies| {
            energies
                .iter()
                .filter(|&&e| e < fermi)
                .count()
                .min(energies.len().saturating_sub(1))
        })
        .unwrap_or(0);

    let mut contour_x = Vec::new();
    let mut contour_y = Vec::new();
    for (j, &ky) in ky_values.iter().enumerate() {
        for i in 0..n_kx.saturating_sub(1) {
            let e0 = samples[j][i].get(band).copied().unwrap_or(f64::NAN);
            let e1 = samples[j][i + 1].get(band).copied().unwrap_or(f64::NAN);
            if !e0.is_finite() || !e1.is_finite() {
                continue;
            }
            let f0 = e0 - fermi;
            let f1 = e1 - fermi;
            if f0 == 0.0 {
                contour_x.push(kx_values[i]);
                contour_y.push(ky);
                continue;
            }
            if f0 * f1 > 0.0 {
                continue;
            }
            let span = kx_values[i + 1] - kx_values[i];
            let energy_at = |x: f64| -> f64 {
                if span == 0.0 {
                    e0
                } else {
                    e0 + (x - kx_values[i]) / span * (e1 - e0)
                }
            };
            let (mut lo, mut hi) = (kx_values[i], kx_values[i + 1]);
            for _ in 0..60 {
                let mid = 0.5 * (lo + hi);
                if (energy_at(mid) - fermi) * f0 <= 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            contour_x.push(0.5 * (lo + hi));
            contour_y.push(ky);
        }
    }
    (contour_x, contour_y, fermi)
}

/// Parse a `MicrokineticNetworkConfig`, falling back to a built-in A ⇌ B demo.
#[cfg(feature = "alpha-kinetics")]
fn network_config(
    network_json: &str,
    n_steps: usize,
) -> Result<sci_form::alpha::kinetics::MicrokineticNetworkConfig, String> {
    use sci_form::alpha::kinetics::{
        ElementaryStep, MicrokineticNetworkConfig, ThermodynamicState,
    };
    let trimmed = network_json.trim();
    if !trimmed.is_empty() && trimmed != "null" {
        let mut config: MicrokineticNetworkConfig = serde_json::from_str(trimmed)
            .map_err(|e| format!("invalid microkinetic network JSON: {e}"))?;
        if config.n_frames < 2 {
            config.n_frames = n_steps.max(2);
        }
        return Ok(config);
    }
    Ok(MicrokineticNetworkConfig {
        species: vec!["A".into(), "B".into()],
        steps: vec![ElementaryStep {
            step_id: "A->B".into(),
            activation_free_energy_ev: 0.5,
            reaction_free_energy_ev: -0.1,
            prefactor_s_inv: None,
        }],
        stoichiometry: vec![vec![(0, -1.0), (1, 1.0)]],
        initial_populations: vec![1.0, 0.0],
        total_time_s: 1e-3,
        n_frames: n_steps.max(2),
        state: ThermodynamicState {
            temperature_k: 298.15,
            pressure_bar: 1.0,
        },
    })
}

/// Parse either the chart wrapper envelope or a bare chart payload.
#[cfg(feature = "alpha-render-bridge")]
fn parse_chart(json_str: &str) -> Result<sci_form::render_bridge::ChartPayload, String> {
    let value: serde_json::Value =
        serde_json::from_str(json_str).map_err(|e| format!("invalid JSON: {e}"))?;
    let chart_value = value.get("data").cloned().unwrap_or(value);
    let chart: sci_form::render_bridge::ChartPayload =
        serde_json::from_value(chart_value).map_err(|e| format!("invalid chart payload: {e}"))?;
    if chart.series.iter().any(|s| s.x.len() != s.y.len()) {
        return Err("chart series must have matching x and y lengths".into());
    }
    Ok(chart)
}

/// Serialize a chart payload to the canonical chart JSON.
#[cfg(all(
    feature = "alpha-render-bridge",
    any(
        feature = "alpha-edl",
        feature = "alpha-kinetics",
        feature = "alpha-periodic-linear",
        feature = "alpha-gsm"
    )
))]
fn chart_json(chart: &sci_form::render_bridge::ChartPayload) -> String {
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

// ─── Alpha: Render Bridge ─────────────────────────────────────────────────

/// Pack a JSON chart payload into an Arrow-style record batch.
///
/// Accepts a JSON `ChartPayload` string (or the `chart-to-json` wrapper).
/// Returns JSON describing the packed batch.
#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn pack_chart_payload_wasm(chart_json: &str) -> String {
    use sci_form::render_bridge::{pack_chart_payload, verify_chart_round_trip};
    let chart = match parse_chart(chart_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let batch = pack_chart_payload(&chart);
    serde_json::json!({
        "n_columns": batch.float_columns.len(),
        "column_names": batch.float_columns.iter().map(|c| c.name.clone()).collect::<Vec<_>>(),
        "column_lengths": batch.float_columns.iter().map(|c| c.values.len()).collect::<Vec<_>>(),
        "round_trip_ok": verify_chart_round_trip(&chart, 0.0)
    })
    .to_string()
}

/// Alias of [`pack_chart_payload`], the name of the Rust root function.
///
/// The `_wasm` suffix existed only to avoid a name clash inside the flat
/// WASM namespace; the canonical name is exported as well (blocker C-H1).
#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn pack_chart_payload(chart_json: &str) -> String {
    pack_chart_payload_wasm(chart_json)
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
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-edl"))]
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
        Ok(profile) => chart_json(&edl_profile_chart(&profile)),
        Err(e) => json_error(&e),
    }
}

/// Build a capacitance scan chart (visualization payload).
///
/// The interfacial charge is obtained by integrating the computed charge
/// density of every EDL profile in the scan.
///
/// Returns JSON ChartPayload with potential vs capacitance series.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-edl"))]
#[wasm_bindgen]
pub fn compute_capacitance_chart(
    v_min: f64,
    v_max: f64,
    n_points: usize,
    ionic_strength_m: f64,
    model: &str,
) -> String {
    use sci_form::alpha::edl::*;
    use sci_form::render_bridge::{ChartPayload, ChartSeries};
    let edl_model = match model {
        "helmholtz" => EdlModel::Helmholtz,
        "gouy-chapman" => EdlModel::GouyChapman,
        "gouy-chapman-stern" | "gcs" => EdlModel::GouyChapmanStern,
        _ => {
            return json_error(&invalid_value(
                "'model'",
                model,
                &["helmholtz", "gouy-chapman", "gouy-chapman-stern", "gcs"],
            ))
        }
    };
    let config = EdlConfig {
        model: edl_model,
        ionic_strength_m,
        ..Default::default()
    };
    let scan = match scan_edl_capacitance(v_min, v_max, n_points, &config) {
        Ok(scan) => scan,
        Err(e) => return json_error(&e),
    };
    let mut potentials = Vec::with_capacity(scan.len());
    let mut capacitances = Vec::with_capacity(scan.len());
    let mut charges = Vec::with_capacity(scan.len());
    for (potential, capacitance) in &scan {
        let profile = match compute_edl_profile(*potential, &config) {
            Ok(profile) => profile,
            Err(e) => return json_error(&e),
        };
        potentials.push(*potential);
        capacitances.push(*capacitance);
        charges.push(integrated_charge_c_per_m2(&profile));
    }
    let chart = ChartPayload {
        title: "EDL Capacitance Scan".to_string(),
        series: vec![
            ChartSeries {
                series_id: "charge_vs_potential".into(),
                label: "Q(μ)".into(),
                x: potentials.clone(),
                y: charges,
                x_unit: "V".into(),
                y_unit: "C/m^2".into(),
            },
            ChartSeries {
                series_id: "capacitance".into(),
                label: "C(μ)".into(),
                x: potentials,
                y: capacitances,
                x_unit: "V".into(),
                y_unit: "F/m^2".into(),
            },
        ],
    };
    chart_json(&chart)
}

/// Integrate the charge density of an EDL profile into a surface charge
/// density in coulombs per square metre.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-edl"))]
fn integrated_charge_c_per_m2(profile: &sci_form::alpha::edl::EdlProfileResult) -> f64 {
    const ANGSTROM_TO_M: f64 = 1.0e-10;
    let mut integral = 0.0_f64;
    for i in 1..profile.distance_axis_angstrom.len() {
        let dx_m = (profile.distance_axis_angstrom[i] - profile.distance_axis_angstrom[i - 1])
            * ANGSTROM_TO_M;
        integral += 0.5
            * (profile.charge_density_c_per_m3[i] + profile.charge_density_c_per_m3[i - 1])
            * dx_m;
    }
    integral
}

/// Build an Arrhenius chart from temperature-rate pairs.
///
/// The chart is the core Arrhenius projection `ln(k)` vs `1000/T`; the
/// supplied forward rates are used verbatim.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-kinetics"))]
#[wasm_bindgen]
pub fn compute_arrhenius_chart(temperatures_json: &str, rates_json: &str) -> String {
    use sci_form::render_bridge::{ChartPayload, ChartSeries};
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
    let chart = ChartPayload {
        title: "Arrhenius Plot".into(),
        series: vec![ChartSeries {
            series_id: "arrhenius".into(),
            label: "ln(k) vs 1000/T".into(),
            x: temps.iter().map(|t| 1000.0 / t).collect(),
            y: rates.iter().map(|r| r.ln()).collect(),
            x_unit: "1000/K".into(),
            y_unit: "ln(s⁻¹)".into(),
        }],
    };
    chart_json(&chart)
}

/// Build a thermal-properties chart from an elementary step.
///
/// The HTST rates and equilibrium constants are computed by the core
/// transition-state evaluator over the requested temperature sweep.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-kinetics"))]
#[wasm_bindgen]
pub fn compute_thermal_prop_chart(
    activation_free_energy_ev: f64,
    reaction_free_energy_ev: f64,
    temperatures_json: &str,
) -> String {
    use sci_form::alpha::kinetics::{evaluate_htst_temperature_sweep, ElementaryStep};
    use sci_form::alpha::render_bridge::arrhenius_chart;
    use sci_form::render_bridge::{ChartPayload, ChartSeries};
    let temperatures: Vec<f64> = match serde_json::from_str(temperatures_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'temperatures'", temperatures_json, &e)),
    };
    let step = ElementaryStep {
        step_id: "thermal".to_string(),
        activation_free_energy_ev,
        reaction_free_energy_ev,
        prefactor_s_inv: None,
    };
    let results = match evaluate_htst_temperature_sweep(&step, &temperatures, 1.0) {
        Ok(results) => results,
        Err(e) => return json_error(&e),
    };
    let base = arrhenius_chart(&results);
    let mut series = base.series.into_iter().collect::<Vec<_>>();
    series.push(ChartSeries {
        series_id: "equilibrium_constant".into(),
        label: "ln(K_eq) vs 1000/T".into(),
        x: results
            .iter()
            .map(|r| 1000.0 / r.state.temperature_k)
            .collect(),
        y: results
            .iter()
            .map(|r| r.equilibrium_constant.ln())
            .collect(),
        x_unit: "1000/K".into(),
        y_unit: "ln".into(),
    });
    let chart = ChartPayload {
        title: "Thermal Properties".into(),
        series,
    };
    chart_json(&chart)
}

/// Build a phase portrait chart from a solved microkinetic network.
///
/// The two axes are the populations of the first two species over the
/// simulated trajectory.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-kinetics"))]
#[wasm_bindgen]
pub fn compute_phase_portrait_chart(network_json: &str, n_steps: usize) -> String {
    use sci_form::alpha::kinetics::solve_microkinetic_network;
    use sci_form::render_bridge::{ChartPayload, ChartSeries};
    let config = match network_config(network_json, n_steps) {
        Ok(config) => config,
        Err(e) => return json_error(&e),
    };
    if config.species.len() < 2 {
        return json_error("a phase portrait needs a network with at least two species");
    }
    let trace = match solve_microkinetic_network(&config) {
        Ok(trace) => trace,
        Err(e) => return json_error(&e),
    };
    let mut x = Vec::with_capacity(trace.frames.len());
    let mut y = Vec::with_capacity(trace.frames.len());
    for frame in &trace.frames {
        let first = frame
            .species
            .iter()
            .find(|s| s.species_id == config.species[0])
            .map(|s| s.population)
            .unwrap_or_default();
        let second = frame
            .species
            .iter()
            .find(|s| s.species_id == config.species[1])
            .map(|s| s.population)
            .unwrap_or_default();
        x.push(first);
        y.push(second);
    }
    let chart = ChartPayload {
        title: "Phase Portrait".into(),
        series: vec![ChartSeries {
            series_id: "phase_portrait".into(),
            label: format!("[{}] vs [{}]", config.species[0], config.species[1]),
            x,
            y,
            x_unit: "arb.".into(),
            y_unit: "arb.".into(),
        }],
    };
    chart_json(&chart)
}

/// Build a band-structure chart from a periodic EHT calculation.
///
/// `elements_json` are atomic numbers, `coords_json` the flat Cartesian
/// coordinates of the unit cell and `lattice_json` the nine lattice-vector
/// components.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
#[wasm_bindgen]
pub fn compute_band_structure_chart(
    elements_json: &str,
    coords_json: &str,
    lattice_json: &str,
    n_kpoints_per_segment: usize,
) -> String {
    use sci_form::alpha::periodic_linear::*;
    use sci_form::alpha::render_bridge::band_structure_chart;
    let elements: Vec<u8> = match serde_json::from_str(elements_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'elements'", elements_json, &e)),
    };
    let coords: Vec<f64> = match serde_json::from_str(coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'coords'", coords_json, &e)),
    };
    let lattice_values: Vec<f64> = match serde_json::from_str(lattice_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'lattice'", lattice_json, &e)),
    };
    let lattice = match parse_lattice(&lattice_values) {
        Ok(lattice) => lattice,
        Err(e) => return json_error(&e),
    };
    let positions = positions_from_flat(&coords);
    let n_electrons = match sci_form::eht::solve_eht(&elements, &positions, None) {
        Ok(result) => result.n_electrons,
        Err(e) => return json_error(&e),
    };
    let config = BandStructureAdapterConfig {
        n_kpoints_per_segment,
        path: BandStructureAdapterConfig::default().path,
    };
    match compute_band_structure_adapter(&elements, &positions, &lattice, n_electrons, &config) {
        Ok(bands) => chart_json(&band_structure_chart(&bands)),
        Err(e) => json_error(&e),
    }
}

/// Compute DOS (density of states) chart.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
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
    chart_json(&periodic_dos_chart(&dos_result))
}

/// Build a k-point path chart from the core high-symmetry path.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
#[wasm_bindgen]
pub fn compute_kpoint_path_chart(n_points_per_segment: usize) -> String {
    use sci_form::alpha::periodic_linear::BandStructureAdapterConfig;
    use sci_form::render_bridge::{ChartPayload, ChartSeries};
    let path = BandStructureAdapterConfig::default().path;
    let (distance, points) = interpolate_path(&path, n_points_per_segment);
    let chart = ChartPayload {
        title: "K-Point Path".into(),
        series: vec![
            ChartSeries {
                series_id: "k_x".into(),
                label: "k_x".into(),
                x: distance.clone(),
                y: points.iter().map(|k| k[0]).collect(),
                x_unit: "1/B".into(),
                y_unit: "1/B".into(),
            },
            ChartSeries {
                series_id: "k_y".into(),
                label: "k_y".into(),
                x: distance.clone(),
                y: points.iter().map(|k| k[1]).collect(),
                x_unit: "1/B".into(),
                y_unit: "1/B".into(),
            },
            ChartSeries {
                series_id: "k_z".into(),
                label: "k_z".into(),
                x: distance,
                y: points.iter().map(|k| k[2]).collect(),
                x_unit: "1/B".into(),
                y_unit: "1/B".into(),
            },
        ],
    };
    chart_json(&chart)
}

/// Build a Fermi-surface slice chart from real-space operator blocks.
///
/// The band energies are evaluated on a regular `n_kx` × `n_ky` mesh in the
/// `kz = 0` plane; the contour where the selected band crosses the Fermi
/// energy is extracted by bisection.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
#[wasm_bindgen]
pub fn compute_fermi_surface_chart(
    blocks_json: &str,
    fermi_energy_ev: Option<f64>,
    n_kx: usize,
    n_ky: usize,
) -> String {
    use sci_form::render_bridge::{ChartPayload, ChartSeries};
    let blocks = match parse_blocks(blocks_json) {
        Ok(blocks) => blocks,
        Err(e) => return json_error(&e),
    };
    let refs = as_core_refs(&blocks);
    let (contour_x, contour_y, fermi) = fermi_contour(&refs, fermi_energy_ev, n_kx, n_ky);
    let chart = ChartPayload {
        title: "Fermi Surface".into(),
        series: vec![ChartSeries {
            series_id: "fermi_contour".into(),
            label: format!("E_F = {fermi:.4} eV"),
            x: contour_x,
            y: contour_y,
            x_unit: "1/B".into(),
            y_unit: "1/B".into(),
        }],
    };
    chart_json(&chart)
}

/// Build a trajectory chart from a GSM reaction path.
///
/// The path is computed by the core growing-string method with the requested
/// energy backend; the chart shows the energy profile and the accumulated
/// path length of the string.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-gsm"))]
#[wasm_bindgen]
pub fn compute_trajectory_chart(
    smiles: &str,
    reactant_coords_json: &str,
    product_coords_json: &str,
    method: &str,
) -> String {
    use sci_form::render_bridge::{ChartPayload, ChartSeries};
    let reactant: Vec<f64> = match serde_json::from_str(reactant_coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'reactant coords'", reactant_coords_json, &e)),
    };
    let product: Vec<f64> = match serde_json::from_str(product_coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'product coords'", product_coords_json, &e)),
    };
    let result = match gsm_path(smiles, &reactant, &product, method) {
        Ok(result) => result,
        Err(e) => return json_error(&e),
    };
    let distances = path_distances(&result.path_coords);
    let index: Vec<f64> = (0..result.path_energies.len()).map(|i| i as f64).collect();
    let chart = ChartPayload {
        title: "GSM Trajectory".into(),
        series: vec![
            ChartSeries {
                series_id: "path_energy".into(),
                label: "E(path)".into(),
                x: index.clone(),
                y: result.path_energies.clone(),
                x_unit: "node".into(),
                y_unit: "kcal/mol".into(),
            },
            ChartSeries {
                series_id: "path_distance".into(),
                label: "path length".into(),
                x: index,
                y: distances,
                x_unit: "node".into(),
                y_unit: "A".into(),
            },
        ],
    };
    let mut payload: serde_json::Value =
        serde_json::from_str(&chart_json(&chart)).unwrap_or(serde_json::Value::Null);
    let trajectory =
        sci_form::render_bridge::gsm_path_trajectory(&result.path_coords, &result.path_energies);
    payload["trajectory"] = serde_json::json!(trajectory
        .iter()
        .map(|frame| serde_json::json!({
            "frame_time": frame.frame_time,
            "coords": frame.coords
        }))
        .collect::<Vec<_>>());
    payload.to_string()
}

/// Build a reaction-coordinate chart from a GSM reaction path.
///
/// The reaction coordinate is the accumulated path length in Å; the chart
/// shows the potential energy along that coordinate.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-gsm"))]
#[wasm_bindgen]
pub fn compute_reaction_coordinate_chart(
    smiles: &str,
    reactant_coords_json: &str,
    product_coords_json: &str,
    method: &str,
) -> String {
    use sci_form::render_bridge::{ChartPayload, ChartSeries};
    let reactant: Vec<f64> = match serde_json::from_str(reactant_coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'reactant coords'", reactant_coords_json, &e)),
    };
    let product: Vec<f64> = match serde_json::from_str(product_coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'product coords'", product_coords_json, &e)),
    };
    let result = match gsm_path(smiles, &reactant, &product, method) {
        Ok(result) => result,
        Err(e) => return json_error(&e),
    };
    let chart = ChartPayload {
        title: "Reaction Coordinate".into(),
        series: vec![ChartSeries {
            series_id: "reaction_coordinate".into(),
            label: "E(rc)".into(),
            x: path_distances(&result.path_coords),
            y: result.path_energies.clone(),
            x_unit: "A".into(),
            y_unit: "kcal/mol".into(),
        }],
    };
    chart_json(&chart)
}

/// Run the core GSM search with the requested energy backend.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-gsm"))]
fn gsm_path(
    smiles: &str,
    reactant: &[f64],
    product: &[f64],
    method: &str,
) -> Result<sci_form::alpha::gsm::GsmResult, String> {
    use sci_form::alpha::gsm::{find_transition_state_with_backend, GsmConfig, GsmEnergyBackend};
    let backend: GsmEnergyBackend = method
        .parse()
        .map_err(|e: String| format!("unknown GSM backend '{method}': {e}"))?;
    find_transition_state_with_backend(smiles, reactant, product, backend, &GsmConfig::default())
}

/// Accumulated path length of the string nodes.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-gsm"))]
fn path_distances(path_coords: &[Vec<f64>]) -> Vec<f64> {
    let mut distances = Vec::with_capacity(path_coords.len());
    let mut accumulated = 0.0_f64;
    let mut previous: Option<&Vec<f64>> = None;
    for coords in path_coords {
        if let Some(prev) = previous {
            accumulated += prev
                .iter()
                .zip(coords.iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();
        }
        previous = Some(coords);
        distances.push(accumulated);
    }
    distances
}

// ─── PERIODIC LINEAR EXPORTS ───────────────────────────────────

/// Compute the k-averaged density of states with the core KPM solver.
///
/// `blocks_json` is a JSON array of real-space operator blocks; the Bloch
/// Hamiltonians are assembled for every k-point of the requested mesh.
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn compute_periodic_dos_wasm(
    blocks_json: &str,
    grid_json: &str,
    centering: &str,
    order: usize,
    e_min: f64,
    e_max: f64,
    n_points: usize,
    fermi_energy_ev: Option<f64>,
) -> String {
    use sci_form::alpha::periodic_linear::*;
    let blocks = match parse_blocks(blocks_json) {
        Ok(blocks) => blocks,
        Err(e) => return json_error(&e),
    };
    let grid: [usize; 3] = match serde_json::from_str(grid_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'grid'", grid_json, &e)),
    };
    let centering = match parse_centering(centering) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let kmesh_config = KMeshConfig { grid, centering };
    let mesh = match monkhorst_pack_mesh(&kmesh_config) {
        Ok(mesh) => mesh,
        Err(e) => return json_error(&e),
    };
    let hamiltonian_refs = as_core_refs(&blocks);
    let (hamiltonians, _overlaps) =
        sci_form::kmesh::assemble_periodic_operators(&hamiltonian_refs, &[], &mesh);
    let config = PeriodicKpmAdapterConfig {
        kmesh: kmesh_config,
        order,
        e_min_ev: e_min,
        e_max_ev: e_max,
        n_points,
        fermi_energy_ev,
        ..Default::default()
    };
    match compute_periodic_kpm_dos_from_operators(&hamiltonians, &config) {
        Ok(result) => serde_json::json!({
            "energies_ev": result.energies_ev,
            "total_dos": result.total_dos,
            "n_kpoints": result.diagnostics.n_kpoints,
            "order": order,
            "e_min": e_min,
            "e_max": e_max,
            "fermi_energy_ev": fermi_energy_ev,
            "band_edges": result.band_edges,
            "diagnostics": result.diagnostics,
        })
        .to_string(),
        Err(e) => json_error(&e),
    }
}

/// Alias of [`compute_periodic_dos`], the name of the Rust root function.
///
/// The `_wasm` suffix existed only to avoid a name clash inside the flat
/// WASM namespace; the canonical name is exported as well (blocker C-H1).
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn compute_periodic_dos(
    blocks_json: &str,
    grid_json: &str,
    centering: &str,
    order: usize,
    e_min: f64,
    e_max: f64,
    n_points: usize,
    fermi_energy_ev: Option<f64>,
) -> String {
    compute_periodic_dos_wasm(
        blocks_json,
        grid_json,
        centering,
        order,
        e_min,
        e_max,
        n_points,
        fermi_energy_ev,
    )
}

/// Solve the periodic generalized eigenproblem with the core RandNLA solver.
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn solve_periodic_randnla_wasm(
    hamiltonian_blocks_json: &str,
    overlap_blocks_json: &str,
    n_electrons: usize,
    grid_json: &str,
    centering: &str,
    sketch_size: Option<usize>,
) -> String {
    use sci_form::alpha::periodic_linear::*;
    let hamiltonian_blocks = match parse_blocks(hamiltonian_blocks_json) {
        Ok(blocks) => blocks,
        Err(e) => return json_error(&e),
    };
    let overlap_blocks = match parse_blocks(overlap_blocks_json) {
        Ok(blocks) => blocks,
        Err(e) => return json_error(&e),
    };
    let grid: [usize; 3] = match serde_json::from_str(grid_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'grid'", grid_json, &e)),
    };
    let centering = match parse_centering(centering) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let kmesh_config = KMeshConfig { grid, centering };
    let mesh = match monkhorst_pack_mesh(&kmesh_config) {
        Ok(mesh) => mesh,
        Err(e) => return json_error(&e),
    };
    let hamiltonian_refs = as_core_refs(&hamiltonian_blocks);
    let overlap_refs = as_core_refs(&overlap_blocks);
    let (hamiltonians, overlaps) =
        sci_form::kmesh::assemble_periodic_operators(&hamiltonian_refs, &overlap_refs, &mesh);
    let config = PeriodicRandNlaAdapterConfig {
        kmesh: kmesh_config,
        sketch_size,
        ..Default::default()
    };
    match sci_form::alpha::periodic_linear::solve_periodic_randnla(
        &hamiltonians,
        &overlaps,
        n_electrons,
        &config,
    ) {
        Ok(result) => serde_json::json!({
            "n_kpoints": result.diagnostics.n_kpoints,
            "sketch_size": sketch_size,
            "homo_energy_ev": result.band_edges.homo_energy_ev,
            "lumo_energy_ev": result.band_edges.lumo_energy_ev,
            "band_gap_ev": result.band_edges.direct_gap_ev,
            "direct_gap_ev": result.band_edges.direct_gap_ev,
            "indirect_gap_ev": result.band_edges.indirect_gap_ev,
            "max_residual": result.diagnostics.max_residual,
            "used_fallback": result.diagnostics.used_fallback,
            "kpoint_results": result.kpoint_results,
        })
        .to_string(),
        Err(e) => json_error(&e),
    }
}

/// Alias of [`solve_periodic_randnla`], the name of the Rust root function.
///
/// The `_wasm` suffix existed only to avoid a name clash inside the flat
/// WASM namespace; the canonical name is exported as well (blocker C-H1).
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn solve_periodic_randnla(
    hamiltonian_blocks_json: &str,
    overlap_blocks_json: &str,
    n_electrons: usize,
    grid_json: &str,
    centering: &str,
    sketch_size: Option<usize>,
) -> String {
    solve_periodic_randnla_wasm(
        hamiltonian_blocks_json,
        overlap_blocks_json,
        n_electrons,
        grid_json,
        centering,
        sketch_size,
    )
}

/// Compute the Bloch phase factor for a lattice translation.
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn bloch_phase_wasm(k_x: f64, k_y: f64, k_z: f64, t_x: i32, t_y: i32, t_z: i32) -> String {
    let (cos_theta, sin_theta) = sci_form::kmesh::bloch_phase([k_x, k_y, k_z], [t_x, t_y, t_z]);
    serde_json::json!({
        "cos": cos_theta,
        "sin": sin_theta
    })
    .to_string()
}

/// Alias of [`bloch_phase`], the name of the Rust root function.
///
/// The `_wasm` suffix existed only to avoid a name clash inside the flat
/// WASM namespace; the canonical name is exported as well (blocker C-H1).
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn bloch_phase(k_x: f64, k_y: f64, k_z: f64, t_x: i32, t_y: i32, t_z: i32) -> String {
    bloch_phase_wasm(k_x, k_y, k_z, t_x, t_y, t_z)
}

/// Build the Bloch Hamiltonian `H(k)` from real-space operator blocks.
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn build_bloch_hamiltonian_wasm(blocks_json: &str, k_x: f64, k_y: f64, k_z: f64) -> String {
    let blocks = match parse_blocks(blocks_json) {
        Ok(blocks) => blocks,
        Err(e) => return json_error(&e),
    };
    let refs = as_core_refs(&blocks);
    let hamiltonian = sci_form::kmesh::build_bloch_hamiltonian(&refs, [k_x, k_y, k_z]);
    let n = hamiltonian.nrows();
    serde_json::json!({
        "n_basis": n,
        "k": [k_x, k_y, k_z],
        "matrix_shape": [n, n],
        "matrix": matrix_rows(&hamiltonian),
        "eigenvalues": symmetric_eigenvalues(&hamiltonian),
    })
    .to_string()
}

/// Alias of [`build_bloch_hamiltonian`], the name of the Rust root function.
///
/// The `_wasm` suffix existed only to avoid a name clash inside the flat
/// WASM namespace; the canonical name is exported as well (blocker C-H1).
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn build_bloch_hamiltonian(blocks_json: &str, k_x: f64, k_y: f64, k_z: f64) -> String {
    build_bloch_hamiltonian_wasm(blocks_json, k_x, k_y, k_z)
}

/// Assemble periodic Hamiltonians and overlaps for every k-point of a mesh.
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn assemble_periodic_operators_wasm(
    hamiltonian_blocks_json: &str,
    overlap_blocks_json: &str,
    grid_json: &str,
    centering: &str,
) -> String {
    use sci_form::alpha::periodic_linear::{monkhorst_pack_mesh, KMeshConfig};
    let hamiltonian_blocks = match parse_blocks(hamiltonian_blocks_json) {
        Ok(blocks) => blocks,
        Err(e) => return json_error(&e),
    };
    let overlap_blocks = if overlap_blocks_json.trim().is_empty() {
        Vec::new()
    } else {
        match parse_blocks(overlap_blocks_json) {
            Ok(blocks) => blocks,
            Err(e) => return json_error(&e),
        }
    };
    let grid: [usize; 3] = match serde_json::from_str(grid_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'grid'", grid_json, &e)),
    };
    let centering = match parse_centering(centering) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let mesh = match monkhorst_pack_mesh(&KMeshConfig { grid, centering }) {
        Ok(mesh) => mesh,
        Err(e) => return json_error(&e),
    };
    let hamiltonian_refs = as_core_refs(&hamiltonian_blocks);
    let overlap_refs = as_core_refs(&overlap_blocks);
    let (hamiltonians, overlaps) =
        sci_form::kmesh::assemble_periodic_operators(&hamiltonian_refs, &overlap_refs, &mesh);
    let n_basis = hamiltonians.first().map(|m| m.nrows()).unwrap_or_default();
    serde_json::json!({
        "n_kpoints": hamiltonians.len(),
        "n_basis": n_basis,
        "assembled": hamiltonians.len() == mesh.points.len(),
        "kpoints": mesh.points.iter().map(|p| serde_json::json!({
            "fractional": p.fractional,
            "index": p.index,
            "weight": p.weight,
        })).collect::<Vec<_>>(),
        "hamiltonians": hamiltonians.iter().map(matrix_rows).collect::<Vec<_>>(),
        "overlaps": overlaps.iter().map(matrix_rows).collect::<Vec<_>>(),
    })
    .to_string()
}

/// Alias of [`assemble_periodic_operators`], the name of the Rust root function.
///
/// The `_wasm` suffix existed only to avoid a name clash inside the flat
/// WASM namespace; the canonical name is exported as well (blocker C-H1).
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn assemble_periodic_operators(
    hamiltonian_blocks_json: &str,
    overlap_blocks_json: &str,
    grid_json: &str,
    centering: &str,
) -> String {
    assemble_periodic_operators_wasm(
        hamiltonian_blocks_json,
        overlap_blocks_json,
        grid_json,
        centering,
    )
}

/// Integrate a scalar observable over the Brillouin zone.
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn bz_integrate_scalar_wasm(values_json: &str, weights_json: &str) -> String {
    let values: Vec<f64> = match serde_json::from_str(values_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'values'", values_json, &e)),
    };
    let weights: Vec<f64> = match serde_json::from_str(weights_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'weights'", weights_json, &e)),
    };
    if values.len() != weights.len() {
        return json_error("values and weights must have the same length");
    }
    let result = sci_form::kmesh::bz_integrate_scalar(&values, &weights);
    serde_json::json!({"result": result}).to_string()
}

/// Alias of [`bz_integrate_scalar`], the name of the Rust root function.
///
/// The `_wasm` suffix existed only to avoid a name clash inside the flat
/// WASM namespace; the canonical name is exported as well (blocker C-H1).
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn bz_integrate_scalar(values_json: &str, weights_json: &str) -> String {
    bz_integrate_scalar_wasm(values_json, weights_json)
}

/// Integrate a vector observable over the Brillouin zone.
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn bz_integrate_vector_wasm(vectors_json: &str, weights_json: &str) -> String {
    let vectors: Vec<Vec<f64>> = match serde_json::from_str(vectors_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'vectors'", vectors_json, &e)),
    };
    let weights: Vec<f64> = match serde_json::from_str(weights_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'weights'", weights_json, &e)),
    };
    if vectors.len() != weights.len() {
        return json_error("vectors and weights must have the same length");
    }
    let result = sci_form::kmesh::bz_integrate_vector(&vectors, &weights);
    serde_json::json!({"result": result}).to_string()
}

/// Alias of [`bz_integrate_vector`], the name of the Rust root function.
///
/// The `_wasm` suffix existed only to avoid a name clash inside the flat
/// WASM namespace; the canonical name is exported as well (blocker C-H1).
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn bz_integrate_vector(vectors_json: &str, weights_json: &str) -> String {
    bz_integrate_vector_wasm(vectors_json, weights_json)
}

/// Check that the electron count fits into the available band/k-point states.
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

/// Alias of [`validate_electron_count`], the name of the Rust root function.
///
/// The `_wasm` suffix existed only to avoid a name clash inside the flat
/// WASM namespace; the canonical name is exported as well (blocker C-H1).
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn validate_electron_count(n_electrons: usize, n_bands: usize, n_kpoints: usize) -> String {
    validate_electron_count_wasm(n_electrons, n_bands, n_kpoints)
}

/// Verify that k-mesh weights sum to 1 within a tolerance.
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn validate_kmesh_weights_wasm(weights_json: &str, tolerance: f64) -> String {
    use sci_form::kmesh::{KMesh, KMeshCentering, KPoint};
    let weights: Vec<f64> = match serde_json::from_str(weights_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'weights'", weights_json, &e)),
    };
    let points: Vec<KPoint> = weights
        .iter()
        .map(|&weight| KPoint {
            fractional: [0.0, 0.0, 0.0],
            index: [0, 0, 0],
            weight,
        })
        .collect();
    let mesh = KMesh {
        total_weight: weights.iter().sum(),
        points,
        grid: [1, 1, 1],
        centering: KMeshCentering::MonkhorstPack,
    };
    serde_json::json!({
        "valid": sci_form::kmesh::validate_kmesh_weights(&mesh, tolerance),
        "total_weight": weights.iter().sum::<f64>()
    })
    .to_string()
}

/// Alias of [`validate_kmesh_weights`], the name of the Rust root function.
///
/// The `_wasm` suffix existed only to avoid a name clash inside the flat
/// WASM namespace; the canonical name is exported as well (blocker C-H1).
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn validate_kmesh_weights(weights_json: &str, tolerance: f64) -> String {
    validate_kmesh_weights_wasm(weights_json, tolerance)
}

// ─── KINETICS EXPORTS ──────────────────────────────────────────

/// Extract convergence diagnostics from a solved microkinetic trace.
///
/// `network_json` is an optional `MicrokineticNetworkConfig` JSON document;
/// when omitted the built-in A ⇌ B demo network is solved.
#[cfg(feature = "alpha-kinetics")]
#[wasm_bindgen]
pub fn extract_kinetics_diagnostics_wasm(network_json: &str, n_steps: usize) -> String {
    use sci_form::alpha::kinetics::*;
    let config = match network_config(network_json, n_steps) {
        Ok(config) => config,
        Err(e) => return json_error(&e),
    };
    let trace = match solve_microkinetic_network(&config) {
        Ok(trace) => trace,
        Err(e) => return json_error(&e),
    };
    let diagnostics = extract_kinetics_diagnostics(&trace);
    serde_json::json!({
        "max_population": diagnostics.max_population,
        "total_population": diagnostics.total_population,
        "mass_conservation_error": diagnostics.mass_conservation_error,
        "max_rate_of_change": diagnostics.max_rate_of_change,
        "is_steady_state": diagnostics.is_steady_state,
    })
    .to_string()
}

/// Alias of [`extract_kinetics_diagnostics`], the name of the Rust root function.
///
/// The `_wasm` suffix existed only to avoid a name clash inside the flat
/// WASM namespace; the canonical name is exported as well (blocker C-H1).
#[cfg(feature = "alpha-kinetics")]
#[wasm_bindgen]
pub fn extract_kinetics_diagnostics(network_json: &str, n_steps: usize) -> String {
    extract_kinetics_diagnostics_wasm(network_json, n_steps)
}

/// Solve a microkinetic network and report the resulting trace.
#[cfg(feature = "alpha-kinetics")]
#[wasm_bindgen]
pub fn solve_microkinetic_network_wasm(network_json: &str, n_steps: usize) -> String {
    use sci_form::alpha::kinetics::*;
    let config = match network_config(network_json, n_steps) {
        Ok(config) => config,
        Err(e) => return json_error(&e),
    };
    let trace = match solve_microkinetic_network(&config) {
        Ok(trace) => trace,
        Err(e) => return json_error(&e),
    };
    let diagnostics = extract_kinetics_diagnostics(&trace);
    let last = trace.frames.last();
    serde_json::json!({
        "n_steps": trace.frames.len(),
        "n_rates": trace.rates.len(),
        "converged": diagnostics.is_steady_state,
        "final_time_s": last.map(|f| f.time_s),
        "final_populations": last
            .map(|f| f.species.iter().map(|s| s.population).collect::<Vec<_>>())
            .unwrap_or_default(),
        "diagnostics": diagnostics,
    })
    .to_string()
}

/// Alias of [`solve_microkinetic_network`], the name of the Rust root function.
///
/// The `_wasm` suffix existed only to avoid a name clash inside the flat
/// WASM namespace; the canonical name is exported as well (blocker C-H1).
#[cfg(feature = "alpha-kinetics")]
#[wasm_bindgen]
pub fn solve_microkinetic_network(network_json: &str, n_steps: usize) -> String {
    solve_microkinetic_network_wasm(network_json, n_steps)
}

/// Solve a microkinetic network and propagate it to steady state.
#[cfg(feature = "alpha-kinetics")]
#[wasm_bindgen]
pub fn solve_microkinetic_steady_state_wasm(
    network_json: &str,
    n_steps: usize,
    tol: f64,
) -> String {
    use sci_form::alpha::kinetics::*;
    let config = match network_config(network_json, n_steps) {
        Ok(config) => config,
        Err(e) => return json_error(&e),
    };
    let trace = match solve_microkinetic_steady_state(&config, tol) {
        Ok(trace) => trace,
        Err(e) => return json_error(&e),
    };
    let diagnostics = extract_kinetics_diagnostics(&trace);
    let last = trace.frames.last();
    serde_json::json!({
        "n_steps": trace.frames.len(),
        "converged": diagnostics.is_steady_state,
        "residual": diagnostics.max_rate_of_change,
        "tol": tol,
        "final_time_s": last.map(|f| f.time_s),
        "final_populations": last
            .map(|f| f.species.iter().map(|s| s.population).collect::<Vec<_>>())
            .unwrap_or_default(),
        "diagnostics": diagnostics,
    })
    .to_string()
}

/// Alias of [`solve_microkinetic_steady_state`], the name of the Rust root function.
///
/// The `_wasm` suffix existed only to avoid a name clash inside the flat
/// WASM namespace; the canonical name is exported as well (blocker C-H1).
#[cfg(feature = "alpha-kinetics")]
#[wasm_bindgen]
pub fn solve_microkinetic_steady_state(network_json: &str, n_steps: usize, tol: f64) -> String {
    solve_microkinetic_steady_state_wasm(network_json, n_steps, tol)
}

/// Analyze one GSM + MBH + HTST coupled elementary step.
#[cfg(feature = "alpha-kinetics")]
#[wasm_bindgen]
pub fn analyze_gsm_mbh_htst_step_wasm(
    smiles: &str,
    reactant_coords_json: &str,
    product_coords_json: &str,
    method: &str,
    temperature_k: f64,
    pressure_bar: f64,
) -> String {
    use sci_form::alpha::gsm::GsmEnergyBackend;
    use sci_form::alpha::kinetics::*;
    let reactant: Vec<f64> = match serde_json::from_str(reactant_coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'reactant coords'", reactant_coords_json, &e)),
    };
    let product: Vec<f64> = match serde_json::from_str(product_coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'product coords'", product_coords_json, &e)),
    };
    let backend: GsmEnergyBackend = match method.parse() {
        Ok(backend) => backend,
        Err(_) => {
            return json_error(&invalid_value(
                "'method'",
                method,
                &[
                    "uff", "mmff94", "pm3", "xtb", "gfn1", "gfn2", "hf3c", "reaxff",
                ],
            ))
        }
    };
    let config = HtstAdapterConfig {
        step_id: "elementary-step".to_string(),
        state: ThermodynamicState {
            temperature_k,
            pressure_bar,
        },
        ..Default::default()
    };
    match analyze_gsm_mbh_htst_step_with_backend(smiles, &reactant, &product, backend, &config) {
        Ok(result) => serde_json::json!({
            "method": backend.as_str(),
            "step_id": result.rate.step_id,
            "converged": result.rate.forward_rate_s_inv.is_finite(),
            "barrier_kcal_mol": result.activation_energy_kcal_mol,
            "reverse_barrier_kcal_mol": result.reverse_barrier_kcal_mol,
            "ts_energy_kcal_mol": result.ts_energy_kcal_mol,
            "forward_rate_s_inv": result.rate.forward_rate_s_inv,
            "reverse_rate_s_inv": result.rate.reverse_rate_s_inv,
            "equilibrium_constant": result.rate.equilibrium_constant,
            "imaginary_frequency_cm1": result.imaginary_frequency_cm1,
            "n_blocks": result.n_blocks,
            "n_flexible": result.n_flexible,
            "n_nodes": result.n_nodes,
            "energy_evaluations": result.energy_evaluations,
        })
        .to_string(),
        Err(e) => json_error(&e),
    }
}

/// Alias of [`analyze_gsm_mbh_htst_step`], the name of the Rust root function.
///
/// The `_wasm` suffix existed only to avoid a name clash inside the flat
/// WASM namespace; the canonical name is exported as well (blocker C-H1).
#[cfg(feature = "alpha-kinetics")]
#[wasm_bindgen]
pub fn analyze_gsm_mbh_htst_step(
    smiles: &str,
    reactant_coords_json: &str,
    product_coords_json: &str,
    method: &str,
    temperature_k: f64,
    pressure_bar: f64,
) -> String {
    analyze_gsm_mbh_htst_step_wasm(
        smiles,
        reactant_coords_json,
        product_coords_json,
        method,
        temperature_k,
        pressure_bar,
    )
}

// ──── RENDER BRIDGE: chart utilities ──────────────────────────────────────

/// Serialize a chart payload to the canonical wrapper JSON.
#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn chart_to_json_wasm(chart_json: &str) -> String {
    let chart = match parse_chart(chart_json) {
        Ok(chart) => chart,
        Err(e) => return json_error(&e),
    };
    serde_json::json!({
        "format": "json",
        "version": "1.0",
        "data": chart,
    })
    .to_string()
}

/// Alias of [`chart_to_json`], the name of the Rust root function.
///
/// The `_wasm` suffix existed only to avoid a name clash inside the flat
/// WASM namespace; the canonical name is exported as well (blocker C-H1).
#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn chart_to_json(chart_json: &str) -> String {
    chart_to_json_wasm(chart_json)
}

/// Deserialize a chart payload (wrapper or bare chart) to chart JSON.
#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn chart_from_json_wasm(json_str: &str) -> String {
    match parse_chart(json_str) {
        Ok(chart) => serde_json::to_string(&chart)
            .unwrap_or_else(|e| serde_json::json!({"error": e.to_string()}).to_string()),
        Err(e) => serde_json::json!({"error": e}).to_string(),
    }
}

/// Alias of [`chart_from_json`], the name of the Rust root function.
///
/// The `_wasm` suffix existed only to avoid a name clash inside the flat
/// WASM namespace; the canonical name is exported as well (blocker C-H1).
#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn chart_from_json(json_str: &str) -> String {
    chart_from_json_wasm(json_str)
}

/// Validate the structure of a chart payload (wrapper or bare chart).
#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn validate_chart_schema_wasm(json_str: &str) -> bool {
    parse_chart(json_str).is_ok()
}

/// Alias of [`validate_chart_schema`], the name of the Rust root function.
///
/// The `_wasm` suffix existed only to avoid a name clash inside the flat
/// WASM namespace; the canonical name is exported as well (blocker C-H1).
#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
pub fn validate_chart_schema(json_str: &str) -> bool {
    validate_chart_schema_wasm(json_str)
}
