//! Spectroscopy WASM bindings: sTDA UV-Vis, IR, NMR.

use crate::helpers::*;
use wasm_bindgen::prelude::*;

/// Compute an sTDA UV-Vis spectrum with proper oscillator strengths.
#[wasm_bindgen]
pub fn compute_stda_uvvis(
    elements_json: &str,
    coords_flat_json: &str,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
    broadening: &str,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let bt = match broadening {
        "lorentzian" | "Lorentzian" => sci_form::reactivity::BroadeningType::Lorentzian,
        _ => sci_form::reactivity::BroadeningType::Gaussian,
    };
    match sci_form::compute_stda_uvvis(&elems, &positions, sigma, e_min, e_max, n_points, bt) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Perform vibrational analysis via numerical Hessian.
#[wasm_bindgen]
pub fn compute_vibrational_analysis(
    elements_json: &str,
    coords_flat_json: &str,
    method: &str,
    step_size: f64,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let step = if step_size > 0.0 {
        Some(step_size)
    } else {
        None
    };
    match sci_form::compute_vibrational_analysis(&elems, &positions, method, step) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Perform vibrational analysis using the fast UFF analytical Hessian path.
#[wasm_bindgen]
pub fn compute_vibrational_analysis_uff(
    smiles: &str,
    elements_json: &str,
    coords_flat_json: &str,
    step_size: f64,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let step = if step_size > 0.0 {
        Some(step_size)
    } else {
        None
    };
    match sci_form::compute_vibrational_analysis_uff(smiles, &elems, &positions, step) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Generate a Lorentzian-broadened IR spectrum from vibrational analysis JSON.
#[wasm_bindgen]
pub fn compute_ir_spectrum(
    analysis_json: &str,
    gamma: f64,
    wn_min: f64,
    wn_max: f64,
    n_points: usize,
) -> String {
    let analysis: sci_form::ir::VibrationalAnalysis = match serde_json::from_str(analysis_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad analysis JSON: {}", e)),
    };
    let result = sci_form::compute_ir_spectrum(&analysis, gamma, wn_min, wn_max, n_points);
    serialize_or_error(&result)
}

/// Predict NMR chemical shifts for all active nuclei.
#[wasm_bindgen]
pub fn predict_nmr_shifts(smiles: &str) -> String {
    match sci_form::predict_nmr_shifts(smiles) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Predict J-coupling constants.
#[wasm_bindgen]
pub fn predict_nmr_couplings(smiles: &str, coords_flat_json: &str) -> String {
    let flat: Vec<f64> = match serde_json::from_str(coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad coords: {}", e)),
    };
    let positions: Vec<[f64; 3]> = flat.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    match sci_form::predict_nmr_couplings(smiles, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Generate a complete NMR spectrum.
///
/// `nucleus`: "1H", "13C", "19F", "31P", "15N", "11B", "29Si", "77Se", "17O", or "33S".
#[wasm_bindgen]
pub fn compute_nmr_spectrum(
    smiles: &str,
    nucleus: &str,
    gamma: f64,
    ppm_min: f64,
    ppm_max: f64,
    n_points: usize,
) -> String {
    match sci_form::compute_nmr_spectrum(smiles, nucleus, gamma, ppm_min, ppm_max, n_points) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Generate a complete NMR spectrum using optional 3D coordinates for ³J couplings.
#[wasm_bindgen]
pub fn compute_nmr_spectrum_with_coords(
    smiles: &str,
    coords_flat_json: &str,
    nucleus: &str,
    gamma: f64,
    ppm_min: f64,
    ppm_max: f64,
    n_points: usize,
) -> String {
    let flat: Vec<f64> = match serde_json::from_str(coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad coords: {}", e)),
    };
    let positions: Vec<[f64; 3]> = flat.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    match sci_form::compute_nmr_spectrum_with_coords(
        smiles, &positions, nucleus, gamma, ppm_min, ppm_max, n_points,
    ) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Generate HOSE codes for all atoms in a molecule.
#[wasm_bindgen]
pub fn compute_hose_codes(smiles: &str, max_radius: usize) -> String {
    match sci_form::compute_hose_codes(smiles, max_radius) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}
