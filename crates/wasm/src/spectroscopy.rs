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

/// Fluorescence from the same sTDA singlets as absorption.
/// `stokes_ev` is an extra red-shift; 0 keeps the vertical emission energy.
#[wasm_bindgen]
pub fn compute_fluorescence(
    elements_json: &str,
    coords_flat_json: &str,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
    stokes_ev: f64,
    refractive_index: f64,
    dark_threshold: f64,
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
    let config = sci_form::spectroscopy::FluorescenceConfig {
        sigma,
        e_min,
        e_max,
        n_points,
        stokes_ev,
        refractive_index,
        dark_threshold,
        broadening: bt,
    };
    match sci_form::compute_fluorescence(&elems, &positions, config) {
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

/// Generate a broadened IR spectrum from an internally computed vibrational analysis.
///
/// `method`: semi-empirical Hessian method (`eht`, `pm3`, `xtb`, …) or `uff`.
/// `broadening`: `gaussian` or `lorentzian`.
/// `smiles`: required when `method` is `uff`, otherwise ignored.
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn compute_ir_spectrum_broadened(
    elements_json: &str,
    coords_flat_json: &str,
    method: &str,
    gamma: f64,
    wn_min: f64,
    wn_max: f64,
    n_points: usize,
    broadening: &str,
    step_size: f64,
    smiles: &str,
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
    let analysis = if method.eq_ignore_ascii_case("uff") {
        if smiles.is_empty() {
            return json_error("UFF IR spectrum requires a smiles argument");
        }
        sci_form::compute_vibrational_analysis_uff(smiles, &elems, &positions, step)
    } else {
        sci_form::compute_vibrational_analysis(&elems, &positions, method, step)
    };
    match analysis {
        Ok(a) => {
            let spectrum = sci_form::compute_ir_spectrum_broadened(
                &a, gamma, wn_min, wn_max, n_points, broadening,
            );
            serialize_or_error(&spectrum)
        }
        Err(e) => json_error(&e),
    }
}

/// Assign IR peaks to functional groups using standard group-frequency tables.
///
/// `frequencies_json` / `intensities_json`: JSON arrays of numbers.
/// `elements_json`: JSON array of atomic numbers.
/// `mode_displacements_json`: JSON array of modes, each an array of `[x,y,z]`
/// triples; pass `""` to skip the normal-mode refinement.
#[wasm_bindgen]
pub fn assign_peaks(
    frequencies_json: &str,
    intensities_json: &str,
    elements_json: &str,
    mode_displacements_json: &str,
) -> String {
    let frequencies: Vec<f64> = match serde_json::from_str(frequencies_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad frequencies: {e}")),
    };
    let intensities: Vec<f64> = match serde_json::from_str(intensities_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad intensities: {e}")),
    };
    let elements: Vec<u8> = match serde_json::from_str(elements_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad elements: {e}")),
    };
    let displacements: Option<Vec<Vec<[f64; 3]>>> = if mode_displacements_json.is_empty() {
        None
    } else {
        match serde_json::from_str::<Vec<Vec<[f64; 3]>>>(mode_displacements_json) {
            Ok(v) => Some(v),
            Err(e) => return json_error(&format!("bad mode displacements: {e}")),
        }
    };
    let result = sci_form::ir::peak_assignment::assign_peaks(
        &frequencies,
        &intensities,
        &elements,
        displacements.as_deref(),
    );
    serialize_or_error(&result)
}

/// Predict J-couplings averaged over a conformer ensemble (Boltzmann weighting).
///
/// `conformer_coords_json`: JSON array of conformers, each a flat xyz array.
/// `energies_kcal_json`: JSON array of conformer energies in kcal/mol.
#[wasm_bindgen]
pub fn compute_ensemble_j_couplings(
    smiles: &str,
    conformer_coords_json: &str,
    energies_kcal_json: &str,
    temperature_k: f64,
) -> String {
    let conformers: Vec<Vec<f64>> = match serde_json::from_str(conformer_coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad conformer coords: {e}")),
    };
    let energies: Vec<f64> = match serde_json::from_str(energies_kcal_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad conformer energies: {e}")),
    };
    match sci_form::compute_ensemble_j_couplings(smiles, &conformers, &energies, temperature_k) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Predict NMR chemical shifts for all active nuclei.
#[wasm_bindgen]
pub fn predict_nmr_shifts(smiles: &str) -> String {
    match sci_form::predict_nmr_shifts(smiles) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Predict NMR chemical shifts for a specific nucleus.
#[wasm_bindgen]
pub fn predict_nmr_shifts_for_nucleus(smiles: &str, nucleus: &str) -> String {
    match sci_form::predict_nmr_shifts_for_nucleus(smiles, nucleus) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Run the public SCF-backed GIAO NMR route for a specific nucleus.
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn compute_giao_nmr(
    elements_json: &str,
    coords_flat_json: &str,
    nucleus: &str,
    charge: i32,
    multiplicity: u32,
    max_scf_iter: usize,
    allow_basis_fallback: bool,
    use_parallel_eri: bool,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };

    let config = sci_form::GiaoNmrConfig {
        charge,
        multiplicity,
        max_scf_iterations: max_scf_iter,
        use_parallel_eri,
        allow_basis_fallback,
    };

    match sci_form::compute_giao_nmr_configured(&elems, &positions, nucleus, &config) {
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
/// `nucleus`: any supported nucleus alias, such as "1H", "13C", "35Cl", "79Br", or "195Pt".
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
