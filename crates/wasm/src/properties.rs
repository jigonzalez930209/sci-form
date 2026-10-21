//! Molecular property calculations WASM bindings.
//! Charges, SASA, population, bond orders, dipole, DOS, RMSD.

use crate::helpers::*;
use wasm_bindgen::prelude::*;

/// Compute Gasteiger-Marsili partial charges from SMILES.
#[wasm_bindgen]
pub fn compute_charges(smiles: &str) -> String {
    match sci_form::compute_charges(smiles) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Compute solvent-accessible surface area.
#[wasm_bindgen]
pub fn compute_sasa(elements: &str, coords_flat: &str, probe_radius: f64) -> String {
    let elems = match parse_elements(elements) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let flat = match parse_flat_coords(coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let pr = if probe_radius <= 0.0 {
        None
    } else {
        Some(probe_radius)
    };
    match sci_form::compute_sasa(&elems, &flat, pr) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Compute Mulliken & Löwdin population analysis.
#[wasm_bindgen]
pub fn compute_population(elements: &str, coords_flat: &str) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_population(&elems, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Compute Wiberg-like and Mayer-like bond orders.
#[wasm_bindgen]
pub fn compute_bond_orders(elements: &str, coords_flat: &str) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_bond_orders(&elems, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Compute molecular dipole moment.
#[wasm_bindgen]
pub fn compute_dipole(elements: &str, coords_flat: &str) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_dipole(&elems, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Compute DOS/PDOS from EHT orbital energies.
#[wasm_bindgen]
pub fn compute_dos(
    elements: &str,
    coords_flat: &str,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_dos(&elems, &positions, sigma, e_min, e_max, n_points) {
        Ok(result) => {
            let energies_json: Vec<String> = result
                .energies
                .iter()
                .map(|v| format!("{:.4}", v))
                .collect();
            let dos_json: Vec<String> = result
                .total_dos
                .iter()
                .map(|v| format!("{:.6}", v))
                .collect();
            format!(
                "{{\"sigma\":{},\"energies\":[{}],\"total_dos\":[{}],\"n_atoms_pdos\":{}}}",
                result.sigma,
                energies_json.join(","),
                dos_json.join(","),
                result.pdos.len()
            )
        }
        Err(e) => json_error(&e),
    }
}

/// Compute Gasteiger-Marsili partial charges with configurable parameters.
///
/// `max_iter`: number of electronegativity-equalization sweeps.
/// `initial_damping`: first damping factor (halved every iteration).
/// `convergence_threshold`: max charge-transfer convergence criterion.
#[wasm_bindgen]
pub fn compute_charges_configured(
    smiles: &str,
    max_iter: usize,
    initial_damping: f64,
    convergence_threshold: f64,
) -> String {
    let config = sci_form::charges::gasteiger::GasteigerConfig {
        max_iter,
        initial_damping,
        convergence_threshold,
    };
    match sci_form::compute_charges_configured(smiles, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Compute DOS from the orbital energies of a chosen electronic-structure method.
///
/// `method`: `eht`, `pm3`, `xtb` (GFN2), `gfn1`, `gfn2` or `hf3c`.
#[wasm_bindgen]
pub fn compute_dos_multimethod(
    elements: &str,
    coords_flat: &str,
    method: &str,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let dos_method = match method.to_ascii_lowercase().as_str() {
        "eht" => sci_form::dos::multi_method::DosMethod::Eht,
        "pm3" => sci_form::dos::multi_method::DosMethod::Pm3,
        "xtb" | "gfn0" => sci_form::dos::multi_method::DosMethod::Xtb,
        "gfn1" => sci_form::dos::multi_method::DosMethod::Gfn1,
        "gfn2" => sci_form::dos::multi_method::DosMethod::Gfn2,
        "hf3c" | "hf-3c" => sci_form::dos::multi_method::DosMethod::Hf3c,
        other => {
            return json_error(&format!(
                "unknown DOS method '{other}' (expected eht, pm3, xtb, gfn1, gfn2 or hf3c)"
            ))
        }
    };
    match sci_form::dos::multi_method::compute_dos_multimethod(
        &elems, &positions, dos_method, sigma, e_min, e_max, n_points,
    ) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Compute RMSD between two coordinate sets after Kabsch alignment.
#[wasm_bindgen]
pub fn compute_rmsd(coords: &str, reference: &str) -> String {
    let c: Vec<f64> = match serde_json::from_str(coords) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad coords: {}", e)),
    };
    let r: Vec<f64> = match serde_json::from_str(reference) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad reference: {}", e)),
    };
    let result = sci_form::alignment::align_coordinates(&c, &r);
    format!("{{\"rmsd\":{:.6}}}", result.rmsd)
}
