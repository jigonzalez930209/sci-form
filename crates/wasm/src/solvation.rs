//! WASM bindings for solvation energy calculations.

use crate::helpers::{json_error, parse_elements_and_positions, serialize_or_error};
use wasm_bindgen::prelude::*;

/// Non-polar solvation energy from SASA with atomic solvation parameters.
///
/// Returns JSON `NonPolarSolvation`.
#[wasm_bindgen]
pub fn compute_nonpolar_solvation(elements: &str, coords_flat: &str, probe_radius: f64) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let r = sci_form::compute_nonpolar_solvation(&elems, &positions, Some(probe_radius));
    serialize_or_error(&r)
}

/// Generalized Born solvation energy (electrostatic + non-polar).
///
/// `charges` is a JSON array of partial charges.
/// Returns JSON `GbSolvation`.
#[wasm_bindgen]
pub fn compute_gb_solvation(
    elements: &str,
    coords_flat: &str,
    charges: &str,
    solvent_dielectric: f64,
    solute_dielectric: f64,
    probe_radius: f64,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let q: Vec<f64> = match serde_json::from_str(charges) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad charges: {e}")),
    };
    let r = sci_form::compute_gb_solvation(
        &elems,
        &positions,
        &q,
        Some(solvent_dielectric),
        Some(solute_dielectric),
        Some(probe_radius),
    );
    serialize_or_error(&r)
}
