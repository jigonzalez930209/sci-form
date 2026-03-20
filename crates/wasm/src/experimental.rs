//! WASM bindings for experimental modules.
//!
//! Each function is gated behind its experimental feature flag.
//! All functions accept/return JSON strings.

use crate::helpers::{json_error, parse_elements_and_positions};
use wasm_bindgen::prelude::*;

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
    let config = sci_form::experimental::eeq::EeqConfig {
        total_charge,
        regularization: 1e-10,
    };
    let r = sci_form::experimental::eeq::compute_eeq_charges(&elems, &pos, &config);
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
    let config = sci_form::experimental::eeq::EeqConfig {
        total_charge,
        regularization: 1e-10,
    };
    let r = sci_form::experimental::eeq::compute_eeq_energy(&elems, &pos, &config);
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
    let config = sci_form::experimental::alpb::AlpbConfig {
        solvent_dielectric,
        probe_radius,
        surface_tension,
    };
    let r = sci_form::experimental::alpb::compute_alpb_solvation(&elems, &pos, &q, &config);
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
    let r = sci_form::experimental::alpb::compute_born_radii(&elems, &pos, probe_radius);
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
    let config = sci_form::experimental::d4::D4Config {
        s6, s8, a1, a2, three_body, s9: 1.0,
    };
    let r = sci_form::experimental::d4::compute_d4_energy(&elems, &pos, &config);
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
    let config = sci_form::experimental::cpm::CpmConfig {
        mu_ev,
        dielectric,
        max_iter: 100,
        charge_tol: 1e-6,
    };
    let r = sci_form::experimental::cpm::compute_cpm_charges(&elems, &pos, &config);
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
    let r = sci_form::experimental::cpm::compute_cpm_surface(
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
