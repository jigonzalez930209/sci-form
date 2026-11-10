//! WASM bindings for EPR, MP2 and XPS.

use crate::helpers::{bad_json, json_error, parse_elements_and_positions, serialize_or_error};
use wasm_bindgen::prelude::*;

/// Isotropic EPR hyperfine coupling. `nucleus = 0` evaluates all known nuclei.
#[wasm_bindgen]
pub fn compute_epr_hyperfine(
    elements: &str,
    coords_flat: &str,
    nucleus: u8,
    charge: i32,
    multiplicity: u32,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(value) => value,
        Err(error) => return json_error(&error),
    };
    let config = sci_form::epr::EprConfig {
        nucleus: if nucleus == 0 { None } else { Some(nucleus) },
        charge,
        multiplicity,
    };
    match sci_form::compute_epr_hyperfine(&elems, &positions, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// Conventional closed-shell MP2 correlation energy.
#[wasm_bindgen]
pub fn compute_mp2_energy(elements: &str, coords_flat: &str) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(value) => value,
        Err(error) => return json_error(&error),
    };
    match sci_form::compute_mp2_energy(&elems, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// Koopmans core-level XPS binding energies.
#[wasm_bindgen]
pub fn compute_xps_binding_energies(
    elements: &str,
    coords_flat: &str,
    elements_filter_json: &str,
    charge: i32,
    multiplicity: u32,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(value) => value,
        Err(error) => return json_error(&error),
    };
    let filter = if elements_filter_json.trim().is_empty() {
        None
    } else {
        match serde_json::from_str::<Vec<u8>>(elements_filter_json) {
            Ok(value) => Some(value),
            Err(error) => {
                return json_error(&bad_json(
                    "'elements_filter_json'",
                    elements_filter_json,
                    &error,
                ))
            }
        }
    };
    let config = sci_form::xps::XpsConfig {
        charge,
        multiplicity,
        elements: filter,
    };
    match sci_form::compute_xps_binding_energies(&elems, &positions, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// GCMC adsorption isotherm from JSON framework, guest and config payloads.
#[wasm_bindgen]
pub fn compute_gcmc_isotherm(framework_json: &str, guest_json: &str, config_json: &str) -> String {
    let framework: sci_form::gcmc::GcmcFramework = match serde_json::from_str(framework_json) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'framework_json'", framework_json, &error)),
    };
    let guest: sci_form::gcmc::GuestModel = match serde_json::from_str(guest_json) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'guest_json'", guest_json, &error)),
    };
    let config: sci_form::gcmc::GcmcConfig = match serde_json::from_str(config_json) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'config_json'", config_json, &error)),
    };
    match sci_form::compute_gcmc_isotherm(&framework, &guest, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}
