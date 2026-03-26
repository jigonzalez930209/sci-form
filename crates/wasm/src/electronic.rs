//! Electronic structure method WASM bindings: PM3, xTB, GFN1, GFN2, HF-3c, ANI.

use crate::helpers::*;
use wasm_bindgen::prelude::*;

/// Run a PM3 semi-empirical calculation.
#[wasm_bindgen]
pub fn compute_pm3(elements_json: &str, coords_flat_json: &str) -> String {
    let (elements, positions) = match parse_elements_and_positions(elements_json, coords_flat_json)
    {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_pm3(&elements, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Run an xTB tight-binding calculation.
#[wasm_bindgen]
pub fn compute_xtb(elements_json: &str, coords_flat_json: &str) -> String {
    let (elements, positions) = match parse_elements_and_positions(elements_json, coords_flat_json)
    {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_xtb(&elements, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Run a GFN1-xTB tight-binding calculation.
#[wasm_bindgen]
pub fn compute_gfn1(elements_json: &str, coords_flat_json: &str) -> String {
    let (elements, positions) = match parse_elements_and_positions(elements_json, coords_flat_json)
    {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::xtb::solve_gfn1(&elements, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Run a GFN2-xTB tight-binding calculation.
#[wasm_bindgen]
pub fn compute_gfn2(elements_json: &str, coords_flat_json: &str) -> String {
    let (elements, positions) = match parse_elements_and_positions(elements_json, coords_flat_json)
    {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::xtb::solve_gfn2(&elements, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Run a HF-3c calculation with default configuration.
#[wasm_bindgen]
pub fn compute_hf3c(elements_json: &str, coords_flat_json: &str) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let config = sci_form::hf::HfConfig::default();
    match sci_form::compute_hf3c(&elems, &positions, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Run a HF-3c calculation with custom configuration.
#[wasm_bindgen]
pub fn compute_hf3c_custom(
    elements_json: &str,
    coords_flat_json: &str,
    max_scf_iter: usize,
    n_cis_states: usize,
    corrections: bool,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let config = sci_form::hf::HfConfig {
        max_scf_iter,
        diis_size: 6,
        n_cis_states,
        corrections,
    };
    match sci_form::compute_hf3c(&elems, &positions, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Compute ANI neural-network potential energy.
#[wasm_bindgen]
pub fn compute_ani(elements_json: &str, coords_flat_json: &str) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_ani(&elems, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}
