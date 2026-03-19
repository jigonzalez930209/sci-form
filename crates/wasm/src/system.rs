//! System info and method comparison WASM bindings.

use crate::helpers::*;
use wasm_bindgen::prelude::*;

/// Library version.
#[wasm_bindgen]
pub fn version() -> String {
    sci_form::version()
}

/// Query operation capabilities from a JSON element array.
#[wasm_bindgen]
pub fn system_capabilities(elements: &str) -> String {
    let elems = match parse_elements(elements) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    serialize_or_error(&sci_form::get_system_capabilities(&elems))
}

/// Build a structured method plan from a JSON element array.
#[wasm_bindgen]
pub fn system_method_plan(elements: &str) -> String {
    let elems = match parse_elements(elements) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    serialize_or_error(&sci_form::get_system_method_plan(&elems))
}

/// Compare available methods on the same system/geometry.
#[wasm_bindgen]
pub fn compare_methods(
    smiles: &str,
    elements: &str,
    coords_flat: &str,
    allow_experimental_eht: bool,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    serialize_or_error(&sci_form::compare_methods(
        smiles,
        &elems,
        &positions,
        allow_experimental_eht,
    ))
}
