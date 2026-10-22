//! System info and method comparison WASM bindings.

use crate::helpers::*;
use wasm_bindgen::prelude::*;

/// Library version.
#[wasm_bindgen]
pub fn version() -> String {
    sci_form::version()
}

/// Query operation capabilities from a JSON element array.
///
/// The Rust root calls this `get_system_capabilities`; the `get_` prefix is
/// not idiomatic and is dropped here (blocker C-M7).  The Rust name is kept
/// as an alias below.
#[wasm_bindgen]
pub fn system_capabilities(elements: &str) -> String {
    let elems = match parse_elements(elements) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    serialize_or_error(&sci_form::system_capabilities(&elems))
}

/// Query operation capabilities from a JSON element array.
///
/// **Alias** of [`system_capabilities`] matching the Rust root name
/// `sci_form::get_system_capabilities`.
#[wasm_bindgen]
#[allow(deprecated)] // this alias exists to keep the old name working
pub fn get_system_capabilities(elements: &str) -> String {
    let elems = match parse_elements(elements) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    serialize_or_error(&sci_form::get_system_capabilities(&elems))
}

/// Build a structured method plan from a JSON element array.
///
/// The Rust root calls this `get_system_method_plan`; the `get_` prefix is
/// not idiomatic and is dropped here (blocker C-M7).  The Rust name is kept
/// as an alias below.
#[wasm_bindgen]
pub fn system_method_plan(elements: &str) -> String {
    let elems = match parse_elements(elements) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    serialize_or_error(&sci_form::system_method_plan(&elems))
}

/// Build a structured method plan from a JSON element array.
///
/// **Alias** of [`system_method_plan`] matching the Rust root name
/// `sci_form::get_system_method_plan`.
#[wasm_bindgen]
#[allow(deprecated)] // this alias exists to keep the old name working
pub fn get_system_method_plan(elements: &str) -> String {
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
