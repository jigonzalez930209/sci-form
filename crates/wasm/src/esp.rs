//! ESP (electrostatic potential) WASM bindings.

use crate::helpers::*;
use wasm_bindgen::prelude::*;

/// Compute the full electrostatic potential grid as JSON.
#[wasm_bindgen]
pub fn compute_esp(elements: &str, coords_flat: &str, spacing: f64, padding: f64) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_esp(&elems, &positions, spacing, padding) {
        Ok(grid) => serialize_or_error(&grid),
        Err(e) => json_error(&e),
    }
}

/// Compute ESP grid values as Float64Array (typed-array transfer).
#[wasm_bindgen]
pub fn compute_esp_grid_typed(
    elements: &str,
    coords_flat: &str,
    spacing: f64,
    padding: f64,
) -> js_sys::Float64Array {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(_) => return js_sys::Float64Array::new_with_length(0),
    };
    match sci_form::compute_esp(&elems, &positions, spacing, padding) {
        Ok(grid) => {
            let arr = js_sys::Float64Array::new_with_length(grid.values.len() as u32);
            arr.copy_from(&grid.values);
            arr
        }
        Err(_) => js_sys::Float64Array::new_with_length(0),
    }
}

/// Get ESP grid metadata (origin, spacing, dimensions) as JSON.
#[wasm_bindgen]
pub fn compute_esp_grid_info(
    elements: &str,
    coords_flat: &str,
    spacing: f64,
    padding: f64,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_esp(&elems, &positions, spacing, padding) {
        Ok(grid) => format!(
            "{{\"origin\":[{:.4},{:.4},{:.4}],\"spacing\":{:.4},\"dims\":[{},{},{}]}}",
            grid.origin[0],
            grid.origin[1],
            grid.origin[2],
            grid.spacing,
            grid.dims[0],
            grid.dims[1],
            grid.dims[2]
        ),
        Err(e) => json_error(&e),
    }
}
