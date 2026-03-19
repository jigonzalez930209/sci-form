//! EHT (Extended Hückel Theory) WASM bindings.

use crate::helpers::*;
use wasm_bindgen::prelude::*;

fn orbital_grid_from_coefficients(
    elements: &str,
    coords_flat: &str,
    coefficients_json: &str,
    mo_index: usize,
    spacing: f64,
) -> Result<sci_form::eht::VolumetricGrid, String> {
    let (elems, positions) = parse_elements_and_positions(elements, coords_flat)?;
    let coefficients: Vec<Vec<f64>> =
        serde_json::from_str(coefficients_json).map_err(|e| format!("bad coefficients: {}", e))?;
    let basis = sci_form::eht::basis::build_basis(&elems, &positions);

    if basis.is_empty() {
        return Err("No basis functions found for orbital evaluation".to_string());
    }
    if mo_index >= coefficients.len() {
        return Err(format!(
            "orbital index {} out of range for {} orbitals",
            mo_index,
            coefficients.len()
        ));
    }
    if coefficients.len() != basis.len() {
        return Err(format!(
            "coefficient row count {} does not match basis size {}",
            coefficients.len(),
            basis.len()
        ));
    }
    if coefficients.iter().any(|row| mo_index >= row.len()) {
        return Err(format!(
            "orbital index {} exceeds coefficient columns",
            mo_index
        ));
    }

    #[cfg(feature = "parallel")]
    {
        Ok(sci_form::eht::evaluate_orbital_on_grid_parallel(
            &basis,
            &coefficients,
            mo_index,
            &positions,
            spacing,
            3.0,
        ))
    }

    #[cfg(not(feature = "parallel"))]
    {
        Ok(sci_form::eht::evaluate_orbital_on_grid(
            &basis,
            &coefficients,
            mo_index,
            &positions,
            spacing,
            3.0,
        ))
    }
}

/// Run an EHT calculation on a molecule.
///
/// - `elements`: JSON array of atomic numbers
/// - `coords_flat`: JSON array of flat xyz coords
/// - `k`: Wolfsberg-Helmholtz constant (0.0 for default 1.75)
#[wasm_bindgen]
pub fn eht_calculate(elements: &str, coords_flat: &str, k: f64) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let k_opt = if k <= 0.0 { None } else { Some(k) };
    match sci_form::eht::solve_eht(&elems, &positions, k_opt) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Run EHT or route to UFF fallback based on support/confidence.
#[wasm_bindgen]
pub fn eht_or_uff_fallback(
    smiles: &str,
    elements: &str,
    coords_flat: &str,
    allow_experimental_eht: bool,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_eht_or_uff_fallback(smiles, &elems, &positions, allow_experimental_eht)
    {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Generate an orbital isosurface mesh from EHT.
#[wasm_bindgen]
pub fn eht_orbital_mesh(
    elements: &str,
    coords_flat: &str,
    mo_index: usize,
    spacing: f64,
    isovalue: f32,
) -> String {
    let result_json = eht_calculate(elements, coords_flat, 0.0);
    let result: sci_form::eht::EhtResult = match serde_json::from_str(&result_json) {
        Ok(result) => result,
        Err(_) => return result_json,
    };
    let coeff_json = match serde_json::to_string(&result.coefficients) {
        Ok(json) => json,
        Err(e) => return json_error(&e.to_string()),
    };
    let grid =
        match orbital_grid_from_coefficients(elements, coords_flat, &coeff_json, mo_index, spacing)
        {
            Ok(grid) => grid,
            Err(e) => return json_error(&e),
        };
    let mesh = sci_form::eht::marching_cubes(&grid, isovalue);
    serialize_or_error(&mesh)
}

/// Query EHT support metadata.
#[wasm_bindgen]
pub fn eht_support(elements: &str) -> String {
    let elems: Vec<u8> = match parse_elements(elements) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    serialize_or_error(&sci_form::get_eht_support(&elems))
}

/// Compute orbital grid as Float32Array (typed-array transfer).
#[wasm_bindgen]
pub fn eht_orbital_grid_typed(
    elements: &str,
    coords_flat: &str,
    mo_index: usize,
    spacing: f64,
) -> js_sys::Float32Array {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(_) => return js_sys::Float32Array::new_with_length(0),
    };
    let result = match sci_form::eht::solve_eht(&elems, &positions, None) {
        Ok(r) => r,
        Err(_) => return js_sys::Float32Array::new_with_length(0),
    };
    let coeff_json = match serde_json::to_string(&result.coefficients) {
        Ok(json) => json,
        Err(_) => return js_sys::Float32Array::new_with_length(0),
    };
    let grid =
        match orbital_grid_from_coefficients(elements, coords_flat, &coeff_json, mo_index, spacing)
        {
            Ok(grid) => grid,
            Err(_) => return js_sys::Float32Array::new_with_length(0),
        };
    let arr = js_sys::Float32Array::new_with_length(grid.values.len() as u32);
    arr.copy_from(&grid.values);
    arr
}

/// Compute orbital grid from precomputed EHT coefficients as Float32Array.
#[wasm_bindgen]
pub fn eht_orbital_grid_from_coefficients_typed(
    elements: &str,
    coords_flat: &str,
    coefficients_json: &str,
    mo_index: usize,
    spacing: f64,
) -> js_sys::Float32Array {
    let grid = match orbital_grid_from_coefficients(
        elements,
        coords_flat,
        coefficients_json,
        mo_index,
        spacing,
    ) {
        Ok(grid) => grid,
        Err(_) => return js_sys::Float32Array::new_with_length(0),
    };
    let arr = js_sys::Float32Array::new_with_length(grid.values.len() as u32);
    arr.copy_from(&grid.values);
    arr
}
