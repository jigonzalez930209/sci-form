//! Orbital mesh generation WASM bindings for all methods.

use crate::helpers::*;
use wasm_bindgen::prelude::*;

/// Compute an orbital isosurface mesh using any implemented method.
///
/// `method`: "eht", "pm3", "xtb", or "hf3c".
/// `mo_index`: molecular orbital index (0-based).
/// `spacing`: grid spacing in Å (e.g. 0.2).
/// `isovalue`: isosurface cutoff (e.g. 0.02).
///
/// Returns JSON OrbitalMeshResult with mesh, grid, orbital energies, HOMO-LUMO gap.
#[wasm_bindgen]
pub fn compute_orbital_mesh(
    elements_json: &str,
    coords_flat_json: &str,
    method: &str,
    mo_index: usize,
    spacing: f64,
    padding: f64,
    isovalue: f32,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_orbital_mesh(
        &elems, &positions, method, mo_index, spacing, padding, isovalue,
    ) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}
