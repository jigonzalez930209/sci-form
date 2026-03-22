//! Orbital mesh generation WASM bindings for all methods.

use crate::helpers::*;
use wasm_bindgen::prelude::*;

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
use crate::eht::compute_eht_orbital_grid_wasm_accelerated;

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
use crate::webgpu::{self, WasmExecutionMode};

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

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
#[wasm_bindgen]
pub async fn compute_orbital_mesh_accelerated(
    elements_json: &str,
    coords_flat_json: &str,
    method: &str,
    mo_index: usize,
    spacing: f64,
    padding: f64,
    isovalue: f32,
    execution_mode: &str,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let mode = match webgpu::parse_execution_mode(execution_mode) {
        Ok(mode) => mode,
        Err(err) => return json_error(&err),
    };

    match compute_orbital_mesh_wasm_accelerated(
        &elems, &positions, method, mo_index, spacing, padding, isovalue, mode,
    )
    .await
    {
        Ok((result, backend, used_gpu, note)) => serde_json::json!({
            "mesh": result.mesh,
            "grid": result.grid,
            "method": result.method,
            "mo_index": result.mo_index,
            "homo_index": result.homo_index,
            "orbital_energies": result.orbital_energies,
            "gap": result.gap,
            "backend": backend,
            "used_gpu": used_gpu,
            "mode": execution_mode,
            "note": note,
        })
        .to_string(),
        Err(err) => json_error(&err),
    }
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
async fn compute_orbital_mesh_wasm_accelerated(
    elements: &[u8],
    positions: &[[f64; 3]],
    method: &str,
    mo_index: usize,
    spacing: f64,
    padding: f64,
    isovalue: f32,
    mode: WasmExecutionMode,
) -> Result<(sci_form::mesh::OrbitalMeshResult, String, bool, String), String> {
    let (mesh_method, homo_index, orbital_energies, gap) =
        mesh_method_metadata(elements, positions, method)?;
    let (grid, backend, used_gpu, grid_note) = compute_eht_orbital_grid_wasm_accelerated(
        elements, positions, mo_index, spacing, padding, mode,
    )
    .await?;
    let mesh = sci_form::eht::marching_cubes(&grid, isovalue);

    Ok((
        sci_form::mesh::OrbitalMeshResult {
            mesh,
            grid,
            method: mesh_method,
            mo_index,
            homo_index,
            orbital_energies,
            gap,
        },
        backend,
        used_gpu,
        format!("{grid_note}; marching cubes remains CPU-side after the accelerated grid stage"),
    ))
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
fn mesh_method_metadata(
    elements: &[u8],
    positions: &[[f64; 3]],
    method: &str,
) -> Result<(sci_form::mesh::MeshMethod, usize, Vec<f64>, f64), String> {
    match method.trim().to_ascii_lowercase().as_str() {
        "eht" => {
            let result = sci_form::eht::solve_eht(elements, positions, None)?;
            Ok((
                sci_form::mesh::MeshMethod::Eht,
                result.homo_index,
                result.energies,
                result.gap,
            ))
        }
        "pm3" => {
            let result = sci_form::compute_pm3(elements, positions)?;
            let homo_index = if result.n_electrons > 0 {
                result.n_electrons / 2 - 1
            } else {
                0
            };
            Ok((
                sci_form::mesh::MeshMethod::Pm3,
                homo_index,
                result.orbital_energies,
                result.gap,
            ))
        }
        "xtb" => {
            let result = sci_form::compute_xtb(elements, positions)?;
            let homo_index = if result.n_electrons > 0 {
                result.n_electrons / 2 - 1
            } else {
                0
            };
            Ok((
                sci_form::mesh::MeshMethod::Xtb,
                homo_index,
                result.orbital_energies,
                result.gap,
            ))
        }
        "hf3c" => {
            let config = sci_form::hf::HfConfig::default();
            let result = sci_form::hf::solve_hf3c(elements, positions, &config)?;
            let n_electrons: usize = elements.iter().map(|&z| z as usize).sum();
            let homo_index = if n_electrons > 0 {
                n_electrons / 2 - 1
            } else {
                0
            };
            let gap = if result.orbital_energies.len() > homo_index + 1 {
                result.orbital_energies[homo_index + 1] - result.orbital_energies[homo_index]
            } else {
                0.0
            };
            Ok((
                sci_form::mesh::MeshMethod::Hf3c,
                homo_index,
                result.orbital_energies,
                gap,
            ))
        }
        other => Err(format!("unsupported mesh method: {other}")),
    }
}
