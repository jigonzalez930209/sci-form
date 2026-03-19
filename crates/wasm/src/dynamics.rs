//! Molecular dynamics and pathway WASM bindings.

use crate::helpers::*;
use wasm_bindgen::prelude::*;

/// Run short NVE molecular dynamics with UFF force field.
///
/// `smiles`: SMILES string for topology.
/// `coords_json`: JSON array of flat xyz coords.
/// `n_steps`: number of MD steps.
/// `dt_fs`: time step in femtoseconds.
/// `seed`: random seed for initial velocities.
///
/// Returns JSON MdTrajectory with frames, energies, temperatures.
#[wasm_bindgen]
pub fn compute_md_trajectory(
    smiles: &str,
    coords_json: &str,
    n_steps: usize,
    dt_fs: f64,
    seed: u32,
) -> String {
    let flat = match parse_flat_coords(coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_md_trajectory(smiles, &flat, n_steps, dt_fs, seed as u64) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Run short NVT molecular dynamics with Berendsen thermostat.
///
/// `target_temp_k`: target temperature in Kelvin.
/// `thermostat_tau_fs`: thermostat coupling time in fs.
#[wasm_bindgen]
pub fn compute_md_trajectory_nvt(
    smiles: &str,
    coords_json: &str,
    n_steps: usize,
    dt_fs: f64,
    seed: u32,
    target_temp_k: f64,
    thermostat_tau_fs: f64,
) -> String {
    let flat = match parse_flat_coords(coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_md_trajectory_nvt(
        smiles,
        &flat,
        n_steps,
        dt_fs,
        seed as u64,
        target_temp_k,
        thermostat_tau_fs,
    ) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Build a simplified NEB path between two geometries.
///
/// `start_coords_json`, `end_coords_json`: JSON flat xyz coords.
/// `n_images`: number of intermediate images.
/// `n_iter`: optimization iterations.
/// `spring_k`: spring constant.
/// `step_size`: optimization step size.
#[wasm_bindgen]
pub fn compute_simplified_neb_path(
    smiles: &str,
    start_coords_json: &str,
    end_coords_json: &str,
    n_images: usize,
    n_iter: usize,
    spring_k: f64,
    step_size: f64,
) -> String {
    let start = match parse_flat_coords(start_coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let end = match parse_flat_coords(end_coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_simplified_neb_path(
        smiles, &start, &end, n_images, n_iter, spring_k, step_size,
    ) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Search conformers by sampling, torsion refinement, and UFF ranking.
///
/// `n_samples`: number of initial conformer samples.
/// `rmsd_threshold`: RMSD cutoff for duplicate filtering (Å).
#[wasm_bindgen]
pub fn search_conformers_with_uff(
    smiles: &str,
    n_samples: usize,
    seed: u32,
    rmsd_threshold: f64,
) -> String {
    match sci_form::search_conformers_with_uff(smiles, n_samples, seed as u64, rmsd_threshold) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}
