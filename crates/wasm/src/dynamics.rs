//! Molecular dynamics and pathway WASM bindings.

use crate::helpers::*;
use wasm_bindgen::prelude::*;

/// Run short NVE molecular dynamics with UFF force field.
///
/// `smiles`: SMILES string for topology.
/// `coords_json`: JSON array of flat xyz coords.
/// `n_steps`: number of MD steps.
/// `dt_fs`: time step in femtoseconds.
/// `seed`: random seed for initial velocities (32-bit; use
/// [`compute_md_trajectory_seed64`] for the full 64-bit range).
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
    compute_md_trajectory_seed64(smiles, coords_json, n_steps, dt_fs, seed as u64)
}

/// Run short NVE molecular dynamics with UFF force field, 64-bit seed.
///
/// `seed` is a JavaScript `BigInt` (e.g. `42n`).
#[wasm_bindgen]
pub fn compute_md_trajectory_seed64(
    smiles: &str,
    coords_json: &str,
    n_steps: usize,
    dt_fs: f64,
    seed: u64,
) -> String {
    let flat = match parse_flat_coords(coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_md_trajectory(smiles, &flat, n_steps, dt_fs, seed) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Run short NVT molecular dynamics with Berendsen thermostat.
///
/// `target_temp_k`: target temperature in Kelvin.
/// `thermostat_tau_fs`: thermostat coupling time in fs.
/// `seed`: 32-bit random seed (use [`compute_md_trajectory_nvt_seed64`] for the
/// full 64-bit range).
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
    compute_md_trajectory_nvt_seed64(
        smiles,
        coords_json,
        n_steps,
        dt_fs,
        seed as u64,
        target_temp_k,
        thermostat_tau_fs,
    )
}

/// Run short NVT molecular dynamics with Berendsen thermostat, 64-bit seed.
///
/// `seed` is a JavaScript `BigInt` (e.g. `42n`).
#[wasm_bindgen]
pub fn compute_md_trajectory_nvt_seed64(
    smiles: &str,
    coords_json: &str,
    n_steps: usize,
    dt_fs: f64,
    seed: u64,
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
        seed,
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
/// `seed`: 32-bit random seed (use [`search_conformers_with_uff_seed64`] for
/// the full 64-bit range).
#[wasm_bindgen]
pub fn search_conformers_with_uff(
    smiles: &str,
    n_samples: usize,
    seed: u32,
    rmsd_threshold: f64,
) -> String {
    search_conformers_with_uff_seed64(smiles, n_samples, seed as u64, rmsd_threshold)
}

/// Search conformers by sampling, torsion refinement, and UFF ranking, 64-bit seed.
///
/// `seed` is a JavaScript `BigInt` (e.g. `42n`).
#[wasm_bindgen]
pub fn search_conformers_with_uff_seed64(
    smiles: &str,
    n_samples: usize,
    seed: u64,
    rmsd_threshold: f64,
) -> String {
    match sci_form::search_conformers_with_uff(smiles, n_samples, seed, rmsd_threshold) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}
