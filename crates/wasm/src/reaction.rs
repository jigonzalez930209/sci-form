//! Reaction-focused WASM bindings — SMIRKS transforms, multi-method NEB,
//! single-point energy queries, and configurable MD for reaction dynamics.

use crate::helpers::*;
use wasm_bindgen::prelude::*;

// ─── SMIRKS ─────────────────────────────────────────────────────────────────

/// Parse a SMIRKS reaction transform string.
///
/// Returns JSON `SmirksTransform`:
/// ```json
/// {
///   "reactant_smarts": ["[C:1](=O)[OH:2]"],
///   "product_smarts":  ["[C:1](=O)[O-:2]"],
///   "atom_map": {1: 1, 2: 2},
///   "bond_changes": [...],
///   "smirks": "[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]"
/// }
/// ```
#[wasm_bindgen]
pub fn parse_smirks(smirks: &str) -> String {
    match sci_form::smirks::parse_smirks(smirks) {
        Ok(transform) => serialize_or_error(&transform),
        Err(e) => json_error(&e),
    }
}

/// Apply a SMIRKS reaction transform to a molecule.
///
/// Returns JSON `SmirksResult`:
/// ```json
/// {
///   "products": ["[C](=O)[O-]"],
///   "atom_mapping": {0: 0, 1: 1},
///   "n_transforms": 1,
///   "success": true,
///   "messages": []
/// }
/// ```
#[wasm_bindgen]
pub fn apply_smirks(smirks: &str, smiles: &str) -> String {
    match sci_form::smirks::apply_smirks(smirks, smiles) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

// ─── Multi-method NEB ───────────────────────────────────────────────────────

/// Build a simplified NEB path with a configurable energy backend.
///
/// `method`: `"uff"`, `"mmff94"`, `"pm3"`, `"xtb"`, `"gfn1"`, `"gfn2"`, `"hf3c"`.
///
/// Returns JSON `NebPathResult` with images, energies, TS info.
#[wasm_bindgen]
pub fn compute_neb_path_with_method(
    smiles: &str,
    start_coords_json: &str,
    end_coords_json: &str,
    n_images: usize,
    n_iter: usize,
    spring_k: f64,
    step_size: f64,
    method: &str,
) -> String {
    let start = match parse_flat_coords(start_coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let end = match parse_flat_coords(end_coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_simplified_neb_path_configurable(
        smiles, &start, &end, n_images, n_iter, spring_k, step_size, method,
    ) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

// ─── Single-point energy query ──────────────────────────────────────────────

/// Evaluate the energy (kcal/mol) of a geometry with a given NEB backend.
///
/// `method`: `"uff"`, `"mmff94"`, `"pm3"`, `"xtb"`, `"gfn1"`, `"gfn2"`, `"hf3c"`.
///
/// Returns JSON `{"energy": <f64>, "method": "<str>"}`.
#[wasm_bindgen]
pub fn neb_single_point_energy(method: &str, smiles: &str, coords_json: &str) -> String {
    let coords = match parse_flat_coords(coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::neb_backend_energy_kcal(method, smiles, &coords) {
        Ok(energy) => serde_json::json!({
            "energy": energy,
            "method": method,
        })
        .to_string(),
        Err(e) => json_error(&e),
    }
}

/// Evaluate the energy (kcal/mol) and gradient of a geometry with a given NEB backend.
///
/// `method`: `"uff"`, `"mmff94"`, `"pm3"`, `"xtb"`, `"gfn1"`, `"gfn2"`, `"hf3c"`.
///
/// Returns JSON `{"energy": <f64>, "gradient": [<f64>...], "method": "<str>"}`.
#[wasm_bindgen]
pub fn neb_energy_and_gradient(method: &str, smiles: &str, coords_json: &str) -> String {
    let coords = match parse_flat_coords(coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::neb_backend_energy_and_gradient(method, smiles, &coords) {
        Ok((energy, gradient)) => serde_json::json!({
            "energy": energy,
            "gradient": gradient,
            "method": method,
        })
        .to_string(),
        Err(e) => json_error(&e),
    }
}

// ─── Multi-method MD ────────────────────────────────────────────────────────

/// Run velocity-Verlet MD with a configurable energy backend.
///
/// `backend`: `"uff"`, `"pm3"`, `"xtb"`.
/// `target_temp_k`: if > 0, enables Berendsen thermostat with 100 fs coupling.
///
/// Returns JSON `MdTrajectory`.
#[wasm_bindgen]
pub fn compute_md_trajectory_with_method(
    smiles: &str,
    coords_json: &str,
    elements_json: &str,
    n_steps: usize,
    dt_fs: f64,
    seed: u32,
    target_temp_k: f64,
    backend: &str,
) -> String {
    let flat = match parse_flat_coords(coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let elements = match parse_elements(elements_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };

    let md_backend = match backend.to_lowercase().as_str() {
        "uff" => sci_form::dynamics::MdBackend::Uff,
        "pm3" => sci_form::dynamics::MdBackend::Pm3,
        "xtb" => sci_form::dynamics::MdBackend::Xtb,
        _ => return json_error(&format!("unsupported MD backend: {}", backend)),
    };

    let thermostat = if target_temp_k > 0.0 {
        Some((target_temp_k, 100.0)) // (target_temp, tau_fs)
    } else {
        None
    };

    match sci_form::dynamics::simulate_velocity_verlet(
        smiles,
        &flat,
        &elements,
        n_steps,
        dt_fs,
        seed as u64,
        thermostat,
        md_backend,
    ) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// List available NEB backend methods for a given set of elements.
///
/// Returns JSON array of `{"method": string, "available": bool, "note": string}`.
#[wasm_bindgen]
pub fn list_neb_backends(elements_json: &str) -> String {
    let elements = match parse_elements(elements_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };

    let has_transition_metals = elements
        .iter()
        .any(|&z| matches!(z, 22 | 24..=30 | 44 | 46 | 47 | 78 | 79));

    let backends = vec![
        serde_json::json!({
            "method": "uff",
            "available": true,
            "note": "Universal Force Field — fastest, broadest element coverage"
        }),
        serde_json::json!({
            "method": "mmff94",
            "available": !has_transition_metals,
            "note": "MMFF94 — organic molecules only"
        }),
        serde_json::json!({
            "method": "pm3",
            "available": true,
            "note": "PM3 semi-empirical SCF"
        }),
        serde_json::json!({
            "method": "xtb",
            "available": true,
            "note": "GFN0-xTB tight-binding"
        }),
        serde_json::json!({
            "method": "gfn1",
            "available": true,
            "note": "GFN1-xTB with D3 dispersion"
        }),
        serde_json::json!({
            "method": "gfn2",
            "available": true,
            "note": "GFN2-xTB with D4 dispersion — most accurate"
        }),
        serde_json::json!({
            "method": "hf3c",
            "available": !has_transition_metals,
            "note": "HF-3c minimal-basis Hartree-Fock — slowest, most accurate"
        }),
    ];

    serde_json::to_string(&backends).unwrap_or_else(|e| json_error(&e.to_string()))
}

// ─── Full reaction dynamics ─────────────────────────────────────────────────

/// Compute a full reaction dynamics path: embed reactant + product fragments,
/// build oriented complexes, run NEB for the reactive zone, and generate
/// approach/departure frames — all energies computed in Rust.
///
/// `reactant_smiles_json`: JSON array of reactant SMILES, e.g. `["CC(=O)O", "N"]`.
/// `product_smiles_json`:  JSON array of product SMILES, e.g. `["CC(=O)N", "O"]`.
/// `method`: energy backend — `"uff"`, `"mmff94"`, `"pm3"`, `"xtb"`, `"gfn1"`, `"gfn2"`, `"hf3c"`.
/// `config_json`: optional JSON [`ReactionDynamicsConfig`]; pass `""` or `"{}"` for defaults.
///
/// Returns JSON [`ReactionDynamicsResult`] with frames, energies, TS info.
#[wasm_bindgen]
pub fn compute_reaction_dynamics(
    reactant_smiles_json: &str,
    product_smiles_json: &str,
    method: &str,
    config_json: &str,
) -> String {
    let reactants: Vec<String> = match serde_json::from_str(reactant_smiles_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad reactant_smiles: {}", e)),
    };
    let products: Vec<String> = match serde_json::from_str(product_smiles_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad product_smiles: {}", e)),
    };

    let config: sci_form::dynamics::ReactionDynamicsConfig = if config_json.is_empty()
        || config_json == "{}"
    {
        sci_form::dynamics::ReactionDynamicsConfig::default()
    } else {
        match serde_json::from_str(config_json) {
            Ok(c) => c,
            Err(e) => return json_error(&format!("bad config: {}", e)),
        }
    };

    let r_refs: Vec<&str> = reactants.iter().map(String::as_str).collect();
    let p_refs: Vec<&str> = products.iter().map(String::as_str).collect();

    match sci_form::compute_reaction_dynamics(&r_refs, &p_refs, method, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}
