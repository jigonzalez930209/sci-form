//! **ALPHA** — WASM bindings for all alpha-tier experimental modules.
//!
//! Import separately from the stable API:
//!
//! ```js
//! import { compute_dft, compute_reaxff_gradient, compute_mlff }
//!   from 'sci-form-wasm/alpha';
//! ```
//!
//! Each function is individually gated; only the ones compiled with the
//! matching feature flag are available in a given build.

#[allow(unused_imports)]
use crate::helpers::{json_error, parse_elements_and_positions, parse_flat_coords};
#[allow(unused_imports)]
use wasm_bindgen::prelude::*;

// ─── A1: Kohn-Sham DFT ──────────────────────────────────────────────────────

/// Run a Kohn-Sham DFT single-point calculation.
///
/// `method`: `"svwn"` (LDA) | `"pbe"` (GGA, default).
///
/// Returns JSON `DftResult`:
/// ```json
/// { "energy": -1.12, "homo_energy": -0.5, "lumo_energy": 0.1, "gap": 0.6,
///   "converged": true, "n_basis": 5, "scf_iterations": 12,
///   "mulliken_charges": [...], "orbital_energies": [...],
///   "nuclear_repulsion": 0.72, "xc_energy": -0.3 }
/// ```
#[cfg(feature = "alpha-dft")]
#[wasm_bindgen]
pub fn alpha_compute_dft(elements_json: &str, coords_flat_json: &str, method: &str) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let dft_method = match method.to_lowercase().as_str() {
        "svwn" | "lda" => sci_form::dft::ks_fock::DftMethod::Svwn,
        _ => sci_form::dft::ks_fock::DftMethod::Pbe,
    };
    let config = sci_form::dft::ks_fock::DftConfig {
        method: dft_method,
        ..Default::default()
    };
    match sci_form::dft::ks_fock::solve_ks_dft(&elems, &positions, &config) {
        Ok(r) => serde_json::to_string(&r).unwrap_or_else(|e| json_error(&e.to_string())),
        Err(e) => json_error(&e),
    }
}

// ─── A2: ReaxFF ─────────────────────────────────────────────────────────────

/// Compute ReaxFF energy and gradient for a molecule.
///
/// Returns JSON:
/// ```json
/// { "energy_kcal_mol": 12.3,
///   "gradient": [gx0, gy0, gz0, gx1, gy1, gz1, ...] }
/// ```
#[cfg(feature = "alpha-reaxff")]
#[wasm_bindgen]
pub fn alpha_compute_reaxff_gradient(elements_json: &str, coords_flat_json: &str) -> String {
    let flat: Vec<f64> = match parse_flat_coords(coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let elems: Vec<u8> = match serde_json::from_str(elements_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad elements: {}", e)),
    };
    let params = sci_form::forcefield::reaxff::params::ReaxffParams::default_chon();
    match sci_form::forcefield::reaxff::gradients::compute_reaxff_gradient(&flat, &elems, &params) {
        Ok((energy, gradient)) => serde_json::json!({
            "energy_kcal_mol": energy,
            "gradient": gradient
        })
        .to_string(),
        Err(e) => json_error(&e),
    }
}

/// Compute ReaxFF energy components for a molecule (no gradient).
///
/// Returns JSON with keys `bonded`, `coulomb`, `van_der_waals`, `total`.
#[cfg(feature = "alpha-reaxff")]
#[wasm_bindgen]
pub fn alpha_compute_reaxff_energy(elements_json: &str, coords_flat_json: &str) -> String {
    let flat: Vec<f64> = match parse_flat_coords(coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let elems: Vec<u8> = match serde_json::from_str(elements_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad elements: {}", e)),
    };
    let positions: Vec<[f64; 3]> = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    let params = sci_form::forcefield::reaxff::params::ReaxffParams::default_chon();
    let charges = sci_form::forcefield::reaxff::eem::solve_eem(
        &elems,
        &positions,
        &sci_form::forcefield::reaxff::eem::default_eem_params(),
    );
    let bo =
        sci_form::forcefield::reaxff::bond_order::compute_bond_orders(&elems, &positions, &params);
    let bonded = sci_form::forcefield::reaxff::energy::compute_bonded_energy(&bo, &elems, &params);
    let nonbonded = sci_form::forcefield::reaxff::nonbonded::compute_nonbonded_energy(
        &elems, &positions, &charges, &params,
    );
    serde_json::json!({
        "bonded": bonded,
        "coulomb": nonbonded.coulomb,
        "van_der_waals": nonbonded.van_der_waals,
        "total": bonded + nonbonded.coulomb + nonbonded.van_der_waals
    })
    .to_string()
}

/// Solve EEM charges for a molecule using the ReaxFF charge model.
///
/// Returns JSON: `{ "charges": [...], "total_charge": 0.0 }`
#[cfg(feature = "alpha-reaxff")]
#[wasm_bindgen]
pub fn alpha_compute_eem_charges(elements_json: &str, coords_flat_json: &str) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let params = sci_form::forcefield::reaxff::eem::default_eem_params();
    let charges = sci_form::forcefield::reaxff::eem::solve_eem(&elems, &positions, &params);
    let total: f64 = charges.iter().sum();
    serde_json::json!({
        "charges": charges,
        "total_charge": total
    })
    .to_string()
}

// ─── A3: MLFF Neural Network Force Field ────────────────────────────────────

/// Compute MLFF energy and forces using pre-configured element networks.
///
/// `config_json` is a JSON object matching `MlffConfig`:
/// ```json
/// {
///   "aev_params": { "radial_cutoff": 6.0, "angular_cutoff": 3.5, ... },
///   "element_nets": {
///     "1":  { "layers": [[...weights...]], "biases": [[...]] },
///     "6":  { ... }
///   }
/// }
/// ```
///
/// Returns JSON: `{ "energy": -1.23, "atomic_energies": [...], "forces": [...] }`
#[cfg(feature = "alpha-mlff")]
#[wasm_bindgen]
pub fn alpha_compute_mlff(
    elements_json: &str,
    coords_flat_json: &str,
    config_json: &str,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let config: sci_form::mlff::MlffConfig = match serde_json::from_str(config_json) {
        Ok(c) => c,
        Err(e) => return json_error(&format!("bad mlff config: {}", e)),
    };
    match sci_form::mlff::compute_mlff(&elems, &positions, &config) {
        Ok(r) => serde_json::to_string(&r).unwrap_or_else(|e| json_error(&e.to_string())),
        Err(e) => json_error(&e),
    }
}

/// Compute Atomic Environment Vectors (AEVs) for each atom.
///
/// Uses default symmetry function parameters unless `params_json` is provided.
///
/// Returns JSON: `{ "n_atoms": 5, "aev_length": 384, "aevs": [[...], [...]] }`
#[cfg(feature = "alpha-mlff")]
#[wasm_bindgen]
pub fn alpha_compute_aevs(
    elements_json: &str,
    coords_flat_json: &str,
    params_json: &str,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let params: sci_form::mlff::SymmetryFunctionParams =
        if params_json.trim().is_empty() || params_json == "null" || params_json == "{}" {
            sci_form::mlff::SymmetryFunctionParams::default()
        } else {
            match serde_json::from_str(params_json) {
                Ok(p) => p,
                Err(e) => return json_error(&format!("bad params: {}", e)),
            }
        };
    let aevs = sci_form::mlff::compute_aevs(&elems, &positions, &params);
    let aev_length = aevs.first().map(|v| v.len()).unwrap_or(0);
    let aev_arrays: Vec<Vec<f64>> = aevs.iter().map(|v| v.as_slice().to_vec()).collect();
    serde_json::json!({
        "n_atoms": elems.len(),
        "aev_length": aev_length,
        "aevs": aev_arrays
    })
    .to_string()
}

// ─── A4: Obara-Saika ERIs ───────────────────────────────────────────────────

/// Compute the Boys function F_n(x).
///
/// Returns JSON: `{ "value": 0.9846 }`
#[cfg(feature = "alpha-obara-saika")]
#[wasm_bindgen]
pub fn alpha_boys_function(n: u32, x: f64) -> String {
    let v = sci_form::scf::obara_saika::boys_function(n as usize, x);
    serde_json::json!({ "value": v }).to_string()
}

/// Compute a single (ss|ss) two-electron repulsion integral using Obara-Saika.
///
/// `shell_a_json` / `shell_b_json` / `shell_c_json` / `shell_d_json` are JSON
/// objects `{ "alpha": 1.0, "center": [x, y, z] }`.
///
/// Returns JSON: `{ "eri": 4.37 }`
#[cfg(feature = "alpha-obara-saika")]
#[wasm_bindgen]
pub fn alpha_eri_ssss(
    shell_a_json: &str,
    shell_b_json: &str,
    shell_c_json: &str,
    shell_d_json: &str,
) -> String {
    #[derive(serde::Deserialize)]
    struct Shell {
        alpha: f64,
        center: [f64; 3],
    }
    let parse_shell = |s: &str| -> Result<Shell, String> {
        serde_json::from_str(s).map_err(|e| format!("bad shell: {}", e))
    };
    let (a, b, c, d) = match (
        parse_shell(shell_a_json),
        parse_shell(shell_b_json),
        parse_shell(shell_c_json),
        parse_shell(shell_d_json),
    ) {
        (Ok(a), Ok(b), Ok(c), Ok(d)) => (a, b, c, d),
        (Err(e), _, _, _) | (_, Err(e), _, _) | (_, _, Err(e), _) | (_, _, _, Err(e)) => {
            return json_error(&e)
        }
    };
    let sp_ab =
        sci_form::scf::obara_saika::ShellPairData::new(a.alpha, a.center, b.alpha, b.center);
    let sp_cd =
        sci_form::scf::obara_saika::ShellPairData::new(c.alpha, c.center, d.alpha, d.center);
    let eri = sci_form::scf::obara_saika::eri_ssss(&sp_ab, &sp_cd);
    serde_json::json!({ "eri": eri }).to_string()
}

/// Compute Schwarz pre-screening bound Q_ab = √〈ab|ab〉.
///
/// Returns JSON: `{ "schwarz_bound": 0.8432 }`.
#[cfg(feature = "alpha-obara-saika")]
#[wasm_bindgen]
pub fn alpha_schwarz_bound(
    alpha: f64,
    center_a_json: &str,
    beta: f64,
    center_b_json: &str,
) -> String {
    let ca: [f64; 3] = match serde_json::from_str(center_a_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad center_a: {}", e)),
    };
    let cb: [f64; 3] = match serde_json::from_str(center_b_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad center_b: {}", e)),
    };
    let v = sci_form::scf::obara_saika::schwarz_bound(alpha, ca, beta, cb);
    serde_json::json!({ "schwarz_bound": v }).to_string()
}

// ─── A5: CGA Motors ─────────────────────────────────────────────────────────

/// Apply a CGA dihedral rotation to a sub-tree of atoms.
///
/// `coords_flat_json` — flat [x0,y0,z0,...] for all atoms.
/// `subtree_indices_json` — atom indices to rotate.
/// `axis_a_json` / `axis_b_json` — the two atoms defining the rotation axis ([x,y,z]).
/// `angle_rad` — rotation angle in radians.
///
/// Returns JSON: `{ "coords": [x0,y0,z0,...] }`
#[cfg(feature = "alpha-cga")]
#[wasm_bindgen]
pub fn alpha_rotate_dihedral_cga(
    coords_flat_json: &str,
    subtree_indices_json: &str,
    axis_a_json: &str,
    axis_b_json: &str,
    angle_rad: f64,
) -> String {
    let flat: Vec<f64> = match parse_flat_coords(coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let indices: Vec<usize> = match serde_json::from_str(subtree_indices_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad indices: {}", e)),
    };
    let a: [f64; 3] = match serde_json::from_str(axis_a_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad axis_a: {}", e)),
    };
    let b: [f64; 3] = match serde_json::from_str(axis_b_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad axis_b: {}", e)),
    };

    let motor = sci_form::alpha::cga::dihedral_motor(a, b, angle_rad);
    let mut coords = flat.clone();
    sci_form::alpha::cga::apply_motor_to_subtree(&mut coords, &indices, &motor);
    serde_json::json!({ "coords": coords }).to_string()
}

/// Refine a single dihedral torsion to a target angle using CGA.
///
/// `smiles` — SMILES string for bond connectivity.
/// `torsion_indices_json` — 4 atom indices [i, j, k, l].
/// `target_angle_rad` — target dihedral angle in radians.
///
/// Returns JSON: `{ "coords": [...] }`
#[cfg(feature = "alpha-cga")]
#[wasm_bindgen]
pub fn alpha_refine_torsion_cga(
    elements_json: &str,
    coords_flat_json: &str,
    smiles: &str,
    torsion_indices_json: &str,
    target_angle_rad: f64,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let torsion: Vec<usize> = match serde_json::from_str(torsion_indices_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad torsion_indices: {}", e)),
    };
    if torsion.len() != 4 {
        return json_error("torsion_indices must have exactly 4 elements");
    }
    // Parse bonds from SMILES
    let mol = match sci_form::parse(smiles) {
        Ok(m) => m,
        Err(e) => return json_error(&format!("parse error: {}", e)),
    };
    let bonds: Vec<(usize, usize)> = mol
        .graph
        .edge_indices()
        .map(|e| {
            let (a, b) = mol.graph.edge_endpoints(e).unwrap();
            (a.index(), b.index())
        })
        .collect();
    let flat: Vec<f64> = positions.iter().flat_map(|p| p.iter().cloned()).collect();
    let new_coords = sci_form::alpha::cga::refine_torsion_cga(
        &elems,
        &flat,
        &bonds,
        &[torsion[0], torsion[1], torsion[2], torsion[3]],
        target_angle_rad,
    );
    serde_json::json!({ "coords": new_coords }).to_string()
}

// ─── A6: Growing String Method ──────────────────────────────────────────────

/// Interpolate a node along the reaction path at parameter t ∈ [0, 1].
///
/// `reactant_json` / `product_json` — flat coordinate arrays.
///
/// Returns JSON: `{ "coords": [...] }`
#[cfg(feature = "alpha-gsm")]
#[wasm_bindgen]
pub fn alpha_gsm_interpolate(reactant_json: &str, product_json: &str, t: f64) -> String {
    let reactant: Vec<f64> = match serde_json::from_str(reactant_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad reactant: {}", e)),
    };
    let product: Vec<f64> = match serde_json::from_str(product_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad product: {}", e)),
    };
    let node = sci_form::alpha::gsm::interpolate_node(&reactant, &product, t);
    serde_json::json!({ "coords": node }).to_string()
}

/// Grow the reaction string and find the transition state node between
/// reactant and product geometries using UFF energy.
///
/// Returns JSON `GsmResult`:
/// ```json
/// { "ts_coords": [...], "ts_energy": 12.3,
///   "path_energies": [...], "n_nodes": 5, "converged": true }
/// ```
#[cfg(feature = "alpha-gsm")]
#[wasm_bindgen]
pub fn alpha_gsm_find_ts(
    smiles: &str,
    reactant_coords_json: &str,
    product_coords_json: &str,
    config_json: &str,
) -> String {
    let reactant: Vec<f64> = match serde_json::from_str(reactant_coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad reactant: {}", e)),
    };
    let product: Vec<f64> = match serde_json::from_str(product_coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad product: {}", e)),
    };
    let config: sci_form::alpha::gsm::GsmConfig =
        if config_json.trim().is_empty() || config_json == "null" || config_json == "{}" {
            sci_form::alpha::gsm::GsmConfig::default()
        } else {
            match serde_json::from_str(config_json) {
                Ok(c) => c,
                Err(e) => return json_error(&format!("bad gsm config: {}", e)),
            }
        };

    // Use UFF as energy evaluator (only topology needed from SMILES)
    let energy_fn = |coords: &[f64]| -> f64 {
        sci_form::compute_uff_energy(smiles, coords).unwrap_or(f64::MAX)
    };

    match sci_form::alpha::gsm::find_transition_state(&reactant, &product, &config, &energy_fn) {
        Ok(r) => serde_json::to_string(&r).unwrap_or_else(|e| json_error(&e.to_string())),
        Err(e) => json_error(&e),
    }
}

// ─── A7: SDR Embedding ──────────────────────────────────────────────────────

/// Embed a molecule from pairwise distance constraints using Semidefinite Relaxation.
///
/// `distance_pairs_json` — array of `[[i, j, d_ij], ...]` where d_ij is in Å.
/// `n_atoms` — total number of atoms.
/// `config_json` — optional SdrConfig JSON or `{}` for defaults.
///
/// Returns JSON `SdrResult`:
/// ```json
/// { "coords": [...], "residual": 0.001, "converged": true, "iterations": 45 }
/// ```
#[cfg(feature = "alpha-sdr")]
#[wasm_bindgen]
pub fn alpha_sdr_embed(distance_pairs_json: &str, n_atoms: usize, config_json: &str) -> String {
    // distance_pairs: [[i, j, d], ...]
    let raw_pairs: Vec<[f64; 3]> = match serde_json::from_str(distance_pairs_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad distance_pairs: {}", e)),
    };
    let distance_pairs: Vec<(usize, usize, f64)> = raw_pairs
        .iter()
        .map(|&[i, j, d]| (i as usize, j as usize, d))
        .collect();

    let config: sci_form::alpha::sdr::SdrConfig =
        if config_json.trim().is_empty() || config_json == "null" || config_json == "{}" {
            sci_form::alpha::sdr::SdrConfig::default()
        } else {
            match serde_json::from_str(config_json) {
                Ok(c) => c,
                Err(e) => return json_error(&format!("bad sdr config: {}", e)),
            }
        };

    let result = sci_form::alpha::sdr::sdr_embed(&distance_pairs, n_atoms, &config);
    serde_json::json!({
        "coords": result.coords,
        "residual": result.residual,
        "converged": result.converged,
        "iterations": result.iterations
    })
    .to_string()
}

// ─── A8: Live MD Simulation (re-exported from dynamics_live module) ──────────
//
// The full `LiveSimulation` JS class is exposed in `crates/wasm/src/dynamics_live.rs`.
// This section adds metadata helpers for the alpha subpath export.

/// Return alpha module metadata.
///
/// Returns JSON: `{ "alpha_modules": ["dft", "reaxff", "mlff", "obara-saika",
///                   "cga", "gsm", "sdr", "dynamics-live"] }`
#[wasm_bindgen]
pub fn alpha_modules_info() -> String {
    serde_json::json!({
        "alpha_modules": [
            "dft",
            "reaxff",
            "mlff",
            "obara-saika",
            "cga",
            "gsm",
            "sdr",
            "dynamics-live"
        ],
        "stability": "alpha — experimental, subject to breaking changes",
        "version": "0.11.2"
    })
    .to_string()
}
