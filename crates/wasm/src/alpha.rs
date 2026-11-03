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
use crate::helpers::{
    bad_json, json_error, parse_elements_and_positions, parse_flat_coords, smiles_error,
};
#[allow(unused_imports)]
use wasm_bindgen::prelude::*;

#[cfg(feature = "alpha-reaxff")]
fn build_reaxff_atom_params(
    elements: &[u8],
) -> Vec<sci_form::forcefield::reaxff::params::ReaxffAtomParams> {
    let params = sci_form::forcefield::reaxff::params::ReaxffParams::default_chon();
    let fallback = params.atom_params.last().cloned().unwrap_or(
        sci_form::forcefield::reaxff::params::ReaxffAtomParams {
            element: 1,
            r_sigma: 1.0,
            r_pi: 0.0,
            r_pipi: 0.0,
            p_bo1: 0.0,
            p_bo2: 1.0,
            p_bo3: 0.0,
            p_bo4: 1.0,
            p_bo5: 0.0,
            p_bo6: 1.0,
            valence: 1.0,
        },
    );

    elements
        .iter()
        .map(|&z| {
            params
                .element_index(z)
                .map(|idx| params.atom_params[idx].clone())
                .unwrap_or_else(|| fallback.clone())
        })
        .collect()
}

fn dihedral_angle(p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], p4: [f64; 3]) -> f64 {
    let b1 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
    let b2 = [p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]];
    let b3 = [p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2]];

    let n1 = [
        b1[1] * b2[2] - b1[2] * b2[1],
        b1[2] * b2[0] - b1[0] * b2[2],
        b1[0] * b2[1] - b1[1] * b2[0],
    ];
    let n2 = [
        b2[1] * b3[2] - b2[2] * b3[1],
        b2[2] * b3[0] - b2[0] * b3[2],
        b2[0] * b3[1] - b2[1] * b3[0],
    ];
    let b2_norm = (b2[0] * b2[0] + b2[1] * b2[1] + b2[2] * b2[2]).sqrt();
    let b2_unit = if b2_norm > 1e-15 {
        [b2[0] / b2_norm, b2[1] / b2_norm, b2[2] / b2_norm]
    } else {
        [0.0, 0.0, 0.0]
    };
    let m1 = [
        n1[1] * b2_unit[2] - n1[2] * b2_unit[1],
        n1[2] * b2_unit[0] - n1[0] * b2_unit[2],
        n1[0] * b2_unit[1] - n1[1] * b2_unit[0],
    ];

    let x = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2];
    let y = m1[0] * n2[0] + m1[1] * n2[1] + m1[2] * n2[2];
    (-y).atan2(-x) + std::f64::consts::PI
}

fn rotation_subtree(
    bond_a: usize,
    bond_b: usize,
    bonds: &[(usize, usize)],
    n_atoms: usize,
) -> Vec<usize> {
    let mut adjacency = vec![Vec::new(); n_atoms];
    for &(i, j) in bonds {
        if i < n_atoms && j < n_atoms {
            adjacency[i].push(j);
            adjacency[j].push(i);
        }
    }

    let mut visited = vec![false; n_atoms];
    visited[bond_a] = true;
    visited[bond_b] = true;
    let mut stack = vec![bond_b];
    let mut subtree = vec![bond_b];

    while let Some(current) = stack.pop() {
        for &neighbor in &adjacency[current] {
            if !visited[neighbor] {
                visited[neighbor] = true;
                subtree.push(neighbor);
                stack.push(neighbor);
            }
        }
    }

    subtree
}

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
    let dft_method: sci_form::dft::ks_fock::DftMethod = match method.parse() {
        Ok(m) => m,
        Err(e) => return json_error(&e),
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
        Err(e) => return json_error(&bad_json("'elements'", elements_json, &e)),
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
        Err(e) => return json_error(&bad_json("'elements'", elements_json, &e)),
    };
    let params = sci_form::forcefield::reaxff::params::ReaxffParams::default_chon();
    let atom_params = build_reaxff_atom_params(&elems);
    let eem_params: Vec<_> = elems
        .iter()
        .map(|&z| sci_form::forcefield::reaxff::eem::default_eem_params(z))
        .collect();
    let charges = match sci_form::forcefield::reaxff::eem::solve_eem(&flat, &eem_params, 0.0) {
        Ok(charges) => charges,
        Err(e) => return json_error(&e),
    };
    let bo = sci_form::forcefield::reaxff::bond_order::compute_bond_orders(
        &flat,
        &atom_params,
        params.cutoff,
    );
    let bonded = sci_form::forcefield::reaxff::energy::compute_bonded_energy(&flat, &bo, &params);
    let (van_der_waals, coulomb) =
        sci_form::forcefield::reaxff::nonbonded::compute_nonbonded_energy(
            &flat,
            &charges,
            &elems,
            params.cutoff,
        );
    serde_json::json!({
        "bonded": bonded.total,
        "coulomb": coulomb,
        "van_der_waals": van_der_waals,
        "total": bonded.total + coulomb + van_der_waals
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
    let flat: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();
    let eem_params: Vec<_> = elems
        .iter()
        .map(|&z| sci_form::forcefield::reaxff::eem::default_eem_params(z))
        .collect();
    let charges = match sci_form::forcefield::reaxff::eem::solve_eem(&flat, &eem_params, 0.0) {
        Ok(charges) => charges,
        Err(e) => return json_error(&e),
    };
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
        Err(e) => return json_error(&bad_json("'mlff config'", config_json, &e)),
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
                Err(e) => return json_error(&bad_json("'params'", params_json, &e)),
            }
        };
    let aevs = sci_form::mlff::compute_aevs(&elems, &positions, &params);
    let aev_length = aevs.first().map(|v| v.len()).unwrap_or(0);
    let aev_arrays: Vec<Vec<f64>> = aevs.iter().map(|v| v.to_vec()).collect();
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
        serde_json::from_str(s).map_err(|e| bad_json("'shell'", s, &e))
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
        Err(e) => return json_error(&bad_json("'center a'", center_a_json, &e)),
    };
    let cb: [f64; 3] = match serde_json::from_str(center_b_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'center b'", center_b_json, &e)),
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
        Err(e) => return json_error(&bad_json("'indices'", subtree_indices_json, &e)),
    };
    let a: [f64; 3] = match serde_json::from_str(axis_a_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'axis a'", axis_a_json, &e)),
    };
    let b: [f64; 3] = match serde_json::from_str(axis_b_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'axis b'", axis_b_json, &e)),
    };

    let motor = sci_form::alpha::cga::dihedral_motor(a, b, angle_rad);
    let coords = sci_form::alpha::cga::apply_motor_to_subtree(&flat, &indices, &motor);
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
        Err(e) => return json_error(&bad_json("'torsion indices'", torsion_indices_json, &e)),
    };
    if torsion.len() != 4 {
        return json_error("torsion_indices must have exactly 4 elements");
    }
    let mol = match sci_form::parse(smiles) {
        Ok(m) => m,
        Err(e) => return json_error(&smiles_error(&e, smiles)),
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
    let i = torsion[0];
    let j = torsion[1];
    let k = torsion[2];
    let l = torsion[3];
    let current = dihedral_angle(
        [flat[i * 3], flat[i * 3 + 1], flat[i * 3 + 2]],
        [flat[j * 3], flat[j * 3 + 1], flat[j * 3 + 2]],
        [flat[k * 3], flat[k * 3 + 1], flat[k * 3 + 2]],
        [flat[l * 3], flat[l * 3 + 1], flat[l * 3 + 2]],
    );
    let mut delta = target_angle_rad - current;
    while delta <= -std::f64::consts::PI {
        delta += 2.0 * std::f64::consts::PI;
    }
    while delta > std::f64::consts::PI {
        delta -= 2.0 * std::f64::consts::PI;
    }
    let subtree = rotation_subtree(j, k, &bonds, elems.len());
    let motor = sci_form::alpha::cga::dihedral_motor(
        [flat[j * 3], flat[j * 3 + 1], flat[j * 3 + 2]],
        [flat[k * 3], flat[k * 3 + 1], flat[k * 3 + 2]],
        delta,
    );
    let new_coords = sci_form::alpha::cga::apply_motor_to_subtree(&flat, &subtree, &motor);
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
        Err(e) => return json_error(&bad_json("'reactant'", reactant_json, &e)),
    };
    let product: Vec<f64> = match serde_json::from_str(product_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'product'", product_json, &e)),
    };
    let node = sci_form::alpha::gsm::interpolate_node(&reactant, &product, t);
    serde_json::json!({ "coords": node }).to_string()
}

/// Plan which GSM backends are available for a molecule.
#[cfg(feature = "alpha-gsm")]
#[wasm_bindgen]
pub fn alpha_gsm_backend_plan(smiles: &str) -> String {
    match sci_form::alpha::gsm::plan_gsm_backends_for_smiles(smiles) {
        Ok(plan) => serde_json::to_string(&plan).unwrap_or_else(|e| json_error(&e.to_string())),
        Err(e) => json_error(&e),
    }
}

/// Compare one or more GSM backends on the same geometry.
///
/// `methods_json` can be `[]` to evaluate the full backend list.
#[cfg(feature = "alpha-gsm")]
#[wasm_bindgen]
pub fn alpha_gsm_compare_backends(
    smiles: &str,
    coords_flat_json: &str,
    methods_json: &str,
) -> String {
    let coords: Vec<f64> = match serde_json::from_str(coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'coords'", coords_flat_json, &e)),
    };
    let methods: Vec<sci_form::alpha::gsm::GsmEnergyBackend> =
        if methods_json.trim().is_empty() || methods_json == "null" {
            Vec::new()
        } else {
            let raw_methods: Vec<String> = match serde_json::from_str(methods_json) {
                Ok(v) => v,
                Err(e) => return json_error(&bad_json("'methods'", methods_json, &e)),
            };
            let mut parsed = Vec::with_capacity(raw_methods.len());
            for method in raw_methods {
                match method.parse::<sci_form::alpha::gsm::GsmEnergyBackend>() {
                    Ok(value) => parsed.push(value),
                    Err(err) => return json_error(&err),
                }
            }
            parsed
        };

    match sci_form::alpha::gsm::compare_gsm_backends(smiles, &coords, &methods) {
        Ok(result) => serde_json::to_string(&result).unwrap_or_else(|e| json_error(&e.to_string())),
        Err(e) => json_error(&e),
    }
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
    alpha_gsm_find_ts_with_method(
        smiles,
        reactant_coords_json,
        product_coords_json,
        config_json,
        "uff",
    )
}

/// Grow the reaction string and find the transition state node between
/// reactant and product geometries using the selected backend.
#[cfg(feature = "alpha-gsm")]
#[wasm_bindgen]
pub fn alpha_gsm_find_ts_with_method(
    smiles: &str,
    reactant_coords_json: &str,
    product_coords_json: &str,
    config_json: &str,
    method: &str,
) -> String {
    let reactant: Vec<f64> = match serde_json::from_str(reactant_coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'reactant'", reactant_coords_json, &e)),
    };
    let product: Vec<f64> = match serde_json::from_str(product_coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'product'", product_coords_json, &e)),
    };
    let config: sci_form::alpha::gsm::GsmConfig =
        if config_json.trim().is_empty() || config_json == "null" || config_json == "{}" {
            sci_form::alpha::gsm::GsmConfig::default()
        } else {
            match serde_json::from_str(config_json) {
                Ok(c) => c,
                Err(e) => return json_error(&bad_json("'gsm config'", config_json, &e)),
            }
        };

    let backend: sci_form::alpha::gsm::GsmEnergyBackend = match method.parse() {
        Ok(value) => value,
        Err(err) => return json_error(&err),
    };
    let r = match sci_form::alpha::gsm::find_transition_state_with_backend(
        smiles, &reactant, &product, backend, &config,
    ) {
        Ok(result) => result,
        Err(err) => return json_error(&err),
    };
    serde_json::json!({
        "backend": backend.as_str(),
        "ts_coords": r.ts_coords,
        "ts_energy": r.ts_energy,
        "activation_energy": r.activation_energy,
        "reverse_barrier": r.reverse_barrier,
        "path_energies": r.path_energies,
        "path_coords": r.path_coords,
        "ts_node_index": r.ts_node_index,
        "n_nodes": r.n_nodes,
        "energy_evaluations": r.energy_evaluations
    })
    .to_string()
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
        Err(e) => return json_error(&bad_json("'distance pairs'", distance_pairs_json, &e)),
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
                Err(e) => return json_error(&bad_json("'sdr config'", config_json, &e)),
            }
        };

    let result = sci_form::alpha::sdr::sdr_embed(n_atoms, &distance_pairs, &config);
    serde_json::json!({
        "coords": result.coords,
        "num_atoms": result.num_atoms,
        "convergence": {
            "iterations": result.convergence.iterations,
            "converged": result.convergence.converged,
            "final_residual": result.convergence.final_residual,
            "neg_eigenvalues_removed": result.convergence.neg_eigenvalues_removed
        },
        "max_distance_error": result.max_distance_error,
        "retries_avoided": result.retries_avoided
    })
    .to_string()
}

// ─── A8: Live MD Simulation (re-exported from dynamics_live module) ──────────
//
// The full `LiveSimulation` JS class is exposed in `crates/wasm/src/dynamics_live.rs`.
// This section adds metadata helpers for the alpha subpath export.

/// Return alpha module metadata.
///
/// The module list is derived from the features the binary was compiled
/// with, so it always describes what the package actually contains.
///
/// Returns JSON: `{ "alpha_modules": ["dft", "reaxff", …],
///                   "stability": "alpha — experimental, …",
///                   "version": "1.0.0" }`
#[wasm_bindgen]
pub fn alpha_modules_info() -> String {
    let mut modules: Vec<&str> = Vec::new();
    if cfg!(feature = "alpha-dft") {
        modules.push("dft");
    }
    if cfg!(feature = "alpha-reaxff") {
        modules.push("reaxff");
    }
    if cfg!(feature = "alpha-mlff") {
        modules.push("mlff");
    }
    if cfg!(feature = "alpha-obara-saika") {
        modules.push("obara-saika");
    }
    if cfg!(feature = "alpha-cga") {
        modules.push("cga");
    }
    if cfg!(feature = "alpha-gsm") {
        modules.push("gsm");
    }
    if cfg!(feature = "alpha-sdr") {
        modules.push("sdr");
    }
    if cfg!(feature = "alpha-dynamics-live") {
        modules.push("dynamics-live");
    }
    if cfg!(feature = "alpha-imd") {
        modules.push("imd");
    }
    if cfg!(feature = "alpha-edl") {
        modules.push("edl");
    }
    if cfg!(feature = "alpha-periodic-linear") {
        modules.push("periodic-linear");
    }
    if cfg!(feature = "alpha-kinetics") {
        modules.push("kinetics");
    }
    if cfg!(feature = "alpha-render-bridge") {
        modules.push("render-bridge");
    }
    serde_json::json!({
        "alpha_modules": modules,
        "stability": "alpha — experimental, subject to breaking changes",
        "version": env!("CARGO_PKG_VERSION"),
    })
    .to_string()
}

// ─── RENDER BRIDGE (charts / trajectories / transport) ─────────────────────
//
// Thin JSON wrappers over the Rust chart generators in
// `sci_form::alpha::render_bridge`. Each one is gated on the feature that
// provides its input type, mirroring the Rust-side `cfg`s.

#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-edl"))]
#[wasm_bindgen]
/// Build an EDL profile chart from an `EdlProfileResult` JSON payload.
///
/// @param profile_json  JSON `EdlProfileResult`
/// @returns JSON chart payload
pub fn alpha_render_edl_profile_chart(profile_json: &str) -> String {
    use sci_form::alpha::edl::EdlProfileResult;
    match serde_json::from_str::<EdlProfileResult>(profile_json) {
        Ok(profile) => {
            let chart = sci_form::alpha::render_bridge::edl_profile_chart(&profile);
            serde_json::to_string(&chart).unwrap_or_else(|e| json_error(&e.to_string()))
        }
        Err(e) => json_error(&bad_json("'EdlProfileResult'", profile_json, &e)),
    }
}

#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-edl"))]
#[wasm_bindgen]
/// Build a capacitance-scan chart from a `CpmEdlScanResult` JSON payload.
///
/// @param scan_json  JSON `CpmEdlScanResult`
/// @returns JSON chart payload
pub fn alpha_render_capacitance_scan_chart(scan_json: &str) -> String {
    use sci_form::alpha::edl::CpmEdlScanResult;
    match serde_json::from_str::<CpmEdlScanResult>(scan_json) {
        Ok(scan) => {
            let chart = sci_form::alpha::render_bridge::capacitance_scan_chart(&scan);
            serde_json::to_string(&chart).unwrap_or_else(|e| json_error(&e.to_string()))
        }
        Err(e) => json_error(&bad_json("'CpmEdlScanResult'", scan_json, &e)),
    }
}

#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
#[wasm_bindgen]
/// Build a periodic-DOS chart from a `PeriodicKpmDosResult` JSON payload.
///
/// @param dos_json  JSON `PeriodicKpmDosResult`
/// @returns JSON chart payload
pub fn alpha_render_periodic_dos_chart(dos_json: &str) -> String {
    use sci_form::alpha::periodic_linear::PeriodicKpmDosResult;
    match serde_json::from_str::<PeriodicKpmDosResult>(dos_json) {
        Ok(dos) => {
            let chart = sci_form::alpha::render_bridge::periodic_dos_chart(&dos);
            serde_json::to_string(&chart).unwrap_or_else(|e| json_error(&e.to_string()))
        }
        Err(e) => json_error(&bad_json("'PeriodicKpmDosResult'", dos_json, &e)),
    }
}

#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
#[wasm_bindgen]
/// Build a band-structure chart from a `BandStructureAdapterResult` JSON payload.
///
/// @param bs_json  JSON `BandStructureAdapterResult`
/// @returns JSON chart payload
pub fn alpha_render_band_structure_chart(bs_json: &str) -> String {
    use sci_form::alpha::periodic_linear::BandStructureAdapterResult;
    match serde_json::from_str::<BandStructureAdapterResult>(bs_json) {
        Ok(bs) => {
            let chart = sci_form::alpha::render_bridge::band_structure_chart(&bs);
            serde_json::to_string(&chart).unwrap_or_else(|e| json_error(&e.to_string()))
        }
        Err(e) => json_error(&bad_json("'BandStructureAdapterResult'", bs_json, &e)),
    }
}

#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-kinetics"))]
#[wasm_bindgen]
/// Build a microkinetic population chart from a `MicrokineticTrace` JSON payload.
///
/// @param trace_json  JSON `MicrokineticTrace`
/// @returns JSON chart payload
pub fn alpha_render_microkinetic_population_chart(trace_json: &str) -> String {
    use sci_form::alpha::kinetics::MicrokineticTrace;
    match serde_json::from_str::<MicrokineticTrace>(trace_json) {
        Ok(trace) => {
            let chart = sci_form::alpha::render_bridge::microkinetic_population_chart(&trace);
            serde_json::to_string(&chart).unwrap_or_else(|e| json_error(&e.to_string()))
        }
        Err(e) => json_error(&bad_json("'MicrokineticTrace'", trace_json, &e)),
    }
}

#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-kinetics"))]
#[wasm_bindgen]
/// Build an Arrhenius chart from a JSON array of `ElementaryRateResult`.
///
/// @param results_json  JSON array of `ElementaryRateResult`
/// @returns JSON chart payload
pub fn alpha_render_arrhenius_chart(results_json: &str) -> String {
    use sci_form::alpha::kinetics::ElementaryRateResult;
    match serde_json::from_str::<Vec<ElementaryRateResult>>(results_json) {
        Ok(results) => {
            let chart = sci_form::alpha::render_bridge::arrhenius_chart(&results);
            serde_json::to_string(&chart).unwrap_or_else(|e| json_error(&e.to_string()))
        }
        Err(e) => json_error(&bad_json("'ElementaryRateResult'", results_json, &e)),
    }
}

#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-kinetics"))]
#[wasm_bindgen]
/// Build GSM path trajectory frames from coords + energies JSON payloads.
///
/// @param coords_json  JSON array of flat coordinate frames
/// @param energies_json  JSON array of frame energies
/// @returns JSON array of trajectory frames
pub fn alpha_render_gsm_path_trajectory(coords_json: &str, energies_json: &str) -> String {
    let coords: Vec<Vec<f64>> = match serde_json::from_str(coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'coords'", coords_json, &e)),
    };
    let energies: Vec<f64> = match serde_json::from_str(energies_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'energies'", energies_json, &e)),
    };
    let frames = sci_form::alpha::render_bridge::gsm_path_trajectory(&coords, &energies);
    serde_json::to_string(&frames).unwrap_or_else(|e| json_error(&e.to_string()))
}

#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
/// Pack a chart payload into the columnar transport format.
///
/// @param chart_json  JSON `ChartPayload`
/// @returns JSON packed batch
pub fn alpha_render_pack_chart_payload(chart_json: &str) -> String {
    use sci_form::alpha::render_bridge::ChartPayload;
    match serde_json::from_str::<ChartPayload>(chart_json) {
        Ok(chart) => {
            let batch = sci_form::alpha::render_bridge::pack_chart_payload(&chart);
            serde_json::to_string(&batch).unwrap_or_else(|e| json_error(&e.to_string()))
        }
        Err(e) => json_error(&bad_json("'ChartPayload'", chart_json, &e)),
    }
}

#[cfg(feature = "alpha-render-bridge")]
#[wasm_bindgen]
/// Verify a chart payload survives a pack/unpack round trip.
///
/// @param chart_json  JSON `ChartPayload`
/// @param tol  Comparison tolerance
/// @returns `true` when the round trip is lossless within `tol`
pub fn alpha_render_verify_chart_round_trip(chart_json: &str, tol: f64) -> bool {
    use sci_form::alpha::render_bridge::ChartPayload;
    match serde_json::from_str::<ChartPayload>(chart_json) {
        Ok(chart) => sci_form::alpha::render_bridge::verify_chart_round_trip(&chart, tol),
        Err(_) => false,
    }
}

// ─── PERIODIC LINEAR (K-meshes, BZ integration, EHT band structure) ────────

#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
/// Build a Monkhorst-Pack (or gamma-centred) k-mesh.
///
/// @param config_json  JSON `KMeshConfig` or `"{}"`
/// @returns JSON k-mesh
pub fn alpha_periodic_monkhorst_pack(config_json: &str) -> String {
    use sci_form::alpha::periodic_linear::{monkhorst_pack_mesh, KMeshConfig};
    match serde_json::from_str::<KMeshConfig>(config_json) {
        Ok(config) => match monkhorst_pack_mesh(&config) {
            Ok(mesh) => serde_json::to_string(&mesh).unwrap_or_else(|e| json_error(&e.to_string())),
            Err(e) => json_error(&e),
        },
        Err(e) => json_error(&bad_json("'KMeshConfig'", config_json, &e)),
    }
}

#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
/// Bloch phase for a fractional k-point and a lattice translation.
///
/// @param k_fractional_json  JSON `[kx, ky, kz]`
/// @param translation_json  JSON `[tx, ty, tz]`
/// @returns JSON Bloch phase
pub fn alpha_periodic_bloch_phase(k_fractional_json: &str, translation_json: &str) -> String {
    use sci_form::alpha::periodic_linear::bloch_phase;
    let k: [f64; 3] = match serde_json::from_str(k_fractional_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'k-point'", k_fractional_json, &e)),
    };
    let r: [i32; 3] = match serde_json::from_str(translation_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'translation'", translation_json, &e)),
    };
    let (re, im) = bloch_phase(k, r);
    serde_json::json!({ "re": re, "im": im }).to_string()
}

#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
/// Brillouin-zone weighted integral of a scalar field.
///
/// @param values_json  JSON array of scalars
/// @param weights_json  JSON array of k-point weights
/// @returns the weighted integral
pub fn alpha_periodic_bz_integrate_scalar(values_json: &str, weights_json: &str) -> f64 {
    use sci_form::alpha::periodic_linear::bz_integrate_scalar;
    let values: Vec<f64> = serde_json::from_str(values_json).unwrap_or_default();
    let weights: Vec<f64> = serde_json::from_str(weights_json).unwrap_or_default();
    bz_integrate_scalar(&values, &weights)
}

#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
/// Brillouin-zone weighted integral of a vector field.
///
/// @param vectors_json  JSON array of `[x, y, z]`
/// @param weights_json  JSON array of k-point weights
/// @returns JSON `[x, y, z]` integral
pub fn alpha_periodic_bz_integrate_vector(vectors_json: &str, weights_json: &str) -> String {
    use sci_form::alpha::periodic_linear::bz_integrate_vector;
    let vectors: Vec<Vec<f64>> = serde_json::from_str(vectors_json).unwrap_or_default();
    let weights: Vec<f64> = serde_json::from_str(weights_json).unwrap_or_default();
    serde_json::to_string(&bz_integrate_vector(&vectors, &weights))
        .unwrap_or_else(|e| json_error(&e.to_string()))
}

#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
/// Validate that the k-mesh weights sum to one.
///
/// @param kmesh_json  JSON `KMesh`
/// @param tol  Tolerance
/// @returns `true` when the weights are normalised within `tol`
pub fn alpha_periodic_validate_kmesh_weights(kmesh_json: &str, tol: f64) -> bool {
    use sci_form::alpha::periodic_linear::{validate_kmesh_weights, KMesh};
    match serde_json::from_str::<KMesh>(kmesh_json) {
        Ok(mesh) => validate_kmesh_weights(&mesh, tol),
        Err(_) => false,
    }
}

/// EHT band structure along a k-path (`eht::band_structure::compute_band_structure`).
#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
pub fn alpha_periodic_band_structure_eht(
    elements_json: &str,
    coords_flat_json: &str,
    lattice_json: &str,
    config_json: &str,
    n_electrons: usize,
) -> String {
    use sci_form::eht::band_structure::{compute_band_structure, BandStructureConfig};
    let (elements, positions) = match parse_elements_and_positions(elements_json, coords_flat_json)
    {
        Ok(v) => v,
        Err(e) => return e,
    };
    let lattice: [[f64; 3]; 3] = match serde_json::from_str(lattice_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'lattice'", lattice_json, &e)),
    };
    let config: BandStructureConfig = match serde_json::from_str(config_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'BandStructureConfig'", config_json, &e)),
    };
    match compute_band_structure(&elements, &positions, &lattice, &config, n_electrons) {
        Ok(bs) => serde_json::to_string(&bs).unwrap_or_else(|e| json_error(&e.to_string())),
        Err(e) => json_error(&e),
    }
}

#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
/// Validate the electron count implied by a periodic DOS result.
///
/// @param dos_json  JSON `PeriodicKpmDosResult`
/// @param fermi_energy_ev  Fermi energy in eV
/// @param expected_electrons  Expected number of electrons
/// @returns the number of electrons found
pub fn alpha_periodic_validate_electron_count(
    dos_json: &str,
    fermi_energy_ev: f64,
    expected_electrons: f64,
) -> f64 {
    use sci_form::alpha::periodic_linear::{validate_electron_count, PeriodicKpmDosResult};
    match serde_json::from_str::<PeriodicKpmDosResult>(dos_json) {
        Ok(dos) => validate_electron_count(&dos, fermi_energy_ev, expected_electrons),
        Err(_) => f64::NAN,
    }
}

/// One real-space operator block: the matrix and its lattice translation.
#[cfg(feature = "alpha-periodic-linear")]
type RealSpaceBlock = (nalgebra::DMatrix<f64>, [i32; 3]);

/// Parse a JSON array of real-space operator blocks into matrices + translations.
#[cfg(feature = "alpha-periodic-linear")]
fn parse_real_space_blocks(blocks_json: &str) -> Result<Vec<RealSpaceBlock>, String> {
    let value: serde_json::Value = serde_json::from_str(blocks_json)
        .map_err(|e| format!("invalid real-space block JSON: {e}"))?;
    let array = value
        .as_array()
        .ok_or_else(|| "real-space blocks must be a JSON array".to_string())?;
    let mut blocks = Vec::with_capacity(array.len());
    for item in array {
        let matrix = item
            .get("matrix")
            .and_then(|m| m.as_array())
            .ok_or_else(|| "each real-space block needs a 'matrix' array".to_string())?;
        let translation = item
            .get("translation")
            .and_then(|t| t.as_array())
            .ok_or_else(|| "each real-space block needs a 'translation' array".to_string())?;
        if translation.len() != 3 {
            return Err("lattice translations need exactly three components".into());
        }
        let mut trans = [0i32; 3];
        for (slot, value) in trans.iter_mut().zip(translation.iter()) {
            *slot = value
                .as_i64()
                .ok_or_else(|| "lattice translations must be integers".to_string())?
                as i32;
        }
        let n = matrix.len();
        if n == 0 {
            return Err("real-space matrices must be non-empty".into());
        }
        let mut flat = Vec::with_capacity(n * n);
        for row in matrix {
            let row = row
                .as_array()
                .ok_or_else(|| "real-space matrix rows must be arrays".to_string())?;
            if row.len() != n {
                return Err("real-space matrices must be square".into());
            }
            for entry in row {
                flat.push(
                    entry
                        .as_f64()
                        .ok_or_else(|| "real-space matrix entries must be numbers".to_string())?,
                );
            }
        }
        blocks.push((nalgebra::DMatrix::from_row_slice(n, n, &flat), trans));
    }
    if blocks.is_empty() {
        return Err("at least one real-space block is required".into());
    }
    Ok(blocks)
}

#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
/// Interpolate the default high-symmetry k-path (`Γ → X → M → Γ`).
///
/// @param n_points_per_segment  Intervals per high-symmetry segment (≥ 1)
/// @returns JSON `KPathStandard`
pub fn alpha_periodic_kpath_standard(n_points_per_segment: usize) -> String {
    use sci_form::alpha::periodic_linear::kpath_standard;
    serde_json::to_string(&kpath_standard(n_points_per_segment))
        .unwrap_or_else(|e| json_error(&e.to_string()))
}

#[cfg(feature = "alpha-periodic-linear")]
#[wasm_bindgen]
/// Trace the Fermi-surface contour of real-space operator blocks in `kz = 0`.
///
/// @param blocks_json  JSON array of `{ "matrix": [[…]], "translation": [i, j, k] }`
/// @param fermi_energy_ev  Fermi energy in eV (`null`/`undefined` = median estimate)
/// @param n_kx  Number of kx samples (≥ 2)
/// @param n_ky  Number of ky samples (≥ 2)
/// @returns JSON `FermiSurfaceContour`
pub fn alpha_periodic_fermi_surface_from_bands(
    blocks_json: &str,
    fermi_energy_ev: Option<f64>,
    n_kx: usize,
    n_ky: usize,
) -> String {
    use sci_form::alpha::periodic_linear::fermi_surface_from_bands;
    let blocks = match parse_real_space_blocks(blocks_json) {
        Ok(blocks) => blocks,
        Err(e) => return json_error(&e),
    };
    let refs: Vec<(&nalgebra::DMatrix<f64>, [i32; 3])> = blocks
        .iter()
        .map(|(matrix, trans)| (matrix, *trans))
        .collect();
    let contour = fermi_surface_from_bands(&refs, fermi_energy_ev, n_kx, n_ky);
    serde_json::to_string(&contour).unwrap_or_else(|e| json_error(&e.to_string()))
}
