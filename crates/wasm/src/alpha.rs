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
        Err(e) => return json_error(&format!("bad torsion_indices: {}", e)),
    };
    if torsion.len() != 4 {
        return json_error("torsion_indices must have exactly 4 elements");
    }
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
        Err(e) => return json_error(&format!("bad reactant: {}", e)),
    };
    let product: Vec<f64> = match serde_json::from_str(product_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad product: {}", e)),
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
        Err(e) => return json_error(&format!("bad coords: {}", e)),
    };
    let methods: Vec<sci_form::alpha::gsm::GsmEnergyBackend> =
        if methods_json.trim().is_empty() || methods_json == "null" {
            Vec::new()
        } else {
            let raw_methods: Vec<String> = match serde_json::from_str(methods_json) {
                Ok(v) => v,
                Err(e) => return json_error(&format!("bad methods: {}", e)),
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
            "dynamics-live",
            "reaction-dynamics"
        ],
        "stability": "alpha — experimental, subject to breaking changes",
        "version": "0.14.4"
    })
    .to_string()
}

// ─── A9: Alpha Reaction Dynamics 3D ─────────────────────────────────────────

/// Compute a full 3D reaction dynamics path using the alpha pipeline.
///
/// Uses CI-NEB with IDPP initialisation, constrained geometry relaxation,
/// orbital/electrostatic approach guidance, and optional SMIRKS atom mapping.
///
/// `reactant_smiles_json`: JSON array of reactant SMILES, e.g. `["[Cl-]", "CBr"]`.
/// `product_smiles_json`:  JSON array of product SMILES, e.g. `["ClC", "[Br-]"]`.
/// `config_json`: optional JSON config; pass `""` or `"{}"` for defaults (GFN2-xTB).
///
/// Returns JSON with frames (approach + NEB + departure), energies, TS info.
#[cfg(feature = "alpha-reaction-dynamics")]
#[wasm_bindgen]
pub fn compute_reaction_dynamics_3d(
    reactant_smiles_json: &str,
    product_smiles_json: &str,
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

    let config: sci_form::alpha::reaction_dynamics::ReactionDynamics3DConfig =
        if config_json.is_empty() || config_json == "{}" {
            sci_form::alpha::reaction_dynamics::ReactionDynamics3DConfig::default()
        } else {
            match serde_json::from_str(config_json) {
                Ok(c) => c,
                Err(e) => return json_error(&format!("bad config: {}", e)),
            }
        };

    let r_refs: Vec<&str> = reactants.iter().map(String::as_str).collect();
    let p_refs: Vec<&str> = products.iter().map(String::as_str).collect();

    match sci_form::alpha::reaction_dynamics::compute_reaction_dynamics_3d(
        &r_refs, &p_refs, &config,
    ) {
        Ok(result) => serde_json::to_string(&result).unwrap_or_else(|e| json_error(&e.to_string())),
        Err(e) => json_error(&e),
    }
}

/// Reactant-only approach PES scan — no product SMILES required.
///
/// Embeds the reactant fragments, identifies reactive sites using Fukui/FMO
/// and electrostatic analysis, then performs a constrained-relaxed approach
/// scan from `far_distance` down to `reactive_distance`.
///
/// This is the correct workflow when only reactants are known:
/// sci-form builds a physically meaningful approach trajectory with real
/// quantum energies without needing to specify the products.
///
/// `reactant_smiles_json`: JSON array of reactant SMILES, e.g. `["CC(=O)O", "N"]`.
/// `config_json`: optional JSON [`ReactionDynamics3DConfig`]; pass `""` for defaults.
///
/// Returns JSON `ReactionDynamics3DResult` with approach frames only
/// (`phase = "approach"`). `reaction_energy_kcal_mol` is `0.0` (no products known).
#[cfg(feature = "alpha-reaction-dynamics")]
#[wasm_bindgen]
pub fn compute_approach_dynamics(reactant_smiles_json: &str, config_json: &str) -> String {
    let reactants: Vec<String> = match serde_json::from_str(reactant_smiles_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad reactant_smiles: {}", e)),
    };

    let mut config: sci_form::alpha::reaction_dynamics::ReactionDynamics3DConfig =
        if config_json.is_empty() || config_json == "{}" {
            sci_form::alpha::reaction_dynamics::ReactionDynamics3DConfig::default()
        } else {
            match serde_json::from_str(config_json) {
                Ok(c) => c,
                Err(e) => return json_error(&format!("bad config: {}", e)),
            }
        };

    // For approach-only, orbital and electrostatic guidance default to true
    // so the user gets a physically meaningful direction without NEB.
    if config_json.is_empty() || config_json == "{}" {
        config.use_orbital_guidance = true;
        config.use_electrostatic_steering = true;
    }

    let r_refs: Vec<&str> = reactants.iter().map(String::as_str).collect();

    match sci_form::alpha::reaction_dynamics::compute_approach_dynamics(&r_refs, &config) {
        Ok(result) => serde_json::to_string(&result).unwrap_or_else(|e| json_error(&e.to_string())),
        Err(e) => json_error(&e),
    }
}
