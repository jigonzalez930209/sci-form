//! **BETA** — WASM bindings for all beta-tier experimental modules.
//!
//! Import separately from the stable API:
//!
//! ```js
//! import { compute_kpm_dos, compute_mbh_frequencies, solve_eht_randnla_wasm }
//!   from 'sci-form-wasm/beta';
//! ```
//!
//! Each function is individually gated; only the ones compiled with the
//! matching feature flag are available in a given build.

#[allow(unused_imports)]
use crate::helpers::{json_error, parse_elements_and_positions, parse_flat_coords};
#[allow(unused_imports)]
use wasm_bindgen::prelude::*;

fn build_eht_matrices(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> (
    Vec<sci_form::eht::basis::AtomicOrbital>,
    nalgebra::DMatrix<f64>,
    nalgebra::DMatrix<f64>,
) {
    let basis = sci_form::eht::basis::build_basis(elements, positions);
    let overlap = sci_form::eht::build_overlap_matrix(&basis);
    let hamiltonian = sci_form::eht::build_hamiltonian(&basis, &overlap, None);
    (basis, hamiltonian, overlap)
}

// ─── B1: KPM — Kernel Polynomial Method ─────────────────────────────────────

/// Compute EHT Density of States using the Kernel Polynomial Method (O(N)).
///
/// Builds the EHT Hamiltonian internally from geometry, then runs KPM.
///
/// `config_json` fields (all optional):
/// - `order` (int, default 100): Chebyshev expansion order
/// - `n_points` (int, default 200): Energy grid points
/// - `e_min` (float): Override spectral minimum
/// - `e_max` (float): Override spectral maximum
/// - `temperature` (float, default 0): Electronic temperature in K
///
/// Returns JSON `KpmDosResult`:
/// ```json
/// { "energies": [...], "total_dos": [...], "e_min": -30.0, "e_max": 5.0,
///   "order": 100 }
/// ```
#[cfg(feature = "beta-kpm")]
#[wasm_bindgen]
pub fn beta_compute_kpm_dos(
    elements_json: &str,
    coords_flat_json: &str,
    config_json: &str,
) -> String {
    use sci_form::beta::kpm::{compute_kpm_dos, KpmConfig};

    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };

    let (_basis, hamiltonian, _overlap) = build_eht_matrices(&elems, &positions);
    if hamiltonian.nrows() == 0 {
        return json_error("empty Hamiltonian");
    }

    // Optional config override
    #[derive(serde::Deserialize, Default)]
    struct KpmConfigInput {
        order: Option<usize>,
        n_points: Option<usize>,
        e_min: Option<f64>,
        e_max: Option<f64>,
        temperature: Option<f64>,
    }
    let user: KpmConfigInput = serde_json::from_str(config_json).unwrap_or_default();

    let config = KpmConfig {
        order: user.order.unwrap_or(100),
        n_vectors: 0, // exact
        seed: 42,
        temperature: user.temperature.unwrap_or(0.0),
    };

    let result = compute_kpm_dos(
        &hamiltonian,
        &config,
        user.e_min.unwrap_or(-30.0),
        user.e_max.unwrap_or(5.0),
        user.n_points.unwrap_or(200),
    );

    serde_json::json!({
        "energies": result.energies,
        "total_dos": result.total_dos,
        "e_min": result.e_min,
        "e_max": result.e_max,
        "order": result.order
    })
    .to_string()
}

/// Compute EHT Mulliken populations using KPM (O(N) stochastic trace).
///
/// Returns JSON: `{ "mulliken_charges": [...], "orbital_populations": [...] }`
#[cfg(feature = "beta-kpm")]
#[wasm_bindgen]
pub fn beta_compute_kpm_mulliken(elements_json: &str, coords_flat_json: &str) -> String {
    use sci_form::beta::kpm::{compute_kpm_mulliken, KpmConfig};

    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let (basis, hamiltonian, overlap) = build_eht_matrices(&elems, &positions);
    let eht = match sci_form::eht::solve_eht(&elems, &positions, None) {
        Ok(r) => r,
        Err(e) => return json_error(&format!("EHT failed: {}", e)),
    };

    let config = KpmConfig::default();
    let nuclear_charges: Vec<f64> = basis
        .iter()
        .map(|orbital| elems[orbital.atom_index] as f64)
        .collect();
    let result = compute_kpm_mulliken(
        &hamiltonian,
        &overlap,
        eht.n_electrons,
        &nuclear_charges,
        &config,
    );

    let mut mulliken_charges = vec![0.0; elems.len()];
    let mut orbital_populations = vec![0.0; basis.len()];
    for (idx, orbital) in basis.iter().enumerate() {
        let charge = result.charges.get(idx).copied().unwrap_or(0.0);
        mulliken_charges[orbital.atom_index] += charge;
        orbital_populations[idx] = nuclear_charges[idx] - charge;
    }

    serde_json::json!({
        "mulliken_charges": mulliken_charges,
        "orbital_populations": orbital_populations
    })
    .to_string()
}

// ─── B2: MBH — Mobile Block Hessian ─────────────────────────────────────────

/// Compute vibrational frequencies using Mobile Block Hessian (reduced DOF).
///
/// Uses UFF as the energy evaluator internally.
///
/// Returns JSON `MbhResult`:
/// ```json
/// { "frequencies": [...], "n_blocks": 3, "n_flexible": 2,
///   "n_dof_reduced": 12, "n_dof_full": 24, "speedup": 2.0 }
/// ```
#[cfg(feature = "beta-mbh")]
#[wasm_bindgen]
pub fn beta_compute_mbh_frequencies(
    elements_json: &str,
    coords_flat_json: &str,
    smiles: &str,
) -> String {
    use sci_form::beta::mbh::hessian::{compute_mbh_frequencies, MbhConfig};

    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };

    // Get ring info from SMILES
    let rings: Vec<(Vec<usize>, bool)> = match sci_form::compute_sssr(smiles) {
        Ok(sssr) => sssr
            .rings
            .iter()
            .map(|r| (r.atoms.clone(), r.is_aromatic))
            .collect(),
        Err(_) => vec![],
    };

    let config = MbhConfig::default();

    // UFF energy evaluator
    let energy_fn =
        |coords: &[f64]| -> f64 { sci_form::compute_uff_energy(smiles, coords).unwrap_or(0.0) };

    let result = compute_mbh_frequencies(&elems, &positions, &rings, &energy_fn, &config);

    serde_json::json!({
        "frequencies": result.frequencies,
        "n_blocks": result.n_blocks,
        "n_flexible": result.n_flexible,
        "n_dof_reduced": result.n_dof_reduced,
        "n_dof_full": result.n_dof_full,
        "speedup": result.speedup
    })
    .to_string()
}

// ─── B3: RandNLA — Randomized EHT ───────────────────────────────────────────

/// Solve EHT with a randomized Nyström eigensolver instead of O(N³) exact diag.
///
/// `config_json` fields (all optional):
/// - `sketch_size` (int): Nyström rank k. Default: √N
/// - `seed` (int, default 42)
/// - `max_error` (float, default 0.001): Max relative residual before fallback
/// - `fallback_enabled` (bool, default true)
///
/// Returns JSON:
/// ```json
/// { "orbital_energies": [...], "homo_index": 4, "homo_energy": -0.5,
///   "lumo_energy": 0.1, "gap": 0.6,
///   "k": 12, "residual_error": 0.0003, "used_fallback": false }
/// ```
#[cfg(feature = "beta-randnla")]
#[wasm_bindgen]
pub fn beta_solve_eht_randnla(
    elements_json: &str,
    coords_flat_json: &str,
    config_json: &str,
) -> String {
    use sci_form::beta::rand_nla::{solve_eht_randnla, RandNlaConfig};

    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let (_basis, hamiltonian, overlap) = build_eht_matrices(&elems, &positions);
    let eht = match sci_form::eht::solve_eht(&elems, &positions, None) {
        Ok(r) => r,
        Err(e) => return json_error(&format!("EHT failed: {}", e)),
    };

    let cfg: RandNlaConfig =
        if config_json.trim().is_empty() || config_json == "null" || config_json == "{}" {
            RandNlaConfig::default()
        } else {
            match serde_json::from_str(config_json) {
                Ok(c) => c,
                Err(e) => return json_error(&format!("bad config: {}", e)),
            }
        };

    let (eigenvalues, _eigenvectors, info) = solve_eht_randnla(&hamiltonian, &overlap, &cfg);

    let energies: Vec<f64> = eigenvalues.iter().cloned().collect();
    let n_occ = eht.n_electrons / 2;
    let homo_idx = n_occ.saturating_sub(1);
    let lumo_idx = n_occ.min(energies.len().saturating_sub(1));
    let homo = energies.get(homo_idx).copied().unwrap_or(0.0);
    let lumo = energies.get(lumo_idx).copied().unwrap_or(0.0);

    serde_json::json!({
        "orbital_energies": energies,
        "homo_index": homo_idx,
        "homo_energy": homo,
        "lumo_energy": lumo,
        "gap": (lumo - homo).max(0.0),
        "k": info.k,
        "residual_error": info.residual_error,
        "used_fallback": info.used_fallback
    })
    .to_string()
}

// ─── B4: Riemannian Geometry ─────────────────────────────────────────────────

/// Compute Riemannian distance between two PSD matrices.
///
/// `matrix_a_json` / `matrix_b_json` — flat row-major arrays with `dim` dimension.
///
/// Returns JSON: `{ "distance": 0.12 }`
#[cfg(feature = "beta-riemannian")]
#[wasm_bindgen]
pub fn beta_psd_distance(matrix_a_json: &str, matrix_b_json: &str, dim: usize) -> String {
    let flat_a: Vec<f64> = match serde_json::from_str(matrix_a_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad matrix_a: {}", e)),
    };
    let flat_b: Vec<f64> = match serde_json::from_str(matrix_b_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad matrix_b: {}", e)),
    };
    if flat_a.len() != dim * dim || flat_b.len() != dim * dim {
        return json_error("matrix flat length must equal dim*dim");
    }
    let a = nalgebra::DMatrix::from_row_slice(dim, dim, &flat_a);
    let b = nalgebra::DMatrix::from_row_slice(dim, dim, &flat_b);
    let d = sci_form::beta::riemannian::psd_distance(&a, &b);
    serde_json::json!({ "distance": d }).to_string()
}

/// Project a symmetric matrix onto the PSD (positive semidefinite) cone.
///
/// `matrix_json` — flat row-major array with `dim` dimension.
///
/// Returns JSON: `{ "matrix": [...flat psd projection...], "was_psd": false }`
#[cfg(feature = "beta-riemannian")]
#[wasm_bindgen]
pub fn beta_psd_projection(matrix_json: &str, dim: usize) -> String {
    let flat: Vec<f64> = match serde_json::from_str(matrix_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad matrix: {}", e)),
    };
    if flat.len() != dim * dim {
        return json_error("matrix flat length must equal dim*dim");
    }
    let m = nalgebra::DMatrix::from_row_slice(dim, dim, &flat);
    let projected = sci_form::beta::riemannian::psd_projection(&m);
    let was_psd = (&m - &projected).norm() < 1e-10;
    serde_json::json!({
        "matrix": projected.as_slice().to_vec(),
        "was_psd": was_psd
    })
    .to_string()
}

// ─── B5: CPM — Constant Potential Method ────────────────────────────────────

/// Compute constant-potential charges using the CPM model.
///
/// `potential` — target electrode potential in eV.
///
/// Returns JSON: `{ "charges": [...], "total_charge": 0.5, "energy": -1.2 }`
#[cfg(feature = "beta-cpm")]
#[wasm_bindgen]
pub fn beta_compute_cpm_charges(
    elements_json: &str,
    coords_flat_json: &str,
    potential: f64,
) -> String {
    use sci_form::beta::cpm::grand_potential::{compute_cpm_charges, CpmConfig};

    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let config = CpmConfig {
        mu_ev: potential,
        ..Default::default()
    };
    let r = compute_cpm_charges(&elems, &positions, &config);
    serde_json::json!({
        "charges": r.charges,
        "total_charge": r.total_charge,
        "energy": r.grand_potential
    })
    .to_string()
}

// ─── Metadata ────────────────────────────────────────────────────────────────

/// Return beta module metadata.
///
/// Returns JSON: `{ "beta_modules": ["kpm", "mbh", "randnla", "riemannian", "cpm"] }`
#[wasm_bindgen]
pub fn beta_modules_info() -> String {
    serde_json::json!({
        "beta_modules": [
            "kpm",
            "mbh",
            "randnla",
            "riemannian",
            "cpm"
        ],
        "stability": "beta — validated, approaching stable",
        "version": "0.11.2"
    })
    .to_string()
}
