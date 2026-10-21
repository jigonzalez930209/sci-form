//! Electronic structure method WASM bindings: PM3, xTB, GFN1, GFN2, HF-3c, ANI.

use crate::helpers::*;
use wasm_bindgen::prelude::*;

/// Run a PM3 semi-empirical calculation.
#[wasm_bindgen]
pub fn compute_pm3(elements_json: &str, coords_flat_json: &str) -> String {
    let (elements, positions) = match parse_elements_and_positions(elements_json, coords_flat_json)
    {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_pm3(&elements, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Run an xTB tight-binding calculation.
#[wasm_bindgen]
pub fn compute_xtb(elements_json: &str, coords_flat_json: &str) -> String {
    let (elements, positions) = match parse_elements_and_positions(elements_json, coords_flat_json)
    {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_xtb(&elements, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Run a GFN1-xTB tight-binding calculation.
#[wasm_bindgen]
pub fn compute_gfn1(elements_json: &str, coords_flat_json: &str) -> String {
    let (elements, positions) = match parse_elements_and_positions(elements_json, coords_flat_json)
    {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::xtb::solve_gfn1(&elements, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Run a GFN2-xTB tight-binding calculation.
#[wasm_bindgen]
pub fn compute_gfn2(elements_json: &str, coords_flat_json: &str) -> String {
    let (elements, positions) = match parse_elements_and_positions(elements_json, coords_flat_json)
    {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::xtb::solve_gfn2(&elements, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Run a HF-3c calculation with default configuration.
#[wasm_bindgen]
pub fn compute_hf3c(elements_json: &str, coords_flat_json: &str) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let config = sci_form::hf::HfConfig::default();
    match sci_form::compute_hf3c(&elems, &positions, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Run a HF-3c calculation with custom configuration.
#[wasm_bindgen]
pub fn compute_hf3c_custom(
    elements_json: &str,
    coords_flat_json: &str,
    max_scf_iter: usize,
    n_cis_states: usize,
    corrections: bool,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let config = sci_form::hf::HfConfig {
        max_scf_iter,
        diis_size: 6,
        n_cis_states,
        corrections,
    };
    match sci_form::compute_hf3c(&elems, &positions, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Compute ANI neural-network potential energy.
#[wasm_bindgen]
pub fn compute_ani(elements_json: &str, coords_flat_json: &str) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_ani(&elems, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Compute CISD excited states from an internal STO-3G RHF calculation.
///
/// Runs the shared SCF engine and the validated two-electron integral engine
/// internally, then feeds the canonical orbitals and ERIs to the CISD driver.
#[wasm_bindgen]
pub fn compute_cisd(elements_json: &str, coords_flat_json: &str, n_states: usize) -> String {
    let (elements, positions) = match parse_elements_and_positions(elements_json, coords_flat_json)
    {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let system = sci_form::scf::types::MolecularSystem::from_angstrom(&elements, &positions, 0, 1);
    let scf_config = sci_form::scf::scf_loop::ScfConfig::default();
    let scf = sci_form::scf::scf_loop::run_scf(&system, &scf_config);
    if !scf.converged {
        return json_error("RHF/STO-3G SCF did not converge; CISD requires a converged reference");
    }

    let basis =
        sci_form::scf::basis::BasisSet::sto3g(&system.atomic_numbers, &system.positions_bohr);
    let n_basis = basis.n_basis;
    let engine = sci_form::scf::two_electron::TwoElectronIntegrals::compute(&basis);
    let n_pair = n_basis * (n_basis + 1) / 2;
    let mut eris = vec![0.0f64; n_pair * (n_pair + 1) / 2];
    for i in 0..n_basis {
        for j in 0..n_basis {
            for k in 0..n_basis {
                for l in 0..n_basis {
                    eris[sci_form::hf::integrals::eri_index(i, j, k, l, n_basis)] =
                        engine.get(i, j, k, l);
                }
            }
        }
    }

    let n_occupied = scf.n_electrons / 2;
    let result = sci_form::hf::compute_cisd(
        &scf.orbital_energies,
        &scf.mo_coefficients,
        &eris,
        n_basis,
        n_occupied,
        n_states,
    );
    serialize_or_error(&result)
}

/// Compute ANI-TM atomic environment vectors (24-element extended set).
///
/// Uses the built-in transition-metal AEV parameter set.
#[wasm_bindgen]
pub fn compute_aevs_tm(elements_json: &str, coords_flat_json: &str) -> String {
    let (elements, positions) = match parse_elements_and_positions(elements_json, coords_flat_json)
    {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let params = sci_form::ani::ani_tm::default_ani_tm_params();
    let aevs = sci_form::ani::ani_tm::compute_aevs_tm(&elements, &positions, &params);
    let aev_length = aevs.first().map(|a| a.len()).unwrap_or(0);
    serialize_or_error(&serde_json::json!({
        "n_atoms": elements.len(),
        "aev_length": aev_length,
        "aevs": aevs,
    }))
}

/// **ALPHA** — Run Kohn-Sham DFT calculation.
///
/// `method`: `"svwn"` (LDA) or `"pbe"` (GGA).
/// Returns JSON DftResult with energy, orbital energies, gap, Mulliken charges.
#[cfg(feature = "alpha-dft")]
#[wasm_bindgen]
pub fn compute_dft(elements_json: &str, coords_flat_json: &str, method: &str) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements_json, coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let dft_method = match method {
        "pbe" => sci_form::dft::DftMethod::Pbe,
        _ => sci_form::dft::DftMethod::Svwn,
    };
    let config = sci_form::dft::ks_fock::DftConfig {
        method: dft_method,
        ..Default::default()
    };
    match sci_form::dft::solve_ks_dft(&elems, &positions, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}
