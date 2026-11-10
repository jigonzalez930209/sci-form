//! WASM bindings for EPR, MP2 and XPS.

use crate::helpers::{bad_json, json_error, parse_elements_and_positions, serialize_or_error};
use wasm_bindgen::prelude::*;

/// Isotropic EPR hyperfine coupling. `nucleus = 0` evaluates all known nuclei.
#[wasm_bindgen]
pub fn compute_epr_hyperfine(
    elements: &str,
    coords_flat: &str,
    nucleus: u8,
    charge: i32,
    multiplicity: u32,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(value) => value,
        Err(error) => return json_error(&error),
    };
    let config = sci_form::epr::EprConfig {
        nucleus: if nucleus == 0 { None } else { Some(nucleus) },
        charge,
        multiplicity,
    };
    match sci_form::compute_epr_hyperfine(&elems, &positions, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// Conventional closed-shell MP2 correlation energy.
#[wasm_bindgen]
pub fn compute_mp2_energy(elements: &str, coords_flat: &str) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(value) => value,
        Err(error) => return json_error(&error),
    };
    match sci_form::compute_mp2_energy(&elems, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// Koopmans core-level XPS binding energies.
#[wasm_bindgen]
pub fn compute_xps_binding_energies(
    elements: &str,
    coords_flat: &str,
    elements_filter_json: &str,
    charge: i32,
    multiplicity: u32,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(value) => value,
        Err(error) => return json_error(&error),
    };
    let filter = if elements_filter_json.trim().is_empty() {
        None
    } else {
        match serde_json::from_str::<Vec<u8>>(elements_filter_json) {
            Ok(value) => Some(value),
            Err(error) => {
                return json_error(&bad_json(
                    "'elements_filter_json'",
                    elements_filter_json,
                    &error,
                ))
            }
        }
    };
    let config = sci_form::xps::XpsConfig {
        charge,
        multiplicity,
        elements: filter,
    };
    match sci_form::compute_xps_binding_energies(&elems, &positions, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// GCMC adsorption isotherm from JSON framework, guest and config payloads.
#[wasm_bindgen]
pub fn compute_gcmc_isotherm(framework_json: &str, guest_json: &str, config_json: &str) -> String {
    let framework: sci_form::gcmc::GcmcFramework = match serde_json::from_str(framework_json) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'framework_json'", framework_json, &error)),
    };
    let guest: sci_form::gcmc::GuestModel = match serde_json::from_str(guest_json) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'guest_json'", guest_json, &error)),
    };
    let config: sci_form::gcmc::GcmcConfig = match serde_json::from_str(config_json) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'config_json'", config_json, &error)),
    };
    match sci_form::compute_gcmc_isotherm(&framework, &guest, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// Rule-based reaction-network discovery from a JSON TemplateNetworkConfig.
#[wasm_bindgen]
pub fn discover_reaction_network(config_json: &str) -> String {
    let config: sci_form::networks::TemplateNetworkConfig = match serde_json::from_str(config_json)
    {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'config_json'", config_json, &error)),
    };
    match config.discover() {
        Ok(network) => serialize_or_error(&network),
        Err(error) => json_error(&error),
    }
}

/// Raman activities from the response polarizability and UFF modes.
#[wasm_bindgen]
pub fn compute_raman_spectrum(
    smiles: &str,
    elements: &str,
    coords_flat: &str,
    step_angstrom: f64,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(value) => value,
        Err(error) => return json_error(&error),
    };
    let config = sci_form::raman::RamanConfig { step_angstrom };
    match sci_form::compute_raman_spectrum(smiles, &elems, &positions, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// VCD rotational strengths with the nuclear charge-current screening model.
#[wasm_bindgen]
pub fn compute_vcd_spectrum(
    smiles: &str,
    elements: &str,
    coords_flat: &str,
    min_frequency_cm1: f64,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(value) => value,
        Err(error) => return json_error(&error),
    };
    let config = sci_form::vcd_roa::VcdRoaConfig {
        raman: sci_form::raman::RamanConfig::default(),
        min_frequency_cm1,
    };
    match sci_form::compute_vcd_spectrum(smiles, &elems, &positions, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// ROA screening intensities from Raman activities and the VCD magnetic term.
#[wasm_bindgen]
pub fn compute_roa_spectrum(
    smiles: &str,
    elements: &str,
    coords_flat: &str,
    min_frequency_cm1: f64,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(value) => value,
        Err(error) => return json_error(&error),
    };
    let config = sci_form::vcd_roa::VcdRoaConfig {
        raman: sci_form::raman::RamanConfig::default(),
        min_frequency_cm1,
    };
    match sci_form::compute_roa_spectrum(smiles, &elems, &positions, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// Mechanical-embedding QM/MM energy from serialized system/config payloads.
#[wasm_bindgen]
pub fn compute_qmmm_energy(system_json: &str, config_json: &str) -> String {
    let system: sci_form::qmmm::QmMmSystem = match serde_json::from_str(system_json) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'system_json'", system_json, &error)),
    };
    let config: sci_form::qmmm::QmMmConfig = match serde_json::from_str(config_json) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'config_json'", config_json, &error)),
    };
    match sci_form::compute_qmmm_energy(&system, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// Gillespie kinetic Monte Carlo from a JSON KmcConfig.
#[wasm_bindgen]
pub fn run_kinetic_monte_carlo(config_json: &str) -> String {
    let config: sci_form::kmc::KmcConfig = match serde_json::from_str(config_json) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'config_json'", config_json, &error)),
    };
    match sci_form::run_kinetic_monte_carlo(&config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// Refine a barrier with GSM from flat coordinate arrays.
#[wasm_bindgen]
pub fn refine_barrier_with_gsm(
    smiles: &str,
    reactant_coords: &str,
    product_coords: &str,
    method: &str,
    config_json: &str,
) -> String {
    let reactant: Vec<f64> = match serde_json::from_str(reactant_coords) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'reactant_coords'", reactant_coords, &error)),
    };
    let product: Vec<f64> = match serde_json::from_str(product_coords) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'product_coords'", product_coords, &error)),
    };
    let config: sci_form::gsm::GsmConfig = if config_json.trim().is_empty() {
        sci_form::gsm::GsmConfig::default()
    } else {
        match serde_json::from_str(config_json) {
            Ok(value) => value,
            Err(error) => return json_error(&bad_json("'config_json'", config_json, &error)),
        }
    };
    match sci_form::refine_barrier_with_gsm(smiles, &reactant, &product, method, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// Well-tempered metadynamics on a UFF molecule from a JSON config.
#[wasm_bindgen]
pub fn run_uff_metadynamics(smiles: &str, coords_flat: &str, config_json: &str) -> String {
    let coords: Vec<f64> = match serde_json::from_str(coords_flat) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'coords_flat'", coords_flat, &error)),
    };
    let config: sci_form::enhanced_sampling::MetaDynamicsConfig =
        match serde_json::from_str(config_json) {
            Ok(value) => value,
            Err(error) => return json_error(&bad_json("'config_json'", config_json, &error)),
        };
    match sci_form::run_uff_metadynamics(smiles, &coords, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// OPES-style enhanced sampling on a UFF molecule from a JSON config.
#[wasm_bindgen]
pub fn run_uff_opes(smiles: &str, coords_flat: &str, config_json: &str) -> String {
    let coords: Vec<f64> = match serde_json::from_str(coords_flat) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'coords_flat'", coords_flat, &error)),
    };
    let config: sci_form::enhanced_sampling::OpesConfig = match serde_json::from_str(config_json) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'config_json'", config_json, &error)),
    };
    match sci_form::run_uff_opes(smiles, &coords, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// Reweight a biased 1D CV trace from JSON.
#[wasm_bindgen]
pub fn reweight_free_energy(trace_json: &str, config_json: &str) -> String {
    let trace: Vec<sci_form::enhanced_sampling::MetaDynamicsPoint> =
        match serde_json::from_str(trace_json) {
            Ok(value) => value,
            Err(error) => return json_error(&bad_json("'trace_json'", trace_json, &error)),
        };
    let config: sci_form::enhanced_sampling::ReweightConfig =
        match serde_json::from_str(config_json) {
            Ok(value) => value,
            Err(error) => return json_error(&bad_json("'config_json'", config_json, &error)),
        };
    match sci_form::enhanced_sampling::reweight_free_energy(&trace, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// Block-averaged reweighting of a correlated 1D CV trace from JSON.
#[wasm_bindgen]
pub fn reweight_free_energy_correlated(trace_json: &str, config_json: &str) -> String {
    let trace: Vec<sci_form::enhanced_sampling::MetaDynamicsPoint> =
        match serde_json::from_str(trace_json) {
            Ok(value) => value,
            Err(error) => return json_error(&bad_json("'trace_json'", trace_json, &error)),
        };
    let config: sci_form::enhanced_sampling::ReweightConfig =
        match serde_json::from_str(config_json) {
            Ok(value) => value,
            Err(error) => return json_error(&bad_json("'config_json'", config_json, &error)),
        };
    match sci_form::enhanced_sampling::reweight_free_energy_correlated(&trace, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// Reweight a 2D OPES trace from JSON.
#[wasm_bindgen]
pub fn reweight_free_energy_2d(trace_json: &str, config_json: &str) -> String {
    let trace: Vec<sci_form::enhanced_sampling_multi::MultiTracePoint> =
        match serde_json::from_str(trace_json) {
            Ok(value) => value,
            Err(error) => return json_error(&bad_json("'trace_json'", trace_json, &error)),
        };
    let config: sci_form::enhanced_sampling_multi::MultiReweightConfig =
        match serde_json::from_str(config_json) {
            Ok(value) => value,
            Err(error) => return json_error(&bad_json("'config_json'", config_json, &error)),
        };
    match sci_form::enhanced_sampling_multi::reweight_free_energy_2d(&trace, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// Solve 3D-RISM from grid arrays and a JSON config.
#[wasm_bindgen]
pub fn solve_3d_rism(beta_u: &str, susceptibility: &str, config_json: &str) -> String {
    let beta_u: Vec<f64> = match serde_json::from_str(beta_u) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'beta_u'", beta_u, &error)),
    };
    let susceptibility: Vec<f64> = match serde_json::from_str(susceptibility) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'susceptibility'", susceptibility, &error)),
    };
    let config: sci_form::rism::RismConfig = match serde_json::from_str(config_json) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'config_json'", config_json, &error)),
    };
    match sci_form::solve_3d_rism(&beta_u, &susceptibility, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// Estimate an alchemical free-energy difference from JSON samples.
#[wasm_bindgen]
pub fn estimate_free_energy(samples_json: &str, estimator: &str, config_json: &str) -> String {
    let samples: Vec<sci_form::free_energy::AlchemicalSample> =
        match serde_json::from_str(samples_json) {
            Ok(value) => value,
            Err(error) => return json_error(&bad_json("'samples_json'", samples_json, &error)),
        };
    let estimator = match estimator.to_ascii_lowercase().as_str() {
        "fep" => sci_form::free_energy::FreeEnergyEstimator::Fep,
        "bar" => sci_form::free_energy::FreeEnergyEstimator::Bar,
        "mbar" => sci_form::free_energy::FreeEnergyEstimator::Mbar,
        other => {
            return json_error(&format!(
                "unknown estimator '{other}'; expected 'fep', 'bar' or 'mbar'"
            ))
        }
    };
    let config: sci_form::free_energy::FreeEnergyConfig = match serde_json::from_str(config_json) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'config_json'", config_json, &error)),
    };
    match sci_form::free_energy::estimate_free_energy(&samples, estimator, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// Fit a split-conformal calibration radius.
#[wasm_bindgen]
pub fn conformal_calibrate(predictions_json: &str, targets_json: &str, coverage: f64) -> String {
    let predictions: Vec<f64> = match serde_json::from_str(predictions_json) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'predictions_json'", predictions_json, &error)),
    };
    let targets: Vec<f64> = match serde_json::from_str(targets_json) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'targets_json'", targets_json, &error)),
    };
    match sci_form::active_learning::ConformalCalibrator::fit(&predictions, &targets, coverage) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// List every NMR-active isotope registered for an element.
#[wasm_bindgen]
pub fn active_nmr_isotopes(element: u8) -> String {
    serialize_or_error(&sci_form::nmr_all::active_isotopes(element))
}

/// Run the all-nuclide NMR driver from a JSON GiaoNmrConfig.
#[wasm_bindgen]
pub fn compute_all_nmr(elements: &str, coords_flat: &str, config_json: &str) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(value) => value,
        Err(error) => return json_error(&error),
    };
    let config: sci_form::GiaoNmrConfig = match serde_json::from_str(config_json) {
        Ok(value) => value,
        Err(error) => return json_error(&bad_json("'config_json'", config_json, &error)),
    };
    match sci_form::nmr_all::compute_all_nmr(&elems, &positions, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}
