//! Python bindings for the v2.x feature set.
//!
//! Every function returns a JSON string so that the nested result types keep
//! their full structure without duplicating dozens of Python classes. The
//! docstrings name the Rust type that the JSON encodes. Configuration payloads
//! are parsed from the same JSON form used by the CLI and WASM bindings.

use pyo3::prelude::*;

use crate::errors::{engine_error, value_error};
use crate::system::coords_to_positions;

fn to_json<T: serde::Serialize>(value: &T) -> PyResult<String> {
    serde_json::to_string(value)
        .map_err(|error| value_error(format!("could not serialize the result: {error}")))
}

fn parse_json<T: serde::de::DeserializeOwned>(input: &str, what: &str) -> PyResult<T> {
    serde_json::from_str(input)
        .map_err(|error| value_error(format!("invalid {what} JSON: {error}")))
}

/// Isotropic EPR hyperfine coupling from the UHF spin density.
///
/// `elements` and `coords` describe the molecule. `nucleus=None` evaluates
/// every atom with a known nuclear g-factor. Returns the JSON encoding of
/// `sci_form_core::epr::EprHyperfineResult`. Coupling is in MHz.
#[pyfunction]
#[pyo3(signature = (elements, coords, nucleus=None, charge=0, multiplicity=2))]
fn compute_epr_hyperfine(
    elements: Vec<u8>,
    coords: Vec<f64>,
    nucleus: Option<u8>,
    charge: i32,
    multiplicity: u32,
) -> PyResult<String> {
    crate::errors::check_coords_len(&coords, elements.len())?;
    let config = sci_form_core::epr::EprConfig {
        nucleus,
        charge,
        multiplicity,
    };
    let result =
        sci_form_core::compute_epr_hyperfine(&elements, &coords_to_positions(&coords), &config)
            .map_err(engine_error)?;
    to_json(&result)
}

/// Conventional closed-shell MP2 correlation energy on RHF/STO-3G orbitals.
///
/// Returns the JSON encoding of `sci_form_core::mp2::Mp2Result`. Energies are
/// in Hartree and kcal/mol.
#[pyfunction]
fn compute_mp2_energy(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<String> {
    crate::errors::check_coords_len(&coords, elements.len())?;
    let result = sci_form_core::compute_mp2_energy(&elements, &coords_to_positions(&coords))
        .map_err(engine_error)?;
    to_json(&result)
}

/// Koopmans core-level XPS binding energies.
///
/// `elements_filter` restricts the output to those atomic numbers. Returns the
/// JSON encoding of `sci_form_core::xps::XpsResult`. Binding energies are in eV.
#[pyfunction]
#[pyo3(signature = (elements, coords, elements_filter=None, charge=0, multiplicity=1))]
fn compute_xps_binding_energies(
    elements: Vec<u8>,
    coords: Vec<f64>,
    elements_filter: Option<Vec<u8>>,
    charge: i32,
    multiplicity: u32,
) -> PyResult<String> {
    crate::errors::check_coords_len(&coords, elements.len())?;
    let config = sci_form_core::xps::XpsConfig {
        charge,
        multiplicity,
        elements: elements_filter,
    };
    let result = sci_form_core::compute_xps_binding_energies(
        &elements,
        &coords_to_positions(&coords),
        &config,
    )
    .map_err(engine_error)?;
    to_json(&result)
}

/// Raman activities from the response polarizability and UFF normal modes.
///
/// Returns the JSON encoding of `sci_form_core::raman::RamanSpectrum`.
/// Frequencies are cm⁻¹ and the tensor is in Å³.
#[pyfunction]
#[pyo3(signature = (smiles, elements, coords, step_angstrom=0.01))]
fn compute_raman_spectrum(
    smiles: &str,
    elements: Vec<u8>,
    coords: Vec<f64>,
    step_angstrom: f64,
) -> PyResult<String> {
    crate::errors::check_coords_len(&coords, elements.len())?;
    let config = sci_form_core::raman::RamanConfig { step_angstrom };
    let result = sci_form_core::compute_raman_spectrum(
        smiles,
        &elements,
        &coords_to_positions(&coords),
        &config,
    )
    .map_err(engine_error)?;
    to_json(&result)
}

/// VCD rotational strengths with the nuclear charge-current screening model.
///
/// Returns the JSON encoding of `sci_form_core::vcd_roa::VcdSpectrum`.
#[pyfunction]
#[pyo3(signature = (smiles, elements, coords, min_frequency_cm1=1.0))]
fn compute_vcd_spectrum(
    smiles: &str,
    elements: Vec<u8>,
    coords: Vec<f64>,
    min_frequency_cm1: f64,
) -> PyResult<String> {
    crate::errors::check_coords_len(&coords, elements.len())?;
    let config = sci_form_core::vcd_roa::VcdRoaConfig {
        raman: sci_form_core::raman::RamanConfig::default(),
        min_frequency_cm1,
    };
    let result = sci_form_core::compute_vcd_spectrum(
        smiles,
        &elements,
        &coords_to_positions(&coords),
        &config,
    )
    .map_err(engine_error)?;
    to_json(&result)
}

/// ROA screening intensities from Raman activities and the VCD magnetic term.
///
/// Returns the JSON encoding of `sci_form_core::vcd_roa::RoaSpectrum`.
#[pyfunction]
#[pyo3(signature = (smiles, elements, coords, min_frequency_cm1=1.0))]
fn compute_roa_spectrum(
    smiles: &str,
    elements: Vec<u8>,
    coords: Vec<f64>,
    min_frequency_cm1: f64,
) -> PyResult<String> {
    crate::errors::check_coords_len(&coords, elements.len())?;
    let config = sci_form_core::vcd_roa::VcdRoaConfig {
        raman: sci_form_core::raman::RamanConfig::default(),
        min_frequency_cm1,
    };
    let result = sci_form_core::compute_roa_spectrum(
        smiles,
        &elements,
        &coords_to_positions(&coords),
        &config,
    )
    .map_err(engine_error)?;
    to_json(&result)
}

/// Mechanical-embedding QM/MM energy.
///
/// `system_json` is the serialized `QmMmSystem`; `config_json` is the
/// serialized `QmMmConfig`. Returns `QmMmResult` as JSON. Energies are kcal/mol.
#[pyfunction]
fn compute_qmmm_energy(system_json: &str, config_json: &str) -> PyResult<String> {
    let system: sci_form_core::qmmm::QmMmSystem = parse_json(system_json, "QM/MM system")?;
    let config: sci_form_core::qmmm::QmMmConfig = parse_json(config_json, "QM/MM config")?;
    let result = sci_form_core::compute_qmmm_energy(&system, &config).map_err(engine_error)?;
    to_json(&result)
}

/// GCMC adsorption isotherm on a periodic framework.
///
/// The three JSON payloads are the serialized `GcmcFramework`, `GuestModel` and
/// `GcmcConfig`. Returns `AdsorptionIsotherm` as JSON.
#[pyfunction]
fn compute_gcmc_isotherm(
    framework_json: &str,
    guest_json: &str,
    config_json: &str,
) -> PyResult<String> {
    let framework: sci_form_core::gcmc::GcmcFramework =
        parse_json(framework_json, "GCMC framework")?;
    let guest: sci_form_core::gcmc::GuestModel = parse_json(guest_json, "GCMC guest")?;
    let config: sci_form_core::gcmc::GcmcConfig = parse_json(config_json, "GCMC config")?;
    let result =
        sci_form_core::compute_gcmc_isotherm(&framework, &guest, &config).map_err(engine_error)?;
    to_json(&result)
}

/// Rule-based reaction-network discovery.
///
/// `config_json` is the serialized `TemplateNetworkConfig`. Returns
/// `ReactionNetwork` as JSON.
#[pyfunction]
fn discover_reaction_network(config_json: &str) -> PyResult<String> {
    let config: sci_form_core::networks::TemplateNetworkConfig =
        parse_json(config_json, "reaction-network config")?;
    let result = config.discover().map_err(engine_error)?;
    to_json(&result)
}

/// Gillespie kinetic Monte Carlo on an elementary-step network.
///
/// `config_json` is the serialized `KmcConfig`. Returns `KmcResult` as JSON.
#[pyfunction]
fn run_kinetic_monte_carlo(config_json: &str) -> PyResult<String> {
    let config: sci_form_core::kmc::KmcConfig = parse_json(config_json, "kMC config")?;
    let result = sci_form_core::run_kinetic_monte_carlo(&config).map_err(engine_error)?;
    to_json(&result)
}

/// Refine a barrier with the existing GSM backend.
///
/// Returns `BarrierRefinementResult` as JSON; the barrier is in eV and
/// kcal/mol.
#[pyfunction]
#[pyo3(signature = (smiles, reactant_coords, product_coords, method="uff", config_json=None))]
fn refine_barrier_with_gsm(
    smiles: &str,
    reactant_coords: Vec<f64>,
    product_coords: Vec<f64>,
    method: &str,
    config_json: Option<&str>,
) -> PyResult<String> {
    let config: sci_form_core::gsm::GsmConfig = match config_json {
        Some(value) => parse_json(value, "GSM config")?,
        None => sci_form_core::gsm::GsmConfig::default(),
    };
    let result = sci_form_core::refine_barrier_with_gsm(
        smiles,
        &reactant_coords,
        &product_coords,
        method,
        &config,
    )
    .map_err(engine_error)?;
    to_json(&result)
}

/// Well-tempered metadynamics on a UFF molecule.
///
/// `config_json` is the serialized `MetaDynamicsConfig`. Returns
/// `MetaDynamicsResult` as JSON, including the bias-derived and reweighted
/// free-energy surfaces.
#[pyfunction]
fn run_uff_metadynamics(smiles: &str, coords: Vec<f64>, config_json: &str) -> PyResult<String> {
    let config: sci_form_core::enhanced_sampling::MetaDynamicsConfig =
        parse_json(config_json, "metadynamics config")?;
    let result =
        sci_form_core::run_uff_metadynamics(smiles, &coords, &config).map_err(engine_error)?;
    to_json(&result)
}

/// OPES-style enhanced sampling on a UFF molecule.
///
/// `config_json` is the serialized `OpesConfig`. Returns `OpesResult` as JSON.
#[pyfunction]
fn run_uff_opes(smiles: &str, coords: Vec<f64>, config_json: &str) -> PyResult<String> {
    let config: sci_form_core::enhanced_sampling::OpesConfig =
        parse_json(config_json, "OPES config")?;
    let result = sci_form_core::run_uff_opes(smiles, &coords, &config).map_err(engine_error)?;
    to_json(&result)
}

/// Reweight a biased 1D CV trace into an unbiased free-energy histogram.
///
/// Returns `FreeEnergySurface` as JSON. Energies are kcal/mol.
#[pyfunction]
fn reweight_free_energy(trace_json: &str, config_json: &str) -> PyResult<String> {
    let trace: Vec<sci_form_core::enhanced_sampling::MetaDynamicsPoint> =
        parse_json(trace_json, "CV trace")?;
    let config: sci_form_core::enhanced_sampling::ReweightConfig =
        parse_json(config_json, "reweighting config")?;
    let result = sci_form_core::enhanced_sampling::reweight_free_energy(&trace, &config)
        .map_err(engine_error)?;
    to_json(&result)
}

/// Block-averaged reweighting of a correlated 1D CV trace.
///
/// Returns `ReweightResult` as JSON with the surface and per-bin errors.
#[pyfunction]
fn reweight_free_energy_correlated(trace_json: &str, config_json: &str) -> PyResult<String> {
    let trace: Vec<sci_form_core::enhanced_sampling::MetaDynamicsPoint> =
        parse_json(trace_json, "CV trace")?;
    let config: sci_form_core::enhanced_sampling::ReweightConfig =
        parse_json(config_json, "reweighting config")?;
    let result = sci_form_core::enhanced_sampling::reweight_free_energy_correlated(&trace, &config)
        .map_err(engine_error)?;
    to_json(&result)
}

/// Reweight a 2D OPES trace into an unbiased surface.
///
/// Returns `Surface2D` as JSON. The trace is a serialized list of
/// `MultiTracePoint`.
#[pyfunction]
fn reweight_free_energy_2d(trace_json: &str, config_json: &str) -> PyResult<String> {
    let trace: Vec<sci_form_core::enhanced_sampling_multi::MultiTracePoint> =
        parse_json(trace_json, "2D trace")?;
    let config: sci_form_core::enhanced_sampling_multi::MultiReweightConfig =
        parse_json(config_json, "2D reweighting config")?;
    let result = sci_form_core::enhanced_sampling_multi::reweight_free_energy_2d(&trace, &config)
        .map_err(engine_error)?;
    to_json(&result)
}

/// Solve the 3D-RISM Ornstein–Zernike equations on a periodic grid.
///
/// `beta_u` and `susceptibility` are row-major `nx × ny × nz` arrays;
/// `config_json` is the serialized `RismConfig`. Returns `RismResult` as JSON.
#[pyfunction]
fn solve_3d_rism(
    beta_u: Vec<f64>,
    susceptibility: Vec<f64>,
    config_json: &str,
) -> PyResult<String> {
    let config: sci_form_core::rism::RismConfig = parse_json(config_json, "RISM config")?;
    let result =
        sci_form_core::solve_3d_rism(&beta_u, &susceptibility, &config).map_err(engine_error)?;
    to_json(&result)
}

/// Estimate an alchemical free-energy difference.
///
/// `samples_json` is a list of `AlchemicalSample`; `estimator` is `"fep"`,
/// `"bar"` or `"mbar"`; `config_json` is a `FreeEnergyConfig`. Returns
/// `FreeEnergyResult` as JSON.
#[pyfunction]
#[pyo3(signature = (samples_json, estimator="mbar", config_json=None))]
fn estimate_free_energy(
    samples_json: &str,
    estimator: &str,
    config_json: Option<&str>,
) -> PyResult<String> {
    let samples: Vec<sci_form_core::free_energy::AlchemicalSample> =
        parse_json(samples_json, "alchemical samples")?;
    let estimator = match estimator.to_ascii_lowercase().as_str() {
        "fep" => sci_form_core::free_energy::FreeEnergyEstimator::Fep,
        "bar" => sci_form_core::free_energy::FreeEnergyEstimator::Bar,
        "mbar" => sci_form_core::free_energy::FreeEnergyEstimator::Mbar,
        other => {
            return Err(value_error(crate::errors::invalid_value(
                "estimator",
                other,
                &["fep", "bar", "mbar"],
            )))
        }
    };
    let config: sci_form_core::free_energy::FreeEnergyConfig = match config_json {
        Some(value) => parse_json(value, "free-energy config")?,
        None => sci_form_core::free_energy::FreeEnergyConfig::default(),
    };
    let result = sci_form_core::free_energy::estimate_free_energy(&samples, estimator, &config)
        .map_err(engine_error)?;
    to_json(&result)
}

/// Fit a split-conformal calibration radius.
///
/// Returns `ConformalCalibrator` as JSON. The radius is in the target unit and
/// the coverage is the requested level.
#[pyfunction]
fn conformal_calibrate(
    predictions: Vec<f64>,
    targets: Vec<f64>,
    coverage: f64,
) -> PyResult<String> {
    let calibrator =
        sci_form_core::active_learning::ConformalCalibrator::fit(&predictions, &targets, coverage)
            .map_err(engine_error)?;
    to_json(&calibrator)
}

/// List every registered NMR-active isotope for an element.
///
/// Returns a JSON list of `NmrIsotopeInfo`.
#[pyfunction]
fn active_nmr_isotopes(element: u8) -> PyResult<String> {
    let isotopes = sci_form_core::nmr_all::active_isotopes(element);
    to_json(&isotopes)
}

/// Run the all-nuclide NMR driver for a molecule.
///
/// `config_json` is the serialized `GiaoNmrConfig`. Returns `AllNmrResult` as
/// JSON; unsupported elements appear in `failures`.
#[pyfunction]
fn compute_all_nmr(elements: Vec<u8>, coords: Vec<f64>, config_json: &str) -> PyResult<String> {
    crate::errors::check_coords_len(&coords, elements.len())?;
    let config: sci_form_core::GiaoNmrConfig = parse_json(config_json, "NMR config")?;
    let result =
        sci_form_core::nmr_all::compute_all_nmr(&elements, &coords_to_positions(&coords), &config)
            .map_err(engine_error)?;
    to_json(&result)
}

/// Restricted PBE0 energy (requires the `alpha-dft` feature).
///
/// Returns the JSON encoding of `sci_form_core::dft::ks_fock::DftResult`.
#[cfg(feature = "alpha-dft")]
#[pyfunction]
#[pyo3(signature = (elements, coords, max_iterations=100))]
fn compute_pbe0_energy(
    elements: Vec<u8>,
    coords: Vec<f64>,
    max_iterations: usize,
) -> PyResult<String> {
    crate::errors::check_coords_len(&coords, elements.len())?;
    let config = sci_form_core::dft::ks_fock::DftConfig {
        method: sci_form_core::dft::ks_fock::DftMethod::Pbe0,
        max_iterations,
        ..Default::default()
    };
    let result = sci_form_core::dft::ks_fock::solve_ks_dft(
        &elements,
        &coords_to_positions(&coords),
        &config,
    )
    .map_err(engine_error)?;
    to_json(&result)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_epr_hyperfine, m)?)?;
    m.add_function(wrap_pyfunction!(compute_mp2_energy, m)?)?;
    m.add_function(wrap_pyfunction!(compute_xps_binding_energies, m)?)?;
    m.add_function(wrap_pyfunction!(compute_raman_spectrum, m)?)?;
    m.add_function(wrap_pyfunction!(compute_vcd_spectrum, m)?)?;
    m.add_function(wrap_pyfunction!(compute_roa_spectrum, m)?)?;
    m.add_function(wrap_pyfunction!(compute_qmmm_energy, m)?)?;
    m.add_function(wrap_pyfunction!(compute_gcmc_isotherm, m)?)?;
    m.add_function(wrap_pyfunction!(discover_reaction_network, m)?)?;
    m.add_function(wrap_pyfunction!(run_kinetic_monte_carlo, m)?)?;
    m.add_function(wrap_pyfunction!(refine_barrier_with_gsm, m)?)?;
    m.add_function(wrap_pyfunction!(run_uff_metadynamics, m)?)?;
    m.add_function(wrap_pyfunction!(run_uff_opes, m)?)?;
    m.add_function(wrap_pyfunction!(reweight_free_energy, m)?)?;
    m.add_function(wrap_pyfunction!(reweight_free_energy_correlated, m)?)?;
    m.add_function(wrap_pyfunction!(reweight_free_energy_2d, m)?)?;
    m.add_function(wrap_pyfunction!(solve_3d_rism, m)?)?;
    m.add_function(wrap_pyfunction!(estimate_free_energy, m)?)?;
    m.add_function(wrap_pyfunction!(conformal_calibrate, m)?)?;
    m.add_function(wrap_pyfunction!(active_nmr_isotopes, m)?)?;
    m.add_function(wrap_pyfunction!(compute_all_nmr, m)?)?;
    #[cfg(feature = "alpha-dft")]
    m.add_function(wrap_pyfunction!(compute_pbe0_energy, m)?)?;
    Ok(())
}
