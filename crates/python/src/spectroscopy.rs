//! sTDA UV-Vis, vibrational analysis, and IR spectrum bindings.

use crate::system::coords_to_positions;
use pyo3::prelude::*;
use pyo3::types::PyDict;

#[pyclass(name = "StdaExcitation")]
#[derive(Clone)]
pub(crate) struct StdaExcitationPy {
    #[pyo3(get)]
    energy_ev: f64,
    #[pyo3(get)]
    wavelength_nm: f64,
    #[pyo3(get)]
    oscillator_strength: f64,
    #[pyo3(get)]
    from_mo: usize,
    #[pyo3(get)]
    to_mo: usize,
    #[pyo3(get)]
    transition_dipole: f64,
}

#[pyclass(name = "StdaUvVisSpectrum")]
#[derive(Clone)]
pub(crate) struct StdaUvVisSpectrumPy {
    #[pyo3(get)]
    energies_ev: Vec<f64>,
    #[pyo3(get)]
    wavelengths_nm: Vec<f64>,
    #[pyo3(get)]
    absorptivity: Vec<f64>,
    #[pyo3(get)]
    excitations: Vec<StdaExcitationPy>,
    #[pyo3(get)]
    sigma: f64,
    #[pyo3(get)]
    broadening: String,
    #[pyo3(get)]
    notes: Vec<String>,
}

#[pyclass(name = "VibrationalMode")]
#[derive(Clone)]
pub(crate) struct VibrationalModePy {
    #[pyo3(get)]
    frequency_cm1: f64,
    #[pyo3(get)]
    ir_intensity: f64,
    #[pyo3(get)]
    displacement: Vec<f64>,
    #[pyo3(get)]
    is_real: bool,
}

#[pyclass(name = "VibrationalAnalysis")]
#[derive(Clone)]
pub(crate) struct VibrationalAnalysisPy {
    #[pyo3(get)]
    n_atoms: usize,
    #[pyo3(get)]
    modes: Vec<VibrationalModePy>,
    #[pyo3(get)]
    n_real_modes: usize,
    #[pyo3(get)]
    zpve_ev: f64,
    #[pyo3(get)]
    method: String,
    #[pyo3(get)]
    notes: Vec<String>,
}

#[pyclass(name = "IrPeak")]
#[derive(Clone)]
pub(crate) struct IrPeakPy {
    #[pyo3(get)]
    frequency_cm1: f64,
    #[pyo3(get)]
    ir_intensity: f64,
    #[pyo3(get)]
    mode_index: usize,
}

#[pyclass(name = "IrSpectrum")]
#[derive(Clone)]
pub(crate) struct IrSpectrumPy {
    #[pyo3(get)]
    wavenumbers: Vec<f64>,
    #[pyo3(get)]
    intensities: Vec<f64>,
    #[pyo3(get)]
    peaks: Vec<IrPeakPy>,
    #[pyo3(get)]
    gamma: f64,
    #[pyo3(get)]
    notes: Vec<String>,
}

#[pyfunction]
#[pyo3(signature = (elements, coords, sigma=0.3, e_min=1.0, e_max=8.0, n_points=500, broadening="gaussian"))]
fn compute_stda_uvvis(
    elements: Vec<u8>,
    coords: Vec<f64>,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
    broadening: &str,
) -> PyResult<StdaUvVisSpectrumPy> {
    if coords.len() != elements.len() * 3 {
        return Err(crate::errors::value_error(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    let positions = coords_to_positions(&coords);
    let bt: sci_form_core::reactivity::BroadeningType =
        broadening.parse().map_err(crate::errors::value_error)?;
    sci_form_core::compute_stda_uvvis(&elements, &positions, sigma, e_min, e_max, n_points, bt)
        .map(|r| StdaUvVisSpectrumPy {
            energies_ev: r.energies_ev,
            wavelengths_nm: r.wavelengths_nm,
            absorptivity: r.absorptivity,
            excitations: r
                .excitations
                .iter()
                .map(|e| StdaExcitationPy {
                    energy_ev: e.energy_ev,
                    wavelength_nm: e.wavelength_nm,
                    oscillator_strength: e.oscillator_strength,
                    from_mo: e.from_mo,
                    to_mo: e.to_mo,
                    transition_dipole: e.transition_dipole,
                })
                .collect(),
            sigma: r.sigma,
            broadening: format!("{:?}", r.broadening),
            notes: r.notes,
        })
        .map_err(crate::errors::engine_error)
}

#[pyclass(name = "FluorescenceSpectrum")]
#[derive(Clone)]
pub(crate) struct FluorescenceSpectrumPy {
    #[pyo3(get)]
    energies_ev: Vec<f64>,
    #[pyo3(get)]
    wavelengths_nm: Vec<f64>,
    #[pyo3(get)]
    intensity: Vec<f64>,
    #[pyo3(get)]
    absorption_energy_ev: f64,
    #[pyo3(get)]
    emission_energy_ev: f64,
    #[pyo3(get)]
    emission_wavelength_nm: f64,
    #[pyo3(get)]
    oscillator_strength: f64,
    #[pyo3(get)]
    radiative_rate_hz: f64,
    #[pyo3(get)]
    radiative_lifetime_ns: f64,
    #[pyo3(get)]
    kasha: bool,
    #[pyo3(get)]
    from_mo: usize,
    #[pyo3(get)]
    to_mo: usize,
    #[pyo3(get)]
    notes: Vec<String>,
}

#[pyfunction]
#[pyo3(signature = (elements, coords, sigma=0.25, e_min=0.5, e_max=8.0, n_points=400, stokes_ev=0.0, refractive_index=1.0, dark_threshold=1e-3, broadening="gaussian", options=None))]
fn compute_fluorescence(
    elements: Vec<u8>,
    coords: Vec<f64>,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
    stokes_ev: f64,
    refractive_index: f64,
    dark_threshold: f64,
    broadening: &str,
    options: Option<&Bound<'_, PyDict>>,
) -> PyResult<FluorescenceSpectrumPy> {
    // The optional `options` mapping overrides the positional arguments
    // (blocker C-M3).
    const KEYS: &[&str] = &[
        "sigma",
        "e_min",
        "e_max",
        "n_points",
        "stokes_ev",
        "refractive_index",
        "dark_threshold",
        "broadening",
    ];
    crate::options::reject_unknown(options, KEYS)?;
    let sigma = crate::options::f64_from(options, "sigma", sigma)?;
    let e_min = crate::options::f64_from(options, "e_min", e_min)?;
    let e_max = crate::options::f64_from(options, "e_max", e_max)?;
    let n_points = crate::options::usize_from(options, "n_points", n_points)?;
    let stokes_ev = crate::options::f64_from(options, "stokes_ev", stokes_ev)?;
    let refractive_index = crate::options::f64_from(options, "refractive_index", refractive_index)?;
    let dark_threshold = crate::options::f64_from(options, "dark_threshold", dark_threshold)?;
    let broadening = &crate::options::str_from(options, "broadening", broadening)?;
    let broadening = broadening.as_str();
    if coords.len() != elements.len() * 3 {
        return Err(crate::errors::value_error(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    let positions = coords_to_positions(&coords);
    let bt: sci_form_core::reactivity::BroadeningType =
        broadening.parse().map_err(crate::errors::value_error)?;
    let config = sci_form_core::spectroscopy::FluorescenceConfig {
        sigma,
        e_min,
        e_max,
        n_points,
        stokes_ev,
        refractive_index,
        dark_threshold,
        broadening: bt,
    };
    sci_form_core::compute_fluorescence(&elements, &positions, config)
        .map(|r| FluorescenceSpectrumPy {
            energies_ev: r.energies_ev,
            wavelengths_nm: r.wavelengths_nm,
            intensity: r.intensity,
            absorption_energy_ev: r.absorption_energy_ev,
            emission_energy_ev: r.emission_energy_ev,
            emission_wavelength_nm: r.emission_wavelength_nm,
            oscillator_strength: r.oscillator_strength,
            radiative_rate_hz: r.radiative_rate_hz,
            radiative_lifetime_ns: r.radiative_lifetime_ns,
            kasha: r.kasha,
            from_mo: r.from_mo,
            to_mo: r.to_mo,
            notes: r.notes,
        })
        .map_err(crate::errors::engine_error)
}

#[pyfunction]
#[pyo3(signature = (elements, coords, method="eht", step_size=None, smiles=None))]
fn compute_vibrational_analysis(
    elements: Vec<u8>,
    coords: Vec<f64>,
    method: &str,
    step_size: Option<f64>,
    smiles: Option<&str>,
) -> PyResult<VibrationalAnalysisPy> {
    if coords.len() != elements.len() * 3 {
        return Err(crate::errors::value_error(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    let positions = coords_to_positions(&coords);
    // The method name is the caller's input: a wrong one is a ValueError,
    // not an engine failure.
    if !method.eq_ignore_ascii_case("uff") {
        method
            .parse::<sci_form_core::ir::HessianMethod>()
            .map_err(crate::errors::value_error)?;
    }

    let analysis = if method.eq_ignore_ascii_case("uff") {
        let smiles = smiles.ok_or_else(|| {
            crate::errors::value_error("UFF vibrational analysis requires smiles=...")
        })?;
        sci_form_core::compute_vibrational_analysis_uff(smiles, &elements, &positions, step_size)
    } else {
        sci_form_core::compute_vibrational_analysis(&elements, &positions, method, step_size)
    };
    analysis
        .map(|r| VibrationalAnalysisPy {
            n_atoms: r.n_atoms,
            modes: r
                .modes
                .iter()
                .map(|m| VibrationalModePy {
                    frequency_cm1: m.frequency_cm1,
                    ir_intensity: m.ir_intensity,
                    displacement: m.displacement.clone(),
                    is_real: m.is_real,
                })
                .collect(),
            n_real_modes: r.n_real_modes,
            zpve_ev: r.zpve_ev,
            method: r.method.clone(),
            notes: r.notes.clone(),
        })
        .map_err(crate::errors::engine_error)
}

#[pyfunction]
#[pyo3(signature = (elements, coords, method="eht", gamma=15.0, wn_min=400.0, wn_max=4000.0, n_points=1000, step_size=None, smiles=None, options=None))]
fn compute_ir_spectrum(
    elements: Vec<u8>,
    coords: Vec<f64>,
    method: &str,
    gamma: f64,
    wn_min: f64,
    wn_max: f64,
    n_points: usize,
    step_size: Option<f64>,
    smiles: Option<&str>,
    options: Option<&Bound<'_, PyDict>>,
) -> PyResult<IrSpectrumPy> {
    // The optional `options` mapping overrides the positional arguments
    // (blocker C-M3).
    const KEYS: &[&str] = &[
        "method",
        "gamma",
        "wn_min",
        "wn_max",
        "n_points",
        "step_size",
        "smiles",
    ];
    crate::options::reject_unknown(options, KEYS)?;
    let method = &crate::options::str_from(options, "method", method)?;
    let method = method.as_str();
    let gamma = crate::options::f64_from(options, "gamma", gamma)?;
    let wn_min = crate::options::f64_from(options, "wn_min", wn_min)?;
    let wn_max = crate::options::f64_from(options, "wn_max", wn_max)?;
    let n_points = crate::options::usize_from(options, "n_points", n_points)?;
    let step_size = crate::options::opt_f64_from(options, "step_size", step_size)?;
    if coords.len() != elements.len() * 3 {
        return Err(crate::errors::value_error(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    let positions = coords_to_positions(&coords);
    let analysis = if method.eq_ignore_ascii_case("uff") {
        let smiles = smiles
            .ok_or_else(|| crate::errors::value_error("UFF IR spectrum requires smiles=..."))?;
        sci_form_core::compute_vibrational_analysis_uff(smiles, &elements, &positions, step_size)
    } else {
        sci_form_core::compute_vibrational_analysis(&elements, &positions, method, step_size)
    }
    .map_err(crate::errors::engine_error)?;
    let spectrum = sci_form_core::compute_ir_spectrum(&analysis, gamma, wn_min, wn_max, n_points);
    Ok(IrSpectrumPy {
        wavenumbers: spectrum.wavenumbers,
        intensities: spectrum.intensities,
        peaks: spectrum
            .peaks
            .iter()
            .map(|p| IrPeakPy {
                frequency_cm1: p.frequency_cm1,
                ir_intensity: p.ir_intensity,
                mode_index: p.mode_index,
            })
            .collect(),
        gamma: spectrum.gamma,
        notes: spectrum.notes,
    })
}

#[pyclass(name = "FunctionalGroupMatch")]
#[derive(Clone)]
pub(crate) struct FunctionalGroupMatchPy {
    #[pyo3(get)]
    group: String,
    #[pyo3(get)]
    vibration_type: String,
    #[pyo3(get)]
    match_quality: f64,
}

#[pyclass(name = "PeakAssignment")]
#[derive(Clone)]
pub(crate) struct PeakAssignmentPy {
    #[pyo3(get)]
    frequency: f64,
    #[pyo3(get)]
    intensity: f64,
    #[pyo3(get)]
    assignments: Vec<FunctionalGroupMatchPy>,
    #[pyo3(get)]
    confidence: f64,
}

#[pyclass(name = "AssignmentResult")]
#[derive(Clone)]
pub(crate) struct AssignmentResultPy {
    #[pyo3(get)]
    assignments: Vec<PeakAssignmentPy>,
    #[pyo3(get)]
    functional_groups: Vec<String>,
}
#[pymethods]
impl AssignmentResultPy {
    fn __repr__(&self) -> String {
        format!(
            "AssignmentResult(n_peaks={}, groups={:?})",
            self.assignments.len(),
            self.functional_groups
        )
    }
}

#[pyfunction]
#[pyo3(signature = (elements, coords, method="eht", gamma=15.0, wn_min=400.0, wn_max=4000.0, n_points=1000, broadening="lorentzian", step_size=None, smiles=None, options=None))]
fn compute_ir_spectrum_broadened(
    elements: Vec<u8>,
    coords: Vec<f64>,
    method: &str,
    gamma: f64,
    wn_min: f64,
    wn_max: f64,
    n_points: usize,
    broadening: &str,
    step_size: Option<f64>,
    smiles: Option<&str>,
    options: Option<&Bound<'_, PyDict>>,
) -> PyResult<IrSpectrumPy> {
    // The optional `options` mapping overrides the positional arguments
    // (blocker C-M3).
    const KEYS: &[&str] = &[
        "method",
        "gamma",
        "wn_min",
        "wn_max",
        "n_points",
        "broadening",
        "step_size",
        "smiles",
    ];
    crate::options::reject_unknown(options, KEYS)?;
    let method = &crate::options::str_from(options, "method", method)?;
    let method = method.as_str();
    let gamma = crate::options::f64_from(options, "gamma", gamma)?;
    let wn_min = crate::options::f64_from(options, "wn_min", wn_min)?;
    let wn_max = crate::options::f64_from(options, "wn_max", wn_max)?;
    let n_points = crate::options::usize_from(options, "n_points", n_points)?;
    let broadening = &crate::options::str_from(options, "broadening", broadening)?;
    let broadening = broadening.as_str();
    let step_size = crate::options::opt_f64_from(options, "step_size", step_size)?;
    if coords.len() != elements.len() * 3 {
        return Err(crate::errors::value_error(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    let positions = coords_to_positions(&coords);
    let analysis = if method.eq_ignore_ascii_case("uff") {
        let smiles = smiles
            .ok_or_else(|| crate::errors::value_error("UFF IR spectrum requires smiles=..."))?;
        sci_form_core::compute_vibrational_analysis_uff(smiles, &elements, &positions, step_size)
    } else {
        sci_form_core::compute_vibrational_analysis(&elements, &positions, method, step_size)
    }
    .map_err(crate::errors::engine_error)?;
    let spectrum = sci_form_core::compute_ir_spectrum_broadened(
        &analysis, gamma, wn_min, wn_max, n_points, broadening,
    )
    .map_err(crate::errors::value_error)?;
    Ok(IrSpectrumPy {
        wavenumbers: spectrum.wavenumbers,
        intensities: spectrum.intensities,
        peaks: spectrum
            .peaks
            .iter()
            .map(|p| IrPeakPy {
                frequency_cm1: p.frequency_cm1,
                ir_intensity: p.ir_intensity,
                mode_index: p.mode_index,
            })
            .collect(),
        gamma: spectrum.gamma,
        notes: spectrum.notes,
    })
}

/// Assign IR peaks to functional groups using standard group-frequency tables.
///
/// `frequencies` and `intensities` come from a vibrational analysis;
/// `mode_displacements` (optional) is one list of per-atom xyz displacements
/// per mode (each mode a list of `[x, y, z]` triples) and refines the
/// assignment with the normal-mode character.
#[pyfunction]
#[pyo3(signature = (frequencies, intensities, elements, mode_displacements=None))]
fn assign_peaks(
    frequencies: Vec<f64>,
    intensities: Vec<f64>,
    elements: Vec<u8>,
    mode_displacements: Option<Vec<Vec<Vec<f64>>>>,
) -> AssignmentResultPy {
    let displacements: Option<Vec<Vec<[f64; 3]>>> = mode_displacements.map(|modes| {
        modes
            .into_iter()
            .map(|mode| {
                mode.into_iter()
                    .map(|xyz| {
                        if xyz.len() != 3 {
                            [0.0, 0.0, 0.0]
                        } else {
                            [xyz[0], xyz[1], xyz[2]]
                        }
                    })
                    .collect::<Vec<[f64; 3]>>()
            })
            .collect()
    });
    let result = sci_form_core::ir::peak_assignment::assign_peaks(
        &frequencies,
        &intensities,
        &elements,
        displacements.as_deref(),
    );
    AssignmentResultPy {
        assignments: result
            .assignments
            .iter()
            .map(|a| PeakAssignmentPy {
                frequency: a.frequency,
                intensity: a.intensity,
                assignments: a
                    .assignments
                    .iter()
                    .map(|m| FunctionalGroupMatchPy {
                        group: m.group.clone(),
                        vibration_type: m.vibration_type.clone(),
                        match_quality: m.match_quality,
                    })
                    .collect(),
                confidence: a.confidence,
            })
            .collect(),
        functional_groups: result.functional_groups,
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_stda_uvvis, m)?)?;
    m.add_function(wrap_pyfunction!(compute_fluorescence, m)?)?;
    m.add_function(wrap_pyfunction!(compute_vibrational_analysis, m)?)?;
    m.add_function(wrap_pyfunction!(compute_ir_spectrum, m)?)?;
    m.add_function(wrap_pyfunction!(compute_ir_spectrum_broadened, m)?)?;
    m.add_function(wrap_pyfunction!(assign_peaks, m)?)?;
    m.add_class::<StdaExcitationPy>()?;
    m.add_class::<StdaUvVisSpectrumPy>()?;
    m.add_class::<FluorescenceSpectrumPy>()?;
    m.add_class::<VibrationalModePy>()?;
    m.add_class::<VibrationalAnalysisPy>()?;
    m.add_class::<IrPeakPy>()?;
    m.add_class::<IrSpectrumPy>()?;
    m.add_class::<FunctionalGroupMatchPy>()?;
    m.add_class::<PeakAssignmentPy>()?;
    m.add_class::<AssignmentResultPy>()?;
    Ok(())
}
