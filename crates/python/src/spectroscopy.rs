//! sTDA UV-Vis, vibrational analysis, and IR spectrum bindings.

use crate::system::coords_to_positions;
use pyo3::prelude::*;

#[pyclass]
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

#[pyclass]
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

#[pyclass]
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

#[pyclass]
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

#[pyclass]
#[derive(Clone)]
pub(crate) struct IrPeakPy {
    #[pyo3(get)]
    frequency_cm1: f64,
    #[pyo3(get)]
    ir_intensity: f64,
    #[pyo3(get)]
    mode_index: usize,
}

#[pyclass]
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
fn stda_uvvis(
    elements: Vec<u8>,
    coords: Vec<f64>,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
    broadening: &str,
) -> PyResult<StdaUvVisSpectrumPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    let positions = coords_to_positions(&coords);
    let bt = match broadening {
        "lorentzian" | "Lorentzian" => sci_form_core::reactivity::BroadeningType::Lorentzian,
        _ => sci_form_core::reactivity::BroadeningType::Gaussian,
    };
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
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
#[pyo3(signature = (elements, coords, method="eht", step_size=None, smiles=None))]
fn vibrational_analysis(
    elements: Vec<u8>,
    coords: Vec<f64>,
    method: &str,
    step_size: Option<f64>,
    smiles: Option<&str>,
) -> PyResult<VibrationalAnalysisPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    let positions = coords_to_positions(&coords);

    let analysis = if method.eq_ignore_ascii_case("uff") {
        let smiles = smiles.ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err("UFF vibrational analysis requires smiles=...")
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
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
#[pyo3(signature = (elements, coords, method="eht", gamma=15.0, wn_min=400.0, wn_max=4000.0, n_points=1000, step_size=None, smiles=None))]
fn ir_spectrum(
    elements: Vec<u8>,
    coords: Vec<f64>,
    method: &str,
    gamma: f64,
    wn_min: f64,
    wn_max: f64,
    n_points: usize,
    step_size: Option<f64>,
    smiles: Option<&str>,
) -> PyResult<IrSpectrumPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    let positions = coords_to_positions(&coords);
    let analysis = if method.eq_ignore_ascii_case("uff") {
        let smiles = smiles.ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err("UFF IR spectrum requires smiles=...")
        })?;
        sci_form_core::compute_vibrational_analysis_uff(smiles, &elements, &positions, step_size)
    } else {
        sci_form_core::compute_vibrational_analysis(&elements, &positions, method, step_size)
    }
    .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
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

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(stda_uvvis, m)?)?;
    m.add_function(wrap_pyfunction!(vibrational_analysis, m)?)?;
    m.add_function(wrap_pyfunction!(ir_spectrum, m)?)?;
    m.add_class::<StdaExcitationPy>()?;
    m.add_class::<StdaUvVisSpectrumPy>()?;
    m.add_class::<VibrationalModePy>()?;
    m.add_class::<VibrationalAnalysisPy>()?;
    m.add_class::<IrPeakPy>()?;
    m.add_class::<IrSpectrumPy>()?;
    Ok(())
}
