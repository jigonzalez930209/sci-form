//! NMR: chemical shifts, J-couplings, spectra, HOSE codes.

use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
pub(crate) struct ChemicalShiftPy {
    #[pyo3(get)]
    atom_index: usize,
    #[pyo3(get)]
    element: u8,
    #[pyo3(get)]
    shift_ppm: f64,
    #[pyo3(get)]
    environment: String,
    #[pyo3(get)]
    confidence: f64,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct NmrShiftResultPy {
    #[pyo3(get)]
    h_shifts: Vec<ChemicalShiftPy>,
    #[pyo3(get)]
    c_shifts: Vec<ChemicalShiftPy>,
    #[pyo3(get)]
    f_shifts: Vec<ChemicalShiftPy>,
    #[pyo3(get)]
    p_shifts: Vec<ChemicalShiftPy>,
    #[pyo3(get)]
    n_shifts: Vec<ChemicalShiftPy>,
    #[pyo3(get)]
    b_shifts: Vec<ChemicalShiftPy>,
    #[pyo3(get)]
    si_shifts: Vec<ChemicalShiftPy>,
    #[pyo3(get)]
    se_shifts: Vec<ChemicalShiftPy>,
    #[pyo3(get)]
    o_shifts: Vec<ChemicalShiftPy>,
    #[pyo3(get)]
    s_shifts: Vec<ChemicalShiftPy>,
    #[pyo3(get)]
    notes: Vec<String>,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct JCouplingPy {
    #[pyo3(get)]
    h1_index: usize,
    #[pyo3(get)]
    h2_index: usize,
    #[pyo3(get)]
    j_hz: f64,
    #[pyo3(get)]
    n_bonds: usize,
    #[pyo3(get)]
    coupling_type: String,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct NmrPeakPy {
    #[pyo3(get)]
    shift_ppm: f64,
    #[pyo3(get)]
    intensity: f64,
    #[pyo3(get)]
    atom_index: usize,
    #[pyo3(get)]
    multiplicity: String,
    #[pyo3(get)]
    environment: String,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct NmrSpectrumPy {
    #[pyo3(get)]
    ppm_axis: Vec<f64>,
    #[pyo3(get)]
    intensities: Vec<f64>,
    #[pyo3(get)]
    peaks: Vec<NmrPeakPy>,
    #[pyo3(get)]
    nucleus: String,
    #[pyo3(get)]
    gamma: f64,
    #[pyo3(get)]
    notes: Vec<String>,
}

fn map_shifts(shifts: &[sci_form_core::nmr::ChemicalShift]) -> Vec<ChemicalShiftPy> {
    shifts
        .iter()
        .map(|s| ChemicalShiftPy {
            atom_index: s.atom_index,
            element: s.element,
            shift_ppm: s.shift_ppm,
            environment: s.environment.clone(),
            confidence: s.confidence,
        })
        .collect()
}

#[pyfunction]
fn nmr_shifts(smiles: &str) -> PyResult<NmrShiftResultPy> {
    sci_form_core::predict_nmr_shifts(smiles)
        .map(|r| NmrShiftResultPy {
            h_shifts: map_shifts(&r.h_shifts),
            c_shifts: map_shifts(&r.c_shifts),
            f_shifts: map_shifts(&r.f_shifts),
            p_shifts: map_shifts(&r.p_shifts),
            n_shifts: map_shifts(&r.n_shifts),
            b_shifts: map_shifts(&r.b_shifts),
            si_shifts: map_shifts(&r.si_shifts),
            se_shifts: map_shifts(&r.se_shifts),
            o_shifts: map_shifts(&r.o_shifts),
            s_shifts: map_shifts(&r.s_shifts),
            notes: r.notes,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
#[pyo3(signature = (smiles, coords=vec![]))]
fn nmr_couplings(smiles: &str, coords: Vec<f64>) -> PyResult<Vec<JCouplingPy>> {
    let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    sci_form_core::predict_nmr_couplings(smiles, &positions)
        .map(|couplings| {
            couplings
                .iter()
                .map(|c| JCouplingPy {
                    h1_index: c.h1_index,
                    h2_index: c.h2_index,
                    j_hz: c.j_hz,
                    n_bonds: c.n_bonds,
                    coupling_type: c.coupling_type.clone(),
                })
                .collect()
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
#[pyo3(signature = (smiles, nucleus="1H", gamma=0.02, ppm_min=0.0, ppm_max=12.0, n_points=1000))]
fn nmr_spectrum(
    smiles: &str,
    nucleus: &str,
    gamma: f64,
    ppm_min: f64,
    ppm_max: f64,
    n_points: usize,
) -> PyResult<NmrSpectrumPy> {
    sci_form_core::compute_nmr_spectrum(smiles, nucleus, gamma, ppm_min, ppm_max, n_points)
        .map(|r| NmrSpectrumPy {
            ppm_axis: r.ppm_axis,
            intensities: r.intensities,
            peaks: r
                .peaks
                .iter()
                .map(|p| NmrPeakPy {
                    shift_ppm: p.shift_ppm,
                    intensity: p.intensity,
                    atom_index: p.atom_index,
                    multiplicity: p.multiplicity.clone(),
                    environment: p.environment.clone(),
                })
                .collect(),
            nucleus: format!("{:?}", r.nucleus),
            gamma: r.gamma,
            notes: r.notes,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
#[pyo3(signature = (smiles, max_radius=2))]
fn hose_codes(smiles: &str, max_radius: usize) -> PyResult<Vec<(usize, u8, String)>> {
    sci_form_core::compute_hose_codes(smiles, max_radius)
        .map(|codes| {
            codes
                .iter()
                .map(|c| (c.atom_index, c.element, c.full_code.clone()))
                .collect()
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(nmr_shifts, m)?)?;
    m.add_function(wrap_pyfunction!(nmr_couplings, m)?)?;
    m.add_function(wrap_pyfunction!(nmr_spectrum, m)?)?;
    m.add_function(wrap_pyfunction!(hose_codes, m)?)?;
    m.add_class::<ChemicalShiftPy>()?;
    m.add_class::<NmrShiftResultPy>()?;
    m.add_class::<JCouplingPy>()?;
    m.add_class::<NmrPeakPy>()?;
    m.add_class::<NmrSpectrumPy>()?;
    Ok(())
}
