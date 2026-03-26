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
    other_shifts: Vec<NucleusShiftSeriesPy>,
    #[pyo3(get)]
    notes: Vec<String>,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct NucleusShiftSeriesPy {
    #[pyo3(get)]
    nucleus: String,
    #[pyo3(get)]
    shifts: Vec<ChemicalShiftPy>,
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

#[pyclass]
#[derive(Clone)]
pub(crate) struct GiaoNmrEntryPy {
    #[pyo3(get)]
    atom_index: usize,
    #[pyo3(get)]
    element: u8,
    #[pyo3(get)]
    tensor: Vec<Vec<f64>>,
    #[pyo3(get)]
    isotropic: f64,
    #[pyo3(get)]
    anisotropy: f64,
    #[pyo3(get)]
    chemical_shift: f64,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct GiaoNmrResultPy {
    #[pyo3(get)]
    nucleus: String,
    #[pyo3(get)]
    target_atomic_number: u8,
    #[pyo3(get)]
    method: String,
    #[pyo3(get)]
    basis_set: String,
    #[pyo3(get)]
    charge: i32,
    #[pyo3(get)]
    multiplicity: u32,
    #[pyo3(get)]
    scf_converged: bool,
    #[pyo3(get)]
    scf_iterations: usize,
    #[pyo3(get)]
    total_energy_hartree: f64,
    #[pyo3(get)]
    n_basis: usize,
    #[pyo3(get)]
    n_target_atoms: usize,
    #[pyo3(get)]
    chemical_shifts: Vec<f64>,
    #[pyo3(get)]
    shieldings: Vec<GiaoNmrEntryPy>,
    #[pyo3(get)]
    fallback_elements: Vec<u8>,
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

fn map_shift_series(series: &[sci_form_core::nmr::NucleusShiftSeries]) -> Vec<NucleusShiftSeriesPy> {
    series
        .iter()
        .map(|entry| NucleusShiftSeriesPy {
            nucleus: entry.nucleus.canonical().to_string(),
            shifts: map_shifts(&entry.shifts),
        })
        .collect()
}

fn map_giao_result(result: sci_form_core::GiaoNmrResult) -> GiaoNmrResultPy {
    GiaoNmrResultPy {
        nucleus: result.nucleus,
        target_atomic_number: result.target_atomic_number,
        method: result.method,
        basis_set: result.basis_set,
        charge: result.charge,
        multiplicity: result.multiplicity,
        scf_converged: result.scf_converged,
        scf_iterations: result.scf_iterations,
        total_energy_hartree: result.total_energy_hartree,
        n_basis: result.n_basis,
        n_target_atoms: result.n_target_atoms,
        chemical_shifts: result.chemical_shifts,
        shieldings: result
            .shieldings
            .into_iter()
            .map(|entry| GiaoNmrEntryPy {
                atom_index: entry.atom_index,
                element: entry.element,
                tensor: entry.tensor.iter().map(|row| row.to_vec()).collect(),
                isotropic: entry.isotropic,
                anisotropy: entry.anisotropy,
                chemical_shift: entry.chemical_shift,
            })
            .collect(),
        fallback_elements: result.fallback_elements,
        notes: result.notes,
    }
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
            other_shifts: map_shift_series(&r.other_shifts),
            notes: r.notes,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
fn nmr_shifts_for_nucleus(smiles: &str, nucleus: &str) -> PyResult<Vec<ChemicalShiftPy>> {
    sci_form_core::predict_nmr_shifts_for_nucleus(smiles, nucleus)
        .map(|shifts| map_shifts(&shifts))
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
#[pyo3(signature = (elements, coords, nucleus="1H", charge=0, multiplicity=1, max_scf_iter=100, allow_basis_fallback=false, use_parallel_eri=false))]
fn giao_nmr(
    elements: Vec<u8>,
    coords: Vec<f64>,
    nucleus: &str,
    charge: i32,
    multiplicity: u32,
    max_scf_iter: usize,
    allow_basis_fallback: bool,
    use_parallel_eri: bool,
) -> PyResult<GiaoNmrResultPy> {
    let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    let config = sci_form_core::GiaoNmrConfig {
        charge,
        multiplicity,
        max_scf_iterations: max_scf_iter,
        use_parallel_eri,
        allow_basis_fallback,
    };

    sci_form_core::compute_giao_nmr_configured(&elements, &positions, nucleus, &config)
        .map(map_giao_result)
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
#[pyo3(signature = (smiles, coords=vec![], nucleus="1H", gamma=0.02, ppm_min=0.0, ppm_max=12.0, n_points=1000))]
fn nmr_spectrum(
    smiles: &str,
    coords: Vec<f64>,
    nucleus: &str,
    gamma: f64,
    ppm_min: f64,
    ppm_max: f64,
    n_points: usize,
) -> PyResult<NmrSpectrumPy> {
    let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    sci_form_core::compute_nmr_spectrum_with_coords(
        smiles, &positions, nucleus, gamma, ppm_min, ppm_max, n_points,
    )
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
    m.add_function(wrap_pyfunction!(nmr_shifts_for_nucleus, m)?)?;
    m.add_function(wrap_pyfunction!(giao_nmr, m)?)?;
    m.add_function(wrap_pyfunction!(nmr_couplings, m)?)?;
    m.add_function(wrap_pyfunction!(nmr_spectrum, m)?)?;
    m.add_function(wrap_pyfunction!(hose_codes, m)?)?;
    m.add_class::<ChemicalShiftPy>()?;
    m.add_class::<NucleusShiftSeriesPy>()?;
    m.add_class::<NmrShiftResultPy>()?;
    m.add_class::<JCouplingPy>()?;
    m.add_class::<NmrPeakPy>()?;
    m.add_class::<NmrSpectrumPy>()?;
    m.add_class::<GiaoNmrEntryPy>()?;
    m.add_class::<GiaoNmrResultPy>()?;
    Ok(())
}
