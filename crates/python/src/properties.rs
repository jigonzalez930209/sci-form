//! Property bindings: charges, SASA, dipole, DOS, RMSD.

use crate::system::coords_to_positions;
use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
pub(crate) struct ChargeResultPy {
    #[pyo3(get)]
    charges: Vec<f64>,
    #[pyo3(get)]
    iterations: usize,
    #[pyo3(get)]
    total_charge: f64,
}
#[pymethods]
impl ChargeResultPy {
    fn __repr__(&self) -> String {
        format!(
            "ChargeResult(n_atoms={}, total_charge={:.4}, iterations={})",
            self.charges.len(),
            self.total_charge,
            self.iterations
        )
    }
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct SasaResultPy {
    #[pyo3(get)]
    total_sasa: f64,
    #[pyo3(get)]
    atom_sasa: Vec<f64>,
    #[pyo3(get)]
    probe_radius: f64,
    #[pyo3(get)]
    num_points: usize,
}
#[pymethods]
impl SasaResultPy {
    fn __repr__(&self) -> String {
        format!(
            "SasaResult(total={:.2} Ų, n_atoms={}, probe={:.2} Å)",
            self.total_sasa,
            self.atom_sasa.len(),
            self.probe_radius
        )
    }
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct DipoleResultPy {
    #[pyo3(get)]
    vector: [f64; 3],
    #[pyo3(get)]
    magnitude: f64,
    #[pyo3(get)]
    unit: String,
}
#[pymethods]
impl DipoleResultPy {
    fn __repr__(&self) -> String {
        format!("DipoleResult({:.3} {})", self.magnitude, self.unit)
    }
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct DosResultPy {
    #[pyo3(get)]
    energies: Vec<f64>,
    #[pyo3(get)]
    total_dos: Vec<f64>,
    #[pyo3(get)]
    sigma: f64,
}
#[pymethods]
impl DosResultPy {
    fn __repr__(&self) -> String {
        format!(
            "DosResult(n_points={}, sigma={:.2} eV)",
            self.energies.len(),
            self.sigma
        )
    }
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct AlignmentResultPy {
    #[pyo3(get)]
    rmsd: f64,
    #[pyo3(get)]
    aligned_coords: Vec<f64>,
}
#[pymethods]
impl AlignmentResultPy {
    fn __repr__(&self) -> String {
        format!("AlignmentResult(rmsd={:.4} Å)", self.rmsd)
    }
}

#[pyfunction]
fn charges(smiles: &str) -> PyResult<ChargeResultPy> {
    sci_form_core::compute_charges(smiles)
        .map(|r| ChargeResultPy {
            charges: r.charges,
            iterations: r.iterations,
            total_charge: r.total_charge,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
#[pyo3(signature = (elements, coords, probe_radius=1.4))]
fn sasa(elements: Vec<u8>, coords: Vec<f64>, probe_radius: f64) -> PyResult<SasaResultPy> {
    let pr = if (probe_radius - 1.4).abs() < 1e-10 {
        None
    } else {
        Some(probe_radius)
    };
    sci_form_core::compute_sasa(&elements, &coords, pr)
        .map(|r| SasaResultPy {
            total_sasa: r.total_sasa,
            atom_sasa: r.atom_sasa,
            probe_radius: r.probe_radius,
            num_points: r.num_points,
        })
        .map_err(pyo3::exceptions::PyValueError::new_err)
}

#[pyfunction]
fn dipole(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<DipoleResultPy> {
    let positions = coords_to_positions(&coords);
    sci_form_core::compute_dipole(&elements, &positions)
        .map(|r| DipoleResultPy {
            vector: r.vector,
            magnitude: r.magnitude,
            unit: r.unit,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
#[pyo3(signature = (elements, coords, sigma=0.3, e_min=-30.0, e_max=5.0, n_points=500))]
fn dos(
    elements: Vec<u8>,
    coords: Vec<f64>,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> PyResult<DosResultPy> {
    let positions = coords_to_positions(&coords);
    sci_form_core::compute_dos(&elements, &positions, sigma, e_min, e_max, n_points)
        .map(|r| DosResultPy {
            energies: r.energies,
            total_dos: r.total_dos,
            sigma: r.sigma,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
fn rmsd(coords: Vec<f64>, reference: Vec<f64>) -> PyResult<AlignmentResultPy> {
    let result = sci_form_core::alignment::align_coordinates(&coords, &reference);
    Ok(AlignmentResultPy {
        rmsd: result.rmsd,
        aligned_coords: result.aligned_coords,
    })
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(charges, m)?)?;
    m.add_function(wrap_pyfunction!(sasa, m)?)?;
    m.add_function(wrap_pyfunction!(dipole, m)?)?;
    m.add_function(wrap_pyfunction!(dos, m)?)?;
    m.add_function(wrap_pyfunction!(rmsd, m)?)?;
    m.add_class::<ChargeResultPy>()?;
    m.add_class::<SasaResultPy>()?;
    m.add_class::<DipoleResultPy>()?;
    m.add_class::<DosResultPy>()?;
    m.add_class::<AlignmentResultPy>()?;
    Ok(())
}
