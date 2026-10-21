//! EHT bindings: calculation, support, orbital mesh.

use crate::system::coords_to_positions;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::PyAny;

#[pyclass]
#[derive(Clone)]
pub(crate) struct EhtResultPy {
    #[pyo3(get)]
    energies: Vec<f64>,
    #[pyo3(get)]
    n_electrons: usize,
    #[pyo3(get)]
    homo_index: usize,
    #[pyo3(get)]
    lumo_index: usize,
    #[pyo3(get)]
    homo_energy: f64,
    #[pyo3(get)]
    lumo_energy: f64,
    #[pyo3(get)]
    gap: f64,
    #[pyo3(get)]
    support_level: String,
    #[pyo3(get)]
    has_transition_metals: bool,
    #[pyo3(get)]
    supported_elements: Vec<u8>,
    #[pyo3(get)]
    provisional_elements: Vec<u8>,
    #[pyo3(get)]
    unsupported_elements: Vec<u8>,
    #[pyo3(get)]
    warnings: Vec<String>,
}

#[pymethods]
impl EhtResultPy {
    fn __repr__(&self) -> String {
        format!(
            "EhtResult(n_mo={}, gap={:.3} eV, HOMO={:.3} eV, LUMO={:.3} eV)",
            self.energies.len(),
            self.gap,
            self.homo_energy,
            self.lumo_energy
        )
    }
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct EhtSupportPy {
    #[pyo3(get)]
    level: String,
    #[pyo3(get)]
    has_transition_metals: bool,
    #[pyo3(get)]
    supported_elements: Vec<u8>,
    #[pyo3(get)]
    provisional_elements: Vec<u8>,
    #[pyo3(get)]
    unsupported_elements: Vec<u8>,
    #[pyo3(get)]
    warnings: Vec<String>,
}

impl From<sci_form_core::eht::EhtSupport> for EhtSupportPy {
    fn from(v: sci_form_core::eht::EhtSupport) -> Self {
        Self {
            level: format!("{:?}", v.level).to_lowercase(),
            has_transition_metals: v.has_transition_metals,
            supported_elements: v.supported_elements,
            provisional_elements: v.provisional_elements,
            unsupported_elements: v.unsupported_elements,
            warnings: v.warnings,
        }
    }
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct EhtGradientPy {
    #[pyo3(get)]
    gradients: Vec<[f64; 3]>,
    #[pyo3(get)]
    rms_gradient: f64,
    #[pyo3(get)]
    max_gradient: f64,
    #[pyo3(get)]
    energy: f64,
}
#[pymethods]
impl EhtGradientPy {
    fn __repr__(&self) -> String {
        format!(
            "EhtGradient(n_atoms={}, rms={:.6} eV/A, E={:.4} eV)",
            self.gradients.len(),
            self.rms_gradient,
            self.energy
        )
    }
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct EhtOptResultPy {
    #[pyo3(get)]
    positions: Vec<[f64; 3]>,
    #[pyo3(get)]
    energy: f64,
    #[pyo3(get)]
    n_steps: usize,
    #[pyo3(get)]
    converged: bool,
    #[pyo3(get)]
    rms_gradient: f64,
    #[pyo3(get)]
    energies: Vec<f64>,
}
#[pymethods]
impl EhtOptResultPy {
    fn __repr__(&self) -> String {
        format!(
            "EhtOptResult(n_atoms={}, steps={}, converged={}, E={:.4} eV)",
            self.positions.len(),
            self.n_steps,
            self.converged,
            self.energy
        )
    }
}

/// Compute EHT energy-surface gradients (eV/Å) by central finite differences.
#[pyfunction]
fn compute_eht_gradient(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<EhtGradientPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "coords length must be 3 * len(elements)",
        ));
    }
    let positions = coords_to_positions(&coords);
    sci_form_core::eht::gradients::compute_eht_gradient(&elements, &positions)
        .map(|g| EhtGradientPy {
            gradients: g.gradients.clone(),
            rms_gradient: g.rms_gradient,
            max_gradient: g.max_gradient,
            energy: g.energy,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

/// Optimize a geometry on the EHT energy surface with the shared BFGS engine.
#[pyfunction]
#[pyo3(signature = (elements, coords, max_steps=200, grad_threshold=0.01, energy_threshold=1e-6, step_size=0.05))]
fn optimize_geometry_eht(
    elements: Vec<u8>,
    coords: Vec<f64>,
    max_steps: usize,
    grad_threshold: f64,
    energy_threshold: f64,
    step_size: f64,
) -> PyResult<EhtOptResultPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "coords length must be 3 * len(elements)",
        ));
    }
    let positions = coords_to_positions(&coords);
    let config = sci_form_core::eht::gradients::EhtOptConfig {
        max_steps,
        grad_threshold,
        energy_threshold,
        step_size,
    };
    sci_form_core::eht::gradients::optimize_geometry_eht(&elements, &positions, Some(config))
        .map(|r| EhtOptResultPy {
            positions: r.positions.clone(),
            energy: r.energy,
            n_steps: r.n_steps,
            converged: r.converged,
            rms_gradient: r.rms_gradient,
            energies: r.energies,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
#[pyo3(signature = (elements, coords, k=1.75))]
fn eht_calculate(elements: Vec<u8>, coords: Vec<f64>, k: f64) -> PyResult<EhtResultPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "coords length must be 3 * len(elements)",
        ));
    }
    let positions = coords_to_positions(&coords);
    let k_opt = if (k - 1.75).abs() < 1e-10 {
        None
    } else {
        Some(k)
    };
    match sci_form_core::eht::solve_eht(&elements, &positions, k_opt) {
        Ok(r) => Ok(EhtResultPy {
            energies: r.energies,
            n_electrons: r.n_electrons,
            homo_index: r.homo_index,
            lumo_index: r.lumo_index,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap: r.gap,
            support_level: format!("{:?}", r.support.level).to_lowercase(),
            has_transition_metals: r.support.has_transition_metals,
            supported_elements: r.support.supported_elements,
            provisional_elements: r.support.provisional_elements,
            unsupported_elements: r.support.unsupported_elements,
            warnings: r.support.warnings,
        }),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

#[pyfunction]
#[pyo3(signature = (smiles, elements, coords, allow_experimental_eht=false))]
fn eht_or_uff_fallback(
    py: Python<'_>,
    smiles: &str,
    elements: Vec<u8>,
    coords: Vec<f64>,
    allow_experimental_eht: bool,
) -> PyResult<Py<PyAny>> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "coords length must be 3 * len(elements)",
        ));
    }
    let positions = coords_to_positions(&coords);
    let result = sci_form_core::compute_eht_or_uff_fallback(
        smiles,
        &elements,
        &positions,
        allow_experimental_eht,
    )
    .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
    let dict = PyDict::new(py);
    match result {
        sci_form_core::ElectronicWorkflowResult::Eht { result } => {
            dict.set_item("mode", "eht")?;
            dict.set_item("energies", result.energies)?;
            dict.set_item("n_electrons", result.n_electrons)?;
            dict.set_item("homo_index", result.homo_index)?;
            dict.set_item("lumo_index", result.lumo_index)?;
            dict.set_item("homo_energy", result.homo_energy)?;
            dict.set_item("lumo_energy", result.lumo_energy)?;
            dict.set_item("gap", result.gap)?;
            dict.set_item(
                "support_level",
                format!("{:?}", result.support.level).to_lowercase(),
            )?;
            dict.set_item("warnings", result.support.warnings)?;
        }
        sci_form_core::ElectronicWorkflowResult::UffFallback {
            energy_kcal_mol,
            reason,
            support,
        } => {
            dict.set_item("mode", "uff_fallback")?;
            dict.set_item("energy_kcal_mol", energy_kcal_mol)?;
            dict.set_item("reason", reason)?;
            dict.set_item(
                "support_level",
                format!("{:?}", support.level).to_lowercase(),
            )?;
            dict.set_item("warnings", support.warnings)?;
        }
    }
    Ok(dict.into())
}

#[pyfunction]
fn eht_support(elements: Vec<u8>) -> EhtSupportPy {
    sci_form_core::get_eht_support(&elements).into()
}

#[pyfunction]
#[pyo3(signature = (elements, coords, mo_index, spacing=0.2, isovalue=0.02))]
fn eht_orbital_mesh(
    py: Python<'_>,
    elements: Vec<u8>,
    coords: Vec<f64>,
    mo_index: usize,
    spacing: f64,
    isovalue: f32,
) -> PyResult<Py<PyAny>> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "coords length must be 3 * len(elements)",
        ));
    }
    let positions = coords_to_positions(&coords);
    let result = sci_form_core::eht::solve_eht(&elements, &positions, None)
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
    let basis = sci_form_core::eht::basis::build_basis(&elements, &positions);
    let grid = sci_form_core::eht::evaluate_orbital_on_grid(
        &basis,
        &result.coefficients,
        mo_index,
        &positions,
        spacing,
        3.0,
    );
    let mesh = sci_form_core::eht::marching_cubes(&grid, isovalue);
    let dict = PyDict::new(py);
    dict.set_item("vertices", mesh.vertices)?;
    dict.set_item("normals", mesh.normals)?;
    dict.set_item("indices", mesh.indices)?;
    dict.set_item("num_triangles", mesh.num_triangles)?;
    dict.set_item("isovalue", mesh.isovalue)?;
    Ok(dict.into())
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(eht_calculate, m)?)?;
    m.add_function(wrap_pyfunction!(compute_eht_gradient, m)?)?;
    m.add_function(wrap_pyfunction!(optimize_geometry_eht, m)?)?;
    m.add_function(wrap_pyfunction!(eht_or_uff_fallback, m)?)?;
    m.add_function(wrap_pyfunction!(eht_support, m)?)?;
    m.add_function(wrap_pyfunction!(eht_orbital_mesh, m)?)?;
    m.add_class::<EhtResultPy>()?;
    m.add_class::<EhtSupportPy>()?;
    m.add_class::<EhtGradientPy>()?;
    m.add_class::<EhtOptResultPy>()?;
    Ok(())
}
