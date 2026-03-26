//! HF-3c, ANI neural potential, and ESP grid bindings.

use crate::system::coords_to_positions;
use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
pub(crate) struct Hf3cResultPy {
    #[pyo3(get)]
    energy: f64,
    #[pyo3(get)]
    hf_energy: f64,
    #[pyo3(get)]
    nuclear_repulsion: f64,
    #[pyo3(get)]
    d3_energy: f64,
    #[pyo3(get)]
    gcp_energy: f64,
    #[pyo3(get)]
    srb_energy: f64,
    #[pyo3(get)]
    orbital_energies: Vec<f64>,
    #[pyo3(get)]
    scf_iterations: usize,
    #[pyo3(get)]
    converged: bool,
    #[pyo3(get)]
    n_basis: usize,
    #[pyo3(get)]
    n_electrons: usize,
    #[pyo3(get)]
    homo_energy: f64,
    #[pyo3(get)]
    lumo_energy: Option<f64>,
    #[pyo3(get)]
    gap: f64,
    #[pyo3(get)]
    mulliken_charges: Vec<f64>,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct AniResultPy {
    #[pyo3(get)]
    energy: f64,
    #[pyo3(get)]
    forces: Vec<(f64, f64, f64)>,
    #[pyo3(get)]
    species: Vec<u8>,
    #[pyo3(get)]
    atomic_energies: Vec<f64>,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct EspGridPy {
    #[pyo3(get)]
    values: Vec<f64>,
    #[pyo3(get)]
    origin: (f64, f64, f64),
    #[pyo3(get)]
    spacing: f64,
    #[pyo3(get)]
    dims: (usize, usize, usize),
}

#[pyfunction]
#[pyo3(signature = (elements, coords, max_scf_iter=300, n_cis_states=5, corrections=true))]
fn hf3c_calculate(
    elements: Vec<u8>,
    coords: Vec<f64>,
    max_scf_iter: usize,
    n_cis_states: usize,
    corrections: bool,
) -> PyResult<Hf3cResultPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    let positions = coords_to_positions(&coords);
    let config = sci_form_core::hf::HfConfig {
        max_scf_iter,
        diis_size: 6,
        n_cis_states,
        corrections,
    };
    sci_form_core::compute_hf3c(&elements, &positions, &config)
        .map(|r| Hf3cResultPy {
            energy: r.energy,
            hf_energy: r.hf_energy,
            nuclear_repulsion: r.nuclear_repulsion,
            d3_energy: r.d3_energy,
            gcp_energy: r.gcp_energy,
            srb_energy: r.srb_energy,
            orbital_energies: r.orbital_energies,
            scf_iterations: r.scf_iterations,
            converged: r.converged,
            n_basis: r.n_basis,
            n_electrons: r.n_electrons,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap: r.gap,
            mulliken_charges: r.mulliken_charges,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
fn ani_calculate(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<AniResultPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    let positions = coords_to_positions(&coords);
    sci_form_core::compute_ani(&elements, &positions)
        .map(|r| AniResultPy {
            energy: r.energy,
            forces: r.forces.iter().map(|f| (f[0], f[1], f[2])).collect(),
            species: r.species,
            atomic_energies: r.atomic_energies,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
#[pyo3(signature = (elements, coords, spacing=0.5, padding=3.0))]
fn esp(elements: Vec<u8>, coords: Vec<f64>, spacing: f64, padding: f64) -> PyResult<EspGridPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    let positions = coords_to_positions(&coords);
    sci_form_core::compute_esp(&elements, &positions, spacing, padding)
        .map(|grid| EspGridPy {
            values: grid.values,
            origin: (grid.origin[0], grid.origin[1], grid.origin[2]),
            spacing: grid.spacing,
            dims: (grid.dims[0], grid.dims[1], grid.dims[2]),
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(hf3c_calculate, m)?)?;
    m.add_function(wrap_pyfunction!(ani_calculate, m)?)?;
    m.add_function(wrap_pyfunction!(esp, m)?)?;
    m.add_class::<Hf3cResultPy>()?;
    m.add_class::<AniResultPy>()?;
    m.add_class::<EspGridPy>()?;
    Ok(())
}
