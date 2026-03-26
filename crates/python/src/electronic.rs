//! PM3 and xTB-family semi-empirical electronic structure bindings.

use crate::system::coords_to_positions;
use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
pub(crate) struct Pm3ResultPy {
    #[pyo3(get)]
    orbital_energies: Vec<f64>,
    #[pyo3(get)]
    electronic_energy: f64,
    #[pyo3(get)]
    nuclear_repulsion: f64,
    #[pyo3(get)]
    total_energy: f64,
    #[pyo3(get)]
    heat_of_formation: f64,
    #[pyo3(get)]
    n_basis: usize,
    #[pyo3(get)]
    n_electrons: usize,
    #[pyo3(get)]
    homo_energy: f64,
    #[pyo3(get)]
    lumo_energy: f64,
    #[pyo3(get)]
    gap: f64,
    #[pyo3(get)]
    mulliken_charges: Vec<f64>,
    #[pyo3(get)]
    scf_iterations: usize,
    #[pyo3(get)]
    converged: bool,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct XtbResultPy {
    #[pyo3(get)]
    orbital_energies: Vec<f64>,
    #[pyo3(get)]
    total_energy: f64,
    #[pyo3(get)]
    repulsion_energy: f64,
    #[pyo3(get)]
    electronic_energy: f64,
    #[pyo3(get)]
    n_basis: usize,
    #[pyo3(get)]
    n_electrons: usize,
    #[pyo3(get)]
    homo_energy: f64,
    #[pyo3(get)]
    lumo_energy: f64,
    #[pyo3(get)]
    gap: f64,
    #[pyo3(get)]
    mulliken_charges: Vec<f64>,
    #[pyo3(get)]
    converged: bool,
    #[pyo3(get)]
    scf_iterations: usize,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct Gfn1ResultPy {
    #[pyo3(get)]
    orbital_energies: Vec<f64>,
    #[pyo3(get)]
    total_energy: f64,
    #[pyo3(get)]
    repulsive_energy: f64,
    #[pyo3(get)]
    dispersion_energy: f64,
    #[pyo3(get)]
    electronic_energy: f64,
    #[pyo3(get)]
    n_basis: usize,
    #[pyo3(get)]
    n_electrons: usize,
    #[pyo3(get)]
    homo_energy: f64,
    #[pyo3(get)]
    lumo_energy: f64,
    #[pyo3(get)]
    gap: f64,
    #[pyo3(get)]
    mulliken_charges: Vec<f64>,
    #[pyo3(get)]
    shell_charges: Vec<Vec<f64>>,
    #[pyo3(get)]
    converged: bool,
    #[pyo3(get)]
    scc_iterations: usize,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct Gfn2ResultPy {
    #[pyo3(get)]
    orbital_energies: Vec<f64>,
    #[pyo3(get)]
    total_energy: f64,
    #[pyo3(get)]
    repulsive_energy: f64,
    #[pyo3(get)]
    dispersion_energy: f64,
    #[pyo3(get)]
    halogen_bond_energy: f64,
    #[pyo3(get)]
    electronic_energy: f64,
    #[pyo3(get)]
    n_basis: usize,
    #[pyo3(get)]
    n_electrons: usize,
    #[pyo3(get)]
    homo_energy: f64,
    #[pyo3(get)]
    lumo_energy: f64,
    #[pyo3(get)]
    gap: f64,
    #[pyo3(get)]
    mulliken_charges: Vec<f64>,
    #[pyo3(get)]
    atomic_dipoles: Vec<(f64, f64, f64)>,
    #[pyo3(get)]
    atomic_quadrupoles: Vec<f64>,
    #[pyo3(get)]
    converged: bool,
    #[pyo3(get)]
    scc_iterations: usize,
}

fn validate_coords(elements: &[u8], coords: &[f64]) -> PyResult<Vec<[f64; 3]>> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    Ok(coords_to_positions(coords))
}

#[pyfunction]
fn pm3_calculate(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<Pm3ResultPy> {
    let positions = validate_coords(&elements, &coords)?;
    sci_form_core::compute_pm3(&elements, &positions)
        .map(|r| Pm3ResultPy {
            orbital_energies: r.orbital_energies,
            electronic_energy: r.electronic_energy,
            nuclear_repulsion: r.nuclear_repulsion,
            total_energy: r.total_energy,
            heat_of_formation: r.heat_of_formation,
            n_basis: r.n_basis,
            n_electrons: r.n_electrons,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap: r.gap,
            mulliken_charges: r.mulliken_charges,
            scf_iterations: r.scf_iterations,
            converged: r.converged,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
fn xtb_calculate(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<XtbResultPy> {
    let positions = validate_coords(&elements, &coords)?;
    sci_form_core::compute_xtb(&elements, &positions)
        .map(|r| XtbResultPy {
            orbital_energies: r.orbital_energies,
            total_energy: r.total_energy,
            repulsion_energy: r.repulsive_energy,
            electronic_energy: r.electronic_energy,
            n_basis: r.n_basis,
            n_electrons: r.n_electrons,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap: r.gap,
            mulliken_charges: r.mulliken_charges,
            converged: r.converged,
            scf_iterations: r.scc_iterations,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
fn gfn1_calculate(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<Gfn1ResultPy> {
    let positions = validate_coords(&elements, &coords)?;
    sci_form_core::xtb::solve_gfn1(&elements, &positions)
        .map(|r| Gfn1ResultPy {
            orbital_energies: r.orbital_energies,
            total_energy: r.total_energy,
            repulsive_energy: r.repulsive_energy,
            dispersion_energy: r.dispersion_energy,
            electronic_energy: r.electronic_energy,
            n_basis: r.n_basis,
            n_electrons: r.n_electrons,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap: r.gap,
            mulliken_charges: r.mulliken_charges,
            shell_charges: r.shell_charges,
            converged: r.converged,
            scc_iterations: r.scc_iterations,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
fn gfn2_calculate(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<Gfn2ResultPy> {
    let positions = validate_coords(&elements, &coords)?;
    sci_form_core::xtb::solve_gfn2(&elements, &positions)
        .map(|r| Gfn2ResultPy {
            orbital_energies: r.orbital_energies,
            total_energy: r.total_energy,
            repulsive_energy: r.repulsive_energy,
            dispersion_energy: r.dispersion_energy,
            halogen_bond_energy: r.halogen_bond_energy,
            electronic_energy: r.electronic_energy,
            n_basis: r.n_basis,
            n_electrons: r.n_electrons,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap: r.gap,
            mulliken_charges: r.mulliken_charges,
            atomic_dipoles: r
                .atomic_dipoles
                .into_iter()
                .map(|dipole| (dipole[0], dipole[1], dipole[2]))
                .collect(),
            atomic_quadrupoles: r.atomic_quadrupoles,
            converged: r.converged,
            scc_iterations: r.scc_iterations,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pm3_calculate, m)?)?;
    m.add_function(wrap_pyfunction!(xtb_calculate, m)?)?;
    m.add_function(wrap_pyfunction!(gfn1_calculate, m)?)?;
    m.add_function(wrap_pyfunction!(gfn2_calculate, m)?)?;
    m.add_class::<Pm3ResultPy>()?;
    m.add_class::<XtbResultPy>()?;
    m.add_class::<Gfn1ResultPy>()?;
    m.add_class::<Gfn2ResultPy>()?;
    Ok(())
}
