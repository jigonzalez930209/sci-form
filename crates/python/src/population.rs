//! Population analysis, bond orders, frontier descriptors.

use crate::system::coords_to_positions;
use pyo3::prelude::*;

#[pyclass(name = "PopulationResult")]
#[derive(Clone)]
pub(crate) struct PopulationResultPy {
    #[pyo3(get)]
    mulliken_charges: Vec<f64>,
    #[pyo3(get)]
    lowdin_charges: Vec<f64>,
    #[pyo3(get)]
    num_atoms: usize,
    #[pyo3(get)]
    total_charge_mulliken: f64,
    #[pyo3(get)]
    total_charge_lowdin: f64,
    #[pyo3(get)]
    charge_conservation_error: f64,
    #[pyo3(get)]
    homo_energy: Option<f64>,
    #[pyo3(get)]
    lumo_energy: Option<f64>,
    #[pyo3(get)]
    gap: Option<f64>,
    #[pyo3(get)]
    warnings: Vec<String>,
}
#[pymethods]
impl PopulationResultPy {
    fn __repr__(&self) -> String {
        format!(
            "PopulationResult(n_atoms={}, total_mulliken={:.4}, conservation_err={:.1e})",
            self.num_atoms, self.total_charge_mulliken, self.charge_conservation_error
        )
    }
}

#[pyclass(name = "BondOrderResult")]
#[derive(Clone)]
pub(crate) struct BondOrderResultPy {
    #[pyo3(get)]
    atom_pairs: Vec<(usize, usize)>,
    #[pyo3(get)]
    distances: Vec<f64>,
    #[pyo3(get)]
    wiberg: Vec<f64>,
    #[pyo3(get)]
    mayer: Vec<f64>,
    #[pyo3(get)]
    wiberg_valence: Vec<f64>,
    #[pyo3(get)]
    mayer_valence: Vec<f64>,
    #[pyo3(get)]
    num_atoms: usize,
}
#[pymethods]
impl BondOrderResultPy {
    fn __repr__(&self) -> String {
        format!(
            "BondOrderResult(n_atoms={}, n_pairs={})",
            self.num_atoms,
            self.atom_pairs.len()
        )
    }
}

#[pyclass(name = "FrontierDescriptors")]
#[derive(Clone)]
pub(crate) struct FrontierDescriptorsPy {
    #[pyo3(get)]
    homo_atom_contributions: Vec<f64>,
    #[pyo3(get)]
    lumo_atom_contributions: Vec<f64>,
    #[pyo3(get)]
    dual_descriptor: Vec<f64>,
    #[pyo3(get)]
    homo_energy: f64,
    #[pyo3(get)]
    lumo_energy: f64,
    #[pyo3(get)]
    gap: f64,
    #[pyo3(get)]
    num_atoms: usize,
}
#[pymethods]
impl FrontierDescriptorsPy {
    fn __repr__(&self) -> String {
        format!(
            "FrontierDescriptors(n_atoms={}, gap={:.3} eV)",
            self.num_atoms, self.gap
        )
    }
}

#[pyfunction]
fn compute_population(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<PopulationResultPy> {
    let positions = coords_to_positions(&coords);
    sci_form_core::compute_population(&elements, &positions)
        .map(|r| PopulationResultPy {
            mulliken_charges: r.mulliken_charges,
            lowdin_charges: r.lowdin_charges,
            num_atoms: r.num_atoms,
            total_charge_mulliken: r.total_charge_mulliken,
            total_charge_lowdin: r.total_charge_lowdin,
            charge_conservation_error: r.charge_conservation_error,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap: r.gap,
            warnings: r.warnings,
        })
        .map_err(crate::errors::engine_error)
}

#[pyfunction]
fn compute_bond_orders(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<BondOrderResultPy> {
    let positions = coords_to_positions(&coords);
    sci_form_core::compute_bond_orders(&elements, &positions)
        .map(|r| BondOrderResultPy {
            atom_pairs: r.bonds.iter().map(|b| (b.atom_i, b.atom_j)).collect(),
            distances: r.bonds.iter().map(|b| b.distance).collect(),
            wiberg: r.bonds.iter().map(|b| b.wiberg).collect(),
            mayer: r.bonds.iter().map(|b| b.mayer).collect(),
            wiberg_valence: r.wiberg_valence,
            mayer_valence: r.mayer_valence,
            num_atoms: r.num_atoms,
        })
        .map_err(crate::errors::engine_error)
}

#[pyfunction]
fn compute_frontier_descriptors(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<FrontierDescriptorsPy> {
    let positions = coords_to_positions(&coords);
    sci_form_core::compute_frontier_descriptors(&elements, &positions)
        .map(|r| FrontierDescriptorsPy {
            homo_atom_contributions: r.homo_atom_contributions,
            lumo_atom_contributions: r.lumo_atom_contributions,
            dual_descriptor: r.dual_descriptor,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap: r.gap,
            num_atoms: r.num_atoms,
        })
        .map_err(crate::errors::engine_error)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_population, m)?)?;
    m.add_function(wrap_pyfunction!(compute_bond_orders, m)?)?;
    m.add_function(wrap_pyfunction!(compute_frontier_descriptors, m)?)?;
    m.add_class::<PopulationResultPy>()?;
    m.add_class::<BondOrderResultPy>()?;
    m.add_class::<FrontierDescriptorsPy>()?;
    Ok(())
}
