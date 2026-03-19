//! UFF, MMFF94, and empirical pKa bindings.

use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
pub(crate) struct UffHeuristicEnergyPy {
    #[pyo3(get)]
    raw_energy_kcal_mol: f64,
    #[pyo3(get)]
    aromatic_stabilization_kcal_mol: f64,
    #[pyo3(get)]
    corrected_energy_kcal_mol: f64,
    #[pyo3(get)]
    aromatic_bond_count: usize,
    #[pyo3(get)]
    notes: Vec<String>,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct EmpiricalPkaSitePy {
    #[pyo3(get)]
    atom_index: usize,
    #[pyo3(get)]
    site_type: String,
    #[pyo3(get)]
    environment: String,
    #[pyo3(get)]
    estimated_pka: f64,
    #[pyo3(get)]
    confidence: f64,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct EmpiricalPkaResultPy {
    #[pyo3(get)]
    acidic_sites: Vec<EmpiricalPkaSitePy>,
    #[pyo3(get)]
    basic_sites: Vec<EmpiricalPkaSitePy>,
    #[pyo3(get)]
    notes: Vec<String>,
}

#[pyfunction]
fn uff_energy(smiles: &str, coords: Vec<f64>) -> PyResult<f64> {
    sci_form_core::compute_uff_energy(smiles, &coords)
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
fn uff_energy_with_aromatic_heuristics(
    smiles: &str,
    coords: Vec<f64>,
) -> PyResult<UffHeuristicEnergyPy> {
    sci_form_core::compute_uff_energy_with_aromatic_heuristics(smiles, &coords)
        .map(|r| UffHeuristicEnergyPy {
            raw_energy_kcal_mol: r.raw_energy_kcal_mol,
            aromatic_stabilization_kcal_mol: r.aromatic_stabilization_kcal_mol,
            corrected_energy_kcal_mol: r.corrected_energy_kcal_mol,
            aromatic_bond_count: r.aromatic_bond_count,
            notes: r.notes,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
fn empirical_pka(smiles: &str) -> PyResult<EmpiricalPkaResultPy> {
    sci_form_core::compute_empirical_pka(smiles)
        .map(|r| {
            let map_site = |s: sci_form_core::reactivity::EmpiricalPkaSite| EmpiricalPkaSitePy {
                atom_index: s.atom_index,
                site_type: s.site_type,
                environment: s.environment,
                estimated_pka: s.estimated_pka,
                confidence: s.confidence,
            };
            EmpiricalPkaResultPy {
                acidic_sites: r.acidic_sites.into_iter().map(map_site).collect(),
                basic_sites: r.basic_sites.into_iter().map(map_site).collect(),
                notes: r.notes,
            }
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
fn mmff94_energy(smiles: &str, coords: Vec<f64>) -> PyResult<f64> {
    sci_form_core::compute_mmff94_energy(smiles, &coords)
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(uff_energy, m)?)?;
    m.add_function(wrap_pyfunction!(uff_energy_with_aromatic_heuristics, m)?)?;
    m.add_function(wrap_pyfunction!(empirical_pka, m)?)?;
    m.add_function(wrap_pyfunction!(mmff94_energy, m)?)?;
    m.add_class::<UffHeuristicEnergyPy>()?;
    m.add_class::<EmpiricalPkaSitePy>()?;
    m.add_class::<EmpiricalPkaResultPy>()?;
    Ok(())
}
