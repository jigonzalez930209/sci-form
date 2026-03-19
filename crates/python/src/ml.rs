//! ML molecular descriptors and property prediction bindings.

use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
pub(crate) struct MolecularDescriptorsPy {
    #[pyo3(get)]
    molecular_weight: f64,
    #[pyo3(get)]
    n_heavy_atoms: usize,
    #[pyo3(get)]
    n_hydrogens: usize,
    #[pyo3(get)]
    n_bonds: usize,
    #[pyo3(get)]
    n_rotatable_bonds: usize,
    #[pyo3(get)]
    n_hbd: usize,
    #[pyo3(get)]
    n_hba: usize,
    #[pyo3(get)]
    fsp3: f64,
    #[pyo3(get)]
    wiener_index: f64,
    #[pyo3(get)]
    n_rings: usize,
    #[pyo3(get)]
    n_aromatic: usize,
    #[pyo3(get)]
    sum_electronegativity: f64,
    #[pyo3(get)]
    sum_polarizability: f64,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct MlPropertyResultPy {
    #[pyo3(get)]
    logp: f64,
    #[pyo3(get)]
    molar_refractivity: f64,
    #[pyo3(get)]
    log_solubility: f64,
    #[pyo3(get)]
    lipinski_violations: u8,
    #[pyo3(get)]
    lipinski_passes: bool,
    #[pyo3(get)]
    druglikeness: f64,
}

fn parse_mol_for_ml(smiles: &str) -> PyResult<(Vec<u8>, Vec<(usize, usize, u8)>)> {
    let mol = sci_form_core::parse(smiles).map_err(pyo3::exceptions::PyValueError::new_err)?;
    let n = mol.graph.node_count();
    let elements: Vec<u8> = (0..n)
        .map(|i| mol.graph[sci_form_core::graph::NodeIndex::new(i)].element)
        .collect();
    let bonds: Vec<(usize, usize, u8)> = mol
        .graph
        .edge_indices()
        .map(|e| {
            let (a, b) = mol.graph.edge_endpoints(e).unwrap();
            let order = match mol.graph[e].order {
                sci_form_core::graph::BondOrder::Single => 1u8,
                sci_form_core::graph::BondOrder::Double => 2,
                sci_form_core::graph::BondOrder::Triple => 3,
                sci_form_core::graph::BondOrder::Aromatic => 2,
                sci_form_core::graph::BondOrder::Unknown => 1,
            };
            (a.index(), b.index(), order)
        })
        .collect();
    Ok((elements, bonds))
}

#[pyfunction]
fn ml_descriptors(smiles: &str) -> PyResult<MolecularDescriptorsPy> {
    let (elements, bonds) = parse_mol_for_ml(smiles)?;
    let desc = sci_form_core::compute_ml_descriptors(&elements, &bonds, &[], &[]);
    Ok(MolecularDescriptorsPy {
        molecular_weight: desc.molecular_weight,
        n_heavy_atoms: desc.n_heavy_atoms,
        n_hydrogens: desc.n_hydrogens,
        n_bonds: desc.n_bonds,
        n_rotatable_bonds: desc.n_rotatable_bonds,
        n_hbd: desc.n_hbd,
        n_hba: desc.n_hba,
        fsp3: desc.fsp3,
        wiener_index: desc.wiener_index,
        n_rings: desc.n_rings,
        n_aromatic: desc.n_aromatic,
        sum_electronegativity: desc.sum_electronegativity,
        sum_polarizability: desc.sum_polarizability,
    })
}

#[pyfunction]
fn ml_predict(smiles: &str) -> PyResult<MlPropertyResultPy> {
    let (elements, bonds) = parse_mol_for_ml(smiles)?;
    let desc = sci_form_core::compute_ml_descriptors(&elements, &bonds, &[], &[]);
    let result = sci_form_core::predict_ml_properties(&desc);
    Ok(MlPropertyResultPy {
        logp: result.logp,
        molar_refractivity: result.molar_refractivity,
        log_solubility: result.log_solubility,
        lipinski_violations: result.lipinski.violations,
        lipinski_passes: result.lipinski.passes,
        druglikeness: result.druglikeness,
    })
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(ml_descriptors, m)?)?;
    m.add_function(wrap_pyfunction!(ml_predict, m)?)?;
    m.add_class::<MolecularDescriptorsPy>()?;
    m.add_class::<MlPropertyResultPy>()?;
    Ok(())
}
