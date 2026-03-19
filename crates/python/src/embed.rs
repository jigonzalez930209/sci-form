//! Conformer embedding bindings.

use pyo3::prelude::*;
use pyo3::types::PyDict;

#[pyclass]
#[derive(Clone)]
pub(crate) struct ConformerResult {
    #[pyo3(get)]
    pub smiles: String,
    #[pyo3(get)]
    pub num_atoms: usize,
    #[pyo3(get)]
    pub coords: Vec<f64>,
    #[pyo3(get)]
    pub elements: Vec<u8>,
    #[pyo3(get)]
    pub bonds: Vec<(usize, usize, String)>,
    #[pyo3(get)]
    pub error: Option<String>,
    #[pyo3(get)]
    pub time_ms: f64,
}

#[pymethods]
impl ConformerResult {
    fn get_positions(&self) -> Vec<(f64, f64, f64)> {
        self.coords
            .chunks_exact(3)
            .map(|c| (c[0], c[1], c[2]))
            .collect()
    }
    fn is_ok(&self) -> bool {
        self.error.is_none()
    }
    fn __repr__(&self) -> String {
        if let Some(ref e) = self.error {
            format!("ConformerResult(smiles='{}', error='{}')", self.smiles, e)
        } else {
            format!(
                "ConformerResult(smiles='{}', atoms={}, time={:.1}ms)",
                self.smiles, self.num_atoms, self.time_ms
            )
        }
    }
}

impl From<sci_form_core::ConformerResult> for ConformerResult {
    fn from(r: sci_form_core::ConformerResult) -> Self {
        ConformerResult {
            smiles: r.smiles,
            num_atoms: r.num_atoms,
            coords: r.coords,
            elements: r.elements,
            bonds: r.bonds,
            error: r.error,
            time_ms: r.time_ms,
        }
    }
}

#[pyfunction]
#[pyo3(signature = (smiles, seed=42))]
fn embed(smiles: &str, seed: u64) -> ConformerResult {
    sci_form_core::embed(smiles, seed).into()
}

#[pyfunction]
#[pyo3(signature = (smiles_list, seed=42, num_threads=0))]
fn embed_batch(smiles_list: Vec<String>, seed: u64, num_threads: usize) -> Vec<ConformerResult> {
    let refs: Vec<&str> = smiles_list.iter().map(|s| s.as_str()).collect();
    let config = sci_form_core::ConformerConfig { seed, num_threads };
    sci_form_core::embed_batch(&refs, &config)
        .into_iter()
        .map(|r| r.into())
        .collect()
}

#[pyfunction]
fn parse(py: Python<'_>, smiles: &str) -> PyResult<PyObject> {
    match sci_form_core::parse(smiles) {
        Ok(mol) => {
            let dict = PyDict::new_bound(py);
            let n = mol.graph.node_count();
            dict.set_item("num_atoms", n)?;
            dict.set_item("num_bonds", mol.graph.edge_count())?;
            let mut atoms_list: Vec<PyObject> = Vec::new();
            for i in 0..n {
                let idx = sci_form_core::graph::NodeIndex::new(i);
                let atom = &mol.graph[idx];
                let adict = PyDict::new_bound(py);
                adict.set_item("element", atom.element)?;
                adict.set_item("hybridization", format!("{:?}", atom.hybridization))?;
                adict.set_item("formal_charge", atom.formal_charge)?;
                atoms_list.push(adict.into());
            }
            dict.set_item("atoms", atoms_list)?;
            Ok(dict.into())
        }
        Err(e) => Err(pyo3::exceptions::PyValueError::new_err(e)),
    }
}

#[pyfunction]
fn version() -> String {
    sci_form_core::version()
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(embed, m)?)?;
    m.add_function(wrap_pyfunction!(embed_batch, m)?)?;
    m.add_function(wrap_pyfunction!(parse, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_class::<ConformerResult>()?;
    Ok(())
}
