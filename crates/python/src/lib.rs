use pyo3::prelude::*;
use pyo3::types::PyDict;

/// A single conformer result returned to Python.
#[pyclass]
#[derive(Clone)]
struct ConformerResult {
    #[pyo3(get)]
    smiles: String,
    #[pyo3(get)]
    num_atoms: usize,
    #[pyo3(get)]
    coords: Vec<f64>,
    #[pyo3(get)]
    elements: Vec<u8>,
    #[pyo3(get)]
    bonds: Vec<(usize, usize, String)>,
    #[pyo3(get)]
    error: Option<String>,
    #[pyo3(get)]
    time_ms: f64,
}

#[pymethods]
impl ConformerResult {
    /// Get coordinates as list of (x, y, z) tuples.
    fn get_positions(&self) -> Vec<(f64, f64, f64)> {
        self.coords
            .chunks_exact(3)
            .map(|c| (c[0], c[1], c[2]))
            .collect()
    }

    /// Check if generation succeeded.
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

/// Generate a 3D conformer from a SMILES string.
///
/// Args:
///     smiles: SMILES string of the molecule.
///     seed: RNG seed for reproducibility (default: 42).
///
/// Returns:
///     ConformerResult with 3D coordinates.
#[pyfunction]
#[pyo3(signature = (smiles, seed=42))]
fn embed(smiles: &str, seed: u64) -> ConformerResult {
    sci_form_core::embed(smiles, seed).into()
}

/// Batch-generate 3D conformers for multiple SMILES in parallel.
///
/// Args:
///     smiles_list: List of SMILES strings.
///     seed: RNG seed (default: 42).
///     num_threads: Number of threads, 0 = auto (default: 0).
///
/// Returns:
///     List of ConformerResult objects.
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

/// Parse a SMILES string without generating 3D coordinates.
///
/// Returns a dict with atom/bond info.
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

/// Get the library version string.
#[pyfunction]
fn version() -> String {
    sci_form_core::version()
}

/// sci_form Python module
#[pymodule]
fn sci_form(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(embed, m)?)?;
    m.add_function(wrap_pyfunction!(embed_batch, m)?)?;
    m.add_function(wrap_pyfunction!(parse, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_class::<ConformerResult>()?;
    Ok(())
}
