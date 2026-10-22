//! Conformer embedding bindings.

use pyo3::prelude::*;

use crate::buffer::F64Buffer;
use crate::errors::{engine_error, smiles_error, value_error};

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

/// One atom of a parsed molecule.
#[pyclass(name = "ParsedAtom")]
#[derive(Clone)]
pub(crate) struct ParsedAtom {
    /// Atomic number.
    #[pyo3(get)]
    pub element: u8,
    /// Hybridization label, e.g. `"SP3"`.
    #[pyo3(get)]
    pub hybridization: String,
    /// Formal charge.
    #[pyo3(get)]
    pub formal_charge: i8,
    /// Number of explicit hydrogens.
    #[pyo3(get)]
    pub explicit_h: u8,
}

#[pymethods]
impl ParsedAtom {
    fn __repr__(&self) -> String {
        format!(
            "ParsedAtom(element={}, hybridization='{}', formal_charge={})",
            self.element, self.hybridization, self.formal_charge
        )
    }
}

/// The result of [`parse`].
///
/// A proper `pyclass` with the same fields the previous plain `dict`
/// exposed, so `result.num_atoms`, `result.num_bonds` and `result.atoms`
/// keep working while `result["num_atoms"]` now raises a `TypeError`
/// instead of silently succeeding (blocker C-H5).
#[pyclass(name = "ParseResult")]
#[derive(Clone)]
pub(crate) struct ParseResult {
    /// Number of atoms (implicit hydrogens included).
    #[pyo3(get)]
    pub num_atoms: usize,
    /// Number of bonds.
    #[pyo3(get)]
    pub num_bonds: usize,
    /// The atoms of the molecule.
    #[pyo3(get)]
    pub atoms: Vec<ParsedAtom>,
}

#[pymethods]
impl ParseResult {
    fn __repr__(&self) -> String {
        format!(
            "ParseResult(num_atoms={}, num_bonds={})",
            self.num_atoms, self.num_bonds
        )
    }
    fn __len__(&self) -> usize {
        self.num_atoms
    }
    /// Keep the mapping-style access of the previous plain-`dict` payload.
    ///
    /// `result["num_atoms"]` still works, so scripts written against the
    /// dictionary shape do not break.
    fn __getitem__(&self, py: Python<'_>, key: &str) -> PyResult<Py<pyo3::PyAny>> {
        match key {
            "num_atoms" => Ok(self.num_atoms.into_pyobject(py)?.to_owned().into()),
            "num_bonds" => Ok(self.num_bonds.into_pyobject(py)?.to_owned().into()),
            "atoms" => Ok(self.atoms.clone().into_pyobject(py)?.to_owned().into()),
            other => Err(pyo3::exceptions::PyKeyError::new_err(format!(
                "'{other}' is not a key of ParseResult. Supported: 'num_atoms', 'num_bonds', 'atoms'."
            ))),
        }
    }
    /// The keys of the payload, so `dict(result)` keeps working.
    fn keys(&self) -> Vec<&'static str> {
        vec!["num_atoms", "num_bonds", "atoms"]
    }
}

/// Generate a 3D conformer for a SMILES string.
///
/// Raises `ValueError` when the SMILES is malformed and `RuntimeError` when
/// the embedding engine fails; a failed embed never returns a result object
/// with `error` set, so the failure cannot be silently ignored
/// (blocker C-M4).
#[pyfunction]
#[pyo3(signature = (smiles, seed=42))]
fn embed(smiles: &str, seed: u64) -> PyResult<ConformerResult> {
    let result = sci_form_core::embed(smiles, seed);
    if let Some(ref e) = result.error {
        return Err(crate::errors::smiles_or_engine_error(
            format!(
                "Embedding failed for SMILES '{}': {}",
                crate::errors::truncate(smiles, 120),
                e
            ),
            smiles,
        ));
    }
    Ok(result.into())
}

/// Batch-embed multiple SMILES strings.
///
/// The batch itself either succeeds or raises; a *per-molecule* failure is
/// reported through the `error` field of the corresponding
/// [`ConformerResult`], because a batch is expected to contain failures.
#[pyfunction]
#[pyo3(signature = (smiles_list, seed=42, num_threads=0))]
fn embed_batch(
    smiles_list: Vec<String>,
    seed: u64,
    num_threads: usize,
) -> PyResult<Vec<ConformerResult>> {
    let refs: Vec<&str> = smiles_list.iter().map(|s| s.as_str()).collect();
    let config = sci_form_core::ConformerConfig { seed, num_threads };
    sci_form_core::embed_batch(&refs, &config)
        .map_err(engine_error)
        .map(|results| results.into_iter().map(|r| r.into()).collect())
}

/// Generate a diverse conformer ensemble and cluster it by Butina RMSD.
///
/// Embeds `n_conformers` times with consecutive seeds, then keeps only the
/// cluster centroids (RMSD threshold `rmsd_cutoff` in Å).
#[pyfunction]
#[pyo3(signature = (smiles, n_conformers=10, rmsd_cutoff=1.0, base_seed=42))]
fn embed_diverse(
    smiles: &str,
    n_conformers: usize,
    rmsd_cutoff: f64,
    base_seed: u64,
) -> Vec<ConformerResult> {
    sci_form_core::embed_diverse(smiles, n_conformers, rmsd_cutoff, base_seed)
        .into_iter()
        .map(|r| r.into())
        .collect()
}

/// Return the canonical SMILES representation of a molecule.
///
/// The output is deterministic regardless of the input atom ordering.
#[pyfunction]
fn to_canonical_smiles(smiles: &str) -> PyResult<String> {
    sci_form_core::to_canonical_smiles(smiles).map_err(|e| value_error(smiles_error(&e, smiles)))
}

/// Parse a SMILES string into a [`ParseResult`].
///
/// Raises `ValueError` with the position of the offending bracket when the
/// SMILES is malformed.
#[pyfunction]
fn parse(smiles: &str) -> PyResult<ParseResult> {
    let mol = sci_form_core::parse(smiles).map_err(|e| value_error(smiles_error(&e, smiles)))?;
    let n = mol.graph.node_count();
    let mut atoms = Vec::with_capacity(n);
    for i in 0..n {
        let idx = sci_form_core::graph::NodeIndex::new(i);
        let atom = &mol.graph[idx];
        atoms.push(ParsedAtom {
            element: atom.element,
            hybridization: format!("{:?}", atom.hybridization),
            formal_charge: atom.formal_charge,
            explicit_h: atom.explicit_h,
        });
    }
    Ok(ParseResult {
        num_atoms: n,
        num_bonds: mol.graph.edge_count(),
        atoms,
    })
}

/// The library version.
#[pyfunction]
fn version() -> String {
    sci_form_core::version()
}

/// Embed a SMILES and return the coordinates **without copying**.
///
/// `coords` is a writable, C-contiguous `float64` buffer of length
/// `3 * num_atoms`.  The coordinates are written straight into the buffer
/// the caller owns, so a full volumetric grid or a long MD trajectory never
/// crosses the Python boundary twice (blocker C-H6).
///
/// The buffer must be exactly `3 * num_atoms` long; otherwise a
/// `ValueError` is raised and the buffer is left untouched.
///
/// ```python
/// import numpy as np, sci_form
/// buf = np.zeros(3 * sci_form.parse("CCO").num_atoms, dtype=np.float64)
/// sci_form.embed_coords_into("CCO", buf, seed=42)
/// ```
#[pyfunction]
#[pyo3(signature = (smiles, coords, seed=42))]
fn embed_coords_into(smiles: &str, coords: &Bound<'_, PyAny>, seed: u64) -> PyResult<usize> {
    let result = sci_form_core::embed(smiles, seed);
    if let Some(ref e) = result.error {
        return Err(crate::errors::smiles_or_engine_error(
            format!(
                "Embedding failed for SMILES '{}': {}",
                crate::errors::truncate(smiles, 120),
                e
            ),
            smiles,
        ));
    }
    let expected = result.num_atoms * 3;
    let mut buffer = F64Buffer::writable_f64(coords)?;
    if buffer.len() != expected {
        return Err(value_error(format!(
            "Invalid 'coords' buffer length {}. Expected {expected} values (3 per atom).",
            buffer.len()
        )));
    }
    buffer.as_mut_slice().copy_from_slice(&result.coords);
    Ok(result.num_atoms)
}

/// Compute the electrostatic potential on a 3D grid **without copying**.
///
/// `out` is a writable, C-contiguous `float64` buffer of length
/// `dims[0] * dims[1] * dims[2]`.  The grid values are written straight
/// into the buffer the caller owns (blocker C-H6).
///
/// `spacing` and `padding` are in Å; the grid is centred on the molecule
/// and spans `padding` Å beyond the bounding box on every side.
#[pyfunction]
#[pyo3(signature = (elements, coords, out, spacing=0.5, padding=3.0))]
fn compute_esp_grid_into(
    elements: Vec<u8>,
    coords: Vec<f64>,
    out: &Bound<'_, PyAny>,
    spacing: f64,
    padding: f64,
) -> PyResult<usize> {
    let positions = crate::system::coords_to_positions(&coords);
    crate::errors::check_coords_len(&coords, elements.len())?;
    let grid = sci_form_core::compute_esp(&elements, &positions, spacing, padding)
        .map_err(engine_error)?;
    let expected = grid.dims[0] * grid.dims[1] * grid.dims[2];
    let mut buffer = F64Buffer::writable_f64(out)?;
    if buffer.len() != expected {
        return Err(value_error(format!(
            "Invalid 'out' buffer length {}. Expected {expected} values (one per grid point).",
            buffer.len()
        )));
    }
    buffer.as_mut_slice().copy_from_slice(&grid.values);
    Ok(expected)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(embed, m)?)?;
    m.add_function(wrap_pyfunction!(embed_batch, m)?)?;
    m.add_function(wrap_pyfunction!(embed_diverse, m)?)?;
    m.add_function(wrap_pyfunction!(to_canonical_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(parse, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(embed_coords_into, m)?)?;
    m.add_function(wrap_pyfunction!(compute_esp_grid_into, m)?)?;
    m.add_class::<ConformerResult>()?;
    m.add_class::<ParsedAtom>()?;
    m.add_class::<ParseResult>()?;
    Ok(())
}
