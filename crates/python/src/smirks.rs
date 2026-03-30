//! SMIRKS reaction transform bindings for Python.

use pyo3::prelude::*;

/// A parsed SMIRKS reaction transform (Python-exposed).
#[pyclass]
#[derive(Clone)]
pub(crate) struct SmirksTransform {
    #[pyo3(get)]
    pub reactant_smarts: Vec<String>,
    #[pyo3(get)]
    pub product_smarts: Vec<String>,
    #[pyo3(get)]
    pub smirks: String,
}

impl From<sci_form_core::smirks::SmirksTransform> for SmirksTransform {
    fn from(t: sci_form_core::smirks::SmirksTransform) -> Self {
        SmirksTransform {
            reactant_smarts: t.reactant_smarts,
            product_smarts: t.product_smarts,
            smirks: t.smirks,
        }
    }
}

#[pymethods]
impl SmirksTransform {
    fn __repr__(&self) -> String {
        format!(
            "SmirksTransform(smirks='{}', reactants={}, products={})",
            self.smirks,
            self.reactant_smarts.len(),
            self.product_smarts.len()
        )
    }
}

/// Result of applying a SMIRKS transform (Python-exposed).
#[pyclass]
#[derive(Clone)]
pub(crate) struct SmirksResult {
    #[pyo3(get)]
    pub products: Vec<String>,
    #[pyo3(get)]
    pub n_transforms: usize,
    #[pyo3(get)]
    pub success: bool,
    #[pyo3(get)]
    pub messages: Vec<String>,
}

impl From<sci_form_core::smirks::SmirksResult> for SmirksResult {
    fn from(r: sci_form_core::smirks::SmirksResult) -> Self {
        SmirksResult {
            products: r.products,
            n_transforms: r.n_transforms,
            success: r.success,
            messages: r.messages,
        }
    }
}

#[pymethods]
impl SmirksResult {
    fn __repr__(&self) -> String {
        if self.success {
            format!(
                "SmirksResult(success=True, products={}, n_transforms={})",
                self.products.len(),
                self.n_transforms
            )
        } else {
            format!(
                "SmirksResult(success=False, messages={:?})",
                self.messages
            )
        }
    }

    fn is_ok(&self) -> bool {
        self.success
    }
}

/// Parse a SMIRKS reaction string.
///
/// Args:
///     smirks: SMIRKS pattern in format "reactant>>product" with atom maps
///
/// Returns:
///     SmirksTransform object containing parsed reactant/product patterns
///
/// Example:
///     >>> transform = parse_smirks("[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]")
///     >>> print(transform.reactant_smarts)
///     ['[C:1](=O)[OH:2]']
#[pyfunction]
fn parse_smirks(smirks: &str) -> PyResult<SmirksTransform> {
    match sci_form_core::smirks::parse_smirks(smirks) {
        Ok(transform) => Ok(transform.into()),
        Err(e) => Err(pyo3::exceptions::PyValueError::new_err(e)),
    }
}

/// Apply a SMIRKS transform to a molecule.
///
/// Args:
///     smirks: SMIRKS pattern string
///     smiles: Input molecule SMILES string
///
/// Returns:
///     SmirksResult with success status and product SMARTS patterns
///
/// Example:
///     >>> result = apply_smirks("[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]", "CC(=O)O")
///     >>> print(result.success)
///     True
///     >>> print(result.n_transforms)
///     1
#[pyfunction]
fn apply_smirks(smirks: &str, smiles: &str) -> PyResult<SmirksResult> {
    match sci_form_core::smirks::apply_smirks(smirks, smiles) {
        Ok(result) => Ok(result.into()),
        Err(e) => Err(pyo3::exceptions::PyValueError::new_err(e)),
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse_smirks, m)?)?;
    m.add_function(wrap_pyfunction!(apply_smirks, m)?)?;
    m.add_class::<SmirksTransform>()?;
    m.add_class::<SmirksResult>()?;
    Ok(())
}
