//! Optional keyword-argument ("options object") support for the Python
//! entry points that take many parameters (blocker C-M3).
//!
//! The Rust core already groups the parameters of the heavy entry points
//! into configuration structs (`HfConfig`, `AlpbConfig`, `GasteigerConfig`,
//! `FluorescenceConfig`, `ScfConfig`, `UhfConfig`).  The Python surface used
//! to flatten them into long positional signatures, which makes a call like
//! `optimize_framework(...)` with 14 arguments unreadable.
//!
//! Every affected function therefore accepts an **optional** trailing
//! `options` mapping.  The positional form keeps working exactly as before;
//! a key present in the mapping overrides the positional value.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

/// Read a `float` from the options mapping.
pub(crate) fn f64_from(
    opts: Option<&Bound<'_, PyDict>>,
    key: &str,
    fallback: f64,
) -> PyResult<f64> {
    match opts.and_then(|d| d.get_item(key).ok().flatten()) {
        None => Ok(fallback),
        Some(value) => value.extract::<f64>().map_err(|e| {
            PyValueError::new_err(format!("Invalid 'options' value for '{key}': {e}"))
        }),
    }
}

/// Read a non-negative `int` from the options mapping.
pub(crate) fn usize_from(
    opts: Option<&Bound<'_, PyDict>>,
    key: &str,
    fallback: usize,
) -> PyResult<usize> {
    match opts.and_then(|d| d.get_item(key).ok().flatten()) {
        None => Ok(fallback),
        Some(value) => {
            let v = value.extract::<i64>().map_err(|e| {
                PyValueError::new_err(format!("Invalid 'options' value for '{key}': {e}"))
            })?;
            if v < 0 {
                return Err(PyValueError::new_err(format!(
                    "Invalid 'options' value for '{key}': {v}. Expected a non-negative integer."
                )));
            }
            Ok(v as usize)
        }
    }
}

/// Read a `str` from the options mapping.
pub(crate) fn str_from(
    opts: Option<&Bound<'_, PyDict>>,
    key: &str,
    fallback: &str,
) -> PyResult<String> {
    match opts.and_then(|d| d.get_item(key).ok().flatten()) {
        None => Ok(fallback.to_string()),
        Some(value) => value.extract::<String>().map_err(|e| {
            PyValueError::new_err(format!("Invalid 'options' value for '{key}': {e}"))
        }),
    }
}

/// Read an optional `float` from the options mapping.
pub(crate) fn opt_f64_from(
    opts: Option<&Bound<'_, PyDict>>,
    key: &str,
    fallback: Option<f64>,
) -> PyResult<Option<f64>> {
    match opts.and_then(|d| d.get_item(key).ok().flatten()) {
        None => Ok(fallback),
        Some(value) => {
            if value.is_none() {
                return Ok(None);
            }
            value.extract::<f64>().map(Some).map_err(|e| {
                PyValueError::new_err(format!("Invalid 'options' value for '{key}': {e}"))
            })
        }
    }
}

/// Read a list of `int` from the options mapping.
pub(crate) fn usize_list_from(
    opts: Option<&Bound<'_, PyDict>>,
    key: &str,
    fallback: Option<Vec<usize>>,
) -> PyResult<Option<Vec<usize>>> {
    match opts.and_then(|d| d.get_item(key).ok().flatten()) {
        None => Ok(fallback),
        Some(value) => {
            if value.is_none() {
                return Ok(None);
            }
            value.extract::<Vec<usize>>().map(Some).map_err(|e| {
                PyValueError::new_err(format!("Invalid 'options' value for '{key}': {e}"))
            })
        }
    }
}

/// Reject option keys that the function does not know about.
///
/// A typo in an options mapping is a silent no-op otherwise, which is exactly
/// the failure mode the options object is meant to remove.
pub(crate) fn reject_unknown(opts: Option<&Bound<'_, PyDict>>, known: &[&str]) -> PyResult<()> {
    let Some(opts) = opts else { return Ok(()) };
    for key in opts.keys() {
        let key = key.extract::<String>()?;
        if !known.contains(&key.as_str()) {
            let list = known
                .iter()
                .map(|k| format!("'{k}'"))
                .collect::<Vec<_>>()
                .join(", ");
            return Err(PyValueError::new_err(format!(
                "Invalid 'options' key '{key}'. Supported: {list}."
            )));
        }
    }
    Ok(())
}
