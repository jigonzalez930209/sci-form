//! Error taxonomy and canonical error-message helpers for the Python bindings.
//!
//! # Taxonomy (blocker C-M6)
//!
//! | Situation | Exception |
//! |-----------|-----------|
//! | The caller's input is malformed (bad SMILES, wrong coordinate count, unknown method) | `ValueError` |
//! | The engine ran but failed (SCF blew up, geometry optimisation failed) | `RuntimeError` |
//! | The engine ran but did not converge in the allowed number of iterations | `NotConvergedError` (a `RuntimeError`) |
//!
//! `NotConvergedError` is exported from the package as
//! `sci_form.NotConvergedError` so a caller can distinguish "the physics did
//! not settle" from "you passed nonsense".
//!
//! # Message style (blocker C-H9)
//!
//! Every message follows
//!
//! ```text
//! <What failed> '<actual value>'. Supported: 'a', 'b', 'c'.
//! ```
//!
//! All literals are quoted, the offending input is always included and the
//! valid values are always listed when there is a fixed set.  Rust internal
//! names (`positions`, `unwrap`, `None`, a `file:line`) never leak into a
//! Python message.

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;

pyo3::create_exception!(
    sci_form,
    NotConvergedError,
    PyRuntimeError,
    "An iterative engine (SCF, geometry optimisation, CPM) did not converge. \
     Subclasses RuntimeError, so `except RuntimeError` still catches it."
);

/// Mark a message as "the caller's input is wrong".
pub(crate) fn value_error(message: impl Into<String>) -> PyErr {
    PyValueError::new_err(message.into())
}

/// Mark a message as "the engine failed".
pub(crate) fn runtime_error(message: impl Into<String>) -> PyErr {
    PyRuntimeError::new_err(message.into())
}

/// Mark a message as "the engine did not converge".
pub(crate) fn not_converged(message: impl Into<String>) -> PyErr {
    NotConvergedError::new_err(message.into())
}

/// Build the canonical "invalid value" error message.
///
/// `<what> '<value>'. Supported: 'a', 'b', 'c'.`
pub(crate) fn invalid_value(what: &str, value: &str, supported: &[&str]) -> String {
    if supported.is_empty() {
        return format!("Invalid {what} '{value}'.");
    }
    let list = supported
        .iter()
        .map(|s| format!("'{s}'"))
        .collect::<Vec<_>>()
        .join(", ");
    format!("Invalid {what} '{value}'. Supported: {list}.")
}

/// Truncate a long string for inclusion in an error message.
pub(crate) fn truncate(value: &str, max_chars: usize) -> String {
    let trimmed = value.trim();
    if trimmed.chars().count() <= max_chars {
        return trimmed.to_string();
    }
    let head: String = trimmed.chars().take(max_chars).collect();
    format!("{head}…")
}

/// Position of the last unclosed `open` bracket in `text`, if any.
fn unclosed_position(text: &str, open: u8, close: u8) -> Option<usize> {
    let mut stack: Vec<usize> = Vec::new();
    for (i, b) in text.bytes().enumerate() {
        if b == open {
            stack.push(i);
        } else if b == close {
            stack.pop();
        }
    }
    stack.last().copied()
}

/// Enrich a SMILES parser error with the offending input and the position.
///
/// The core parser reports bare messages such as `Expected ')'`.  The
/// bindings rewrite them into the canonical style and add the position of
/// the offending bracket — the single most common user error in the library
/// (blocker C-H9).
pub(crate) fn smiles_error(message: &str, smiles: &str) -> String {
    let fragment = truncate(smiles, 120);
    let detail = match message {
        "Expected ')'" => match unclosed_position(smiles, b'(', b')') {
            Some(pos) => format!("expected ')' to close the '(' at position {pos}"),
            None => "expected ')'".to_string(),
        },
        "Expected ']'" => match unclosed_position(smiles, b'[', b']') {
            Some(pos) => format!("expected ']' to close the '[' at position {pos}"),
            None => "expected ']'".to_string(),
        },
        "Unexpected end of SMILES input" => "unexpected end of input".to_string(),
        other => other.to_string(),
    };
    format!("Invalid 'smiles' '{fragment}'. {detail}")
}

/// Classify a core error string and turn it into the right Python exception.
///
/// * A message that reports a non-convergence becomes [`NotConvergedError`].
/// * Everything else becomes a plain [`RuntimeError`].
///
/// The caller decides separately whether the failure was caused by malformed
/// input; that check needs the *arguments*, not the message.
pub(crate) fn engine_error(message: impl Into<String>) -> PyErr {
    let message = message.into();
    let lower = message.to_ascii_lowercase();
    if lower.contains("did not converge")
        || lower.contains("not converged")
        || lower.contains("failed to converge")
        || lower.contains("convergence criterion")
        || lower.contains("max_iterations")
        || lower.contains("maximum number of iterations")
    {
        not_converged(message)
    } else {
        runtime_error(message)
    }
}

/// Canonical message for a coordinate array whose length is wrong.
pub(crate) fn coords_length_error(what: &str, got: usize, expected: usize) -> String {
    format!("Invalid {what} length {got}. Expected {expected} values (3 per atom).")
}

/// Canonical message for a coordinate array that is not a multiple of three.
pub(crate) fn coords_triplet_error(what: &str, got: usize) -> String {
    format!("Invalid {what} length {got}. Expected a multiple of 3 (x, y, z triplets).")
}

/// Validate that a flat coordinate array has `3 * n_atoms` entries.
///
/// Returns a `ValueError` describing the mismatch.
pub(crate) fn check_coords_len(coords: &[f64], n_atoms: usize) -> Result<(), PyErr> {
    let expected = n_atoms * 3;
    if coords.len() != expected {
        return Err(value_error(coords_length_error(
            "'coords'",
            coords.len(),
            expected,
        )));
    }
    if !coords.len().is_multiple_of(3) {
        return Err(value_error(coords_triplet_error("'coords'", coords.len())));
    }
    Ok(())
}

/// Canonical message for a partial-charge array whose length is wrong.
pub(crate) fn charges_length_error(got: usize, expected: usize) -> String {
    format!("Invalid 'charges' length {got}. Expected {expected} values (one per atom).")
}

/// Validate that a partial-charge array has one entry per atom.
///
/// Returns a `ValueError` describing the mismatch.
pub(crate) fn check_charges_len(charges: &[f64], n_atoms: usize) -> Result<(), PyErr> {
    if charges.len() != n_atoms {
        return Err(value_error(charges_length_error(charges.len(), n_atoms)));
    }
    Ok(())
}

/// Whether a core error message comes from the SMILES parser.
///
/// The core reports a handful of bare messages (`Expected ')'`, `Unknown
/// element: …`, …); they are all *malformed input*, not engine failures.
pub(crate) fn is_parse_error(message: &str) -> bool {
    message.starts_with("Expected '")
        || message.starts_with("Unexpected character")
        || message.starts_with("Unexpected end of SMILES")
        || message.starts_with("Unknown element")
        || message.starts_with("Unknown bracket element")
        || message.starts_with("Invalid SMILES")
        || message.contains("SMILES")
}

/// Classify an error raised by a function that takes a SMILES string.
///
/// A parse failure becomes a `ValueError` enriched with the position of the
/// offending bracket; anything else is an engine failure.
pub(crate) fn smiles_or_engine_error(message: impl Into<String>, smiles: &str) -> PyErr {
    let message = message.into();
    if is_parse_error(&message) {
        value_error(smiles_error(&message, smiles))
    } else {
        engine_error(message)
    }
}

/// Register the error taxonomy on the package.
pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("NotConvergedError", m.py().get_type::<NotConvergedError>())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invalid_value_lists_the_supported_values() {
        let msg = invalid_value("'method'", "foo", &["eht", "pm3"]);
        assert_eq!(msg, "Invalid 'method' 'foo'. Supported: 'eht', 'pm3'.");
    }

    #[test]
    fn smiles_error_adds_the_position() {
        let msg = smiles_error("Expected ')'", "C(CC");
        assert!(msg.contains("position 1"), "got: {msg}");
        assert!(msg.contains("C(CC"), "got: {msg}");
    }

    #[test]
    fn convergence_failures_are_detected() {
        Python::attach(|py| {
            let err = engine_error("PM3 did not converge at this geometry");
            assert!(err.is_instance_of::<NotConvergedError>(py));
            let err = engine_error("something else went wrong");
            assert!(err.is_instance_of::<PyRuntimeError>(py));
            assert!(!err.is_instance_of::<NotConvergedError>(py));
        });
    }
}
