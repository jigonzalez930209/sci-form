//! Typed error type for the public API.
//!
//! [`SciFormError`] is the 2.0 error surface.  It classifies failures so that
//! callers can react to the *kind* of failure instead of matching substrings:
//!
//! ```rust
//! use sci_form::{parse, SciFormError};
//!
//! match parse("Q") {
//!     Ok(_) => unreachable!(),
//!     Err(error) if error.is_parse() => {
//!         assert!(!error.to_string().is_empty());
//!     }
//!     Err(other) => panic!("unexpected error kind: {other}"),
//! }
//! ```
//!
//! Internal modules still return `String` in a few legacy paths; the
//! [`From<String>`] implementation classifies those messages and keeps the
//! migration incremental without breaking `?` conversions.

use std::fmt;
use std::ops::Deref;

/// Classified library error.
///
/// The enum is `#[non_exhaustive]` so new categories can be added without a
/// breaking release.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum SciFormError {
    /// The input string could not be parsed (SMILES, CIF, JSON, …).
    Parse(String),
    /// The input was well-formed but failed validation.
    Validation(String),
    /// A molecular geometry is invalid (mismatched lengths, non-finite values).
    Geometry(String),
    /// A numerical method failed (SCF, EEM, decomposition, non-convergence).
    Numerical(String),
    /// The requested feature or element is not supported.
    Unsupported(String),
    /// Any other failure.
    Other(String),
}

impl SciFormError {
    /// Construct a parse error.
    pub fn parse(message: impl Into<String>) -> Self {
        Self::Parse(message.into())
    }

    /// Construct a validation error.
    pub fn validation(message: impl Into<String>) -> Self {
        Self::Validation(message.into())
    }

    /// Construct a geometry error.
    pub fn geometry(message: impl Into<String>) -> Self {
        Self::Geometry(message.into())
    }

    /// Construct a numerical error.
    pub fn numerical(message: impl Into<String>) -> Self {
        Self::Numerical(message.into())
    }

    /// Construct an unsupported-feature error.
    pub fn unsupported(message: impl Into<String>) -> Self {
        Self::Unsupported(message.into())
    }

    /// Construct a generic error.
    pub fn other(message: impl Into<String>) -> Self {
        Self::Other(message.into())
    }

    /// The human-readable message.
    pub fn message(&self) -> &str {
        match self {
            Self::Parse(message)
            | Self::Validation(message)
            | Self::Geometry(message)
            | Self::Numerical(message)
            | Self::Unsupported(message)
            | Self::Other(message) => message,
        }
    }

    /// True for [`SciFormError::Parse`].
    pub fn is_parse(&self) -> bool {
        matches!(self, Self::Parse(_))
    }

    /// True for [`SciFormError::Validation`].
    pub fn is_validation(&self) -> bool {
        matches!(self, Self::Validation(_))
    }

    /// True for [`SciFormError::Geometry`].
    pub fn is_geometry(&self) -> bool {
        matches!(self, Self::Geometry(_))
    }

    /// True for [`SciFormError::Numerical`].
    pub fn is_numerical(&self) -> bool {
        matches!(self, Self::Numerical(_))
    }

    /// True for [`SciFormError::Unsupported`].
    pub fn is_unsupported(&self) -> bool {
        matches!(self, Self::Unsupported(_))
    }
}

impl fmt::Display for SciFormError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.message())
    }
}

impl std::error::Error for SciFormError {}

impl Deref for SciFormError {
    type Target = str;

    fn deref(&self) -> &Self::Target {
        self.message()
    }
}

impl AsRef<str> for SciFormError {
    fn as_ref(&self) -> &str {
        self.message()
    }
}

impl From<String> for SciFormError {
    fn from(message: String) -> Self {
        classify(&message)
    }
}

impl From<&str> for SciFormError {
    fn from(message: &str) -> Self {
        classify(message)
    }
}

impl From<SciFormError> for String {
    fn from(error: SciFormError) -> Self {
        error.to_string()
    }
}

/// Best-effort classification of a legacy message.
fn classify(message: &str) -> SciFormError {
    let lower = message.to_ascii_lowercase();
    if message.starts_with("Expected '")
        || message.starts_with("Unexpected character")
        || message.starts_with("Unexpected end")
        || message.starts_with("Unknown element")
        || message.starts_with("Unknown bracket element")
        || message.starts_with("Invalid SMILES")
        || lower.contains("smiles")
        || lower.contains("parse")
        || lower.contains("invalid json")
    {
        SciFormError::Parse(message.to_string())
    } else if lower.contains("not supported")
        || lower.contains("unsupported")
        || lower.contains("no parameter")
        || lower.contains("no mmff")
    {
        SciFormError::Unsupported(message.to_string())
    } else if lower.contains("did not converge")
        || lower.contains("not converged")
        || lower.contains("failed to converge")
        || lower.contains("singular")
        || lower.contains("numerical")
    {
        SciFormError::Numerical(message.to_string())
    } else if lower.contains("coords length")
        || lower.contains("coordinates")
        || lower.contains("non-finite")
        || lower.contains("nan")
    {
        SciFormError::Geometry(message.to_string())
    } else if lower.contains("must be")
        || lower.contains("invalid")
        || lower.contains("requires")
        || lower.contains("mismatch")
    {
        SciFormError::Validation(message.to_string())
    } else {
        SciFormError::Other(message.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classifies_parse_messages() {
        let error = SciFormError::from("Unexpected end of SMILES input");
        assert!(error.is_parse());
        assert_eq!(error.message(), "Unexpected end of SMILES input");
    }

    #[test]
    fn classifies_numerical_and_validation() {
        assert!(SciFormError::from("SCF did not converge").is_numerical());
        assert!(SciFormError::from("coords length 4 != 3 * atoms 2").is_geometry());
        assert!(SciFormError::from("n_samples must be > 0").is_validation());
    }

    #[test]
    fn converts_back_to_string() {
        let error = SciFormError::validation("bad input");
        let message: String = error.into();
        assert_eq!(message, "bad input");
    }
}
