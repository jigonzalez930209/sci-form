//! Backwards-compatible re-export of the core [`crate::sdr`] module.
//!
//! SDR was promoted from alpha to core (F83 vs numpy + CVXPY).
//! Prefer `sci_form::sdr`; this path remains so existing `alpha::sdr`
//! imports keep compiling.

pub use crate::sdr::*;
