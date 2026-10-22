//! Backwards-compatible re-export of the core [`crate::sdr`] module.
//!
//! SDR was promoted from alpha to core (F83 vs numpy + CVXPY).
//! Prefer `sci_form::sdr`; this path remains so existing `alpha::sdr`
//! imports keep compiling.
//!
//! # Stability tier
//!
//! Core (promoted). This alias is gated behind the `alpha-sdr` feature flag
//! (alias: `experimental-sdr`); the validation gate that covered the module
//! is `tests/alpha_sdr_validation.rs` (F83).

pub use crate::sdr::*;
