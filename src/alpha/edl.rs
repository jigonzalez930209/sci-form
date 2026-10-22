//! Backwards-compatible re-export of the core [`crate::edl`] module.
//!
//! EDL was promoted from alpha to core (F86 vs analytic Grahame/Stern).
//! Prefer `sci_form::edl`; this path remains so existing `alpha::edl`
//! imports keep compiling.
//!
//! # Stability tier
//!
//! Core (promoted). This alias is gated behind the `alpha-edl` feature flag
//! (which also pulls in `experimental-cpm`, `experimental-eeq` and
//! `experimental-alpb`); the validation gate that covered the module is
//! `tests/alpha_edl_validation.rs` (F86).

pub use crate::edl::*;
