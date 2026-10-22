//! Backwards-compatible re-export of the core [`crate::cga`] module.
//!
//! CGA was promoted from alpha to core (F81 vs Rodrigues + `clifford`).
//! Prefer `sci_form::cga`; this path remains so existing `alpha::cga`
//! imports keep compiling.
//!
//! # Stability tier
//!
//! Core (promoted). This alias is gated behind the `alpha-cga` feature flag
//! (alias: `experimental-cga`); the validation gate that covered the module
//! is `tests/alpha_cga_validation.rs` (F81).

pub use crate::cga::*;
