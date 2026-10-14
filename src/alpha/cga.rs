//! Backwards-compatible re-export of the core [`crate::cga`] module.
//!
//! CGA was promoted from alpha to core (F81 vs Rodrigues + `clifford`).
//! Prefer `sci_form::cga`; this path remains so existing `alpha::cga`
//! imports keep compiling.

pub use crate::cga::*;
