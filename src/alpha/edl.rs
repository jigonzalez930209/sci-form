//! Backwards-compatible re-export of the core [`crate::edl`] module.
//!
//! EDL was promoted from alpha to core (F86 vs analytic Grahame/Stern).
//! Prefer `sci_form::edl`; this path remains so existing `alpha::edl`
//! imports keep compiling.

pub use crate::edl::*;
