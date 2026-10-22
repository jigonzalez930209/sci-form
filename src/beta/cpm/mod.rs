//! Backwards-compatible re-export of the core [`crate::cpm`] module.
//!
//! CPM was promoted from beta to core (F95 vs LAMMPS `electrode/conp`).
//! Prefer `sci_form::cpm`; this path remains so existing `beta::cpm`
//! imports keep compiling.
//!
//! # Stability tier
//!
//! Core (promoted). This alias is gated behind the `beta-cpm` feature flag
//! (alias: `experimental-cpm`); the validation gate that covered the module
//! is `tests/experimental/test_cpm.rs` (F95).

pub use crate::cpm::grand_potential;
pub use crate::cpm::surface;
pub use crate::cpm::{compute_cpm_charges, compute_cpm_surface, CpmConfig, CpmResult, CpmSurface};
