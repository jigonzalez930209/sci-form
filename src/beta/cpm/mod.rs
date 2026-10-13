//! Backwards-compatible re-export of the core [`crate::cpm`] module.
//!
//! CPM was promoted from beta to core (F95 vs LAMMPS `electrode/conp`).
//! Prefer `sci_form::cpm`; this path remains so existing `beta::cpm`
//! imports keep compiling.

pub use crate::cpm::grand_potential;
pub use crate::cpm::surface;
pub use crate::cpm::{
    compute_cpm_charges, compute_cpm_surface, CpmConfig, CpmResult, CpmSurface,
};
