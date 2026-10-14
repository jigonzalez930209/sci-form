//! Backwards-compatible re-export of the core [`crate::riemannian`] module.
//!
//! Riemannian PSD geometry was promoted from beta to core (F94 vs numpy +
//! scipy). Prefer `sci_form::riemannian`; this path remains so existing
//! `beta::riemannian` imports keep compiling.
//!
//! Documented gap: [`crate::riemannian::psd_distance`] is the Frobenius
//! proxy, not the affine-invariant geodesic distance.

pub use crate::riemannian::*;
