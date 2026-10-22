//! Backwards-compatible re-export of the core [`crate::riemannian`] module.
//!
//! Riemannian PSD geometry was promoted from beta to core (F94 vs numpy +
//! scipy). Prefer `sci_form::riemannian`; this path remains so existing
//! `beta::riemannian` imports keep compiling.
//!
//! Documented gap: [`crate::riemannian::psd_distance`] is the Frobenius
//! proxy, not the affine-invariant geodesic distance.
//!
//! # Stability tier
//!
//! Core (promoted). This alias is gated behind the `beta-riemannian` feature
//! flag (alias: `experimental-riemannian`); the validation gate that covered
//! the module is `tests/beta_riemannian_validation.rs` (F94).

pub use crate::riemannian::*;
