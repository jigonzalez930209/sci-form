//! Riemannian Optimization for ETKDG — core since 2026-10-14 (F94).
//!
//! Replaces Euclidean BFGS with Riemannian L-BFGS over the manifold
//! of fixed-rank PSD matrices, eliminating negative eigenvalues by design.
//! Validated vs numpy + scipy L-BFGS-B. Documented gap: `psd_distance` is
//! the Frobenius proxy, not the affine-invariant geodesic distance.

mod lbfgs;
mod manifold;

pub use lbfgs::{DistanceConstraint, RiemannianConfig, RiemannianLbfgs, RiemannianResult};
pub use manifold::{psd_distance, psd_projection, psd_retraction, tangent_projection, PsdManifold};
