//! Riemannian Optimization for ETKDG — Track E3
//!
//! Replaces Euclidean BFGS with Riemannian L-BFGS over the manifold
//! of fixed-rank PSD matrices, eliminating negative eigenvalues by design.

mod lbfgs;
mod manifold;

pub use lbfgs::{DistanceConstraint, RiemannianConfig, RiemannianLbfgs, RiemannianResult};
pub use manifold::{psd_distance, psd_projection, psd_retraction, tangent_projection, PsdManifold};
