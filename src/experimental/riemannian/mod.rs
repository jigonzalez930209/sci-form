//! Riemannian Optimization for ETKDG — Track E3
//!
//! Replaces Euclidean BFGS with Riemannian L-BFGS over the manifold
//! of fixed-rank PSD matrices, eliminating negative eigenvalues by design.

mod manifold;
mod lbfgs;

pub use manifold::{PsdManifold, psd_retraction, psd_projection, tangent_projection};
pub use lbfgs::{RiemannianLbfgs, RiemannianConfig, RiemannianResult, DistanceConstraint};
