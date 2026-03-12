//! Bounds distance force field for initial embedding optimization.
//! Includes bounds violation energy/gradient, chiral enforcement, and BFGS optimizer.

pub mod energy;
pub mod lbfgs;
pub mod bfgs;

pub use energy::*;
pub use lbfgs::*;
pub use bfgs::*;

/// A chiral constraint set matching RDKit's ChiralViolationContribs.
pub struct ChiralSet {
    pub center: usize,
    pub neighbors: [usize; 4],
    pub lower_vol: f32,
    pub upper_vol: f32,
}

/// Pre-computed torsion constraint for 4D embedding.
pub struct EmbedTorsion {
    pub idx: [usize; 4],
    pub weight: f32,
    pub n_fold: f32,
    pub gamma: f32,
}
