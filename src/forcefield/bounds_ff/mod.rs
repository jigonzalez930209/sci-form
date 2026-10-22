//! Bounds-matrix force field used during distance-geometry embedding.
//!
//! This is the restraint field that pulls a randomly sampled metric matrix
//! back into the distance bounds derived from the molecular graph. It follows
//! RDKit's `BoundsForceField` design (`ForceFieldHelpers/..`): a sum of
//! distance-violation, chiral-volume and (in 4D) torsion terms whose gradients
//! drive a quasi-Newton minimiser.
//!
//! The submodules are:
//!
//! - [`energy`] — the bounds-violation energy and its analytical gradient,
//! - [`lbfgs`] — the L-BFGS minimiser used for the embedding refinement,
//! - [`bfgs`] — a BFGS variant kept for comparison and tests.

/// BFGS minimiser for the bounds force field (kept for comparison and tests).
pub mod bfgs;
/// Bounds-violation, chiral and 4D-torsion energy terms and their gradients.
pub mod energy;
/// L-BFGS minimiser for the bounds force field.
pub mod lbfgs;

pub use bfgs::*;
pub use energy::*;
pub use lbfgs::*;

/// A chiral constraint set matching RDKit's ChiralViolationContribs.
///
/// The set pins the signed volume of the tetrahedron formed by the stereocentre
/// and its three neighbours inside `[lower_vol, upper_vol]`, which forces the
/// embedding to reproduce the requested handedness.
pub struct ChiralSet {
    /// Index of the stereocentre (the apex of the tetrahedron).
    pub center: usize,
    /// Indices of the three (plus one padding) neighbours of the stereocentre.
    pub neighbors: [usize; 4],
    /// Lower bound on the signed chiral volume in Å³.
    pub lower_vol: f32,
    /// Upper bound on the signed chiral volume in Å³.
    pub upper_vol: f32,
}

/// Pre-computed torsion constraint for 4D embedding.
///
/// In 4D distance geometry the embedding carries an extra dimension, and a
/// torsion term keeps the projected 3D dihedral consistent with the requested
/// value: `E = w · (1 + cos(n_fold·φ - γ))` with `φ` the dihedral angle in
/// radians.
pub struct EmbedTorsion {
    /// The four atom indices `[i, j, k, l]` defining the dihedral.
    pub idx: [usize; 4],
    /// Energy amplitude `w` of the torsion restraint in kcal/mol.
    pub weight: f32,
    /// Periodicity `n_fold` of the restraint (dimensionless).
    pub n_fold: f32,
    /// Phase offset `γ` of the restraint in radians.
    pub gamma: f32,
}
