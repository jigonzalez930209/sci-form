//! ETKDG 3D force field matching RDKit's `construct3DForceField`.
//!
//! Unlike a plain UFF-style field, this one uses **flat-bottom distance
//! constraints**: the penalty is zero while the distance stays inside
//! `[min_len, max_len]` and quadratic outside it, instead of a harmonic bond
//! stretch plus angle bend. That mirrors RDKit's
//! `ForceFieldHelpers/UFF/...` ETKDG variant and keeps the minimisation
//! well-conditioned for strained geometries.
//!
//! The submodules are:
//!
//! - [`builder`] — assembles the term lists from a molecule and a bounds matrix,
//! - [`energy`] — the total energy and its analytical gradient,
//! - [`gradient`] — gradient helpers for the individual term families,
//! - [`optimizer`] — the BFGS minimiser used for the 3D refinement.

pub mod builder;
/// Total energy and analytical gradient of the ETKDG 3D force field.
pub mod energy;
/// Analytical gradient helpers for the individual term families.
pub mod gradient;
/// BFGS minimiser for the ETKDG 3D refinement stage.
pub mod optimizer;

pub use builder::*;
pub use energy::*;
pub use gradient::*;
pub use optimizer::*;

const KNOWN_DIST_TOL: f64 = 0.01;
const KNOWN_DIST_K: f64 = 100.0;

/// A pre-computed flat-bottom distance constraint.
///
/// The penalty is `½·k·(d - min_len)²` for `d < min_len` and
/// `½·k·(d - max_len)²` for `d > max_len`, and zero in between.
#[derive(Clone)]
pub struct DistConstraint {
    /// Index of the first atom of the constrained pair.
    pub i: usize,
    /// Index of the second atom of the constrained pair.
    pub j: usize,
    /// Lower distance limit in Å (inclusive, no penalty).
    pub min_len: f64,
    /// Upper distance limit in Å (inclusive, no penalty).
    pub max_len: f64,
    /// Force constant in kcal/mol/Å².
    pub k: f64,
}

/// Pre-computed M6 torsion contribution (Fourier series).
///
/// The energy is `Σ_n v_n·(1 + sign_n·cos(n·φ))` over `n = 1…6`, with `φ`
/// the dihedral angle in radians.
#[derive(Clone)]
pub struct M6TorsionContrib {
    /// Index of the first terminal atom of the dihedral.
    pub i: usize,
    /// Index of the second atom of the dihedral.
    pub j: usize,
    /// Index of the third atom of the dihedral.
    pub k: usize,
    /// Index of the fourth atom of the dihedral.
    pub l: usize,
    /// Phase signs `s₁ … s₆` of the six cosine terms (dimensionless, each `±1`).
    pub signs: [f64; 6],
    /// Amplitudes `V₁ … V₆` of the six cosine terms, in kcal/mol.
    pub v: [f64; 6],
}

/// Pre-computed UFF inversion contribution (out-of-plane)
/// Matches RDKit's InversionContrib exactly:
///   E = K * (C0 + C1*sinY + C2*cos2W)
/// where Y is the Wilson out-of-plane angle.
/// For SP2 C/N/O: C0=1, C1=-1, C2=0 → E = K*(1-sinY)
#[derive(Clone)]
pub struct UFFInversionContrib {
    /// Neighbor atom 1 (idx1 in RDKit)
    pub at1: usize,
    /// Central atom (idx2 in RDKit)
    pub at2: usize,
    /// Neighbor atom 2 (idx3 in RDKit)
    pub at3: usize,
    /// Neighbor atom 3 (idx4 in RDKit)
    pub at4: usize,
    /// Inversion force constant `K` in kcal/mol/rad².
    pub force_constant: f64,
    /// Constant term `C0` (dimensionless).
    pub c0: f64,
    /// `sinY` coefficient `C1` (dimensionless).
    pub c1: f64,
    /// `cos2W` coefficient `C2` (dimensionless).
    pub c2: f64,
}

/// Flat-bottom angle constraint (degrees) for SP (linear) centers.
/// E = k * angleTerm^2, where angleTerm = clamp(angle - [min, max], 0)
#[derive(Clone)]
pub struct AngleConstraint {
    /// Index of the first outer atom of the angle.
    pub i: usize,
    /// Index of the central atom of the angle.
    pub j: usize,
    /// Index of the second outer atom of the angle.
    pub k: usize,
    /// Lower limit of the accepted angle window in **degrees**.
    pub min_deg: f64,
    /// Upper limit of the accepted angle window in **degrees**.
    pub max_deg: f64,
    /// Force constant of the penalty in kcal/mol/rad².
    pub force_k: f64,
}

/// Pre-computed 3D ETKDG force field terms.
/// Distance constraints are split into three groups matching RDKit's construct3DForceField
/// ordering: 1-2 bonds, 1-3 angles, and long-range. This ensures identical
/// floating-point accumulation order in energy/gradient, which is critical
/// for BFGS trajectory reproducibility.
pub struct Etkdg3DFF {
    /// Flat-bottom constraints for 1-2 (bonded) pairs.
    pub dist_12: Vec<DistConstraint>,
    /// Flat-bottom constraints for 1-3 (angle) pairs.
    pub dist_13: Vec<DistConstraint>,
    /// Flat-bottom constraints for long-range (1-4 and beyond) pairs.
    pub dist_long: Vec<DistConstraint>,
    /// Flat-bottom angle-window constraints (used for linear centres).
    pub angle_constraints: Vec<AngleConstraint>,
    /// M6 empirical torsion contributions.
    pub torsion_contribs: Vec<M6TorsionContrib>,
    /// UFF out-of-plane (inversion) contributions.
    pub inversion_contribs: Vec<UFFInversionContrib>,
    /// Weight of the generic out-of-plane volume penalty (kcal/mol/Å⁶).
    pub oop_k: f64,
    /// Weight of the generic (non-M6) torsion term in kcal/mol.
    pub torsion_k_omega: f64,
    /// Scaling factor applied to the distance-constraint force constants
    /// (dimensionless).
    pub bounds_force_scaling: f64,
    /// Whether the M6 empirical torsion terms are used (otherwise the generic
    /// torsion term with [`Self::torsion_k_omega`] is used).
    pub use_m6_torsions: bool,
}
