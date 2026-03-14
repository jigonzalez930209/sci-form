//! ETKDG 3D force field matching RDKit's construct3DForceField.
//! Uses flat-bottom distance constraints (zero gradient within tolerance)
//! instead of harmonic bond stretch + angle bend.

pub mod builder;
pub mod energy;
pub mod gradient;
pub mod optimizer;

pub use builder::*;
pub use energy::*;
pub use gradient::*;
pub use optimizer::*;

const KNOWN_DIST_TOL: f64 = 0.01;
const KNOWN_DIST_K: f64 = 100.0;

/// A pre-computed flat-bottom distance constraint
#[derive(Clone)]
pub struct DistConstraint {
    pub i: usize,
    pub j: usize,
    pub min_len: f64,
    pub max_len: f64,
    pub k: f64,
}

/// Pre-computed M6 torsion contribution (Fourier series)
#[derive(Clone)]
pub struct M6TorsionContrib {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub l: usize,
    pub signs: [f64; 6],
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
    pub force_constant: f64,
    pub c0: f64,
    pub c1: f64,
    pub c2: f64,
}

/// Flat-bottom angle constraint (degrees) for SP (linear) centers.
/// E = k * angleTerm^2, where angleTerm = clamp(angle - [min,max], 0)
#[derive(Clone)]
pub struct AngleConstraint {
    pub i: usize, // outer atom 1
    pub j: usize, // central atom
    pub k: usize, // outer atom 2
    pub min_deg: f64,
    pub max_deg: f64,
    pub force_k: f64,
}

/// Pre-computed 3D ETKDG force field terms.
/// Distance constraints are split into three groups matching RDKit's construct3DForceField
/// ordering: 1-2 bonds, 1-3 angles, and long-range. This ensures identical
/// floating-point accumulation order in energy/gradient, which is critical
/// for BFGS trajectory reproducibility.
pub struct Etkdg3DFF {
    pub dist_12: Vec<DistConstraint>,
    pub dist_13: Vec<DistConstraint>,
    pub dist_long: Vec<DistConstraint>,
    pub angle_constraints: Vec<AngleConstraint>,
    pub torsion_contribs: Vec<M6TorsionContrib>,
    pub inversion_contribs: Vec<UFFInversionContrib>,
    pub oop_k: f64,
    pub torsion_k_omega: f64,
    pub bounds_force_scaling: f64,
    pub use_m6_torsions: bool,
}
