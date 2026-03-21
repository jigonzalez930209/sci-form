//! Physical and mathematical constants for quantum chemistry.
//!
//! All values follow CODATA 2018 recommended values.

/// Conversion: 1 Å = 1.8897259886 Bohr.
pub const ANGSTROM_TO_BOHR: f64 = 1.8897259886;

/// Conversion: 1 Bohr = 0.529177 Å.
pub const BOHR_TO_ANGSTROM: f64 = 0.52917721067;

/// 1 Hartree = 27.211386 eV.
pub const HARTREE_TO_EV: f64 = 27.211386245988;

/// 1 eV = 0.0367493 Hartree.
pub const EV_TO_HARTREE: f64 = 1.0 / HARTREE_TO_EV;

/// 1 Hartree = 627.5095 kcal/mol.
pub const HARTREE_TO_KCAL: f64 = 627.5094740631;

/// Default SCF max iterations.
pub const SCF_MAX_ITER: usize = 128;

/// Default SCF energy convergence threshold (Hartree).
pub const SCF_ENERGY_THRESHOLD: f64 = 1e-8;

/// Default SCF density convergence threshold.
pub const SCF_DENSITY_THRESHOLD: f64 = 1e-6;

/// Default DIIS subspace size.
pub const DIIS_SUBSPACE_SIZE: usize = 8;
