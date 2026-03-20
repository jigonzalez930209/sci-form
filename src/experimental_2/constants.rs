//! Physical and mathematical constants used across the experimental engine.
//!
//! All values follow CODATA 2018 recommended values where applicable.

/// Conversion: 1 Å = 1.8897259886 Bohr
pub const ANGSTROM_TO_BOHR: f64 = 1.8897259886;

/// Conversion: 1 Bohr = 0.529177 Å
pub const BOHR_TO_ANGSTROM: f64 = 0.52917721067;

/// 1 Hartree = 27.211386 eV
pub const HARTREE_TO_EV: f64 = 27.211386245988;

/// 1 eV = 0.0367493 Hartree
pub const EV_TO_HARTREE: f64 = 1.0 / HARTREE_TO_EV;

/// 1 Hartree = 627.5095 kcal/mol
pub const HARTREE_TO_KCAL: f64 = 627.5094740631;

/// 1 eV = 23.0605 kcal/mol
pub const EV_TO_KCAL: f64 = 23.060547830619;

/// 1 eV = 8065.54 cm⁻¹
pub const EV_TO_CM1: f64 = 8065.544005;

/// Pi
pub const PI: f64 = std::f64::consts::PI;

/// 2π
pub const TWO_PI: f64 = 2.0 * PI;

/// Boltzmann constant in eV/K
pub const KB_EV: f64 = 8.617333262e-5;

/// Boltzmann constant in Hartree/K
pub const KB_HARTREE: f64 = 3.166811563e-6;

/// Planck constant in eV·s
pub const HBAR_EV_S: f64 = 6.582119569e-16;

/// Speed of light in cm/s
pub const C_CM_S: f64 = 2.99792458e10;

/// Atomic mass unit in kg
pub const AMU_KG: f64 = 1.66053906660e-27;

/// Electron mass in AMU
pub const ELECTRON_MASS_AMU: f64 = 5.4857990943e-4;

/// Bohr magneton in eV/T
pub const BOHR_MAGNETON: f64 = 5.7883818012e-5;

/// Nuclear magneton in eV/T
pub const NUCLEAR_MAGNETON: f64 = 3.15245125844e-8;

/// Fine structure constant
pub const ALPHA_FINE: f64 = 7.2973525693e-3;

/// Maximum angular momentum quantum number supported
pub const L_MAX: usize = 2; // s, p, d

/// Threshold for negligible overlap integrals
pub const INTEGRAL_THRESHOLD: f64 = 1.0e-12;

/// SCF default energy convergence (Hartree)
pub const SCF_ENERGY_THRESHOLD: f64 = 1.0e-8;

/// SCF default density convergence
pub const SCF_DENSITY_THRESHOLD: f64 = 1.0e-6;

/// SCF maximum iterations
pub const SCF_MAX_ITER: usize = 200;

/// DIIS default subspace size
pub const DIIS_SUBSPACE_SIZE: usize = 8;

/// Finite difference step for numerical Hessian (Bohr)
pub const HESSIAN_STEP_BOHR: f64 = 0.005;

/// Wolfsberg-Helmholtz constant K
pub const K_WOLFSBERG_HELMHOLTZ: f64 = 1.75;
