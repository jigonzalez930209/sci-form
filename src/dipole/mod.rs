//! Molecular dipole moments from EHT results.

#[allow(clippy::module_inception)]
pub mod dipole;
pub use dipole::{compute_dipole, compute_dipole_from_eht, DipoleResult};
