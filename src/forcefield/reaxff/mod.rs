//! ReaxFF reactive force field (core since 2026-10-13, F78 vs LAMMPS `reaxff`).
//!
//! Implements continuous bond-order based reactive force fields that allow
//! bond breaking and formation during molecular dynamics.
//!
//! # Modules
//!
//! - [`bond_order`] — Continuous bond-order calculation from distances
//! - [`taper`] — Hermite tapering polynomial for smooth cutoff
//! - [`energy`] — Bonded energy terms as functions of continuous BO
//! - [`nonbonded`] — van der Waals + shielded Coulomb
//! - [`eem`] — Electronegativity Equalization Method for charges
//! - [`params`] — ReaxFF parameter parsing (`ffield.reax.cho` pair table)
//! - [`gradients`] — Numerical (central-difference) gradients
//!
//! `alpha-reaxff` survives as a no-op feature alias.

pub mod bond_order;
pub mod eem;
pub mod energy;
pub mod gradients;
pub mod nonbonded;
pub mod params;
pub mod taper;

pub use bond_order::compute_bond_orders;
pub use energy::ReaxffEnergy;
pub use gradients::compute_reaxff_gradient;
pub use params::ReaxffParams;
