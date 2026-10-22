//! EEQ (Electronegativity Equilibration) Charges — Core Module
//!
//! Provides geometry-dependent partial atomic charges as an alternative to
//! Gasteiger-Marsili. EEQ accounts for 3D structure through coordination
//! numbers and damped Coulomb interactions.
//!
//! ## When to use EEQ vs Gasteiger
//!
//! - **Gasteiger**: Fast, topology-only (no 3D needed). Good for quick screening.
//! - **EEQ**: Geometry-dependent, more accurate for charged/polar systems.
//!   Requires 3D coordinates. Better for solvation, electrostatics.
//!
//! ## Reference
//!
//! Based on the charge equilibration scheme used in GFN-FF
//! (Grimme group, Spicher & Grimme, Angew. Chem. Int. Ed. 2020).
//!
//! # Units
//!
//! Positions and the EEQ charge radii in **Å**, the electronegativity χ and
//! the chemical hardness η in **eV**, and the partial charges in **e**
//! (elementary charge).
//!
//! # Stability tier
//!
//! Core (promoted to core: topology-free partial charges, O(N²); the EEQ2019
//! port was validated against dftd4, F46).

mod charges;
mod eeq2019;
mod eeq2019_data;
mod eeqbc;
mod eeqbc_data;
mod energy;

pub use charges::{
    compute_eeq_charges, fractional_coordination, get_eeq_params, EeqChargeResult, EeqConfig,
    EeqParams,
};
pub use eeq2019::compute_eeq2019_charges;
pub use eeqbc::{
    compute_eeqbc_charges, eeqbc_coordination_numbers, eeqbc_local_charges, EEQBC_MAX_Z,
};
pub use energy::{compute_eeq_energy, compute_eeq_gradient, EeqEnergyResult};
