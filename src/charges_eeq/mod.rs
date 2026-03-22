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

mod charges;
mod energy;

pub use charges::{
    compute_eeq_charges, fractional_coordination, get_eeq_params, EeqChargeResult, EeqConfig,
    EeqParams,
};
pub use energy::{compute_eeq_energy, compute_eeq_gradient, EeqEnergyResult};
