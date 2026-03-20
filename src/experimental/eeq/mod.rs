//! Dynamic EEQ Force Field — Track E5
//!
//! Geometry-dependent electronegativity equalization charges with
//! Coulomb damping and D3-BJ dispersion correction.

mod charges;
mod energy;

pub use charges::{
    compute_eeq_charges, fractional_coordination, EeqChargeResult, EeqConfig, EeqParams,
};
pub use energy::{compute_eeq_energy, compute_eeq_gradient, EeqEnergyResult};
