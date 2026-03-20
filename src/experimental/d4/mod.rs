//! Grimme D4 Dispersion Correction — Track E7
//!
//! Geometry-dependent C6/C8 dispersion with Becke-Johnson damping,
//! supporting EHT-D4 and UFF-D4 augmentation.

mod params;
mod dispersion;

pub use params::{d4_coordination_number, get_c6_reference, D4Params};
pub use dispersion::{compute_d4_energy, compute_d4_gradient, D4Config, D4Result};
