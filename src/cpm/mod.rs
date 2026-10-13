//! E10: Constant Potential Method (CPM)
//!
//! Couples molecular charge equilibration to a virtual electrode potential μ,
//! enabling electrochemical property predictions: quantum capacitance,
//! charge-potential curves, and conformer ranking under applied voltage.

pub mod grand_potential;
pub mod surface;

pub use grand_potential::{compute_cpm_charges, CpmConfig, CpmResult};
pub use surface::{compute_cpm_surface, CpmSurface};
