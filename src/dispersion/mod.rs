//! DFT-D4 Dispersion Correction — Core Module
//!
//! Geometry-dependent London dispersion with Becke-Johnson damping.
//! Provides C6/C8 two-body and optional ATM three-body corrections.
//!
//! ## Relationship to other methods
//!
//! - **HF-3c**: Uses D3-BJ internally; D4 is the next-generation replacement
//! - **PM3**: Can be augmented with D4 for non-covalent interactions
//! - **xTB**: GFN-xTB includes D4 natively; standalone D4 useful for EHT-D4
//! - **UFF/MMFF94**: D4 improves long-range dispersion accuracy
//!
//! ## Reference
//!
//! Caldeweyher, Ehlert, Hansen, Bauer, Spicher, Grimme,
//! J. Chem. Phys. 150, 154122 (2019).

#[allow(clippy::module_inception)]
mod dispersion;
mod params;

pub use dispersion::{compute_d4_energy, compute_d4_gradient, D4Config, D4Result};
pub use params::{
    c8_from_c6, d4_coordination_number, dynamic_c6, get_c6_reference, get_d4_params, D4Params,
};
