//! Partial-charge solvers for molecular systems.
//!
//! Currently implements:
//! - **Gasteiger-Marsili** iterative electronegativity equalization
//!   (J. Gasteiger & M. Marsili, *Tetrahedron* **36**, 3219–3228, 1980)

pub mod gasteiger;

pub use gasteiger::{
    gasteiger_marsili_charges, gasteiger_marsili_charges_configured, GasteigerConfig,
    GasteigerParams,
};
