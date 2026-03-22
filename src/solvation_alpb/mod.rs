//! ALPB (Analytical Linearized Poisson-Boltzmann) Solvation — Core Module
//!
//! Advanced implicit solvation combining Generalized Born electrostatics
//! with the ALPB correction and non-polar SASA cavitation.
//!
//! ## Comparison with existing solvation methods
//!
//! | Method | Electrostatics | Non-polar | Use case |
//! |--------|---------------|-----------|----------|
//! | `nonpolar_solvation` | None | SASA·γ | Hydrophobic effects only |
//! | `gb_solvation` | Still-HCT GB | SASA·γ | General polar solvation |
//! | **`alpb_solvation`** | ALPB-corrected GB | SASA·γ | Most accurate; xTB-compatible |
//!
//! ALPB improves on standard GB by applying the Klamt correction factor
//! that accounts for charge penetration effects in the solute cavity.
//!
//! ## Reference
//!
//! Ehlert, Stahn, Spicher, Grimme, J. Chem. Theory Comput. 17, 4250 (2021).

mod born;
mod solvation;

pub(crate) use born::intrinsic_radius;
pub use born::{compute_born_radii, gb_kernel, AlpbBornRadii};
pub use solvation::{compute_alpb_solvation, AlpbConfig, AlpbResult};
