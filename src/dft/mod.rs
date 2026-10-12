//! **ALPHA** — Kohn-Sham Density Functional Theory.
//!
//! Implements:
//! - Molecular integration grid (Euler-Maclaurin radial + Lebedev angular + Becke partitioning)
//! - LDA functional: SVWN (Slater exchange + VWN-5 correlation)
//! - GGA functional: PBE exchange-correlation
//! - V_XC matrix assembly and KS-DFT Fock builder
//!
//! # Feature gate
//!
//! Requires `alpha-dft`:
//! ```toml
//! sci-form = { version = "0.11", features = ["alpha-dft"] }
//! ```

pub mod becke;
pub mod ecp;
pub mod functionals;
pub mod gradient;
pub mod grid;
pub mod ks_fock;
pub mod lebedev;
mod valence_basis;
pub mod vxc_matrix;

pub use gradient::ks_nuclear_gradient;
pub use grid::MolecularGrid;
pub use ks_fock::{solve_ks_dft, DftConfig, DftMethod, DftResult};
pub use valence_basis::DftBasisKind;
