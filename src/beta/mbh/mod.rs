//! E9: Mobile Block Hessian (MBH)
//!
//! Reduces vibrational analysis cost by treating rigid groups (aromatic rings,
//! clusters) as rigid bodies with 6 DOF each, enabling real-time IR spectra
//! for large molecules.
//!
//! The method follows the Mobile Block Hessian of Ghysels, Vanduyfhuys and
//! Van Speybroeck (J. Chem. Phys. 2007, 126, 241104), which projects the
//! full Hessian onto a reduced basis of rigid blocks and flexible atoms.
//!
//! # Stability tier
//!
//! **Beta** — tested tier, approaching core readiness. Gated behind the
//! `beta-mbh` feature flag (alias: `experimental-mbh`). The validation gate
//! is `tests/beta_mbh_validation.rs` (F92): the block decomposition is
//! compared against the closed-form stretch of diatomics and the benzene
//! block structure. Promotion criteria are listed in
//! [`crate::experimental_status`].
//!
//! # Units
//!
//! Coordinates in **Å**, the finite-difference Hessian step in **Å**,
//! force constants in **kcal/(mol·Å²)**, frequencies in **cm⁻¹** and
//! energies in **kcal/mol**.

pub mod blocks;
pub mod hessian;

pub use blocks::{build_projection_matrix, detect_rigid_blocks, BlockDecomposition, RigidBlock};
pub use hessian::{compute_mbh_frequencies, MbhConfig, MbhResult};
