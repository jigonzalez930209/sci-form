//! RandNLA for EHT — Track E2
//!
//! Replaces the O(N³) EHT diagonalization with randomized Nyström
//! approximation for O(N² k) cost, enabling larger molecular systems.
//!
//! # Stability tier
//!
//! **Beta** — tested tier, approaching core readiness. Gated behind the
//! `beta-randnla` feature flag (alias: `experimental-randnla`). The
//! validation gate is `tests/beta_randnla_validation.rs` (F93): the
//! rank-N/4 approximation error is characterised against the exact
//! solution. The "<1 % error at rank = N/4" promotion criterion is **not**
//! yet met; see [`crate::experimental_status`].
//!
//! # Units
//!
//! Hamiltonian and overlap matrices in **eV** (the EHT convention of this
//! crate); the sketch and Nyström factors are dimensionless.

mod nystrom;
mod solver;

pub use nystrom::{GaussianSketch, NystromApprox};
pub use solver::{solve_eht_randnla, RandNlaConfig, RandNlaInfo};
