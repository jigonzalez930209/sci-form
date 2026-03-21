//! RandNLA for EHT — Track E2
//!
//! Replaces the O(N³) EHT diagonalization with randomized Nyström
//! approximation for O(N² k) cost, enabling larger molecular systems.

mod nystrom;
mod solver;

pub use nystrom::{GaussianSketch, NystromApprox};
pub use solver::{solve_eht_randnla, RandNlaConfig, RandNlaInfo};
