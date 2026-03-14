//! Molecular alignment and RMSD calculation.

pub mod kabsch;
pub use kabsch::{align_coordinates, compute_rmsd, AlignmentResult};
