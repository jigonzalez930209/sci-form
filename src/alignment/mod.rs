//! Molecular alignment and RMSD calculation.

pub mod kabsch;
pub use kabsch::{compute_rmsd, align_coordinates, AlignmentResult};
