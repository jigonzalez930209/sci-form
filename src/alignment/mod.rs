//! Molecular alignment and RMSD calculation.

pub mod kabsch;
pub use kabsch::{
    align_coordinates, align_quaternion, compute_rmsd, try_align_coordinates, try_align_quaternion,
    try_compute_rmsd, AlignmentResult,
};
