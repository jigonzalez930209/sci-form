//! E9: Mobile Block Hessian (MBH)
//!
//! Reduces vibrational analysis cost by treating rigid groups (aromatic rings,
//! clusters) as rigid bodies with 6 DOF each, enabling real-time IR spectra
//! for large molecules.

pub mod blocks;
pub mod hessian;

pub use blocks::{
    detect_rigid_blocks, build_projection_matrix,
    RigidBlock, BlockDecomposition,
};
pub use hessian::{
    compute_mbh_frequencies, MbhConfig, MbhResult,
};
