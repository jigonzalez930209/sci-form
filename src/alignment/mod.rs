//! Molecular alignment and RMSD calculation.
//!
//! Superposition of two structures onto each other by the optimal proper
//! rotation and translation, following the Kabsch algorithm of
//! Kabsch, *Acta Cryst.* **1976**, A32, 922, with a quaternion-based
//! alternative after Coutsias, Dill & Church, *Protein Sci.* **2004**, 13, 275.
//!
//! Coordinates are flat `[x0, y0, z0, x1, ...]` arrays in Å, and the RMSD
//! that comes back is in Å.

pub mod kabsch;
pub use kabsch::{
    align_coordinates, align_quaternion, compute_rmsd, try_align_coordinates, try_align_quaternion,
    try_compute_rmsd, AlignmentResult,
};
