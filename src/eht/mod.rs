//! Extended Hückel Theory (EHT) module.
//!
//! Implements a semiempirical electronic-structure calculation pipeline:
//! - Phase B1: EHT parameters, Slater-type orbitals, STO-nG Gaussian expansions
//! - Phase B2: Overlap matrix S and Hamiltonian matrix H
//! - Phase B3: Generalized eigenproblem solver (Löwdin orthogonalization)
//! - Phase B4: 3D volumetric mapping of molecular orbitals
//! - Phase B5: Output structures for rendering (raw volumes + Marching Cubes)

pub mod basis;
pub mod hamiltonian;
pub mod marching_cubes;
pub mod overlap;
pub mod params;
pub mod solver;
pub mod volume;

pub use basis::{AtomicOrbital, GaussianPrimitive, SlaterOrbital};
pub use hamiltonian::build_hamiltonian;
pub use overlap::build_overlap_matrix;
pub use params::{EhtParams, OrbitalDef};
pub use solver::{EhtResult, solve_eht};
pub use volume::{VolumetricGrid, evaluate_orbital_on_grid};
pub use marching_cubes::{IsosurfaceMesh, marching_cubes};
