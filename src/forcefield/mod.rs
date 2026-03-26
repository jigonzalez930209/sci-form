//! Molecular force fields: UFF, MMFF94, ETKDG refinement, and minimization.
//!
//! # Submodules
//!
//! - [`atom_typer`] — Element + topology → UFF atom type assignment.
//! - [`params`] — UFF parameter tables (bond radii, angles, torsions).
//! - [`uff`] / [`mmff94`] — Energy term implementations per force field.
//! - [`energy`] — Unified total-energy evaluation.
//! - [`gradients`] — Analytical gradient computation.
//! - [`minimizer`] — L-BFGS energy minimizer.
//! - [`bounds_ff`] — Bounds-matrix force field for distance geometry.
//! - [`etkdg_3d`] — ETKDG 3D refinement with torsion terms.
//! - [`dg_terms`] — Distance-geometry specific energy contributions.
//! - [`torsion_scan`] — Systematic torsion angle scanning.

pub mod atom_typer;
pub mod bounds_ff;
pub mod builder;
pub mod dg_terms;
pub mod energy;
pub mod etkdg_3d;
pub mod etkdg_lite;
pub mod gradients;
pub mod minimizer;
pub mod mmff94;
pub mod params;
pub mod torsion_scan;
pub mod traits;
pub mod uff;

pub use bounds_ff::*;
pub use dg_terms::*;
pub use energy::*;
pub use etkdg_3d::*;
pub use gradients::*;
pub use minimizer::*;
pub use torsion_scan::*;
pub use traits::*;
