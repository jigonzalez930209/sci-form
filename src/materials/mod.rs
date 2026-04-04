//! Crystallographic materials: unit cells, space groups, framework assembly,
//! and periodic geometry optimization.
//!
//! - [`cell`] — Unit cell construction from lattice parameters (a, b, c, α, β, γ).
//! - [`space_groups`] — All 230 ITC space groups with symmetry operations.
//! - [`sbu`] — Secondary Building Unit definitions for framework assembly.
//! - [`assembly`] — MOF/COF framework assembly from topology + SBU + cell.
//! - [`geometry_opt`] — BFGS and steepest-descent optimization with PBC.

pub mod assembly;
pub mod cell;
pub mod cif;
pub mod geometry_opt;
pub mod sbu;
pub mod space_groups;

pub use assembly::*;
pub use cell::*;
pub use cif::*;
pub use sbu::*;
