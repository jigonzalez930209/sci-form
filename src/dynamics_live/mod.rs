//! **ALPHA** — Zero-copy dynamic simulation state for real-time molecular dynamics.
//!
//! This module provides:
//! - [`DynamicAtom`] — atom with position, velocity, force vectors
//! - [`LiveMolecularSystem`] — full simulation state with flat buffer for WASM pointer export
//! - [`SimulationBox`] — periodic boundary conditions (minimum-image convention)
//! - [`NeighborList`] — Verlet neighbor list for O(N) force evaluation with PBC
//! - [`BerendsenBarostat`] — Berendsen NPT pressure coupling
//! - [`crate::imd::SteeringForce`] — IMD harmonic steering potential
//!
//! # Feature gate
//!
//! Requires `alpha-dynamics-live`:
//! ```toml
//! sci-form = { version = "0.11", features = ["alpha-dynamics-live"] }
//! ```
//!
//! # Stability tier
//!
//! **Alpha** — early research tier. The validation gate is
//! `tests/alpha_dynamics_validation.rs` (F84): the live forces are compared
//! against RDKit's UFF implementation, and the NVE drift and time-reversibility
//! are checked. Promotion criteria are listed in
//! [`crate::experimental_status`].
//!
//! # Units
//!
//! Positions and lengths in **Å**, velocities in **Å/fs**, forces in
//! **kcal/(mol·Å)**, masses in **amu**, energies in **kcal/mol**, time in
//! **fs**, temperature in **K** and pressure in **bar**.

pub mod barostat;
pub mod langevin;
pub mod neighbor_list;
pub mod pbc;
pub mod state;

/// Steering lives in core [`crate::imd`]; re-exported here for compatibility.
pub use crate::imd::*;
pub use barostat::*;
pub use langevin::*;
pub use neighbor_list::*;
pub use pbc::*;
pub use state::*;
