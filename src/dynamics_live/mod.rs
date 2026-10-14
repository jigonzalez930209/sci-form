//! **ALPHA** — Zero-copy dynamic simulation state for real-time molecular dynamics.
//!
//! This module provides:
//! - [`DynamicAtom`] — atom with position, velocity, force vectors
//! - [`MolecularSystem`] — full simulation state with flat buffer for WASM pointer export
//! - [`SimulationBox`] — periodic boundary conditions (minimum-image convention)
//! - [`NeighborList`] — Verlet neighbor list for O(N) force evaluation with PBC
//! - [`Barostat`] — Berendsen NPT pressure coupling
//! - [`SteeringForce`] — IMD harmonic steering potential
//!
//! # Feature gate
//!
//! Requires `alpha-dynamics-live`:
//! ```toml
//! sci-form = { version = "0.11", features = ["alpha-dynamics-live"] }
//! ```

pub mod barostat;
pub mod langevin;
pub mod neighbor_list;
pub mod pbc;
pub mod state;

pub use barostat::*;
pub use langevin::*;
pub use neighbor_list::*;
pub use pbc::*;
pub use state::*;
//! Steering lives in core [`crate::imd`]; re-exported here for compatibility.
pub use crate::imd::*;
