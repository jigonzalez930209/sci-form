//! Density of States (DOS) and Projected DOS (PDOS).

#[allow(clippy::module_inception)]
pub mod dos;
pub use dos::{compute_dos, compute_pdos, DosResult};
