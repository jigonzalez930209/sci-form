//! Density of States (DOS) and Projected DOS (PDOS).

#[allow(clippy::module_inception)]
pub mod dos;
pub use dos::{compute_dos, compute_pdos, DosResult};
#[cfg(feature = "parallel")]
pub use dos::{compute_dos_parallel, compute_pdos_parallel};
