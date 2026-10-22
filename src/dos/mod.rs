//! Density of States (DOS) and Projected DOS (PDOS).
//!
//! The total DOS is built by Gaussian-smearing the orbital energies of an
//! electronic-structure calculation, and the atom-projected DOS (PDOS) by
//! weighting each orbital with its Mulliken populations.
//!
//! - [`dos`] — DOS/PDOS from Extended Hückel orbital energies.
//! - [`multi_method`] — the same curves from EHT, PM3, GFN0/GFN1/GFN2-xTB
//!   or HF-3c orbital energies, selected with [`multi_method::DosMethod`].
//!
//! Energies are on the eV scale and the DOS is in states/eV.

#[allow(clippy::module_inception)]
pub mod dos;
pub mod multi_method;
pub use dos::{compute_dos, compute_pdos, DosResult};
#[cfg(feature = "parallel")]
pub use dos::{compute_dos_parallel, compute_pdos_parallel};
pub use multi_method::{compute_dos_multimethod, DosMethod, MultiMethodDosResult};
