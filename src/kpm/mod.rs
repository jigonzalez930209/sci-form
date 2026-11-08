//! Kernel Polynomial Method (KPM) — Track E4
//!
//! O(N) computation of DOS, density matrix, and Mulliken charges via
//! Chebyshev polynomial expansion, avoiding diagonalization entirely.
//!
//! The method follows the Kernel Polynomial Method of restraining the
//! spectral function of a Hermitian matrix (Lin, Saad and Yang, *SIAM J.
//! Sci. Comput.* 2011, 33, 549), as used for O(N) electronic-structure
//! scaling.
//!
//! # Stability tier
//!
//! **Core** — promoted from `beta-kpm` on 2026-11-07 after clearing the
//! three technical criteria: sparse CSR support, <0.1 eV spectral accuracy
//! against exact diagonalization, and a 5000-site scale demonstration. The
//! old `sci_form::beta::kpm` path is kept as an alias.
//!
//! Validation gate: `tests/beta_kpm_validation.rs` (F91) compares the
//! Chebyshev moments against exact values, the DOS against
//! `numpy.linalg.eigh`, and the sparse CSR path against the dense path.
//!
//! # Units
//!
//! Energies, spectral bounds and the Fermi level are in **eV**, the
//! electronic temperature in **eV**, the density of states in
//! **states/eV** and the Mulliken charges in **e** (elementary charge).

mod chebyshev;
mod density;
mod sparse;

pub use chebyshev::{estimate_spectral_bounds, jackson_kernel, rescale_matrix, ChebyshevExpansion};
pub use density::{
    compute_kpm_dos, compute_kpm_mulliken, fermi_dirac_coefficients, KpmConfig, KpmDosResult,
    KpmMullikenResult,
};
pub use sparse::{compute_kpm_dos_sparse, SparseMatrix};
