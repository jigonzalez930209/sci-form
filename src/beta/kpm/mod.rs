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
//! **Beta** — tested tier, approaching core readiness. Gated behind the
//! `beta-kpm` feature flag. The validation gate is
//! `tests/beta_kpm_validation.rs` (F91): the Chebyshev moments are compared
//! against exact values and the DOS against `numpy.linalg.eigh`. Promotion
//! criteria are listed in [`crate::experimental_status`].
//!
//! # Units
//!
//! Energies, spectral bounds and the Fermi level are in **eV**, the
//! electronic temperature in **eV**, the density of states in
//! **states/eV** and the Mulliken charges in **e** (elementary charge).

mod chebyshev;
mod density;

pub use chebyshev::{estimate_spectral_bounds, jackson_kernel, rescale_matrix, ChebyshevExpansion};
pub use density::{
    compute_kpm_dos, compute_kpm_mulliken, fermi_dirac_coefficients, KpmConfig, KpmDosResult,
    KpmMullikenResult,
};
