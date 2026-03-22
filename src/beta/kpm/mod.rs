//! Kernel Polynomial Method (KPM) — Track E4
//!
//! O(N) computation of DOS, density matrix, and Mulliken charges via
//! Chebyshev polynomial expansion, avoiding diagonalization entirely.

mod chebyshev;
mod density;

pub use chebyshev::{estimate_spectral_bounds, jackson_kernel, rescale_matrix, ChebyshevExpansion};
pub use density::{
    compute_kpm_dos, compute_kpm_mulliken, fermi_dirac_coefficients, KpmConfig, KpmDosResult,
    KpmMullikenResult,
};
