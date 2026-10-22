//! PM3 (Parameterized Model 3) semi-empirical method.
//!
//! Implements the NDDO (Neglect of Diatomic Differential Overlap) approximation
//! with Stewart's PM3 parameterization for common organic elements.
//!
//! Two solvers are available:
//!
//! - [`faithful`] — a term-by-term port of MOPAC 23.2.5 (analytic s/p
//!   Slater overlaps, the 22 local two-electron integrals of `reppd`, the
//!   NDDO Fock build of `fock1`/`fock2`, and the MOPAC heat-of-formation
//!   convention). It reproduces MOPAC 23.2.5 energies to better than
//!   0.02 kcal/mol (validation F42).
//! - [`solver`] — the crate's own NDDO implementation.
//!
//! [`solver::solve_pm3`] runs the faithful port for molecules whose elements
//! all have `Z <= 17` and falls back to the crate's own NDDO path for
//! heavier elements or if the faithful path reports an error; the `notes`
//! field of the result says which path produced the numbers.
//!
//! Transition metals are evaluated with the separate PM3(tm) parameterization
//! (Cundari et al., J. Chem. Inf. Comput. Sci. 38 (1998): 941); use
//! [`params::pm3_method`] to query which set an element uses.
//!
//! Energies are reported in eV, except the heat of formation, which is in
//! kcal/mol.
//!
//! References: Stewart, J. J. P. "Optimization of Parameters for Semiempirical Methods I.
//! Method," J. Comput. Chem. 10 (1989): 209–220; Cundari et al. (1998).

pub mod faithful;
#[cfg(feature = "experimental-gpu")]
pub mod gpu;
pub mod gradients;
pub mod params;
mod params_table;
pub mod solver;

pub use gradients::{compute_pm3_gradient, Pm3GradientResult};
pub use params::{get_pm3_params, is_pm3_supported, pm3_method, Pm3Method, Pm3Params};
pub use solver::{solve_pm3, Pm3Result};

/// Temporary probe for the faithful PM3 solver (validation harness).
///
/// # Errors
///
/// Propagates the errors of [`faithful::solve`]: an invalid geometry, an
/// element without PM3 parameters, an element with `Z > 17`, or an odd
/// electron count.
#[doc(hidden)]
pub fn faithful_probe(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<faithful::FaithfulResult, String> {
    faithful::solve(elements, positions)
}
