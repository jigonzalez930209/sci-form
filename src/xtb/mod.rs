//! GFN-xTB tight-binding methods (GFN0, GFN1 and GFN2).
//!
//! Self-consistent-charge tight-binding methods of Grimme's GFN family,
//! covering H–Rn with valence-shell parameters:
//!
//! - [`gfn0`] — GFN0-xTB: the simplest member of the family, used for
//!   large systems and for qualitative orbital pictures.
//! - [`gfn1`] — GFN1-xTB: adds shell-resolved electrostatics and a
//!   classical halogen-bond correction.
//! - [`gfn2`] — GFN2-xTB: adds multipole electrostatics (dipole and
//!   quadrupole shells) and the charge-dependent D4 dispersion correction.
//!
//! Energies are reported in eV throughout this module.
//!
//! References:
//! - Grimme, S.; Bannwarth, C.; Shushkov, P. "A Robust and Accurate
//!   Tight-Binding Quantum Chemical Method for Structures, Vibrational
//!   Frequencies, and Noncovalent Interactions of Large Molecular Systems
//!   Parametrized for All spd-Block Elements (Z = 1–86)." *JCTC* 13 (2017): 1989.
//!   (GFN0-xTB.)
//! - Grimme, S.; Bannwarth, C.; Shushkov, P. "GFN1-xTB: A Self-Consistent
//!   Tight-Binding Method Using the Generalised Koopmans Approximations for
//!   Large-Scale Molecular Systems." *JCTC* 15 (2019): 1652. (GFN1-xTB.)
//! - Bannwarth, C.; Ehlert, S.; Grimme, S. "GFN2-xTB: A Robust and Accurate
//!   Self-Consistent Tight-Binding Quantum Chemical Method for Large
//!   Molecular Systems." *JCTC* 15 (2019): 1652. (GFN2-xTB.)

pub mod broyden;
pub mod d4;
pub mod d4_data;
pub mod gfn0;
pub mod gfn0_params;
pub mod gfn1;
pub mod gfn1_params;
pub mod gfn2;
pub mod gfn2_params;
#[cfg(feature = "experimental-gpu")]
pub mod gpu;
pub mod gradients;
pub mod params;
pub mod solver;
pub mod spin_params;
pub mod sto_overlap;
pub mod uks;

pub use gfn0::{solve_gfn0, Gfn0Result};
pub use gfn1::{solve_gfn1, Gfn1Result};
pub use gfn2::{solve_gfn2, Gfn2Result};
pub use gradients::{compute_xtb_gradient, XtbGradientResult};
pub use params::{get_xtb_params, is_xtb_supported, XtbParams};
pub use solver::{solve_xtb, XtbResult};
pub use uks::{solve_gfn1_unrestricted, solve_gfn2_unrestricted, Gfn1UksResult, UksResult};
