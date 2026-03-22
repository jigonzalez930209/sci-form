//! GFN0-xTB-inspired tight-binding method.
//!
//! A simplified self-consistent tight-binding approach inspired by
//! Grimme's GFN family, covering H–Rn with valence-shell parameters.
//!
//! Reference: Grimme, S.; Bannwarth, C.; Shushkov, P. "A Robust and Accurate
//! Tight-Binding Quantum Chemical Method for Structures, Vibrational
//! Frequencies, and Noncovalent Interactions of Large Molecular Systems
//! Parametrized for All spd-Block Elements (Z = 1–86)." *JCTC* 13 (2017): 1989.

pub mod gfn1;
pub mod gfn2;
#[cfg(feature = "experimental-gpu")]
pub mod gpu;
pub mod gradients;
pub mod params;
pub mod solver;

pub use gfn1::{solve_gfn1, Gfn1Result};
pub use gfn2::{solve_gfn2, Gfn2Result};
pub use gradients::{compute_xtb_gradient, XtbGradientResult};
pub use params::{get_xtb_params, is_xtb_supported, XtbParams};
pub use solver::{solve_xtb, XtbResult};
