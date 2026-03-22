//! PM3 (Parameterized Model 3) semi-empirical method.
//!
//! Implements the NDDO (Neglect of Diatomic Differential Overlap) approximation
//! with Stewart's PM3 parameterization for common organic elements.
//!
//! Reference: Stewart, J. J. P. "Optimization of Parameters for Semiempirical Methods I.
//! Method," J. Comput. Chem. 10 (1989): 209–220.

#[cfg(feature = "experimental-gpu")]
pub mod gpu;
pub mod gradients;
pub mod params;
pub mod solver;

pub use gradients::{compute_pm3_gradient, Pm3GradientResult};
pub use params::{get_pm3_params, is_pm3_supported, Pm3Params};
pub use solver::{solve_pm3, Pm3Result};
