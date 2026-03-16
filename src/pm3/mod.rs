//! PM3 (Parameterized Model 3) semi-empirical method.
//!
//! Implements the NDDO (Neglect of Diatomic Differential Overlap) approximation
//! with Stewart's PM3 parameterization for common organic elements.
//!
//! Reference: Stewart, J. J. P. "Optimization of Parameters for Semiempirical Methods I.
//! Method," J. Comput. Chem. 10 (1989): 209–220.

pub mod params;
pub mod solver;

pub use params::{Pm3Params, get_pm3_params, is_pm3_supported};
pub use solver::{Pm3Result, solve_pm3};
