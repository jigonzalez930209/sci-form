//! PM3 (Parameterized Model 3) semi-empirical method.
//!
//! Implements the NDDO (Neglect of Diatomic Differential Overlap) approximation
//! with Stewart's PM3 parameterization for common organic elements.
//!
//! Transition metals are evaluated with the separate PM3(tm) parameterization
//! (Cundari et al., J. Chem. Inf. Comput. Sci. 38 (1998): 941); use
//! [`params::pm3_method`] to query which set an element uses. Both
//! parameterizations share the NDDO core, whose MOPAC parity is a documented
//! gap (validation F42).
//!
//! References: Stewart, J. J. P. "Optimization of Parameters for Semiempirical Methods I.
//! Method," J. Comput. Chem. 10 (1989): 209–220; Cundari et al. (1998).

#[cfg(feature = "experimental-gpu")]
pub mod gpu;
pub mod gradients;
pub mod params;
pub mod solver;

pub use gradients::{compute_pm3_gradient, Pm3GradientResult};
pub use params::{get_pm3_params, is_pm3_supported, pm3_method, Pm3Method, Pm3Params};
pub use solver::{solve_pm3, Pm3Result};
