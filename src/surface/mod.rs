//! Surface analysis: solvent-accessible surface area and related properties.
//!
//! Implements:
//! - **Shrake-Rupley** SASA algorithm
//!   (A. Shrake & J.A. Rupley, *J. Mol. Biol.* **79**, 351–371, 1973)

pub mod sasa;

#[cfg(feature = "parallel")]
pub use sasa::compute_sasa_parallel;
pub use sasa::{compute_sasa, SasaResult};
