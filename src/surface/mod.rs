//! Surface analysis: solvent-accessible surface area and related properties.
//!
//! Implements:
//! - **Shrake-Rupley** SASA algorithm
//!   (A. Shrake & J.A. Rupley, *J. Mol. Biol.* **79**, 351–371, 1973)
//!
//! ## Radii and fallbacks
//!
//! The atomic radii are Bondi values for the common elements and a documented
//! 1.70 Å fallback otherwise; see [`sasa::vdw_radius`] and
//! [`sasa::vdw_radius_source`] (validation F28 compares against FreeSASA with
//! the same radii/probe/point set).
//!
//! ## Connolly (solvent-excluded) surface
//!
//! The Shrake-Rupley implementation returns the solvent-*accessible* surface
//! (SAS, probe centre). A true Connolly solvent-*excluded* surface (SES) needs
//! the reentrant region between the inflated spheres; a mesh-level SES is not
//! implemented. Assessment: the SAS area is what the solvation models (SASA
//! non-polar term, GB/ALPB) consume, and the SES area can be approximated from
//! the SAS by the standard probe-offset relation when needed, so a dedicated
//! Connolly tessellation is deferred until a consumer requires the mesh.

pub mod sasa;

#[cfg(feature = "parallel")]
pub use sasa::compute_sasa_parallel;
pub use sasa::{compute_sasa, SasaResult, VDW_FALLBACK_RADIUS};
