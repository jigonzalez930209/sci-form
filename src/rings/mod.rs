//! Smallest Set of Smallest Rings (SSSR) perception and ECFP fingerprints.
//!
//! - SSSR: Identifies the fundamental cycle basis of the molecular graph
//! - ECFP: Extended-Connectivity Fingerprints (Morgan algorithm)

pub mod ecfp;
pub mod sssr;

pub use ecfp::{compute_ecfp, compute_tanimoto, ECFPFingerprint};
pub use sssr::{compute_sssr, RingInfo};
