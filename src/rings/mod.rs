//! Smallest Set of Smallest Rings (SSSR) perception and ECFP fingerprints.
//!
//! - SSSR: Identifies the fundamental cycle basis of the molecular graph
//! - Aromaticity: Hückel perception for Kekulé input (RDKit parity)
//! - ECFP: Extended-Connectivity Fingerprints (Morgan algorithm)

pub mod aromaticity;
pub mod ecfp;
pub mod sssr;

pub use aromaticity::{perceive_aromaticity, perceived_aromaticity};
pub use ecfp::{compute_ecfp, compute_tanimoto, ECFPFingerprint};
pub use sssr::{compute_sssr, RingInfo};
