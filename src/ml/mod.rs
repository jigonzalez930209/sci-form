//! Lightweight ML property proxies.
//!
//! Descriptor-based property prediction using pre-fitted linear models.
//! These provide fast estimates when full quantum-chemical calculations
//! are too expensive.

pub mod descriptors;
pub mod models;

pub use descriptors::{compute_descriptors, MolecularDescriptors};
pub use models::{predict_properties, MlPropertyResult};
