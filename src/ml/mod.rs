//! Lightweight ML property proxies.
//!
//! Descriptor-based property prediction using pre-fitted linear models.
//! These provide fast estimates when full quantum-chemical calculations
//! are too expensive.

pub mod advanced_models;
pub mod descriptors;
pub mod ensemble;
pub mod getaway;
pub mod models;
pub mod pharmacophore;
pub mod rdf_descriptors;
pub mod whim;

#[cfg(feature = "alpha-mlff")]
pub mod inference;
#[cfg(feature = "alpha-mlff")]
pub mod mlff;
#[cfg(feature = "alpha-mlff")]
pub mod symmetry_functions;

pub use descriptors::{
    compute_3d_descriptors, compute_descriptors, Descriptors3D, MolecularDescriptors,
};
pub use ensemble::{compute_tpsa, predict_ensemble, EnsembleResult, VeberResult};
pub use models::{predict_properties, MlPropertyResult, PredictionUncertainty};
pub use pharmacophore::{
    compute_pharmacophore_fingerprint, detect_features, pharmacophore_tanimoto, PharmFeature,
    PharmFeatureType, PharmacophoreFingerprint,
};
