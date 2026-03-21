//! E8: Semidefinite Relaxation Embedding (SDR)
//!
//! Convex optimization approach to conformer embedding that guarantees
//! positive semidefinite Gram matrices, eliminating retry loops from
//! negative eigenvalues.

pub mod projections;
pub mod embedding;

pub use projections::{
    project_psd, project_distances, alternating_projections,
    SdrConfig, SdrConvergence,
};
pub use embedding::{
    sdr_embed, warm_start_gram, extract_coordinates,
    SdrResult,
};
