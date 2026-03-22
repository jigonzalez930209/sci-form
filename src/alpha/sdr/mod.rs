//! E8: Semidefinite Relaxation Embedding (SDR)
//!
//! Convex optimization approach to conformer embedding that guarantees
//! positive semidefinite Gram matrices, eliminating retry loops from
//! negative eigenvalues.

pub mod embedding;
pub mod projections;

pub use embedding::{extract_coordinates, sdr_embed, warm_start_gram, SdrResult};
pub use projections::{
    alternating_projections, project_distances, project_psd, SdrConfig, SdrConvergence,
};
