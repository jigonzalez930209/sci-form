//! E8: Semidefinite Relaxation Embedding (SDR) — core since 2026-10-14 (F83).
//!
//! Convex optimization approach to conformer embedding that guarantees
//! positive semidefinite Gram matrices, eliminating retry loops from
//! negative eigenvalues. Validated vs numpy + CVXPY/SCS; `alpha-sdr`
//! survives as a no-op feature alias.

pub mod embedding;
pub mod projections;

pub use embedding::{extract_coordinates, sdr_embed, warm_start_gram, SdrResult};
pub use projections::{
    alternating_projections, project_distances, project_psd, svt_step, SdrConfig, SdrConvergence,
};
