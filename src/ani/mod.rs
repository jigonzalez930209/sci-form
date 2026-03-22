//! ANI (Accurate NeurAl networK engINe) machine-learning potentials.
//!
//! Implements the ANI potential family for fast energy and force prediction
//! with near-DFT accuracy using neural-network inference on Behler-Parrinello
//! atomic environment vectors (AEVs).
//!
//! Architecture: positions → neighbor list → AEVs → per-element NN → total energy.
//! Forces are computed via analytical backpropagation through the AEV pipeline.

pub mod aev;
pub mod aev_params;
pub mod ani_tm;
pub mod api;
pub mod cutoff;
pub mod gradients;
pub mod neighbor;
pub mod nn;
pub mod weights;

pub use ani_tm::{compute_aevs_tm, is_ani_tm_supported, AniTmResult};
pub use api::{compute_ani, AniConfig, AniResult};
