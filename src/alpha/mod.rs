//! Alpha-tier compatibility aliases and remaining early-research modules.
//!
//! Most of the original alpha tier has been promoted to core after its
//! cross-library gates passed. The paths below remain available as aliases;
//! the remaining gated research modules are DFT, MLFF, reaction dynamics and
//! the render-bridge domain builders.
//!
//! ## Modules
//!
//! | Module | Purpose | Status |
//! |--------|---------|--------|
//! | CGA | Conformal Geometric Algebra | core `crate::cga`; `alpha::cga` alias |
//! | GSM | Growing String Method | core `crate::gsm`; `alpha::gsm` alias |
//! | SDR | Spectral dimensionality reduction | core `crate::sdr`; `alpha::sdr` alias |
//! | EDL | Electrochemical double layer | core `crate::edl`; `alpha::edl` alias |
//! | PeriodicLinear | Periodic adapters (KPM/RandNLA) | core `crate::periodic_linear`; alias |
//! | Kinetics | HTST and microkinetic networks | core `crate::kinetics`; alias |
//! | RenderBridge | Chart builders | core payloads + gated domain builders |
//! | ReactionDynamics | Full 3D reaction path | gated; self-consistency only |
//!
//! ## Usage
//!
//! Core modules need no feature flag. Compatibility aliases require their
//! historical flag:
//! ```toml
//! [dependencies]
//! sci-form = { version = "1.1", features = ["alpha-cga"] }
//! ```

#[cfg(feature = "alpha-cga")]
pub mod cga;

#[cfg(feature = "alpha-edl")]
pub mod edl;

/// Backwards-compatible path for the promoted core GSM module.
#[cfg(feature = "alpha-gsm")]
pub use crate::gsm;

/// Backwards-compatible path for the promoted core kinetics module.
#[cfg(feature = "alpha-kinetics")]
pub use crate::kinetics;

/// Backwards-compatible path for the promoted core periodic-linear module.
#[cfg(feature = "alpha-periodic-linear")]
pub use crate::periodic_linear;

#[cfg(feature = "alpha-render-bridge")]
pub mod render_bridge;

#[cfg(feature = "alpha-sdr")]
pub mod sdr;

#[cfg(feature = "alpha-reaction-dynamics")]
pub mod reaction_dynamics;
