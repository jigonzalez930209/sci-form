//! Beta-tier compatibility aliases and the remaining gated modules.
//!
//! KPM, CPM, RandNLA and Riemannian are core; `beta::*` paths and the
//! `beta-*` feature flags remain as compatibility aliases. Multi-Basin
//! Hopping (MBH) is the only module still gated as beta research.
//!
//! ## Modules
//!
//! | Module | Purpose | Status |
//! |--------|---------|--------|
//! | KPM | Kernel Polynomial Method | core `crate::kpm`; alias |
//! | CPM | Constant Potential Method | core `crate::cpm`; alias |
//! | RandNLA | Randomized eigensolvers | core `crate::rand_nla`; alias |
//! | Riemannian | PSD manifold geometry | core `crate::riemannian`; alias |
//! | MBH | Multi-Basin Hopping | gated (`beta-mbh`) |
//!
//! ## Usage
//!
//! Core modules need no feature flag; MBH keeps `beta-mbh`:
//! ```toml
//! [dependencies]
//! sci-form = { version = "2.0", features = ["beta-mbh"] }
//! ```

/// Backwards-compatible path for the promoted core KPM module.
pub use crate::kpm;

#[cfg(feature = "beta-mbh")]
pub mod mbh;

#[cfg(feature = "beta-cpm")]
pub mod cpm;

/// Backwards-compatible path for the promoted core RandNLA module.
pub use crate::rand_nla;

#[cfg(feature = "beta-riemannian")]
pub mod riemannian;
