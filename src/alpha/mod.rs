//! Alpha-stage experimental modules — early research, proof-of-concept.
//!
//! These modules are at the earliest stage of development. They demonstrate
//! algorithmic feasibility but lack:
//! - Production-quality validation
//! - Integration with core APIs
//! - Performance benchmarks against established methods
//! - Comprehensive edge-case handling
//!
//! ## Modules
//!
//! | Module | Purpose | Promotion criteria |
//! |--------|---------|-------------------|
//! | CGA | Conformal Geometric Algebra for molecular geometry | Benchmark vs ETKDG, symmetry-constrained assembly |
//! | GSM | Growing String Method for reaction pathways | Couple with PM3/xTB gradients, reproduce known barriers |
//! | SDR | Spectral Dimensionality Reduction for chemical space | Validate SAR preservation, adaptive bandwidth |
//! | EDL | Electrochemical double-layer contracts and compact-layer CPU reference | Couple CPM/EEQ/ALPB into validated interface scans |
//! | PeriodicLinear | Periodic linear-scaling contracts, k-meshes, and Bloch helpers | k-mesh DOS and low-rank refinement on periodic operators |
//! | Kinetics | HTST and microkinetic contracts with CPU rate references | Complete GSM -> MBH -> HTST and network validation |
//! | RenderBridge | Browser- and GPU-ready rendering payload contracts | Typed-array/columnar export with round-trip benchmarks |
//!
//! ## Usage
//!
//! Each module requires its own feature flag:
//! ```toml
//! [dependencies]
//! sci-form = { version = "0.9", features = ["alpha-cga"] }
//! ```

#[cfg(feature = "alpha-cga")]
pub mod cga;

#[cfg(feature = "alpha-edl")]
pub mod edl;

#[cfg(feature = "alpha-gsm")]
pub mod gsm;

#[cfg(feature = "alpha-kinetics")]
pub mod kinetics;

#[cfg(feature = "alpha-periodic-linear")]
pub mod periodic_linear;

#[cfg(feature = "alpha-render-bridge")]
pub mod render_bridge;

#[cfg(feature = "alpha-sdr")]
pub mod sdr;
