//! Beta-stage experimental modules — tested, approaching core readiness.
//!
//! These modules have passed initial validation and demonstrate clear
//! scientific value, but need additional work before core promotion:
//! - More extensive numerical validation
//! - Integration with existing core APIs
//! - Performance benchmarks and optimization
//! - Edge-case coverage
//!
//! ## Modules
//!
//! | Module | Purpose | Promotion criteria |
//! |--------|---------|-------------------|
//! | KPM | Kernel Polynomial Method (O(N) DOS) | Wall-clock advantage >2000 atoms, 0.1 eV accuracy |
//! | MBH | Multi-Basin Hopping (global optimization) | Integrate with UFF/MMFF94, reproduce known minima |
//! | CPM | Continuous Phase Methods (crystal surfaces) | Validate vs DFT reference, integrate with materials |
//! | RandNLA | Randomized NLA (fast eigensolvers) | Wire into EHT/SCF, <1% error at rank=N/4 |
//! | Riemannian | Riemannian geometry optimization | Benchmark vs Cartesian L-BFGS, manifold constraints |
//!
//! ## Usage
//!
//! Each module requires its own feature flag:
//! ```toml
//! [dependencies]
//! sci-form = { version = "0.9", features = ["beta-kpm"] }
//! ```
//!
//! ## GPU acceleration candidates
//!
//! Several beta modules are candidates for GPU acceleration:
//! - **KPM**: Chebyshev matrix-vector products → compute shader
//! - **RandNLA**: Random projection and Nyström → parallel GPU
//! - **CPM**: Phase-field PDE solver → grid-parallel compute

#[cfg(feature = "beta-kpm")]
pub mod kpm;

#[cfg(feature = "beta-mbh")]
pub mod mbh;

#[cfg(feature = "beta-cpm")]
pub mod cpm;

#[cfg(feature = "beta-randnla")]
pub mod rand_nla;

#[cfg(feature = "beta-riemannian")]
pub mod riemannian;
