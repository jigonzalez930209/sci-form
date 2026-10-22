// Several wasm-bindgen entry points mirror the Python/Rust surface and take more
// than seven arguments (coordinates plus per-method parameters).
#![allow(clippy::too_many_arguments)]

//! sci-form WASM bindings.
//!
//! All public functions are split into focused modules.
//!
//! # Conventions
//!
//! * Every function returns a **JSON string**.  On success it is the
//!   documented payload; on failure it is the canonical error payload
//!   `{"error": "<message>"}`, which always parses (blocker C-C1).
//! * The function names match the Rust root functions, so the same capability
//!   has the same name in Rust, Python, WASM and the CLI (blocker C-H1).
//! * Seeds are 64-bit.  `wasm-bindgen` maps a `u64` parameter to a JavaScript
//!   `BigInt`, so every seed-taking entry point has a `*_seed64` variant that
//!   takes the full range; the 32-bit variants are kept for backwards
//!   compatibility and document their narrower range (blocker C-H7).
//! * The `alpha` and `beta` modules are published under the
//!   `sci-form-wasm/alpha` and `sci-form-wasm/beta` sub-paths; the
//!   hand-written declarations live in `crates/wasm/js/{alpha,beta}`.

#[cfg(feature = "parallel")]
pub use wasm_bindgen_rayon::init_thread_pool;

mod helpers;

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
pub mod webgpu;

pub mod dynamics;
pub mod eht;
pub mod electronic;
pub mod embed;
pub mod esp;
pub mod experimental;
pub mod forcefield;
pub mod materials;
pub mod mesh;
pub mod ml;
pub mod properties;
pub mod reaction;
pub mod reactivity;
pub mod rings;
pub mod solvation;
pub mod spectroscopy;
pub mod stereo;
pub mod system;
pub mod transport;

pub mod alpha;
pub mod beta;

#[cfg(feature = "alpha-dynamics-live")]
pub mod dynamics_live;
