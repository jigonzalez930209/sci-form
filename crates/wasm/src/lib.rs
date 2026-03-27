//! sci-form WASM bindings.
//!
//! All public functions are split into focused modules.

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
