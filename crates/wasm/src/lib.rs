//! sci-form WASM bindings.
//!
//! All public functions are split into focused modules.

#[cfg(feature = "parallel")]
pub use wasm_bindgen_rayon::init_thread_pool;

mod helpers;

pub mod dynamics;
pub mod eht;
pub mod electronic;
pub mod embed;
pub mod esp;
pub mod forcefield;
pub mod materials;
pub mod mesh;
pub mod ml;
pub mod properties;
pub mod reactivity;
pub mod spectroscopy;
pub mod system;
pub mod transport;
