//! E11: Growing String Method (GSM)
//!
//! Reaction path and transition state finder driven by force-field gradients.
//! Grows a string of images from reactant and product geometries toward
//! the transition state without requiring an initial path guess.

pub mod saddle;
pub mod string;

pub use saddle::{find_transition_state, refine_saddle, GsmResult};
pub use string::{grow_string, interpolate_node, GsmConfig, GsmPath};
