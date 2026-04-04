//! E11: Growing String Method (GSM)
//!
//! Reaction path and transition state finder driven by force-field gradients.
//! Grows a string of images from reactant and product geometries toward
//! the transition state without requiring an initial path guess.

pub mod backend;
pub mod saddle;
pub mod string;

pub use backend::{
    compare_gsm_backends, evaluate_gsm_backend, find_transition_state_with_backend,
    plan_gsm_backends, plan_gsm_backends_for_smiles, GsmBackendCapability, GsmBackendEvaluation,
    GsmBackendStatus, GsmEnergyBackend,
};
pub use saddle::{find_transition_state, refine_saddle, GsmResult};
pub use string::{grow_string, interpolate_node, GsmConfig, GsmPath};
