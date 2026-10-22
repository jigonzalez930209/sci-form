//! E11: Growing String Method (GSM)
//!
//! Reaction path and transition state finder driven by force-field gradients.
//! Grows a string of images from reactant and product geometries toward
//! the transition state without requiring an initial path guess.
//!
//! The method follows the Growing String Method of Behler and Carpenter
//! (J. Chem. Phys. 1992, 96, 6468): a string of images is grown from both
//! ends, re-parametrised by arc length at every iteration and relaxed
//! perpendicular to the string until the highest image converges to the
//! transition state.
//!
//! # Stability tier
//!
//! **Alpha** — early research tier. Gated behind the `alpha-gsm` feature
//! flag (alias: `experimental-gsm`). The validation gate is
//! `tests/alpha_gsm_validation.rs` (F82): the transition state is compared
//! against ASE's climbing-image NEB on the Müller–Brown potential, and the
//! promotion criteria are listed in [`crate::experimental_status`].
//!
//! # Units
//!
//! Coordinates are flat `3N` arrays in **Å**, energies are in **kcal/mol**
//! and gradients in **kcal/(mol·Å)**.

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
