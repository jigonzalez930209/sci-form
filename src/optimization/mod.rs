//! Generic numerical optimizers used internally by force fields and geometry refinement.
//!
//! Currently provides a BFGS quasi-Newton minimizer ([`bfgs::RustBfgsEngine`])
//! with backtracking line search, used by conformer refinement and bounds-based
//! minimization.

pub mod bfgs;
