//! Phase 4: Spectroscopy Module
//!
//! Implements spectroscopic property calculations from converged SCF results:
//!
//! - **UV-Vis**: Simplified Tamm-Dancoff Approximation (sTDA)
//! - **IR**: Semi-numerical Hessian → normal modes → IR intensities
//! - **NMR**: Gauge-Including Atomic Orbitals (GIAO) for chemical shifts
//! - **Transition dipoles**: Electric dipole transition moments

pub mod giao_nmr;
pub mod hessian;
pub mod ir_intensities;
pub mod stda_uvvis;
pub mod transition_dipoles;
