//! Advanced Spectroscopy — Core Module
//!
//! Contains spectroscopic methods promoted from experimental:
//! - **sTDA UV-Vis**: Simplified Tamm-Dancoff approximation for electronic excitations
//! - **GIAO NMR**: Gauge-including atomic orbitals for chemical shifts
//!
//! These complement the existing topology-based NMR (HOSE codes) and
//! numerical IR vibrational analysis in the core.
//!
//! ## sTDA UV-Vis
//!
//! Requires a converged SCF result (from HF-3c, PM3+basis, or EHT).
//! Provides excitation energies, wavelengths, and oscillator strengths.
//!
//! ## GIAO NMR
//!
//! Quantum-mechanical NMR chemical shifts. More accurate than topological
//! (HOSE code) NMR but requires 3D coordinates and a SCF calculation.
//! The user can choose fast topological NMR or slower quantum NMR.

mod giao_nmr;
pub mod hessian;
pub mod ir_intensities;
mod stda_uvvis;
pub mod transition_dipoles;
mod types;

pub use giao_nmr::{
    compute_nmr_shieldings, compute_nmr_shieldings_for_nucleus, shieldings_to_shifts,
};
pub use stda_uvvis::{compute_stda, StdaConfig};
pub use types::{
    NmrShieldingResult, ScfInput, ShieldingTensor, SpectroscopyResult, TransitionInfo,
};
