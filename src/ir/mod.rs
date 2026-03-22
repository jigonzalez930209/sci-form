//! IR Spectroscopy module: numerical Hessian, vibrational frequencies, and IR intensities.
//!
//! Implements Phase D2 of the spectroscopy roadmap:
//! - Numerical Hessian via central finite differences (6N energy evaluations)
//! - Mass-weighted Hessian diagonalization for normal modes and frequencies
//! - IR intensities from numerical dipole derivatives along normal modes
//! - Lorentzian-broadened IR spectrum generation

pub mod hessian;
pub mod peak_assignment;
pub mod vibrations;

pub use hessian::{compute_numerical_hessian, compute_uff_analytical_hessian, HessianMethod};
pub use peak_assignment::{assign_peaks, AssignmentResult, PeakAssignment};
pub use vibrations::{
    compute_ir_spectrum, compute_ir_spectrum_with_broadening, compute_vibrational_analysis,
    compute_vibrational_analysis_uff, BroadeningType, IrPeak, IrSpectrum, Thermochemistry,
    VibrationalAnalysis, VibrationalMode,
};
