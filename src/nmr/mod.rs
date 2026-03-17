//! NMR Spectroscopy module: chemical shift prediction and spectrum generation.
//!
//! Implements Phase D3 of the spectroscopy roadmap:
//! - HOSE code generation for atomic environment characterization
//! - Empirical chemical shift prediction for ¹H and ¹³C
//! - J-coupling estimation via Karplus equation
//! - Lorentzian-broadened NMR spectrum generation

pub mod coupling;
pub mod hose;
pub mod shifts;
pub mod spectrum;

pub use coupling::{predict_j_couplings, JCoupling};
pub use hose::HoseCode;
pub use shifts::{predict_chemical_shifts, ChemicalShift, NmrShiftResult};
pub use spectrum::{compute_nmr_spectrum, NmrNucleus, NmrPeak, NmrSpectrum};
