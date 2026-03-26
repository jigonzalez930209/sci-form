//! NMR Spectroscopy module: chemical shift prediction and spectrum generation.
//!
//! Implements Phase D3 of the spectroscopy roadmap:
//! - HOSE code generation for atomic environment characterization
//! - Empirical chemical shift prediction for ¹H and ¹³C
//! - J-coupling estimation via Karplus equation
//! - Lorentzian-broadened NMR spectrum generation

pub mod coupling;
pub mod hose;
pub mod nucleus;
pub mod shifts;
pub mod spectrum;

#[cfg(feature = "parallel")]
pub use coupling::ensemble_averaged_j_couplings_parallel;
pub use coupling::{ensemble_averaged_j_couplings, predict_j_couplings, JCoupling, KarplusParams};
pub use hose::{HoseCode, HoseShiftLookup};
pub use nucleus::NmrNucleus;
pub use shifts::{
    predict_chemical_shifts, predict_chemical_shifts_for_nucleus, ChemicalShift, NmrShiftResult,
    NucleusShiftSeries,
};
pub use spectrum::{compute_nmr_spectrum, NmrPeak, NmrSpectrum, PeakIntegration};
