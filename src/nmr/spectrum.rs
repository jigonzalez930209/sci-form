//! NMR spectrum generation with Lorentzian broadening.
//!
//! Combines predicted chemical shifts and J-couplings to produce
//! a synthetic 1D NMR spectrum with Lorentzian line shapes.

use serde::{Deserialize, Serialize};

use super::coupling::JCoupling;
use super::shifts::{ChemicalShift, NmrShiftResult};

/// Target NMR nucleus.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum NmrNucleus {
    /// ¹H NMR
    H1,
    /// ¹³C NMR
    C13,
    /// ¹⁹F NMR
    F19,
    /// ³¹P NMR
    P31,
    /// ¹⁵N NMR
    N15,
    /// ¹¹B NMR
    B11,
    /// ²⁹Si NMR
    Si29,
    /// ⁷⁷Se NMR
    Se77,
    /// ¹⁷O NMR
    O17,
    /// ³³S NMR
    S33,
}

/// A single peak in the NMR spectrum.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NmrPeak {
    /// Chemical shift in ppm.
    pub shift_ppm: f64,
    /// Relative intensity (integration).
    pub intensity: f64,
    /// Atom index.
    pub atom_index: usize,
    /// Multiplicity label (s, d, t, q, m).
    pub multiplicity: String,
    /// Environment description.
    pub environment: String,
}

/// Complete NMR spectrum.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NmrSpectrum {
    /// Chemical shift axis (ppm), typically decreasing.
    pub ppm_axis: Vec<f64>,
    /// Intensity axis.
    pub intensities: Vec<f64>,
    /// Identified peaks.
    pub peaks: Vec<NmrPeak>,
    /// Nucleus type.
    pub nucleus: NmrNucleus,
    /// Line width parameter (ppm).
    pub gamma: f64,
    /// Notes.
    pub notes: Vec<String>,
}

/// Lorentzian line shape: L(x) = (1/π) · γ / [(x - x₀)² + γ²]
fn lorentzian(x: f64, x0: f64, gamma: f64) -> f64 {
    let g = gamma.max(1e-8);
    let dx = x - x0;
    g / (std::f64::consts::PI * (dx * dx + g * g))
}

/// Determine multiplicity from J-coupling count.
fn multiplicity_label(n_couplings: usize) -> &'static str {
    match n_couplings {
        0 => "s", // singlet
        1 => "d", // doublet
        2 => "t", // triplet
        3 => "q", // quartet
        _ => "m", // multiplet
    }
}

/// Generate an NMR spectrum from shift predictions and J-couplings.
///
/// `shifts`: predicted chemical shifts
/// `couplings`: predicted J-coupling constants
/// `nucleus`: which nucleus to display
/// `gamma`: Lorentzian line width in ppm (typically 0.01–0.05 for ¹H, 0.5–2.0 for ¹³C)
/// `ppm_min`, `ppm_max`: spectral window
/// `n_points`: grid resolution
pub fn compute_nmr_spectrum(
    shifts: &NmrShiftResult,
    couplings: &[JCoupling],
    nucleus: NmrNucleus,
    gamma: f64,
    ppm_min: f64,
    ppm_max: f64,
    n_points: usize,
) -> NmrSpectrum {
    let n_points = n_points.max(2);
    let step = (ppm_max - ppm_min) / (n_points as f64 - 1.0);
    // NMR convention: ppm axis runs from high to low
    let ppm_axis: Vec<f64> = (0..n_points).map(|i| ppm_max - step * i as f64).collect();
    let mut intensities = vec![0.0; n_points];

    let active_shifts: &[ChemicalShift] = match nucleus {
        NmrNucleus::H1 => &shifts.h_shifts,
        NmrNucleus::C13 => &shifts.c_shifts,
        NmrNucleus::F19 => &shifts.f_shifts,
        NmrNucleus::P31 => &shifts.p_shifts,
        NmrNucleus::N15 => &shifts.n_shifts,
        NmrNucleus::B11 => &shifts.b_shifts,
        NmrNucleus::Si29 => &shifts.si_shifts,
        NmrNucleus::Se77 => &shifts.se_shifts,
        NmrNucleus::O17 => &shifts.o_shifts,
        NmrNucleus::S33 => &shifts.s_shifts,
    };

    let mut peaks = Vec::with_capacity(active_shifts.len());

    for shift in active_shifts {
        // Count J-couplings involving this atom
        let n_j = if matches!(nucleus, NmrNucleus::H1) {
            couplings
                .iter()
                .filter(|c| c.h1_index == shift.atom_index || c.h2_index == shift.atom_index)
                .filter(|c| c.n_bonds == 3) // only vicinal for splitting
                .count()
        } else {
            0 // Non-¹H nuclei typically shown as singlets in decoupled spectra
        };

        let mult = multiplicity_label(n_j);
        let intensity = 1.0; // Each H contributes equally

        peaks.push(NmrPeak {
            shift_ppm: shift.shift_ppm,
            intensity,
            atom_index: shift.atom_index,
            multiplicity: mult.to_string(),
            environment: shift.environment.clone(),
        });

        // Simple splitting pattern (first-order only)
        if n_j == 0 || !matches!(nucleus, NmrNucleus::H1) {
            // Singlet or non-¹H: single Lorentzian
            for (idx, &ppm) in ppm_axis.iter().enumerate() {
                intensities[idx] += intensity * lorentzian(ppm, shift.shift_ppm, gamma);
            }
        } else {
            // For ¹H with coupling: generate split pattern
            // Average J-coupling for this atom
            let avg_j: f64 = couplings
                .iter()
                .filter(|c| c.h1_index == shift.atom_index || c.h2_index == shift.atom_index)
                .filter(|c| c.n_bonds == 3)
                .map(|c| c.j_hz)
                .sum::<f64>()
                / n_j.max(1) as f64;

            // Convert J from Hz to ppm (assume 400 MHz spectrometer)
            let j_ppm = avg_j / 400.0;

            // Generate Pascal's triangle splitting
            let n_lines = n_j + 1;
            let coeffs = pascal_row(n_j);
            let total: f64 = coeffs.iter().sum::<f64>();

            for (k, &coeff) in coeffs.iter().enumerate() {
                let offset = (k as f64 - n_j as f64 / 2.0) * j_ppm;
                let line_ppm = shift.shift_ppm + offset;
                let line_intensity = intensity * coeff / total;

                for (idx, &ppm) in ppm_axis.iter().enumerate() {
                    intensities[idx] += line_intensity * lorentzian(ppm, line_ppm, gamma);
                }
            }

            // Update peak multiplicity
            if let Some(p) = peaks.last_mut() {
                p.multiplicity = format!("{} (n+1={}, J≈{:.1} Hz)", mult, n_lines, avg_j);
            }
        }
    }

    let nucleus_label = match nucleus {
        NmrNucleus::H1 => "¹H",
        NmrNucleus::C13 => "¹³C",
        NmrNucleus::F19 => "¹⁹F",
        NmrNucleus::P31 => "³¹P",
        NmrNucleus::N15 => "¹⁵N",
        NmrNucleus::B11 => "¹¹B",
        NmrNucleus::Si29 => "²⁹Si",
        NmrNucleus::Se77 => "⁷⁷Se",
        NmrNucleus::O17 => "¹⁷O",
        NmrNucleus::S33 => "³³S",
    };

    NmrSpectrum {
        ppm_axis,
        intensities,
        peaks,
        nucleus,
        gamma,
        notes: vec![
            format!(
                "{} NMR spectrum with Lorentzian broadening (γ = {} ppm).",
                nucleus_label, gamma
            ),
            "Chemical shifts from empirical additivity rules; J-couplings from Karplus equation.".to_string(),
            "First-order splitting only; higher-order effects (roofing, strong coupling) not modeled.".to_string(),
        ],
    }
}

/// Pascal's triangle row n: [C(n,0), C(n,1), ..., C(n,n)].
fn pascal_row(n: usize) -> Vec<f64> {
    let mut row = vec![1.0];
    for k in 1..=n {
        let prev = row[k - 1];
        row.push(prev * (n - k + 1) as f64 / k as f64);
    }
    row
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lorentzian_normalization() {
        // Integral of Lorentzian should be 1 (approximately)
        let gamma = 0.02;
        let n = 10000;
        let dx = 20.0 / n as f64;
        let integral: f64 = (0..n)
            .map(|i| {
                let x = -10.0 + i as f64 * dx;
                lorentzian(x, 0.0, gamma) * dx
            })
            .sum();
        assert!(
            (integral - 1.0).abs() < 0.01,
            "Lorentzian integral = {}, expected ~1.0",
            integral
        );
    }

    #[test]
    fn test_pascal_row() {
        assert_eq!(pascal_row(0), vec![1.0]);
        assert_eq!(pascal_row(1), vec![1.0, 1.0]);
        assert_eq!(pascal_row(2), vec![1.0, 2.0, 1.0]);
        assert_eq!(pascal_row(3), vec![1.0, 3.0, 3.0, 1.0]);
    }

    #[test]
    fn test_nmr_spectrum_h1_ethanol() {
        let mol = crate::graph::Molecule::from_smiles("CCO").unwrap();
        let shifts = super::super::shifts::predict_chemical_shifts(&mol);
        let couplings = super::super::coupling::predict_j_couplings(&mol, &[]);

        let spectrum =
            compute_nmr_spectrum(&shifts, &couplings, NmrNucleus::H1, 0.02, 0.0, 12.0, 1000);

        assert_eq!(spectrum.ppm_axis.len(), 1000);
        assert_eq!(spectrum.intensities.len(), 1000);
        assert!(!spectrum.peaks.is_empty());

        // Verify ppm axis goes from high to low (NMR convention)
        assert!(spectrum.ppm_axis[0] > spectrum.ppm_axis[999]);
    }

    #[test]
    fn test_nmr_spectrum_c13_benzene() {
        let mol = crate::graph::Molecule::from_smiles("c1ccccc1").unwrap();
        let shifts = super::super::shifts::predict_chemical_shifts(&mol);
        let couplings = super::super::coupling::predict_j_couplings(&mol, &[]);

        let spectrum =
            compute_nmr_spectrum(&shifts, &couplings, NmrNucleus::C13, 1.0, 0.0, 220.0, 1000);

        assert!(!spectrum.peaks.is_empty(), "Benzene should have ¹³C peaks");

        // All aromatic C should cluster near 128 ppm
        for peak in &spectrum.peaks {
            assert!(
                (peak.shift_ppm - 128.0).abs() < 20.0,
                "Benzene ¹³C peak at {} ppm should be near 128",
                peak.shift_ppm
            );
        }
    }

    #[test]
    fn test_multiplicity_labels() {
        assert_eq!(multiplicity_label(0), "s");
        assert_eq!(multiplicity_label(1), "d");
        assert_eq!(multiplicity_label(2), "t");
        assert_eq!(multiplicity_label(3), "q");
        assert_eq!(multiplicity_label(4), "m");
    }
}
