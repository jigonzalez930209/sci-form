//! NMR spectrum generation with Lorentzian broadening.
//!
//! Combines predicted chemical shifts and J-couplings to produce
//! a synthetic 1D NMR spectrum with Lorentzian line shapes.

use serde::{Deserialize, Serialize};

use super::coupling::JCoupling;
use super::nucleus::NmrNucleus;
use super::shifts::{ChemicalShift, NmrShiftResult};

#[derive(Debug, Clone)]
struct ShiftGroup {
    shift_ppm: f64,
    environment: String,
    atom_indices: Vec<usize>,
    confidence: f64,
    exchangeable: bool,
}

#[derive(Debug, Clone)]
struct CouplingGroup {
    n_equivalent_neighbors: usize,
    average_j_hz: f64,
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
    /// All equivalent atom indices contributing to this signal.
    pub atom_indices: Vec<usize>,
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
    /// Peak group integrations (relative areas, normalized to largest = 1.0).
    pub integrations: Vec<PeakIntegration>,
    /// Notes.
    pub notes: Vec<String>,
}

/// Integration result for a group of equivalent peaks.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PeakIntegration {
    /// Center of the peak group (ppm).
    pub center_ppm: f64,
    /// Integration bounds (ppm_low, ppm_high).
    pub bounds: (f64, f64),
    /// Raw area (sum of intensities × step).
    pub raw_area: f64,
    /// Relative area (normalized so largest group = 1.0).
    pub relative_area: f64,
    /// Number of equivalent atoms contributing.
    pub n_atoms: usize,
}

/// Default FWHM (full-width at half-maximum) for each nucleus type, in Hz.
fn default_fwhm_hz(nucleus: NmrNucleus) -> f64 {
    nucleus.default_fwhm_hz()
}

/// Default spectrometer frequency for each nucleus (MHz).
fn default_frequency_mhz(nucleus: NmrNucleus) -> f64 {
    nucleus.default_frequency_mhz()
}

/// Lorentzian line shape: L(x) = (1/π) · γ / [(x - x₀)² + γ²]
fn lorentzian(x: f64, x0: f64, gamma: f64) -> f64 {
    let g = gamma.max(1e-8);
    let dx = x - x0;
    g / (std::f64::consts::PI * (dx * dx + g * g))
}

/// Determine multiplicity from J-coupling count.
#[cfg(test)]
fn multiplicity_label(n_couplings: usize) -> &'static str {
    match n_couplings {
        0 => "s", // singlet
        1 => "d", // doublet
        2 => "t", // triplet
        3 => "q", // quartet
        _ => "m", // multiplet
    }
}

fn shift_group_tolerance(nucleus: NmrNucleus) -> f64 {
    match nucleus {
        NmrNucleus::H1 | NmrNucleus::H2 | NmrNucleus::H3 => 1e-3,
        NmrNucleus::C13 => 1e-2,
        _ => 5e-3,
    }
}

fn is_exchangeable_environment(environment: &str) -> bool {
    environment.contains("O-H") || environment.contains("N-H") || environment.contains("S-H")
}

fn build_shift_groups(active_shifts: &[ChemicalShift], nucleus: NmrNucleus) -> Vec<ShiftGroup> {
    let tolerance = shift_group_tolerance(nucleus);
    let mut groups: Vec<ShiftGroup> = Vec::new();

    for shift in active_shifts {
        if let Some(group) = groups.iter_mut().find(|group| {
            group.environment == shift.environment
                && (group.shift_ppm - shift.shift_ppm).abs() <= tolerance
        }) {
            let count = group.atom_indices.len() as f64;
            group.shift_ppm = (group.shift_ppm * count + shift.shift_ppm) / (count + 1.0);
            group.confidence = group.confidence.min(shift.confidence);
            group.atom_indices.push(shift.atom_index);
        } else {
            groups.push(ShiftGroup {
                shift_ppm: shift.shift_ppm,
                environment: shift.environment.clone(),
                atom_indices: vec![shift.atom_index],
                confidence: shift.confidence,
                exchangeable: is_exchangeable_environment(&shift.environment),
            });
        }
    }

    groups.sort_by(|left, right| {
        left.shift_ppm
            .partial_cmp(&right.shift_ppm)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    groups
}

fn build_coupling_groups(
    source_group: &ShiftGroup,
    all_groups: &[ShiftGroup],
    couplings: &[JCoupling],
) -> Vec<CouplingGroup> {
    if source_group.exchangeable {
        return Vec::new();
    }

    let representative = match source_group.atom_indices.first() {
        Some(index) => *index,
        None => return Vec::new(),
    };

    let mut groups = Vec::new();

    for target_group in all_groups {
        if std::ptr::eq(source_group, target_group) || target_group.exchangeable {
            continue;
        }

        let target_couplings: Vec<f64> = couplings
            .iter()
            .filter(|coupling| coupling.n_bonds == 3)
            .filter_map(|coupling| {
                if (coupling.h1_index == representative
                    && target_group.atom_indices.contains(&coupling.h2_index))
                    || (coupling.h2_index == representative
                        && target_group.atom_indices.contains(&coupling.h1_index))
                {
                    Some(coupling.j_hz.abs())
                } else {
                    None
                }
            })
            .filter(|j_hz| *j_hz >= 0.5)
            .collect();

        if target_couplings.is_empty() {
            continue;
        }

        let average_j_hz = target_couplings.iter().sum::<f64>() / target_couplings.len() as f64;
        groups.push(CouplingGroup {
            n_equivalent_neighbors: target_couplings.len(),
            average_j_hz,
        });
    }

    groups.sort_by(|left, right| {
        right
            .average_j_hz
            .partial_cmp(&left.average_j_hz)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    groups
}

fn merge_nearby_lines(lines: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
    if lines.is_empty() {
        return lines;
    }

    let mut sorted = lines;
    sorted.sort_by(|left, right| {
        left.0
            .partial_cmp(&right.0)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut merged: Vec<(f64, f64)> = Vec::new();
    for (position, intensity) in sorted {
        if let Some((prev_position, prev_intensity)) = merged.last_mut() {
            if (position - *prev_position).abs() <= 1e-6 {
                let total = *prev_intensity + intensity;
                *prev_position = (*prev_position * *prev_intensity + position * intensity) / total;
                *prev_intensity = total;
                continue;
            }
        }
        merged.push((position, intensity));
    }

    merged
}

fn split_lines(
    center_ppm: f64,
    total_intensity: f64,
    coupling_groups: &[CouplingGroup],
    spectrometer_frequency_mhz: f64,
) -> Vec<(f64, f64)> {
    let mut lines = vec![(center_ppm, total_intensity)];

    for group in coupling_groups {
        let j_ppm = group.average_j_hz / spectrometer_frequency_mhz.max(1e-6);
        let coeffs = pascal_row(group.n_equivalent_neighbors);
        let coeff_sum = coeffs.iter().sum::<f64>().max(1e-12);
        let mut split = Vec::with_capacity(lines.len() * coeffs.len());

        for (line_center, line_intensity) in &lines {
            for (index, coeff) in coeffs.iter().enumerate() {
                let offset = (index as f64 - group.n_equivalent_neighbors as f64 / 2.0) * j_ppm;
                split.push((line_center + offset, line_intensity * coeff / coeff_sum));
            }
        }

        lines = merge_nearby_lines(split);
    }

    lines
}

fn multiplicity_code(n_equivalent_neighbors: usize) -> &'static str {
    match n_equivalent_neighbors {
        0 => "s",
        1 => "d",
        2 => "t",
        3 => "q",
        4 => "quint",
        _ => "m",
    }
}

fn format_multiplicity(coupling_groups: &[CouplingGroup]) -> String {
    if coupling_groups.is_empty() {
        return "s".to_string();
    }

    let mut codes: Vec<&str> = coupling_groups
        .iter()
        .map(|group| multiplicity_code(group.n_equivalent_neighbors))
        .collect();
    let base = if codes.len() == 1 {
        codes.remove(0).to_string()
    } else if codes.len() == 2 && codes.iter().all(|code| *code != "m" && *code != "quint") {
        codes.join("")
    } else {
        "m".to_string()
    };

    let j_values = coupling_groups
        .iter()
        .map(|group| format!("{:.1}", group.average_j_hz))
        .collect::<Vec<_>>()
        .join(", ");

    format!("{} (J≈{} Hz)", base, j_values)
}

/// Generate an NMR spectrum from shift predictions and J-couplings.
///
/// `shifts`: predicted chemical shifts
/// `couplings`: predicted J-coupling constants
/// `nucleus`: which nucleus to display
/// `gamma`: Lorentzian line width in ppm (if 0.0, uses nucleus-specific FWHM default)
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
    compute_nmr_spectrum_for_shifts(
        shifts.shifts_for_nucleus(nucleus),
        couplings,
        nucleus,
        gamma,
        ppm_min,
        ppm_max,
        n_points,
    )
}

pub fn compute_nmr_spectrum_for_shifts(
    active_shifts: &[ChemicalShift],
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

    // Use nucleus-specific FWHM if gamma is 0 or very small
    let effective_gamma = if gamma > 1e-6 {
        gamma
    } else {
        default_fwhm_hz(nucleus) / default_frequency_mhz(nucleus)
    };

    let shift_groups = build_shift_groups(active_shifts, nucleus);
    let spectrometer_frequency_mhz = default_frequency_mhz(nucleus);
    let mut peaks = Vec::with_capacity(shift_groups.len());

    for group in &shift_groups {
        let peak_intensity = group.atom_indices.len() as f64;
        let coupling_groups = if matches!(nucleus, NmrNucleus::H1) {
            build_coupling_groups(group, &shift_groups, couplings)
        } else {
            Vec::new()
        };

        peaks.push(NmrPeak {
            shift_ppm: group.shift_ppm,
            intensity: peak_intensity,
            atom_index: group.atom_indices[0],
            atom_indices: group.atom_indices.clone(),
            multiplicity: format_multiplicity(&coupling_groups),
            environment: group.environment.clone(),
        });

        let lines = if matches!(nucleus, NmrNucleus::H1) {
            split_lines(
                group.shift_ppm,
                peak_intensity,
                &coupling_groups,
                spectrometer_frequency_mhz,
            )
        } else {
            vec![(group.shift_ppm, peak_intensity)]
        };

        for (line_ppm, line_intensity) in lines {
            for (idx, &ppm) in ppm_axis.iter().enumerate() {
                intensities[idx] += line_intensity * lorentzian(ppm, line_ppm, effective_gamma);
            }
        }
    }

    // Compute integrations for each peak group
    let integration_width = effective_gamma * 10.0; // integrate ±10γ around each peak
    let integrations =
        compute_integrations(&peaks, &ppm_axis, &intensities, step, integration_width);

    let mut notes = vec![
        format!(
            "{} NMR spectrum with Lorentzian broadening (γ = {:.4} ppm, FWHM = {:.1} Hz).",
            nucleus.unicode_label(),
            effective_gamma,
            effective_gamma * default_frequency_mhz(nucleus)
        ),
        "Chemical shifts come from the fast empirical inference layer. ¹H uses explicit vicinal J-coupling splitting; other nuclei are rendered as screening-level singlets unless a dedicated model is available.".to_string(),
        "Equivalent nuclei are grouped before rendering, and exchangeable O-H/N-H/S-H couplings are suppressed in the default 1H spectrum.".to_string(),
        "First-order splitting uses explicit coupling groups; higher-order effects (roofing, strong coupling) are not modeled.".to_string(),
    ];
    if nucleus.is_quadrupolar() {
        notes.push(
            "Quadrupolar nucleus selected: linewidths and positions are approximate relative indicators, not quantitative simulations of relaxation or isotope abundance.".to_string(),
        );
    }

    NmrSpectrum {
        ppm_axis,
        intensities,
        peaks,
        nucleus,
        gamma: effective_gamma,
        integrations,
        notes,
    }
}

/// Compute peak integrations for groups of peaks in the spectrum.
fn compute_integrations(
    peaks: &[NmrPeak],
    ppm_axis: &[f64],
    intensities: &[f64],
    step: f64,
    integration_width: f64,
) -> Vec<PeakIntegration> {
    if peaks.is_empty() || ppm_axis.is_empty() {
        return Vec::new();
    }

    let mut integrations: Vec<PeakIntegration> = peaks
        .iter()
        .map(|peak| {
            let low = peak.shift_ppm - integration_width;
            let high = peak.shift_ppm + integration_width;

            let raw_area: f64 = ppm_axis
                .iter()
                .zip(intensities.iter())
                .filter(|(&ppm, _)| ppm >= low && ppm <= high)
                .map(|(_, &intensity)| intensity * step)
                .sum();

            PeakIntegration {
                center_ppm: peak.shift_ppm,
                bounds: (low, high),
                raw_area,
                relative_area: 0.0,
                n_atoms: peak.intensity.round().max(1.0) as usize,
            }
        })
        .collect();

    // Normalize relative areas using the equivalent-atom count to keep areas chemically meaningful.
    let max_atoms = integrations
        .iter()
        .map(|integration| integration.n_atoms)
        .max()
        .unwrap_or(1);
    let max_area = integrations
        .iter()
        .map(|i| i.raw_area)
        .fold(0.0f64, f64::max);

    if max_atoms > 0 {
        for int in &mut integrations {
            int.relative_area = int.n_atoms as f64 / max_atoms as f64;
            if max_area <= 1e-30 {
                int.raw_area = int.n_atoms as f64;
            }
        }
    }

    integrations
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

    #[test]
    fn test_shift_grouping_and_integrations_for_ethanol() {
        let spectrum = crate::compute_nmr_spectrum("CCO", "1H", 0.02, 0.0, 12.0, 1000).unwrap();

        assert_eq!(
            spectrum.peaks.len(),
            3,
            "ethanol should collapse to CH3, CH2, and OH groups"
        );

        let mut atom_counts: Vec<usize> = spectrum
            .integrations
            .iter()
            .map(|integration| integration.n_atoms)
            .collect();
        atom_counts.sort_unstable();
        assert_eq!(atom_counts, vec![1, 2, 3]);

        assert!(
            spectrum
                .peaks
                .iter()
                .any(|peak| peak.environment.contains("methyl")
                    && peak.multiplicity.starts_with('t')),
            "methyl group should appear as a triplet-like signal: {:?}",
            spectrum
                .peaks
                .iter()
                .map(|peak| (&peak.environment, &peak.multiplicity))
                .collect::<Vec<_>>()
        );
        assert!(
            spectrum
                .peaks
                .iter()
                .any(|peak| peak.environment.contains("methylene")
                    && peak.multiplicity.starts_with('q')),
            "methylene group should appear as a quartet-like signal: {:?}",
            spectrum
                .peaks
                .iter()
                .map(|peak| (&peak.environment, &peak.multiplicity))
                .collect::<Vec<_>>()
        );
        assert!(
            spectrum
                .peaks
                .iter()
                .any(|peak| peak.environment.contains("O-H") && peak.multiplicity == "s"),
            "exchangeable alcohol proton should default to a singlet: {:?}",
            spectrum
                .peaks
                .iter()
                .map(|peak| (&peak.environment, &peak.multiplicity))
                .collect::<Vec<_>>()
        );
    }
}
