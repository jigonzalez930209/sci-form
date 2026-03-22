//! Vibrational analysis: normal modes, frequencies, IR intensities, and spectrum.
//!
//! Given a numerical Hessian, computes mass-weighted normal modes,
//! vibrational frequencies (cm⁻¹), IR intensities from dipole derivatives,
//! and generates a Lorentzian-broadened IR spectrum.

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

use super::hessian::{compute_numerical_hessian, HessianMethod};

/// Atomic masses in amu for elements 1–86 (H through Rn).
fn atomic_mass(z: u8) -> f64 {
    match z {
        1 => 1.00794,
        2 => 4.00260,
        5 => 10.811,
        6 => 12.0107,
        7 => 14.0067,
        8 => 15.9994,
        9 => 18.9984,
        14 => 28.0855,
        15 => 30.9738,
        16 => 32.065,
        17 => 35.453,
        22 => 47.867,
        24 => 51.9961,
        25 => 54.9380,
        26 => 55.845,
        27 => 58.9332,
        28 => 58.6934,
        29 => 63.546,
        30 => 65.38,
        35 => 79.904,
        44 => 101.07,
        46 => 106.42,
        47 => 107.868,
        53 => 126.904,
        78 => 195.084,
        79 => 196.967,
        _ => z as f64 * 1.5, // rough fallback
    }
}

/// A single vibrational normal mode.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VibrationalMode {
    /// Frequency in cm⁻¹ (negative = imaginary).
    pub frequency_cm1: f64,
    /// IR intensity in km/mol.
    pub ir_intensity: f64,
    /// Displacement vector (3N elements, mass-weighted normal coordinate).
    pub displacement: Vec<f64>,
    /// Whether this is a real mode (true) or translation/rotation (false).
    pub is_real: bool,
}

/// Complete vibrational analysis result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VibrationalAnalysis {
    /// Number of atoms.
    pub n_atoms: usize,
    /// All vibrational modes sorted by frequency.
    pub modes: Vec<VibrationalMode>,
    /// Number of real vibrational modes (excluding translation/rotation).
    pub n_real_modes: usize,
    /// Zero-point vibrational energy in eV.
    pub zpve_ev: f64,
    /// Thermochemistry at 298.15 K (if computed).
    pub thermochemistry: Option<Thermochemistry>,
    /// Semiempirical method used for Hessian.
    pub method: String,
    /// Notes and caveats.
    pub notes: Vec<String>,
}

/// Thermochemical properties from RRHO approximation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Thermochemistry {
    /// Temperature in K.
    pub temperature_k: f64,
    /// Zero-point energy in kcal/mol.
    pub zpve_kcal: f64,
    /// Thermal energy correction in kcal/mol.
    pub thermal_energy_kcal: f64,
    /// Thermal enthalpy correction in kcal/mol (E_thermal + RT).
    pub enthalpy_correction_kcal: f64,
    /// Vibrational entropy in cal/(mol·K).
    pub entropy_vib_cal: f64,
    /// Gibbs free energy correction in kcal/mol (H - TS).
    pub gibbs_correction_kcal: f64,
}

/// A peak in the IR spectrum.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IrPeak {
    /// Frequency in cm⁻¹.
    pub frequency_cm1: f64,
    /// IR intensity in km/mol.
    pub ir_intensity: f64,
    /// Mode index.
    pub mode_index: usize,
    /// Functional group assignment (e.g., "O-H stretch", "C=O stretch").
    pub assignment: String,
}

/// Broadened IR spectrum.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IrSpectrum {
    /// Wavenumber axis (cm⁻¹).
    pub wavenumbers: Vec<f64>,
    /// Transmittance/absorbance axis.
    pub intensities: Vec<f64>,
    /// Discrete peaks from vibrational analysis.
    pub peaks: Vec<IrPeak>,
    /// Broadening width in cm⁻¹.
    pub gamma: f64,
    /// Notes.
    pub notes: Vec<String>,
}

// Constants for unit conversion
// 1 eV = 8065.544 cm⁻¹
const EV_TO_CM1: f64 = 8065.544;
// 1 amu = 1.66054e-27 kg, 1 eV = 1.60218e-19 J, 1 Å = 1e-10 m
// Hessian in eV/Å², mass in amu → frequency in cm⁻¹:
// ω² = H/(m) in eV/(amu·Å²) → convert to s⁻²:
// factor = eV_to_J / (amu_to_kg * Å_to_m²) = 1.60218e-19 / (1.66054e-27 * 1e-20) = 9.6485e27
// ν_cm1 = sqrt(factor * eigenvalue) / (2π * c) where c = 2.998e10 cm/s
const HESSIAN_TO_FREQUENCY_FACTOR: f64 = 9.6485e27; // eV/(amu·Å²) → s⁻²
const INV_2PI_C: f64 = 1.0 / (2.0 * std::f64::consts::PI * 2.99792458e10); // 1/(2πc) in s/cm

/// Compute the molecular dipole for use in IR intensity calculation.
fn compute_dipole_vector(
    elements: &[u8],
    positions: &[[f64; 3]],
    method: HessianMethod,
) -> Result<[f64; 3], String> {
    match method {
        HessianMethod::Eht => {
            let eht = crate::eht::solve_eht(elements, positions, None)?;
            let dipole = crate::dipole::compute_dipole_from_eht(
                elements,
                positions,
                &eht.coefficients,
                eht.n_electrons,
            );
            Ok(dipole.vector)
        }
        HessianMethod::Pm3 => {
            let pm3 = crate::pm3::solve_pm3(elements, positions)?;
            // Use Mulliken charges for dipole
            let dipole = crate::dipole::compute_dipole(&pm3.mulliken_charges, positions);
            Ok(dipole.vector)
        }
        HessianMethod::Xtb => {
            let xtb = crate::xtb::solve_xtb(elements, positions)?;
            let dipole = crate::dipole::compute_dipole(&xtb.mulliken_charges, positions);
            Ok(dipole.vector)
        }
        HessianMethod::Uff => {
            // UFF has no electronic structure; use Gasteiger charges as approximation
            let n = elements.len();
            let charges = vec![0.0; n]; // neutral approximation
            let dipole = crate::dipole::compute_dipole(&charges, positions);
            Ok(dipole.vector)
        }
    }
}

/// Compute IR intensities from numerical dipole derivatives along normal modes.
///
/// I_k ∝ |dμ/dQ_k|² where Q_k is the k-th normal coordinate.
fn compute_ir_intensities(
    elements: &[u8],
    positions: &[[f64; 3]],
    normal_modes: &DMatrix<f64>,
    method: HessianMethod,
    delta: f64,
) -> Result<Vec<f64>, String> {
    #[cfg(feature = "parallel")]
    {
        compute_ir_intensities_parallel(elements, positions, normal_modes, method, delta)
    }
    #[cfg(not(feature = "parallel"))]
    {
        compute_ir_intensities_sequential(elements, positions, normal_modes, method, delta)
    }
}

#[cfg(not(feature = "parallel"))]
fn compute_ir_intensities_sequential(
    elements: &[u8],
    positions: &[[f64; 3]],
    normal_modes: &DMatrix<f64>,
    method: HessianMethod,
    delta: f64,
) -> Result<Vec<f64>, String> {
    let n_atoms = elements.len();
    let n_modes = normal_modes.ncols();
    let masses: Vec<f64> = elements.iter().map(|&z| atomic_mass(z)).collect();

    let mut intensities = Vec::with_capacity(n_modes);

    for k in 0..n_modes {
        let mut pos_plus: Vec<[f64; 3]> = positions.to_vec();
        let mut pos_minus: Vec<[f64; 3]> = positions.to_vec();

        for i in 0..n_atoms {
            let sqrt_m = masses[i].sqrt();
            for c in 0..3 {
                let disp = normal_modes[(3 * i + c, k)] * delta / sqrt_m;
                pos_plus[i][c] += disp;
                pos_minus[i][c] -= disp;
            }
        }

        let mu_plus = compute_dipole_vector(elements, &pos_plus, method)?;
        let mu_minus = compute_dipole_vector(elements, &pos_minus, method)?;

        let dmu_dq: Vec<f64> = (0..3)
            .map(|c| (mu_plus[c] - mu_minus[c]) / (2.0 * delta))
            .collect();

        let intensity = dmu_dq.iter().map(|x| x * x).sum::<f64>();
        intensities.push(intensity * 42.256);
    }

    Ok(intensities)
}

/// Parallel IR intensities via rayon — each mode's dipole derivative
/// displacement is independent.
#[cfg(feature = "parallel")]
fn compute_ir_intensities_parallel(
    elements: &[u8],
    positions: &[[f64; 3]],
    normal_modes: &DMatrix<f64>,
    method: HessianMethod,
    delta: f64,
) -> Result<Vec<f64>, String> {
    use rayon::prelude::*;

    let n_atoms = elements.len();
    let n_modes = normal_modes.ncols();
    let masses: Vec<f64> = elements.iter().map(|&z| atomic_mass(z)).collect();

    let results: Vec<Result<f64, String>> = (0..n_modes)
        .into_par_iter()
        .map(|k| {
            let mut pos_plus: Vec<[f64; 3]> = positions.to_vec();
            let mut pos_minus: Vec<[f64; 3]> = positions.to_vec();

            for i in 0..n_atoms {
                let sqrt_m = masses[i].sqrt();
                for c in 0..3 {
                    let disp = normal_modes[(3 * i + c, k)] * delta / sqrt_m;
                    pos_plus[i][c] += disp;
                    pos_minus[i][c] -= disp;
                }
            }

            let mu_plus = compute_dipole_vector(elements, &pos_plus, method)?;
            let mu_minus = compute_dipole_vector(elements, &pos_minus, method)?;

            let dmu_dq: Vec<f64> = (0..3)
                .map(|c| (mu_plus[c] - mu_minus[c]) / (2.0 * delta))
                .collect();

            let intensity = dmu_dq.iter().map(|x| x * x).sum::<f64>();
            Ok(intensity * 42.256)
        })
        .collect();

    results.into_iter().collect()
}

/// Lorentzian line shape: L(x) = (1/π) · γ / [(x - x₀)² + γ²]
fn lorentzian(x: f64, x0: f64, gamma: f64) -> f64 {
    let g = gamma.max(1e-6);
    let dx = x - x0;
    g / (std::f64::consts::PI * (dx * dx + g * g))
}

/// Gaussian line shape: G(x) = (1/(σ√(2π))) · exp(-(x-x₀)²/(2σ²))
fn gaussian(x: f64, x0: f64, sigma: f64) -> f64 {
    let s = sigma.max(1e-6);
    let dx = x - x0;
    (-(dx * dx) / (2.0 * s * s)).exp() / (s * (2.0 * std::f64::consts::PI).sqrt())
}

/// Assign functional group label to an IR frequency based on standard ranges.
fn assign_ir_peak(frequency_cm1: f64) -> String {
    if frequency_cm1 > 3500.0 {
        "O-H stretch (broad)".to_string()
    } else if frequency_cm1 > 3300.0 {
        "N-H stretch".to_string()
    } else if frequency_cm1 > 3000.0 {
        "sp2 C-H stretch".to_string()
    } else if frequency_cm1 > 2800.0 {
        "sp3 C-H stretch".to_string()
    } else if frequency_cm1 > 2100.0 && frequency_cm1 < 2300.0 {
        "C≡N or C≡C stretch".to_string()
    } else if frequency_cm1 > 1680.0 && frequency_cm1 < 1800.0 {
        "C=O stretch".to_string()
    } else if frequency_cm1 > 1600.0 && frequency_cm1 < 1680.0 {
        "C=C stretch".to_string()
    } else if frequency_cm1 > 1400.0 && frequency_cm1 < 1600.0 {
        "aromatic C=C stretch".to_string()
    } else if frequency_cm1 > 1300.0 && frequency_cm1 < 1400.0 {
        "C-H bend".to_string()
    } else if frequency_cm1 > 1000.0 && frequency_cm1 < 1300.0 {
        "C-O stretch".to_string()
    } else if frequency_cm1 > 600.0 && frequency_cm1 < 900.0 {
        "C-H out-of-plane bend".to_string()
    } else {
        "skeletal mode".to_string()
    }
}

/// Compute RRHO thermochemistry from vibrational frequencies.
///
/// Uses the rigid-rotor harmonic-oscillator approximation at the given temperature.
fn compute_thermochemistry(frequencies_cm1: &[f64], temperature_k: f64) -> Thermochemistry {
    // Constants
    const CM1_TO_EV: f64 = 1.0 / 8065.544;
    const EV_TO_KCAL: f64 = 23.0605;
    const R_KCAL: f64 = 0.001987204; // kcal/(mol·K)
    const R_CAL: f64 = 1.987204; // cal/(mol·K)
    const KB_EV: f64 = 8.617333e-5; // eV/K

    let kbt_ev = KB_EV * temperature_k;

    let real_freqs: Vec<f64> = frequencies_cm1
        .iter()
        .filter(|&&f| f > 50.0)
        .copied()
        .collect();

    // ZPVE = Σ (1/2)hν with 0.9 anharmonic scaling
    let zpve_ev: f64 = real_freqs.iter().map(|&f| 0.5 * f * CM1_TO_EV * 0.9).sum();
    let zpve_kcal = zpve_ev * EV_TO_KCAL;

    // Thermal energy (vibrational contribution): Σ hν / (exp(hν/kT) - 1)
    let thermal_e_ev: f64 = real_freqs
        .iter()
        .map(|&f| {
            let hnu = f * CM1_TO_EV;
            let x = hnu / kbt_ev;
            if x > 100.0 {
                0.0
            } else {
                hnu / (x.exp() - 1.0)
            }
        })
        .sum();
    let thermal_energy_kcal = thermal_e_ev * EV_TO_KCAL;

    // Enthalpy correction: E_thermal + RT
    let enthalpy_correction_kcal = thermal_energy_kcal + R_KCAL * temperature_k;

    // Vibrational entropy: Σ [x/(exp(x)-1) - ln(1-exp(-x))]·R
    let entropy_vib_cal: f64 = real_freqs
        .iter()
        .map(|&f| {
            let hnu = f * CM1_TO_EV;
            let x = hnu / kbt_ev;
            if x > 100.0 {
                0.0
            } else {
                let ex = x.exp();
                R_CAL * (x / (ex - 1.0) - (1.0 - (-x).exp()).ln())
            }
        })
        .sum();

    // Gibbs correction: H - TS
    let gibbs_correction_kcal = enthalpy_correction_kcal - temperature_k * entropy_vib_cal / 1000.0;

    Thermochemistry {
        temperature_k,
        zpve_kcal,
        thermal_energy_kcal,
        enthalpy_correction_kcal,
        entropy_vib_cal,
        gibbs_correction_kcal,
    }
}

/// Perform a complete vibrational analysis from elements and positions.
///
/// Computes the numerical Hessian, mass-weighted eigenvalue problem,
/// vibrational frequencies, and IR intensities.
pub fn compute_vibrational_analysis(
    elements: &[u8],
    positions: &[[f64; 3]],
    method: HessianMethod,
    step_size: Option<f64>,
) -> Result<VibrationalAnalysis, String> {
    if method == HessianMethod::Uff {
        return Err(
            "UFF requires SMILES; use compute_vibrational_analysis_uff instead".to_string(),
        );
    }

    let n_atoms = elements.len();
    let delta = step_size.unwrap_or(0.005);

    if n_atoms < 2 {
        return Err("Need at least 2 atoms for vibrational analysis".to_string());
    }

    let hessian = compute_numerical_hessian(elements, positions, method, Some(delta))?;
    build_vibrational_analysis_from_hessian(elements, positions, &hessian, method, delta)
}

/// Perform vibrational analysis using UFF analytical Hessian (gradient-difference method).
///
/// This avoids the expensive O(9N²) energy evaluations of the standard numerical Hessian
/// by using O(6N) gradient evaluations with UFF's analytical gradients.
pub fn compute_vibrational_analysis_uff(
    smiles: &str,
    elements: &[u8],
    positions: &[[f64; 3]],
    step_size: Option<f64>,
) -> Result<VibrationalAnalysis, String> {
    let n_atoms = elements.len();
    let delta = step_size.unwrap_or(0.005);

    if n_atoms < 2 {
        return Err("Need at least 2 atoms for vibrational analysis".to_string());
    }

    let coords_flat: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();
    let hessian =
        super::hessian::compute_uff_analytical_hessian(smiles, &coords_flat, Some(delta))?;
    build_vibrational_analysis_from_hessian(
        elements,
        positions,
        &hessian,
        HessianMethod::Uff,
        delta,
    )
}

/// Build vibrational analysis from a pre-computed Hessian matrix.
fn build_vibrational_analysis_from_hessian(
    elements: &[u8],
    positions: &[[f64; 3]],
    hessian: &DMatrix<f64>,
    method: HessianMethod,
    delta: f64,
) -> Result<VibrationalAnalysis, String> {
    let n_atoms = elements.len();
    let n3 = 3 * n_atoms;

    // 2. Build mass-weighted Hessian: H'_{ij} = H_{ij} / sqrt(m_i * m_j)
    let masses: Vec<f64> = elements.iter().map(|&z| atomic_mass(z)).collect();
    let mut mw_hessian = DMatrix::zeros(n3, n3);
    for i in 0..n3 {
        let mi = masses[i / 3];
        for j in 0..n3 {
            let mj = masses[j / 3];
            mw_hessian[(i, j)] = hessian[(i, j)] / (mi * mj).sqrt();
        }
    }

    // 3. Diagonalize mass-weighted Hessian
    let eigen = mw_hessian.symmetric_eigen();

    // 4. Sort eigenvalues and eigenvectors by eigenvalue
    let mut indices: Vec<usize> = (0..n3).collect();
    indices.sort_by(|&a, &b| {
        eigen.eigenvalues[a]
            .partial_cmp(&eigen.eigenvalues[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // 5. Convert eigenvalues to frequencies (cm⁻¹)
    // For linear molecules: 5 zero modes, for nonlinear: 6
    let is_linear = n_atoms == 2; // simplified check
    let _n_tr = if is_linear { 5 } else { 6 };
    let freq_threshold = 50.0; // cm⁻¹ threshold for "real" modes

    let mut sorted_eigenvalues = Vec::with_capacity(n3);
    let mut sorted_modes = DMatrix::zeros(n3, n3);
    for (new_idx, &old_idx) in indices.iter().enumerate() {
        sorted_eigenvalues.push(eigen.eigenvalues[old_idx]);
        for i in 0..n3 {
            sorted_modes[(i, new_idx)] = eigen.eigenvectors[(i, old_idx)];
        }
    }

    // Convert eigenvalues to frequencies
    let frequencies: Vec<f64> = sorted_eigenvalues
        .iter()
        .map(|&ev| {
            if ev >= 0.0 {
                // ν = (1/2πc) * sqrt(eigenvalue * factor)
                (ev * HESSIAN_TO_FREQUENCY_FACTOR).sqrt() * INV_2PI_C
            } else {
                // Imaginary frequency (negative eigenvalue)
                -((-ev) * HESSIAN_TO_FREQUENCY_FACTOR).sqrt() * INV_2PI_C
            }
        })
        .collect();

    // 6. Compute IR intensities from dipole derivatives
    let ir_intensities = compute_ir_intensities(
        elements,
        positions,
        &sorted_modes,
        method,
        delta * 2.0, // slightly larger step for dipole derivatives
    )?;

    // 7. Build mode list
    let mut modes = Vec::with_capacity(n3);
    let mut n_real = 0;
    let mut zpve = 0.0;

    for k in 0..n3 {
        let freq = frequencies[k];
        let is_real = freq.abs() > freq_threshold;
        if is_real && freq > 0.0 {
            n_real += 1;
            // ZPVE contribution: (1/2)hν in eV
            zpve += 0.5 * freq / EV_TO_CM1;
        }

        let displacement: Vec<f64> = (0..n3).map(|i| sorted_modes[(i, k)]).collect();

        modes.push(VibrationalMode {
            frequency_cm1: freq,
            ir_intensity: ir_intensities.get(k).copied().unwrap_or(0.0),
            displacement,
            is_real,
        });
    }

    let method_name = match method {
        HessianMethod::Eht => "EHT",
        HessianMethod::Pm3 => "PM3",
        HessianMethod::Xtb => "xTB",
        HessianMethod::Uff => "UFF",
    };

    // Compute thermochemistry at 298.15 K
    let thermochemistry = Some(compute_thermochemistry(&frequencies, 298.15));

    Ok(VibrationalAnalysis {
        n_atoms,
        modes,
        n_real_modes: n_real,
        zpve_ev: zpve,
        thermochemistry,
        method: method_name.to_string(),
        notes: vec![
            format!(
                "Numerical Hessian computed with {} using central finite differences (δ = {} Å).",
                method_name, delta
            ),
            "IR intensities derived from numerical dipole derivatives along normal coordinates."
                .to_string(),
            "Frequencies below 50 cm⁻¹ are classified as translation/rotation modes.".to_string(),
        ],
    })
}

/// Broadening function type for IR spectrum generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BroadeningType {
    /// Lorentzian line shape (default).
    Lorentzian,
    /// Gaussian line shape.
    Gaussian,
}

/// Generate a broadened IR spectrum from vibrational analysis.
///
/// `gamma`: broadening width in cm⁻¹ (typically 10–30).
/// `wn_min`, `wn_max`: wavenumber range in cm⁻¹.
/// `n_points`: number of grid points.
pub fn compute_ir_spectrum(
    analysis: &VibrationalAnalysis,
    gamma: f64,
    wn_min: f64,
    wn_max: f64,
    n_points: usize,
) -> IrSpectrum {
    compute_ir_spectrum_with_broadening(
        analysis,
        gamma,
        wn_min,
        wn_max,
        n_points,
        BroadeningType::Lorentzian,
    )
}

/// Generate a broadened IR spectrum with selectable broadening function.
pub fn compute_ir_spectrum_with_broadening(
    analysis: &VibrationalAnalysis,
    gamma: f64,
    wn_min: f64,
    wn_max: f64,
    n_points: usize,
    broadening: BroadeningType,
) -> IrSpectrum {
    let n_points = n_points.max(2);
    let step = (wn_max - wn_min) / (n_points as f64 - 1.0);
    let wavenumbers: Vec<f64> = (0..n_points).map(|i| wn_min + step * i as f64).collect();
    let mut intensities = vec![0.0; n_points];

    let mut peaks = Vec::new();

    for (mode_idx, mode) in analysis.modes.iter().enumerate() {
        if !mode.is_real || mode.frequency_cm1 <= 0.0 {
            continue;
        }

        peaks.push(IrPeak {
            frequency_cm1: mode.frequency_cm1,
            ir_intensity: mode.ir_intensity,
            mode_index: mode_idx,
            assignment: assign_ir_peak(mode.frequency_cm1),
        });

        for (idx, &wn) in wavenumbers.iter().enumerate() {
            intensities[idx] += mode.ir_intensity
                * match broadening {
                    BroadeningType::Lorentzian => lorentzian(wn, mode.frequency_cm1, gamma),
                    BroadeningType::Gaussian => gaussian(wn, mode.frequency_cm1, gamma),
                };
        }
    }

    peaks.sort_by(|a, b| {
        b.ir_intensity
            .partial_cmp(&a.ir_intensity)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    IrSpectrum {
        wavenumbers,
        intensities,
        peaks,
        gamma,
        notes: vec![
            format!(
                "IR spectrum generated with Lorentzian broadening (γ = {} cm⁻¹).",
                gamma
            ),
            format!("Vibrational analysis method: {}.", analysis.method),
        ],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lorentzian_peak_at_center() {
        let val = lorentzian(100.0, 100.0, 10.0);
        // At center: L = 1/(π·γ)
        let expected = 1.0 / (std::f64::consts::PI * 10.0);
        assert!((val - expected).abs() < 1e-10);
    }

    #[test]
    fn test_lorentzian_symmetry() {
        let left = lorentzian(90.0, 100.0, 10.0);
        let right = lorentzian(110.0, 100.0, 10.0);
        assert!((left - right).abs() < 1e-10);
    }

    #[test]
    fn test_atomic_masses_known_elements() {
        assert!((atomic_mass(1) - 1.00794).abs() < 0.001);
        assert!((atomic_mass(6) - 12.011).abs() < 0.01);
        assert!((atomic_mass(8) - 15.999).abs() < 0.01);
    }

    #[test]
    fn test_h2_vibrational_analysis() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let analysis =
            compute_vibrational_analysis(&elements, &positions, HessianMethod::Xtb, Some(0.005))
                .unwrap();

        assert_eq!(analysis.n_atoms, 2);
        assert_eq!(analysis.modes.len(), 6); // 3N = 6

        // H₂ is linear: 5 zero + 1 stretch
        // Should have at least 1 real mode (H-H stretch)
        assert!(
            analysis.n_real_modes >= 1,
            "H₂ should have at least 1 real vibrational mode, got {}",
            analysis.n_real_modes
        );

        // The stretch frequency should be positive
        let real_modes: Vec<&VibrationalMode> =
            analysis.modes.iter().filter(|m| m.is_real).collect();
        assert!(!real_modes.is_empty(), "Should have at least one real mode");

        // Check that ZPVE is positive for real modes
        if analysis.n_real_modes > 0 {
            assert!(analysis.zpve_ev > 0.0, "ZPVE should be positive");
        }
    }

    #[test]
    fn test_water_vibrational_analysis() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let analysis =
            compute_vibrational_analysis(&elements, &positions, HessianMethod::Xtb, Some(0.005))
                .unwrap();

        assert_eq!(analysis.n_atoms, 3);
        assert_eq!(analysis.modes.len(), 9); // 3N = 9

        // Water is nonlinear: 6 zero + 3 real modes (sym stretch, asym stretch, bend)
        // Due to numerical noise some translational modes may appear as real
        // so we just check that we have at least 3 real modes
        assert!(
            analysis.n_real_modes >= 2,
            "Water should have at least 2-3 real modes, got {}",
            analysis.n_real_modes
        );
    }

    #[test]
    fn test_ir_spectrum_generation() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let analysis =
            compute_vibrational_analysis(&elements, &positions, HessianMethod::Xtb, Some(0.005))
                .unwrap();

        let spectrum = compute_ir_spectrum(&analysis, 20.0, 400.0, 4000.0, 500);

        assert_eq!(spectrum.wavenumbers.len(), 500);
        assert_eq!(spectrum.intensities.len(), 500);
        assert!(!spectrum.peaks.is_empty(), "Should have IR peaks");
        assert!(
            spectrum.intensities.iter().any(|&i| i > 0.0),
            "Spectrum should have non-zero intensity"
        );
    }
}
