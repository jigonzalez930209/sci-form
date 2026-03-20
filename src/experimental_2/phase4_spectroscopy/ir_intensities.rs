//! IR intensities from dipole derivatives.
//!
//! IR absorption intensity for normal mode k is proportional to:
//!
//!   I_k = |∂μ/∂Q_k|²
//!
//! where μ is the molecular dipole moment and Q_k is the normal coordinate.
//!
//! The dipole derivative is computed numerically:
//!   ∂μ/∂Q_k ≈ [μ(Q_k + δ) - μ(Q_k - δ)] / (2δ)
//!
//! Units: km/mol (standard IR intensity units).

use crate::experimental_2::types::{IrResult, MolecularSystem};
use crate::experimental_2::phase3_scf_engine::scf_loop::{run_scf, ScfConfig};

use super::hessian::{atomic_mass, HessianResult};

/// Configuration for IR intensity calculation.
#[derive(Debug, Clone)]
pub struct IrConfig {
    /// Displacement step for numerical dipole derivative (Bohr).
    pub dipole_step: f64,
    /// SCF configuration for energy evaluations.
    pub scf_config: ScfConfig,
    /// Broadening width for spectral convolution (cm⁻¹).
    pub broadening: f64,
    /// Minimum frequency to include (cm⁻¹).
    pub freq_min: f64,
    /// Maximum frequency to include (cm⁻¹).
    pub freq_max: f64,
}

impl Default for IrConfig {
    fn default() -> Self {
        Self {
            dipole_step: 1e-3,
            scf_config: ScfConfig::default(),
            broadening: 20.0,
            freq_min: 0.0,
            freq_max: 4000.0,
        }
    }
}

/// IR peak data.
#[derive(Debug, Clone)]
pub struct IrPeak {
    /// Frequency in cm⁻¹.
    pub frequency: f64,
    /// Intensity in km/mol.
    pub intensity: f64,
    /// Normal mode index.
    pub mode_index: usize,
    /// Assignment label (empty if not assigned).
    pub assignment: String,
}

/// Compute IR spectrum from Hessian results and SCF data.
///
/// Steps:
/// 1. For each normal mode Q_k, displace geometry along ±Q_k
/// 2. Compute dipole moment at displaced geometries
/// 3. Numerical derivative gives ∂μ/∂Q_k
/// 4. Intensity = |∂μ/∂Q_k|²
pub fn compute_ir_intensities(
    system: &MolecularSystem,
    hessian_result: &HessianResult,
    config: &IrConfig,
) -> IrResult {
    let n_atoms = system.n_atoms();
    let n_coords = n_atoms * 3;

    // Number of vibrational modes (3N - 6 for nonlinear, 3N - 5 for linear)
    let n_tr = if is_linear(system) { 5 } else { 6 };
    let _n_vib = if n_coords > n_tr { n_coords - n_tr } else { 0 };

    let masses: Vec<f64> = system
        .atomic_numbers
        .iter()
        .flat_map(|&z| {
            let m = atomic_mass(z);
            vec![m, m, m]
        })
        .collect();

    let mut peaks = Vec::new();

    // Skip first n_tr modes (translations + rotations)
    for mode_idx in n_tr..n_coords {
        let freq = hessian_result.frequencies[mode_idx];

        // Skip imaginary or near-zero frequencies
        if freq.abs() < config.freq_min || freq > config.freq_max {
            continue;
        }

        // Get normal mode displacement vector (in mass-weighted coordinates)
        let q_mw: Vec<f64> = (0..n_coords)
            .map(|i| hessian_result.normal_modes[(i, mode_idx)])
            .collect();

        // Convert to Cartesian displacements: δR_i = Q_i / sqrt(m_i)
        let q_cart: Vec<f64> = q_mw
            .iter()
            .enumerate()
            .map(|(i, &q)| q / masses[i].sqrt())
            .collect();

        // Displace along +Q and -Q
        let mut sys_plus = system.clone();
        let mut sys_minus = system.clone();

        for atom in 0..n_atoms {
            for coord in 0..3 {
                let idx = atom * 3 + coord;
                sys_plus.positions_bohr[atom][coord] += config.dipole_step * q_cart[idx];
                sys_minus.positions_bohr[atom][coord] -= config.dipole_step * q_cart[idx];
            }
        }

        // Compute dipole moments at displaced geometries
        let scf_plus = run_scf(&sys_plus, &config.scf_config);
        let scf_minus = run_scf(&sys_minus, &config.scf_config);

        let dipole_plus = compute_dipole_from_scf(&scf_plus, &sys_plus);
        let dipole_minus = compute_dipole_from_scf(&scf_minus, &sys_minus);

        // Numerical derivative: ∂μ/∂Q = [μ(+δ) - μ(-δ)] / (2δ)
        let dmu_dq: [f64; 3] = [
            (dipole_plus[0] - dipole_minus[0]) / (2.0 * config.dipole_step),
            (dipole_plus[1] - dipole_minus[1]) / (2.0 * config.dipole_step),
            (dipole_plus[2] - dipole_minus[2]) / (2.0 * config.dipole_step),
        ];

        // Intensity = |∂μ/∂Q|² (convert to km/mol)
        let intensity_au = dmu_dq[0].powi(2) + dmu_dq[1].powi(2) + dmu_dq[2].powi(2);
        // Conversion: 1 (e·a₀)²/(amu·a₀²) → km/mol
        let au_to_km_mol = 42.256; // approximate conversion
        let intensity = intensity_au * au_to_km_mol;

        peaks.push(IrPeak {
            frequency: freq,
            intensity,
            mode_index: mode_idx,
            assignment: String::new(),
        });
    }

    // Sort by frequency
    peaks.sort_by(|a, b| a.frequency.partial_cmp(&b.frequency).unwrap());

    let frequencies: Vec<f64> = peaks.iter().map(|p| p.frequency).collect();
    let intensities: Vec<f64> = peaks.iter().map(|p| p.intensity).collect();

    IrResult {
        frequencies,
        intensities,
        n_modes: peaks.len(),
    }
}

/// Compute electric dipole moment from SCF density matrix.
///
///   μ_k = -Σ_μν P_μν <μ|r_k|ν> + Σ_A Z_A R_{A,k}
///
/// Simplified: uses Mulliken charges for quick dipole estimate.
fn compute_dipole_from_scf(
    scf: &crate::experimental_2::types::ScfResult,
    system: &MolecularSystem,
) -> [f64; 3] {
    let n_atoms = system.n_atoms();
    let mut dipole = [0.0; 3];

    // Nuclear contribution
    for a in 0..n_atoms {
        let z = system.atomic_numbers[a] as f64;
        for k in 0..3 {
            dipole[k] += z * system.positions_bohr[a][k];
        }
    }

    // Electronic contribution (using Mulliken charges as approximation)
    for a in 0..n_atoms {
        let q = scf.mulliken_charges[a]; // partial charge on atom A
        for k in 0..3 {
            dipole[k] -= (system.atomic_numbers[a] as f64 - q) * system.positions_bohr[a][k];
        }
    }

    dipole
}

/// Simple linearity check for a molecule.
fn is_linear(system: &MolecularSystem) -> bool {
    if system.n_atoms() <= 2 {
        return true;
    }

    // Check if all atoms are collinear
    let p0 = system.positions_bohr[0];
    let p1 = system.positions_bohr[1];

    let v01 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let v01_len = (v01[0] * v01[0] + v01[1] * v01[1] + v01[2] * v01[2]).sqrt();

    if v01_len < 1e-10 {
        return false; // Atoms are on top of each other
    }

    for i in 2..system.n_atoms() {
        let pi = system.positions_bohr[i];
        let v0i = [pi[0] - p0[0], pi[1] - p0[1], pi[2] - p0[2]];

        // Cross product v01 × v0i
        let cross = [
            v01[1] * v0i[2] - v01[2] * v0i[1],
            v01[2] * v0i[0] - v01[0] * v0i[2],
            v01[0] * v0i[1] - v01[1] * v0i[0],
        ];

        let cross_len = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();

        // If cross product magnitude relative to vectors is > threshold, not linear
        let v0i_len = (v0i[0] * v0i[0] + v0i[1] * v0i[1] + v0i[2] * v0i[2]).sqrt();
        if v0i_len > 1e-10 && cross_len / (v01_len * v0i_len) > 1e-3 {
            return false;
        }
    }

    true
}

/// Broaden IR spectrum with Lorentzian line shape.
///
///   L(ω; ω₀, γ) = γ / [π · ((ω - ω₀)² + γ²)]
pub fn broaden_spectrum(
    peaks: &[IrPeak],
    gamma: f64,
    wn_min: f64,
    wn_max: f64,
    n_points: usize,
) -> (Vec<f64>, Vec<f64>) {
    let step = (wn_max - wn_min) / (n_points - 1) as f64;

    let wavenumbers: Vec<f64> = (0..n_points)
        .map(|i| wn_min + i as f64 * step)
        .collect();

    let spectrum: Vec<f64> = wavenumbers
        .iter()
        .map(|&wn| {
            peaks
                .iter()
                .map(|peak| {
                    let dw = wn - peak.frequency;
                    peak.intensity * gamma / (std::f64::consts::PI * (dw * dw + gamma * gamma))
                })
                .sum()
        })
        .collect();

    (wavenumbers, spectrum)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_linear_diatomic() {
        let system = MolecularSystem {
            atomic_numbers: vec![1, 1],
            positions_bohr: vec![[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]],
            charge: 0,
            multiplicity: 1,
        };
        assert!(is_linear(&system));
    }

    #[test]
    fn test_is_linear_triatomic_linear() {
        let system = MolecularSystem {
            atomic_numbers: vec![8, 6, 8],
            positions_bohr: vec![
                [-2.2, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [2.2, 0.0, 0.0],
            ],
            charge: 0,
            multiplicity: 1,
        };
        assert!(is_linear(&system));
    }

    #[test]
    fn test_is_not_linear_water() {
        let system = MolecularSystem {
            atomic_numbers: vec![8, 1, 1],
            positions_bohr: vec![
                [0.0, 0.0, 0.0],
                [1.43, 1.11, 0.0],
                [-1.43, 1.11, 0.0],
            ],
            charge: 0,
            multiplicity: 1,
        };
        assert!(!is_linear(&system));
    }

    #[test]
    fn test_lorentzian_broadening() {
        let peaks = vec![IrPeak {
            frequency: 1000.0,
            intensity: 100.0,
            mode_index: 0,
            assignment: String::new(),
        }];

        let (wn, spectrum) = broaden_spectrum(&peaks, 10.0, 500.0, 1500.0, 1001);

        // Maximum should be at or near 1000 cm⁻¹
        let max_idx = spectrum
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap()
            .0;

        assert!((wn[max_idx] - 1000.0).abs() < 2.0, "Peak should be near 1000 cm⁻¹");
    }
}
