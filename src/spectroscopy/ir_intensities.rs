//! IR intensities from dipole derivatives.
//!
//! IR absorption intensity for normal mode k is proportional to:
//!
//!   I_k = |∂μ/∂Q_k|²
//!
//! where μ is the molecular dipole moment and Q_k is the normal coordinate.

use crate::scf::scf_loop::{run_scf, ScfConfig};
use crate::scf::types::{IrResult, MolecularSystem, ScfResult};

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
    /// Assignment label.
    pub assignment: String,
}

/// Compute IR spectrum from Hessian results and SCF data.
pub fn compute_ir_intensities(
    system: &MolecularSystem,
    hessian_result: &HessianResult,
    config: &IrConfig,
) -> IrResult {
    let n_atoms = system.n_atoms();
    let n_coords = n_atoms * 3;

    let n_tr = if is_linear(system) { 5 } else { 6 };

    let masses: Vec<f64> = system
        .atomic_numbers
        .iter()
        .flat_map(|&z| {
            let m = atomic_mass(z);
            vec![m, m, m]
        })
        .collect();

    let mut peaks = Vec::new();

    for mode_idx in n_tr..n_coords {
        let freq = hessian_result.frequencies[mode_idx];

        if freq.abs() < config.freq_min || freq > config.freq_max {
            continue;
        }

        let q_mw: Vec<f64> = (0..n_coords)
            .map(|i| hessian_result.normal_modes[(i, mode_idx)])
            .collect();

        let q_cart: Vec<f64> = q_mw
            .iter()
            .enumerate()
            .map(|(i, &q)| q / masses[i].sqrt())
            .collect();

        let mut sys_plus = system.clone();
        let mut sys_minus = system.clone();

        for atom in 0..n_atoms {
            for coord in 0..3 {
                let idx = atom * 3 + coord;
                sys_plus.positions_bohr[atom][coord] += config.dipole_step * q_cart[idx];
                sys_minus.positions_bohr[atom][coord] -= config.dipole_step * q_cart[idx];
            }
        }

        let scf_plus = run_scf(&sys_plus, &config.scf_config);
        let scf_minus = run_scf(&sys_minus, &config.scf_config);

        let dipole_plus = compute_dipole_from_scf(&scf_plus, &sys_plus);
        let dipole_minus = compute_dipole_from_scf(&scf_minus, &sys_minus);

        let dmu_dq: [f64; 3] = [
            (dipole_plus[0] - dipole_minus[0]) / (2.0 * config.dipole_step),
            (dipole_plus[1] - dipole_minus[1]) / (2.0 * config.dipole_step),
            (dipole_plus[2] - dipole_minus[2]) / (2.0 * config.dipole_step),
        ];

        let intensity_au = dmu_dq[0].powi(2) + dmu_dq[1].powi(2) + dmu_dq[2].powi(2);
        let au_to_km_mol = 42.256;
        let intensity = intensity_au * au_to_km_mol;

        peaks.push(IrPeak {
            frequency: freq,
            intensity,
            mode_index: mode_idx,
            assignment: String::new(),
        });
    }

    peaks.sort_by(|a, b| a.frequency.partial_cmp(&b.frequency).unwrap());

    let frequencies: Vec<f64> = peaks.iter().map(|p| p.frequency).collect();
    let intensities: Vec<f64> = peaks.iter().map(|p| p.intensity).collect();

    IrResult {
        frequencies,
        intensities,
        n_modes: peaks.len(),
    }
}

/// Compute electric dipole moment from SCF density matrix (Mulliken approximation).
fn compute_dipole_from_scf(scf: &ScfResult, system: &MolecularSystem) -> [f64; 3] {
    let n_atoms = system.n_atoms();
    let mut dipole = [0.0; 3];

    for a in 0..n_atoms {
        let z = system.atomic_numbers[a] as f64;
        for k in 0..3 {
            dipole[k] += z * system.positions_bohr[a][k];
        }
    }

    for a in 0..n_atoms {
        let q = scf.mulliken_charges[a];
        for k in 0..3 {
            dipole[k] -= (system.atomic_numbers[a] as f64 - q) * system.positions_bohr[a][k];
        }
    }

    dipole
}

/// Simple linearity check.
fn is_linear(system: &MolecularSystem) -> bool {
    if system.n_atoms() <= 2 {
        return true;
    }

    let p0 = system.positions_bohr[0];
    let p1 = system.positions_bohr[1];

    let v01 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let v01_len = (v01[0] * v01[0] + v01[1] * v01[1] + v01[2] * v01[2]).sqrt();

    if v01_len < 1e-10 {
        return false;
    }

    for i in 2..system.n_atoms() {
        let pi = system.positions_bohr[i];
        let v0i = [pi[0] - p0[0], pi[1] - p0[1], pi[2] - p0[2]];

        let cross = [
            v01[1] * v0i[2] - v01[2] * v0i[1],
            v01[2] * v0i[0] - v01[0] * v0i[2],
            v01[0] * v0i[1] - v01[1] * v0i[0],
        ];

        let cross_len = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
        let v0i_len = (v0i[0] * v0i[0] + v0i[1] * v0i[1] + v0i[2] * v0i[2]).sqrt();
        if v0i_len > 1e-10 && cross_len / (v01_len * v0i_len) > 1e-3 {
            return false;
        }
    }

    true
}
