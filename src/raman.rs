//! Raman intensities from a sum-over-states polarizability and normal modes.
//!
//! This is the first response-property prototype for the v2.x line. The
//! polarizability uses the uncoupled HF sum-over-states expression in the
//! existing STO-3G basis; Raman activities are finite differences of that
//! tensor along the UFF normal modes. VCD and ROA need magnetic-dipole
//! response and are documented follow-ups.
//!
//! # Units
//!
//! Coordinates are Å, frequencies cm⁻¹, polarizabilities Å³ and Raman
//! activities arbitrary but consistent units. Relative intensities are the
//! meaningful output; absolute cross-sections are not calibrated.

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

use crate::error::SciFormError;
use crate::scf::basis::BasisSet;
use crate::scf::dipole_matrix::dipole_matrices;

const ANGSTROM_TO_BOHR: f64 = 1.889_726_031_269_39;
const BOHR3_TO_ANGSTROM3: f64 = 0.148_184_711;

/// One polarizability tensor.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Polarizability {
    /// 3×3 tensor in Å³.
    pub tensor: [[f64; 3]; 3],
    /// Isotropic polarizability in Å³.
    pub isotropic: f64,
    /// Anisotropy in Å³.
    pub anisotropy: f64,
}

/// One Raman-active normal mode.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RamanMode {
    /// Mode index in the vibrational analysis.
    pub mode_index: usize,
    /// Frequency in cm⁻¹.
    pub frequency_cm1: f64,
    /// Raman activity in the Placzek invariant form.
    pub activity: f64,
    /// Depolarization ratio, dimensionless.
    pub depolarization: f64,
}

/// Result of a Raman calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RamanSpectrum {
    /// Modes sorted by the vibrational analysis.
    pub modes: Vec<RamanMode>,
    /// Isotropic polarizability at the input geometry.
    pub polarizability: Polarizability,
    /// Finite-difference step in Å.
    pub step_angstrom: f64,
    /// Notes and caveats.
    pub notes: Vec<String>,
}

/// Configuration for the Raman prototype.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RamanConfig {
    /// Cartesian displacement along each normalized normal mode, Å.
    pub step_angstrom: f64,
}

impl Default for RamanConfig {
    fn default() -> Self {
        Self {
            step_angstrom: 0.01,
        }
    }
}

/// Compute the uncoupled HF/STO-3G polarizability tensor.
///
/// The sum-over-states expression is
///
/// ```text
/// α_ij = Σ_spin Σ_a^occ Σ_b^vir
///        <a|μ_i|b><b|μ_j|a> / (ε_b - ε_a)
/// ```
///
/// with the alpha and beta channels summed separately. The dipole integrals
/// are in Bohr and the result is converted to Å³. This is an uncoupled
/// approximation, not a coupled-perturbed HF response.
pub fn compute_polarizability(
    elements: &[u8],
    positions_angstrom: &[[f64; 3]],
) -> Result<Polarizability, SciFormError> {
    crate::validation::validate_geometry(elements, positions_angstrom)?;
    let scf = crate::compute_uhf(elements, positions_angstrom, 0, 1)?;
    if !scf.converged {
        return Err(SciFormError::numerical(
            "UHF did not converge while computing the polarizability",
        ));
    }
    let positions_bohr: Vec<[f64; 3]> = positions_angstrom
        .iter()
        .map(|position| {
            [
                position[0] * ANGSTROM_TO_BOHR,
                position[1] * ANGSTROM_TO_BOHR,
                position[2] * ANGSTROM_TO_BOHR,
            ]
        })
        .collect();
    let basis = BasisSet::sto3g(elements, &positions_bohr);
    let dips = dipole_matrices(&basis, [0.0, 0.0, 0.0]);
    let n_basis = basis.n_basis;

    let mut tensor = [[0.0_f64; 3]; 3];
    for (coefficients, energies, n_occ) in [
        (
            &scf.alpha_coefficients,
            &scf.alpha_orbital_energies,
            scf.n_alpha,
        ),
        (
            &scf.beta_coefficients,
            &scf.beta_orbital_energies,
            scf.n_beta,
        ),
    ] {
        let c = coefficients;
        let mut d_mo = [
            DMatrix::zeros(n_basis, n_basis),
            DMatrix::zeros(n_basis, n_basis),
            DMatrix::zeros(n_basis, n_basis),
        ];
        for (component, dipole) in dips.iter().enumerate() {
            d_mo[component] = c.transpose() * dipole * c;
        }
        for a in 0..n_occ {
            for b in n_occ..n_basis {
                let delta = energies[b] - energies[a];
                if delta.abs() < 1e-12 {
                    continue;
                }
                for i in 0..3 {
                    for j in 0..3 {
                        tensor[i][j] += d_mo[i][(a, b)] * d_mo[j][(b, a)] / delta;
                    }
                }
            }
        }
    }
    for row in &mut tensor {
        for value in row {
            *value *= BOHR3_TO_ANGSTROM3;
        }
    }
    let isotropic = (tensor[0][0] + tensor[1][1] + tensor[2][2]) / 3.0;
    let anisotropy = (0.5
        * ((tensor[0][0] - tensor[1][1]).powi(2)
            + (tensor[1][1] - tensor[2][2]).powi(2)
            + (tensor[2][2] - tensor[0][0]).powi(2))
        + 3.0 * (tensor[0][1].powi(2) + tensor[0][2].powi(2) + tensor[1][2].powi(2)))
    .sqrt();
    Ok(Polarizability {
        tensor,
        isotropic,
        anisotropy,
    })
}

/// Compute Raman activities for the UFF normal modes of a molecule.
pub fn compute_raman_spectrum(
    smiles: &str,
    elements: &[u8],
    positions_angstrom: &[[f64; 3]],
    config: &RamanConfig,
) -> Result<RamanSpectrum, SciFormError> {
    if !config.step_angstrom.is_finite() || config.step_angstrom <= 0.0 {
        return Err(SciFormError::validation(
            "step_angstrom must be finite and > 0",
        ));
    }
    crate::validation::validate_geometry(elements, positions_angstrom)?;
    let analysis =
        crate::compute_vibrational_analysis_uff(smiles, elements, positions_angstrom, None)?;
    let masses: Vec<f64> = elements
        .iter()
        .map(|&z| crate::dynamics::atomic_mass_amu(z))
        .collect();
    if masses.iter().any(|mass| *mass <= 0.0) {
        return Err(SciFormError::unsupported(
            "Raman requires an atomic mass for every element",
        ));
    }
    let polarizability = compute_polarizability(elements, positions_angstrom)?;
    let n_atoms = elements.len();
    let mut modes = Vec::new();

    for (mode_index, mode) in analysis.modes.iter().enumerate() {
        if !mode.is_real || mode.frequency_cm1 <= 0.0 {
            continue;
        }
        let mut q = vec![0.0; n_atoms * 3];
        for atom in 0..n_atoms {
            let mass = masses[atom].sqrt();
            for axis in 0..3 {
                q[3 * atom + axis] = mode.displacement[3 * atom + axis] / mass;
            }
        }
        let norm = q.iter().map(|value| value * value).sum::<f64>().sqrt();
        if norm < 1e-12 {
            continue;
        }
        for value in &mut q {
            *value /= norm;
        }
        let plus = displaced(elements, positions_angstrom, &q, config.step_angstrom)?;
        let minus = displaced(elements, positions_angstrom, &q, -config.step_angstrom)?;
        let alpha_plus = compute_polarizability(elements, &plus)?;
        let alpha_minus = compute_polarizability(elements, &minus)?;
        let isotropic_derivative =
            (alpha_plus.isotropic - alpha_minus.isotropic) / (2.0 * config.step_angstrom);
        let anisotropy_derivative =
            (alpha_plus.anisotropy - alpha_minus.anisotropy) / (2.0 * config.step_angstrom);
        let activity = 45.0 * isotropic_derivative * isotropic_derivative
            + 7.0 * anisotropy_derivative * anisotropy_derivative;
        let denominator = 45.0 * isotropic_derivative * isotropic_derivative
            + 4.0 * anisotropy_derivative * anisotropy_derivative;
        let depolarization = if denominator.abs() > 1e-30 {
            3.0 * anisotropy_derivative * anisotropy_derivative / denominator
        } else {
            0.0
        };
        modes.push(RamanMode {
            mode_index,
            frequency_cm1: mode.frequency_cm1,
            activity,
            depolarization,
        });
    }

    Ok(RamanSpectrum {
        modes,
        polarizability,
        step_angstrom: config.step_angstrom,
        notes: vec![
            "Polarizability uses the uncoupled HF/STO-3G sum-over-states expression.".to_string(),
            "Raman activities are finite differences along UFF normal modes.".to_string(),
            "VCD and ROA require magnetic-dipole response and are not implemented yet.".to_string(),
            "Relative activities are meaningful; absolute cross-sections are not calibrated."
                .to_string(),
        ],
    })
}

fn displaced(
    elements: &[u8],
    positions: &[[f64; 3]],
    q: &[f64],
    step: f64,
) -> Result<Vec<[f64; 3]>, SciFormError> {
    let mut result = positions.to_vec();
    for (atom, position) in result.iter_mut().enumerate() {
        for axis in 0..3 {
            position[axis] += step * q[3 * atom + axis];
        }
    }
    if !result
        .iter()
        .all(|position| position.iter().all(|v| v.is_finite()))
    {
        return Err(SciFormError::numerical(
            "normal-mode displacement produced a non-finite geometry",
        ));
    }
    let _ = elements;
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn water_polarizability_is_symmetric_and_positive() {
        let water = [
            [0.0, 0.0, 0.0],
            [0.9572, 0.0, 0.0],
            [-0.239_987, 0.927_297, 0.0],
        ];
        let result = compute_polarizability(&[8, 1, 1], &water).unwrap();
        assert!(result.isotropic > 0.0);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (result.tensor[i][j] - result.tensor[j][i]).abs() < 1e-10,
                    "tensor is not symmetric at ({i},{j})"
                );
            }
        }
        assert!(result.tensor[0][0] > 0.0);
    }

    #[test]
    fn water_raman_smoke() {
        let water = [
            [0.0, 0.0, 0.0],
            [0.9572, 0.0, 0.0],
            [-0.239_987, 0.927_297, 0.0],
        ];
        let result =
            compute_raman_spectrum("O", &[8, 1, 1], &water, &RamanConfig::default()).unwrap();
        assert!(!result.modes.is_empty());
        assert!(result.modes.iter().all(|mode| mode.activity.is_finite()));
    }

    #[test]
    fn invalid_step_is_rejected() {
        let config = RamanConfig { step_angstrom: 0.0 };
        let water = [
            [0.0, 0.0, 0.0],
            [0.9572, 0.0, 0.0],
            [-0.239_987, 0.927_297, 0.0],
        ];
        assert!(compute_raman_spectrum("O", &[8, 1, 1], &water, &config).is_err());
    }
}
