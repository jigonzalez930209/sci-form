//! Vibrational circular dichroism (VCD) and Raman optical activity (ROA)
//! screening models.
//!
//! The first implementation uses a **nuclear charge-current** approximation for
//! the magnetic dipole derivative and the existing UFF normal modes. VCD
//! rotational strengths are `μ'·m`; ROA circular intensities weight the Raman
//! activity by the same magnetic term. This is a qualitative screening model,
//! not a coupled-perturbed magnetic response.
//!
//! # Units
//!
//! Frequencies are cm⁻¹, charges are elementary charge, positions Å. VCD
//! rotational strengths and ROA circular intensities are arbitrary but
//! consistent units; relative signs and magnitudes are the meaningful output.

use serde::{Deserialize, Serialize};

use crate::error::SciFormError;
use crate::raman::{RamanConfig, RamanSpectrum};

/// One VCD-active mode.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VcdMode {
    /// Mode index in the vibrational analysis.
    pub mode_index: usize,
    /// Frequency in cm⁻¹.
    pub frequency_cm1: f64,
    /// Electric-dipole derivative vector.
    pub electric_dipole_derivative: [f64; 3],
    /// Magnetic-dipole derivative vector.
    pub magnetic_dipole_derivative: [f64; 3],
    /// Rotational strength `μ'·m` in the screening model.
    pub rotational_strength: f64,
}

/// Result of a VCD calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VcdSpectrum {
    /// VCD modes.
    pub modes: Vec<VcdMode>,
    /// Notes and caveats.
    pub notes: Vec<String>,
}

/// One ROA-active mode.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RoaMode {
    /// Mode index in the vibrational analysis.
    pub mode_index: usize,
    /// Frequency in cm⁻¹.
    pub frequency_cm1: f64,
    /// Raman activity from the existing Raman prototype.
    pub raman_activity: f64,
    /// Screening circular intensity difference.
    pub circular_intensity: f64,
}

/// Result of an ROA calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RoaSpectrum {
    /// ROA modes.
    pub modes: Vec<RoaMode>,
    /// Raman spectrum used for the activities.
    pub raman: RamanSpectrum,
    /// Notes and caveats.
    pub notes: Vec<String>,
}

/// Configuration for the VCD/ROA prototypes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VcdRoaConfig {
    /// Raman configuration used for the ROA activities.
    pub raman: RamanConfig,
    /// Keep only modes above this frequency in cm⁻¹.
    pub min_frequency_cm1: f64,
}

impl Default for VcdRoaConfig {
    fn default() -> Self {
        Self {
            raman: RamanConfig::default(),
            min_frequency_cm1: 1.0,
        }
    }
}

/// Compute VCD rotational strengths for the UFF normal modes.
pub fn compute_vcd_spectrum(
    smiles: &str,
    elements: &[u8],
    positions_angstrom: &[[f64; 3]],
    config: &VcdRoaConfig,
) -> Result<VcdSpectrum, SciFormError> {
    crate::validation::validate_geometry(elements, positions_angstrom)?;
    if !config.min_frequency_cm1.is_finite() || config.min_frequency_cm1 < 0.0 {
        return Err(SciFormError::validation(
            "min_frequency_cm1 must be finite and >= 0",
        ));
    }
    let analysis =
        crate::compute_vibrational_analysis_uff(smiles, elements, positions_angstrom, None)?;
    let charges = crate::compute_charges(smiles)?.charges;
    if charges.len() != elements.len() {
        return Err(SciFormError::numerical(
            "charge count does not match the atom count",
        ));
    }
    let masses: Vec<f64> = elements
        .iter()
        .map(|&z| crate::dynamics::atomic_mass_amu(z))
        .collect();
    let mut modes = Vec::new();

    for (mode_index, mode) in analysis.modes.iter().enumerate() {
        if !mode.is_real
            || mode.frequency_cm1 <= config.min_frequency_cm1
            || mode.displacement.len() != elements.len() * 3
        {
            continue;
        }
        let displacement = normalized_cartesian_displacement(&mode.displacement, &masses);
        let mut electric = [0.0; 3];
        let mut magnetic = [0.0; 3];
        for atom in 0..elements.len() {
            let e = [
                displacement[3 * atom],
                displacement[3 * atom + 1],
                displacement[3 * atom + 2],
            ];
            let position = positions_angstrom[atom];
            let charge = charges[atom];
            for axis in 0..3 {
                electric[axis] += charge * e[axis];
            }
            let cross = [
                position[1] * e[2] - position[2] * e[1],
                position[2] * e[0] - position[0] * e[2],
                position[0] * e[1] - position[1] * e[0],
            ];
            for axis in 0..3 {
                magnetic[axis] += 0.5 * charge * cross[axis];
            }
        }
        let rotational_strength = dot(electric, magnetic);
        modes.push(VcdMode {
            mode_index,
            frequency_cm1: mode.frequency_cm1,
            electric_dipole_derivative: electric,
            magnetic_dipole_derivative: magnetic,
            rotational_strength,
        });
    }

    Ok(VcdSpectrum {
        modes,
        notes: vec![
            "Nuclear charge-current approximation for the magnetic dipole derivative.".to_string(),
            "No coupled-perturbed electronic magnetic response is included.".to_string(),
            "Relative signs and magnitudes are qualitative screening values.".to_string(),
        ],
    })
}

/// Compute ROA screening intensities from the Raman activities and the VCD
/// magnetic term.
pub fn compute_roa_spectrum(
    smiles: &str,
    elements: &[u8],
    positions_angstrom: &[[f64; 3]],
    config: &VcdRoaConfig,
) -> Result<RoaSpectrum, SciFormError> {
    let vcd = compute_vcd_spectrum(smiles, elements, positions_angstrom, config)?;
    let raman = crate::compute_raman_spectrum(smiles, elements, positions_angstrom, &config.raman)?;
    let by_index: std::collections::HashMap<usize, f64> = raman
        .modes
        .iter()
        .map(|mode| (mode.mode_index, mode.activity))
        .collect();
    let modes = vcd
        .modes
        .iter()
        .filter_map(|mode| {
            let activity = *by_index.get(&mode.mode_index)?;
            Some(RoaMode {
                mode_index: mode.mode_index,
                frequency_cm1: mode.frequency_cm1,
                raman_activity: activity,
                circular_intensity: mode.rotational_strength * activity,
            })
        })
        .collect();

    Ok(RoaSpectrum {
        modes,
        raman,
        notes: vec![
            "ROA circular intensity is a screening product of the Raman activity and the nuclear magnetic term."
                .to_string(),
            "This is not a full electric-dipole/magnetic-dipole response calculation.".to_string(),
            "Use relative signs and trends, not absolute ROA cross-sections.".to_string(),
        ],
    })
}

fn normalized_cartesian_displacement(displacement: &[f64], masses: &[f64]) -> Vec<f64> {
    let mut q = vec![0.0; displacement.len()];
    for atom in 0..masses.len() {
        let scale = masses[atom].sqrt();
        for axis in 0..3 {
            q[3 * atom + axis] = displacement[3 * atom + axis] / scale;
        }
    }
    let norm = q.iter().map(|value| value * value).sum::<f64>().sqrt();
    if norm > 1e-12 {
        for value in &mut q {
            *value /= norm;
        }
    }
    q
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chiral_molecule_has_nonzero_vcd_strengths() {
        let conformer = crate::embed("C[C@H](O)Cl", 42);
        assert!(conformer.error.is_none(), "embed failed");
        let positions: Vec<[f64; 3]> = conformer
            .coords
            .as_chunks::<3>()
            .0
            .iter()
            .map(|chunk| [chunk[0], chunk[1], chunk[2]])
            .collect();
        let result = compute_vcd_spectrum(
            "C[C@H](O)Cl",
            &conformer.elements,
            &positions,
            &VcdRoaConfig::default(),
        )
        .unwrap();
        assert!(!result.modes.is_empty());
        assert!(result
            .modes
            .iter()
            .any(|mode| mode.rotational_strength.abs() > 1e-8));
    }

    #[test]
    #[ignore = "heavy: polarizability finite differences per mode"]
    fn roa_smoke_is_finite() {
        let conformer = crate::embed("CCO", 42);
        let positions: Vec<[f64; 3]> = conformer
            .coords
            .as_chunks::<3>()
            .0
            .iter()
            .map(|chunk| [chunk[0], chunk[1], chunk[2]])
            .collect();
        let result = compute_roa_spectrum(
            "CCO",
            &conformer.elements,
            &positions,
            &VcdRoaConfig::default(),
        )
        .unwrap();
        assert!(result
            .modes
            .iter()
            .all(|mode| mode.circular_intensity.is_finite()));
    }

    #[test]
    fn invalid_geometry_is_rejected() {
        assert!(compute_vcd_spectrum(
            "CCO",
            &[6, 6, 8],
            &[[0.0, 0.0, 0.0]],
            &VcdRoaConfig::default()
        )
        .is_err());
    }
}
