//! ALPB Solvation Energy — Core
//!
//! Generalized Born electrostatic solvation with ALPB correction
//! and non-polar SASA contribution.
//!
//! ## Units
//! - Positions and Born radii: Ångström (Å)
//! - Charges: elementary charges (e)
//! - Energies (result): kcal/mol
//! - Surface tension: kcal/(mol·Å²)
//! - Coulomb constant: 332.063 kcal·Å/e² = e²·Nₐ/(4πε₀)

use super::born::{compute_born_radii, gb_kernel};
use serde::{Deserialize, Serialize};

/// Configuration for ALPB solvation.
///
/// Default: water (`ε = 78.5`), standard probe radius (1.4 Å),
/// surface tension 0.005 kcal/(mol·Å²).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlpbConfig {
    /// Solvent dielectric constant.
    pub solvent_dielectric: f64,
    /// Probe radius for SASA (Å).
    pub probe_radius: f64,
    /// Surface tension for non-polar term (kcal/(mol·Å²)).
    pub surface_tension: f64,
}

impl Default for AlpbConfig {
    fn default() -> Self {
        Self {
            solvent_dielectric: 78.5,
            probe_radius: 1.4,
            surface_tension: 0.005,
        }
    }
}

/// Result of ALPB solvation calculation.
///
/// Energies are in kcal/mol. Born radii are in Å.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlpbResult {
    /// Electrostatic solvation energy (kcal/mol).
    pub electrostatic_energy: f64,
    /// Non-polar solvation energy (kcal/mol).
    pub nonpolar_energy: f64,
    /// Total solvation energy (kcal/mol).
    pub total_energy: f64,
    /// Born radii (Å).
    pub born_radii: Vec<f64>,
    /// ALPB correction factor.
    pub alpb_factor: f64,
}

/// Compute ALPB solvation energy.
///
/// ALPB applies a Klamt correction to standard GB:
/// E_ALPB = E_GB · (ε-1)/(ε + A·x)  where A = 0.571412
pub fn compute_alpb_solvation(
    elements: &[u8],
    positions: &[[f64; 3]],
    charges: &[f64],
    config: &AlpbConfig,
) -> AlpbResult {
    let n = elements.len();
    let born = compute_born_radii(elements, positions, config.probe_radius);

    let eps = config.solvent_dielectric;
    let prefactor = -0.5 * (1.0 - 1.0 / eps);

    let mut e_gb = 0.0;
    let mut e_coulomb = 0.0;

    for i in 0..n {
        for j in i..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r_ij = (dx * dx + dy * dy + dz * dz).sqrt();

            let f = gb_kernel(r_ij, born.radii[i], born.radii[j]);
            let factor = if i == j { 1.0 } else { 2.0 };

            e_gb += factor * charges[i] * charges[j] / f;

            if i != j && r_ij > 1e-10 {
                e_coulomb += factor * charges[i] * charges[j] / r_ij;
            }
        }
    }

    e_gb *= prefactor;

    // ALPB correction (Klamt scheme)
    let klamt_a = 0.571412;
    let (alpb_factor, e_alpb) = if (eps - 1.0).abs() < 1e-12 {
        (0.0, 0.0)
    } else {
        let x = if e_coulomb.abs() > 1e-15 {
            e_gb / (prefactor * e_coulomb)
        } else {
            1.0
        };
        let af = (eps - 1.0) / (eps + klamt_a * x.abs().max(0.01));
        let e = e_gb * af / ((eps - 1.0) / eps);
        (af, e)
    };

    // Coulomb constant in kcal·Å/e²: e²·Nₐ / (4π ε₀) ≈ 332.06
    let coulomb_kcal_ang = 332.063;
    let e_electrostatic = e_alpb * coulomb_kcal_ang;

    // Non-polar SASA term — use the full Shrake-Rupley implementation
    let sasa_result = crate::surface::sasa::compute_sasa(elements, positions, Some(config.probe_radius), None);
    let sasa_total = sasa_result.total_sasa;

    let e_nonpolar = config.surface_tension * sasa_total;

    AlpbResult {
        electrostatic_energy: e_electrostatic,
        nonpolar_energy: e_nonpolar,
        total_energy: e_electrostatic + e_nonpolar,
        born_radii: born.radii,
        alpb_factor,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alpb_water_in_water() {
        let elements = [8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let charges = [-0.834, 0.417, 0.417];
        let config = AlpbConfig::default();
        let result = compute_alpb_solvation(&elements, &positions, &charges, &config);
        assert!(result.total_energy.is_finite());
        assert!(
            result.electrostatic_energy < 0.0 || result.electrostatic_energy.abs() < 1.0,
            "Electrostatic = {}",
            result.electrostatic_energy
        );
    }

    #[test]
    fn test_alpb_vacuum_dielectric() {
        let elements = [8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let charges = [-0.834, 0.417, 0.417];
        let config = AlpbConfig {
            solvent_dielectric: 1.0,
            ..Default::default()
        };
        let result = compute_alpb_solvation(&elements, &positions, &charges, &config);
        assert!(
            result.electrostatic_energy.abs() < 1e-10,
            "Vacuum solvation should be 0: {}",
            result.electrostatic_energy
        );
    }

    #[test]
    fn test_alpb_nonpolar_positive() {
        let elements = [6, 1, 1, 1, 1];
        let positions = [
            [0.0, 0.0, 0.0],
            [0.63, 0.63, 0.63],
            [-0.63, -0.63, 0.63],
            [-0.63, 0.63, -0.63],
            [0.63, -0.63, -0.63],
        ];
        let charges = [0.0, 0.0, 0.0, 0.0, 0.0];
        let config = AlpbConfig::default();
        let result = compute_alpb_solvation(&elements, &positions, &charges, &config);
        assert!(result.nonpolar_energy >= 0.0);
    }

    #[test]
    fn test_alpb_factor_bounded() {
        let elements = [8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let charges = [-0.834, 0.417, 0.417];
        let config = AlpbConfig::default();
        let result = compute_alpb_solvation(&elements, &positions, &charges, &config);
        assert!(
            result.alpb_factor >= 0.0 && result.alpb_factor <= 1.0,
            "ALPB factor should be in [0,1]: {}",
            result.alpb_factor
        );
    }
}
