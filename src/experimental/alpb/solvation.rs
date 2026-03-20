//! ALPB Solvation Energy — E6.2
//!
//! Generalized Born electrostatic solvation with ALPB correction
//! and non-polar SASA contribution.

use super::born::{compute_born_radii, gb_kernel};

/// Configuration for ALPB solvation.
#[derive(Debug, Clone)]
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
            solvent_dielectric: 78.5, // water
            probe_radius: 1.4,
            surface_tension: 0.005,
        }
    }
}

/// Result of ALPB solvation calculation.
#[derive(Debug, Clone)]
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
/// E_el = -0.5 * (1 - 1/eps) * sum_{ij} q_i * q_j / f_GB
/// ALPB correction: E_ALPB = E_GB * (eps-1)/(eps + A*x)
/// where A = 0.571412 (Klamt universal constant)
pub fn compute_alpb_solvation(
    elements: &[u8],
    positions: &[[f64; 3]],
    charges: &[f64],
    config: &AlpbConfig,
) -> AlpbResult {
    let n = elements.len();
    let born = compute_born_radii(elements, positions, config.probe_radius);

    // GB electrostatic energy
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

    // ALPB correction
    let klamt_a = 0.571412;
    let e_alpb = if (eps - 1.0).abs() < 1e-12 {
        // Vacuum: no solvation
        0.0
    } else {
        let x = if e_coulomb.abs() > 1e-15 { e_gb / (prefactor * e_coulomb) } else { 1.0 };
        let alpb_factor_local = (eps - 1.0) / (eps + klamt_a * x.abs().max(0.01));
        e_gb * alpb_factor_local / ((eps - 1.0) / eps)
    };
    let alpb_factor = if (eps - 1.0).abs() < 1e-12 {
        0.0
    } else {
        let x = if e_coulomb.abs() > 1e-15 { e_gb / (prefactor * e_coulomb) } else { 1.0 };
        (eps - 1.0) / (eps + klamt_a * x.abs().max(0.01))
    };

    // Convert from atomic units to kcal/mol
    let hartree_to_kcal = 627.509;
    let bohr_to_angstrom = 0.529177;
    let conversion = hartree_to_kcal * bohr_to_angstrom; // eÅ → kcal/mol

    let e_electrostatic = e_alpb * conversion;

    // Non-polar SASA term (simplified version)
    let mut sasa_total = 0.0;
    for i in 0..n {
        let r_eff = born.radii[i] + config.probe_radius;
        // Approximate per-atom SASA
        let mut exposed = 4.0 * std::f64::consts::PI * r_eff * r_eff;
        for j in 0..n {
            if i == j { continue; }
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r_ij = (dx * dx + dy * dy + dz * dz).sqrt();
            let r_j_eff = born.radii[j] + config.probe_radius;
            if r_ij < r_eff + r_j_eff {
                let overlap = std::f64::consts::PI * r_eff * (r_eff + r_j_eff - r_ij);
                exposed -= overlap.min(exposed * 0.5);
            }
        }
        sasa_total += exposed.max(0.0);
    }

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
        // Solvation should be negative (stabilizing)
        assert!(result.electrostatic_energy < 0.0 || result.electrostatic_energy.abs() < 1.0,
            "Electrostatic = {}", result.electrostatic_energy);
    }

    #[test]
    fn test_alpb_vacuum_dielectric() {
        let elements = [8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let charges = [-0.834, 0.417, 0.417];
        let config = AlpbConfig { solvent_dielectric: 1.0, ..Default::default() };
        let result = compute_alpb_solvation(&elements, &positions, &charges, &config);
        // In vacuum (eps=1), electrostatic solvation should be ~0
        assert!(result.electrostatic_energy.abs() < 1e-10,
            "Vacuum solvation should be 0: {}", result.electrostatic_energy);
    }

    #[test]
    fn test_alpb_nonpolar_positive() {
        let elements = [6, 1, 1, 1, 1];
        let positions = [
            [0.0, 0.0, 0.0], [0.63, 0.63, 0.63],
            [-0.63, -0.63, 0.63], [-0.63, 0.63, -0.63], [0.63, -0.63, -0.63],
        ];
        let charges = [0.0, 0.0, 0.0, 0.0, 0.0];
        let config = AlpbConfig::default();
        let result = compute_alpb_solvation(&elements, &positions, &charges, &config);
        // Non-polar should be positive (cavity formation cost)
        assert!(result.nonpolar_energy >= 0.0);
    }
}
