//! Implicit solvation modeling: ASP non-polar solvation and Generalized Born (GB).
//!
//! Implements the Still Generalized Born model with Hawkins-Cramer-Truhlar (HCT)
//! descreening for effective Born radii, combined with non-polar SASA-based solvation.

use serde::{Deserialize, Serialize};

/// Atomic solvation parameters (ASP) for non-polar solvation energy in cal/(mol·Å²).
/// Values from Wesson & Eisenberg, Protein Science 1992.
fn atomic_solvation_parameter(z: u8) -> f64 {
    match z {
        1 => 7.0,    // H: hydrophobic
        6 => 12.0,   // C: hydrocarbon
        7 => -116.0, // N: polar
        8 => -166.0, // O: polar
        9 => -5.0,   // F: slightly polar
        15 => -20.0, // P: polar
        16 => -32.0, // S: slightly polar
        17 => 18.0,  // Cl: hydrophobic
        35 => 22.0,  // Br: hydrophobic
        53 => 28.0,  // I: hydrophobic
        _ => 0.0,
    }
}

/// Intrinsic Born radii (Å) used as initial radii for GB calculations.
/// Based on Bondi vdW radii with scaling factor 0.8.
fn intrinsic_born_radius(z: u8) -> f64 {
    let vdw = match z {
        1 => 1.20,
        5 => 1.92,
        6 => 1.70,
        7 => 1.55,
        8 => 1.52,
        9 => 1.47,
        14 => 2.10,
        15 => 1.80,
        16 => 1.80,
        17 => 1.75,
        35 => 1.85,
        53 => 1.98,
        _ => 1.70,
    };
    vdw * 0.8 // HCT scaling
}

/// HCT descreening parameter by element.
fn hct_descreening_scale(z: u8) -> f64 {
    match z {
        1 => 0.85,
        6 => 0.72,
        7 => 0.79,
        8 => 0.85,
        9 => 0.88,
        15 => 0.86,
        16 => 0.96,
        _ => 0.80,
    }
}

/// Result of non-polar solvation energy calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NonPolarSolvation {
    /// Total non-polar solvation energy in kcal/mol.
    pub energy_kcal_mol: f64,
    /// Per-atom solvation contributions in kcal/mol.
    pub atom_contributions: Vec<f64>,
    /// Per-atom SASA values used (Å²).
    pub atom_sasa: Vec<f64>,
    /// Total SASA (Å²).
    pub total_sasa: f64,
}

/// Result of Generalized Born electrostatic solvation calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GbSolvation {
    /// Electrostatic solvation energy in kcal/mol.
    pub electrostatic_energy_kcal_mol: f64,
    /// Non-polar solvation energy in kcal/mol.
    pub nonpolar_energy_kcal_mol: f64,
    /// Total solvation energy (electrostatic + non-polar) in kcal/mol.
    pub total_energy_kcal_mol: f64,
    /// Effective Born radii for each atom (Å).
    pub born_radii: Vec<f64>,
    /// Partial charges used.
    pub charges: Vec<f64>,
    /// Solvent dielectric constant used.
    pub solvent_dielectric: f64,
    /// Solute dielectric constant used.
    pub solute_dielectric: f64,
}

/// Compute non-polar solvation energy from per-atom SASA and atomic solvation parameters.
///
/// ΔG_np = Σ_i σ_i · A_i where σ_i is the ASP and A_i is the SASA of atom i.
pub fn compute_nonpolar_solvation(
    elements: &[u8],
    positions: &[[f64; 3]],
    probe_radius: Option<f64>,
) -> NonPolarSolvation {
    let sasa = crate::surface::sasa::compute_sasa(elements, positions, probe_radius, Some(960));

    let mut atom_contributions = Vec::with_capacity(elements.len());
    for (i, &z) in elements.iter().enumerate() {
        let asp = atomic_solvation_parameter(z);
        // Convert cal/(mol·Å²) to kcal/(mol·Å²)
        let contrib = asp * sasa.atom_sasa[i] / 1000.0;
        atom_contributions.push(contrib);
    }

    let energy_kcal_mol: f64 = atom_contributions.iter().sum();

    NonPolarSolvation {
        energy_kcal_mol,
        atom_contributions,
        atom_sasa: sasa.atom_sasa,
        total_sasa: sasa.total_sasa,
    }
}

/// Compute effective Born radii using the Hawkins-Cramer-Truhlar (HCT) pairwise descreening model.
///
/// R_i^eff = 1 / (1/ρ_i - I_i) where I_i is the descreening integral.
pub fn compute_born_radii(elements: &[u8], positions: &[[f64; 3]]) -> Vec<f64> {
    let n = elements.len();
    let mut born_radii = Vec::with_capacity(n);

    for i in 0..n {
        let rho_i = intrinsic_born_radius(elements[i]);
        let mut integral = 0.0;

        for j in 0..n {
            if i == j {
                continue;
            }

            let rho_j = intrinsic_born_radius(elements[j]);
            let scale_j = hct_descreening_scale(elements[j]);
            let scaled_rj = rho_j * scale_j;

            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let rij = (dx * dx + dy * dy + dz * dz).sqrt();

            if rij > rho_i + scaled_rj {
                // Atom j is completely outside atom i's sphere
                let term = 0.5
                    * (1.0 / (rij - scaled_rj) - 1.0 / (rij + scaled_rj)
                        + scaled_rj / (rij * rij - scaled_rj * scaled_rj)
                            * (rij / (2.0 * (rij * rij - scaled_rj * scaled_rj).abs().max(1e-10))
                                + 0.5 * (1.0 / rij).ln().exp() * 0.0));
                // Simplified HCT integral: I_ij ≈ ½(1/(r-s) - 1/(r+s)) + s/(r²-s²)·ln(r/s)/(2r)
                let ljr = if rij > scaled_rj && scaled_rj > 1e-10 {
                    (rij / scaled_rj).ln()
                } else {
                    0.0
                };
                let denom1 = (rij - scaled_rj).max(1e-10);
                let denom2 = rij + scaled_rj;
                let denom3 = (rij * rij - scaled_rj * scaled_rj).abs().max(1e-10);
                let _ = term;
                integral += 0.5 * (1.0 / denom1 - 1.0 / denom2)
                    + scaled_rj * ljr / (2.0 * rij * denom3.max(1e-10));
            } else if rij + rho_i > scaled_rj {
                // Partial overlap
                let denom = (rij - scaled_rj).abs().max(1e-10);
                integral += 0.5 * (1.0 / denom - 1.0 / (rij + scaled_rj));
            }
            // If fully enclosed, skip (integral contribution is handled differently)
        }

        let inv_r = 1.0 / rho_i - integral;
        let born_r = if inv_r > 1e-10 { 1.0 / inv_r } else { 50.0 }; // cap at 50 Å
        born_radii.push(born_r.max(rho_i)); // Born radius should not be smaller than intrinsic
    }

    born_radii
}

/// Compute Generalized Born electrostatic solvation energy using the Still equation.
///
/// ΔG_GB = -½(1/ε_in - 1/ε_out) Σ_{i,j} q_i·q_j / f_GB(r_ij, R_i, R_j)
///
/// where f_GB = sqrt(r_ij² + R_i·R_j·exp(-r_ij²/(4·R_i·R_j)))
pub fn compute_gb_solvation(
    elements: &[u8],
    positions: &[[f64; 3]],
    charges: &[f64],
    solvent_dielectric: Option<f64>,
    solute_dielectric: Option<f64>,
    probe_radius: Option<f64>,
) -> GbSolvation {
    let n = elements.len();
    let eps_out = solvent_dielectric.unwrap_or(78.5); // water at 25°C
    let eps_in = solute_dielectric.unwrap_or(1.0);

    let born_radii = compute_born_radii(elements, positions);

    // Electrostatic GB energy
    let prefactor = -332.05 * 0.5 * (1.0 / eps_in - 1.0 / eps_out); // 332.05 = e²/(4πε₀) in kcal·Å/mol
    let mut elec_energy = 0.0;

    for i in 0..n {
        for j in i..n {
            let qi = charges[i];
            let qj = charges[j];
            if qi.abs() < 1e-12 && qj.abs() < 1e-12 {
                continue;
            }

            let rij_sq = if i == j {
                0.0
            } else {
                let dx = positions[i][0] - positions[j][0];
                let dy = positions[i][1] - positions[j][1];
                let dz = positions[i][2] - positions[j][2];
                dx * dx + dy * dy + dz * dz
            };

            let ri_rj = born_radii[i] * born_radii[j];
            let f_gb = (rij_sq + ri_rj * (-rij_sq / (4.0 * ri_rj).max(1e-10)).exp()).sqrt();

            let factor = if i == j { 1.0 } else { 2.0 }; // count (i,j) and (j,i)
            elec_energy += factor * prefactor * qi * qj / f_gb;
        }
    }

    // Non-polar contribution
    let nonpolar = compute_nonpolar_solvation(elements, positions, probe_radius);

    GbSolvation {
        electrostatic_energy_kcal_mol: elec_energy,
        nonpolar_energy_kcal_mol: nonpolar.energy_kcal_mol,
        total_energy_kcal_mol: elec_energy + nonpolar.energy_kcal_mol,
        born_radii,
        charges: charges.to_vec(),
        solvent_dielectric: eps_out,
        solute_dielectric: eps_in,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nonpolar_solvation_methane() {
        // Single carbon → positive ASP → positive solvation
        let elements = vec![6u8];
        let positions = vec![[0.0, 0.0, 0.0]];
        let result = compute_nonpolar_solvation(&elements, &positions, None);
        assert!(
            result.energy_kcal_mol > 0.0,
            "Carbon ASP should be positive"
        );
        assert!(result.total_sasa > 0.0);
    }

    #[test]
    fn test_born_radii_positive() {
        let elements = vec![8u8, 1, 1];
        let positions = vec![[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let radii = compute_born_radii(&elements, &positions);
        assert_eq!(radii.len(), 3);
        for r in &radii {
            assert!(*r > 0.0, "Born radius should be positive, got {}", r);
            assert!(*r <= 50.0, "Born radius should be <= 50 Å, got {}", r);
        }
    }

    #[test]
    fn test_gb_solvation_water() {
        let elements = vec![8u8, 1, 1];
        let positions = vec![[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let charges = vec![-0.834, 0.417, 0.417]; // TIP3P charges
        let result = compute_gb_solvation(&elements, &positions, &charges, None, None, None);
        // Water should have negative electrostatic solvation (stabilizing)
        assert!(
            result.electrostatic_energy_kcal_mol < 0.0,
            "Water GB energy should be negative, got {}",
            result.electrostatic_energy_kcal_mol
        );
    }

    #[test]
    fn test_neutral_molecule_near_zero() {
        // All-zero charges → near-zero electrostatic solvation
        let elements = vec![6u8, 1, 1, 1, 1];
        let positions = vec![
            [0.0, 0.0, 0.0],
            [1.09, 0.0, 0.0],
            [-0.36, 1.03, 0.0],
            [-0.36, -0.52, 0.89],
            [-0.36, -0.52, -0.89],
        ];
        let charges = vec![0.0, 0.0, 0.0, 0.0, 0.0];
        let result = compute_gb_solvation(&elements, &positions, &charges, None, None, None);
        assert!(
            result.electrostatic_energy_kcal_mol.abs() < 1e-6,
            "Zero-charge should give zero electrostatic solvation"
        );
    }
}
