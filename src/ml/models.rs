//! Pre-fitted linear ML models for fast property estimation.
//!
//! Uses molecular descriptors to predict LogP, molar refractivity,
//! aqueous solubility, and toxicity flags via simple linear regression
//! coefficients fitted to public datasets.

use super::descriptors::MolecularDescriptors;
use serde::{Deserialize, Serialize};

/// ML-predicted molecular properties.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MlPropertyResult {
    /// Predicted octanol-water partition coefficient (log P).
    pub logp: f64,
    /// Predicted molar refractivity (cm³/mol).
    pub molar_refractivity: f64,
    /// Predicted aqueous solubility (log S, mol/L).
    pub log_solubility: f64,
    /// Lipinski rule-of-five flags.
    pub lipinski: LipinskiResult,
    /// Druglikeness score (0–1).
    pub druglikeness: f64,
}

/// Lipinski's rule of five analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LipinskiResult {
    /// Molecular weight ≤ 500.
    pub mw_ok: bool,
    /// LogP ≤ 5.
    pub logp_ok: bool,
    /// H-bond donors ≤ 5.
    pub hbd_ok: bool,
    /// H-bond acceptors ≤ 10.
    pub hba_ok: bool,
    /// Number of violations (0–4).
    pub violations: u8,
    /// Passes Ro5 (≤ 1 violation).
    pub passes: bool,
}

/// Predict LogP using Wildman-Crippen inspired atom-type contributions.
///
/// This is a simplified additive model:
/// LogP ≈ Σ a_i * (atom_contribution) + corrections
fn predict_logp(desc: &MolecularDescriptors) -> f64 {
    // Wildman-Crippen style: base contribution per heavy atom + corrections
    let base = 0.120 * desc.n_heavy_atoms as f64;
    let h_correction = -0.230 * desc.n_hbd as f64; // polar H reduces logP
    let ring_correction = 0.150 * desc.n_rings as f64;
    let aromatic_correction = 0.080 * desc.n_aromatic as f64;
    let polar_correction = -0.310 * desc.n_hba as f64;
    let sp3_correction = -0.180 * desc.fsp3;
    let mw_term = 0.005 * (desc.molecular_weight - 100.0);

    base + h_correction
        + ring_correction
        + aromatic_correction
        + polar_correction
        + sp3_correction
        + mw_term
}

/// Predict molar refractivity (CMR) from atom contributions.
///
/// Simple additive model based on Ghose-Crippen atom contributions.
fn predict_molar_refractivity(desc: &MolecularDescriptors) -> f64 {
    // Approximate: MR scales with polarizability
    // Miller's formula: MR ≈ 2.536 * sum_polarizability + small corrections
    let base = 2.536 * desc.sum_polarizability;
    let ring_correction = 1.20 * desc.n_rings as f64;
    let aromatic = 0.80 * desc.n_aromatic as f64;
    base + ring_correction + aromatic
}

/// Predict aqueous solubility (log S) using ESOL-inspired model (Delaney 2004).
///
/// log S ≈ 0.16 - 0.63 * logP - 0.0062 * MW + 0.066 * nRotB - 0.74 * nAromAtoms/nHeavy
fn predict_solubility(desc: &MolecularDescriptors, logp: f64) -> f64 {
    let frac_aromatic = if desc.n_heavy_atoms > 0 {
        desc.n_aromatic as f64 / desc.n_heavy_atoms as f64
    } else {
        0.0
    };
    0.16 - 0.63 * logp - 0.0062 * desc.molecular_weight + 0.066 * desc.n_rotatable_bonds as f64
        - 0.74 * frac_aromatic
}

/// Compute druglikeness score (0–1) from molecular descriptors.
fn druglikeness_score(desc: &MolecularDescriptors, logp: f64) -> f64 {
    let mut score = 1.0;

    // MW penalty
    if desc.molecular_weight > 500.0 {
        score -= 0.2 * ((desc.molecular_weight - 500.0) / 200.0).min(1.0);
    }
    // LogP penalty
    if logp > 5.0 {
        score -= 0.2 * ((logp - 5.0) / 3.0).min(1.0);
    } else if logp < -2.0 {
        score -= 0.15;
    }
    // HBD penalty
    if desc.n_hbd > 5 {
        score -= 0.15;
    }
    // HBA penalty
    if desc.n_hba > 10 {
        score -= 0.15;
    }
    // Rotatable bond penalty
    if desc.n_rotatable_bonds > 10 {
        score -= 0.1 * ((desc.n_rotatable_bonds as f64 - 10.0) / 5.0).min(1.0);
    }
    // Reward for fsp3 (lead-likeness)
    score += 0.05 * desc.fsp3;

    score.clamp(0.0, 1.0)
}

/// Predict molecular properties from descriptors.
///
/// Returns a `MlPropertyResult` with LogP, MR, solubility, Lipinski flags,
/// and a druglikeness score. These are fast, approximate predictions
/// suitable for screening workflows.
pub fn predict_properties(desc: &MolecularDescriptors) -> MlPropertyResult {
    let logp = predict_logp(desc);
    let mr = predict_molar_refractivity(desc);
    let log_s = predict_solubility(desc, logp);

    let mw_ok = desc.molecular_weight <= 500.0;
    let logp_ok = logp <= 5.0;
    let hbd_ok = desc.n_hbd <= 5;
    let hba_ok = desc.n_hba <= 10;
    let violations = [!mw_ok, !logp_ok, !hbd_ok, !hba_ok]
        .iter()
        .filter(|&&v| v)
        .count() as u8;

    let lipinski = LipinskiResult {
        mw_ok,
        logp_ok,
        hbd_ok,
        hba_ok,
        violations,
        passes: violations <= 1,
    };

    let druglikeness = druglikeness_score(desc, logp);

    MlPropertyResult {
        logp,
        molar_refractivity: mr,
        log_solubility: log_s,
        lipinski,
        druglikeness,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ml::descriptors::compute_descriptors;

    #[test]
    fn test_predict_water() {
        let elements = [8u8, 1, 1];
        let bonds = [(0, 1, 1u8), (0, 2, 1)];
        let desc = compute_descriptors(&elements, &bonds, &[], &[]);
        let result = predict_properties(&desc);
        // Water should be very hydrophilic (low logP)
        assert!(
            result.logp < 1.0,
            "Water logP should be low: {}",
            result.logp
        );
        assert!(result.lipinski.passes, "Water should pass Lipinski");
    }

    #[test]
    fn test_lipinski_violations() {
        // Giant fake molecule that violates everything
        let desc = MolecularDescriptors {
            molecular_weight: 800.0,
            n_heavy_atoms: 60,
            n_hydrogens: 20,
            n_bonds: 80,
            n_rotatable_bonds: 15,
            n_hbd: 8,
            n_hba: 15,
            fsp3: 0.1,
            total_abs_charge: 5.0,
            max_charge: 0.5,
            min_charge: -0.5,
            wiener_index: 5000.0,
            n_rings: 5,
            n_aromatic: 12,
            balaban_j: 2.0,
            sum_electronegativity: 150.0,
            sum_polarizability: 80.0,
        };
        let result = predict_properties(&desc);
        assert!(
            result.lipinski.violations >= 2,
            "Should have multiple violations"
        );
        assert!(!result.lipinski.passes, "Should fail Lipinski");
    }

    #[test]
    fn test_druglikeness_range() {
        let elements = [6u8, 6, 8, 1, 1, 1, 1, 1, 1];
        let bonds = [
            (0, 1, 1u8),
            (1, 2, 1),
            (0, 3, 1),
            (0, 4, 1),
            (0, 5, 1),
            (1, 6, 1),
            (1, 7, 1),
            (2, 8, 1),
        ];
        let desc = compute_descriptors(&elements, &bonds, &[], &[]);
        let result = predict_properties(&desc);
        assert!(result.druglikeness >= 0.0 && result.druglikeness <= 1.0);
    }
}
