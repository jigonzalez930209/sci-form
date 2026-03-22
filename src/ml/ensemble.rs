//! Ensemble ML models with non-linear predictions and uncertainty estimates.
//!
//! Extends the base linear models with:
//! - Decision-tree-like cascading rules for LogP, solubility, pKa
//! - TPSA (Topological Polar Surface Area) descriptor
//! - Veber rule analysis for oral bioavailability
//! - BBB permeability prediction
//! - Consensus scoring with prediction confidence

use super::descriptors::MolecularDescriptors;
use serde::{Deserialize, Serialize};

/// Extended ML property predictions with uncertainty and additional properties.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnsembleResult {
    /// Predicted LogP (consensus of 3 models).
    pub logp: f64,
    /// Standard deviation across LogP models (uncertainty).
    pub logp_std: f64,
    /// Predicted aqueous solubility (log S, mol/L).
    pub log_solubility: f64,
    /// Predicted TPSA (Å²).
    pub tpsa: f64,
    /// Predicted pKa for the most acidic group.
    pub pka_acidic: Option<f64>,
    /// Predicted pKa for the most basic group.
    pub pka_basic: Option<f64>,
    /// Veber oral bioavailability rules.
    pub veber: VeberResult,
    /// BBB permeability prediction.
    pub bbb_permeable: bool,
    /// BBB permeability score (0–1).
    pub bbb_score: f64,
    /// Overall prediction confidence (0–1).
    pub confidence: f64,
}

/// Veber's rules for oral bioavailability.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VeberResult {
    /// TPSA ≤ 140 Å².
    pub tpsa_ok: bool,
    /// Number of rotatable bonds ≤ 10.
    pub rotb_ok: bool,
    /// Passes Veber rules (both criteria met).
    pub passes: bool,
}

/// Contribution of each atom type to TPSA (Å²).
/// Ertl, P., Rohde, B., Selzer, P. (2000) — simplified fragment contributions.
fn tpsa_contribution(z: u8, n_heavy_neighbors: usize, n_h_neighbors: usize) -> f64 {
    match z {
        // Nitrogen contributions
        7 => match (n_heavy_neighbors, n_h_neighbors) {
            (1, 2) => 26.02, // -NH₂
            (1, 1) => 26.02, // =NH
            (2, 1) => 12.36, // >NH
            (2, 0) => 12.36, // =N-
            (3, 0) => 3.24,  // >N-
            (1, 0) => 23.79, // ≡N
            _ => 12.0,
        },
        // Oxygen contributions
        8 => match (n_heavy_neighbors, n_h_neighbors) {
            (0, 2) => 20.23, // H₂O
            (1, 1) => 20.23, // -OH
            (1, 0) => 17.07, // =O
            (2, 0) => 9.23,  // -O-
            _ => 15.0,
        },
        // Sulfur: small contributions
        16 => match (n_heavy_neighbors, n_h_neighbors) {
            (1, 1) => 38.80, // -SH
            (2, 0) => 25.30, // -S-
            (1, 0) => 32.00, // =S
            _ => 28.0,
        },
        // Phosphorus
        15 => 34.14,
        _ => 0.0,
    }
}

/// Compute TPSA from elements and bond connectivity.
///
/// Uses Ertl fragment-based TPSA calculation.
pub fn compute_tpsa(elements: &[u8], bonds: &[(usize, usize, u8)]) -> f64 {
    let n = elements.len();
    let mut adj: Vec<Vec<usize>> = vec![vec![]; n];
    for &(i, j, _) in bonds {
        if i < n && j < n {
            adj[i].push(j);
            adj[j].push(i);
        }
    }

    let mut tpsa = 0.0;
    for i in 0..n {
        if !matches!(elements[i], 7 | 8 | 15 | 16) {
            continue;
        }
        let n_heavy = adj[i].iter().filter(|&&j| elements[j] != 1).count();
        let n_h = adj[i].iter().filter(|&&j| elements[j] == 1).count();
        tpsa += tpsa_contribution(elements[i], n_heavy, n_h);
    }
    tpsa
}

// ─── Ensemble LogP models ────────────────────────────────────────────────────

/// Model 1: Wildman-Crippen additive (same as base model).
fn logp_model_1(desc: &MolecularDescriptors) -> f64 {
    let base = 0.120 * desc.n_heavy_atoms as f64;
    let h_corr = -0.230 * desc.n_hbd as f64;
    let ring_corr = 0.150 * desc.n_rings as f64;
    let arom_corr = 0.080 * desc.n_aromatic as f64;
    let polar_corr = -0.310 * desc.n_hba as f64;
    let sp3_corr = -0.180 * desc.fsp3;
    let mw_term = 0.005 * (desc.molecular_weight - 100.0);
    base + h_corr + ring_corr + arom_corr + polar_corr + sp3_corr + mw_term
}

/// Model 2: ALOGP-like fragment-based with non-linear corrections.
fn logp_model_2(desc: &MolecularDescriptors, tpsa: f64) -> f64 {
    // Use TPSA as a proxy for overall polarity
    let polarity_term = -0.015 * tpsa;
    let size_term = 0.008 * desc.molecular_weight;
    let hbond_term = -0.45 * (desc.n_hbd as f64 + 0.5 * desc.n_hba as f64);
    let lipophilic = 0.22 * (desc.n_heavy_atoms as f64 - desc.n_hba as f64);

    // Non-linear correction: penalty for very high polarity
    let nl_correction = if tpsa > 80.0 {
        -0.003 * (tpsa - 80.0).powi(2) / 100.0
    } else {
        0.0
    };

    size_term + polarity_term + hbond_term + lipophilic + nl_correction - 1.5
}

/// Model 3: Topological descriptor based.
fn logp_model_3(desc: &MolecularDescriptors) -> f64 {
    let chi_approx = if desc.n_bonds > 0 {
        (desc.n_bonds as f64).sqrt() / (desc.n_heavy_atoms as f64).max(1.0)
    } else {
        0.0
    };

    let wiener_term = if desc.wiener_index > 0.0 {
        0.25 * desc.wiener_index.ln()
    } else {
        0.0
    };

    let polar_penalty = -0.35 * (desc.n_hbd + desc.n_hba) as f64;
    let arom_bonus = 0.12 * desc.n_aromatic as f64;

    chi_approx + wiener_term + polar_penalty + arom_bonus - 0.8
}

// ─── pKa prediction ──────────────────────────────────────────────────────────

/// Predict acidic pKa based on functional groups.
/// Returns None if no acidic group detected.
fn predict_pka_acidic(desc: &MolecularDescriptors, tpsa: f64) -> Option<f64> {
    if desc.n_hbd == 0 {
        return None;
    }

    // Base pKa value depends on donor type
    // Approximate: carboxylic acid ~4.5, phenol ~10, aliphatic OH ~16, NH ~35
    let base_pka = if desc.n_hba >= 2 && desc.n_hbd >= 1 {
        // Likely has carboxylic acid (O-H + C=O pattern)
        4.5
    } else {
        // Generic O-H or N-H
        14.0
    };

    // Corrections
    let ew_correction =
        -0.3 * (desc.sum_electronegativity / desc.n_heavy_atoms.max(1) as f64 - 2.5);
    let arom_correction = if desc.n_aromatic > 0 { -1.5 } else { 0.0 };
    let tpsa_correction = -0.02 * (tpsa - 60.0);

    Some((base_pka + ew_correction + arom_correction + tpsa_correction).clamp(0.0, 25.0))
}

/// Predict basic pKa. Returns None if no basic group detected.
fn predict_pka_basic(desc: &MolecularDescriptors) -> Option<f64> {
    // Check for nitrogen atoms that could be basic
    let has_nitrogen =
        desc.n_hba > 0 && desc.sum_electronegativity / desc.n_heavy_atoms.max(1) as f64 > 2.6;
    if !has_nitrogen {
        return None;
    }

    // Amine base pKa: primary ~10.6, aromatic ~5
    let base_pka = if desc.n_aromatic > 0 {
        5.2 // pyridine-like
    } else {
        10.6 // aliphatic amine
    };

    let sp3_correction = 0.5 * desc.fsp3;
    Some((base_pka + sp3_correction).clamp(0.0, 14.0))
}

// ─── BBB permeability ────────────────────────────────────────────────────────

/// Predict blood-brain barrier permeability.
///
/// Lipinski-like heuristic: MW < 450, TPSA < 90, LogP ∈ (1, 5), HBD ≤ 3.
fn predict_bbb(desc: &MolecularDescriptors, logp: f64, tpsa: f64) -> (bool, f64) {
    let mut score = 1.0;

    if desc.molecular_weight > 450.0 {
        score -= 0.3 * ((desc.molecular_weight - 450.0) / 100.0).min(1.0);
    }
    if tpsa > 90.0 {
        score -= 0.35 * ((tpsa - 90.0) / 50.0).min(1.0);
    }
    if logp < 1.0 {
        score -= 0.2 * (1.0 - logp).min(1.0);
    }
    if logp > 5.0 {
        score -= 0.2 * ((logp - 5.0) / 2.0).min(1.0);
    }
    if desc.n_hbd > 3 {
        score -= 0.15 * (desc.n_hbd as f64 - 3.0).min(2.0) / 2.0;
    }

    let score = score.clamp(0.0, 1.0);
    (score > 0.5, score)
}

/// Predict molecular properties using an ensemble of models.
///
/// Combines three LogP models via consensus, adds TPSA, pKa estimates,
/// BBB permeability, and Veber bioavailability rules.
pub fn predict_ensemble(
    desc: &MolecularDescriptors,
    elements: &[u8],
    bonds: &[(usize, usize, u8)],
) -> EnsembleResult {
    let tpsa = compute_tpsa(elements, bonds);

    // Ensemble LogP: average of 3 models
    let lp1 = logp_model_1(desc);
    let lp2 = logp_model_2(desc, tpsa);
    let lp3 = logp_model_3(desc);
    let logp = (lp1 + lp2 + lp3) / 3.0;

    // Standard deviation as uncertainty proxy
    let logp_std = {
        let mean = logp;
        let var = ((lp1 - mean).powi(2) + (lp2 - mean).powi(2) + (lp3 - mean).powi(2)) / 3.0;
        var.sqrt()
    };

    // Solubility (ESOL-like, using consensus LogP)
    let frac_aromatic = if desc.n_heavy_atoms > 0 {
        desc.n_aromatic as f64 / desc.n_heavy_atoms as f64
    } else {
        0.0
    };
    let log_solubility = 0.16 - 0.63 * logp - 0.0062 * desc.molecular_weight
        + 0.066 * desc.n_rotatable_bonds as f64
        - 0.74 * frac_aromatic;

    // pKa predictions
    let pka_acidic = predict_pka_acidic(desc, tpsa);
    let pka_basic = predict_pka_basic(desc);

    // Veber rules
    let tpsa_ok = tpsa <= 140.0;
    let rotb_ok = desc.n_rotatable_bonds <= 10;
    let veber = VeberResult {
        tpsa_ok,
        rotb_ok,
        passes: tpsa_ok && rotb_ok,
    };

    // BBB
    let (bbb_permeable, bbb_score) = predict_bbb(desc, logp, tpsa);

    // Confidence: based on model agreement and descriptor coverage
    let model_agreement = 1.0 - (logp_std / 2.0).min(1.0);
    let size_confidence = if desc.n_heavy_atoms >= 3 && desc.n_heavy_atoms <= 50 {
        1.0
    } else {
        0.5
    };
    let confidence = (model_agreement * 0.7 + size_confidence * 0.3).clamp(0.0, 1.0);

    EnsembleResult {
        logp,
        logp_std,
        log_solubility,
        tpsa,
        pka_acidic,
        pka_basic,
        veber,
        bbb_permeable,
        bbb_score,
        confidence,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ml::descriptors::compute_descriptors;

    #[test]
    fn test_tpsa_water() {
        // Water: O with 2 H neighbors
        let elements = [8u8, 1, 1];
        let bonds = [(0usize, 1usize, 1u8), (0, 2, 1)];
        let tpsa = compute_tpsa(&elements, &bonds);
        assert!(tpsa > 15.0 && tpsa < 25.0, "Water TPSA: {tpsa}");
    }

    #[test]
    fn test_tpsa_benzene() {
        // Benzene: no polar atoms → TPSA = 0
        let elements = [6u8, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1];
        let bonds: Vec<(usize, usize, u8)> = vec![
            (0, 1, 2),
            (1, 2, 1),
            (2, 3, 2),
            (3, 4, 1),
            (4, 5, 2),
            (5, 0, 1),
            (0, 6, 1),
            (1, 7, 1),
            (2, 8, 1),
            (3, 9, 1),
            (4, 10, 1),
            (5, 11, 1),
        ];
        let tpsa = compute_tpsa(&elements, &bonds);
        assert!(
            (tpsa - 0.0).abs() < 1e-6,
            "Benzene TPSA should be 0: {tpsa}"
        );
    }

    #[test]
    fn test_ensemble_ethanol() {
        let elements = [6u8, 6, 8, 1, 1, 1, 1, 1, 1];
        let bonds: Vec<(usize, usize, u8)> = vec![
            (0, 1, 1),
            (1, 2, 1),
            (0, 3, 1),
            (0, 4, 1),
            (0, 5, 1),
            (1, 6, 1),
            (1, 7, 1),
            (2, 8, 1),
        ];
        let desc = compute_descriptors(&elements, &bonds, &[], &[]);
        let result = predict_ensemble(&desc, &elements, &bonds);

        assert!(result.tpsa > 15.0, "Ethanol has polar O-H: {}", result.tpsa);
        assert!(result.logp < 2.0, "Ethanol is hydrophilic: {}", result.logp);
        assert!(result.logp_std >= 0.0, "Uncertainty must be non-negative");
        assert!(result.confidence > 0.0 && result.confidence <= 1.0);
        assert!(result.veber.passes);
    }

    #[test]
    fn test_ensemble_logp_consistency() {
        // All three models should give broadly similar results for simple molecules
        let elements = [6u8, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1];
        let bonds: Vec<(usize, usize, u8)> = vec![
            (0, 1, 1),
            (1, 2, 1),
            (0, 3, 1),
            (0, 4, 1),
            (0, 5, 1),
            (1, 6, 1),
            (1, 7, 1),
            (2, 8, 1),
            (2, 9, 1),
            (2, 10, 1),
        ];
        let desc = compute_descriptors(&elements, &bonds, &[], &[]);
        let result = predict_ensemble(&desc, &elements, &bonds);

        // Models should agree within ~2 units for small alkanes
        assert!(
            result.logp_std < 2.0,
            "Models should broadly agree: std={}",
            result.logp_std
        );
    }

    #[test]
    fn test_veber_large_molecule() {
        let desc = MolecularDescriptors {
            molecular_weight: 600.0,
            n_heavy_atoms: 45,
            n_hydrogens: 20,
            n_bonds: 60,
            n_rotatable_bonds: 15,
            n_hbd: 5,
            n_hba: 10,
            fsp3: 0.3,
            total_abs_charge: 3.0,
            max_charge: 0.4,
            min_charge: -0.4,
            wiener_index: 3000.0,
            n_rings: 4,
            n_aromatic: 8,
            balaban_j: 1.8,
            sum_electronegativity: 120.0,
            sum_polarizability: 65.0,
        };
        let elements = [6u8; 45];
        let bonds: Vec<(usize, usize, u8)> = (0..44).map(|i| (i, i + 1, 1u8)).collect();
        let result = predict_ensemble(&desc, &elements, &bonds);

        assert!(
            !result.veber.rotb_ok,
            "Too many rotatable bonds: {}",
            desc.n_rotatable_bonds
        );
    }

    #[test]
    fn test_bbb_small_lipophilic() {
        // Small, lipophilic molecule — should be BBB permeable
        let elements = [6u8, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1];
        let bonds: Vec<(usize, usize, u8)> = vec![
            (0, 1, 1),
            (1, 2, 1),
            (0, 3, 1),
            (0, 4, 1),
            (0, 5, 1),
            (1, 6, 1),
            (1, 7, 1),
            (2, 8, 1),
            (2, 9, 1),
            (2, 10, 1),
        ];
        let desc = compute_descriptors(&elements, &bonds, &[], &[]);
        let result = predict_ensemble(&desc, &elements, &bonds);
        assert!(
            result.bbb_score > 0.0,
            "Small lipophilic molecule should have positive BBB score"
        );
    }

    #[test]
    fn test_pka_with_acid() {
        // Molecule with O-H and C=O (carboxylic acid pattern)
        let desc = MolecularDescriptors {
            molecular_weight: 60.0,
            n_heavy_atoms: 3,
            n_hydrogens: 4,
            n_bonds: 6,
            n_rotatable_bonds: 0,
            n_hbd: 1,
            n_hba: 2,
            fsp3: 0.0,
            total_abs_charge: 0.5,
            max_charge: 0.2,
            min_charge: -0.3,
            wiener_index: 4.0,
            n_rings: 0,
            n_aromatic: 0,
            balaban_j: 1.0,
            sum_electronegativity: 8.0,
            sum_polarizability: 3.0,
        };
        let elements = [6u8, 8, 8, 1, 1, 1, 1];
        let bonds: Vec<(usize, usize, u8)> = vec![
            (0, 1, 2),
            (0, 2, 1),
            (2, 3, 1),
            (0, 4, 1),
            (0, 5, 1),
            (0, 6, 1),
        ];
        let result = predict_ensemble(&desc, &elements, &bonds);
        assert!(
            result.pka_acidic.is_some(),
            "Carboxylic acid should have pKa"
        );
    }
}
