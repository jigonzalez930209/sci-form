//! QM/MM with mechanical embedding.
//!
//! The first v2.2 multiscale prototype combines a quantum region with a
//! classical point-charge environment. The current embedding is **mechanical**:
//! the QM calculation is performed in vacuum, the QM Mulliken charges are then
//! coupled to the MM point charges, and MM-MM point-charge interactions are
//! added. Electrostatic embedding (the MM charges entering the QM Hamiltonian)
//! is the documented follow-up.
//!
//! # Units
//!
//! Positions in Å, charges in elementary charge, energies in kcal/mol. The
//! Coulomb constant is `332.06371 kcal·Å/(mol·e²)`.

use serde::{Deserialize, Serialize};

use crate::error::SciFormError;

const HARTREE_TO_KCAL_MOL: f64 = 627.509_474_063_1;
const EV_TO_KCAL_MOL: f64 = 23.060_548_867;
const COULOMB_KCAL_ANGSTROM: f64 = 332.063_71;

/// Electronic-structure backend for the QM region.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum QmMethod {
    /// Closed-shell UHF/STO-3G. Supports charged QM regions.
    Uhf,
    /// GFN2-xTB. Neutral closed-shell QM regions only.
    Xtb,
}

/// QM/MM system definition.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QmMmSystem {
    /// QM atomic numbers.
    pub qm_elements: Vec<u8>,
    /// QM Cartesian coordinates in Å.
    pub qm_positions: Vec<[f64; 3]>,
    /// MM atomic numbers (kept for a future Lennard-Jones term).
    pub mm_elements: Vec<u8>,
    /// MM Cartesian coordinates in Å.
    pub mm_positions: Vec<[f64; 3]>,
    /// MM point charges in elementary charge.
    pub mm_charges: Vec<f64>,
}

/// QM/MM configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QmMmConfig {
    /// QM method.
    pub method: QmMethod,
    /// Total QM charge.
    pub qm_charge: i32,
    /// QM spin multiplicity.
    pub qm_multiplicity: u32,
    /// Include the MM-MM point-charge interaction.
    pub include_mm_mm: bool,
}

impl Default for QmMmConfig {
    fn default() -> Self {
        Self {
            method: QmMethod::Uhf,
            qm_charge: 0,
            qm_multiplicity: 1,
            include_mm_mm: true,
        }
    }
}

/// QM/MM energy decomposition.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QmMmResult {
    /// QM energy in kcal/mol.
    pub qm_energy_kcal_mol: f64,
    /// QM-MM point-charge interaction in kcal/mol.
    pub qm_mm_energy_kcal_mol: f64,
    /// MM-MM point-charge interaction in kcal/mol.
    pub mm_mm_energy_kcal_mol: f64,
    /// Total energy in kcal/mol.
    pub total_energy_kcal_mol: f64,
    /// QM Mulliken charges in elementary charge.
    pub qm_charges: Vec<f64>,
    /// Number of QM atoms.
    pub n_qm: usize,
    /// Number of MM atoms.
    pub n_mm: usize,
    /// Notes and caveats.
    pub notes: Vec<String>,
}

impl QmMmSystem {
    fn validate(&self) -> Result<(), SciFormError> {
        if self.qm_elements.is_empty() || self.qm_positions.is_empty() {
            return Err(SciFormError::validation(
                "the QM region must contain at least one atom",
            ));
        }
        if self.qm_elements.len() != self.qm_positions.len() {
            return Err(SciFormError::geometry(
                "qm_elements and qm_positions must have the same length",
            ));
        }
        if self.mm_elements.len() != self.mm_positions.len()
            || self.mm_positions.len() != self.mm_charges.len()
        {
            return Err(SciFormError::geometry(
                "mm_elements, mm_positions and mm_charges must have the same length",
            ));
        }
        if !self
            .qm_positions
            .iter()
            .chain(&self.mm_positions)
            .all(|position| position.iter().all(|value| value.is_finite()))
        {
            return Err(SciFormError::geometry("all coordinates must be finite"));
        }
        if self.mm_charges.iter().any(|charge| !charge.is_finite()) {
            return Err(SciFormError::validation("MM charges must be finite"));
        }
        Ok(())
    }
}

/// Compute a mechanical-embedding QM/MM energy.
pub fn compute_qmmm_energy(
    system: &QmMmSystem,
    config: &QmMmConfig,
) -> Result<QmMmResult, SciFormError> {
    system.validate()?;
    let (qm_energy_kcal, qm_charges) = match config.method {
        QmMethod::Uhf => {
            let scf = crate::compute_uhf(
                &system.qm_elements,
                &system.qm_positions,
                config.qm_charge,
                config.qm_multiplicity,
            )?;
            if !scf.converged {
                return Err(SciFormError::numerical("QM region SCF did not converge"));
            }
            (scf.total_energy * HARTREE_TO_KCAL_MOL, scf.mulliken_charges)
        }
        QmMethod::Xtb => {
            if config.qm_charge != 0 {
                return Err(SciFormError::unsupported(
                    "the xTB QM backend supports neutral regions only; use UHF for charged QM regions",
                ));
            }
            let result = crate::compute_xtb(&system.qm_elements, &system.qm_positions)?;
            (
                result.total_energy * EV_TO_KCAL_MOL,
                result.mulliken_charges,
            )
        }
    };

    let mut qm_mm = 0.0;
    for (qm_index, qm_charge) in qm_charges.iter().enumerate() {
        for (mm_index, mm_charge) in system.mm_charges.iter().enumerate() {
            let r = distance(system.qm_positions[qm_index], system.mm_positions[mm_index]);
            if r > 1e-10 {
                qm_mm += COULOMB_KCAL_ANGSTROM * qm_charge * mm_charge / r;
            }
        }
    }

    let mut mm_mm = 0.0;
    if config.include_mm_mm {
        for i in 0..system.mm_charges.len() {
            for j in (i + 1)..system.mm_charges.len() {
                let r = distance(system.mm_positions[i], system.mm_positions[j]);
                if r > 1e-10 {
                    mm_mm +=
                        COULOMB_KCAL_ANGSTROM * system.mm_charges[i] * system.mm_charges[j] / r;
                }
            }
        }
    }

    let total = qm_energy_kcal + qm_mm + mm_mm;
    Ok(QmMmResult {
        qm_energy_kcal_mol: qm_energy_kcal,
        qm_mm_energy_kcal_mol: qm_mm,
        mm_mm_energy_kcal_mol: mm_mm,
        total_energy_kcal_mol: total,
        qm_charges,
        n_qm: system.qm_elements.len(),
        n_mm: system.mm_charges.len(),
        notes: vec![
            "Mechanical embedding: the QM region is computed in vacuum and coupled through Mulliken charges."
                .to_string(),
            "MM-MM uses point-charge Coulomb only; no Lennard-Jones or bonded MM terms are included."
                .to_string(),
            "Electrostatic embedding and link atoms are follow-ups.".to_string(),
        ],
    })
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn water_system() -> QmMmSystem {
        QmMmSystem {
            qm_elements: vec![8, 1, 1],
            qm_positions: vec![
                [0.0, 0.0, 0.0],
                [0.9572, 0.0, 0.0],
                [-0.239_987, 0.927_297, 0.0],
            ],
            mm_elements: vec![8],
            mm_positions: vec![[0.0, 0.0, 5.0]],
            mm_charges: vec![-0.5],
        }
    }

    #[test]
    fn empty_mm_region_matches_qm_energy() {
        let mut system = water_system();
        system.mm_elements.clear();
        system.mm_positions.clear();
        system.mm_charges.clear();
        let result = compute_qmmm_energy(&system, &QmMmConfig::default()).unwrap();
        assert!((result.qm_mm_energy_kcal_mol).abs() < 1e-12);
        assert!((result.mm_mm_energy_kcal_mol).abs() < 1e-12);
        assert!((result.total_energy_kcal_mol - result.qm_energy_kcal_mol).abs() < 1e-12);
    }

    #[test]
    fn qm_mm_coulomb_matches_manual_sum() {
        let system = water_system();
        let result = compute_qmmm_energy(&system, &QmMmConfig::default()).unwrap();
        let manual: f64 = result
            .qm_charges
            .iter()
            .enumerate()
            .map(|(index, charge)| {
                let distance = distance(system.qm_positions[index], system.mm_positions[0]);
                COULOMB_KCAL_ANGSTROM * charge * system.mm_charges[0] / distance
            })
            .sum();
        assert!((result.qm_mm_energy_kcal_mol - manual).abs() < 1e-8);
    }

    #[test]
    fn mismatched_mm_lengths_are_rejected() {
        let mut system = water_system();
        system.mm_charges.clear();
        assert!(compute_qmmm_energy(&system, &QmMmConfig::default()).is_err());
    }
}
