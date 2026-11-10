//! Isotropic EPR hyperfine coupling from the UHF spin density.
//!
//! This is the first EPR prototype: it evaluates the Fermi-contact spin density
//! at each selected nucleus and converts it to MHz. The g-tensor, anisotropic
//! hyperfine tensor, relativistic corrections and paramagnetic shifts are
//! documented follow-ups.
//!
//! # Units
//!
//! Spin density is in electrons per Bohr³, hyperfine coupling in MHz and
//! positions in Å. The conversion uses the proton as the reference:
//! `A(H) = 1420.405751 MHz` for one unpaired electron in a hydrogen 1s orbital.
//!
//! # Reference
//!
//! Fermi, *Z. Phys.* **1930**, 60, 320; Neese, *Chem. Rev.* **2012**, 112, 108.

use serde::{Deserialize, Serialize};

use crate::error::SciFormError;
use crate::scf::basis::BasisSet;

const ANGSTROM_TO_BOHR: f64 = 1.889_726_031_269_39;
const PROTON_G_FACTOR: f64 = 5.585_694_689_3;
const PROTON_HYPERFINE_MHZ: f64 = 1_420.405_751;

/// Configuration for an isotropic hyperfine calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EprConfig {
    /// Target atomic number. `None` evaluates every atom with a known g-factor.
    pub nucleus: Option<u8>,
    /// Total molecular charge.
    pub charge: i32,
    /// Spin multiplicity; the default is a doublet.
    pub multiplicity: u32,
}

impl Default for EprConfig {
    fn default() -> Self {
        Self {
            nucleus: None,
            charge: 0,
            multiplicity: 2,
        }
    }
}

/// One nucleus with its isotropic hyperfine coupling.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HyperfineEntry {
    /// Atom index in the input arrays.
    pub atom_index: usize,
    /// Atomic number.
    pub element: u8,
    /// Nuclear g-factor used.
    pub g_factor: f64,
    /// Fermi-contact spin density at the nucleus, electrons/Bohr³.
    pub spin_density_au: f64,
    /// Isotropic hyperfine coupling in MHz.
    pub a_iso_mhz: f64,
}

/// Result of an isotropic hyperfine calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EprHyperfineResult {
    /// One entry per selected nucleus.
    pub entries: Vec<HyperfineEntry>,
    /// Spin multiplicity used.
    pub multiplicity: u32,
    /// Notes and caveats.
    pub notes: Vec<String>,
}

/// Nuclear g-factors for common NMR/EPR-active nuclei.
pub fn nuclear_g_factor(z: u8) -> Option<f64> {
    Some(match z {
        1 => 5.585_694_689_3,
        2 => 0.857_438_230,
        3 => 5.708_3,
        5 => 1.792_4,
        6 => 1.404_274_8,
        7 => 0.403_760_7,
        8 => -0.757_516,
        9 => 5.257_736,
        11 => 1.478_4,
        13 => 1.456_6,
        14 => -1.110_58,
        15 => 2.263_20,
        16 => 0.429_1,
        17 => 0.547_9,
        19 => 0.261_0,
        23 => 1.471_0,
        25 => 1.381_3,
        26 => 0.180_9,
        27 => 1.317_0,
        29 => 1.482_0,
        31 => 1.344_0,
        33 => 0.959_0,
        35 => 0.456_0,
        47 => -0.113_0,
        53 => 1.125_0,
        55 => 1.101_0,
        78 => 1.055_7,
        82 => 1.175_0,
        80 => 0.849_0,
        _ => return None,
    })
}

/// Compute isotropic hyperfine coupling from the UHF spin density.
pub fn compute_epr_hyperfine(
    elements: &[u8],
    positions_angstrom: &[[f64; 3]],
    config: &EprConfig,
) -> Result<EprHyperfineResult, SciFormError> {
    crate::validation::validate_geometry(elements, positions_angstrom)?;
    if config.multiplicity == 0 {
        return Err(SciFormError::validation("multiplicity must be >= 1"));
    }
    let scf = crate::compute_uhf(
        elements,
        positions_angstrom,
        config.charge,
        config.multiplicity,
    )?;
    if !scf.converged {
        return Err(SciFormError::numerical(
            "UHF did not converge for the hyperfine calculation",
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
    let spin_density = &scf.alpha_density - &scf.beta_density;
    let mut entries = Vec::new();
    let mut skipped = Vec::new();

    for (atom_index, &element) in elements.iter().enumerate() {
        if let Some(target) = config.nucleus {
            if element != target {
                continue;
            }
        }
        let Some(g_factor) = nuclear_g_factor(element) else {
            skipped.push(element);
            continue;
        };
        let center = positions_bohr[atom_index];
        let mut density = 0.0;
        for i in 0..basis.n_basis {
            let value_i = basis.functions[i].evaluate(center[0], center[1], center[2]);
            if value_i == 0.0 {
                continue;
            }
            for j in 0..basis.n_basis {
                let value_j = basis.functions[j].evaluate(center[0], center[1], center[2]);
                if value_j != 0.0 {
                    density += spin_density[(i, j)] * value_i * value_j;
                }
            }
        }
        let a_iso = PROTON_HYPERFINE_MHZ
            * (g_factor / PROTON_G_FACTOR)
            * (density / proton_reference_density());
        entries.push(HyperfineEntry {
            atom_index,
            element,
            g_factor,
            spin_density_au: density,
            a_iso_mhz: a_iso,
        });
    }

    let mut notes = vec![
        "Isotropic Fermi-contact term only; the anisotropic tensor and g-tensor are not implemented."
            .to_string(),
        "The conversion uses the proton 1s reference and scales by the nuclear g-factor."
            .to_string(),
        "Spin polarization is taken from the UHF/STO-3G density difference.".to_string(),
    ];
    if !skipped.is_empty() {
        notes.push(format!(
            "No nuclear g-factor table entry for element(s) {skipped:?}; they were skipped."
        ));
    }
    if entries.is_empty() {
        return Err(SciFormError::unsupported(
            "no selected nucleus has a known g-factor",
        ));
    }

    Ok(EprHyperfineResult {
        entries,
        multiplicity: config.multiplicity,
        notes,
    })
}

fn proton_reference_density() -> f64 {
    let basis = BasisSet::sto3g(&[1], &[[0.0, 0.0, 0.0]]);
    basis
        .functions
        .iter()
        .map(|function| function.evaluate(0.0, 0.0, 0.0))
        .sum::<f64>()
        .powi(2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hydrogen_atom_hyperfine_is_near_1420_mhz() {
        let result = compute_epr_hyperfine(
            &[1],
            &[[0.0, 0.0, 0.0]],
            &EprConfig {
                nucleus: Some(1),
                multiplicity: 2,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(result.entries.len(), 1);
        let value = result.entries[0].a_iso_mhz;
        assert!((value - 1_420.405_751).abs() < 5.0, "A_iso = {value} MHz");
    }

    #[test]
    fn closed_shell_water_has_near_zero_spin_density() {
        let water = [
            [0.0, 0.0, 0.0],
            [0.9572, 0.0, 0.0],
            [-0.239_987, 0.927_297, 0.0],
        ];
        let result = compute_epr_hyperfine(
            &[8, 1, 1],
            &water,
            &EprConfig {
                nucleus: Some(1),
                multiplicity: 1,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(result
            .entries
            .iter()
            .all(|entry| entry.a_iso_mhz.abs() < 1.0));
    }

    #[test]
    fn unknown_g_factor_is_reported() {
        let result = compute_epr_hyperfine(
            &[1],
            &[[0.0, 0.0, 0.0]],
            &EprConfig {
                nucleus: Some(99),
                multiplicity: 2,
                ..Default::default()
            },
        );
        assert!(result.is_err());
    }
}
