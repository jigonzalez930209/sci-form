//! Core-level XPS binding energies in the Koopmans approximation.
//!
//! The first v2.3 core-spectroscopy prototype identifies the most localized
//! core-like s orbital on each atom and reports its UHF/STO-3G orbital energy as
//! a binding-energy estimate. This is a Koopmans screening value, not a ΔSCF
//! core-hole calculation and not XANES.
//!
//! # Units
//!
//! Orbital energies in Hartree, binding energies in eV, positions in Å.

use serde::{Deserialize, Serialize};

use crate::error::SciFormError;
use crate::scf::basis::BasisSet;

const ANGSTROM_TO_BOHR: f64 = 1.889_726_031_269_39;
const HARTREE_TO_EV: f64 = 27.211_386_245_988;

/// Configuration for a Koopmans XPS estimate.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XpsConfig {
    /// Total molecular charge.
    pub charge: i32,
    /// Spin multiplicity.
    pub multiplicity: u32,
    /// Restrict the output to these elements. `None` reports every atom.
    pub elements: Option<Vec<u8>>,
}

impl Default for XpsConfig {
    fn default() -> Self {
        Self {
            charge: 0,
            multiplicity: 1,
            elements: None,
        }
    }
}

/// One core-level binding estimate.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XpsBinding {
    /// Atom index.
    pub atom_index: usize,
    /// Atomic number.
    pub element: u8,
    /// Orbital energy used, Hartree.
    pub orbital_energy_hartree: f64,
    /// Koopmans binding energy, eV.
    pub binding_energy_ev: f64,
    /// Fraction of the MO localized on the atom's s functions.
    pub localization: f64,
}

/// Result of a Koopmans XPS calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XpsResult {
    /// Binding estimates.
    pub bindings: Vec<XpsBinding>,
    /// Notes and caveats.
    pub notes: Vec<String>,
}

/// Compute Koopmans core binding energies for every selected atom.
pub fn compute_xps_binding_energies(
    elements: &[u8],
    positions_angstrom: &[[f64; 3]],
    config: &XpsConfig,
) -> Result<XpsResult, SciFormError> {
    crate::validation::validate_geometry(elements, positions_angstrom)?;
    let scf = crate::compute_uhf(
        elements,
        positions_angstrom,
        config.charge,
        config.multiplicity.max(1),
    )?;
    if !scf.converged {
        return Err(SciFormError::numerical(
            "UHF did not converge for the XPS calculation",
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
    let n_basis = basis.n_basis;
    let n_mo = scf.alpha_orbital_energies.len();
    let mut bindings = Vec::new();

    for (atom_index, &element) in elements.iter().enumerate() {
        if let Some(allowed) = &config.elements {
            if !allowed.contains(&element) {
                continue;
            }
        }
        let s_functions: Vec<usize> = (0..n_basis)
            .filter(|&index| {
                basis.functions[index].atom_index == atom_index
                    && basis.functions[index].l_total == 0
            })
            .collect();
        if s_functions.is_empty() {
            continue;
        }
        let mut best_mo = None;
        let mut best_weight: f64 = 0.0;
        for mo in 0..n_mo.min(scf.n_alpha.max(1)) {
            let weight: f64 = s_functions
                .iter()
                .map(|&mu| scf.alpha_coefficients[(mu, mo)] * scf.alpha_coefficients[(mu, mo)])
                .sum();
            if weight > best_weight {
                best_weight = weight;
                best_mo = Some(mo);
            }
        }
        let Some(mo) = best_mo else {
            continue;
        };
        let energy = scf.alpha_orbital_energies[mo];
        bindings.push(XpsBinding {
            atom_index,
            element,
            orbital_energy_hartree: energy,
            binding_energy_ev: -energy * HARTREE_TO_EV,
            localization: best_weight,
        });
    }

    if bindings.is_empty() {
        return Err(SciFormError::unsupported(
            "no selected atom has a core-like s orbital in the STO-3G basis",
        ));
    }

    Ok(XpsResult {
        bindings,
        notes: vec![
            "Koopmans approximation: binding energy = -orbital energy, no core-hole relaxation."
                .to_string(),
            "UHF/STO-3G only; absolute values are screening estimates, not benchmarked XPS."
                .to_string(),
            "XANES and transition-potential core excitations are not implemented.".to_string(),
            "The selected orbital is the most atom-localized s orbital, not a symmetry-labelled core state."
                .to_string(),
        ],
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn methane_carbon_binding_is_in_the_core_range() {
        let methane = [
            [0.0, 0.0, 0.0],
            [0.629, 0.629, 0.629],
            [0.629, -0.629, -0.629],
            [-0.629, 0.629, -0.629],
            [-0.629, -0.629, 0.629],
        ];
        let result = compute_xps_binding_energies(
            &[6, 1, 1, 1, 1],
            &methane,
            &XpsConfig {
                elements: Some(vec![6]),
                ..XpsConfig::default()
            },
        )
        .unwrap();
        assert_eq!(result.bindings.len(), 1);
        let binding = result.bindings[0].binding_energy_ev;
        assert!(
            (200.0..400.0).contains(&binding),
            "C 1s binding = {binding} eV"
        );
    }

    #[test]
    fn invalid_geometry_is_rejected() {
        assert!(compute_xps_binding_energies(
            &[6],
            &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            &XpsConfig::default()
        )
        .is_err());
    }
}
