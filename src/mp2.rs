//! Conventional closed-shell MP2 correlation energy.
//!
//! The first v2.3 correlation prototype transforms the existing STO-3G ERIs to
//! the UHF/RHF orbital basis and evaluates the closed-shell MP2 expression. It
//! is a conventional four-index MP2, not the RI/local variant: the RI
//! approximation and local PNO correlation are follow-ups.
//!
//! # Units
//!
//! Orbital energies and correlation energy are in Hartree; the result also
//! reports kcal/mol. Positions are Å.

use serde::{Deserialize, Serialize};

use crate::error::SciFormError;
use crate::scf::basis::BasisSet;
use crate::scf::two_electron::TwoElectronIntegrals;

const ANGSTROM_TO_BOHR: f64 = 1.889_726_031_269_39;
const HARTREE_TO_KCAL_MOL: f64 = 627.509_474_063_1;

/// Result of a conventional MP2 calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Mp2Result {
    /// HF total energy in Hartree.
    pub hf_energy_hartree: f64,
    /// MP2 correlation energy in Hartree.
    pub correlation_energy_hartree: f64,
    /// MP2 correlation energy in kcal/mol.
    pub correlation_energy_kcal_mol: f64,
    /// Total HF + MP2 energy in Hartree.
    pub total_energy_hartree: f64,
    /// Number of basis functions.
    pub n_basis: usize,
    /// Number of occupied orbitals.
    pub n_occupied: usize,
    /// Number of virtual orbitals.
    pub n_virtual: usize,
    /// Notes and caveats.
    pub notes: Vec<String>,
}

/// Compute the conventional closed-shell MP2 correlation energy.
///
/// The implementation uses the RHF/STO-3G orbitals from the existing UHF path
/// (singlet, neutral). The four-index transform is performed one index at a
/// time, so the memory is `O(N² N_occ N_vir)` rather than `O(N⁴)`.
pub fn compute_mp2_energy(
    elements: &[u8],
    positions_angstrom: &[[f64; 3]],
) -> Result<Mp2Result, SciFormError> {
    crate::validation::validate_geometry(elements, positions_angstrom)?;
    let scf = crate::compute_uhf(elements, positions_angstrom, 0, 1)?;
    if !scf.converged {
        return Err(SciFormError::numerical(
            "UHF did not converge for the MP2 calculation",
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
    let n = basis.n_basis;
    let n_occ = scf.n_alpha;
    let n_vir = n.saturating_sub(n_occ);
    if n_occ == 0 || n_vir == 0 {
        return Err(SciFormError::unsupported(
            "MP2 requires both occupied and virtual orbitals",
        ));
    }
    let eri = TwoElectronIntegrals::compute(&basis);
    let c = &scf.alpha_coefficients;

    // eri_1[mu, nu, lambda, b] = sum_sigma (mu nu | lambda sigma) C[sigma, b]
    let mut eri_1 = vec![0.0f64; n * n * n * n_vir];
    for mu in 0..n {
        for nu in 0..n {
            for lambda in 0..n {
                for b in 0..n_vir {
                    let mut value = 0.0;
                    for sigma in 0..n {
                        value += eri.get(mu, nu, lambda, sigma) * c[(sigma, n_occ + b)];
                    }
                    eri_1[((mu * n + nu) * n + lambda) * n_vir + b] = value;
                }
            }
        }
    }

    // eri_2[mu, nu, j, b] = sum_lambda eri_1 C[lambda, j]
    let mut eri_2 = vec![0.0f64; n * n * n_occ * n_vir];
    for mu in 0..n {
        for nu in 0..n {
            for j in 0..n_occ {
                for b in 0..n_vir {
                    let mut value = 0.0;
                    for lambda in 0..n {
                        value += eri_1[((mu * n + nu) * n + lambda) * n_vir + b] * c[(lambda, j)];
                    }
                    eri_2[((mu * n + nu) * n_occ + j) * n_vir + b] = value;
                }
            }
        }
    }

    // eri_3[mu, a, j, b] = sum_nu eri_2 C[nu, a]
    let mut eri_3 = vec![0.0f64; n * n_vir * n_occ * n_vir];
    for mu in 0..n {
        for a in 0..n_vir {
            for j in 0..n_occ {
                for b in 0..n_vir {
                    let mut value = 0.0;
                    for nu in 0..n {
                        value +=
                            eri_2[((mu * n + nu) * n_occ + j) * n_vir + b] * c[(nu, n_occ + a)];
                    }
                    eri_3[((mu * n_vir + a) * n_occ + j) * n_vir + b] = value;
                }
            }
        }
    }

    // eri_4[i, a, j, b] = sum_mu eri_3 C[mu, i]
    let mut eri_4 = vec![0.0f64; n_occ * n_vir * n_occ * n_vir];
    for i in 0..n_occ {
        for a in 0..n_vir {
            for j in 0..n_occ {
                for b in 0..n_vir {
                    let mut value = 0.0;
                    for mu in 0..n {
                        value += eri_3[((mu * n_vir + a) * n_occ + j) * n_vir + b] * c[(mu, i)];
                    }
                    eri_4[((i * n_vir + a) * n_occ + j) * n_vir + b] = value;
                }
            }
        }
    }

    let mut correlation = 0.0;
    for i in 0..n_occ {
        for j in 0..n_occ {
            for a in 0..n_vir {
                for b in 0..n_vir {
                    let iajb = eri_4[((i * n_vir + a) * n_occ + j) * n_vir + b];
                    let ibja = eri_4[((i * n_vir + b) * n_occ + j) * n_vir + a];
                    let denominator = scf.alpha_orbital_energies[i] + scf.alpha_orbital_energies[j]
                        - scf.alpha_orbital_energies[n_occ + a]
                        - scf.alpha_orbital_energies[n_occ + b];
                    if denominator.abs() > 1e-12 {
                        correlation += (2.0 * iajb - ibja) * iajb / denominator;
                    }
                }
            }
        }
    }

    if !correlation.is_finite() {
        return Err(SciFormError::numerical(
            "MP2 correlation energy is not finite",
        ));
    }

    Ok(Mp2Result {
        hf_energy_hartree: scf.total_energy,
        correlation_energy_hartree: correlation,
        correlation_energy_kcal_mol: correlation * HARTREE_TO_KCAL_MOL,
        total_energy_hartree: scf.total_energy + correlation,
        n_basis: n,
        n_occupied: n_occ,
        n_virtual: n_vir,
        notes: vec![
            "Conventional four-index MP2, closed-shell RHF/STO-3G.".to_string(),
            "The RI approximation and local PNO correlation are not implemented yet.".to_string(),
            "No frozen-core option; all occupied orbitals are correlated.".to_string(),
        ],
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn water_mp2_correlation_is_negative_and_finite() {
        let water = [
            [0.0, 0.0, 0.0],
            [0.9572, 0.0, 0.0],
            [-0.239_987, 0.927_297, 0.0],
        ];
        let result = compute_mp2_energy(&[8, 1, 1], &water).unwrap();
        assert!(result.correlation_energy_hartree < 0.0);
        assert!(result.correlation_energy_hartree.is_finite());
        assert!(result.n_occupied > 0 && result.n_virtual > 0);
        assert!(
            (result.total_energy_hartree
                - (result.hf_energy_hartree + result.correlation_energy_hartree))
                .abs()
                < 1e-12
        );
    }

    #[test]
    fn invalid_geometry_is_rejected() {
        assert!(compute_mp2_energy(&[8], &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]).is_err());
    }
}
