//! Multi-method DOS: use orbital energies from EHT, PM3, xTB, or HF-3c.
//!
//! Provides a unified API for computing DOS from any supported electronic
//! structure method, choosing the appropriate level of theory via `DosMethod`.

use super::dos::compute_dos;
use serde::{Deserialize, Serialize};

/// Method for obtaining orbital energies for DOS computation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DosMethod {
    /// Extended Hückel Theory (fastest, qualitative).
    Eht,
    /// PM3 semi-empirical (good for organic molecules).
    Pm3,
    /// GFN-xTB tight-binding (good for organometallics).
    Xtb,
    /// GFN1-xTB (improved dispersion).
    Gfn1,
    /// GFN2-xTB (multipole electrostatics).
    Gfn2,
    /// HF-3c minimal Hartree-Fock.
    Hf3c,
}

/// Multi-method DOS result extending base DosResult.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiMethodDosResult {
    /// Energy grid (eV).
    pub energies: Vec<f64>,
    /// Total DOS (states/eV).
    pub total_dos: Vec<f64>,
    /// Method used.
    pub method: DosMethod,
    /// HOMO energy (eV).
    pub homo_energy: f64,
    /// LUMO energy (eV).
    pub lumo_energy: f64,
    /// HOMO-LUMO gap (eV).
    pub gap: f64,
    /// Orbital energies used (eV).
    pub orbital_energies: Vec<f64>,
    /// Smearing width (eV).
    pub sigma: f64,
}

/// Compute DOS using a specified electronic structure method.
///
/// Runs the chosen method to obtain orbital energies, then broadens
/// them with Gaussian smearing to produce the DOS.
pub fn compute_dos_multimethod(
    elements: &[u8],
    positions: &[[f64; 3]],
    method: DosMethod,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> Result<MultiMethodDosResult, String> {
    // Get orbital energies from the specified method
    let (orbital_energies, homo_energy, lumo_energy, gap) = match method {
        DosMethod::Eht => {
            let eht_result = crate::eht::solve_eht(elements, positions, None)?;
            let homo = eht_result.homo_energy;
            let lumo = eht_result.lumo_energy;
            (eht_result.energies, homo, lumo, eht_result.gap)
        }
        DosMethod::Pm3 => {
            let pm3 = crate::compute_pm3(elements, positions)?;
            let homo = pm3.homo_energy;
            let lumo = pm3.lumo_energy;
            (pm3.orbital_energies, homo, lumo, pm3.gap)
        }
        DosMethod::Xtb => {
            let xtb = crate::xtb::solve_xtb(elements, positions)?;
            let homo = xtb.homo_energy;
            let lumo = xtb.lumo_energy;
            (xtb.orbital_energies, homo, lumo, xtb.gap)
        }
        DosMethod::Gfn1 => {
            let gfn1 = crate::xtb::gfn1::solve_gfn1(elements, positions)?;
            let homo = gfn1.homo_energy;
            let lumo = gfn1.lumo_energy;
            (gfn1.orbital_energies, homo, lumo, gfn1.gap)
        }
        DosMethod::Gfn2 => {
            let gfn2 = crate::xtb::gfn2::solve_gfn2(elements, positions)?;
            let homo = gfn2.homo_energy;
            let lumo = gfn2.lumo_energy;
            (gfn2.orbital_energies, homo, lumo, gfn2.gap)
        }
        DosMethod::Hf3c => {
            let config = crate::hf::HfConfig::default();
            let hf = crate::hf::solve_hf3c(elements, positions, &config)?;
            // Estimate n_occ from element electron count
            let total_e: usize = elements.iter().map(|&z| z as usize).sum();
            let n_occ = total_e / 2;
            let homo = if n_occ > 0 && n_occ <= hf.orbital_energies.len() {
                hf.orbital_energies[n_occ - 1] * 27.2114 // Hartree → eV
            } else {
                0.0
            };
            let lumo = if n_occ < hf.orbital_energies.len() {
                hf.orbital_energies[n_occ] * 27.2114
            } else {
                0.0
            };
            let gap = lumo - homo;
            // Convert all orbital energies to eV
            let oe_ev: Vec<f64> = hf.orbital_energies.iter().map(|e| e * 27.2114).collect();
            (oe_ev, homo, lumo, gap)
        }
    };

    // Compute DOS with Gaussian broadening
    let dos = compute_dos(&orbital_energies, sigma, e_min, e_max, n_points);

    Ok(MultiMethodDosResult {
        energies: dos.energies,
        total_dos: dos.total_dos,
        method,
        homo_energy,
        lumo_energy,
        gap,
        orbital_energies,
        sigma,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dos_method_eht() {
        let elements = vec![8u8, 1, 1];
        let positions = vec![
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469],
        ];
        let result =
            compute_dos_multimethod(&elements, &positions, DosMethod::Eht, 0.3, -30.0, 5.0, 100);
        assert!(result.is_ok());
        let r = result.unwrap();
        assert_eq!(r.energies.len(), 100);
        assert!(r.gap > 0.0);
    }
}
