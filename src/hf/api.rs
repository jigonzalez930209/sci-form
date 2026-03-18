//! Public API for HF-3c composite method.
//!
//! Ties together: SCF (Hartree-Fock) + D3 + gCP + SRB corrections,
//! plus optional CIS excited-state calculation for UV-Vis spectroscopy.

use super::basis::{build_sto3g_basis, ANG_TO_BOHR};
use super::cis::{compute_cis, CisResult};
use super::d3::compute_d3_energy;
use super::fock::nuclear_repulsion;
use super::gcp::compute_gcp;
use super::integrals::compute_eris;
use super::nuclear::compute_nuclear_matrix;
use super::overlap_kin::{compute_kinetic_matrix, compute_overlap_matrix};
use super::scf::{solve_scf, ScfConfig};
use super::srb::compute_srb;
use serde::{Deserialize, Serialize};

/// Configuration for HF-3c calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HfConfig {
    /// Maximum SCF iterations.
    pub max_scf_iter: usize,
    /// DIIS subspace size.
    pub diis_size: usize,
    /// Number of CIS excited states to compute (0 = skip CIS).
    pub n_cis_states: usize,
    /// Include empirical corrections (D3, gCP, SRB).
    pub corrections: bool,
}

impl Default for HfConfig {
    fn default() -> Self {
        HfConfig {
            max_scf_iter: 100,
            diis_size: 6,
            n_cis_states: 5,
            corrections: true,
        }
    }
}

/// Result of an HF-3c calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Hf3cResult {
    /// Total HF-3c energy (Hartree).
    pub energy: f64,
    /// Pure HF electronic energy.
    pub hf_energy: f64,
    /// Nuclear repulsion energy.
    pub nuclear_repulsion: f64,
    /// D3 dispersion correction energy.
    pub d3_energy: f64,
    /// gCP BSSE correction energy.
    pub gcp_energy: f64,
    /// SRB short-range correction energy.
    pub srb_energy: f64,
    /// Orbital energies (sorted).
    pub orbital_energies: Vec<f64>,
    /// Number of SCF iterations.
    pub scf_iterations: usize,
    /// Whether SCF converged.
    pub converged: bool,
    /// CIS excitation results (if requested).
    pub cis: Option<CisResult>,
}

/// Run a complete HF-3c calculation.
pub fn solve_hf3c(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &HfConfig,
) -> Result<Hf3cResult, String> {
    if elements.len() != positions.len() {
        return Err("elements/positions length mismatch".to_string());
    }
    if elements.is_empty() {
        return Err("empty molecule".to_string());
    }

    // Convert positions to Bohr
    let pos_bohr: Vec<[f64; 3]> = positions
        .iter()
        .map(|p| [p[0] * ANG_TO_BOHR, p[1] * ANG_TO_BOHR, p[2] * ANG_TO_BOHR])
        .collect();

    // Build basis set
    let basis = build_sto3g_basis(elements, positions);
    let n_basis = basis.n_basis();

    // Count electrons
    let n_electrons: usize = elements.iter().map(|&z| z as usize).sum();

    // One-electron integrals
    let s_mat = compute_overlap_matrix(&basis);
    let t_mat = compute_kinetic_matrix(&basis);
    let v_mat = compute_nuclear_matrix(&basis, elements, &pos_bohr);
    let h_core = &t_mat + &v_mat;

    // Two-electron integrals
    let eris = compute_eris(&basis);

    // SCF
    let scf_config = ScfConfig {
        max_iter: config.max_scf_iter,
        diis_size: config.diis_size,
        ..ScfConfig::default()
    };
    let scf_result = solve_scf(&h_core, &s_mat, &eris, n_electrons, &scf_config);

    // Nuclear repulsion
    let e_nuc = nuclear_repulsion(elements, &pos_bohr);

    // Empirical corrections
    let (d3_e, gcp_e, srb_e) = if config.corrections {
        (
            compute_d3_energy(elements, &pos_bohr).energy,
            compute_gcp(elements, &pos_bohr),
            compute_srb(elements, &pos_bohr),
        )
    } else {
        (0.0, 0.0, 0.0)
    };

    let total = scf_result.energy + e_nuc + d3_e + gcp_e + srb_e;

    // CIS excited states
    let cis = if config.n_cis_states > 0 && scf_result.converged {
        let n_occ = n_electrons / 2;
        Some(compute_cis(
            &scf_result.orbital_energies,
            &scf_result.coefficients,
            &eris,
            n_basis,
            n_occ,
            config.n_cis_states,
        ))
    } else {
        None
    };

    Ok(Hf3cResult {
        energy: total,
        hf_energy: scf_result.energy + e_nuc,
        nuclear_repulsion: e_nuc,
        d3_energy: d3_e,
        gcp_energy: gcp_e,
        srb_energy: srb_e,
        orbital_energies: scf_result.orbital_energies,
        scf_iterations: scf_result.iterations,
        converged: scf_result.converged,
        cis,
    })
}

/// Batch-process multiple HF-3c calculations in parallel.
#[cfg(feature = "parallel")]
pub fn solve_hf3c_batch(
    molecules: &[(&[u8], &[[f64; 3]])],
    config: &HfConfig,
) -> Vec<Result<Hf3cResult, String>> {
    use rayon::prelude::*;
    molecules
        .par_iter()
        .map(|(els, pos)| solve_hf3c(els, pos, config))
        .collect()
}

/// Batch-process multiple HF-3c calculations sequentially.
#[cfg(not(feature = "parallel"))]
pub fn solve_hf3c_batch(
    molecules: &[(&[u8], &[[f64; 3]])],
    config: &HfConfig,
) -> Vec<Result<Hf3cResult, String>> {
    molecules
        .iter()
        .map(|(els, pos)| solve_hf3c(els, pos, config))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_h2_hf3c() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]];
        let config = HfConfig {
            n_cis_states: 0,
            ..Default::default()
        };
        let result = solve_hf3c(&elements, &positions, &config).unwrap();
        assert!(result.energy.is_finite(), "Energy should be finite");
        assert!(result.energy < 0.0, "H2 total energy should be negative");
    }

    #[test]
    fn test_water_hf3c() {
        let elements = [8u8, 1, 1];
        let positions = [
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469],
        ];
        let result = solve_hf3c(&elements, &positions, &HfConfig::default()).unwrap();
        assert!(result.energy.is_finite());
        assert!(result.orbital_energies.len() == 7); // 7 basis functions
    }
}
