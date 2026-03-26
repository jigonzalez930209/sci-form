//! Public API for HF-3c composite method.
//!
//! Ties together: SCF (Hartree-Fock) + D3 + gCP + SRB corrections,
//! plus optional CIS excited-state calculation for UV-Vis spectroscopy.

use super::basis::{build_sto3g_basis, ANG_TO_BOHR};
use super::cis::{compute_cis_with_dipole, CisResult};
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
            max_scf_iter: 300,
            diis_size: 6,
            n_cis_states: 5,
            corrections: true,
        }
    }
}

/// Result of an HF-3c calculation.
///
/// Energy breakdown: `energy = hf_energy + d3_energy + gcp_energy + srb_energy`
/// where `hf_energy = electronic + nuclear_repulsion`.
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
    /// Orbital energies (sorted, eV).
    pub orbital_energies: Vec<f64>,
    /// Number of SCF iterations.
    pub scf_iterations: usize,
    /// Whether SCF converged.
    pub converged: bool,
    /// CIS excitation results (if requested).
    pub cis: Option<CisResult>,
    /// Number of basis functions.
    pub n_basis: usize,
    /// Number of electrons.
    pub n_electrons: usize,
    /// HOMO energy (eV).
    pub homo_energy: f64,
    /// LUMO energy (eV), if available.
    pub lumo_energy: Option<f64>,
    /// HOMO–LUMO gap (eV).
    pub gap: f64,
    /// Mulliken charges per atom.
    pub mulliken_charges: Vec<f64>,
}

#[cfg(feature = "experimental-gpu")]
fn hf_basis_to_gpu_basis(basis: &super::basis::BasisSet) -> crate::scf::basis::BasisSet {
    use crate::scf::basis::{
        BasisFunction as GpuBasisFunction, BasisSet as GpuBasisSet,
        ContractedShell as GpuContractedShell, GaussianPrimitive,
    };

    let mut functions = Vec::new();
    let mut shells = Vec::new();
    let mut function_to_atom = Vec::new();

    for shell in &basis.shells {
        let primitives: Vec<GaussianPrimitive> = shell
            .exponents
            .iter()
            .zip(shell.coefficients.iter())
            .map(|(&alpha, &coefficient)| GaussianPrimitive { alpha, coefficient })
            .collect();

        let l = match shell.shell_type {
            super::basis::ShellType::S => 0,
            super::basis::ShellType::P => 1,
        };

        shells.push(GpuContractedShell {
            atom_index: shell.center_idx,
            center: shell.center,
            l,
            primitives: primitives.clone(),
        });

        match shell.shell_type {
            super::basis::ShellType::S => {
                functions.push(GpuBasisFunction {
                    atom_index: shell.center_idx,
                    center: shell.center,
                    angular: [0, 0, 0],
                    l_total: 0,
                    primitives: primitives.clone(),
                });
                function_to_atom.push(shell.center_idx);
            }
            super::basis::ShellType::P => {
                for angular in [[1, 0, 0], [0, 1, 0], [0, 0, 1]] {
                    functions.push(GpuBasisFunction {
                        atom_index: shell.center_idx,
                        center: shell.center,
                        angular,
                        l_total: 1,
                        primitives: primitives.clone(),
                    });
                    function_to_atom.push(shell.center_idx);
                }
            }
        }
    }

    let n_basis = functions.len();
    GpuBasisSet {
        functions,
        shells,
        n_basis,
        function_to_atom,
    }
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

    #[cfg(feature = "experimental-gpu")]
    let gpu_eris_full = if n_basis >= 4 {
        // Memory check: full N⁴ tensor requires n_basis⁴ × 8 bytes.
        // Cap at ~512 MB to prevent OOM.
        let n4 = (n_basis as u64)
            .saturating_mul(n_basis as u64)
            .saturating_mul(n_basis as u64)
            .saturating_mul(n_basis as u64);
        let mem_bytes = n4.saturating_mul(8);
        let max_mem: u64 = 512 * 1024 * 1024; // 512 MB

        if mem_bytes > max_mem {
            None // Too large for dense tensor; fall back to CPU packed ERIs
        } else if let Ok(ctx) = crate::gpu::context::GpuContext::try_create() {
            let gpu_basis = hf_basis_to_gpu_basis(&basis);
            crate::gpu::two_electron_gpu::compute_eris_gpu(&ctx, &gpu_basis)
                .ok()
                .map(|gpu_eris| {
                    let cap = n_basis * n_basis * n_basis * n_basis;
                    let mut full = Vec::with_capacity(cap);
                    for mu in 0..n_basis {
                        for nu in 0..n_basis {
                            for lam in 0..n_basis {
                                for sig in 0..n_basis {
                                    full.push(gpu_eris.get(mu, nu, lam, sig));
                                }
                            }
                        }
                    }
                    full
                })
        } else {
            None
        }
    } else {
        None
    };

    #[cfg(not(feature = "experimental-gpu"))]
    let gpu_eris_full: Option<Vec<f64>> = None;

    // SCF
    let scf_config = ScfConfig {
        max_iter: config.max_scf_iter,
        diis_size: config.diis_size,
        ..ScfConfig::default()
    };
    let scf_result = solve_scf(
        &h_core,
        &s_mat,
        &eris,
        gpu_eris_full.as_deref(),
        n_electrons,
        &scf_config,
    );

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
        let ao_map = super::basis::ao_to_atom_map(&basis);
        Some(compute_cis_with_dipole(
            &scf_result.orbital_energies,
            &scf_result.coefficients,
            &eris,
            n_basis,
            n_occ,
            config.n_cis_states,
            Some(&pos_bohr),
            Some(&ao_map),
        ))
    } else {
        None
    };

    // Extract HOMO/LUMO from orbital energies
    let n_occ = n_electrons / 2;
    let homo_energy = if n_occ > 0 && n_occ <= scf_result.orbital_energies.len() {
        scf_result.orbital_energies[n_occ - 1]
    } else {
        0.0
    };
    let lumo_energy = if n_occ < scf_result.orbital_energies.len() {
        Some(scf_result.orbital_energies[n_occ])
    } else {
        None
    };
    let gap = lumo_energy.map_or(0.0, |l| l - homo_energy);

    // Mulliken charges from converged density
    let mulliken_charges = if scf_result.converged {
        let ps = &scf_result.density * &s_mat;
        let ao_to_atom = super::basis::ao_to_atom_map(&basis);
        let mut charges = vec![0.0_f64; elements.len()];
        for mu in 0..n_basis {
            charges[ao_to_atom[mu]] += ps[(mu, mu)];
        }
        charges
            .iter()
            .enumerate()
            .map(|(i, &pop)| elements[i] as f64 - pop)
            .collect()
    } else {
        vec![0.0; elements.len()]
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
        n_basis,
        n_electrons,
        homo_energy,
        lumo_energy,
        gap,
        mulliken_charges,
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
    fn test_default_hf3c_iteration_budget() {
        assert_eq!(HfConfig::default().max_scf_iter, 300);
    }

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
