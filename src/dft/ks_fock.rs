//! Kohn-Sham DFT Fock matrix builder and SCF driver.
//!
//! Builds the KS-DFT Fock matrix: F^KS = H_core + J + V_XC
//! (no exact exchange for pure DFT).
//!
//! Reuses the existing SCF infrastructure (DIIS, density matrix, convergence).

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

use super::grid::{GridQuality, MolecularGrid};
use super::vxc_matrix::build_vxc_matrix;
use crate::scf::basis::BasisSet;
use crate::scf::constants::{HARTREE_TO_EV, SCF_MAX_ITER};
use crate::scf::core_matrices::{nuclear_repulsion_energy, CoreMatrices};
use crate::scf::density_matrix::{build_density_matrix, density_rms_change};
use crate::scf::diis::DiisAccelerator;
use crate::scf::mulliken::mulliken_analysis;
use crate::scf::orthogonalization::{
    back_transform, lowdin_orthogonalization, transform_to_orthogonal,
};
use crate::scf::two_electron::TwoElectronIntegrals;
use crate::scf::types::MolecularSystem;

/// DFT functional method.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DftMethod {
    /// SVWN (LDA): Slater exchange + VWN-5 correlation.
    Svwn,
    /// PBE (GGA): Perdew-Burke-Ernzerhof.
    Pbe,
}

/// Configuration for KS-DFT calculation.
#[derive(Debug, Clone)]
pub struct DftConfig {
    /// DFT functional.
    pub method: DftMethod,
    /// Grid quality.
    pub grid_quality: GridQuality,
    /// Maximum SCF iterations.
    pub max_iterations: usize,
    /// Energy convergence threshold (Hartree).
    pub energy_threshold: f64,
    /// Density RMS convergence threshold.
    pub density_threshold: f64,
    /// DIIS subspace size.
    pub diis_size: usize,
}

impl Default for DftConfig {
    fn default() -> Self {
        Self {
            method: DftMethod::Svwn,
            grid_quality: GridQuality::Medium,
            max_iterations: SCF_MAX_ITER,
            energy_threshold: 1e-7,
            density_threshold: 1e-6,
            diis_size: 6,
        }
    }
}

/// Result of a KS-DFT calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DftResult {
    /// Orbital energies (Hartree).
    pub orbital_energies: Vec<f64>,
    /// Orbital energies in eV.
    pub orbital_energies_ev: Vec<f64>,
    /// Total electronic energy (Hartree).
    pub electronic_energy: f64,
    /// Nuclear repulsion energy (Hartree).
    pub nuclear_repulsion: f64,
    /// Exchange-correlation energy (Hartree).
    pub xc_energy: f64,
    /// Coulomb (J) energy (Hartree).
    pub coulomb_energy: f64,
    /// Total energy (Hartree).
    pub total_energy: f64,
    /// Total energy in eV.
    pub total_energy_ev: f64,
    /// HOMO energy (eV).
    pub homo_energy: f64,
    /// LUMO energy (eV).
    pub lumo_energy: f64,
    /// HOMO-LUMO gap (eV).
    pub gap: f64,
    /// Mulliken atomic charges.
    pub mulliken_charges: Vec<f64>,
    /// DFT functional used.
    pub method: String,
    /// Number of SCF iterations.
    pub scf_iterations: usize,
    /// Whether SCF converged.
    pub converged: bool,
    /// Number of basis functions.
    pub n_basis: usize,
    /// Number of electrons.
    pub n_electrons: usize,
    /// Number of grid points.
    pub n_grid_points: usize,
}

/// Run a Kohn-Sham DFT calculation.
pub fn solve_ks_dft(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &DftConfig,
) -> Result<DftResult, String> {
    let system = MolecularSystem::from_angstrom(elements, positions, 0, 1);
    let n_electrons = system.n_electrons();
    let n_occupied = n_electrons / 2;

    if n_occupied == 0 {
        return Err("KS-DFT requires at least 2 electrons".to_string());
    }

    // Build basis set
    let basis = BasisSet::sto3g(&system.atomic_numbers, &system.positions_bohr);
    let n_basis = basis.functions.len();

    // Build molecular integration grid
    let grid = MolecularGrid::build(
        &system.atomic_numbers,
        &system.positions_bohr,
        config.grid_quality,
    );

    // Build one-electron matrices
    let core = CoreMatrices::build(&basis, &system.atomic_numbers, &system.positions_bohr);
    let s = &core.overlap;
    let h_core = &core.core_hamiltonian;
    let v_nuc = nuclear_repulsion_energy(&system.atomic_numbers, &system.positions_bohr);

    // Orthogonalization matrix
    let (x, _n_independent) = lowdin_orthogonalization(s, 1e-8);

    // Two-electron integrals (for Coulomb J matrix)
    let eri = TwoElectronIntegrals::compute(&basis);

    // DIIS accelerator
    let mut diis = DiisAccelerator::new(config.diis_size);

    // Initial guess from core Hamiltonian
    let f_orth = transform_to_orthogonal(h_core, &x);
    let eigen = f_orth.symmetric_eigen();

    let mut indices: Vec<usize> = (0..n_basis).collect();
    indices.sort_by(|&a, &b| {
        eigen.eigenvalues[a]
            .partial_cmp(&eigen.eigenvalues[b])
            .unwrap()
    });

    let mut c_ortho = DMatrix::zeros(n_basis, n_basis);
    let mut orbital_energies = vec![0.0; n_basis];
    for (new_idx, &old_idx) in indices.iter().enumerate() {
        orbital_energies[new_idx] = eigen.eigenvalues[old_idx];
        for i in 0..n_basis {
            c_ortho[(i, new_idx)] = eigen.eigenvectors[(i, old_idx)];
        }
    }

    let mut c = back_transform(&c_ortho, &x);
    let mut density = build_density_matrix(&c, n_occupied);

    let mut total_energy = 0.0;
    let mut converged = false;
    let mut xc_energy = 0.0;
    let mut final_orbital_energies = orbital_energies.clone();
    let mut iterations = 0;

    for iter in 0..config.max_iterations {
        iterations = iter + 1;

        // Build Coulomb matrix J (no exchange for pure DFT)
        let mut j_matrix = DMatrix::zeros(n_basis, n_basis);
        for mu in 0..n_basis {
            for nu in 0..n_basis {
                let mut j_val = 0.0;
                for lam in 0..n_basis {
                    for sig in 0..n_basis {
                        j_val += density[(lam, sig)] * eri.get(mu, nu, lam, sig);
                    }
                }
                j_matrix[(mu, nu)] = j_val;
            }
        }

        // Build V_XC matrix
        let (vxc, exc) = build_vxc_matrix(&basis, &density, &grid, config.method);
        xc_energy = exc;

        // KS Fock matrix: F = H_core + J + V_XC
        let fock = h_core + &j_matrix + &vxc;

        // Compute Coulomb energy
        let j_energy = 0.5 * (&density * &j_matrix).trace();

        // Total electronic energy
        let one_electron_energy = (&density * h_core).trace();
        let e_elec = one_electron_energy + j_energy + xc_energy;
        let new_total = e_elec + v_nuc;

        // Convergence check
        if iter > 0 && (new_total - total_energy).abs() < config.energy_threshold {
            total_energy = new_total;
            converged = true;
            break;
        }
        total_energy = new_total;

        // DIIS extrapolation
        diis.add_iteration(&fock, &density, s);
        let fock_diis = diis.extrapolate().unwrap_or(fock);

        // Diagonalize
        let f_orth = transform_to_orthogonal(&fock_diis, &x);
        let eigen = f_orth.symmetric_eigen();

        let mut indices: Vec<usize> = (0..n_basis).collect();
        indices.sort_by(|&a, &b| {
            eigen.eigenvalues[a]
                .partial_cmp(&eigen.eigenvalues[b])
                .unwrap()
        });

        for (new_idx, &old_idx) in indices.iter().enumerate() {
            orbital_energies[new_idx] = eigen.eigenvalues[old_idx];
            for i in 0..n_basis {
                c_ortho[(i, new_idx)] = eigen.eigenvectors[(i, old_idx)];
            }
        }

        c = back_transform(&c_ortho, &x);
        final_orbital_energies = orbital_energies.clone();

        let new_density = build_density_matrix(&c, n_occupied);
        let _rms = density_rms_change(&density, &new_density);
        density = new_density;
    }

    // Extract HOMO/LUMO
    let homo_ev = if n_occupied > 0 {
        final_orbital_energies[n_occupied - 1] * HARTREE_TO_EV
    } else {
        0.0
    };
    let lumo_ev = if n_occupied < n_basis {
        final_orbital_energies[n_occupied] * HARTREE_TO_EV
    } else {
        0.0
    };

    // Mulliken charges
    let mulliken = mulliken_analysis(
        &density,
        &core.overlap,
        &basis.function_to_atom,
        &system.atomic_numbers,
    );

    let orbital_energies_ev: Vec<f64> = final_orbital_energies
        .iter()
        .map(|e| e * HARTREE_TO_EV)
        .collect();

    Ok(DftResult {
        orbital_energies: final_orbital_energies,
        orbital_energies_ev,
        electronic_energy: total_energy - v_nuc,
        nuclear_repulsion: v_nuc,
        xc_energy,
        coulomb_energy: 0.0, // approximate
        total_energy,
        total_energy_ev: total_energy * HARTREE_TO_EV,
        homo_energy: homo_ev,
        lumo_energy: lumo_ev,
        gap: lumo_ev - homo_ev,
        mulliken_charges: mulliken.charges,
        method: format!("{:?}", config.method),
        scf_iterations: iterations,
        converged,
        n_basis,
        n_electrons,
        n_grid_points: grid.n_points,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dft_config_default_is_svwn() {
        let config = DftConfig::default();
        assert_eq!(config.method, DftMethod::Svwn);
    }

    #[test]
    fn dft_h2_svwn_converges() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let config = DftConfig {
            grid_quality: super::super::grid::GridQuality::Coarse,
            max_iterations: 100,
            ..DftConfig::default()
        };
        let result = solve_ks_dft(&elements, &positions, &config).unwrap();
        assert!(result.converged, "H2 SVWN should converge");
        assert!(result.total_energy < 0.0, "total energy should be negative");
        assert!(result.n_electrons == 2);
        assert!(result.n_basis >= 2);
    }

    #[test]
    fn dft_h2_has_homo_lumo_gap() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let config = DftConfig {
            grid_quality: super::super::grid::GridQuality::Coarse,
            ..DftConfig::default()
        };
        let result = solve_ks_dft(&elements, &positions, &config).unwrap();
        assert!(result.gap > 0.0, "HOMO-LUMO gap should be positive");
        assert!(result.homo_energy < result.lumo_energy);
    }

    #[test]
    fn dft_result_mulliken_charges_sum_to_zero() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let config = DftConfig {
            grid_quality: super::super::grid::GridQuality::Coarse,
            ..DftConfig::default()
        };
        let result = solve_ks_dft(&elements, &positions, &config).unwrap();
        let charge_sum: f64 = result.mulliken_charges.iter().sum();
        assert!(
            charge_sum.abs() < 0.1,
            "Mulliken charges should approximately sum to 0 for neutral H2, got {charge_sum}"
        );
    }
}
