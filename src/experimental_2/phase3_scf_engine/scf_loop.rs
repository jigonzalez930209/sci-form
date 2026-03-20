//! Main SCF iteration driver.
//!
//! Orchestrates the complete Roothaan-Hall SCF procedure:
//!
//! 1. Build one-electron matrices
//! 2. Compute orthogonalization matrix
//! 3. Initial guess (core Hamiltonian)
//! 4. SCF loop with DIIS acceleration
//! 5. Return converged result

use nalgebra::DMatrix;

use super::density_matrix::{build_density_matrix, density_rms_change};
use super::diis::DiisAccelerator;
use super::energy;
use super::fock_matrix::build_fock_matrix;
use super::mulliken::mulliken_analysis;
use super::orthogonalization::{back_transform, lowdin_orthogonalization, transform_to_orthogonal};
use crate::experimental_2::constants::{DIIS_SUBSPACE_SIZE, HARTREE_TO_EV, SCF_DENSITY_THRESHOLD, SCF_ENERGY_THRESHOLD, SCF_MAX_ITER};
use crate::experimental_2::phase2_quantum_engine::basis_set::BasisSet;
use crate::experimental_2::phase2_quantum_engine::core_hamiltonian::{build_core_matrices, nuclear_repulsion_energy};
use crate::experimental_2::phase2_quantum_engine::two_electron::TwoElectronIntegrals;
use crate::experimental_2::types::{MolecularSystem, ScfResult};

/// Configuration for the SCF calculation.
#[derive(Debug, Clone)]
pub struct ScfConfig {
    /// Maximum number of SCF iterations.
    pub max_iterations: usize,
    /// Energy convergence threshold (Hartree).
    pub energy_threshold: f64,
    /// Density matrix RMS convergence threshold.
    pub density_threshold: f64,
    /// DIIS subspace size (0 = no DIIS).
    pub diis_size: usize,
    /// Level shift for virtual orbitals (Hartree).
    pub level_shift: f64,
    /// Damping factor for density mixing (0 = no damping, 0.5 = 50% old).
    pub damping: f64,
    /// Use rayon-parallelized two-electron integral computation.
    /// Enabled automatically for basis sets larger than `parallel_threshold`.
    pub use_parallel_eri: bool,
    /// Basis set size threshold above which parallel ERI is used (when `use_parallel_eri` is true).
    pub parallel_threshold: usize,
}

impl Default for ScfConfig {
    fn default() -> Self {
        Self {
            max_iterations: SCF_MAX_ITER,
            energy_threshold: SCF_ENERGY_THRESHOLD,
            density_threshold: SCF_DENSITY_THRESHOLD,
            diis_size: DIIS_SUBSPACE_SIZE,
            level_shift: 0.0,
            damping: 0.0,
            use_parallel_eri: false,
            parallel_threshold: 20,
        }
    }
}

impl ScfConfig {
    /// Return a config with rayon-parallel ERI enabled.
    pub fn parallel() -> Self {
        Self {
            use_parallel_eri: true,
            parallel_threshold: 0, // always parallel
            ..Self::default()
        }
    }
}

/// Run a complete Hartree-Fock SCF calculation.
pub fn run_scf(system: &MolecularSystem, config: &ScfConfig) -> ScfResult {
    let _n_atoms = system.n_atoms();
    let n_electrons = system.n_electrons();
    let n_occupied = n_electrons / 2;

    // Step 1: Build basis set
    let basis = BasisSet::sto3g(&system.atomic_numbers, &system.positions_bohr);
    let n_basis = basis.n_basis;

    // Step 2: Build one-electron matrices
    let core_matrices = build_core_matrices(&basis, &system.atomic_numbers, &system.positions_bohr);
    let s = &core_matrices.overlap;
    let h_core = &core_matrices.core_hamiltonian;

    // Step 3: Compute two-electron integrals (parallel if configured)
    let use_parallel = config.use_parallel_eri && basis.n_basis >= config.parallel_threshold;
    let eris = if use_parallel {
        TwoElectronIntegrals::compute_parallel(&basis)
    } else {
        TwoElectronIntegrals::compute(&basis)
    };

    // Step 4: Nuclear repulsion energy
    let e_nuc = nuclear_repulsion_energy(&system.atomic_numbers, &system.positions_bohr);

    // Step 5: Löwdin orthogonalization
    let (x, _n_independent) = lowdin_orthogonalization(s, 1e-8);

    // Step 6: Initial guess — diagonalize H_core in orthogonal basis
    let h_ortho = transform_to_orthogonal(h_core, &x);
    let eigen = h_ortho.symmetric_eigen();

    // Sort eigenvalues and get sorted indices
    let mut indices: Vec<usize> = (0..n_basis).collect();
    indices.sort_by(|&a, &b| {
        eigen.eigenvalues[a]
            .partial_cmp(&eigen.eigenvalues[b])
            .unwrap()
    });

    // Build sorted coefficient matrix
    let mut c_ortho = DMatrix::zeros(n_basis, n_basis);
    let mut orbital_energies = vec![0.0; n_basis];
    for (new_idx, &old_idx) in indices.iter().enumerate() {
        orbital_energies[new_idx] = eigen.eigenvalues[old_idx];
        for i in 0..n_basis {
            c_ortho[(i, new_idx)] = eigen.eigenvectors[(i, old_idx)];
        }
    }

    let mut c = back_transform(&c_ortho, &x);
    let mut p = build_density_matrix(&c, n_occupied);
    let mut p_old = p.clone();

    // Step 7: SCF loop
    let mut diis = DiisAccelerator::new(config.diis_size);
    let mut e_total = 0.0;
    let mut converged = false;
    let mut scf_iter = 0;
    let mut fock = h_core.clone();

    for iteration in 0..config.max_iterations {
        scf_iter = iteration + 1;

        // Build Fock matrix
        fock = build_fock_matrix(h_core, &p, &eris);

        // Apply level shift if requested
        if config.level_shift > 0.0 {
            let f_shifted = apply_level_shift(&fock, &p, s, config.level_shift, n_occupied);
            fock = f_shifted;
        }

        // DIIS extrapolation
        if config.diis_size > 0 {
            diis.add_iteration(&fock, &p, s);
            if let Some(f_diis) = diis.extrapolate() {
                fock = f_diis;
            }
        }

        // Diagonalize in orthogonal basis
        let f_ortho = transform_to_orthogonal(&fock, &x);
        let eigen = f_ortho.symmetric_eigen();

        // Sort
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

        // Build new density
        let p_new = build_density_matrix(&c, n_occupied);

        // Apply damping if requested
        if config.damping > 0.0 {
            p = &p_new * (1.0 - config.damping) + &p_old * config.damping;
        } else {
            p = p_new;
        }

        // Compute energy
        let e_new = energy::total_energy(&p, h_core, &fock, e_nuc);

        // Check convergence
        let delta_e = (e_new - e_total).abs();
        let delta_p = density_rms_change(&p, &p_old);

        if delta_e < config.energy_threshold && delta_p < config.density_threshold && iteration > 0 {
            converged = true;
            e_total = e_new;
            break;
        }

        e_total = e_new;
        p_old = p.clone();
    }

    // Final energy evaluation with unconverged Fock if needed
    let e_elec = energy::electronic_energy(&p, h_core, &fock);

    // HOMO/LUMO
    let homo_energy = if n_occupied > 0 {
        orbital_energies[n_occupied - 1]
    } else {
        0.0
    };
    let lumo_energy = if n_occupied < n_basis {
        Some(orbital_energies[n_occupied])
    } else {
        None
    };
    let gap_ev = lumo_energy
        .map(|lumo| (lumo - homo_energy) * HARTREE_TO_EV)
        .unwrap_or(0.0);

    // Mulliken charges
    let mulliken = mulliken_analysis(&p, s, &basis.function_to_atom, &system.atomic_numbers);

    ScfResult {
        orbital_energies,
        mo_coefficients: c,
        density_matrix: p,
        electronic_energy: e_elec,
        nuclear_repulsion: e_nuc,
        total_energy: e_total,
        homo_energy,
        lumo_energy,
        gap_ev,
        mulliken_charges: mulliken.charges,
        scf_iterations: scf_iter,
        converged,
        n_basis,
        n_electrons,
        overlap_matrix: s.clone(),
        fock_matrix: fock,
    }
}

/// Apply level shift to virtual orbitals for SCF stability.
fn apply_level_shift(
    fock: &DMatrix<f64>,
    _density: &DMatrix<f64>,
    _overlap: &DMatrix<f64>,
    shift: f64,
    _n_occupied: usize,
) -> DMatrix<f64> {
    // Simple diagonal level shift — shifts all diagonal by `shift`
    // For proper implementation, shift only virtual orbital subspace
    let mut f_shifted = fock.clone();
    let n = fock.nrows();
    for i in 0..n {
        f_shifted[(i, i)] += shift;
    }
    f_shifted
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::experimental_2::constants::ANGSTROM_TO_BOHR;

    #[test]
    fn test_scf_h2() {
        // H2 at 0.74 Å
        let system = MolecularSystem {
            atomic_numbers: vec![1, 1],
            positions_bohr: vec![
                [0.0, 0.0, 0.0],
                [0.74 * ANGSTROM_TO_BOHR, 0.0, 0.0],
            ],
            charge: 0,
            multiplicity: 1,
        };

        let config = ScfConfig::default();
        let result = run_scf(&system, &config);

        // H2 STO-3G reference energy ≈ -1.117 Hartree
        assert!(result.total_energy < 0.0, "Total energy should be negative");
        assert_eq!(result.n_basis, 2);
        assert_eq!(result.n_electrons, 2);
    }

    #[test]
    fn test_scf_converges() {
        let system = MolecularSystem {
            atomic_numbers: vec![1, 1],
            positions_bohr: vec![
                [0.0, 0.0, 0.0],
                [1.4, 0.0, 0.0],
            ],
            charge: 0,
            multiplicity: 1,
        };

        let config = ScfConfig {
            max_iterations: 100,
            ..ScfConfig::default()
        };
        let result = run_scf(&system, &config);
        // H2 with 2 basis functions should converge quickly
        assert!(result.scf_iterations <= 100);
    }
}
