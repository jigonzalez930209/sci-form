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

use super::basis::BasisSet;
use super::constants::{
    DIIS_SUBSPACE_SIZE, HARTREE_TO_EV, SCF_DENSITY_THRESHOLD, SCF_ENERGY_THRESHOLD, SCF_MAX_ITER,
};
use super::core_matrices::{nuclear_repulsion_energy, CoreMatrices};
use super::density_matrix::{build_density_matrix, density_rms_change};
use super::diis::DiisAccelerator;
use super::energy;
use super::fock_matrix::build_fock_matrix;
use super::mulliken::mulliken_analysis;
use super::orthogonalization::{back_transform, lowdin_orthogonalization, transform_to_orthogonal};
use super::two_electron::TwoElectronIntegrals;
use super::types::{MolecularSystem, ScfResult};

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
    pub use_parallel_eri: bool,
    /// Basis set size threshold for parallel ERI.
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
            parallel_threshold: 0,
            ..Self::default()
        }
    }
}

/// Run a complete Hartree-Fock SCF calculation.
pub fn run_scf(system: &MolecularSystem, config: &ScfConfig) -> ScfResult {
    let n_electrons = system.n_electrons();
    let n_occupied = n_electrons / 2;

    // Build basis set
    let basis = BasisSet::sto3g(&system.atomic_numbers, &system.positions_bohr);
    let n_basis = basis.n_basis;

    // Build one-electron matrices
    let core_matrices = CoreMatrices::build(&basis, &system.atomic_numbers, &system.positions_bohr);
    let s = &core_matrices.overlap;
    let h_core = &core_matrices.core_hamiltonian;

    // Two-electron integrals
    let use_parallel = config.use_parallel_eri && basis.n_basis >= config.parallel_threshold;
    let eris = if use_parallel {
        #[cfg(feature = "parallel")]
        {
            TwoElectronIntegrals::compute_parallel(&basis)
        }
        #[cfg(not(feature = "parallel"))]
        {
            TwoElectronIntegrals::compute(&basis)
        }
    } else {
        TwoElectronIntegrals::compute(&basis)
    };

    // Nuclear repulsion
    let e_nuc = nuclear_repulsion_energy(&system.atomic_numbers, &system.positions_bohr);

    // Löwdin orthogonalization
    let (x, _n_independent) = lowdin_orthogonalization(s, 1e-8);

    // Initial guess — diagonalize H_core in orthogonal basis
    let h_ortho = transform_to_orthogonal(h_core, &x);
    let eigen = h_ortho.symmetric_eigen();

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
    let mut p = build_density_matrix(&c, n_occupied);
    let mut p_old = p.clone();

    // SCF loop
    let mut diis = DiisAccelerator::new(config.diis_size);
    let mut e_total = 0.0;
    let mut converged = false;
    let mut scf_iter = 0;
    let mut fock = h_core.clone();

    for iteration in 0..config.max_iterations {
        scf_iter = iteration + 1;

        fock = build_fock_matrix(h_core, &p, &eris);

        // Level shift
        if config.level_shift > 0.0 {
            let n = fock.nrows();
            for i in 0..n {
                fock[(i, i)] += config.level_shift;
            }
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
        let p_new = build_density_matrix(&c, n_occupied);

        if config.damping > 0.0 {
            p = &p_new * (1.0 - config.damping) + &p_old * config.damping;
        } else {
            p = p_new;
        }

        let e_new = energy::total_energy(&p, h_core, &fock, e_nuc);
        let delta_e = (e_new - e_total).abs();
        let delta_p = density_rms_change(&p, &p_old);

        if delta_e < config.energy_threshold && delta_p < config.density_threshold && iteration > 0
        {
            converged = true;
            e_total = e_new;
            break;
        }

        e_total = e_new;
        p_old = p.clone();
    }

    let e_elec = energy::electronic_energy(&p, h_core, &fock);

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

#[cfg(test)]
mod tests {
    use super::super::constants::ANGSTROM_TO_BOHR;
    use super::*;

    #[test]
    fn test_scf_h2() {
        let system = MolecularSystem {
            atomic_numbers: vec![1, 1],
            positions_bohr: vec![[0.0, 0.0, 0.0], [0.74 * ANGSTROM_TO_BOHR, 0.0, 0.0]],
            charge: 0,
            multiplicity: 1,
        };

        let config = ScfConfig::default();
        let result = run_scf(&system, &config);

        assert!(result.total_energy < 0.0, "Total energy should be negative");
        assert_eq!(result.n_basis, 2);
        assert_eq!(result.n_electrons, 2);
    }

    #[test]
    fn test_scf_converges() {
        let system = MolecularSystem {
            atomic_numbers: vec![1, 1],
            positions_bohr: vec![[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]],
            charge: 0,
            multiplicity: 1,
        };

        let result = run_scf(&system, &ScfConfig::default());
        assert!(result.converged, "SCF should converge for H2");
        assert!(result.scf_iterations < 50);
    }
}
