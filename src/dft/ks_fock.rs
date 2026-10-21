//! Kohn-Sham DFT Fock matrix builder and SCF driver.
//!
//! Builds the KS-DFT Fock matrix: F^KS = H_core + J + V_XC
//! (no exact exchange for pure DFT).
//!
//! Reuses the existing SCF infrastructure (DIIS, density matrix, convergence).

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

use super::grid::{GridQuality, MolecularGrid};
use super::valence_basis::DftBasisKind;
use super::vxc_matrix::{build_vxc_matrix, build_vxc_uks};
use crate::scf::constants::{HARTREE_TO_EV, SCF_MAX_ITER};
use crate::scf::core_matrices::{nuclear_repulsion_energy_charges, CoreMatrices};
use crate::scf::density_matrix::{build_density_matrix, build_spin_density};
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
    pub charge: i32,
    pub multiplicity: u32,
    /// Force UKS even for singlets.
    pub unrestricted: bool,
    pub basis: DftBasisKind,
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
            charge: 0,
            multiplicity: 1,
            unrestricted: false,
            basis: DftBasisKind::Minimal,
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
    pub n_alpha: usize,
    pub n_beta: usize,
    pub unrestricted: bool,
    /// AO density, row-major n_basis² (total density).
    pub density: Vec<f64>,
    /// MO coefficients, row-major n_basis² (alpha / RKS).
    pub coefficients: Vec<f64>,
}

/// Run a Kohn-Sham DFT calculation.
pub fn solve_ks_dft(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &DftConfig,
) -> Result<DftResult, String> {
    let system = MolecularSystem::from_angstrom(
        elements,
        positions,
        config.charge,
        config.multiplicity.max(1),
    );
    let n_core: usize = system
        .atomic_numbers
        .iter()
        .map(|&z| super::ecp::n_core(z) as usize)
        .sum();
    let n_electrons = system.n_electrons().saturating_sub(n_core);
    let n_unpaired = (config.multiplicity.max(1) as usize).saturating_sub(1);
    let uks = config.unrestricted || n_unpaired > 0 || n_electrons % 2 == 1;
    if n_electrons == 0 {
        return Err("KS-DFT requires at least 1 valence electron".to_string());
    }
    if uks {
        return solve_uks(&system, config, n_electrons, n_unpaired);
    }
    let n_occupied = n_electrons / 2;
    if n_occupied == 0 {
        return Err("closed-shell KS-DFT needs an even electron count".to_string());
    }

    let basis = super::valence_basis::ks_basis_kind(
        &system.atomic_numbers,
        &system.positions_bohr,
        config.basis,
    )?;
    let n_basis = basis.functions.len();
    let z_eff = super::ecp::effective_charges(&system.atomic_numbers);

    // Build molecular integration grid
    let grid = MolecularGrid::build(
        &system.atomic_numbers,
        &system.positions_bohr,
        config.grid_quality,
    );

    // Build one-electron matrices
    let mut core = CoreMatrices::build_with_charges(&basis, &z_eff, &system.positions_bohr);
    super::ecp::add_local_ecp(&mut core.nuclear, &basis, &system.atomic_numbers);
    core.core_hamiltonian = &core.kinetic + &core.nuclear;
    let s = &core.overlap;
    let h_core = &core.core_hamiltonian;
    let v_nuc = nuclear_repulsion_energy_charges(&z_eff, &system.positions_bohr);

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
    // Coulomb (Hartree) energy of the converged density, in Hartree.
    let mut j_energy = 0.0;
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

        // Compute Coulomb (Hartree) energy
        j_energy = 0.5 * (&density * &j_matrix).trace();

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
        coulomb_energy: j_energy,
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
        n_alpha: n_occupied,
        n_beta: n_occupied,
        unrestricted: false,
        density: flatten_matrix(&density),
        coefficients: flatten_matrix(&c),
    })
}

fn flatten_matrix(m: &DMatrix<f64>) -> Vec<f64> {
    let n = m.nrows();
    let mut out = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            out.push(m[(i, j)]);
        }
    }
    out
}

fn coulomb_matrix(density: &DMatrix<f64>, eri: &TwoElectronIntegrals, n: usize) -> DMatrix<f64> {
    let mut j_matrix = DMatrix::zeros(n, n);
    for mu in 0..n {
        for nu in 0..n {
            let mut j_val = 0.0;
            for lam in 0..n {
                for sig in 0..n {
                    j_val += density[(lam, sig)] * eri.get(mu, nu, lam, sig);
                }
            }
            j_matrix[(mu, nu)] = j_val;
        }
    }
    j_matrix
}

fn solve_uks(
    system: &MolecularSystem,
    config: &DftConfig,
    n_electrons: usize,
    n_unpaired: usize,
) -> Result<DftResult, String> {
    let n_alpha = (n_electrons + n_unpaired) / 2;
    let n_beta = n_electrons - n_alpha;
    if n_alpha == 0 {
        return Err("UKS needs at least one alpha electron".to_string());
    }
    let basis = super::valence_basis::ks_basis_kind(
        &system.atomic_numbers,
        &system.positions_bohr,
        config.basis,
    )?;
    let n_basis = basis.functions.len();
    let z_eff = super::ecp::effective_charges(&system.atomic_numbers);
    let grid = MolecularGrid::build(
        &system.atomic_numbers,
        &system.positions_bohr,
        config.grid_quality,
    );
    let mut core = CoreMatrices::build_with_charges(&basis, &z_eff, &system.positions_bohr);
    super::ecp::add_local_ecp(&mut core.nuclear, &basis, &system.atomic_numbers);
    core.core_hamiltonian = &core.kinetic + &core.nuclear;
    let v_nuc = nuclear_repulsion_energy_charges(&z_eff, &system.positions_bohr);
    let (x, _) = lowdin_orthogonalization(&core.overlap, 1e-8);
    let eri = TwoElectronIntegrals::compute(&basis);
    let mut diis_a = DiisAccelerator::new(config.diis_size);
    let mut diis_b = DiisAccelerator::new(config.diis_size);

    let f_orth = transform_to_orthogonal(&core.core_hamiltonian, &x);
    let eigen = f_orth.symmetric_eigen();
    let mut indices: Vec<usize> = (0..n_basis).collect();
    indices.sort_by(|&a, &b| {
        eigen.eigenvalues[a]
            .partial_cmp(&eigen.eigenvalues[b])
            .unwrap()
    });
    let mut c_ortho = DMatrix::zeros(n_basis, n_basis);
    for (new_idx, &old_idx) in indices.iter().enumerate() {
        for i in 0..n_basis {
            c_ortho[(i, new_idx)] = eigen.eigenvectors[(i, old_idx)];
        }
    }
    let mut c_a = back_transform(&c_ortho, &x);
    let mut c_b = c_a.clone();
    let mut p_a = build_spin_density(&c_a, n_alpha);
    let mut p_b = build_spin_density(&c_b, n_beta);
    let mut total_energy = 0.0;
    let mut xc_energy = 0.0;
    // Coulomb (Hartree) energy of the converged total density, in Hartree.
    let mut j_energy = 0.0;
    let mut converged = false;
    let mut iterations = 0;
    let mut e_a = vec![0.0; n_basis];
    let mut e_b = vec![0.0; n_basis];

    for iter in 0..config.max_iterations {
        iterations = iter + 1;
        let p_tot = &p_a + &p_b;
        let j_matrix = coulomb_matrix(&p_tot, &eri, n_basis);
        let (vxa, vxb, exc) = build_vxc_uks(&basis, &p_a, &p_b, &grid, config.method);
        xc_energy = exc;
        let fa = &core.core_hamiltonian + &j_matrix + &vxa;
        let fb = &core.core_hamiltonian + &j_matrix + &vxb;
        j_energy = 0.5 * (&p_tot * &j_matrix).trace();
        let e1 = (&p_tot * &core.core_hamiltonian).trace();
        let new_total = e1 + j_energy + xc_energy + v_nuc;
        if iter > 0 && (new_total - total_energy).abs() < config.energy_threshold {
            total_energy = new_total;
            converged = true;
            break;
        }
        total_energy = new_total;
        diis_a.add_iteration(&fa, &p_a, &core.overlap);
        diis_b.add_iteration(&fb, &p_b, &core.overlap);
        let fa = diis_a.extrapolate().unwrap_or(fa);
        let fb = diis_b.extrapolate().unwrap_or(fb);
        for (fock, c, energies, n_occ) in [
            (&fa, &mut c_a, &mut e_a, n_alpha),
            (&fb, &mut c_b, &mut e_b, n_beta),
        ] {
            let _ = n_occ;
            let f_orth = transform_to_orthogonal(fock, &x);
            let eigen = f_orth.symmetric_eigen();
            let mut indices: Vec<usize> = (0..n_basis).collect();
            indices.sort_by(|&a, &b| {
                eigen.eigenvalues[a]
                    .partial_cmp(&eigen.eigenvalues[b])
                    .unwrap()
            });
            let mut c_o = DMatrix::zeros(n_basis, n_basis);
            for (new_idx, &old_idx) in indices.iter().enumerate() {
                energies[new_idx] = eigen.eigenvalues[old_idx];
                for i in 0..n_basis {
                    c_o[(i, new_idx)] = eigen.eigenvectors[(i, old_idx)];
                }
            }
            *c = back_transform(&c_o, &x);
        }
        p_a = build_spin_density(&c_a, n_alpha);
        p_b = build_spin_density(&c_b, n_beta);
    }

    let homo = if n_beta > 0 {
        e_a[n_alpha - 1].max(e_b[n_beta - 1])
    } else {
        e_a[n_alpha - 1]
    };
    let lumo_candidates = [
        if n_alpha < n_basis {
            Some(e_a[n_alpha])
        } else {
            None
        },
        if n_beta < n_basis {
            Some(e_b[n_beta])
        } else {
            None
        },
    ];
    let lumo = lumo_candidates
        .into_iter()
        .flatten()
        .fold(f64::INFINITY, f64::min);
    let lumo = if lumo.is_finite() { lumo } else { 0.0 };
    let p_tot = &p_a + &p_b;
    let mulliken = mulliken_analysis(
        &p_tot,
        &core.overlap,
        &basis.function_to_atom,
        &system.atomic_numbers,
    );
    let homo_ev = homo * HARTREE_TO_EV;
    let lumo_ev = lumo * HARTREE_TO_EV;
    Ok(DftResult {
        orbital_energies: e_a.clone(),
        orbital_energies_ev: e_a.iter().map(|e| e * HARTREE_TO_EV).collect(),
        electronic_energy: total_energy - v_nuc,
        nuclear_repulsion: v_nuc,
        xc_energy,
        coulomb_energy: j_energy,
        total_energy,
        total_energy_ev: total_energy * HARTREE_TO_EV,
        homo_energy: homo_ev,
        lumo_energy: lumo_ev,
        gap: lumo_ev - homo_ev,
        mulliken_charges: mulliken.charges,
        method: format!("{:?} UKS", config.method),
        scf_iterations: iterations,
        converged,
        n_basis,
        n_electrons,
        n_grid_points: grid.n_points,
        n_alpha,
        n_beta,
        unrestricted: true,
        density: flatten_matrix(&p_tot),
        coefficients: flatten_matrix(&c_a),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The Coulomb (Hartree) energy must be a real quantity, not the 0.0
    /// placeholder it used to be (audit P2-41).
    #[test]
    fn test_coulomb_energy_is_reported() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]];
        let config = DftConfig::default();
        let r = solve_ks_dft(&elements, &positions, &config).expect("KS-DFT");
        assert!(
            r.coulomb_energy.is_finite(),
            "coulomb energy is not finite: {}",
            r.coulomb_energy
        );
        assert!(
            r.coulomb_energy > 0.0,
            "coulomb energy should be positive for a bound system: {}",
            r.coulomb_energy
        );
    }

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
    fn dft_lih_uses_a_valence_basis() {
        let elements = [3u8, 1];
        let positions = [[0.0, 0.0, 0.0], [1.60, 0.0, 0.0]];
        let config = DftConfig {
            grid_quality: super::super::grid::GridQuality::Coarse,
            max_iterations: 40,
            ..DftConfig::default()
        };
        let result = solve_ks_dft(&elements, &positions, &config).unwrap();
        assert!(result.n_basis > 2, "Li must not collapse to a single 1s");
        assert!(result.total_energy.is_finite());
        assert_eq!(result.n_electrons, 4);
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

    #[test]
    fn dft_zn_atom_uses_d_shell_coulomb() {
        let elements = [30u8];
        let positions = [[0.0, 0.0, 0.0]];
        let config = DftConfig {
            grid_quality: super::super::grid::GridQuality::Coarse,
            max_iterations: 40,
            ..DftConfig::default()
        };
        let result = solve_ks_dft(&elements, &positions, &config).unwrap();
        assert!(
            result.n_basis >= 19,
            "Zn all-electron minimal basis includes 3d (got {})",
            result.n_basis
        );
        assert_eq!(result.n_electrons, 30);
        assert!(result.total_energy.is_finite() && result.total_energy < 0.0);
        assert!(result.nuclear_repulsion.abs() < 1e-12);
    }

    #[test]
    fn uks_hydrogen_atom_doublet() {
        let config = DftConfig {
            multiplicity: 2,
            grid_quality: super::super::grid::GridQuality::Coarse,
            max_iterations: 30,
            ..DftConfig::default()
        };
        let result = solve_ks_dft(&[1u8], &[[0.0, 0.0, 0.0]], &config).unwrap();
        assert!(result.unrestricted);
        assert_eq!(result.n_alpha, 1);
        assert_eq!(result.n_beta, 0);
        assert_eq!(result.n_electrons, 1);
        assert!(result.total_energy.is_finite() && result.total_energy < 0.0);
    }

    #[test]
    fn ecp_hbr_does_not_reject_bromine() {
        let config = DftConfig {
            grid_quality: super::super::grid::GridQuality::Coarse,
            max_iterations: 20,
            ..DftConfig::default()
        };
        let result =
            solve_ks_dft(&[35u8, 1], &[[0.0, 0.0, 0.0], [1.41, 0.0, 0.0]], &config).unwrap();
        assert!(result.n_electrons < 36, "Br core electrons must be removed");
        assert!(result.n_basis > 2);
        assert!(result.total_energy.is_finite());
    }

    #[test]
    fn split_valence_h2_has_more_aos_than_sto3g() {
        let config = DftConfig {
            basis: super::super::valence_basis::DftBasisKind::SplitValence,
            grid_quality: super::super::grid::GridQuality::Coarse,
            max_iterations: 30,
            ..DftConfig::default()
        };
        let result =
            solve_ks_dft(&[1u8, 1], &[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]], &config).unwrap();
        assert_eq!(result.n_basis, 4);
        assert!(result.total_energy.is_finite() && result.total_energy < 0.0);
    }
}
