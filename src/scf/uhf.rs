//! Unrestricted Hartree-Fock (UHF) and Restricted Open-shell HF (ROHF).
//!
//! UHF uses separate alpha and beta Fock matrices and density matrices,
//! allowing different spatial orbitals for different spins.
//!
//! **UHF equations:**
//!   F^α = H⁰ + Σ_{λσ} [P^T_{λσ}(μν|λσ) - P^α_{λσ}(μλ|νσ)]
//!   F^β = H⁰ + Σ_{λσ} [P^T_{λσ}(μν|λσ) - P^β_{λσ}(μλ|νσ)]
//!
//! where P^T = P^α + P^β is the total density.

use nalgebra::DMatrix;

use super::basis::BasisSet;
use super::constants::{
    DIIS_SUBSPACE_SIZE, HARTREE_TO_EV, SCF_DENSITY_THRESHOLD, SCF_ENERGY_THRESHOLD, SCF_MAX_ITER,
};
use super::core_matrices::{nuclear_repulsion_energy, CoreMatrices};
use super::diis::DiisAccelerator;
use super::orthogonalization::{back_transform, lowdin_orthogonalization, transform_to_orthogonal};
use super::two_electron::TwoElectronIntegrals;
use super::types::MolecularSystem;

/// Result of an Unrestricted Hartree-Fock (UHF) calculation.
#[derive(Debug, Clone)]
pub struct UhfResult {
    /// Alpha orbital energies (Hartree, ascending).
    pub alpha_orbital_energies: Vec<f64>,
    /// Beta orbital energies (Hartree, ascending).
    pub beta_orbital_energies: Vec<f64>,
    /// Alpha MO coefficients (n_basis × n_basis).
    pub alpha_coefficients: DMatrix<f64>,
    /// Beta MO coefficients (n_basis × n_basis).
    pub beta_coefficients: DMatrix<f64>,
    /// Alpha density matrix.
    pub alpha_density: DMatrix<f64>,
    /// Beta density matrix.
    pub beta_density: DMatrix<f64>,
    /// Total density matrix (P^α + P^β).
    pub total_density: DMatrix<f64>,
    /// Electronic energy (Hartree).
    pub electronic_energy: f64,
    /// Nuclear repulsion energy (Hartree).
    pub nuclear_repulsion: f64,
    /// Total energy (Hartree).
    pub total_energy: f64,
    /// HOMO energy (Hartree) — highest among alpha/beta.
    pub homo_energy: f64,
    /// LUMO energy (Hartree) — lowest among alpha/beta.
    pub lumo_energy: Option<f64>,
    /// HOMO-LUMO gap (eV).
    pub gap_ev: f64,
    /// Mulliken charges per atom.
    pub mulliken_charges: Vec<f64>,
    /// Number of SCF iterations.
    pub scf_iterations: usize,
    /// Whether SCF converged.
    pub converged: bool,
    /// Number of basis functions.
    pub n_basis: usize,
    /// Number of alpha electrons.
    pub n_alpha: usize,
    /// Number of beta electrons.
    pub n_beta: usize,
    /// Expectation value <S²> (ideally S(S+1) for clean spin state).
    pub s2_expectation: f64,
    /// Spin contamination: <S²> - S(S+1).
    pub spin_contamination: f64,
    /// Overlap matrix.
    pub overlap_matrix: DMatrix<f64>,
    /// Alpha Fock matrix at convergence.
    pub alpha_fock: DMatrix<f64>,
    /// Beta Fock matrix at convergence.
    pub beta_fock: DMatrix<f64>,
}

/// Configuration for UHF/ROHF calculations.
#[derive(Debug, Clone)]
pub struct UhfConfig {
    pub max_iterations: usize,
    pub energy_threshold: f64,
    pub density_threshold: f64,
    pub diis_size: usize,
    pub level_shift: f64,
    pub damping: f64,
    pub use_parallel_eri: bool,
    /// If true, apply ROHF constraint after UHF convergence.
    pub rohf: bool,
}

impl Default for UhfConfig {
    fn default() -> Self {
        Self {
            max_iterations: SCF_MAX_ITER,
            energy_threshold: SCF_ENERGY_THRESHOLD,
            density_threshold: SCF_DENSITY_THRESHOLD,
            diis_size: DIIS_SUBSPACE_SIZE,
            level_shift: 0.0,
            damping: 0.0,
            use_parallel_eri: false,
            rohf: false,
        }
    }
}

/// Build single-spin density matrix: P^σ_μν = Σ_{k ∈ occ_σ} C^σ_μk · C^σ_νk
fn build_spin_density(coefficients: &DMatrix<f64>, n_occ: usize) -> DMatrix<f64> {
    let n = coefficients.nrows();
    let mut p = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..=i {
            let mut val = 0.0;
            for k in 0..n_occ {
                val += coefficients[(i, k)] * coefficients[(j, k)];
            }
            p[(i, j)] = val;
            p[(j, i)] = val;
        }
    }
    p
}

/// Build UHF Fock matrix for one spin channel.
///
/// F^σ_μν = H⁰_μν + Σ_{λσ} [P^T_{λσ}(μν|λσ) − P^σ_{λσ}(μλ|νσ)]
fn build_uhf_fock(
    h_core: &DMatrix<f64>,
    p_total: &DMatrix<f64>,
    p_spin: &DMatrix<f64>,
    eris: &TwoElectronIntegrals,
) -> DMatrix<f64> {
    let n = h_core.nrows();
    let mut fock = h_core.clone();

    for mu in 0..n {
        for nu in 0..=mu {
            let mut g = 0.0;
            for lam in 0..n {
                for sig in 0..n {
                    // Coulomb from total density
                    let j = p_total[(lam, sig)] * eris.get(mu, nu, lam, sig);
                    // Exchange from same-spin density only
                    let k = p_spin[(lam, sig)] * eris.get(mu, lam, nu, sig);
                    g += j - k;
                }
            }
            fock[(mu, nu)] += g;
            if mu != nu {
                fock[(nu, mu)] += g;
            }
        }
    }
    fock
}

/// Compute UHF electronic energy.
///
/// E_elec = 0.5 · Tr[P^T · H⁰ + P^α · F^α + P^β · F^β]
fn uhf_electronic_energy(
    p_alpha: &DMatrix<f64>,
    p_beta: &DMatrix<f64>,
    h_core: &DMatrix<f64>,
    f_alpha: &DMatrix<f64>,
    f_beta: &DMatrix<f64>,
) -> f64 {
    let n = h_core.nrows();
    let mut e = 0.0;
    for i in 0..n {
        for j in 0..n {
            let p_t = p_alpha[(i, j)] + p_beta[(i, j)];
            e += p_t * h_core[(i, j)]
                + p_alpha[(i, j)] * f_alpha[(i, j)]
                + p_beta[(i, j)] * f_beta[(i, j)];
        }
    }
    0.5 * e
}

/// Compute <S²> expectation value for a UHF wavefunction.
///
/// <S²> = S_z(S_z+1) + n_β − Σ_{i∈α,j∈β} |⟨ψ^α_i|ψ^β_j⟩|²
fn compute_s2(
    c_alpha: &DMatrix<f64>,
    c_beta: &DMatrix<f64>,
    overlap: &DMatrix<f64>,
    n_alpha: usize,
    n_beta: usize,
) -> f64 {
    let sz = 0.5 * (n_alpha as f64 - n_beta as f64);
    let mut overlap_sum = 0.0;

    // Compute S^α_β = C^α^T · S · C^β for occupied orbitals
    for i in 0..n_alpha {
        for j in 0..n_beta {
            let mut s_ij = 0.0;
            let n = overlap.nrows();
            for mu in 0..n {
                for nu in 0..n {
                    s_ij += c_alpha[(mu, i)] * overlap[(mu, nu)] * c_beta[(nu, j)];
                }
            }
            overlap_sum += s_ij * s_ij;
        }
    }

    sz * (sz + 1.0) + n_beta as f64 - overlap_sum
}

/// RMS change between two density matrices.
fn density_rms(a: &DMatrix<f64>, b: &DMatrix<f64>) -> f64 {
    let n = a.nrows();
    let mut sum = 0.0;
    for i in 0..n {
        for j in 0..n {
            let d = a[(i, j)] - b[(i, j)];
            sum += d * d;
        }
    }
    (sum / (n * n) as f64).sqrt()
}

/// Diagonalize a Fock matrix in orthogonal basis, return sorted
/// (coefficients_AO, orbital_energies).
fn diag_fock(fock: &DMatrix<f64>, x: &DMatrix<f64>, n_basis: usize) -> (DMatrix<f64>, Vec<f64>) {
    let f_ortho = transform_to_orthogonal(fock, x);
    let eigen = f_ortho.symmetric_eigen();

    let mut indices: Vec<usize> = (0..n_basis).collect();
    indices.sort_by(|&a, &b| {
        eigen.eigenvalues[a]
            .partial_cmp(&eigen.eigenvalues[b])
            .unwrap()
    });

    let mut c_ortho = DMatrix::zeros(n_basis, n_basis);
    let mut energies = vec![0.0; n_basis];
    for (new_idx, &old_idx) in indices.iter().enumerate() {
        energies[new_idx] = eigen.eigenvalues[old_idx];
        for i in 0..n_basis {
            c_ortho[(i, new_idx)] = eigen.eigenvectors[(i, old_idx)];
        }
    }

    let c = back_transform(&c_ortho, x);
    (c, energies)
}

/// Run an Unrestricted Hartree-Fock (UHF) calculation.
///
/// Supports open-shell systems with arbitrary spin multiplicity.
/// Set `multiplicity = 1` for closed-shell (equivalent to RHF),
/// `multiplicity = 2` for doublet (radical), `multiplicity = 3` for
/// triplet, etc.
pub fn run_uhf(system: &MolecularSystem, config: &UhfConfig) -> UhfResult {
    let n_electrons = system.n_electrons();
    let n_unpaired = (system.multiplicity as usize).saturating_sub(1);
    let n_alpha = (n_electrons + n_unpaired) / 2;
    let n_beta = n_electrons - n_alpha;

    // Build basis
    let basis = BasisSet::sto3g(&system.atomic_numbers, &system.positions_bohr);
    let n_basis = basis.n_basis;

    // One-electron matrices
    let core = CoreMatrices::build(&basis, &system.atomic_numbers, &system.positions_bohr);
    let s = &core.overlap;
    let h_core = &core.core_hamiltonian;
    let e_nuc = nuclear_repulsion_energy(&system.atomic_numbers, &system.positions_bohr);

    // Two-electron integrals
    let eris = {
        #[cfg(feature = "parallel")]
        {
            if config.use_parallel_eri {
                TwoElectronIntegrals::compute_parallel(&basis)
            } else {
                TwoElectronIntegrals::compute(&basis)
            }
        }
        #[cfg(not(feature = "parallel"))]
        {
            TwoElectronIntegrals::compute(&basis)
        }
    };

    // Löwdin orthogonalization
    let (x, _) = lowdin_orthogonalization(s, 1e-8);

    // Initial guess: diagonalize H_core for both spins
    let (mut c_alpha, mut e_alpha) = diag_fock(h_core, &x, n_basis);
    let (mut c_beta, mut e_beta) = diag_fock(h_core, &x, n_basis);

    let mut p_alpha = build_spin_density(&c_alpha, n_alpha);
    let mut p_beta = build_spin_density(&c_beta, n_beta);
    let mut p_total = &p_alpha + &p_beta;

    let mut p_alpha_old = p_alpha.clone();
    let mut p_beta_old = p_beta.clone();

    // SCF loop
    let mut diis_alpha = DiisAccelerator::new(config.diis_size);
    let mut diis_beta = DiisAccelerator::new(config.diis_size);
    let mut e_total = 0.0;
    let mut converged = false;
    let mut scf_iter = 0;
    let mut f_alpha = h_core.clone();
    let mut f_beta = h_core.clone();

    for iteration in 0..config.max_iterations {
        scf_iter = iteration + 1;

        // Build Fock matrices
        f_alpha = build_uhf_fock(h_core, &p_total, &p_alpha, &eris);
        f_beta = build_uhf_fock(h_core, &p_total, &p_beta, &eris);

        // Level shift
        if config.level_shift > 0.0 {
            for i in 0..n_basis {
                f_alpha[(i, i)] += config.level_shift;
                f_beta[(i, i)] += config.level_shift;
            }
        }

        // DIIS
        if config.diis_size > 0 {
            diis_alpha.add_iteration(&f_alpha, &p_alpha, s);
            diis_beta.add_iteration(&f_beta, &p_beta, s);
            if let Some(fa) = diis_alpha.extrapolate() {
                f_alpha = fa;
            }
            if let Some(fb) = diis_beta.extrapolate() {
                f_beta = fb;
            }
        }

        // Diagonalize
        let (ca_new, ea_new) = diag_fock(&f_alpha, &x, n_basis);
        let (cb_new, eb_new) = diag_fock(&f_beta, &x, n_basis);
        c_alpha = ca_new;
        c_beta = cb_new;
        e_alpha = ea_new;
        e_beta = eb_new;

        let pa_new = build_spin_density(&c_alpha, n_alpha);
        let pb_new = build_spin_density(&c_beta, n_beta);

        if config.damping > 0.0 {
            p_alpha = &pa_new * (1.0 - config.damping) + &p_alpha_old * config.damping;
            p_beta = &pb_new * (1.0 - config.damping) + &p_beta_old * config.damping;
        } else {
            p_alpha = pa_new;
            p_beta = pb_new;
        }
        p_total = &p_alpha + &p_beta;

        let e_elec = uhf_electronic_energy(&p_alpha, &p_beta, h_core, &f_alpha, &f_beta);
        let e_new = e_elec + e_nuc;

        let delta_e = (e_new - e_total).abs();
        let delta_pa = density_rms(&p_alpha, &p_alpha_old);
        let delta_pb = density_rms(&p_beta, &p_beta_old);
        let delta_p = delta_pa.max(delta_pb);

        if delta_e < config.energy_threshold && delta_p < config.density_threshold && iteration > 0
        {
            converged = true;
            e_total = e_new;
            break;
        }

        e_total = e_new;
        p_alpha_old = p_alpha.clone();
        p_beta_old = p_beta.clone();
    }

    // ROHF: apply canonical ROHF projection if requested
    if config.rohf {
        // Effective Fock = (F^α + F^β)/2 with coupling operator
        let f_eff = (&f_alpha + &f_beta) * 0.5;
        let (c_rohf, e_rohf) = diag_fock(&f_eff, &x, n_basis);
        // Use same spatial orbitals for both spins
        p_alpha = build_spin_density(&c_rohf, n_alpha);
        p_beta = build_spin_density(&c_rohf, n_beta);
        p_total = &p_alpha + &p_beta;
        c_alpha = c_rohf.clone();
        c_beta = c_rohf;
        e_alpha = e_rohf.clone();
        e_beta = e_rohf;
    }

    let e_elec = uhf_electronic_energy(&p_alpha, &p_beta, h_core, &f_alpha, &f_beta);

    // HOMO/LUMO
    let homo_a = if n_alpha > 0 {
        e_alpha[n_alpha - 1]
    } else {
        f64::NEG_INFINITY
    };
    let homo_b = if n_beta > 0 {
        e_beta[n_beta - 1]
    } else {
        f64::NEG_INFINITY
    };
    let homo = homo_a.max(homo_b);

    let lumo_a = if n_alpha < n_basis {
        Some(e_alpha[n_alpha])
    } else {
        None
    };
    let lumo_b = if n_beta < n_basis {
        Some(e_beta[n_beta])
    } else {
        None
    };
    let lumo = match (lumo_a, lumo_b) {
        (Some(a), Some(b)) => Some(a.min(b)),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        _ => None,
    };
    let gap = lumo.map(|l| (l - homo) * HARTREE_TO_EV).unwrap_or(0.0);

    // <S²>
    let s2 = compute_s2(&c_alpha, &c_beta, s, n_alpha, n_beta);
    let s_exact = 0.5 * n_unpaired as f64;
    let s2_exact = s_exact * (s_exact + 1.0);
    let contamination = s2 - s2_exact;

    // Mulliken from total density
    let mulliken = super::mulliken::mulliken_analysis(
        &p_total,
        s,
        &basis.function_to_atom,
        &system.atomic_numbers,
    );

    UhfResult {
        alpha_orbital_energies: e_alpha,
        beta_orbital_energies: e_beta,
        alpha_coefficients: c_alpha,
        beta_coefficients: c_beta,
        alpha_density: p_alpha,
        beta_density: p_beta,
        total_density: p_total,
        electronic_energy: e_elec,
        nuclear_repulsion: e_nuc,
        total_energy: e_total,
        homo_energy: homo,
        lumo_energy: lumo,
        gap_ev: gap,
        mulliken_charges: mulliken.charges,
        scf_iterations: scf_iter,
        converged,
        n_basis,
        n_alpha,
        n_beta,
        s2_expectation: s2,
        spin_contamination: contamination,
        overlap_matrix: s.clone(),
        alpha_fock: f_alpha,
        beta_fock: f_beta,
    }
}

/// Convenience: run UHF with parallel ERI.
pub fn run_uhf_parallel(system: &MolecularSystem) -> UhfResult {
    run_uhf(
        system,
        &UhfConfig {
            use_parallel_eri: true,
            ..UhfConfig::default()
        },
    )
}

/// Convenience: run ROHF (restricted open-shell).
pub fn run_rohf(system: &MolecularSystem) -> UhfResult {
    run_uhf(
        system,
        &UhfConfig {
            rohf: true,
            use_parallel_eri: true,
            ..UhfConfig::default()
        },
    )
}
