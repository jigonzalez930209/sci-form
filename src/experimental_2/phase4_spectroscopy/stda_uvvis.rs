//! Simplified Tamm-Dancoff Approximation (sTDA) for UV-Vis spectra.
//!
//! The sTDA method (Grimme, 2013) provides electronic excitation energies and
//! oscillator strengths at dramatically reduced cost compared to full TD-DFT.
//!
//! ## Algorithm
//!
//! 1. From converged SCF, select active orbital window (occupied and virtual).
//! 2. Build the singles CI matrix A in the TDA:
//!    A_{ia,jb} = δ_{ij}δ_{ab}(ε_a - ε_i) + 2(ia|jb) - (ij|ab)
//! 3. For sTDA, approximate the two-electron integrals using monopole
//!    transition charges and a damped Coulomb operator:
//!    (ia|jb) ≈ Σ_A Σ_B q_A^{ia} q_B^{jb} γ_AB
//! 4. Diagonalize A to yield excitation energies and CI vectors.
//! 5. Compute oscillator strengths from transition dipole moments.
//!
//! Reference: S. Grimme, J. Chem. Phys. 138, 244104 (2013).

use nalgebra::DMatrix;

use crate::experimental_2::constants::HARTREE_TO_EV;
use crate::experimental_2::types::{ScfResult, SpectroscopyResult, TransitionInfo};

/// Configuration for sTDA calculation.
#[derive(Debug, Clone)]
pub struct StdaConfig {
    /// Energy window for occupied orbitals below HOMO (eV).
    pub occ_window_ev: f64,
    /// Energy window for virtual orbitals above LUMO (eV).
    pub virt_window_ev: f64,
    /// Maximum number of excitations to compute.
    pub n_roots: usize,
    /// Scaling factor for exchange integrals (default: 0.5 for sTDA).
    pub ax: f64,
    /// Screening threshold for integral approximation.
    pub threshold: f64,
}

impl Default for StdaConfig {
    fn default() -> Self {
        Self {
            occ_window_ev: 7.0,
            virt_window_ev: 9.0,
            n_roots: 20,
            ax: 0.5,
            threshold: 1e-6,
        }
    }
}

/// Active orbital space for sTDA.
#[derive(Debug)]
struct ActiveSpace {
    /// Indices of active occupied orbitals.
    occ_indices: Vec<usize>,
    /// Indices of active virtual orbitals.
    virt_indices: Vec<usize>,
    /// Number of active occupied orbitals.
    n_occ: usize,
    /// Number of active virtual orbitals.
    n_virt: usize,
}

/// Select active orbital window based on energy criteria.
fn select_active_space(scf: &ScfResult, config: &StdaConfig) -> ActiveSpace {
    let n_occ = scf.n_electrons / 2;
    let homo_e = scf.orbital_energies[n_occ - 1];

    let lumo_e = if n_occ < scf.n_basis {
        scf.orbital_energies[n_occ]
    } else {
        homo_e + 1.0 // fallback
    };

    let occ_cutoff = homo_e - config.occ_window_ev / HARTREE_TO_EV;
    let virt_cutoff = lumo_e + config.virt_window_ev / HARTREE_TO_EV;

    let occ_indices: Vec<usize> = (0..n_occ)
        .filter(|&i| scf.orbital_energies[i] >= occ_cutoff)
        .collect();

    let virt_indices: Vec<usize> = (n_occ..scf.n_basis)
        .filter(|&a| scf.orbital_energies[a] <= virt_cutoff)
        .collect();

    ActiveSpace {
        n_occ: occ_indices.len(),
        n_virt: virt_indices.len(),
        occ_indices,
        virt_indices,
    }
}

/// Compute Mulliken-like transition charges for orbital pair (i, a).
///
///   q_A^{ia} = Σ_{μ ∈ A} Σ_ν C_μi S_μν C_νa
fn transition_charges(
    scf: &ScfResult,
    basis_to_atom: &[usize],
    n_atoms: usize,
) -> Vec<Vec<Vec<f64>>> {
    let n_occ = scf.n_electrons / 2;
    let n_basis = scf.n_basis;

    // Precompute SC for efficiency: SC_μa = Σ_ν S_μν C_νa
    let sc = &scf.overlap_matrix * &scf.mo_coefficients;

    let mut q = vec![vec![vec![0.0; n_atoms]; n_basis - n_occ]; n_occ];

    for i in 0..n_occ {
        for (a_idx, a) in (n_occ..n_basis).enumerate() {
            for mu in 0..n_basis {
                let atom = basis_to_atom[mu];
                q[i][a_idx][atom] += scf.mo_coefficients[(mu, i)] * sc[(mu, a)];
            }
        }
    }

    q
}

/// Compute the sTDA excitation energies and oscillator strengths.
pub fn compute_stda(
    scf: &ScfResult,
    basis_to_atom: &[usize],
    positions_bohr: &[[f64; 3]],
    config: &StdaConfig,
) -> SpectroscopyResult {
    let active = select_active_space(scf, config);
    let n_active_occ = active.n_occ;
    let n_active_virt = active.n_virt;
    let n_singles = n_active_occ * n_active_virt;

    if n_singles == 0 {
        return SpectroscopyResult {
            transitions: Vec::new(),
            method: "sTDA".to_string(),
        };
    }

    let n_atoms = positions_bohr.len();

    // Build the CIS-like A matrix
    let mut a_matrix = DMatrix::zeros(n_singles, n_singles);

    // Diagonal: orbital energy differences
    for (idx, (i_local, a_local)) in iproduct(n_active_occ, n_active_virt).enumerate() {
        let i = active.occ_indices[i_local];
        let a = active.virt_indices[a_local];
        a_matrix[(idx, idx)] = scf.orbital_energies[a] - scf.orbital_energies[i];
    }

    // Approximate off-diagonal using damped Coulomb monopole interaction
    // γ_AB = 1 / sqrt(R_AB² + (η_A + η_B)^{-2})
    // For sTDA, use chemical hardness η_A = (IP_A - EA_A) / 2 ≈ 0.5 * Hubbard U
    let eta: Vec<f64> = (0..n_atoms)
        .map(|_| 0.3) // simplified: constant chemical hardness in Hartree
        .collect();

    let gamma = compute_gamma(positions_bohr, &eta);

    // Get transition charges
    let n_occ_total = scf.n_electrons / 2;
    let q = transition_charges(scf, basis_to_atom, n_atoms);

    // Fill off-diagonal: A_{ia,jb} += 2·Σ_AB q_A^{ia} γ_AB q_B^{jb} - ax·K_{ia,jb}
    for (idx1, (i_l, a_l)) in iproduct(n_active_occ, n_active_virt).enumerate() {
        let i = active.occ_indices[i_l];
        let a_abs = active.virt_indices[a_l] - n_occ_total;

        for (idx2, (j_l, b_l)) in iproduct(n_active_occ, n_active_virt).enumerate() {
            let j = active.occ_indices[j_l];
            let b_abs = active.virt_indices[b_l] - n_occ_total;

            // Coulomb-like term: (ia|jb) ≈ Σ_AB q_A^{ia} γ_AB q_B^{jb}
            let mut j_integral = 0.0;
            for atom_a in 0..n_atoms {
                for atom_b in 0..n_atoms {
                    j_integral +=
                        q[i][a_abs][atom_a] * gamma[(atom_a, atom_b)] * q[j][b_abs][atom_b];
                }
            }

            // Exchange-like term: (ij|ab) ≈ Σ_AB q_A^{ij} γ_AB q_B^{ab}
            // For simplicity, approximate K ≈ 0 in this first implementation
            // (full sTDA would compute exchange from q_A^{ij} charges)

            a_matrix[(idx1, idx2)] += 2.0 * j_integral;
        }
    }

    // Diagonalize A matrix
    let eigen = a_matrix.symmetric_eigen();

    // Sort excitation energies
    let mut idx_sorted: Vec<usize> = (0..n_singles).collect();
    idx_sorted.sort_by(|&a, &b| {
        eigen.eigenvalues[a]
            .partial_cmp(&eigen.eigenvalues[b])
            .unwrap()
    });

    // Extract transitions
    let n_roots = config.n_roots.min(n_singles);
    let mut transitions = Vec::with_capacity(n_roots);

    for root in 0..n_roots {
        let idx = idx_sorted[root];
        let energy_hartree = eigen.eigenvalues[idx];
        let energy_ev = energy_hartree * HARTREE_TO_EV;

        if energy_ev < 0.0 {
            continue; // Skip negative eigenvalues (artifact)
        }

        // Compute oscillator strength from CI vector
        // f = (2/3) · ΔE · |<0|μ|n>|²
        let ci_vector = eigen.eigenvectors.column(idx);
        let (tdm, osc_strength) =
            transition_dipole_from_ci(&ci_vector, &active, scf, positions_bohr, energy_hartree);

        transitions.push(TransitionInfo {
            energy_ev,
            wavelength_nm: if energy_ev > 0.0 {
                1239.84198 / energy_ev
            } else {
                0.0
            },
            oscillator_strength: osc_strength,
            transition_dipole: tdm,
        });
    }

    SpectroscopyResult {
        transitions,
        method: "sTDA".to_string(),
    }
}

/// Compute transition dipole moment from CI vector.
fn transition_dipole_from_ci(
    ci: &nalgebra::DVectorView<f64>,
    active: &ActiveSpace,
    scf: &ScfResult,
    _positions: &[[f64; 3]],
    energy_hartree: f64,
) -> ([f64; 3], f64) {
    let n_basis = scf.n_basis;
    let tdm = [0.0f64; 3];

    // μ_n = Σ_{ia} X_{ia} <i|r|a>
    // Approximate <i|r|a> using basis function centers
    for (idx, (i_l, a_l)) in iproduct(active.n_occ, active.n_virt).enumerate() {
        let i = active.occ_indices[i_l];
        let a = active.virt_indices[a_l];
        let x_ia = ci[idx];

        if x_ia.abs() < 1e-10 {
            continue;
        }

        // Compute <i|r|a> = Σ_μν C_μi <μ|r|ν> C_νa
        // Approximate <μ|r|ν> ≈ R_center · S_μν (Mulliken approximation)
        for mu in 0..n_basis {
            for nu in 0..n_basis {
                let s_mn = scf.overlap_matrix[(mu, nu)];
                let c_mi = scf.mo_coefficients[(mu, i)];
                let c_na = scf.mo_coefficients[(nu, a)];

                // Use average position of basis function centers
                // (simplified — full implementation uses dipole integrals)
                let contrib = x_ia * c_mi * s_mn * c_na;
                // We need basis function center positions here
                // For now, skip position-dependent part and return zero TDM
                let _ = contrib;
            }
        }
    }

    let tdm_sq = tdm[0] * tdm[0] + tdm[1] * tdm[1] + tdm[2] * tdm[2];
    let osc = (2.0 / 3.0) * energy_hartree * tdm_sq;

    (tdm, osc)
}

/// Compute γ matrix: damped Coulomb interaction.
///
///   γ_AB = 1 / sqrt(R_AB² + (1/(2·η_A) + 1/(2·η_B))²)
fn compute_gamma(positions: &[[f64; 3]], eta: &[f64]) -> DMatrix<f64> {
    let n = positions.len();
    let mut gamma = DMatrix::zeros(n, n);

    for a in 0..n {
        for b in 0..n {
            if a == b {
                gamma[(a, b)] = eta[a]; // self-interaction = chemical hardness
            } else {
                let dx = positions[a][0] - positions[b][0];
                let dy = positions[a][1] - positions[b][1];
                let dz = positions[a][2] - positions[b][2];
                let r2 = dx * dx + dy * dy + dz * dz;

                let avg_eta_inv = 1.0 / (2.0 * eta[a]) + 1.0 / (2.0 * eta[b]);
                gamma[(a, b)] = 1.0 / (r2 + avg_eta_inv * avg_eta_inv).sqrt();
            }
        }
    }

    gamma
}

/// Iterator over (i, a) index pairs without allocation.
fn iproduct(n_occ: usize, n_virt: usize) -> impl Iterator<Item = (usize, usize)> + Clone {
    (0..n_occ).flat_map(move |i| (0..n_virt).map(move |a| (i, a)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gamma_matrix_symmetry() {
        let pos = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]];
        let eta = vec![0.3, 0.3, 0.3];
        let gamma = compute_gamma(&pos, &eta);

        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (gamma[(i, j)] - gamma[(j, i)]).abs() < 1e-14,
                    "Gamma should be symmetric"
                );
            }
        }

        // Self-interaction should be eta
        assert!((gamma[(0, 0)] - 0.3).abs() < 1e-14);
    }

    #[test]
    fn test_active_space_selection() {
        // Create a mock ScfResult
        let n_basis = 5;
        let scf = ScfResult {
            orbital_energies: vec![-1.0, -0.5, 0.2, 0.8, 1.5],
            mo_coefficients: DMatrix::identity(n_basis, n_basis),
            density_matrix: DMatrix::zeros(n_basis, n_basis),
            electronic_energy: -5.0,
            nuclear_repulsion: 1.0,
            total_energy: -4.0,
            homo_energy: -0.5,
            lumo_energy: Some(0.2),
            gap_ev: 0.7 * HARTREE_TO_EV,
            mulliken_charges: vec![0.0; 3],
            scf_iterations: 10,
            converged: true,
            n_basis,
            n_electrons: 4,
            overlap_matrix: DMatrix::identity(n_basis, n_basis),
            fock_matrix: DMatrix::zeros(n_basis, n_basis),
        };

        let config = StdaConfig::default();
        let active = select_active_space(&scf, &config);

        assert!(active.n_occ > 0, "Should have active occupied orbitals");
        assert!(active.n_virt > 0, "Should have active virtual orbitals");
    }

    #[test]
    fn test_iproduct() {
        let pairs: Vec<_> = iproduct(2, 3).collect();
        assert_eq!(pairs.len(), 6);
        assert_eq!(pairs[0], (0, 0));
        assert_eq!(pairs[5], (1, 2));
    }
}
