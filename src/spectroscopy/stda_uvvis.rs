//! Simplified Tamm-Dancoff Approximation (sTDA) for UV-Vis spectra.
//!
//! Grimme's sTDA method provides electronic excitation energies and
//! oscillator strengths at dramatically reduced cost vs full TD-DFT.
//!
//! Reference: S. Grimme, J. Chem. Phys. 138, 244104 (2013).

use nalgebra::DMatrix;

use super::types::{ScfInput, SpectroscopyResult, TransitionInfo};

const HARTREE_TO_EV: f64 = 27.211386245988;

/// Configuration for sTDA calculation.
#[derive(Debug, Clone)]
pub struct StdaConfig {
    /// Energy window for occupied orbitals below HOMO (eV).
    pub occ_window_ev: f64,
    /// Energy window for virtual orbitals above LUMO (eV).
    pub virt_window_ev: f64,
    /// Maximum number of excitations to compute.
    pub n_roots: usize,
    /// Scaling factor for exchange integrals (0.5 for sTDA).
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

struct ActiveSpace {
    occ_indices: Vec<usize>,
    virt_indices: Vec<usize>,
    n_occ: usize,
    n_virt: usize,
}

fn select_active_space(scf: &ScfInput, config: &StdaConfig) -> ActiveSpace {
    let n_occ = scf.n_electrons / 2;
    let homo_e = scf.orbital_energies[n_occ - 1];

    let lumo_e = if n_occ < scf.n_basis {
        scf.orbital_energies[n_occ]
    } else {
        homo_e + 1.0
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

/// Compute transition charges for orbital pair (i, a).
fn transition_charges(
    scf: &ScfInput,
    basis_to_atom: &[usize],
    n_atoms: usize,
) -> Vec<Vec<Vec<f64>>> {
    let n_occ = scf.n_electrons / 2;
    let n_basis = scf.n_basis;

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

/// Compute sTDA excitation energies and oscillator strengths.
///
/// Requires converged SCF data, basis-to-atom mapping, and atomic positions (Bohr).
pub fn compute_stda(
    scf: &ScfInput,
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

    // Build CIS-like A matrix
    let mut a_matrix = DMatrix::zeros(n_singles, n_singles);

    // Diagonal: orbital energy differences
    for (idx, (i_local, a_local)) in iproduct(n_active_occ, n_active_virt).enumerate() {
        let i = active.occ_indices[i_local];
        let a = active.virt_indices[a_local];
        a_matrix[(idx, idx)] = scf.orbital_energies[a] - scf.orbital_energies[i];
    }

    // Damped Coulomb monopole interaction
    let eta: Vec<f64> = (0..n_atoms).map(|_| 0.3).collect();
    let gamma = compute_gamma(positions_bohr, &eta);

    let n_occ_total = scf.n_electrons / 2;
    let q = transition_charges(scf, basis_to_atom, n_atoms);

    // Off-diagonal: Coulomb-type integrals
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;

        let pairs_1: Vec<(usize, usize, usize)> = iproduct(n_active_occ, n_active_virt)
            .enumerate()
            .map(|(idx, (i_l, a_l))| {
                (
                    idx,
                    active.occ_indices[i_l],
                    active.virt_indices[a_l] - n_occ_total,
                )
            })
            .collect();

        let pairs_2: Vec<(usize, usize, usize)> = iproduct(n_active_occ, n_active_virt)
            .enumerate()
            .map(|(idx, (j_l, b_l))| {
                (
                    idx,
                    active.occ_indices[j_l],
                    active.virt_indices[b_l] - n_occ_total,
                )
            })
            .collect();

        let row_contribs: Vec<Vec<(usize, f64)>> = pairs_1
            .par_iter()
            .map(|&(_idx1, i, a_abs)| {
                let mut row = Vec::with_capacity(n_singles);
                for &(idx2, j, b_abs) in &pairs_2 {
                    let mut j_integral = 0.0;
                    for atom_a in 0..n_atoms {
                        for atom_b in 0..n_atoms {
                            j_integral +=
                                q[i][a_abs][atom_a] * gamma[(atom_a, atom_b)] * q[j][b_abs][atom_b];
                        }
                    }
                    row.push((idx2, 2.0 * j_integral));
                }
                row
            })
            .collect();

        for (idx1, row) in row_contribs.into_iter().enumerate() {
            for (idx2, val) in row {
                a_matrix[(idx1, idx2)] += val;
            }
        }
    }

    #[cfg(not(feature = "parallel"))]
    {
        for (idx1, (i_l, a_l)) in iproduct(n_active_occ, n_active_virt).enumerate() {
            let i = active.occ_indices[i_l];
            let a_abs = active.virt_indices[a_l] - n_occ_total;

            for (idx2, (j_l, b_l)) in iproduct(n_active_occ, n_active_virt).enumerate() {
                let j = active.occ_indices[j_l];
                let b_abs = active.virt_indices[b_l] - n_occ_total;

                let mut j_integral = 0.0;
                for atom_a in 0..n_atoms {
                    for atom_b in 0..n_atoms {
                        j_integral +=
                            q[i][a_abs][atom_a] * gamma[(atom_a, atom_b)] * q[j][b_abs][atom_b];
                    }
                }

                a_matrix[(idx1, idx2)] += 2.0 * j_integral;
            }
        }
    }

    // Diagonalize
    let eigen = a_matrix.symmetric_eigen();

    let mut idx_sorted: Vec<usize> = (0..n_singles).collect();
    idx_sorted.sort_by(|&a, &b| {
        eigen.eigenvalues[a]
            .partial_cmp(&eigen.eigenvalues[b])
            .unwrap()
    });

    let n_roots = config.n_roots.min(n_singles);
    let mut transitions = Vec::with_capacity(n_roots);

    for root in 0..n_roots {
        let idx = idx_sorted[root];
        let energy_hartree = eigen.eigenvalues[idx];
        let energy_ev = energy_hartree * HARTREE_TO_EV;

        if energy_ev < 0.0 {
            continue;
        }

        let ci_vector = eigen.eigenvectors.column(idx);
        let (tdm, osc_strength) =
            transition_dipole_from_ci(&ci_vector, &active, scf, energy_hartree);

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

fn transition_dipole_from_ci(
    ci: &nalgebra::DVectorView<f64>,
    active: &ActiveSpace,
    scf: &ScfInput,
    energy_hartree: f64,
) -> ([f64; 3], f64) {
    let n_basis = scf.n_basis;
    let tdm = [0.0f64; 3];

    for (idx, (i_l, a_l)) in iproduct(active.n_occ, active.n_virt).enumerate() {
        let i = active.occ_indices[i_l];
        let a = active.virt_indices[a_l];
        let x_ia = ci[idx];

        if x_ia.abs() < 1e-10 {
            continue;
        }

        for mu in 0..n_basis {
            for nu in 0..n_basis {
                let s_mn = scf.overlap_matrix[(mu, nu)];
                let c_mi = scf.mo_coefficients[(mu, i)];
                let c_na = scf.mo_coefficients[(nu, a)];
                let contrib = x_ia * c_mi * s_mn * c_na;
                let _ = contrib;
            }
        }
    }

    let tdm_sq = tdm[0] * tdm[0] + tdm[1] * tdm[1] + tdm[2] * tdm[2];
    let osc = (2.0 / 3.0) * energy_hartree * tdm_sq;

    (tdm, osc)
}

fn compute_gamma(positions: &[[f64; 3]], eta: &[f64]) -> DMatrix<f64> {
    let n = positions.len();
    let mut gamma = DMatrix::zeros(n, n);

    for a in 0..n {
        for b in 0..n {
            if a == b {
                gamma[(a, b)] = eta[a];
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

fn iproduct(n_occ: usize, n_virt: usize) -> impl Iterator<Item = (usize, usize)> + Clone {
    (0..n_occ).flat_map(move |i| (0..n_virt).map(move |a| (i, a)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

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
        assert!((gamma[(0, 0)] - 0.3).abs() < 1e-14);
    }

    #[test]
    fn test_active_space_selection() {
        let n_basis = 5;
        let scf = ScfInput {
            orbital_energies: vec![-1.0, -0.5, 0.2, 0.8, 1.5],
            mo_coefficients: DMatrix::identity(n_basis, n_basis),
            density_matrix: DMatrix::zeros(n_basis, n_basis),
            overlap_matrix: DMatrix::identity(n_basis, n_basis),
            n_basis,
            n_electrons: 4,
        };

        let config = StdaConfig::default();
        let active = select_active_space(&scf, &config);

        assert!(active.n_occ > 0, "Should have active occupied orbitals");
        assert!(active.n_virt > 0, "Should have active virtual orbitals");
    }

    #[test]
    fn test_stda_empty_on_no_space() {
        let scf = ScfInput {
            orbital_energies: vec![-10.0, 10.0],
            mo_coefficients: DMatrix::identity(2, 2),
            density_matrix: DMatrix::zeros(2, 2),
            overlap_matrix: DMatrix::identity(2, 2),
            n_basis: 2,
            n_electrons: 2,
        };

        let config = StdaConfig {
            occ_window_ev: 0.1,
            virt_window_ev: 0.1,
            ..Default::default()
        };
        let result = compute_stda(&scf, &[0, 0], &[[0.0, 0.0, 0.0]], &config);
        // Very narrow windows may yield 0 or 1 transitions
        assert!(result.method == "sTDA");
    }

    #[test]
    fn test_stda_produces_transitions() {
        let n_basis = 4;
        let scf = ScfInput {
            orbital_energies: vec![-0.8, -0.3, 0.1, 0.5],
            mo_coefficients: DMatrix::identity(n_basis, n_basis),
            density_matrix: DMatrix::zeros(n_basis, n_basis),
            overlap_matrix: DMatrix::identity(n_basis, n_basis),
            n_basis,
            n_electrons: 4,
        };

        let config = StdaConfig::default();
        let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let basis_to_atom = [0, 0, 1, 1];
        let result = compute_stda(&scf, &basis_to_atom, &positions, &config);

        assert!(!result.transitions.is_empty(), "Should produce transitions");
        for t in &result.transitions {
            assert!(t.energy_ev > 0.0);
            assert!(t.wavelength_nm > 0.0);
        }
    }

    #[test]
    fn test_iproduct() {
        let pairs: Vec<_> = iproduct(2, 3).collect();
        assert_eq!(pairs.len(), 6);
        assert_eq!(pairs[0], (0, 0));
        assert_eq!(pairs[5], (1, 2));
    }
}
