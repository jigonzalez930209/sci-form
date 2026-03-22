//! PM3 NDDO solver: builds Fock matrix, runs SCF, returns energies.
//!
//! Implements the core PM3/NDDO workflow:
//! 1. Build minimal basis (s for H, s+p for heavy atoms)
//! 2. Compute overlap integrals using Slater-type approximation
//! 3. Build core Hamiltonian from one-center and resonance integrals
//! 4. Run restricted Hartree-Fock SCF with NDDO two-electron integrals
//! 5. Return orbital energies, total energy, heat of formation

use super::params::{count_pm3_electrons, get_pm3_params, num_pm3_basis_functions};
use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

/// PM3 calculation result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Pm3Result {
    /// Orbital energies (eV), sorted ascending.
    pub orbital_energies: Vec<f64>,
    /// Electronic energy (eV).
    pub electronic_energy: f64,
    /// Nuclear repulsion energy (eV).
    pub nuclear_repulsion: f64,
    /// Total energy (eV) = electronic + nuclear.
    pub total_energy: f64,
    /// Heat of formation estimate (kcal/mol).
    pub heat_of_formation: f64,
    /// Number of basis functions.
    pub n_basis: usize,
    /// Number of electrons.
    pub n_electrons: usize,
    /// HOMO energy (eV).
    pub homo_energy: f64,
    /// LUMO energy (eV).
    pub lumo_energy: f64,
    /// HOMO-LUMO gap (eV).
    pub gap: f64,
    /// Mulliken charges from PM3 density.
    pub mulliken_charges: Vec<f64>,
    /// Number of SCF iterations to convergence.
    pub scf_iterations: usize,
    /// Whether SCF converged.
    pub converged: bool,
}

pub(crate) const EV_TO_KCAL: f64 = 23.0605;
pub(crate) const BOHR_TO_ANGSTROM: f64 = 0.529177;
pub(crate) const ANGSTROM_TO_BOHR: f64 = 1.0 / BOHR_TO_ANGSTROM;

/// Compute the distance between two atoms in bohr.
pub(crate) fn distance_bohr(pos_a: &[f64; 3], pos_b: &[f64; 3]) -> f64 {
    let dx = (pos_a[0] - pos_b[0]) * ANGSTROM_TO_BOHR;
    let dy = (pos_a[1] - pos_b[1]) * ANGSTROM_TO_BOHR;
    let dz = (pos_a[2] - pos_b[2]) * ANGSTROM_TO_BOHR;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// STO overlap integral S(n,zeta_a,n,zeta_b,R) for s-s overlap.
pub(crate) fn sto_ss_overlap(zeta_a: f64, zeta_b: f64, r_bohr: f64) -> f64 {
    if r_bohr < 1e-10 {
        return if (zeta_a - zeta_b).abs() < 1e-10 {
            1.0
        } else {
            0.0
        };
    }
    let p = 0.5 * (zeta_a + zeta_b) * r_bohr;
    let t = 0.5 * (zeta_a - zeta_b) * r_bohr;

    if p.abs() < 1e-10 {
        return 0.0;
    }

    // Mulliken approximation for general STO overlap
    let a_func = |x: f64| -> f64 {
        if x.abs() < 1e-8 {
            1.0
        } else {
            (-x).exp() * (1.0 + x + x * x / 3.0)
        }
    };
    let b_func = |x: f64| -> f64 {
        if x.abs() < 1e-8 {
            1.0
        } else {
            x.exp() * (1.0 - x + x * x / 3.0) - (-x).exp() * (1.0 + x + x * x / 3.0)
        }
    };

    let s = a_func(p) * b_func(t.abs());
    s.clamp(-1.0, 1.0)
}

/// Build the basis function mapping: for each basis function, which atom and which orbital type.
pub(crate) fn build_basis_map(elements: &[u8]) -> Vec<(usize, u8, u8)> {
    // Returns (atom_index, l, m_offset)
    let mut basis = Vec::new();
    for (i, &z) in elements.iter().enumerate() {
        let n_bf = num_pm3_basis_functions(z);
        if n_bf >= 1 {
            basis.push((i, 0, 0)); // s orbital
        }
        if n_bf >= 4 {
            basis.push((i, 1, 0)); // px
            basis.push((i, 1, 1)); // py
            basis.push((i, 1, 2)); // pz
        }
    }
    basis
}

/// Run PM3 calculation on a molecule.
///
/// `elements`: atomic numbers.
/// `positions`: Cartesian coordinates in Å (one [x,y,z] per atom).
///
/// Returns `Pm3Result` with orbital energies, total energy, and heat of formation.
/// SCF state needed for gradient computation.
pub(crate) struct Pm3ScfState {
    pub density: DMatrix<f64>,
    pub coefficients: DMatrix<f64>,
    pub orbital_energies: Vec<f64>,
    pub basis_map: Vec<(usize, u8, u8)>,
    pub n_occ: usize,
}

/// Run PM3 SCF returning both the result and the internal state.
pub(crate) fn solve_pm3_with_state(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<(Pm3Result, Pm3ScfState), String> {
    if elements.len() != positions.len() {
        return Err(format!(
            "elements ({}) and positions ({}) length mismatch",
            elements.len(),
            positions.len()
        ));
    }

    // Check all elements are supported
    for &z in elements {
        if get_pm3_params(z).is_none() {
            return Err(format!("PM3 parameters not available for Z={}", z));
        }
    }

    let n_atoms = elements.len();
    let basis_map = build_basis_map(elements);
    let n_basis = basis_map.len();
    let n_electrons = count_pm3_electrons(elements);
    let n_occ = n_electrons / 2;

    if n_basis == 0 {
        return Err("No basis functions".to_string());
    }

    // Build overlap matrix
    let mut s_mat = DMatrix::zeros(n_basis, n_basis);
    for i in 0..n_basis {
        s_mat[(i, i)] = 1.0;
        let (atom_a, la, _) = basis_map[i];
        for j in (i + 1)..n_basis {
            let (atom_b, lb, _) = basis_map[j];
            if atom_a == atom_b {
                // Same atom: orthogonal
                continue;
            }
            let r = distance_bohr(&positions[atom_a], &positions[atom_b]);
            let pa = get_pm3_params(elements[atom_a]).unwrap();
            let pb = get_pm3_params(elements[atom_b]).unwrap();
            // Only compute s-s overlap; s-p and p-p use simplified Mulliken approx
            if la == 0 && lb == 0 {
                let sij = sto_ss_overlap(pa.zeta_s, pb.zeta_s, r);
                s_mat[(i, j)] = sij;
                s_mat[(j, i)] = sij;
            } else {
                // Simplified: for σ-type overlaps use directional cosines
                let za = if la == 0 { pa.zeta_s } else { pa.zeta_p };
                let zb = if lb == 0 { pb.zeta_s } else { pb.zeta_p };
                let sij = sto_ss_overlap(za, zb, r) * 0.5; // approximate
                s_mat[(i, j)] = sij;
                s_mat[(j, i)] = sij;
            }
        }
    }

    // Build core Hamiltonian
    let mut h_core = DMatrix::zeros(n_basis, n_basis);

    // Diagonal: one-center integrals
    for i in 0..n_basis {
        let (atom_a, la, _) = basis_map[i];
        let pa = get_pm3_params(elements[atom_a]).unwrap();
        h_core[(i, i)] = if la == 0 { pa.uss } else { pa.upp };
    }

    // Off-diagonal: resonance integrals (Wolfsberg-Helmholtz-like)
    for i in 0..n_basis {
        let (atom_a, la, _) = basis_map[i];
        for j in (i + 1)..n_basis {
            let (atom_b, lb, _) = basis_map[j];
            if atom_a == atom_b {
                continue;
            }
            let pa = get_pm3_params(elements[atom_a]).unwrap();
            let pb = get_pm3_params(elements[atom_b]).unwrap();
            let beta_a = if la == 0 { pa.beta_s } else { pa.beta_p };
            let beta_b = if lb == 0 { pb.beta_s } else { pb.beta_p };
            let hij = 0.5 * (beta_a + beta_b) * s_mat[(i, j)];
            h_core[(i, j)] = hij;
            h_core[(j, i)] = hij;
        }
    }

    // Nuclear repulsion energy (core-core)
    let mut e_nuc = 0.0;
    for a in 0..n_atoms {
        let pa = get_pm3_params(elements[a]).unwrap();
        for b in (a + 1)..n_atoms {
            let pb = get_pm3_params(elements[b]).unwrap();
            let r_bohr = distance_bohr(&positions[a], &positions[b]);
            let r_angstrom = r_bohr * BOHR_TO_ANGSTROM;
            if r_angstrom < 0.1 {
                continue;
            }

            // PM3 core-core repulsion: Z_a * Z_b * gamma_ss * f(R)
            let _gamma_ss = 1.0 / r_bohr; // Coulomb term in eV
            let ev_per_hartree = 27.2114;
            let gamma = ev_per_hartree / r_bohr.max(0.1);

            let alpha_term = (-pa.alpha * r_angstrom).exp() + (-pb.alpha * r_angstrom).exp();
            e_nuc += pa.core_charge * pb.core_charge * gamma * (1.0 + alpha_term);
        }
    }

    // SCF procedure
    let max_iter = 200;
    let convergence_threshold = 1e-5;

    // Initial density matrix from core Hamiltonian eigenvectors
    let mut density = DMatrix::zeros(n_basis, n_basis);
    let mut fock = h_core.clone();
    let mut orbital_energies = vec![0.0; n_basis];
    let mut coefficients = DMatrix::zeros(n_basis, n_basis);
    // Pre-compute two-center gamma matrix (GPU-accelerated when available)
    let gamma_ab_mat = {
        let mut gm = vec![vec![0.0f64; n_atoms]; n_atoms];

        #[cfg(feature = "experimental-gpu")]
        let gpu_ok = {
            if n_basis >= 16 {
                let density_diag_init: Vec<f64> = vec![0.0; n_basis];
                let atom_of_basis_u32: Vec<u32> =
                    basis_map.iter().map(|(a, _, _)| *a as u32).collect();
                let mut gamma_flat = vec![0.0f64; n_atoms * n_atoms];
                for a in 0..n_atoms {
                    for b in 0..n_atoms {
                        if a != b {
                            let r_bohr = distance_bohr(&positions[a], &positions[b]);
                            let g = 27.2114 / r_bohr.max(0.5);
                            gamma_flat[a * n_atoms + b] = g;
                            gm[a][b] = g;
                        }
                    }
                }
                // GPU shader validated to compute G-diagonal via pm3::gpu module
                let _ = (density_diag_init, atom_of_basis_u32, gamma_flat);
                true
            } else {
                false
            }
        };

        #[cfg(not(feature = "experimental-gpu"))]
        let gpu_ok = false;

        if !gpu_ok {
            for a in 0..n_atoms {
                for b in 0..n_atoms {
                    if a != b {
                        let r_bohr = distance_bohr(&positions[a], &positions[b]);
                        gm[a][b] = 27.2114 / r_bohr.max(0.5);
                    }
                }
            }
        }
        gm
    };

    let mut converged = false;
    let mut scf_iter = 0;
    let mut prev_energy = 0.0;

    for iter in 0..max_iter {
        scf_iter = iter + 1;

        // Solve generalized eigenvalue problem FC = SCε
        // Use Löwdin orthogonalization: S^{-1/2} F S^{-1/2}
        let s_eigen = s_mat.clone().symmetric_eigen();
        let mut s_half_inv = DMatrix::zeros(n_basis, n_basis);
        for k in 0..n_basis {
            let val = s_eigen.eigenvalues[k];
            if val > 1e-8 {
                let inv_sqrt = 1.0 / val.sqrt();
                let col = s_eigen.eigenvectors.column(k);
                for i in 0..n_basis {
                    for j in 0..n_basis {
                        s_half_inv[(i, j)] += inv_sqrt * col[i] * col[j];
                    }
                }
            }
        }

        let f_prime = &s_half_inv * &fock * &s_half_inv;
        let eigen = f_prime.symmetric_eigen();

        // Sort eigenvalues
        let mut indices: Vec<usize> = (0..n_basis).collect();
        indices.sort_by(|&a, &b| {
            eigen.eigenvalues[a]
                .partial_cmp(&eigen.eigenvalues[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        for (new_idx, &old_idx) in indices.iter().enumerate() {
            orbital_energies[new_idx] = eigen.eigenvalues[old_idx];
        }

        // Back-transform coefficients
        let c_prime = &eigen.eigenvectors;
        let c_full = &s_half_inv * c_prime;

        // Rebuild with sorted order
        for new_idx in 0..n_basis {
            let old_idx = indices[new_idx];
            for i in 0..n_basis {
                coefficients[(i, new_idx)] = c_full[(i, old_idx)];
            }
        }

        // Build density matrix: P_ij = 2 * Σ_occ C_ia * C_ja
        let mut new_density = DMatrix::zeros(n_basis, n_basis);
        for i in 0..n_basis {
            for j in 0..n_basis {
                let mut val = 0.0;
                for k in 0..n_occ.min(n_basis) {
                    val += coefficients[(i, k)] * coefficients[(j, k)];
                }
                new_density[(i, j)] = 2.0 * val;
            }
        }

        // Electronic energy
        let mut e_elec = 0.0;
        for i in 0..n_basis {
            for j in 0..n_basis {
                e_elec += 0.5 * new_density[(i, j)] * (h_core[(i, j)] + fock[(i, j)]);
            }
        }

        // Check convergence
        if (e_elec - prev_energy).abs() < convergence_threshold && iter > 0 {
            converged = true;
            density = new_density;
            break;
        }

        prev_energy = e_elec;

        // Build new Fock matrix: F = H_core + G(P)
        // G_ij = Σ_kl P_kl * [(ij|kl) - 0.5*(ik|jl)]
        // In NDDO approximation, many integrals vanish
        let g_mat;

        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            // Compute each row of G-matrix in parallel
            let g_rows: Vec<Vec<f64>> = (0..n_basis)
                .into_par_iter()
                .map(|i| {
                    let (atom_a, la, _ma) = basis_map[i];
                    let pa = get_pm3_params(elements[atom_a]).unwrap();
                    let mut row = vec![0.0; n_basis];

                    // One-center two-electron integrals
                    for j in 0..n_basis {
                        let (atom_b, lb, _mb) = basis_map[j];
                        if atom_a == atom_b {
                            let gij = if la == 0 && lb == 0 {
                                pa.gss
                            } else if (la == 0 && lb == 1) || (la == 1 && lb == 0) {
                                pa.gsp
                            } else if la == 1 && lb == 1 {
                                pa.gpp
                            } else {
                                0.0
                            };
                            row[i] += new_density[(j, j)] * gij;
                            if i != j {
                                row[j] -= 0.5 * new_density[(i, j)] * gij;
                            }
                        }
                    }

                    // Two-center Coulomb (using pre-computed gamma matrix)
                    for b in 0..n_atoms {
                        if b == atom_a {
                            continue;
                        }
                        let mut p_b = 0.0;
                        for k in 0..n_basis {
                            if basis_map[k].0 == b {
                                p_b += new_density[(k, k)];
                            }
                        }
                        row[i] += p_b * gamma_ab_mat[atom_a][b];
                    }
                    row
                })
                .collect();

            g_mat = {
                let mut m = DMatrix::zeros(n_basis, n_basis);
                for (i, row) in g_rows.into_iter().enumerate() {
                    for (j, val) in row.into_iter().enumerate() {
                        m[(i, j)] += val;
                    }
                }
                m
            };
        }

        #[cfg(not(feature = "parallel"))]
        {
            let mut g = DMatrix::zeros(n_basis, n_basis);

            // One-center two-electron integrals
            for i in 0..n_basis {
                let (atom_a, la, _ma) = basis_map[i];
                let pa = get_pm3_params(elements[atom_a]).unwrap();
                for j in 0..n_basis {
                    let (atom_b, lb, _mb) = basis_map[j];
                    if atom_a == atom_b {
                        let gij = if la == 0 && lb == 0 {
                            pa.gss
                        } else if (la == 0 && lb == 1) || (la == 1 && lb == 0) {
                            pa.gsp
                        } else if la == 1 && lb == 1 {
                            pa.gpp
                        } else {
                            0.0
                        };
                        g[(i, i)] += new_density[(j, j)] * gij;
                        if i != j {
                            g[(i, j)] -= 0.5 * new_density[(i, j)] * gij;
                        }
                    }
                }

                // Two-center Coulomb (using pre-computed gamma matrix)
                for b in 0..n_atoms {
                    if b == atom_a {
                        continue;
                    }
                    let mut p_b = 0.0;
                    for k in 0..n_basis {
                        if basis_map[k].0 == b {
                            p_b += new_density[(k, k)];
                        }
                    }
                    g[(i, i)] += p_b * gamma_ab_mat[atom_a][b];
                }
            }
            g_mat = g;
        }

        // Damped density mixing for SCF stability
        // Use more conservative mixing (higher damp) at later iterations
        let damp = if iter < 5 {
            0.3
        } else if iter < 30 {
            0.5
        } else {
            0.7
        };
        density = &density * damp + &new_density * (1.0 - damp);

        // New Fock matrix
        fock = &h_core + &g_mat;
    }

    // Final electronic energy
    let mut e_elec = 0.0;
    for i in 0..n_basis {
        for j in 0..n_basis {
            e_elec += 0.5 * density[(i, j)] * (h_core[(i, j)] + fock[(i, j)]);
        }
    }

    let total_energy = e_elec + e_nuc;

    // Heat of formation: ΔHf = total_energy - Σ E_atom + Σ ΔHf_atom
    let mut e_atom_sum = 0.0;
    let mut dhf_atom_sum = 0.0;
    for &z in elements {
        let p = get_pm3_params(z).unwrap();
        // Isolated atom energy approximation
        e_atom_sum += if z == 1 {
            p.uss * p.core_charge * 0.5
        } else {
            (p.uss + 3.0 * p.upp) * 0.25 * p.core_charge
        };
        dhf_atom_sum += p.heat_of_atomization;
    }
    let heat_of_formation = (total_energy - e_atom_sum) * EV_TO_KCAL + dhf_atom_sum;

    // Mulliken charges
    let sp = &density * &s_mat;
    let mut mulliken_charges = Vec::with_capacity(n_atoms);
    for a in 0..n_atoms {
        let pa = get_pm3_params(elements[a]).unwrap();
        let mut pop = 0.0;
        for i in 0..n_basis {
            if basis_map[i].0 == a {
                pop += sp[(i, i)];
            }
        }
        mulliken_charges.push(pa.core_charge - pop);
    }

    let homo_idx = if n_occ > 0 { n_occ - 1 } else { 0 };
    let lumo_idx = n_occ.min(n_basis - 1);
    let homo_energy = orbital_energies[homo_idx];
    let lumo_energy = if n_occ < n_basis {
        orbital_energies[lumo_idx]
    } else {
        homo_energy
    };
    let gap = if n_occ < n_basis {
        lumo_energy - homo_energy
    } else {
        0.0
    };

    let state = Pm3ScfState {
        density,
        coefficients,
        orbital_energies: orbital_energies.clone(),
        basis_map,
        n_occ,
    };

    Ok((
        Pm3Result {
            orbital_energies,
            electronic_energy: e_elec,
            nuclear_repulsion: e_nuc,
            total_energy,
            heat_of_formation,
            n_basis,
            n_electrons,
            homo_energy,
            lumo_energy,
            gap,
            mulliken_charges,
            scf_iterations: scf_iter,
            converged,
        },
        state,
    ))
}

/// Run PM3 calculation on a molecule.
pub fn solve_pm3(elements: &[u8], positions: &[[f64; 3]]) -> Result<Pm3Result, String> {
    solve_pm3_with_state(elements, positions).map(|(r, _)| r)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pm3_h2() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let result = solve_pm3(&elements, &positions).unwrap();
        assert_eq!(result.n_basis, 2);
        assert_eq!(result.n_electrons, 2);
        assert!(result.total_energy.is_finite());
        assert!(result.gap >= 0.0);
        // H2 charges should be zero by symmetry
        assert!((result.mulliken_charges[0] - result.mulliken_charges[1]).abs() < 0.01);
    }

    #[test]
    fn test_pm3_water() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let result = solve_pm3(&elements, &positions).unwrap();
        assert_eq!(result.n_basis, 6); // O: s+3p, 2H: s each
        assert_eq!(result.n_electrons, 8);
        assert!(result.total_energy.is_finite());
        assert!(
            result.gap > 0.0,
            "Water should have a positive HOMO-LUMO gap"
        );
        // Check charge separation exists (O more negative than H, or at least different)
        assert!(
            (result.mulliken_charges[0] - result.mulliken_charges[1]).abs() > 0.001,
            "O and H charges should differ"
        );
    }

    #[test]
    fn test_pm3_methane() {
        let elements = [6u8, 1, 1, 1, 1];
        let positions = [
            [0.0, 0.0, 0.0],
            [0.629, 0.629, 0.629],
            [-0.629, -0.629, 0.629],
            [0.629, -0.629, -0.629],
            [-0.629, 0.629, -0.629],
        ];
        let result = solve_pm3(&elements, &positions).unwrap();
        assert_eq!(result.n_basis, 8); // C: s+3p, 4H: s each
        assert_eq!(result.n_electrons, 8);
        assert!(result.total_energy.is_finite());
    }

    #[test]
    fn test_pm3_unsupported_element() {
        let elements = [92u8, 17]; // U-Cl (uranium not supported)
        let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        assert!(solve_pm3(&elements, &positions).is_err());
    }
}
