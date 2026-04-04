//! GFN0-xTB-inspired tight-binding solver.
//!
//! Implements a charge-self-consistent tight-binding scheme with
//! repulsive pair potentials and Mulliken charge analysis.

use super::params::{count_xtb_electrons, get_xtb_params, num_xtb_basis_functions};
use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

/// Result of an xTB tight-binding calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XtbResult {
    /// Orbital energies (eV), sorted ascending.
    pub orbital_energies: Vec<f64>,
    /// Electronic energy (eV).
    pub electronic_energy: f64,
    /// Repulsive energy (eV).
    pub repulsive_energy: f64,
    /// Total energy (eV) = electronic + repulsive.
    pub total_energy: f64,
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
    /// Mulliken charges from TB density.
    pub mulliken_charges: Vec<f64>,
    /// Number of SCC iterations.
    pub scc_iterations: usize,
    /// Whether SCC converged.
    pub converged: bool,
}

pub(crate) const ANGSTROM_TO_BOHR: f64 = 1.889_725_988_6;
pub(crate) const EV_PER_HARTREE: f64 = 27.211_385_05;

/// Compute distance in bohr between two atoms.
pub(crate) fn distance_bohr(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = (a[0] - b[0]) * ANGSTROM_TO_BOHR;
    let dy = (a[1] - b[1]) * ANGSTROM_TO_BOHR;
    let dz = (a[2] - b[2]) * ANGSTROM_TO_BOHR;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Compute damped STO overlap integral (s-s approximation).
pub(crate) fn sto_overlap(zeta_a: f64, zeta_b: f64, r_bohr: f64) -> f64 {
    if r_bohr < 1e-10 {
        return if (zeta_a - zeta_b).abs() < 1e-10 {
            1.0
        } else {
            0.0
        };
    }
    let p = 0.5 * (zeta_a + zeta_b) * r_bohr;
    (-p).exp() * (1.0 + p + p * p / 3.0)
}

/// Build basis map: (atom_index, l_quantum, m_offset).
pub(crate) fn build_basis_map(elements: &[u8]) -> Vec<(usize, u8, u8)> {
    let mut basis = Vec::new();
    for (i, &z) in elements.iter().enumerate() {
        let n = num_xtb_basis_functions(z);
        if n >= 1 {
            basis.push((i, 0, 0));
        } // s
        if n >= 4 {
            basis.push((i, 1, 0)); // px
            basis.push((i, 1, 1)); // py
            basis.push((i, 1, 2)); // pz
        }
        if n >= 9 {
            for m in 0..5u8 {
                basis.push((i, 2, m));
            } // d orbitals
        }
    }
    basis
}

/// Run an xTB tight-binding calculation.
///
/// `elements`: atomic numbers.
/// `positions`: Cartesian coordinates in Å.
/// SCF state for xTB gradient computation and GFN1 shell SCC.
pub(crate) struct XtbScfState {
    pub density: DMatrix<f64>,
    pub coefficients: DMatrix<f64>,
    pub orbital_energies: Vec<f64>,
    pub basis_map: Vec<(usize, u8, u8)>,
    pub n_occ: usize,
    pub charges: Vec<f64>,
    pub h_diag: Vec<f64>,
    pub overlap: DMatrix<f64>,
    pub hamiltonian: DMatrix<f64>,
    pub s_half_inv: DMatrix<f64>,
}

/// Run xTB calculation returning both result and internal SCF state.
pub(crate) fn solve_xtb_with_state(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<(XtbResult, XtbScfState), String> {
    if elements.len() != positions.len() {
        return Err(format!(
            "elements ({}) and positions ({}) length mismatch",
            elements.len(),
            positions.len()
        ));
    }

    for &z in elements {
        if get_xtb_params(z).is_none() {
            return Err(format!("xTB parameters not available for Z={}", z));
        }
    }

    let n_atoms = elements.len();
    let basis_map = build_basis_map(elements);
    let n_basis = basis_map.len();
    let n_electrons = count_xtb_electrons(elements);
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
                continue;
            }
            let r = distance_bohr(&positions[atom_a], &positions[atom_b]);
            let pa = get_xtb_params(elements[atom_a]).unwrap();
            let pb = get_xtb_params(elements[atom_b]).unwrap();
            let za = match la {
                0 => pa.zeta_s,
                1 => pa.zeta_p,
                _ => pa.zeta_d,
            };
            let zb = match lb {
                0 => pb.zeta_s,
                1 => pb.zeta_p,
                _ => pb.zeta_d,
            };
            if za < 1e-10 || zb < 1e-10 {
                continue;
            }
            // Shell-dependent overlap scaling factors from Grimme's GFN0 parametrization.
            // s-s: full overlap; s-p: reduced due to angular mismatch; p-p: further reduced.
            // d-orbital scaling follows similar attenuation pattern.
            let scale = match (la, lb) {
                (0, 0) => 1.0,           // s-s
                (0, 1) | (1, 0) => 0.65, // s-p (angular mismatch)
                (1, 1) => 0.55,          // p-p σ approximation
                (0, 2) | (2, 0) => 0.40, // s-d
                (1, 2) | (2, 1) => 0.35, // p-d
                (2, 2) => 0.30,          // d-d
                _ => 0.5,
            };
            let sij = sto_overlap(za, zb, r) * scale;
            s_mat[(i, j)] = sij;
            s_mat[(j, i)] = sij;
        }
    }

    // Build Hamiltonian: H_ii = level energy, H_ij = Wolfsberg-Helmholtz
    let mut h_mat = DMatrix::zeros(n_basis, n_basis);
    for i in 0..n_basis {
        let (atom_a, la, _) = basis_map[i];
        let pa = get_xtb_params(elements[atom_a]).unwrap();
        h_mat[(i, i)] = match la {
            0 => pa.h_s,
            1 => pa.h_p,
            _ => pa.h_d,
        };
    }
    for i in 0..n_basis {
        for j in (i + 1)..n_basis {
            let (atom_a, _, _) = basis_map[i];
            let (atom_b, _, _) = basis_map[j];
            if atom_a == atom_b {
                continue;
            }
            let k_wh = 1.75;
            let hij = 0.5 * k_wh * s_mat[(i, j)] * (h_mat[(i, i)] + h_mat[(j, j)]);
            h_mat[(i, j)] = hij;
            h_mat[(j, i)] = hij;
        }
    }

    // Repulsive energy: pair potential with coordination-number damping.
    // First pass: compute coordination numbers for each atom.
    let coord_numbers: Vec<f64> = (0..n_atoms)
        .map(|a| {
            let pa = get_xtb_params(elements[a]).unwrap();
            let mut cn = 0.0;
            for b in 0..n_atoms {
                if b == a {
                    continue;
                }
                let pb = get_xtb_params(elements[b]).unwrap();
                let dx = positions[a][0] - positions[b][0];
                let dy = positions[a][1] - positions[b][1];
                let dz = positions[a][2] - positions[b][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                let r_ref = pa.r_cov + pb.r_cov;
                // Fermi-type counting function
                cn += 1.0 / (1.0 + (-16.0 * (r_ref / r - 1.0)).exp());
            }
            cn
        })
        .collect();

    let mut e_rep = 0.0;
    for a in 0..n_atoms {
        let pa = get_xtb_params(elements[a]).unwrap();
        for b in (a + 1)..n_atoms {
            let pb = get_xtb_params(elements[b]).unwrap();
            let r_ang = {
                let dx = positions[a][0] - positions[b][0];
                let dy = positions[a][1] - positions[b][1];
                let dz = positions[a][2] - positions[b][2];
                (dx * dx + dy * dy + dz * dz).sqrt()
            };
            if r_ang < 0.1 {
                continue;
            }
            let r_ref = pa.r_cov + pb.r_cov;
            // Short-range repulsive with coordination-number dependent scaling.
            // Effective Z is reduced for highly-coordinated atoms.
            let alpha = 6.0;
            let cn_a = coord_numbers[a];
            let cn_b = coord_numbers[b];
            let z_eff_a = (pa.n_valence as f64) / (1.0 + 0.1 * cn_a);
            let z_eff_b = (pb.n_valence as f64) / (1.0 + 0.1 * cn_b);
            e_rep += z_eff_a * z_eff_b * EV_PER_HARTREE / (r_ang * ANGSTROM_TO_BOHR)
                * (-alpha * (r_ang / r_ref - 1.0)).exp();
        }
    }

    // SCC (self-consistent charges) loop
    let max_iter = 250;
    let convergence = 1e-6;
    let mut charges = vec![0.0f64; n_atoms];
    let mut orbital_energies = vec![0.0; n_basis];
    let mut coefficients = DMatrix::zeros(n_basis, n_basis);
    let mut converged = false;
    let mut scc_iter = 0;
    let mut prev_e_elec = 0.0;

    // Broyden mixer for SCC convergence acceleration (same scheme as GFN2)
    let mut mixer = super::broyden::BroydenMixer::new(n_atoms, 15, 0.4);

    // Löwdin S^{-1/2}
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

    // Pre-compute atom-pair gamma matrix (GPU-accelerated when available)
    let gamma_atoms = {
        let mut gm = vec![vec![0.0f64; n_atoms]; n_atoms];

        #[cfg(feature = "experimental-gpu")]
        let gpu_ok = {
            let eta_vec: Vec<f64> = (0..n_atoms)
                .map(|a| get_xtb_params(elements[a]).unwrap().eta)
                .collect();
            let pos_bohr: Vec<[f64; 3]> = positions
                .iter()
                .map(|p| {
                    [
                        p[0] * ANGSTROM_TO_BOHR,
                        p[1] * ANGSTROM_TO_BOHR,
                        p[2] * ANGSTROM_TO_BOHR,
                    ]
                })
                .collect();
            if n_atoms >= 8 {
                if let Ok(ctx) = crate::gpu::context::GpuContext::try_create() {
                    if let Ok(gpu_gamma) =
                        super::gpu::build_xtb_gamma_gpu(&ctx, &eta_vec, &pos_bohr)
                    {
                        for a in 0..n_atoms {
                            for b in 0..n_atoms {
                                gm[a][b] = gpu_gamma[(a, b)];
                            }
                        }
                        true
                    } else {
                        false
                    }
                } else {
                    false
                }
            } else {
                false
            }
        };

        #[cfg(not(feature = "experimental-gpu"))]
        let gpu_ok = false;

        if !gpu_ok {
            for a in 0..n_atoms {
                let pa = get_xtb_params(elements[a]).unwrap();
                gm[a][a] = pa.eta; // self-interaction (eV)
                for b in (a + 1)..n_atoms {
                    let pb = get_xtb_params(elements[b]).unwrap();
                    let r_bohr = distance_bohr(&positions[a], &positions[b]);
                    // Klopman-Ohno gamma in consistent units (Hartree/bohr).
                    // Convert η from eV → Hartree, compute γ in Hartree, convert back to eV.
                    let eta_a_ha = pa.eta / EV_PER_HARTREE;
                    let eta_b_ha = pb.eta / EV_PER_HARTREE;
                    let eta_avg_ha = 0.5 * (eta_a_ha + eta_b_ha);
                    let gamma_ha =
                        1.0 / (r_bohr.powi(2) + eta_avg_ha.powi(-2)).sqrt();
                    let gamma = gamma_ha * EV_PER_HARTREE;
                    gm[a][b] = gamma;
                    gm[b][a] = gamma;
                }
            }
        }
        gm
    };

    for iter in 0..max_iter {
        scc_iter = iter + 1;

        // Store current charges in mixer before SCC step
        mixer.set(&charges);

        // Build charge-shifted Hamiltonian using pre-computed gamma matrix.
        // Diagonal-only SCC shift: H_μμ -= V_A where V_A = Σ_B γ(A,B) * q_B.
        let mut h_scc = h_mat.clone();
        for i in 0..n_basis {
            let atom_a = basis_map[i].0;
            let mut shift = 0.0;
            for b in 0..n_atoms {
                shift += gamma_atoms[atom_a][b] * charges[b];
            }
            h_scc[(i, i)] -= shift;
        }

        // Solve HC = SCε via Löwdin
        let f_prime = &s_half_inv * &h_scc * &s_half_inv;
        let eigen = f_prime.symmetric_eigen();

        let mut indices: Vec<usize> = (0..n_basis).collect();
        indices.sort_by(|&a, &b| {
            eigen.eigenvalues[a]
                .partial_cmp(&eigen.eigenvalues[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        for (new_idx, &old_idx) in indices.iter().enumerate() {
            orbital_energies[new_idx] = eigen.eigenvalues[old_idx];
        }

        let c_prime = &eigen.eigenvectors;
        let c_full = &s_half_inv * c_prime;
        for new_idx in 0..n_basis {
            let old_idx = indices[new_idx];
            for i in 0..n_basis {
                coefficients[(i, new_idx)] = c_full[(i, old_idx)];
            }
        }

        // Build density matrix
        let mut density = DMatrix::zeros(n_basis, n_basis);
        for i in 0..n_basis {
            for j in 0..n_basis {
                let mut val = 0.0;
                for k in 0..n_occ.min(n_basis) {
                    val += coefficients[(i, k)] * coefficients[(j, k)];
                }
                density[(i, j)] = 2.0 * val;
            }
        }

        // Mulliken charges
        let ps = &density * &s_mat;
        let mut new_charges = Vec::with_capacity(n_atoms);
        for a in 0..n_atoms {
            let pa = get_xtb_params(elements[a]).unwrap();
            let mut pop = 0.0;
            for i in 0..n_basis {
                if basis_map[i].0 == a {
                    pop += ps[(i, i)];
                }
            }
            new_charges.push(pa.n_valence as f64 - pop);
        }

        // Electronic energy
        let mut e_elec = 0.0;
        for i in 0..n_basis {
            for j in 0..n_basis {
                e_elec += 0.5 * density[(i, j)] * (h_mat[(i, j)] + h_scc[(i, j)]);
            }
        }

        // Convergence: energy change below threshold.
        let de = (e_elec - prev_e_elec).abs();
        if de < convergence && iter > 0 {
            converged = true;
            prev_e_elec = e_elec;
            charges = new_charges;
            break;
        }
        prev_e_elec = e_elec;

        // Broyden mixing for SCC convergence (replaces simple linear damping)
        mixer.diff(&new_charges);
        if iter > 0 {
            let _ = mixer.next();
        }
        mixer.get(&mut charges);
    }

    // Final electronic energy from last density
    let e_elec = prev_e_elec;
    let total_energy = e_elec + e_rep;

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

    // Save diagonal Hamiltonian for gradient
    let h_diag: Vec<f64> = (0..n_basis).map(|i| h_mat[(i, i)]).collect();

    let state = XtbScfState {
        density: {
            // Rebuild final density from coefficients
            let mut d = DMatrix::zeros(n_basis, n_basis);
            for i in 0..n_basis {
                for j in 0..n_basis {
                    let mut val = 0.0;
                    for k in 0..n_occ.min(n_basis) {
                        val += coefficients[(i, k)] * coefficients[(j, k)];
                    }
                    d[(i, j)] = 2.0 * val;
                }
            }
            d
        },
        coefficients: coefficients.clone(),
        orbital_energies: orbital_energies.clone(),
        basis_map,
        n_occ,
        charges: charges.clone(),
        h_diag,
        overlap: s_mat,
        hamiltonian: h_mat,
        s_half_inv,
    };

    Ok((
        XtbResult {
            orbital_energies,
            electronic_energy: e_elec,
            repulsive_energy: e_rep,
            total_energy,
            n_basis,
            n_electrons,
            homo_energy,
            lumo_energy,
            gap,
            mulliken_charges: charges,
            scc_iterations: scc_iter,
            converged,
        },
        state,
    ))
}

/// Run an xTB tight-binding calculation.
pub fn solve_xtb(elements: &[u8], positions: &[[f64; 3]]) -> Result<XtbResult, String> {
    solve_xtb_with_state(elements, positions).map(|(r, _)| r)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_xtb_h2() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let result = solve_xtb(&elements, &positions).unwrap();
        assert_eq!(result.n_basis, 2);
        assert_eq!(result.n_electrons, 2);
        assert!(result.total_energy.is_finite());
        assert!(result.gap >= 0.0);
    }

    #[test]
    fn test_xtb_water() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let result = solve_xtb(&elements, &positions).unwrap();
        assert_eq!(result.n_basis, 6);
        assert_eq!(result.n_electrons, 8);
        assert!(result.total_energy.is_finite());
        assert!(result.gap > 0.0, "Water should have a positive gap");
    }

    #[test]
    fn test_xtb_ferrocene_atom() {
        // Just Fe atom — should work with s+p+d
        let elements = [26u8];
        let positions = [[0.0, 0.0, 0.0]];
        let result = solve_xtb(&elements, &positions).unwrap();
        assert_eq!(result.n_basis, 9); // s+p+d
        assert_eq!(result.n_electrons, 8);
    }

    #[test]
    fn test_xtb_unsupported() {
        let elements = [92u8]; // uranium
        let positions = [[0.0, 0.0, 0.0]];
        assert!(solve_xtb(&elements, &positions).is_err());
    }
}
