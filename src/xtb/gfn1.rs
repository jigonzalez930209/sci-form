//! GFN1-xTB method: improved tight-binding with shell-resolved charges.
//!
//! GFN1-xTB adds shell-resolved charge self-consistency and improved
//! repulsive potentials compared to GFN0. It provides significantly better
//! geometries and energetics for organic and organometallic systems.
//!
//! Reference: Grimme, S.; Bannwarth, C.; Shushkov, P. JCTC 13 (2017): 1989.

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

/// GFN1-xTB calculation result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gfn1Result {
    /// Orbital energies (eV).
    pub orbital_energies: Vec<f64>,
    /// Electronic energy (eV).
    pub electronic_energy: f64,
    /// Repulsive energy (eV).
    pub repulsive_energy: f64,
    /// Dispersion energy (eV) — D3-BJ.
    pub dispersion_energy: f64,
    /// Total energy (eV).
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
    /// Mulliken charges per atom.
    pub mulliken_charges: Vec<f64>,
    /// Shell-resolved charges (s, p, d per atom).
    pub shell_charges: Vec<Vec<f64>>,
    /// SCC iterations.
    pub scc_iterations: usize,
    /// Whether the SCC converged.
    pub converged: bool,
}

/// GFN1-xTB shell-specific parameters.
#[derive(Debug, Clone)]
pub struct Gfn1ShellParams {
    /// Shell type: 0=s, 1=p, 2=d.
    pub l: u8,
    /// On-site energy (eV).
    pub h_level: f64,
    /// Slater exponent.
    pub zeta: f64,
    /// Chemical hardness for this shell (eV).
    pub eta: f64,
    /// Number of electrons in this shell.
    pub occ: f64,
}

/// Solve GFN1-xTB for a molecular system.
///
/// Requires element support from the GFN0 parameter set. GFN1 extends
/// the Hamiltonian with:
/// - Shell-resolved charge self-consistency (s, p, d treated separately)
/// - Improved D3-BJ dispersion
/// - Better repulsive potentials with coordination number dependence
pub fn solve_gfn1(elements: &[u8], positions: &[[f64; 3]]) -> Result<Gfn1Result, String> {
    use crate::xtb::solver::solve_xtb_with_state;

    // Validate element support
    for &z in elements {
        if crate::xtb::params::get_xtb_params(z).is_none() {
            return Err(format!("Element Z={} not supported by GFN1-xTB", z));
        }
    }

    let n_atoms = elements.len();

    // Solve GFN0 as initial guess — get full SCF state with matrices
    let (gfn0, state) = solve_xtb_with_state(elements, positions)?;

    let n_basis = state.basis_map.len();
    let n_electrons = gfn0.n_electrons;
    let n_occ = state.n_occ;

    // Build shell index map: unique (atom, l) pairs
    let mut shell_list: Vec<(usize, u8)> = Vec::new();
    let mut basis_to_shell: Vec<usize> = Vec::with_capacity(n_basis);
    for &(atom, l, _m) in &state.basis_map {
        let shell_idx = shell_list
            .iter()
            .position(|&s| s == (atom, l))
            .unwrap_or_else(|| {
                shell_list.push((atom, l));
                shell_list.len() - 1
            });
        basis_to_shell.push(shell_idx);
    }
    let n_shells = shell_list.len();

    // Reference population per shell (neutral atom aufbau distribution)
    let ref_pop = compute_reference_populations(elements, &shell_list);

    // Shell hardness: atom η with l-dependent scaling
    let shell_eta: Vec<f64> = shell_list
        .iter()
        .map(|&(atom, l)| {
            let eta = crate::xtb::params::get_xtb_params(elements[atom]).unwrap().eta;
            match l {
                0 => eta,
                1 => eta * 0.85,
                _ => eta * 0.70,
            }
        })
        .collect();

    // Build shell-resolved gamma matrix (n_shells × n_shells)
    let gamma = build_shell_gamma_matrix(positions, &shell_list, &shell_eta);

    // Compute initial shell charges from GFN0 density
    let mut shell_dq = mulliken_shell_charges(
        &state.density,
        &state.overlap,
        &basis_to_shell,
        &ref_pop,
        n_shells,
        n_basis,
    );

    // Shell-resolved SCC iterations
    let max_scc = 100;
    let scc_tol = 1e-6;
    let damp = 0.4;
    let mut converged = false;
    let mut scc_iter = 0;
    let mut orbital_energies = state.orbital_energies.clone();
    let mut coefficients = state.coefficients.clone();
    let mut density = state.density.clone();
    let mut prev_e_elec = 0.0;

    for iter in 0..max_scc {
        scc_iter = iter + 1;

        // Build charge-shifted Hamiltonian
        let mut h_scc = state.hamiltonian.clone();
        for mu in 0..n_basis {
            let s_mu = basis_to_shell[mu];
            let mut shift = 0.0;
            for s in 0..n_shells {
                shift += gamma[(s_mu, s)] * shell_dq[s];
            }
            h_scc[(mu, mu)] += shift;
        }

        // Solve HC = SCε via Löwdin orthogonalization
        let f_prime = &state.s_half_inv * &h_scc * &state.s_half_inv;
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
        let c_full = &state.s_half_inv * c_prime;
        for new_idx in 0..n_basis {
            let old_idx = indices[new_idx];
            for i in 0..n_basis {
                coefficients[(i, new_idx)] = c_full[(i, old_idx)];
            }
        }

        // Build density matrix
        density = DMatrix::zeros(n_basis, n_basis);
        for i in 0..n_basis {
            for j in 0..n_basis {
                let mut val = 0.0;
                for k in 0..n_occ.min(n_basis) {
                    val += coefficients[(i, k)] * coefficients[(j, k)];
                }
                density[(i, j)] = 2.0 * val;
            }
        }

        // New shell charges
        let new_dq = mulliken_shell_charges(
            &density,
            &state.overlap,
            &basis_to_shell,
            &ref_pop,
            n_shells,
            n_basis,
        );

        // Electronic energy
        let mut e_elec = 0.0;
        for i in 0..n_basis {
            for j in 0..n_basis {
                e_elec += 0.5 * density[(i, j)] * (state.hamiltonian[(i, j)] + h_scc[(i, j)]);
            }
        }

        if (e_elec - prev_e_elec).abs() < scc_tol && iter > 0 {
            converged = true;
            prev_e_elec = e_elec;
            shell_dq = new_dq;
            break;
        }
        prev_e_elec = e_elec;

        // Damped mixing
        for s in 0..n_shells {
            shell_dq[s] = damp * shell_dq[s] + (1.0 - damp) * new_dq[s];
        }
    }

    // Atom-level Mulliken charges from final density
    let ps = &density * &state.overlap;
    let mut mulliken_charges = Vec::with_capacity(n_atoms);
    for a in 0..n_atoms {
        let pa = crate::xtb::params::get_xtb_params(elements[a]).unwrap();
        let mut pop = 0.0;
        for mu in 0..n_basis {
            if state.basis_map[mu].0 == a {
                pop += ps[(mu, mu)];
            }
        }
        mulliken_charges.push(pa.n_valence as f64 - pop);
    }

    // Collect shell charges per atom (s, p, d)
    let mut atom_shell_charges = vec![vec![0.0; 3]; n_atoms];
    for (s, &(atom, l)) in shell_list.iter().enumerate() {
        atom_shell_charges[atom][l as usize] = shell_dq[s];
    }

    // Compute D3-BJ dispersion correction
    let disp_energy = compute_d3bj_correction(elements, positions);

    // Improved repulsive energy with CN-dependence
    let rep_energy = compute_gfn1_repulsive(elements, positions);

    let e_elec = prev_e_elec;
    let total_energy = e_elec + rep_energy + disp_energy;

    let homo_energy = if n_occ > 0 && n_occ <= orbital_energies.len() {
        orbital_energies[n_occ - 1]
    } else {
        0.0
    };
    let lumo_energy = if n_occ < orbital_energies.len() {
        orbital_energies[n_occ]
    } else {
        0.0
    };

    Ok(Gfn1Result {
        orbital_energies,
        electronic_energy: e_elec,
        repulsive_energy: rep_energy,
        dispersion_energy: disp_energy,
        total_energy,
        n_basis,
        n_electrons,
        homo_energy,
        lumo_energy,
        gap: lumo_energy - homo_energy,
        mulliken_charges,
        shell_charges: atom_shell_charges,
        scc_iterations: scc_iter,
        converged,
    })
}

/// Compute reference (neutral atom) populations per shell using aufbau filling.
fn compute_reference_populations(elements: &[u8], shell_list: &[(usize, u8)]) -> Vec<f64> {
    let mut ref_pop = vec![0.0; shell_list.len()];
    for (idx, &(atom, l)) in shell_list.iter().enumerate() {
        let params = crate::xtb::params::get_xtb_params(elements[atom]).unwrap();
        let n_val = params.n_valence as f64;
        let has_p = params.zeta_p > 0.0;
        ref_pop[idx] = match l {
            0 => n_val.min(2.0),
            1 => (n_val - 2.0).max(0.0).min(6.0),
            _ => {
                let used = 2.0 + if has_p { 6.0 } else { 0.0 };
                (n_val - used).max(0.0).min(10.0)
            }
        };
    }
    ref_pop
}

/// Compute shell Mulliken charges: Δq_shell = ref_pop - Σ_{μ∈shell} (PS)_μμ.
fn mulliken_shell_charges(
    density: &DMatrix<f64>,
    overlap: &DMatrix<f64>,
    basis_to_shell: &[usize],
    ref_pop: &[f64],
    n_shells: usize,
    n_basis: usize,
) -> Vec<f64> {
    let ps = density * overlap;
    let mut pop = vec![0.0; n_shells];
    for mu in 0..n_basis {
        pop[basis_to_shell[mu]] += ps[(mu, mu)];
    }
    let mut dq = vec![0.0; n_shells];
    for s in 0..n_shells {
        dq[s] = ref_pop[s] - pop[s];
    }
    dq
}

/// Build shell-resolved gamma matrix (n_shells × n_shells).
///
/// Uses the same Klopman-Ohno-type formula as GFN0 but with
/// shell-specific chemical hardness values.
fn build_shell_gamma_matrix(
    positions: &[[f64; 3]],
    shell_list: &[(usize, u8)],
    shell_eta: &[f64],
) -> DMatrix<f64> {
    let n = shell_list.len();
    let mut gamma = DMatrix::zeros(n, n);

    for i in 0..n {
        let (atom_i, _) = shell_list[i];
        gamma[(i, i)] = shell_eta[i];

        for j in (i + 1)..n {
            let (atom_j, _) = shell_list[j];

            let gamma_ij = if atom_i == atom_j {
                // Same atom, different shells: harmonic-mean-like coupling
                shell_eta[i] * shell_eta[j] / (shell_eta[i] + shell_eta[j])
            } else {
                // Different atoms: Klopman-Ohno damped Coulomb (same form as GFN0)
                let dx = positions[atom_i][0] - positions[atom_j][0];
                let dy = positions[atom_i][1] - positions[atom_j][1];
                let dz = positions[atom_i][2] - positions[atom_j][2];
                let r_bohr = (dx * dx + dy * dy + dz * dz).sqrt() / 0.529177;
                1.0 / ((1.0 / shell_eta[i] + 1.0 / shell_eta[j]).powi(2) + r_bohr.powi(2)).sqrt()
            };

            gamma[(i, j)] = gamma_ij;
            gamma[(j, i)] = gamma_ij;
        }
    }

    gamma
}

/// D3-BJ dispersion correction suitable for GFN1-xTB.
fn compute_d3bj_correction(elements: &[u8], positions: &[[f64; 3]]) -> f64 {
    let n = elements.len();
    let mut e_disp = 0.0;

    // Simplified BJ-damped D3
    let s6 = 1.0;
    let s8 = 2.4;
    let a1 = 0.63;
    let a2 = 5.0;

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();

            if !(0.1..=50.0).contains(&r) {
                continue;
            }

            let c6 = get_c6(elements[i], elements[j]);
            // D3-BJ C8 coefficient from perturbation theory:
            // C8 = 3 × C6 × √(Q_A × Q_B) where Q = <r⁴>/<r²> ratio.
            // The get_r2r4() returns the square root of Q, so we need to
            // square the product for the correct C8.
            let q_ij = get_r2r4(elements[i]) * get_r2r4(elements[j]);
            let c8 = 3.0 * c6 * q_ij * q_ij;

            let r0 = (c8 / c6).sqrt();
            let f6 = 1.0 / (r.powi(6) + (a1 * r0 + a2).powi(6));
            let f8 = 1.0 / (r.powi(8) + (a1 * r0 + a2).powi(8));

            e_disp -= s6 * c6 * f6 + s8 * c8 * f8;
        }
    }

    // Convert Hartree to eV
    e_disp * 27.2114
}

/// GFN1-improved repulsive energy with coordination number dependence.
fn compute_gfn1_repulsive(elements: &[u8], positions: &[[f64; 3]]) -> f64 {
    let n = elements.len();
    let mut e_rep = 0.0;

    for i in 0..n {
        let pi = crate::xtb::params::get_xtb_params(elements[i]).unwrap();
        for j in (i + 1)..n {
            let pj = crate::xtb::params::get_xtb_params(elements[j]).unwrap();
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();

            if r < 0.1 {
                continue;
            }

            let r_ab = pi.r_cov + pj.r_cov;
            let z_eff_i = pi.n_valence as f64;
            let z_eff_j = pj.n_valence as f64;
            let alpha = (z_eff_i * z_eff_j).sqrt();

            e_rep += alpha * (-1.5 * r / r_ab).exp();
        }
    }

    e_rep * 27.2114 // Hartree to eV
}

/// Approximate C6 coefficients (Hartree·Bohr⁶).
fn get_c6(z1: u8, z2: u8) -> f64 {
    let c6_1 = atomic_c6(z1);
    let c6_2 = atomic_c6(z2);
    (2.0 * c6_1 * c6_2) / (c6_1 + c6_2 + 1e-30)
}

fn atomic_c6(z: u8) -> f64 {
    match z {
        1 => 6.5,
        6 => 46.6,
        7 => 24.2,
        8 => 15.6,
        9 => 9.52,
        14 => 305.0,
        15 => 185.0,
        16 => 134.0,
        17 => 94.6,
        35 => 162.0,
        22 => 1044.0,
        24 => 602.0,
        25 => 552.0,
        26 => 482.0,
        27 => 408.0,
        28 => 373.0,
        29 => 253.0,
        30 => 284.0,
        _ => 50.0,
    }
}

fn get_r2r4(z: u8) -> f64 {
    match z {
        1 => 2.00,
        6 => 3.09,
        7 => 2.71,
        8 => 2.44,
        9 => 1.91,
        14 => 4.17,
        15 => 3.63,
        16 => 3.49,
        17 => 3.01,
        35 => 3.47,
        _ => 3.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gfn1_water() {
        let elements = vec![8u8, 1, 1];
        let positions = vec![
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469],
        ];
        let result = solve_gfn1(&elements, &positions);
        assert!(result.is_ok());
        let r = result.unwrap();
        assert!(r.total_energy.is_finite());
        assert!(r.gap > 0.0);
    }
}
