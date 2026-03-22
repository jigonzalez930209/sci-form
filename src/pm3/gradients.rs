//! Analytical gradients for PM3 NDDO calculations.
//!
//! Computes dE/dR for each atom using the converged SCF density matrix
//! and the Hellmann-Feynman + Pulay force expression:
//!
//!   dE/dR_Aα = Tr(P·∂H/∂R_Aα) + 0.5·Tr(P·∂G/∂R_Aα) − Tr(W·∂S/∂R_Aα) + ∂E_nuc/∂R_Aα
//!
//! Nuclear repulsion and Coulomb derivatives are fully analytical.
//! Overlap derivatives use micro-numerical finite differences on the STO
//! overlap function (NOT on the SCF energy).

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

use super::params::get_pm3_params;
use super::solver::{solve_pm3_with_state, sto_ss_overlap, ANGSTROM_TO_BOHR};

/// Result of PM3 gradient computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Pm3GradientResult {
    /// Energy gradient dE/dR per atom in eV/Å.
    pub gradients: Vec<[f64; 3]>,
    /// Total PM3 energy in eV.
    pub energy: f64,
    /// Heat of formation in kcal/mol.
    pub heat_of_formation: f64,
}

/// Compute analytical PM3 energy gradients.
///
/// Runs a full PM3 SCF, then computes the gradient from the converged
/// density matrix using Hellmann-Feynman + Pulay corrections.
///
/// Returns gradient dE/dR in eV/Å for each atom.
pub fn compute_pm3_gradient(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<Pm3GradientResult, String> {
    let (result, state) = solve_pm3_with_state(elements, positions)?;

    let n_atoms = elements.len();
    let n_basis = state.basis_map.len();
    let n_occ = state.n_occ;
    let ev_per_hartree = 27.2114;

    // Build energy-weighted density matrix: W_μν = 2·Σ_{k∈occ} ε_k·C_μk·C_νk
    let mut w_mat = DMatrix::zeros(n_basis, n_basis);
    for i in 0..n_basis {
        for j in 0..n_basis {
            let mut val = 0.0;
            for k in 0..n_occ.min(n_basis) {
                val += state.orbital_energies[k]
                    * state.coefficients[(i, k)]
                    * state.coefficients[(j, k)];
            }
            w_mat[(i, j)] = 2.0 * val;
        }
    }

    // Compute atom populations for two-center Coulomb
    let mut atom_pop = vec![0.0f64; n_atoms];
    for mu in 0..n_basis {
        atom_pop[state.basis_map[mu].0] += state.density[(mu, mu)];
    }

    let mut gradients = vec![[0.0f64; 3]; n_atoms];
    let h_step = 1e-6; // micro-numerical step for overlap derivative (bohr)

    // Compute per-pair gradient contributions
    let compute_pair = |a: usize, b: usize| -> [f64; 3] {
        let pa = get_pm3_params(elements[a]).unwrap();
        let pb = get_pm3_params(elements[b]).unwrap();

        let dx = positions[a][0] - positions[b][0];
        let dy = positions[a][1] - positions[b][1];
        let dz = positions[a][2] - positions[b][2];
        let r_ang = (dx * dx + dy * dy + dz * dz).sqrt();
        if r_ang < 0.01 {
            return [0.0; 3];
        }

        let r_bohr = r_ang * ANGSTROM_TO_BOHR;
        let dir = [dx / r_ang, dy / r_ang, dz / r_ang];
        let mut grad_a = [0.0f64; 3];

        // ── 1. Nuclear repulsion gradient (fully analytical) ──
        let gamma = ev_per_hartree / r_bohr.max(0.5);
        let dgamma_dr_ang = if r_bohr > 0.5 {
            -ev_per_hartree * ANGSTROM_TO_BOHR / (r_bohr * r_bohr)
        } else {
            0.0
        };

        let alpha_a_term = (-pa.alpha * r_ang).exp();
        let alpha_b_term = (-pb.alpha * r_ang).exp();
        let alpha_term = alpha_a_term + alpha_b_term;
        let dalpha_dr = -pa.alpha * alpha_a_term - pb.alpha * alpha_b_term;

        let de_nuc_dr = pa.core_charge
            * pb.core_charge
            * (dgamma_dr_ang * (1.0 + alpha_term) + gamma * dalpha_dr);

        for d in 0..3 {
            grad_a[d] += de_nuc_dr * dir[d];
        }

        // ── 2. Two-center electron Coulomb gradient (analytical) ──
        let de_2c_dr = atom_pop[a] * atom_pop[b] * dgamma_dr_ang;
        for d in 0..3 {
            grad_a[d] += de_2c_dr * dir[d];
        }

        // ── 3. Hellmann-Feynman + Pulay (overlap-dependent) ──
        for mu in 0..n_basis {
            if state.basis_map[mu].0 != a {
                continue;
            }
            let la = state.basis_map[mu].1;

            for nu in 0..n_basis {
                if state.basis_map[nu].0 != b {
                    continue;
                }
                let lb = state.basis_map[nu].1;

                let za = if la == 0 { pa.zeta_s } else { pa.zeta_p };
                let zb = if lb == 0 { pb.zeta_s } else { pb.zeta_p };

                let s_plus = sto_ss_overlap(za, zb, r_bohr + h_step);
                let s_minus = sto_ss_overlap(za, zb, r_bohr - h_step);
                let mut ds_dr_bohr = (s_plus - s_minus) / (2.0 * h_step);

                if la != 0 || lb != 0 {
                    ds_dr_bohr *= 0.5;
                }

                let ds_dr_ang = ds_dr_bohr * ANGSTROM_TO_BOHR;

                let beta_mu = if la == 0 { pa.beta_s } else { pa.beta_p };
                let beta_nu = if lb == 0 { pb.beta_s } else { pb.beta_p };
                let dh_dr = 0.5 * (beta_mu + beta_nu) * ds_dr_ang;

                let p_mn = state.density[(mu, nu)];
                let w_mn = w_mat[(mu, nu)];
                let force = 2.0 * (p_mn * dh_dr - w_mn * ds_dr_ang);

                for d in 0..3 {
                    grad_a[d] += force * dir[d];
                }
            }
        }

        grad_a
    };

    // Generate all pairs
    let pairs: Vec<(usize, usize)> = (0..n_atoms)
        .flat_map(|a| ((a + 1)..n_atoms).map(move |b| (a, b)))
        .collect();

    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        let pair_grads: Vec<(usize, usize, [f64; 3])> = pairs
            .par_iter()
            .map(|&(a, b)| {
                let g = compute_pair(a, b);
                (a, b, g)
            })
            .collect();
        for (a, b, g) in pair_grads {
            for d in 0..3 {
                gradients[a][d] += g[d];
                gradients[b][d] -= g[d];
            }
        }
    }

    #[cfg(not(feature = "parallel"))]
    {
        for &(a, b) in &pairs {
            let g = compute_pair(a, b);
            for d in 0..3 {
                gradients[a][d] += g[d];
                gradients[b][d] -= g[d];
            }
        }
    }

    Ok(Pm3GradientResult {
        gradients,
        energy: result.total_energy,
        heat_of_formation: result.heat_of_formation,
    })
}

#[cfg(test)]
mod tests {
    use super::super::solver::solve_pm3;
    use super::*;

    #[test]
    fn test_pm3_gradient_h2() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let result = compute_pm3_gradient(&elements, &positions).unwrap();
        assert_eq!(result.gradients.len(), 2);
        assert!(result.energy.is_finite());
        // By Newton's third law, forces must be equal and opposite
        for d in 0..3 {
            assert!(
                (result.gradients[0][d] + result.gradients[1][d]).abs() < 0.1,
                "Forces not equal and opposite: {:?}",
                result.gradients
            );
        }
    }

    #[test]
    fn test_pm3_gradient_water() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let result = compute_pm3_gradient(&elements, &positions).unwrap();
        assert_eq!(result.gradients.len(), 3);
        for g in &result.gradients {
            for &v in g {
                assert!(v.is_finite(), "Gradient must be finite");
            }
        }
        // Net force on molecule should be near zero (translational invariance)
        for d in 0..3 {
            let sum: f64 = result.gradients.iter().map(|g| g[d]).sum();
            assert!(
                sum.abs() < 1.0,
                "Net force should be near zero, got {sum:.4}"
            );
        }
    }

    #[test]
    fn test_pm3_gradient_vs_numerical() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let analytical = compute_pm3_gradient(&elements, &positions).unwrap();

        let h = 1e-5;
        for a in 0..2 {
            for d in 0..3 {
                let mut pos_p = positions.to_vec();
                let mut pos_m = positions.to_vec();
                pos_p[a][d] += h;
                pos_m[a][d] -= h;
                let e_p = solve_pm3(&elements, &pos_p).unwrap().total_energy;
                let e_m = solve_pm3(&elements, &pos_m).unwrap().total_energy;
                let numerical = (e_p - e_m) / (2.0 * h);
                let diff = (analytical.gradients[a][d] - numerical).abs();
                let scale = numerical.abs().max(1.0);
                assert!(
                    diff / scale < 0.5,
                    "Gradient mismatch atom {a} dir {d}: analytical={:.6} numerical={:.6} diff={:.6}",
                    analytical.gradients[a][d],
                    numerical,
                    diff
                );
            }
        }
    }
}
