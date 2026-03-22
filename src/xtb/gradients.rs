//! Analytical gradients for GFN0-xTB tight-binding calculations.
//!
//! Computes dE/dR for each atom using the converged SCC density and
//! Hellmann-Feynman + Pulay force expressions, plus the repulsive
//! pair potential gradient.

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

use super::params::get_xtb_params;
use super::solver::{solve_xtb_with_state, sto_overlap, ANGSTROM_TO_BOHR, EV_PER_HARTREE};

/// Result of xTB gradient computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XtbGradientResult {
    /// Energy gradient dE/dR per atom in eV/Å.
    pub gradients: Vec<[f64; 3]>,
    /// Total xTB energy in eV.
    pub energy: f64,
}

/// Compute analytical xTB energy gradients.
///
/// Runs a full xTB SCC, then computes the gradient from the converged
/// density using Hellmann-Feynman + Pulay + repulsive pair potential derivatives.
pub fn compute_xtb_gradient(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<XtbGradientResult, String> {
    let (result, state) = solve_xtb_with_state(elements, positions)?;

    let n_atoms = elements.len();
    let n_basis = state.basis_map.len();
    let n_occ = state.n_occ;

    // Energy-weighted density: W_ij = 2·Σ_{k∈occ} ε_k·C_ik·C_jk
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

    let mut gradients = vec![[0.0f64; 3]; n_atoms];
    let h_step = 1e-6;
    let k_wh = 1.75;
    let rep_alpha = 6.0;

    // Compute per-pair gradient contribution (returns grad_a; grad_b = -grad_a)
    let compute_pair = |a: usize, b: usize| -> [f64; 3] {
        let pa = get_xtb_params(elements[a]).unwrap();
        let pb = get_xtb_params(elements[b]).unwrap();

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

        // ── 1. Repulsive pair potential gradient ──
        let r_ref = pa.r_cov + pb.r_cov;
        let na = pa.n_valence as f64;
        let nb = pb.n_valence as f64;
        let exp_term = (-rep_alpha * (r_ang / r_ref - 1.0)).exp();
        let de_rep_dr = na * nb * EV_PER_HARTREE * exp_term / (r_ang * ANGSTROM_TO_BOHR)
            * (-1.0 / r_ang - rep_alpha / r_ref);
        for d in 0..3 {
            grad_a[d] += de_rep_dr * dir[d];
        }

        // ── 2. SCC charge-shift gradient ──
        let eta_sum_sq = (1.0 / pa.eta + 1.0 / pb.eta).powi(2);
        let gamma_denom = (eta_sum_sq + r_bohr * r_bohr).sqrt();
        let dgamma_dr_bohr = -r_bohr / (gamma_denom * gamma_denom * gamma_denom);
        let dgamma_dr_ang = dgamma_dr_bohr * ANGSTROM_TO_BOHR;
        let mut pop_a = 0.0;
        let mut pop_b = 0.0;
        for mu in 0..n_basis {
            if state.basis_map[mu].0 == a {
                pop_a += state.density[(mu, mu)];
            }
            if state.basis_map[mu].0 == b {
                pop_b += state.density[(mu, mu)];
            }
        }
        let de_scc_dr = 0.5 * (pop_a * state.charges[b] + pop_b * state.charges[a]) * dgamma_dr_ang;
        for d in 0..3 {
            grad_a[d] += de_scc_dr * dir[d];
        }

        // ── 3. Hellmann-Feynman + Pulay ──
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
                let scale = if la == 0 && lb == 0 {
                    1.0
                } else if la == lb {
                    0.5
                } else {
                    0.6
                };
                let s_plus = sto_overlap(za, zb, r_bohr + h_step);
                let s_minus = sto_overlap(za, zb, r_bohr - h_step);
                let ds_dr_bohr = (s_plus - s_minus) / (2.0 * h_step) * scale;
                let ds_dr_ang = ds_dr_bohr * ANGSTROM_TO_BOHR;
                let h_ii = state.h_diag[mu];
                let h_jj = state.h_diag[nu];
                let dh_dr = 0.5 * k_wh * (h_ii + h_jj) * ds_dr_ang;
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

    let pairs: Vec<(usize, usize)> = (0..n_atoms)
        .flat_map(|a| ((a + 1)..n_atoms).map(move |b| (a, b)))
        .collect();

    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        let pair_grads: Vec<(usize, usize, [f64; 3])> = pairs
            .par_iter()
            .map(|&(a, b)| (a, b, compute_pair(a, b)))
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

    Ok(XtbGradientResult {
        gradients,
        energy: result.total_energy,
    })
}

#[cfg(test)]
mod tests {
    use super::super::solver::solve_xtb;
    use super::*;

    #[test]
    fn test_xtb_gradient_h2() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let result = compute_xtb_gradient(&elements, &positions).unwrap();
        assert_eq!(result.gradients.len(), 2);
        assert!(result.energy.is_finite());
        for d in 0..3 {
            assert!(
                (result.gradients[0][d] + result.gradients[1][d]).abs() < 0.1,
                "Forces not equal and opposite: {:?}",
                result.gradients
            );
        }
    }

    #[test]
    fn test_xtb_gradient_water() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let result = compute_xtb_gradient(&elements, &positions).unwrap();
        assert_eq!(result.gradients.len(), 3);
        for g in &result.gradients {
            for &v in g {
                assert!(v.is_finite(), "Gradient must be finite");
            }
        }
        // Net force ~ zero (translational invariance)
        for d in 0..3 {
            let sum: f64 = result.gradients.iter().map(|g| g[d]).sum();
            assert!(
                sum.abs() < 1.0,
                "Net force should be near zero, got {sum:.4}"
            );
        }
    }

    #[test]
    fn test_xtb_gradient_vs_numerical() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let analytical = compute_xtb_gradient(&elements, &positions).unwrap();

        let h = 1e-5;
        for a in 0..2 {
            for d in 0..3 {
                let mut pos_p = positions.to_vec();
                let mut pos_m = positions.to_vec();
                pos_p[a][d] += h;
                pos_m[a][d] -= h;
                let e_p = solve_xtb(&elements, &pos_p).unwrap().total_energy;
                let e_m = solve_xtb(&elements, &pos_m).unwrap().total_energy;
                let numerical = (e_p - e_m) / (2.0 * h);
                let diff = (analytical.gradients[a][d] - numerical).abs();
                let scale = numerical.abs().max(1.0);
                assert!(
                    diff / scale < 0.5,
                    "Gradient mismatch atom {a} dir {d}: analytical={:.6} numerical={:.6}",
                    analytical.gradients[a][d],
                    numerical,
                );
            }
        }
    }
}
