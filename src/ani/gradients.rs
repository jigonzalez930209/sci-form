//! Analytical gradient computation (forces) for ANI potentials.
//!
//! Forces = −∇E are computed by backpropagating through:
//! 1. NN output → NN input (dE/dAEV via `FeedForwardNet::backward`)
//! 2. AEV → Cartesian coordinates (chain rule through the symmetry functions)
//!
//! The Cartesian gradient of the radial AEV component is:
//! $$\frac{\partial G_R}{\partial \vec{r}_i} = \sum_j \left[
//!   e^{-\eta(R_{ij}-R_s)^2} \left( -2\eta(R_{ij}-R_s) f_c + f_c' \right)
//! \right] \frac{\vec{r}_{ij}}{R_{ij}}$$
//!
//! The angular component additionally depends on the bond angle θ and on the
//! smoothed neighbour distance $\bar r = (R_{ij} + R_{ik})/2$; both derivatives
//! are propagated explicitly so that the returned forces are the exact
//! negative gradient of the total ANI energy (verified against central
//! finite differences in the unit tests below).

use super::aev::compute_aevs;
use super::aev_params::{species_index, AevParams, N_SPECIES};
use super::cutoff::{cosine_cutoff, cosine_cutoff_deriv};
use super::neighbor::NeighborPair;
use super::nn::FeedForwardNet;
use nalgebra::DVector;
use std::collections::HashMap;

/// Small guard used to avoid division by zero for (near-)collinear triplets.
const SIN_THETA_EPS: f64 = 1e-10;

/// Compute atomic forces from ANI potential.
///
/// Returns forces as `Vec<[f64; 3]>` for each atom (in energy-gradient units).
pub fn compute_forces(
    elements: &[u8],
    positions: &[[f64; 3]],
    neighbors: &[NeighborPair],
    params: &AevParams,
    models: &HashMap<u8, FeedForwardNet>,
) -> Vec<[f64; 3]> {
    let n = elements.len();

    // Forward pass: reuse the exact same AEV implementation as the energy
    // path so the network sees an identical input.
    let aevs = compute_aevs(elements, positions, neighbors, params);

    // dE_i/dAEV_i for every atom.
    let mut aev_grads: Vec<Vec<f64>> = vec![vec![0.0; params.total_aev_length()]; n];
    for i in 0..n {
        if species_index(elements[i]).is_none() {
            continue;
        }
        let net = match models.get(&elements[i]) {
            Some(net) => net,
            None => continue,
        };
        let input = DVector::from_vec(aevs[i].clone());
        aev_grads[i] = net.backward(&input).as_slice().to_vec();
    }

    // Build per-atom neighbor lists (both directions).
    let mut atom_neighbors: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for np in neighbors {
        let d = np.dist_sq.sqrt();
        atom_neighbors[np.i].push((np.j, d));
        atom_neighbors[np.j].push((np.i, d));
    }

    // Accumulate dE/dR (the negative gradient of the total energy).
    let mut grad_e = vec![[0.0f64; 3]; n];

    for i in 0..n {
        if species_index(elements[i]).is_none() {
            continue;
        }
        accumulate_radial_gradient(
            i,
            elements,
            positions,
            &atom_neighbors[i],
            params,
            &aev_grads[i],
            &mut grad_e,
        );
        accumulate_angular_gradient(
            i,
            elements,
            positions,
            &atom_neighbors[i],
            params,
            &aev_grads[i],
            &mut grad_e,
        );
    }

    // F = -dE/dR
    for f in &mut grad_e {
        f[0] = -f[0];
        f[1] = -f[1];
        f[2] = -f[2];
    }
    grad_e
}

/// Accumulate dE/dR contributions from the radial symmetry functions.
fn accumulate_radial_gradient(
    i: usize,
    elements: &[u8],
    positions: &[[f64; 3]],
    neighbors_i: &[(usize, f64)],
    params: &AevParams,
    aev_grad: &[f64],
    grad_e: &mut [[f64; 3]],
) {
    let rad_len = params.radial_length();

    for &(j, rij) in neighbors_i {
        if rij >= params.radial_cutoff || rij < 1e-12 {
            continue;
        }
        let sj = match species_index(elements[j]) {
            Some(s) => s,
            None => continue,
        };
        let fc = cosine_cutoff(rij, params.radial_cutoff);
        let dfc = cosine_cutoff_deriv(rij, params.radial_cutoff);
        let offset = sj * rad_len;
        let rij_inv = 1.0 / rij;

        // Direction unit vector i→j
        let dir = [
            (positions[j][0] - positions[i][0]) * rij_inv,
            (positions[j][1] - positions[i][1]) * rij_inv,
            (positions[j][2] - positions[i][2]) * rij_inv,
        ];

        let mut k = 0;
        for eta in &params.radial_eta {
            for rs in &params.radial_rs {
                let dr = rij - rs;
                let gauss = (-eta * dr * dr).exp();
                // d(AEV_k)/d(rij) = gauss * dfc + fc * gauss * (-2η·dr)
                let daev_dr = gauss * (dfc + fc * (-2.0 * eta * dr));
                let chain = aev_grad[offset + k] * daev_dr;
                for d in 0..3 {
                    grad_e[i][d] -= chain * dir[d];
                    grad_e[j][d] += chain * dir[d];
                }
                k += 1;
            }
        }
    }
}

/// Accumulate dE/dR contributions from the angular symmetry functions.
fn accumulate_angular_gradient(
    i: usize,
    elements: &[u8],
    positions: &[[f64; 3]],
    neighbors_i: &[(usize, f64)],
    params: &AevParams,
    aev_grad: &[f64],
    grad_e: &mut [[f64; 3]],
) {
    let rad_total = N_SPECIES * params.radial_length();
    let ang_len = params.angular_length();

    // Same neighbour filtering as `compute_aevs`.
    let ang_neighbors: Vec<(usize, f64)> = neighbors_i
        .iter()
        .filter(|&&(_, d)| d < params.angular_cutoff)
        .copied()
        .collect();

    for a in 0..ang_neighbors.len() {
        let (j, rij) = ang_neighbors[a];
        let sj = match species_index(elements[j]) {
            Some(s) => s,
            None => continue,
        };
        let fc_ij = cosine_cutoff(rij, params.angular_cutoff);
        let dfc_ij = cosine_cutoff_deriv(rij, params.angular_cutoff);

        for b in (a + 1)..ang_neighbors.len() {
            let (k, rik) = ang_neighbors[b];
            let sk = match species_index(elements[k]) {
                Some(s) => s,
                None => continue,
            };
            let fc_ik = cosine_cutoff(rik, params.angular_cutoff);
            let dfc_ik = cosine_cutoff_deriv(rik, params.angular_cutoff);

            let (s_lo, s_hi) = if sj <= sk { (sj, sk) } else { (sk, sj) };
            let pair_idx = s_lo * N_SPECIES - s_lo * (s_lo + 1) / 2 + s_hi;
            let offset = rad_total + pair_idx * ang_len;

            let Some((theta, dtheta_dj, dtheta_dk)) = angle_with_gradients(positions, i, j, k)
            else {
                continue;
            };

            let r_avg = 0.5 * (rij + rik);

            // Unit vectors i→j and i→k.
            let dir_j = [
                (positions[j][0] - positions[i][0]) / rij,
                (positions[j][1] - positions[i][1]) / rij,
                (positions[j][2] - positions[i][2]) / rij,
            ];
            let dir_k = [
                (positions[k][0] - positions[i][0]) / rik,
                (positions[k][1] - positions[i][1]) / rik,
                (positions[k][2] - positions[i][2]) / rik,
            ];

            let mut m = 0;
            for eta in &params.angular_eta {
                for rs in &params.angular_rs {
                    for zeta in &params.angular_zeta {
                        for theta_s in &params.angular_theta_s {
                            let cos_term = 1.0 + (theta - theta_s).cos();
                            let prefactor = 2.0f64.powf(1.0 - zeta);
                            let angular = prefactor * cos_term.powf(*zeta);
                            let radial = (-eta * (r_avg - rs).powi(2)).exp();
                            let common = angular * radial;

                            // Partial derivatives (θ, Rij, Rik, r̄ held independent).
                            // d/dθ of (1 + cos(θ - θs))^ζ
                            let d_g_d_theta = if cos_term > 0.0 {
                                prefactor
                                    * zeta
                                    * cos_term.powf(zeta - 1.0)
                                    * (-(theta - theta_s).sin())
                                    * radial
                                    * fc_ij
                                    * fc_ik
                            } else {
                                0.0
                            };
                            // d/dRij (fc only)
                            let d_g_d_rij = common * dfc_ij * fc_ik;
                            // d/dRik (fc only)
                            let d_g_d_rik = common * fc_ij * dfc_ik;
                            // d/dr̄ of exp(-η (r̄ - Rs)²) — the cutoffs also
                            // multiply this partial derivative.
                            let d_g_d_ravg = common * (-2.0 * eta * (r_avg - rs)) * fc_ij * fc_ik;

                            // Chain onto Cartesian coordinates.
                            for d in 0..3 {
                                let dg_dj = d_g_d_theta * dtheta_dj[d]
                                    + d_g_d_rij * dir_j[d]
                                    + d_g_d_ravg * 0.5 * dir_j[d];
                                let dg_dk = d_g_d_theta * dtheta_dk[d]
                                    + d_g_d_rik * dir_k[d]
                                    + d_g_d_ravg * 0.5 * dir_k[d];
                                let dg_di = -(dg_dj + dg_dk);

                                let chain = aev_grad[offset + m];
                                grad_e[i][d] += chain * dg_di;
                                grad_e[j][d] += chain * dg_dj;
                                grad_e[k][d] += chain * dg_dk;
                            }
                            m += 1;
                        }
                    }
                }
            }
        }
    }
}

/// Bond angle ∠jik together with its Cartesian gradients.
///
/// Returns `None` for (near-)collinear triplets where the angle gradient is
/// undefined.
fn angle_with_gradients(
    positions: &[[f64; 3]],
    i: usize,
    j: usize,
    k: usize,
) -> Option<(f64, [f64; 3], [f64; 3])> {
    let u = [
        positions[j][0] - positions[i][0],
        positions[j][1] - positions[i][1],
        positions[j][2] - positions[i][2],
    ];
    let w = [
        positions[k][0] - positions[i][0],
        positions[k][1] - positions[i][1],
        positions[k][2] - positions[i][2],
    ];
    let a2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    let b2 = w[0] * w[0] + w[1] * w[1] + w[2] * w[2];
    if a2 < 1e-24 || b2 < 1e-24 {
        return None;
    }
    let a = a2.sqrt();
    let b = b2.sqrt();
    let c = u[0] * w[0] + u[1] * w[1] + u[2] * w[2];
    let cos_theta = (c / (a * b)).clamp(-1.0, 1.0);
    let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
    if sin_theta < SIN_THETA_EPS {
        return None;
    }
    let theta = cos_theta.acos();

    // Total Cartesian derivatives of cosθ.
    //   cosθ = (u·w) / (a·b)  with  u = rj - ri,  w = rk - ri,
    //                             a = |u|,      b = |w|
    //   ∂cosθ/∂rj = w/(ab) - c·u/(a³b)     (c = u·w)
    //   ∂cosθ/∂rk = u/(ab) - c·w/(ab³)
    // and ∂cosθ/∂ri = -(∂cosθ/∂rj + ∂cosθ/∂rk) by translation invariance.
    let dcos_dj = [
        (w[0] / (a * b)) - (c * u[0] / (a * a * a * b)),
        (w[1] / (a * b)) - (c * u[1] / (a * a * a * b)),
        (w[2] / (a * b)) - (c * u[2] / (a * a * a * b)),
    ];
    let dcos_dk = [
        (u[0] / (a * b)) - (c * w[0] / (a * b * b * b)),
        (u[1] / (a * b)) - (c * w[1] / (a * b * b * b)),
        (u[2] / (a * b)) - (c * w[2] / (a * b * b * b)),
    ];

    // dθ = -dcosθ / sinθ
    let inv_sin = -1.0 / sin_theta;
    let dtheta_dj = [
        dcos_dj[0] * inv_sin,
        dcos_dj[1] * inv_sin,
        dcos_dj[2] * inv_sin,
    ];
    let dtheta_dk = [
        dcos_dk[0] * inv_sin,
        dcos_dk[1] * inv_sin,
        dcos_dk[2] * inv_sin,
    ];

    Some((theta, dtheta_dj, dtheta_dk))
}

#[cfg(test)]
mod tests {
    use super::super::aev_params::default_ani2x_params;
    use super::super::neighbor::CellList;
    use super::super::weights::make_test_model;
    use super::*;

    fn water() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![8u8, 1, 1],
            vec![
                [0.0, 0.0, 0.117],
                [0.0, 0.757, -0.469],
                [0.0, -0.757, -0.469],
            ],
        )
    }

    fn models_for(aev_len: usize) -> HashMap<u8, FeedForwardNet> {
        let mut models = HashMap::new();
        models.insert(8u8, make_test_model(aev_len));
        models.insert(1u8, make_test_model(aev_len));
        models
    }

    #[test]
    fn test_forces_sum_near_zero() {
        // Newton's third law: total force should be ~zero
        let (elements, positions) = water();
        let params = default_ani2x_params();
        let cl = CellList::new(&positions, params.radial_cutoff);
        let neighbors = cl.find_neighbors(&positions);
        let models = models_for(params.total_aev_length());

        let forces = compute_forces(&elements, &positions, &neighbors, &params, &models);

        let total = [
            forces.iter().map(|f| f[0]).sum::<f64>(),
            forces.iter().map(|f| f[1]).sum::<f64>(),
            forces.iter().map(|f| f[2]).sum::<f64>(),
        ];
        let magnitude = (total[0] * total[0] + total[1] * total[1] + total[2] * total[2]).sqrt();
        assert!(
            magnitude < 1e-6,
            "Total force should be ~0 (Newton's 3rd law), got {magnitude}"
        );
    }

    #[test]
    fn test_forces_match_finite_difference_energy() {
        // The analytic forces must equal -dE/dR of the total ANI energy.
        use super::super::api::{compute_ani, AniConfig};

        let (elements, positions) = water();
        let params = default_ani2x_params();
        let models = models_for(params.total_aev_length());

        let cl = CellList::new(&positions, params.radial_cutoff);
        let neighbors = cl.find_neighbors(&positions);
        let forces = compute_forces(&elements, &positions, &neighbors, &params, &models);

        let energy_at = |pos: &[[f64; 3]]| -> f64 {
            let cfg = AniConfig {
                cutoff: params.radial_cutoff,
                compute_forces: false,
                output_aevs: false,
            };
            compute_ani(&elements, pos, &cfg, &models).unwrap().energy
        };

        let h = 1e-5;
        for atom in 0..positions.len() {
            for d in 0..3 {
                let mut plus = positions.clone();
                let mut minus = positions.clone();
                plus[atom][d] += h;
                minus[atom][d] -= h;
                let d_e = (energy_at(&plus) - energy_at(&minus)) / (2.0 * h);
                // F = -dE/dR
                let analytic = forces[atom][d];
                let err = (analytic + d_e).abs();
                assert!(
                    err < 1e-6,
                    "force mismatch atom {atom} dim {d}: analytic {analytic:.8e} vs \
                     finite-difference {:.8e} (err {err:.3e})",
                    -d_e
                );
            }
        }
    }

    #[test]
    fn test_forces_match_finite_difference_methane() {
        // A second, low-symmetry geometry: the analytic gradient must still
        // reproduce the numerical gradient of the total energy.
        use super::super::api::{compute_ani, AniConfig};

        let elements = vec![6u8, 1, 1, 1, 1];
        let positions = vec![
            [0.0, 0.0, 0.0],
            [0.63, 0.63, 0.63],
            [-0.63, 0.63, 0.63],
            [-0.63, -0.63, 0.63],
            [0.63, -0.63, 0.63],
        ];
        let params = default_ani2x_params();
        let mut models = HashMap::new();
        models.insert(6u8, make_test_model(params.total_aev_length()));
        models.insert(1u8, make_test_model(params.total_aev_length()));

        let cl = CellList::new(&positions, params.radial_cutoff);
        let neighbors = cl.find_neighbors(&positions);
        let forces = compute_forces(&elements, &positions, &neighbors, &params, &models);

        let energy_at = |pos: &[[f64; 3]]| -> f64 {
            let cfg = AniConfig {
                cutoff: params.radial_cutoff,
                compute_forces: false,
                output_aevs: false,
            };
            compute_ani(&elements, pos, &cfg, &models).unwrap().energy
        };

        let h = 1e-5;
        for atom in 0..positions.len() {
            for d in 0..3 {
                let mut plus = positions.clone();
                let mut minus = positions.clone();
                plus[atom][d] += h;
                minus[atom][d] -= h;
                let d_e = (energy_at(&plus) - energy_at(&minus)) / (2.0 * h);
                let analytic = forces[atom][d];
                let err = (analytic + d_e).abs();
                assert!(
                    err < 1e-6,
                    "force mismatch atom {atom} dim {d}: analytic {analytic:.8e} vs \
                     finite-difference {:.8e} (err {err:.3e})",
                    -d_e
                );
            }
        }
    }
}
