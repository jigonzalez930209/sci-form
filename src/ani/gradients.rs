//! Analytical gradient computation (forces) for ANI potentials.
//!
//! Forces = −∇E are computed by backpropagating through:
//! 1. NN output → NN input (dE/dAEV via nn.backward)
//! 2. AEV → Cartesian coordinates (chain rule through symmetry functions)
//!
//! The Cartesian gradient of the radial AEV component is:
//! $$\frac{\partial G_R}{\partial \vec{r}_i} = \sum_j \left[
//!   e^{-\eta(R_{ij}-R_s)^2} \left( -2\eta(R_{ij}-R_s) f_c + f_c' \right)
//! \right] \frac{\vec{r}_{ij}}{R_{ij}}$$

use super::aev_params::{species_index, AevParams};
use super::cutoff::{cosine_cutoff, cosine_cutoff_deriv};
use super::neighbor::NeighborPair;
use super::nn::FeedForwardNet;
use nalgebra::DVector;
use std::collections::HashMap;

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
    let aev_len = params.total_aev_length();

    // Build per-atom neighbor lists
    let mut atom_neighbors: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for np in neighbors {
        let d = np.dist_sq.sqrt();
        atom_neighbors[np.i].push((np.j, d));
        atom_neighbors[np.j].push((np.i, d));
    }

    // Forward pass: compute AEVs and NN gradients (dE/dAEV)
    let mut aev_grads: Vec<Vec<f64>> = vec![vec![0.0; aev_len]; n];
    for i in 0..n {
        let si = match species_index(elements[i]) {
            Some(s) => s,
            None => continue,
        };
        let _ = si;
        let aev = compute_atom_aev(i, elements, positions, &atom_neighbors[i], params);
        let net = match models.get(&elements[i]) {
            Some(net) => net,
            None => continue,
        };
        let input = DVector::from_vec(aev);
        let grad = net.backward(&input);
        aev_grads[i] = grad.as_slice().to_vec();
    }

    // Backprop: AEV gradients → Cartesian forces
    let mut forces = vec![[0.0f64; 3]; n];
    for i in 0..n {
        if species_index(elements[i]).is_none() {
            continue;
        }
        backprop_radial_forces(
            i,
            elements,
            positions,
            &atom_neighbors[i],
            params,
            &aev_grads[i],
            &mut forces,
        );
    }

    // Negate: F = -∇E
    for f in &mut forces {
        f[0] = -f[0];
        f[1] = -f[1];
        f[2] = -f[2];
    }
    forces
}

fn compute_atom_aev(
    _i: usize,
    elements: &[u8],
    _positions: &[[f64; 3]],
    neighbors_i: &[(usize, f64)],
    params: &AevParams,
) -> Vec<f64> {
    let aev_len = params.total_aev_length();
    let mut aev = vec![0.0f64; aev_len];
    let rad_len = params.radial_length();

    for &(j, rij) in neighbors_i {
        if rij >= params.radial_cutoff {
            continue;
        }
        let sj = match species_index(elements[j]) {
            Some(s) => s,
            None => continue,
        };
        let fc = cosine_cutoff(rij, params.radial_cutoff);
        let offset = sj * rad_len;
        let mut k = 0;
        for eta in &params.radial_eta {
            for rs in &params.radial_rs {
                let dr = rij - rs;
                aev[offset + k] += (-eta * dr * dr).exp() * fc;
                k += 1;
            }
        }
    }
    aev
}

fn backprop_radial_forces(
    i: usize,
    elements: &[u8],
    positions: &[[f64; 3]],
    neighbors_i: &[(usize, f64)],
    params: &AevParams,
    aev_grad: &[f64],
    forces: &mut [[f64; 3]],
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
                    forces[i][d] -= chain * dir[d];
                    forces[j][d] += chain * dir[d];
                }
                k += 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::aev_params::default_ani2x_params;
    use super::super::neighbor::CellList;
    use super::super::weights::make_test_model;
    use super::*;

    #[test]
    fn test_forces_sum_near_zero() {
        // Newton's third law: total force should be ~zero
        let elements = [8u8, 1, 1];
        let positions = [
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469],
        ];
        let params = default_ani2x_params();
        let cl = CellList::new(&positions, params.radial_cutoff);
        let neighbors = cl.find_neighbors(&positions);

        let aev_len = params.total_aev_length();
        let mut models = HashMap::new();
        models.insert(8u8, make_test_model(aev_len));
        models.insert(1u8, make_test_model(aev_len));

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
}
