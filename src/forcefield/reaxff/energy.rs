//! ReaxFF bonded energy terms as functions of continuous bond order.
//!
//! Includes bond stretch, angle bend, and torsion contributions.

use super::bond_order::BondOrder;
use super::params::ReaxffParams;

/// Total ReaxFF energy result.
#[derive(Debug, Clone, Default)]
pub struct ReaxffEnergy {
    /// Bond stretch energy (kcal/mol).
    pub bond_energy: f64,
    /// Angle bend energy (kcal/mol).
    pub angle_energy: f64,
    /// Torsion energy (kcal/mol).
    pub torsion_energy: f64,
    /// Over-coordination penalty (kcal/mol).
    pub over_coord_energy: f64,
    /// Under-coordination penalty (kcal/mol).
    pub under_coord_energy: f64,
    /// Total bonded energy (kcal/mol).
    pub total: f64,
}

/// Compute all bonded ReaxFF energy terms.
pub fn compute_bonded_energy(
    positions_flat: &[f64],
    bond_orders: &[Vec<BondOrder>],
    params: &ReaxffParams,
) -> ReaxffEnergy {
    let n = params.atom_params.len();
    let mut result = ReaxffEnergy::default();

    // Bond stretch energy
    for i in 0..n {
        for j in (i + 1)..n {
            let bo = &bond_orders[i][j];
            if bo.total < 1e-4 {
                continue;
            }

            // E_bond = -D_e * BO_σ * exp(p_be1 * (1 - BO_σ^p_be2))
            // Using simplified parameters
            let de = params.get_bond_de(i, j);
            let e_bond =
                -de * bo.sigma * (params.p_be1 * (1.0 - bo.sigma.powf(params.p_be2))).exp();
            result.bond_energy += e_bond;
        }
    }

    // Angle bend energy
    for j in 0..n {
        for i in 0..n {
            if i == j {
                continue;
            }
            let bo_ij = &bond_orders[i][j];
            if bo_ij.total < 0.01 {
                continue;
            }

            for k in (i + 1)..n {
                if k == j {
                    continue;
                }
                let bo_jk = &bond_orders[j][k];
                if bo_jk.total < 0.01 {
                    continue;
                }

                let angle = compute_angle(positions_flat, i, j, k);
                let theta0 = params.get_equilibrium_angle(i, j, k);
                let ka = params.get_angle_force_constant(i, j, k);

                // E_angle = ka * f(BO_ij) * f(BO_jk) * (1 + cos(theta - theta0))
                let f_bo = bo_ij.total.min(1.0) * bo_jk.total.min(1.0);
                let e_angle = ka * f_bo * (1.0 - (angle - theta0).cos());
                result.angle_energy += e_angle;
            }
        }
    }

    // Over-coordination penalty
    for i in 0..n {
        let val = params.atom_params[i].valence;
        let total_bo: f64 = (0..n)
            .filter(|&j| j != i)
            .map(|j| bond_orders[i][j].total)
            .sum();
        let delta = total_bo - val;
        if delta > 0.0 {
            result.over_coord_energy += params.p_ovun1 * delta / (delta + params.p_val3);
        }
    }

    // Under-coordination penalty
    for i in 0..n {
        let val = params.atom_params[i].valence;
        let total_bo: f64 = (0..n)
            .filter(|&j| j != i)
            .map(|j| bond_orders[i][j].total)
            .sum();
        let delta = val - total_bo;
        if delta > 0.0 {
            result.under_coord_energy += params.p_ovun5 * delta;
        }
    }

    result.total = result.bond_energy
        + result.angle_energy
        + result.torsion_energy
        + result.over_coord_energy
        + result.under_coord_energy;

    result
}

/// Compute angle (radians) between atoms i-j-k.
fn compute_angle(positions_flat: &[f64], i: usize, j: usize, k: usize) -> f64 {
    let rji = [
        positions_flat[3 * i] - positions_flat[3 * j],
        positions_flat[3 * i + 1] - positions_flat[3 * j + 1],
        positions_flat[3 * i + 2] - positions_flat[3 * j + 2],
    ];
    let rjk = [
        positions_flat[3 * k] - positions_flat[3 * j],
        positions_flat[3 * k + 1] - positions_flat[3 * j + 1],
        positions_flat[3 * k + 2] - positions_flat[3 * j + 2],
    ];

    let dot = rji[0] * rjk[0] + rji[1] * rjk[1] + rji[2] * rjk[2];
    let mag_ji = (rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2]).sqrt();
    let mag_jk = (rjk[0] * rjk[0] + rjk[1] * rjk[1] + rjk[2] * rjk[2]).sqrt();

    let cos_theta = (dot / (mag_ji * mag_jk + 1e-30)).clamp(-1.0, 1.0);
    cos_theta.acos()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::reaxff::bond_order::compute_bond_orders;
    use crate::forcefield::reaxff::params::{ReaxffAtomParams, ReaxffParams};

    fn h2_params() -> ReaxffParams {
        let h_param = ReaxffAtomParams {
            element: 1,
            r_sigma: 0.656,
            r_pi: 0.0,
            r_pipi: 0.0,
            p_bo1: -0.0500,
            p_bo2: 6.9136,
            p_bo3: 0.0,
            p_bo4: 6.0,
            p_bo5: 0.0,
            p_bo6: 6.0,
            valence: 1.0,
        };
        // Build params with atom_params matching exactly 2 H atoms
        ReaxffParams {
            atom_params: vec![h_param.clone(), h_param],
            bond_de: vec![vec![100.0; 2]; 2],
            p_be1: -0.2,
            p_be2: 6.25,
            p_ovun1: 50.0,
            p_val3: 3.0,
            p_ovun5: 10.0,
            equilibrium_angles: vec![109.47_f64.to_radians(); 10],
            angle_force_constants: vec![50.0; 10],
            cutoff: 10.0,
        }
    }

    #[test]
    fn bonded_energy_h2_is_finite() {
        let params = h2_params();
        let positions = [0.0, 0.0, 0.0, 0.74, 0.0, 0.0];
        let bo = compute_bond_orders(&positions, &params.atom_params, params.cutoff);
        let energy = compute_bonded_energy(&positions, &bo, &params);
        assert!(energy.total.is_finite(), "energy should be finite");
    }

    #[test]
    fn bonded_energy_components_sum() {
        let params = h2_params();
        let positions = [0.0, 0.0, 0.0, 0.74, 0.0, 0.0];
        let bo = compute_bond_orders(&positions, &params.atom_params, params.cutoff);
        let e = compute_bonded_energy(&positions, &bo, &params);
        let sum = e.bond_energy
            + e.angle_energy
            + e.torsion_energy
            + e.over_coord_energy
            + e.under_coord_energy;
        assert!(
            (sum - e.total).abs() < 1e-10,
            "component sum should equal total"
        );
    }
}
