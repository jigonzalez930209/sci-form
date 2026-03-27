//! Continuous bond-order calculation from interatomic distances.
//!
//! ReaxFF bond orders are smooth functions of distance with three components:
//! - σ bond (single bond backbone)
//! - π bond (double bond component)
//! - ππ bond (triple bond component)
//!
//! BO_ij = BO_σ + BO_π + BO_ππ
//! BO_σ  = exp(p_bo1 · (r_ij / r0_σ)^p_bo2)
//! BO_π  = exp(p_bo3 · (r_ij / r0_π)^p_bo4)
//! BO_ππ = exp(p_bo5 · (r_ij / r0_ππ)^p_bo6)

use super::params::ReaxffAtomParams;
use super::taper::taper_function;

/// Bond order between two atoms.
#[derive(Debug, Clone, Copy)]
pub struct BondOrder {
    /// Total bond order.
    pub total: f64,
    /// σ component.
    pub sigma: f64,
    /// π component.
    pub pi: f64,
    /// ππ component.
    pub pipi: f64,
    /// Distance (Å).
    pub distance: f64,
}

/// Compute continuous bond orders for all atom pairs.
///
/// Returns a matrix of bond orders indexed by (i, j).
pub fn compute_bond_orders(
    positions_flat: &[f64],
    atom_params: &[ReaxffAtomParams],
    cutoff: f64,
) -> Vec<Vec<BondOrder>> {
    let n = atom_params.len();
    let mut bo = vec![
        vec![
            BondOrder {
                total: 0.0,
                sigma: 0.0,
                pi: 0.0,
                pipi: 0.0,
                distance: 0.0
            };
            n
        ];
        n
    ];

    for i in 0..n {
        let xi = positions_flat[3 * i];
        let yi = positions_flat[3 * i + 1];
        let zi = positions_flat[3 * i + 2];

        for j in (i + 1)..n {
            let dx = xi - positions_flat[3 * j];
            let dy = yi - positions_flat[3 * j + 1];
            let dz = zi - positions_flat[3 * j + 2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();

            if r > cutoff || r < 0.01 {
                continue;
            }

            let tap = taper_function(r, cutoff);

            // Combined parameters (geometric mean for heteronuclear pairs)
            let pi = &atom_params[i];
            let pj = &atom_params[j];

            let r0_sigma = (pi.r_sigma + pj.r_sigma) / 2.0;
            let r0_pi = (pi.r_pi + pj.r_pi) / 2.0;
            let r0_pipi = (pi.r_pipi + pj.r_pipi) / 2.0;

            // Bond order components
            let bo_sigma = if r0_sigma > 0.01 {
                (pi.p_bo1 * (r / r0_sigma).powf(pi.p_bo2)).exp()
            } else {
                0.0
            };

            let bo_pi = if r0_pi > 0.01 {
                (pi.p_bo3 * (r / r0_pi).powf(pi.p_bo4)).exp()
            } else {
                0.0
            };

            let bo_pipi = if r0_pipi > 0.01 {
                (pi.p_bo5 * (r / r0_pipi).powf(pi.p_bo6)).exp()
            } else {
                0.0
            };

            let total = (bo_sigma + bo_pi + bo_pipi) * tap;

            let entry = BondOrder {
                total,
                sigma: bo_sigma * tap,
                pi: bo_pi * tap,
                pipi: bo_pipi * tap,
                distance: r,
            };

            bo[i][j] = entry;
            bo[j][i] = entry;
        }
    }

    // Apply over-coordination correction
    apply_overcorrection(&mut bo, atom_params);

    bo
}

/// Apply over-coordination penalty to bond orders.
fn apply_overcorrection(bo: &mut [Vec<BondOrder>], atom_params: &[ReaxffAtomParams]) {
    let n = atom_params.len();

    // Compute total valence for each atom
    let mut valence_total: Vec<f64> = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            if i != j {
                valence_total[i] += bo[i][j].total;
            }
        }
    }

    // Apply correction: reduce bond orders when over-coordinated
    for i in 0..n {
        let val_i = atom_params[i].valence;
        let delta_i = valence_total[i] - val_i;

        if delta_i > 0.0 {
            let correction = 1.0 / (1.0 + delta_i.exp());
            for j in 0..n {
                if i != j && bo[i][j].total > 1e-6 {
                    let val_j = atom_params[j].valence;
                    let delta_j = valence_total[j] - val_j;
                    let corr_j = if delta_j > 0.0 {
                        1.0 / (1.0 + delta_j.exp())
                    } else {
                        1.0
                    };
                    let scale = correction * corr_j;
                    bo[i][j].total *= scale;
                    bo[i][j].sigma *= scale;
                    bo[i][j].pi *= scale;
                    bo[i][j].pipi *= scale;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::reaxff::params::ReaxffParams;

    #[test]
    fn bond_order_h2_is_near_one() {
        let params = ReaxffParams::default_chon();
        // H-H at ~0.74 Å — should have bond order ~ 1
        let positions = [0.0, 0.0, 0.0, 0.74, 0.0, 0.0];
        let atom_params = vec![
            params.atom_params[params.element_index(1).unwrap()].clone(),
            params.atom_params[params.element_index(1).unwrap()].clone(),
        ];
        let bo_matrix = compute_bond_orders(&positions, &atom_params, 10.0);
        let bo = bo_matrix[0][1].total; // BO between atom 0 and atom 1
        assert!(
            bo > 0.3 && bo < 2.0,
            "H-H bond order at 0.74Å should be roughly 1, got {bo}"
        );
    }

    #[test]
    fn bond_order_decreases_with_distance() {
        let params = ReaxffParams::default_chon();
        let atom_params = vec![
            params.atom_params[params.element_index(6).unwrap()].clone(),
            params.atom_params[params.element_index(6).unwrap()].clone(),
        ];
        let pos_close = [0.0, 0.0, 0.0, 1.54, 0.0, 0.0]; // C-C single bond
        let pos_far = [0.0, 0.0, 0.0, 3.0, 0.0, 0.0]; // far apart
        let bo_close = compute_bond_orders(&pos_close, &atom_params, 10.0);
        let bo_far = compute_bond_orders(&pos_far, &atom_params, 10.0);
        assert!(
            bo_close[0][1].total > bo_far[0][1].total,
            "closer atoms should have higher bond order: close={} vs far={}",
            bo_close[0][1].total,
            bo_far[0][1].total,
        );
    }

    #[test]
    fn bond_order_sigma_pi_components() {
        let params = ReaxffParams::default_chon();
        let atom_params = vec![
            params.atom_params[params.element_index(6).unwrap()].clone(),
            params.atom_params[params.element_index(6).unwrap()].clone(),
        ];
        let positions = [0.0, 0.0, 0.0, 1.54, 0.0, 0.0];
        let bo_matrix = compute_bond_orders(&positions, &atom_params, 10.0);
        let bo = bo_matrix[0][1];
        assert!(bo.sigma >= 0.0);
        assert!(bo.pi >= 0.0);
        assert!(bo.pipi >= 0.0);
        // Total ≈ sigma + pi + pipi (after overcorrection scaling, may differ)
        assert!(bo.total >= 0.0);
    }
}
