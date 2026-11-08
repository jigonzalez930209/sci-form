//! Analytical gradients for the ReaxFF energy model.
//!
//! The energy is the one defined by [`super::energy::compute_bonded_energy`]
//! plus [`super::nonbonded::compute_nonbonded_energy_excluding`] and the EEM
//! charges from [`super::eem::solve_eem`].  Every term is differentiated
//! explicitly:
//!
//! * bond orders — chain rule through the exponential/taper form and the
//!   over-coordination scaling;
//! * bond, angle, over- and under-coordination energies — through the bond
//!   orders and the angle geometry;
//! * buffered van der Waals and shielded Coulomb — radial derivatives;
//! * EEM charges — `dq = -A⁻¹ (dA) q` using the cached inverse of the
//!   augmented charge-equilibration system.
//!
//! The implementation is validated against central differences in
//! `tests/alpha_reaxff_validation.rs`.

use super::eem::{default_eem_params, solve_eem_with_inverse};
use super::params::{ReaxffAtomParams, ReaxffParams};
use super::taper::{taper_derivative, taper_function};

const BO_CUT: f64 = 0.3;
const VDW_GAMMA: f64 = 1.0;
const COUL_GAMMA: f64 = 0.5;
const COULOMB_FACTOR: f64 = 332.0638;

/// Per-pair bond-order values and radial derivatives.
struct BondOrderData {
    /// Corrected total bond order matrix (the one the energy uses).
    total: Vec<Vec<f64>>,
    /// Corrected σ bond order matrix.
    sigma: Vec<Vec<f64>>,
    /// Raw total bond order before over-coordination scaling.
    raw_total: Vec<Vec<f64>>,
    /// Raw σ bond order before over-coordination scaling.
    raw_sigma: Vec<Vec<f64>>,
    /// Radial derivative of the raw total bond order.
    d_total: Vec<Vec<f64>>,
    /// Radial derivative of the raw σ bond order.
    d_sigma: Vec<Vec<f64>>,
    /// Over-coordination sigmoid for each atom.
    correction: Vec<f64>,
    /// Whether the atom has positive over-coordination.
    active: Vec<bool>,
    /// Multiplicative factor applied to entries in row `i`.
    row_factor: Vec<f64>,
    /// Derivative of the sigmoid with respect to `Δ_i`.
    d_correction: Vec<f64>,
}

fn atom_params_for(params: &ReaxffParams, elements: &[u8]) -> Vec<ReaxffAtomParams> {
    elements
        .iter()
        .map(|&z| {
            params
                .element_index(z)
                .map(|index| params.atom_params[index].clone())
                .unwrap_or_else(|| params.atom_params.last().cloned().unwrap())
        })
        .collect()
}

fn local_params(params: &ReaxffParams, atom_params: &[ReaxffAtomParams]) -> ReaxffParams {
    let n = atom_params.len();
    ReaxffParams {
        atom_params: atom_params.to_vec(),
        bond_de: (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        let gi = params.element_index(atom_params[i].element);
                        let gj = params.element_index(atom_params[j].element);
                        match (gi, gj) {
                            (Some(a), Some(b)) => params.get_bond_de(a, b),
                            _ => 100.0,
                        }
                    })
                    .collect()
            })
            .collect(),
        ..params.clone()
    }
}

/// Value and radial derivative of one exponential/taper bond-order component.
fn component(r: f64, tap: f64, dtap: f64, b1: f64, b2: f64, r0: f64) -> (f64, f64) {
    if r0 <= 0.01 || r <= 0.0 {
        return (0.0, 0.0);
    }
    let x = r / r0;
    let power = x.powf(b2);
    let exponential = (b1 * power).exp();
    let value = exponential * tap;
    let dvalue = value * b1 * b2 * power / r + exponential * dtap;
    (value, dvalue)
}

fn build_bond_orders(
    positions: &[f64],
    atom_params: &[ReaxffAtomParams],
    params: &ReaxffParams,
    cutoff: f64,
) -> BondOrderData {
    let n = atom_params.len();
    let mut raw_total = vec![vec![0.0; n]; n];
    let mut raw_sigma = vec![vec![0.0; n]; n];
    let mut d_total = vec![vec![0.0; n]; n];
    let mut d_sigma = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[3 * i] - positions[3 * j];
            let dy = positions[3 * i + 1] - positions[3 * j + 1];
            let dz = positions[3 * i + 2] - positions[3 * j + 2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r > cutoff || r < 0.01 {
                continue;
            }
            let tap = taper_function(r, cutoff);
            let dtap = taper_derivative(r, cutoff);
            let pi = &atom_params[i];
            let pj = &atom_params[j];

            let (r0_sigma, b1, b2) = match params.pair_params(pi.element, pj.element) {
                Some(pair) => (pair.r_sigma, pair.p_bo1, pair.p_bo2),
                None => ((pi.r_sigma + pj.r_sigma) / 2.0, pi.p_bo1, pi.p_bo2),
            };
            let (sigma, dsigma) = component(r, tap, dtap, b1, b2, r0_sigma);
            let (pi_part, dpi) =
                component(r, tap, dtap, pi.p_bo3, pi.p_bo4, (pi.r_pi + pj.r_pi) / 2.0);
            let (pipi, dpipi) = component(
                r,
                tap,
                dtap,
                pi.p_bo5,
                pi.p_bo6,
                (pi.r_pipi + pj.r_pipi) / 2.0,
            );

            let total = sigma + pi_part + pipi;
            let dtotal = dsigma + dpi + dpipi;
            raw_sigma[i][j] = sigma;
            raw_sigma[j][i] = sigma;
            raw_total[i][j] = total;
            raw_total[j][i] = total;
            d_sigma[i][j] = dsigma;
            d_sigma[j][i] = dsigma;
            d_total[i][j] = dtotal;
            d_total[j][i] = dtotal;
        }
    }

    let mut correction = vec![1.0; n];
    let mut active = vec![false; n];
    let mut d_correction = vec![0.0; n];
    for i in 0..n {
        let valence: f64 = (0..n).filter(|&j| j != i).map(|j| raw_total[i][j]).sum();
        let delta = valence - atom_params[i].valence;
        if delta > 0.0 {
            active[i] = true;
            let c = if delta > 700.0 {
                0.0
            } else {
                1.0 / (1.0 + delta.exp())
            };
            correction[i] = c;
            d_correction[i] = -c * (1.0 - c);
        }
    }

    let mut total = vec![vec![0.0; n]; n];
    let mut sigma = vec![vec![0.0; n]; n];
    let row_factor: Vec<f64> = (0..n)
        .map(|i| if active[i] { correction[i] } else { 1.0 })
        .collect();
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            if active[i] {
                let scale = correction[i] * row_factor[j];
                total[i][j] = raw_total[i][j] * scale;
                sigma[i][j] = raw_sigma[i][j] * scale;
            } else {
                total[i][j] = raw_total[i][j];
                sigma[i][j] = raw_sigma[i][j];
            }
        }
    }

    BondOrderData {
        total,
        sigma,
        raw_total,
        raw_sigma,
        d_total,
        d_sigma,
        correction,
        active,
        row_factor,
        d_correction,
    }
}

fn angle_geometry(
    positions: &[f64],
    i: usize,
    j: usize,
    k: usize,
) -> (f64, f64, [f64; 3], [f64; 3], f64, f64) {
    let rji = [
        positions[3 * i] - positions[3 * j],
        positions[3 * i + 1] - positions[3 * j + 1],
        positions[3 * i + 2] - positions[3 * j + 2],
    ];
    let rjk = [
        positions[3 * k] - positions[3 * j],
        positions[3 * k + 1] - positions[3 * j + 1],
        positions[3 * k + 2] - positions[3 * j + 2],
    ];
    let d_ji = (rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2])
        .sqrt()
        .max(1e-12);
    let d_jk = (rjk[0] * rjk[0] + rjk[1] * rjk[1] + rjk[2] * rjk[2])
        .sqrt()
        .max(1e-12);
    let cos_theta =
        ((rji[0] * rjk[0] + rji[1] * rjk[1] + rji[2] * rjk[2]) / (d_ji * d_jk)).clamp(-1.0, 1.0);
    let theta = cos_theta.acos();
    (theta, cos_theta, rji, rjk, d_ji, d_jk)
}

fn add_angle_gradient(
    coords: &[f64],
    grad: &mut [f64],
    i: usize,
    j: usize,
    k: usize,
    de_dtheta: f64,
) {
    let (_theta, cos_theta, rji, rjk, d_ji, d_jk) = angle_geometry(coords, i, j, k);
    let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1e-12);
    let pref = de_dtheta * (-1.0 / sin_theta);
    for dim in 0..3 {
        let term_i = pref * (rjk[dim] / (d_ji * d_jk) - cos_theta * rji[dim] / (d_ji * d_ji));
        let term_k = pref * (rji[dim] / (d_ji * d_jk) - cos_theta * rjk[dim] / (d_jk * d_jk));
        grad[3 * i + dim] += term_i;
        grad[3 * k + dim] += term_k;
        grad[3 * j + dim] -= term_i + term_k;
    }
}

/// Compute the ReaxFF energy and its analytical Cartesian gradient.
pub fn compute_reaxff_gradient_analytic(
    positions: &[f64],
    elements: &[u8],
    params: &ReaxffParams,
) -> Result<(f64, Vec<f64>), String> {
    let n = elements.len();
    if positions.len() != 3 * n {
        return Err(format!(
            "positions length {} != 3 * atoms {n}",
            positions.len()
        ));
    }
    let atom_params = atom_params_for(params, elements);
    let local = local_params(params, &atom_params);
    let bo = build_bond_orders(positions, &atom_params, params, params.cutoff);

    let mut grad = vec![0.0f64; 3 * n];
    let mut gt = vec![vec![0.0f64; n]; n];
    let mut gs = vec![vec![0.0f64; n]; n];
    let mut energy = 0.0;

    // ── bond stretch ──
    for i in 0..n {
        for j in (i + 1)..n {
            if bo.total[i][j] < 1e-4 {
                continue;
            }
            let s = bo.sigma[i][j];
            let (de, be1, be2) =
                match local.pair_params(atom_params[i].element, atom_params[j].element) {
                    Some(pair) => (pair.de_sigma, pair.p_be1, pair.p_be2),
                    None => (local.get_bond_de(i, j), local.p_be1, local.p_be2),
                };
            let exponent = (be1 * (1.0 - s.powf(be2))).exp();
            energy += -de * s * exponent;
            let de_ds = -de * exponent * (1.0 - be1 * be2 * s.powf(be2));
            gs[i][j] += de_ds;
        }
    }

    // ── angle bend ──
    for j in 0..n {
        for i in 0..n {
            if i == j {
                continue;
            }
            let bo_ij = bo.total[i][j];
            if bo_ij < 0.01 {
                continue;
            }
            for k in (i + 1)..n {
                if k == j {
                    continue;
                }
                let bo_jk = bo.total[j][k];
                if bo_jk < 0.01 {
                    continue;
                }
                let (theta, cos_theta, _, _, _, _) = angle_geometry(positions, i, j, k);
                let theta0 = local.get_equilibrium_angle(i, j, k);
                let ka = local.get_angle_force_constant(i, j, k);
                let clamp_ij = bo_ij.min(1.0);
                let clamp_jk = bo_jk.min(1.0);
                let d_ij = if bo_ij < 1.0 { 1.0 } else { 0.0 };
                let d_jk = if bo_jk < 1.0 { 1.0 } else { 0.0 };
                let shape = 1.0 - (theta - theta0).cos();
                let force = ka * clamp_ij * clamp_jk;
                energy += force * shape;

                gt[i][j] += ka * d_ij * clamp_jk * shape;
                gt[j][k] += ka * clamp_ij * d_jk * shape;

                let de_dtheta = force * (theta - theta0).sin();
                add_angle_gradient(positions, &mut grad, i, j, k, de_dtheta);
                let _ = cos_theta;
            }
        }
    }

    // ── over- and under-coordination ──
    let mut delta = vec![0.0f64; n];
    for i in 0..n {
        let total_bo: f64 = (0..n).filter(|&j| j != i).map(|j| bo.total[i][j]).sum();
        delta[i] = total_bo - atom_params[i].valence;
        if delta[i] > 0.0 {
            energy += local.p_ovun1 * delta[i] / (delta[i] + local.p_val3);
            let de = local.p_ovun1 * local.p_val3 / (delta[i] + local.p_val3).powi(2);
            for j in 0..n {
                if j != i {
                    gt[i][j] += de;
                }
            }
        }
        let under = atom_params[i].valence - total_bo;
        if under > 0.0 {
            energy += local.p_ovun5 * under;
            for j in 0..n {
                if j != i {
                    gt[i][j] -= local.p_ovun5;
                }
            }
        }
    }

    // ── non-bonded: van der Waals + shielded Coulomb with EEM charges ──
    let eem_params: Vec<_> = elements.iter().map(|&z| default_eem_params(z)).collect();
    let eem = solve_eem_with_inverse(positions, &eem_params, 0.0)?;
    let charges = &eem.charges;
    let inverse = &eem.inverse;
    let q_aug = &eem.augmented;

    // ∂E_coulomb/∂q_k for the fixed-geometry Coulomb energy.
    let mut de_dq = vec![0.0f64; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[3 * i] - positions[3 * j];
            let dy = positions[3 * i + 1] - positions[3 * j + 1];
            let dz = positions[3 * i + 2] - positions[3 * j + 2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r > params.cutoff || r < 0.01 {
                continue;
            }
            let tap = taper_function(r, params.cutoff);
            let r_shield = (r.powi(3) + COUL_GAMMA.powi(3)).powf(1.0 / 3.0);
            let coefficient = COULOMB_FACTOR * tap / r_shield;
            de_dq[i] += coefficient * charges[j];
            de_dq[j] += coefficient * charges[i];
        }
    }

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[3 * i] - positions[3 * j];
            let dy = positions[3 * i + 1] - positions[3 * j + 1];
            let dz = positions[3 * i + 2] - positions[3 * j + 2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r > params.cutoff || r < 0.01 {
                continue;
            }
            let tap = taper_function(r, params.cutoff);
            let dtap = taper_derivative(r, params.cutoff);

            // Shielded Coulomb.
            let r_shield = (r.powi(3) + COUL_GAMMA.powi(3)).powf(1.0 / 3.0);
            let e_coul = COULOMB_FACTOR * charges[i] * charges[j] / r_shield;
            energy += tap * e_coul;
            let d_e_coul = -COULOMB_FACTOR * charges[i] * charges[j] * r.powi(2)
                / (r_shield * r_shield * r_shield * r_shield)
                * tap
                + dtap * e_coul;

            // EEM charge response: dq = -A⁻¹ dA q.
            let gamma_ij = (eem_params[i].gamma + eem_params[j].gamma) / 2.0;
            let shielded = 1.0 / (r.powi(3) + gamma_ij.powi(3)).powf(1.0 / 3.0);
            let coulomb = 14.4 * shielded;
            let d_coulomb = -coulomb * r * r / (r.powi(3) + gamma_ij.powi(3));
            let d_q = -d_coulomb * (q_aug[j] * inverse.column(i) + q_aug[i] * inverse.column(j));
            let mut charge_response = 0.0;
            for k in 0..n {
                charge_response += de_dq[k] * d_q[k];
            }
            let radial = d_e_coul + charge_response;

            // van der Waals (excluded for strongly σ-bonded pairs).
            let bonded = bo.sigma[i][j] > BO_CUT;
            let mut radial_vdw = 0.0;
            if !bonded {
                let (eps, r0, _) = vdw_params(elements[i], elements[j]);
                let r_shielded = (r.powi(7) + VDW_GAMMA.powi(7)).powf(1.0 / 7.0);
                let sigma = r0 / 2.0f64.powf(1.0 / 6.0);
                let sr6 = (sigma / r_shielded).powi(6);
                let e_lj = 4.0 * eps * (sr6 * sr6 - sr6);
                energy += tap * e_lj;
                let dr_shielded = r.powi(6) / r_shielded.powi(6);
                let de_lj = 4.0 * eps * 6.0 * sr6 * (1.0 - 2.0 * sr6) / r_shielded;
                radial_vdw = tap * de_lj * dr_shielded + dtap * e_lj;
            }

            let radial = radial + radial_vdw;
            let inv_r = 1.0 / r;
            for dim in 0..3 {
                let component = radial * [dx, dy, dz][dim] * inv_r;
                grad[3 * i + dim] += component;
                grad[3 * j + dim] -= component;
            }
        }
    }

    // ── bond-order chain rule ──
    for a in 0..n {
        for b in (a + 1)..n {
            let d_total = bo.d_total[a][b];
            let d_sigma = bo.d_sigma[a][b];
            if d_total == 0.0 && d_sigma == 0.0 {
                continue;
            }

            let scale_ab = if bo.active[a] {
                bo.correction[a] * bo.row_factor[b]
            } else {
                1.0
            };
            let scale_ba = if bo.active[b] {
                bo.correction[b] * bo.row_factor[a]
            } else {
                1.0
            };

            let mut de_dr = d_total * (gt[a][b] * scale_ab + gt[b][a] * scale_ba)
                + d_sigma * (gs[a][b] * scale_ab + gs[b][a] * scale_ba);

            // Derivative of the over-coordination scaling.
            let mut row_a = 0.0;
            let mut col_a = 0.0;
            let mut row_b = 0.0;
            let mut col_b = 0.0;
            for j in 0..n {
                if j != a {
                    row_a += gt[a][j] * bo.raw_total[a][j] * bo.row_factor[j]
                        + gs[a][j] * bo.raw_sigma[a][j] * bo.row_factor[j];
                    col_a += if bo.active[j] {
                        gt[j][a] * bo.raw_total[j][a] * bo.correction[j]
                            + gs[j][a] * bo.raw_sigma[j][a] * bo.correction[j]
                    } else {
                        0.0
                    };
                }
                if j != b {
                    row_b += gt[b][j] * bo.raw_total[b][j] * bo.row_factor[j]
                        + gs[b][j] * bo.raw_sigma[b][j] * bo.row_factor[j];
                    col_b += if bo.active[j] {
                        gt[j][b] * bo.raw_total[j][b] * bo.correction[j]
                            + gs[j][b] * bo.raw_sigma[j][b] * bo.correction[j]
                    } else {
                        0.0
                    };
                }
            }
            let u_a = if bo.active[a] {
                bo.d_correction[a] * d_total
            } else {
                0.0
            };
            let u_b = if bo.active[b] {
                bo.d_correction[b] * d_total
            } else {
                0.0
            };
            de_dr += u_a * (row_a + col_a) + u_b * (row_b + col_b);

            let dx = positions[3 * a] - positions[3 * b];
            let dy = positions[3 * a + 1] - positions[3 * b + 1];
            let dz = positions[3 * a + 2] - positions[3 * b + 2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-12);
            let pref = de_dr / r;
            grad[3 * a] += pref * dx;
            grad[3 * a + 1] += pref * dy;
            grad[3 * a + 2] += pref * dz;
            grad[3 * b] -= pref * dx;
            grad[3 * b + 1] -= pref * dy;
            grad[3 * b + 2] -= pref * dz;
        }
    }

    Ok((energy, grad))
}

fn vdw_params(z1: u8, z2: u8) -> (f64, f64, f64) {
    let (e1, r1) = element_vdw(z1);
    let (e2, r2) = element_vdw(z2);
    ((e1 * e2).sqrt(), (r1 + r2) / 2.0, VDW_GAMMA)
}

fn element_vdw(z: u8) -> (f64, f64) {
    match z {
        1 => (0.0200, 3.195),
        6 => (0.0951, 3.851),
        7 => (0.0774, 3.660),
        8 => (0.0957, 3.500),
        9 => (0.0725, 3.364),
        16 => (0.2500, 4.035),
        17 => (0.2270, 3.947),
        _ => (0.0500, 3.500),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::reaxff::bond_order::compute_bond_orders_pairs;

    #[test]
    fn bond_orders_match_reference_implementation() {
        let params = ReaxffParams::default_chon();
        let elements = [6u8, 1, 1, 1, 8];
        let positions = [
            0.0, 0.0, 0.0, 1.09, 0.0, 0.0, -0.36, 1.03, 0.0, -0.36, -0.52, 0.89, 1.43, -0.10, 1.02,
        ];
        let atom_params = atom_params_for(&params, &elements);
        let ours = build_bond_orders(&positions, &atom_params, &params, params.cutoff);
        let reference = compute_bond_orders_pairs(&positions, &atom_params, &params, params.cutoff);
        for i in 0..elements.len() {
            for j in 0..elements.len() {
                assert!(
                    (ours.total[i][j] - reference[i][j].total).abs() < 1e-12,
                    "total[{i}][{j}]: {} vs {}",
                    ours.total[i][j],
                    reference[i][j].total
                );
                assert!(
                    (ours.sigma[i][j] - reference[i][j].sigma).abs() < 1e-12,
                    "sigma[{i}][{j}]: {} vs {}",
                    ours.sigma[i][j],
                    reference[i][j].sigma
                );
            }
        }
    }
}
