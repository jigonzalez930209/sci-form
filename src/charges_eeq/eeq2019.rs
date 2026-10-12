//! EEQ2019 charge model (Grimme group, multicharge reference port).
//!
//! Electronegativity equilibration with coordination-number-dependent
//! electronegativities, erf-damped Coulomb interactions and Gaussian charge
//! widths, solved with a total-charge constraint:
//!
//! ```text
//!   A_ii = eta_i + sqrt(2/pi) / rad_i
//!   A_ij = erf(sqrt(r^2 / (rad_i^2 + rad_j^2))) / r
//!   x_i  = -chi_i + kcnchi_i * CN_i / sqrt(CN_i)
//!   A q = x  (+ Lagrange row/column of ones, x[n] = total charge)
//! ```
//!
//! Atomic units throughout (positions in Bohr, energies in Hartree).
//! Reference parameters from multicharge `eeq2019` (LGPL-3.0-or-later).
//! Validated against dftd4 `get_properties` partial charges (F46).

use nalgebra::{DMatrix, DVector};

use super::eeq2019_data::{EEQ2019_CHI, EEQ2019_ETA, EEQ2019_KCNCHI, EEQ2019_RAD};

/// Steepness of the erf CN counting function (multicharge default).
const EEQ_CN_KCN: f64 = 7.5;
/// Smooth CN cutoff (multicharge default; identity-like below ~6).
const EEQ_CN_MAX: f64 = 8.0;
/// Real-space CN cutoff in Bohr.
const EEQ_CN_CUTOFF: f64 = 25.0;
/// Bohr per Angstrom.
const BOHR_PER_ANGSTROM: f64 = 1.8897259886;

fn log_cn_cut(cn: f64, cn_max: f64) -> f64 {
    (1.0 + cn_max.exp()).ln() - (1.0 + (cn_max - cn).exp()).ln()
}

fn erf_approx(x: f64) -> f64 {
    // Abramowitz & Stegun 7.1.26 (same accuracy class as the reference).
    let (sign, x) = if x >= 0.0 { (1.0, x) } else { (-1.0, -x) };
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = 0.254829592 * t - 0.284496736 * t * t + 1.421413741 * t * t * t
        - 1.453152027 * t * t * t * t
        + 1.061405429 * t * t * t * t * t;
    sign * (1.0 - poly * (-x * x).exp())
}

/// Fractional coordination numbers (plain erf, D3 covalent radii, log-cut).
///
/// Note: like the reference implementation, the radii are the D3 set, i.e.
/// 4/3 x the 2009 values (mctc `get_covalent_rad` returns the D3 set).
fn eeq2019_cn(elements: &[u8], positions_bohr: &[[f64; 3]]) -> Vec<f64> {
    use super::eeq2019_data::EEQ2019_COVRAD_ANG;
    let n = elements.len();
    let mut cn = vec![0.0; n];
    for i in 0..n {
        let zi = elements[i] as usize;
        if zi == 0 || zi > EEQ2019_COVRAD_ANG.len() {
            continue;
        }
        for j in 0..n {
            if i == j {
                continue;
            }
            let zj = elements[j] as usize;
            if zj == 0 || zj > EEQ2019_COVRAD_ANG.len() {
                continue;
            }
            let dx = positions_bohr[i][0] - positions_bohr[j][0];
            let dy = positions_bohr[i][1] - positions_bohr[j][1];
            let dz = positions_bohr[i][2] - positions_bohr[j][2];
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 > EEQ_CN_CUTOFF * EEQ_CN_CUTOFF || r2 < 1e-12 {
                continue;
            }
            let r = r2.sqrt();
            let rc = (EEQ2019_COVRAD_ANG[zi - 1] + EEQ2019_COVRAD_ANG[zj - 1])
                * BOHR_PER_ANGSTROM
                * (4.0 / 3.0);
            cn[i] += 0.5 * (1.0 + erf_approx(-EEQ_CN_KCN * (r - rc) / rc));
        }
    }
    cn.iter().map(|&c| log_cn_cut(c, EEQ_CN_MAX)).collect()
}

fn eeq_param(table: &[f64; 86], z: u8) -> Option<f64> {
    let zi = z as usize;
    if zi == 0 || zi > table.len() {
        return None;
    }
    Some(table[zi - 1])
}

/// EEQ2019 partial charges.
///
/// `positions_bohr`: Cartesian coordinates in Bohr. Returns `None` when an
/// element outside the parameter table is present.
pub fn compute_eeq2019_charges(
    elements: &[u8],
    positions_bohr: &[[f64; 3]],
    total_charge: f64,
) -> Option<Vec<f64>> {
    let n = elements.len();
    if n == 0 {
        return Some(vec![]);
    }
    let mut chi = Vec::with_capacity(n);
    let mut eta = Vec::with_capacity(n);
    let mut rad = Vec::with_capacity(n);
    let mut kcnchi = Vec::with_capacity(n);
    for &z in elements {
        chi.push(eeq_param(&EEQ2019_CHI, z)?);
        eta.push(eeq_param(&EEQ2019_ETA, z)?);
        rad.push(eeq_param(&EEQ2019_RAD, z)?);
        kcnchi.push(eeq_param(&EEQ2019_KCNCHI, z)?);
    }

    let cn = eeq2019_cn(elements, positions_bohr);

    // Coulomb matrix with Gaussian charge widths.
    let mut amat = DMatrix::zeros(n + 1, n + 1);
    for i in 0..n {
        for j in 0..i {
            let r2 = (0..3)
                .map(|d| {
                    let dd = positions_bohr[j][d] - positions_bohr[i][d];
                    dd * dd
                })
                .sum::<f64>();
            let gam = 1.0 / (rad[i] * rad[i] + rad[j] * rad[j]);
            let tmp = erf_approx((r2 * gam).sqrt()) / r2.sqrt();
            amat[(j, i)] += tmp;
            amat[(i, j)] += tmp;
        }
        amat[(i, i)] += eta[i] + (2.0 / std::f64::consts::PI).sqrt() / rad[i];
    }
    for i in 0..=n {
        amat[(n, i)] = 1.0;
        amat[(i, n)] = 1.0;
    }
    amat[(n, n)] = 0.0;

    // Right-hand side: CN-dependent electronegativities.
    let mut xvec = DVector::zeros(n + 1);
    for (i, x) in xvec.iter_mut().enumerate().take(n) {
        let c = cn[i].max(1e-14);
        *x = -chi[i] + kcnchi[i] * cn[i] / c.sqrt();
    }
    xvec[n] = total_charge;

    amat.lu()
        .solve(&xvec)
        .map(|sol| sol.rows(0, n).iter().copied().collect())
}
