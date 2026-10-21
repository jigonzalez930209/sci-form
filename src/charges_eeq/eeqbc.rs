//! EEQ-BC 2025 (bond-capacitor electronegativity equilibration) charge model.
//!
//! Faithful port of multicharge's `eeqbc` model (Froitzheim, Mueller, Hansen
//! and Grimme, *J. Chem. Phys.* **162**, 214109 (2025), DOI 10.1063/5.0268978),
//! including the `new_eeqbc2025_model` wrapper constants.
//!
//! Differences from EEQ2019 ([`super::compute_eeq2019_charges`]):
//!
//! - the Coulomb interaction is screened by a *bond-capacitance matrix* `C`
//!   built from pairwise capacitors `sqrt(cap_i·cap_j)·½(1+erf(arg))` with
//!   `arg = -kbc·(r−r_vdw)/r_vdw` and pair-specific van-der-Waals distances;
//! - the right-hand side is `C·x` (not just `x`);
//! - the CN entering the electronegativity is the *directed*,
//!   electronegativity-weighted CN (local charges `qloc`), which also scales
//!   the hardness through `kqeta`.
//!
//! ```text
//!   qloc_i = CN^en_i + Q/N
//!   C_ii   = Σ_j c_ij,  C_ij = −c_ij,  C_NN = 1
//!   x_i    = C · [−chi + kcnchi·CN + kqchi·qloc, Q]
//!   A_ij   = erf(r·γ)/r · C_ij,  γ² = 1/(rad_i² + rad_j²)
//!   A_ii   = (eta + kqeta·qloc + sqrt(2/π)/rad_i)·C_ii + 1
//!   A q = x   (+ Lagrange row/column of ones, A_NN = 0)
//! ```
//!
//! Atomic units throughout (positions in Bohr, energies in Hartree).
//! Parameters from multicharge `eeqbc2025` (Apache-2.0) and mctc-lib
//! (Apache-2.0), see `scripts/validation/gen_eeqbc_params.py`.
//! Validated against the multicharge CLI (`--model eeqbc`, gate F47).

use nalgebra::{DMatrix, DVector};

use super::eeqbc_data::{
    EEQBC_AVG_CN, EEQBC_CAP, EEQBC_CHI, EEQBC_COVRAD, EEQBC_ETA, EEQBC_KCNCHI, EEQBC_KQCHI,
    EEQBC_KQETA, EEQBC_PAULING_EN, EEQBC_RAD, EEQBC_VDW_PAIRS,
};

/// Steepness of the erf CN counting function (`new_eeqbc2025_model`).
const EEQBC_CN_KCN: f64 = 2.0;
/// Exponent of the distance normalization in the CN counting function.
const EEQBC_CN_NORM_EXP: f64 = 0.75;
/// Real-space CN cutoff in Bohr.
const EEQBC_CN_CUTOFF: f64 = 25.0;
/// CN scaling of the Gaussian charge width.
const EEQBC_KCNRAD: f64 = 0.14;
/// Steepness of the bond-capacitor switch.
const EEQBC_KBC: f64 = 0.60;
/// Pauling electronegativity of fluorine (normalization reference).
const PAULING_F: f64 = 3.98;
/// Number of elements covered by the EEQ-BC parameter set.
pub const EEQBC_MAX_Z: usize = 103;

/// `erf(x)` to near machine precision via the regularized incomplete gamma
/// function `erf(x) = sign(x)·P(1/2, x²)` (series / Lentz continued fraction).
///
/// The EEQ2019 module uses the cheaper Abramowitz & Stegun 7.1.26 rational
/// approximation (~1.5e-7); EEQ-BC targets 1e-9 charge agreement with the
/// reference, so the exact path is used here.
fn erf_exact(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let t = x * x;
    let p = if t < 1.5 {
        gamma_p_series(0.5, t)
    } else {
        1.0 - gamma_q_cf(0.5, t)
    };
    sign * p
}

/// ln Γ(1/2) = ln(√π).
const LN_GAMMA_HALF: f64 = 0.572_364_942_924_700_1;

/// Series representation of the regularized lower incomplete gamma `P(a, x)`
/// (converges quickly for `x < a + 1`).
fn gamma_p_series(a: f64, x: f64) -> f64 {
    let mut ap = a;
    let mut sum = 1.0 / a;
    let mut del = sum;
    for _ in 0..500 {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if del.abs() < sum.abs() * 1e-16 {
            break;
        }
    }
    sum * (-x + a * x.ln() - LN_GAMMA_HALF).exp()
}

/// Continued fraction representation of the regularized upper incomplete
/// gamma `Q(a, x)` with the modified Lentz algorithm (for `x >= a + 1`).
fn gamma_q_cf(a: f64, x: f64) -> f64 {
    const TINY: f64 = 1e-300;
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / TINY;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..500 {
        let an = -(i as f64) * (i as f64 - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < TINY {
            d = TINY;
        }
        c = b + an / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-16 {
            break;
        }
    }
    h * (-x + a * x.ln() - LN_GAMMA_HALF).exp()
}

/// erf counting function: `½(1 + erf(−kcn·(r−rc)/rc^norm_exp))`.
fn cn_count(r: f64, rc: f64) -> f64 {
    0.5 * (1.0 + erf_exact(-EEQBC_CN_KCN * (r - rc) / rc.powf(EEQBC_CN_NORM_EXP)))
}

/// Electronegativity weight of the directed CN: Pauling EN (normalized to
/// fluorine) with the actinide overrides applied by `new_eeqbc2025_model`.
fn eeqbc_en(z: u8) -> f64 {
    let z = z as usize;
    let mut en = EEQBC_PAULING_EN[z - 1];
    if z >= 90 {
        en = 1.30;
    }
    if z == 87 {
        en = 0.80;
    }
    if z == 89 {
        en = 1.00;
    }
    if z == 90 || z == 91 || z == 92 || z == 95 {
        en = 1.10;
    }
    if z == 93 || z == 94 || z == 97 || z == 103 {
        en = 1.20;
    }
    en / PAULING_F
}

/// Pairwise van-der-Waals radius (Bohr) from the packed mctc-lib table.
fn eeqbc_rvdw(zi: u8, zj: u8) -> f64 {
    let (lo, hi) = if zi <= zj {
        (zi as usize, zj as usize)
    } else {
        (zj as usize, zi as usize)
    };
    EEQBC_VDW_PAIRS[lo - 1 + hi * (hi - 1) / 2]
}

/// Bond capacitance between two atoms: `sqrt(capi·capj)·½(1+erf(arg))`.
fn bond_capacitance(r1: f64, rvdw: f64, capi: f64, capj: f64) -> f64 {
    let arg = -EEQBC_KBC * (r1 - rvdw) / rvdw;
    (capi * capj).sqrt() * 0.5 * (1.0 + erf_exact(arg))
}

fn pair_distance(pos: &[[f64; 3]], i: usize, j: usize) -> f64 {
    let dx = pos[i][0] - pos[j][0];
    let dy = pos[i][1] - pos[j][1];
    let dz = pos[i][2] - pos[j][2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Plain (undirected) erf coordination numbers.
pub fn eeqbc_coordination_numbers(elements: &[u8], positions_bohr: &[[f64; 3]]) -> Vec<f64> {
    let n = elements.len();
    let mut cn = vec![0.0; n];
    for i in 0..n {
        let zi = elements[i] as usize;
        if zi == 0 || zi > EEQBC_MAX_Z {
            continue;
        }
        for j in 0..i {
            let zj = elements[j] as usize;
            if zj == 0 || zj > EEQBC_MAX_Z {
                continue;
            }
            let r = pair_distance(positions_bohr, i, j);
            if !(1e-6..=EEQBC_CN_CUTOFF).contains(&r) {
                continue;
            }
            let rc = EEQBC_COVRAD[zi - 1] + EEQBC_COVRAD[zj - 1];
            let c = cn_count(r, rc);
            cn[i] += c;
            cn[j] += c;
        }
    }
    cn
}

/// Directed, electronegativity-weighted CN used for the local charges
/// (`qloc`): contributions are added to the more electronegative partner and
/// removed from the less electronegative one.
pub fn eeqbc_local_charges(
    elements: &[u8],
    positions_bohr: &[[f64; 3]],
    total_charge: f64,
) -> Vec<f64> {
    let n = elements.len();
    let mut qloc = vec![0.0; n];
    for i in 0..n {
        let zi = elements[i] as usize;
        if zi == 0 || zi > EEQBC_MAX_Z {
            continue;
        }
        for j in 0..i {
            let zj = elements[j] as usize;
            if zj == 0 || zj > EEQBC_MAX_Z {
                continue;
            }
            let r = pair_distance(positions_bohr, i, j);
            if !(1e-6..=EEQBC_CN_CUTOFF).contains(&r) {
                continue;
            }
            let rc = EEQBC_COVRAD[zi - 1] + EEQBC_COVRAD[zj - 1];
            let den = eeqbc_en(elements[j]) - eeqbc_en(elements[i]);
            let c = den * cn_count(r, rc);
            qloc[i] += c;
            qloc[j] -= c;
        }
    }
    let share = total_charge / n as f64;
    for q in qloc.iter_mut() {
        *q += share;
    }
    qloc
}

/// EEQ-BC 2025 partial charges.
///
/// `positions_bohr`: Cartesian coordinates in Bohr. Returns `None` when an
/// element outside the parameter table (Z = 1..103) is present.
pub fn compute_eeqbc_charges(
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
    let mut kqchi = Vec::with_capacity(n);
    let mut kqeta = Vec::with_capacity(n);
    let mut cap = Vec::with_capacity(n);
    let mut avg_cn = Vec::with_capacity(n);
    for &z in elements {
        let zi = z as usize;
        if zi == 0 || zi > EEQBC_MAX_Z {
            return None;
        }
        chi.push(EEQBC_CHI[zi - 1]);
        eta.push(EEQBC_ETA[zi - 1]);
        rad.push(EEQBC_RAD[zi - 1]);
        kcnchi.push(EEQBC_KCNCHI[zi - 1]);
        kqchi.push(EEQBC_KQCHI[zi - 1]);
        kqeta.push(EEQBC_KQETA[zi - 1]);
        cap.push(EEQBC_CAP[zi - 1]);
        avg_cn.push(EEQBC_AVG_CN[zi - 1]);
    }

    let cn = eeqbc_coordination_numbers(elements, positions_bohr);
    let qloc = eeqbc_local_charges(elements, positions_bohr, total_charge);

    // Bond-capacitance matrix (n + 1, constraint row/column included).
    let mut cmat = DMatrix::zeros(n + 1, n + 1);
    for i in 0..n {
        let capi = cap[i];
        for j in 0..i {
            let r1 = pair_distance(positions_bohr, i, j);
            let rvdw = eeqbc_rvdw(elements[i], elements[j]);
            let c = bond_capacitance(r1, rvdw, capi, cap[j]);
            cmat[(j, i)] -= c;
            cmat[(i, j)] -= c;
            cmat[(i, i)] += c;
            cmat[(j, j)] += c;
        }
    }
    cmat[(n, n)] = 1.0;

    // Right-hand side: xvec = C · xtmp.
    let mut xtmp = DVector::zeros(n + 1);
    for i in 0..n {
        xtmp[i] = -chi[i] + kcnchi[i] * cn[i] + kqchi[i] * qloc[i];
    }
    xtmp[n] = total_charge;
    let xvec = &cmat * &xtmp;

    // Effective (CN-compressed) charge widths.
    let rad_eff: Vec<f64> = (0..n)
        .map(|i| rad[i] * (1.0 - EEQBC_KCNRAD * cn[i] / avg_cn[i].powf(EEQBC_CN_NORM_EXP)))
        .collect();

    // Coulomb matrix with the capacitance screening.
    let mut amat = DMatrix::zeros(n + 1, n + 1);
    for i in 0..n {
        for j in 0..i {
            let r = pair_distance(positions_bohr, i, j);
            let r2 = r * r;
            let gam2 = 1.0 / (rad_eff[i] * rad_eff[i] + rad_eff[j] * rad_eff[j]);
            let tmp = erf_exact((r2 * gam2).sqrt()) / r2.sqrt() * cmat[(j, i)];
            amat[(j, i)] = tmp;
            amat[(i, j)] = tmp;
        }
        let tmp = eta[i] + kqeta[i] * qloc[i] + (2.0 / std::f64::consts::PI).sqrt() / rad_eff[i];
        amat[(i, i)] += tmp * cmat[(i, i)] + 1.0;
    }
    for i in 0..=n {
        amat[(n, i)] = 1.0;
        amat[(i, n)] = 1.0;
    }
    amat[(n, n)] = 0.0;

    amat.lu()
        .solve(&xvec)
        .map(|sol| sol.rows(0, n).iter().copied().collect())
}
