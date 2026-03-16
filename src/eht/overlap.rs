//! Overlap matrix S construction using STO-3G Gaussian integrals.
//!
//! S_ij = ∫ φ_i(r) φ_j(r) dr
//! Computed analytically using the Gaussian product theorem.

use nalgebra::DMatrix;
use std::f64::consts::PI;

use super::basis::{orbital_cartesian_terms, AtomicOrbital};

/// Build the symmetric overlap matrix S for the given molecular basis.
pub fn build_overlap_matrix(basis: &[AtomicOrbital]) -> DMatrix<f64> {
    let n = basis.len();
    let mut s = DMatrix::zeros(n, n);

    for i in 0..n {
        s[(i, i)] = 1.0; // S_ii = 1 by definition for normalized orbitals
        for j in (i + 1)..n {
            let sij = overlap_integral(&basis[i], &basis[j]);
            s[(i, j)] = sij;
            s[(j, i)] = sij;
        }
    }

    s
}

/// Compute the overlap integral between two atomic orbitals using
/// their STO-3G contracted Gaussian expansions.
fn overlap_integral(a: &AtomicOrbital, b: &AtomicOrbital) -> f64 {
    let dx = a.center[0] - b.center[0];
    let dy = a.center[1] - b.center[1];
    let dz = a.center[2] - b.center[2];
    let r2_ab = dx * dx + dy * dy + dz * dz;

    let mut s = 0.0;

    let terms_a = orbital_cartesian_terms(a.l, a.m);
    let terms_b = orbital_cartesian_terms(b.l, b.m);

    for ga in &a.gaussians {
        for gb in &b.gaussians {
            let alpha = ga.alpha;
            let beta = gb.alpha;
            let gamma = alpha + beta;

            let px = (alpha * a.center[0] + beta * b.center[0]) / gamma;
            let py = (alpha * a.center[1] + beta * b.center[1]) / gamma;
            let pz = (alpha * a.center[2] + beta * b.center[2]) / gamma;

            let pa = [px - a.center[0], py - a.center[1], pz - a.center[2]];
            let pb = [px - b.center[0], py - b.center[1], pz - b.center[2]];
            let pre = (-alpha * beta / gamma * r2_ab).exp();

            let mut prim_sum = 0.0;
            for (coef_a, [lax, lay, laz]) in &terms_a {
                for (coef_b, [lbx, lby, lbz]) in &terms_b {
                    let ix = overlap_1d(*lax as usize, *lbx as usize, pa[0], pb[0], gamma);
                    let iy = overlap_1d(*lay as usize, *lby as usize, pa[1], pb[1], gamma);
                    let iz = overlap_1d(*laz as usize, *lbz as usize, pa[2], pb[2], gamma);

                    let norm_a = cart_norm(alpha, *lax as i32, *lay as i32, *laz as i32);
                    let norm_b = cart_norm(beta, *lbx as i32, *lby as i32, *lbz as i32);
                    prim_sum += coef_a * coef_b * norm_a * norm_b * ix * iy * iz;
                }
            }

            s += ga.coeff * gb.coeff * pre * prim_sum;
        }
    }

    s
}

fn odd_double_factorial(n: i32) -> f64 {
    if n <= 0 {
        return 1.0;
    }
    let mut acc = 1.0;
    let mut k = n;
    while k > 0 {
        acc *= k as f64;
        k -= 2;
    }
    acc
}

fn cart_norm(alpha: f64, lx: i32, ly: i32, lz: i32) -> f64 {
    let lsum = lx + ly + lz;
    let pref = (2.0 * alpha / PI).powf(0.75);
    let ang = (4.0 * alpha).powf(lsum as f64 / 2.0);
    let denom = (odd_double_factorial(2 * lx - 1)
        * odd_double_factorial(2 * ly - 1)
        * odd_double_factorial(2 * lz - 1))
    .sqrt();
    pref * ang / denom
}

fn overlap_1d(la: usize, lb: usize, pa: f64, pb: f64, gamma: f64) -> f64 {
    let mut e = vec![vec![0.0f64; lb + 1]; la + 1];
    e[0][0] = (PI / gamma).sqrt();

    for i in 1..=la {
        let term1 = pa * e[i - 1][0];
        let term2 = if i > 1 {
            (i as f64 - 1.0) * e[i - 2][0] / (2.0 * gamma)
        } else {
            0.0
        };
        e[i][0] = term1 + term2;
    }

    for j in 1..=lb {
        let term1 = pb * e[0][j - 1];
        let term2 = if j > 1 {
            (j as f64 - 1.0) * e[0][j - 2] / (2.0 * gamma)
        } else {
            0.0
        };
        e[0][j] = term1 + term2;
    }

    for i in 1..=la {
        for j in 1..=lb {
            let t1 = pb * e[i][j - 1];
            let t2 = if j > 1 {
                (j as f64 - 1.0) * e[i][j - 2] / (2.0 * gamma)
            } else {
                0.0
            };
            let t3 = (i as f64) * e[i - 1][j - 1] / (2.0 * gamma);
            e[i][j] = t1 + t2 + t3;
        }
    }

    e[la][lb]
}

#[cfg(test)]
mod tests {
    use super::super::basis::build_basis;
    use super::*;

    #[test]
    fn test_overlap_diagonal_is_one() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let basis = build_basis(&elements, &positions);
        let s = build_overlap_matrix(&basis);
        for i in 0..s.nrows() {
            assert!(
                (s[(i, i)] - 1.0).abs() < 1e-10,
                "S[{},{}] = {}",
                i,
                i,
                s[(i, i)]
            );
        }
    }

    #[test]
    fn test_overlap_symmetry() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let basis = build_basis(&elements, &positions);
        let s = build_overlap_matrix(&basis);
        for i in 0..s.nrows() {
            for j in 0..s.ncols() {
                assert!(
                    (s[(i, j)] - s[(j, i)]).abs() < 1e-14,
                    "S not symmetric at ({},{})",
                    i,
                    j,
                );
            }
        }
    }

    #[test]
    fn test_overlap_off_diagonal_range() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let basis = build_basis(&elements, &positions);
        let s = build_overlap_matrix(&basis);
        // All off-diagonal elements should have |S_ij| < 1
        for i in 0..s.nrows() {
            for j in 0..s.ncols() {
                if i != j {
                    assert!(
                        s[(i, j)].abs() <= 1.0 + 1e-10,
                        "S[{},{}] = {} exceeds 1.0",
                        i,
                        j,
                        s[(i, j)]
                    );
                }
            }
        }
    }

    #[test]
    fn test_overlap_h2_ss() {
        // Two H atoms at 0.74 Å: their 1s-1s overlap should be ~0.6-0.8
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let basis = build_basis(&elements, &positions);
        let s = build_overlap_matrix(&basis);
        assert_eq!(s.nrows(), 2);
        let s12 = s[(0, 1)];
        assert!(s12 > 0.3 && s12 < 0.9, "S_12 for H2 = {}", s12);
    }

    #[test]
    fn test_metal_overlap_matrix_is_finite() {
        let elements = [26u8, 17, 17];
        let positions = [[0.0, 0.0, 0.0], [2.2, 0.0, 0.0], [-2.2, 0.0, 0.0]];
        let basis = build_basis(&elements, &positions);
        let s = build_overlap_matrix(&basis);

        for i in 0..s.nrows() {
            for j in 0..s.ncols() {
                assert!(s[(i, j)].is_finite(), "S[{},{}] is not finite", i, j);
            }
        }
    }
}
