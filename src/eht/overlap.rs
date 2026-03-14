//! Overlap matrix S construction using STO-3G Gaussian integrals.
//!
//! S_ij = ∫ φ_i(r) φ_j(r) dr
//! Computed analytically using the Gaussian product theorem.

use nalgebra::DMatrix;
use std::f64::consts::PI;

use super::basis::AtomicOrbital;

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

    for ga in &a.gaussians {
        for gb in &b.gaussians {
            let alpha = ga.alpha;
            let beta = gb.alpha;
            let gamma = alpha + beta;

            // Gaussian product center
            let px = (alpha * a.center[0] + beta * b.center[0]) / gamma;
            let py = (alpha * a.center[1] + beta * b.center[1]) / gamma;
            let pz = (alpha * a.center[2] + beta * b.center[2]) / gamma;

            // Pre-exponential factor from Gaussian product theorem
            let k_ab = (-alpha * beta / gamma * r2_ab).exp();

            let contrib = match (a.l, a.m, b.l, b.m) {
                // (s|s) overlap
                (0, 0, 0, 0) => {
                    let norm_a = (2.0 * alpha / PI).powf(0.75);
                    let norm_b = (2.0 * beta / PI).powf(0.75);
                    let s_prim = (PI / gamma).powf(1.5) * k_ab;
                    norm_a * norm_b * s_prim
                }
                // (s|p) overlap
                (0, 0, 1, m) | (1, m, 0, 0) => {
                    let (s_orb, p_orb, s_alpha, p_alpha) = if a.l == 0 {
                        (a, b, alpha, beta)
                    } else {
                        (b, a, beta, alpha)
                    };
                    let s_a = s_alpha;
                    let p_a = p_alpha;
                    let g = s_a + p_a;

                    let norm_s = (2.0 * s_a / PI).powf(0.75);
                    let norm_p = (128.0 * p_a.powi(5) / (PI * PI * PI)).powf(0.25);

                    // Product center
                    let pc = [
                        (s_a * s_orb.center[0] + p_a * p_orb.center[0]) / g,
                        (s_a * s_orb.center[1] + p_a * p_orb.center[1]) / g,
                        (s_a * s_orb.center[2] + p_a * p_orb.center[2]) / g,
                    ];

                    // Component index from m: -1→x, 0→y, 1→z
                    let comp_idx = match m {
                        -1 => 0usize,
                        0 => 1,
                        1 => 2,
                        _ => 0,
                    };

                    // (P - B_p) where B_p is the p-orbital center
                    let pb = pc[comp_idx] - p_orb.center[comp_idx];

                    let base = (PI / g).powf(1.5) * k_ab;
                    norm_s * norm_p * pb * base
                }
                // (p|p) overlap
                (1, m_a, 1, m_b) => {
                    let norm_a = (128.0 * alpha.powi(5) / (PI * PI * PI)).powf(0.25);
                    let norm_b = (128.0 * beta.powi(5) / (PI * PI * PI)).powf(0.25);

                    let comp_a = match m_a {
                        -1 => 0usize,
                        0 => 1,
                        1 => 2,
                        _ => 0,
                    };
                    let comp_b = match m_b {
                        -1 => 0usize,
                        0 => 1,
                        1 => 2,
                        _ => 0,
                    };

                    // PA_i and PB_j
                    let pa_i = px - a.center[comp_a];
                    let pb_j;
                    // For p-p: if same direction use shift product + 1/(2γ),
                    // else just shift product
                    if comp_a == comp_b {
                        let shift_a = match comp_a {
                            0 => px - a.center[0],
                            1 => py - a.center[1],
                            2 => pz - a.center[2],
                            _ => 0.0,
                        };
                        let shift_b = match comp_b {
                            0 => px - b.center[0],
                            1 => py - b.center[1],
                            2 => pz - b.center[2],
                            _ => 0.0,
                        };
                        let base = (PI / gamma).powf(1.5) * k_ab;
                        let val = shift_a * shift_b + 0.5 / gamma;
                        norm_a * norm_b * val * base
                    } else {
                        pb_j = match comp_b {
                            0 => px - b.center[0],
                            1 => py - b.center[1],
                            2 => pz - b.center[2],
                            _ => 0.0,
                        };
                        let base = (PI / gamma).powf(1.5) * k_ab;
                        norm_a * norm_b * pa_i * pb_j * base
                    }
                }
                _ => 0.0,
            };

            s += ga.coeff * gb.coeff * contrib;
        }
    }

    s
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
}
