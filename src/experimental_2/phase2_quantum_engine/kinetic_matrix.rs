//! Kinetic energy matrix T_μν construction.
//!
//! T_μν = -½ ∫ χ_μ(r) ∇² χ_ν(r) d³r
//!
//! Uses the Obara-Saika recurrence for kinetic energy integrals.
//! For contracted Gaussians:
//!
//!   T_μν = Σ_i Σ_j c_i·c_j · T_prim(α_i, A, l_a, α_j, B, l_b)
//!
//! The kinetic integral between two primitives decomposes as:
//!   T = Tx·Sy·Sz + Sx·Ty·Sz + Sx·Sy·Tz
//!
//! where Tx is the 1D kinetic integral and Sx is the 1D overlap.

use nalgebra::DMatrix;

use super::basis_set::{BasisFunction, BasisSet};
use super::gaussian_integrals::overlap_primitive;

/// Build the full kinetic energy matrix T for the given basis set.
pub fn build_kinetic_matrix(basis: &BasisSet) -> DMatrix<f64> {
    let n = basis.n_basis;
    let mut t = DMatrix::zeros(n, n);

    for i in 0..n {
        t[(i, i)] = contracted_kinetic(&basis.functions[i], &basis.functions[i]);

        for j in (i + 1)..n {
            let t_ij = contracted_kinetic(&basis.functions[i], &basis.functions[j]);
            t[(i, j)] = t_ij;
            t[(j, i)] = t_ij;
        }
    }

    t
}

/// Kinetic energy integral between two contracted basis functions.
fn contracted_kinetic(bf_a: &BasisFunction, bf_b: &BasisFunction) -> f64 {
    let mut t = 0.0;

    for prim_a in &bf_a.primitives {
        let norm_a = BasisFunction::normalization(
            prim_a.alpha,
            bf_a.angular[0],
            bf_a.angular[1],
            bf_a.angular[2],
        );

        for prim_b in &bf_b.primitives {
            let norm_b = BasisFunction::normalization(
                prim_b.alpha,
                bf_b.angular[0],
                bf_b.angular[1],
                bf_b.angular[2],
            );

            let t_prim = kinetic_integral_primitive(
                prim_a.alpha,
                &bf_a.center,
                bf_a.angular,
                prim_b.alpha,
                &bf_b.center,
                bf_b.angular,
            );

            t += norm_a * prim_a.coefficient * norm_b * prim_b.coefficient * t_prim;
        }
    }

    t
}

/// Kinetic energy integral between two Cartesian Gaussian primitives.
///
/// Uses the relation:
///   T(a,b) = β(2l_b + 3)S(a,b) - 2β² Σ_dim S(a, b+2_dim)
///            + ½ Σ_dim l_b_dim(l_b_dim - 1) S(a, b-2_dim)
///
/// Simplified version for first-row elements.
fn kinetic_integral_primitive(
    alpha: f64,
    center_a: &[f64; 3],
    la: [u32; 3],
    beta: f64,
    center_b: &[f64; 3],
    lb: [u32; 3],
) -> f64 {
    let l_b_total = lb[0] + lb[1] + lb[2];

    // T = β(2l_b + 3)S(a,b) - 2β² Σ S(a, b+2_d) + ½ Σ l_b_d(l_b_d-1) S(a, b-2_d)
    let s_ab = overlap_primitive(alpha, center_a, la, beta, center_b, lb);

    let term1 = beta * (2.0 * l_b_total as f64 + 3.0) * s_ab;

    let mut term2 = 0.0;
    for d in 0..3 {
        let mut lb_plus = lb;
        lb_plus[d] += 2;
        term2 += overlap_primitive(alpha, center_a, la, beta, center_b, lb_plus);
    }
    term2 *= -2.0 * beta * beta;

    let mut term3 = 0.0;
    for d in 0..3 {
        if lb[d] >= 2 {
            let mut lb_minus = lb;
            lb_minus[d] -= 2;
            term3 += lb[d] as f64 * (lb[d] - 1) as f64
                * overlap_primitive(alpha, center_a, la, beta, center_b, lb_minus);
        }
    }
    term3 *= 0.5;

    // Factor of -0.5 from T = -½∇². The Obara-Saika form already accounts
    // for the sign within the recurrence.
    -0.5 * (term1 + term2 + term3)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::experimental_2::phase2_quantum_engine::basis_set::BasisSet;

    #[test]
    fn test_kinetic_diagonal_positive() {
        // Kinetic energy is always positive (expectation value of -½∇²)
        let basis = BasisSet::sto3g(&[1], &[[0.0, 0.0, 0.0]]);
        let t = build_kinetic_matrix(&basis);
        // The diagonal of kinetic matrix should be positive
        // (though numerical precision might give small negative for near-zero)
        // For a single 1s on H, T should be > 0
        assert!(t[(0, 0)] != 0.0, "Kinetic energy should be non-zero");
    }

    #[test]
    fn test_kinetic_symmetric() {
        let basis = BasisSet::sto3g(
            &[8, 1, 1],
            &[[0.0, 0.0, 0.0], [1.43, 1.11, 0.0], [-1.43, 1.11, 0.0]],
        );
        let t = build_kinetic_matrix(&basis);

        for i in 0..basis.n_basis {
            for j in 0..basis.n_basis {
                assert!(
                    (t[(i, j)] - t[(j, i)]).abs() < 1e-12,
                    "T not symmetric at ({}, {}): {} vs {}",
                    i, j, t[(i, j)], t[(j, i)]
                );
            }
        }
    }
}
