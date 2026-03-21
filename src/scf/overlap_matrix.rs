//! Overlap matrix S_μν construction from contracted Gaussian basis.
//!
//! S_μν = ∫ χ_μ(r) χ_ν(r) d³r

use nalgebra::DMatrix;

use super::basis::{BasisFunction, BasisSet};
use super::gaussian_integrals::overlap_primitive;

/// Build the full overlap matrix S for the given basis set.
pub fn build_overlap_matrix(basis: &BasisSet) -> DMatrix<f64> {
    let n = basis.n_basis;
    let mut s = DMatrix::zeros(n, n);

    for i in 0..n {
        s[(i, i)] = contracted_overlap(&basis.functions[i], &basis.functions[i]);

        for j in (i + 1)..n {
            let s_ij = contracted_overlap(&basis.functions[i], &basis.functions[j]);
            s[(i, j)] = s_ij;
            s[(j, i)] = s_ij;
        }
    }

    s
}

/// Compute overlap between two contracted Gaussian basis functions.
fn contracted_overlap(bf_a: &BasisFunction, bf_b: &BasisFunction) -> f64 {
    let mut s = 0.0;

    for prim_a in &bf_a.primitives {
        let norm_a = BasisFunction::normalization(
            prim_a.alpha, bf_a.angular[0], bf_a.angular[1], bf_a.angular[2],
        );

        for prim_b in &bf_b.primitives {
            let norm_b = BasisFunction::normalization(
                prim_b.alpha, bf_b.angular[0], bf_b.angular[1], bf_b.angular[2],
            );

            let s_prim = overlap_primitive(
                prim_a.alpha, &bf_a.center, bf_a.angular,
                prim_b.alpha, &bf_b.center, bf_b.angular,
            );

            s += norm_a * prim_a.coefficient * norm_b * prim_b.coefficient * s_prim;
        }
    }

    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_overlap_diagonal_positive() {
        let basis = BasisSet::sto3g(&[1, 1], &[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]);
        let s = build_overlap_matrix(&basis);

        for i in 0..basis.n_basis {
            assert!(s[(i, i)] > 0.0, "S[{0},{0}] = {1} should be > 0", i, s[(i, i)]);
        }
    }

    #[test]
    fn test_overlap_symmetric() {
        let basis = BasisSet::sto3g(
            &[8, 1, 1],
            &[[0.0, 0.0, 0.0], [1.43, 1.11, 0.0], [-1.43, 1.11, 0.0]],
        );
        let s = build_overlap_matrix(&basis);

        for i in 0..basis.n_basis {
            for j in 0..basis.n_basis {
                assert!(
                    (s[(i, j)] - s[(j, i)]).abs() < 1e-14,
                    "S not symmetric at ({}, {})", i, j
                );
            }
        }
    }

    #[test]
    fn test_overlap_h2_reasonable() {
        let basis = BasisSet::sto3g(&[1, 1], &[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]);
        let s = build_overlap_matrix(&basis);
        assert!(s[(0, 1)] > 0.0 && s[(0, 1)] < 1.0);
    }
}
