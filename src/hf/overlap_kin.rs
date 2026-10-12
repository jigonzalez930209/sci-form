//! Overlap (S) and kinetic energy (T) integral evaluator.
//!
//! Uses 1D Cartesian factorization and analytical derivatives.

use super::basis::BasisSet;
use nalgebra::DMatrix;

pub fn compute_overlap_matrix(basis: &BasisSet) -> DMatrix<f64> {
    // Use the cross-validated `scf` engine (fichas F36/F37).
    let scf_basis = super::basis::to_scf_basis(basis);
    crate::scf::overlap_matrix::build_overlap_matrix(&scf_basis)
}

pub fn compute_kinetic_matrix(basis: &BasisSet) -> DMatrix<f64> {
    let scf_basis = super::basis::to_scf_basis(basis);
    crate::scf::kinetic_matrix::build_kinetic_matrix(&scf_basis)
}

#[cfg(test)]
mod tests {
    use super::super::basis::build_sto3g_basis;
    use super::*;

    #[test]
    fn test_overlap_diagonal_one() {
        let basis = build_sto3g_basis(&[1], &[[0.0, 0.0, 0.0]]);
        let s = compute_overlap_matrix(&basis);
        assert!((s[(0, 0)] - 1.0).abs() < 1e-10, "S_00 = {}", s[(0, 0)]);
    }

    #[test]
    fn test_overlap_symmetric() {
        let basis = build_sto3g_basis(
            &[8, 1, 1],
            &[[0.0, 0.0, 0.0], [0.0, 0.0, 1.8], [0.0, 1.8, 0.0]],
        );
        let s = compute_overlap_matrix(&basis);
        for i in 0..s.nrows() {
            for j in 0..s.nrows() {
                assert!((s[(i, j)] - s[(j, i)]).abs() < 1e-10);
            }
        }
    }
}
