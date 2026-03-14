//! Hamiltonian matrix H construction using EHT (Extended Hückel Theory).
//!
//! Diagonal: H_ii = VSIP of orbital i
//! Off-diagonal: Wolfsberg-Helmholtz approximation:
//!   H_ij = 0.5 * K * S_ij * (H_ii + H_jj)

use nalgebra::DMatrix;

use super::basis::AtomicOrbital;
use super::params::DEFAULT_K;

/// Build the EHT Hamiltonian matrix.
///
/// - `basis`: molecular orbital basis set
/// - `s_matrix`: pre-computed overlap matrix
/// - `k`: Wolfsberg-Helmholtz constant (pass `None` for default 1.75)
pub fn build_hamiltonian(
    basis: &[AtomicOrbital],
    s_matrix: &DMatrix<f64>,
    k: Option<f64>,
) -> DMatrix<f64> {
    let k_val = k.unwrap_or(DEFAULT_K);
    let n = basis.len();
    let mut h = DMatrix::zeros(n, n);

    // Diagonal: H_ii = VSIP
    for i in 0..n {
        h[(i, i)] = basis[i].vsip;
    }

    // Off-diagonal: Wolfsberg-Helmholtz
    for i in 0..n {
        for j in (i + 1)..n {
            let hij = 0.5 * k_val * s_matrix[(i, j)] * (h[(i, i)] + h[(j, j)]);
            h[(i, j)] = hij;
            h[(j, i)] = hij;
        }
    }

    h
}

#[cfg(test)]
mod tests {
    use super::super::basis::build_basis;
    use super::super::overlap::build_overlap_matrix;
    use super::*;

    #[test]
    fn test_hamiltonian_diagonal_vsip() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let basis = build_basis(&elements, &positions);
        let s = build_overlap_matrix(&basis);
        let h = build_hamiltonian(&basis, &s, None);

        // Diagonal should be the VSIP values
        for (i, orb) in basis.iter().enumerate() {
            assert!(
                (h[(i, i)] - orb.vsip).abs() < 1e-12,
                "H[{},{}] = {}, expected VSIP = {}",
                i,
                i,
                h[(i, i)],
                orb.vsip,
            );
        }
    }

    #[test]
    fn test_hamiltonian_symmetry() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let basis = build_basis(&elements, &positions);
        let s = build_overlap_matrix(&basis);
        let h = build_hamiltonian(&basis, &s, None);
        for i in 0..h.nrows() {
            for j in 0..h.ncols() {
                assert!(
                    (h[(i, j)] - h[(j, i)]).abs() < 1e-14,
                    "H not symmetric at ({},{})",
                    i,
                    j,
                );
            }
        }
    }

    #[test]
    fn test_wolfsberg_helmholtz_formula() {
        // Verify off-diagonal is exactly 0.5 * K * S_ij * (H_ii + H_jj)
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let basis = build_basis(&elements, &positions);
        let s = build_overlap_matrix(&basis);
        let h = build_hamiltonian(&basis, &s, None);

        let expected = 0.5 * DEFAULT_K * s[(0, 1)] * (basis[0].vsip + basis[1].vsip);
        assert!(
            (h[(0, 1)] - expected).abs() < 1e-12,
            "H_01 = {}, expected = {}",
            h[(0, 1)],
            expected,
        );
    }

    #[test]
    fn test_custom_k_constant() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let basis = build_basis(&elements, &positions);
        let s = build_overlap_matrix(&basis);
        let h = build_hamiltonian(&basis, &s, Some(2.0));

        let expected = 0.5 * 2.0 * s[(0, 1)] * (basis[0].vsip + basis[1].vsip);
        assert!((h[(0, 1)] - expected).abs() < 1e-12);
    }
}
