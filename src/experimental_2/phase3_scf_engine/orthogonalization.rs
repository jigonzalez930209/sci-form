//! Löwdin and canonical orthogonalization of the overlap matrix.
//!
//! The generalized eigenvalue problem Fc = Scε is converted to a
//! standard eigenvalue problem by orthogonalizing the basis:
//!
//!   F' = X†FX,  F'c' = c'ε,  c = Xc'
//!
//! where X = S^{-1/2} is the (symmetric) orthogonalization matrix.
//!
//! # Löwdin Orthogonalization
//!
//! S = UΛU†  →  S^{-1/2} = UΛ^{-1/2}U†
//!
//! This preserves the character of the original basis functions
//! as much as possible (least-distortion orthogonalization).

use nalgebra::DMatrix;

/// Compute the symmetric orthogonalization matrix X = S^{-1/2}.
///
/// Uses eigendecomposition of S:
///   S = U Λ U^T
///   X = U Λ^{-1/2} U^T
///
/// Eigenvalues below `threshold` are discarded (linear dependency).
pub fn lowdin_orthogonalization(s: &DMatrix<f64>, threshold: f64) -> (DMatrix<f64>, usize) {
    let n = s.nrows();
    let eigen = s.clone().symmetric_eigen();

    // Count valid eigenvalues (above threshold)
    let n_independent = eigen.eigenvalues.iter().filter(|&&e| e > threshold).count();

    // Build X = U Λ^{-1/2} U^T
    let mut x = DMatrix::zeros(n, n);

    for k in 0..n {
        let eigenval = eigen.eigenvalues[k];
        if eigenval > threshold {
            let inv_sqrt = 1.0 / eigenval.sqrt();
            // X += inv_sqrt * u_k * u_k^T
            for i in 0..n {
                for j in 0..n {
                    x[(i, j)] += inv_sqrt * eigen.eigenvectors[(i, k)] * eigen.eigenvectors[(j, k)];
                }
            }
        }
    }

    (x, n_independent)
}

/// Transform a matrix to the orthogonal basis: M' = X† M X.
pub fn transform_to_orthogonal(m: &DMatrix<f64>, x: &DMatrix<f64>) -> DMatrix<f64> {
    x.transpose() * m * x
}

/// Back-transform coefficients from orthogonal to original basis: C = X C'.
pub fn back_transform(c_prime: &DMatrix<f64>, x: &DMatrix<f64>) -> DMatrix<f64> {
    x * c_prime
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lowdin_identity() {
        // S = I → X = I
        let s = DMatrix::identity(3, 3);
        let (x, n_indep) = lowdin_orthogonalization(&s, 1e-10);
        assert_eq!(n_indep, 3);

        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((x[(i, j)] - expected).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_lowdin_produces_orthonormal_basis() {
        // Create a simple 2×2 overlap matrix
        let mut s = DMatrix::identity(2, 2);
        s[(0, 1)] = 0.5;
        s[(1, 0)] = 0.5;

        let (x, _) = lowdin_orthogonalization(&s, 1e-10);

        // X†SX should equal I
        let identity_check = x.transpose() * &s * &x;
        for i in 0..2 {
            for j in 0..2 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (identity_check[(i, j)] - expected).abs() < 1e-10,
                    "X†SX[{},{}] = {}, expected {}",
                    i, j, identity_check[(i, j)], expected
                );
            }
        }
    }
}
