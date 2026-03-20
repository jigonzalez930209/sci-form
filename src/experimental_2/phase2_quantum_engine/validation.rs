//! Validation utilities for comparing CPU and GPU results.
//!
//! Provides tools to verify that GPU-computed matrices match their
//! CPU reference implementations within acceptable tolerances.

use nalgebra::DMatrix;

/// Result of a matrix comparison.
#[derive(Debug, Clone)]
pub struct ComparisonResult {
    /// Maximum absolute difference between elements.
    pub max_abs_error: f64,
    /// Mean absolute error.
    pub mean_abs_error: f64,
    /// Root mean square error.
    pub rms_error: f64,
    /// Whether all elements match within the tolerance.
    pub passed: bool,
    /// Number of elements that exceed the tolerance.
    pub n_failures: usize,
    /// Total number of elements compared.
    pub n_elements: usize,
}

/// Compare two matrices element-wise with configurable tolerance.
pub fn compare_matrices(a: &DMatrix<f64>, b: &DMatrix<f64>, tolerance: f64) -> ComparisonResult {
    assert_eq!(a.nrows(), b.nrows(), "Row count mismatch");
    assert_eq!(a.ncols(), b.ncols(), "Column count mismatch");

    let n = a.nrows() * a.ncols();
    let mut max_err = 0.0f64;
    let mut sum_err = 0.0;
    let mut sum_sq_err = 0.0;
    let mut n_fail = 0;

    for i in 0..a.nrows() {
        for j in 0..a.ncols() {
            let diff = (a[(i, j)] - b[(i, j)]).abs();
            max_err = max_err.max(diff);
            sum_err += diff;
            sum_sq_err += diff * diff;
            if diff > tolerance {
                n_fail += 1;
            }
        }
    }

    ComparisonResult {
        max_abs_error: max_err,
        mean_abs_error: sum_err / n as f64,
        rms_error: (sum_sq_err / n as f64).sqrt(),
        passed: n_fail == 0,
        n_failures: n_fail,
        n_elements: n,
    }
}

/// Check that a matrix is symmetric within some tolerance.
pub fn check_symmetry(m: &DMatrix<f64>, tolerance: f64) -> bool {
    if m.nrows() != m.ncols() {
        return false;
    }
    for i in 0..m.nrows() {
        for j in (i + 1)..m.ncols() {
            if (m[(i, j)] - m[(j, i)]).abs() > tolerance {
                return false;
            }
        }
    }
    true
}

/// Check that a matrix is positive definite (all eigenvalues > 0).
pub fn check_positive_definite(m: &DMatrix<f64>) -> bool {
    let eigen = m.clone().symmetric_eigen();
    eigen.eigenvalues.iter().all(|&e| e > -1e-10)
}

/// Verify overlap matrix properties:
/// 1. Symmetric
/// 2. Positive definite
/// 3. Diagonal elements > 0
pub fn validate_overlap_matrix(s: &DMatrix<f64>) -> Vec<String> {
    let mut issues = Vec::new();

    if !check_symmetry(s, 1e-12) {
        issues.push("Overlap matrix is not symmetric".to_string());
    }

    if !check_positive_definite(s) {
        issues.push("Overlap matrix is not positive definite".to_string());
    }

    for i in 0..s.nrows() {
        if s[(i, i)] <= 0.0 {
            issues.push(format!("S[{},{}] = {} is not positive", i, i, s[(i, i)]));
        }
    }

    issues
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compare_identical_matrices() {
        let a = DMatrix::identity(3, 3);
        let b = DMatrix::identity(3, 3);
        let result = compare_matrices(&a, &b, 1e-10);
        assert!(result.passed);
        assert_eq!(result.n_failures, 0);
        assert!(result.max_abs_error < 1e-15);
    }

    #[test]
    fn test_compare_different_matrices() {
        let a = DMatrix::identity(3, 3);
        let mut b = DMatrix::identity(3, 3);
        b[(0, 0)] = 1.1;
        let result = compare_matrices(&a, &b, 0.05);
        assert!(!result.passed);
        assert_eq!(result.n_failures, 1);
    }

    #[test]
    fn test_symmetry_check() {
        let mut m = DMatrix::zeros(3, 3);
        m[(0, 1)] = 1.0;
        m[(1, 0)] = 1.0;
        assert!(check_symmetry(&m, 1e-10));

        m[(0, 1)] = 1.1;
        assert!(!check_symmetry(&m, 0.05));
    }

    #[test]
    fn test_positive_definite() {
        let m = DMatrix::identity(3, 3);
        assert!(check_positive_definite(&m));
    }
}
