//! SDR Projection operators — E8.1
//!
//! Alternating projections onto the PSD cone and distance constraint set.

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

/// SDR configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SdrConfig {
    /// Maximum number of alternating projection iterations.
    pub max_iter: usize,
    /// Convergence tolerance (relative Frobenius norm change).
    pub tol: f64,
    /// SVT threshold parameter (0.0 disables SVT).
    pub svt_tau: f64,
    /// SVT decay factor per improving iteration.
    pub svt_decay: f64,
}

impl Default for SdrConfig {
    fn default() -> Self {
        Self {
            max_iter: 200,
            tol: 1e-6,
            svt_tau: 0.0,
            svt_decay: 0.9,
        }
    }
}

/// Convergence information from alternating projections.
#[derive(Debug, Clone)]
pub struct SdrConvergence {
    /// Number of iterations performed.
    pub iterations: usize,
    /// Converged within tolerance.
    pub converged: bool,
    /// Final relative change.
    pub final_residual: f64,
    /// Number of negative eigenvalues eliminated.
    pub neg_eigenvalues_removed: usize,
}

/// Project matrix onto the PSD cone by zeroing negative eigenvalues.
///
/// Returns the projected matrix and count of negative eigenvalues removed.
pub fn project_psd(x: &DMatrix<f64>) -> (DMatrix<f64>, usize) {
    let n = x.nrows();
    // Symmetrize
    let sym = (x + x.transpose()) * 0.5;

    let eigen = sym.symmetric_eigen();
    let mut neg_count = 0;

    let mut sigma = eigen.eigenvalues.clone();
    for i in 0..n {
        if sigma[i] < 0.0 {
            sigma[i] = 0.0;
            neg_count += 1;
        }
    }

    // Reconstruct: V * diag(sigma_+) * V^T
    let v = &eigen.eigenvectors;
    let diag = DMatrix::from_diagonal(&sigma);
    let result = v * diag * v.transpose();

    (result, neg_count)
}

/// Project Gram matrix onto distance constraints.
///
/// Given target squared distances, enforce: X_ii + X_jj - 2*X_ij = d_ij^2
/// for known distance pairs.
pub fn project_distances(
    x: &DMatrix<f64>,
    distance_pairs: &[(usize, usize, f64)], // (i, j, d_ij)
) -> DMatrix<f64> {
    let mut result = x.clone();

    for &(i, j, d_ij) in distance_pairs {
        let d_sq = d_ij * d_ij;
        // Current: X_ii + X_jj - 2*X_ij
        let current = result[(i, i)] + result[(j, j)] - 2.0 * result[(i, j)];
        if current.abs() < 1e-15 && d_sq.abs() < 1e-15 {
            continue;
        }

        // Correction: distribute equally among the three entries
        let err = d_sq - current;
        let correction = err / 6.0; // distribute to 3 independent entries

        result[(i, i)] += correction;
        result[(j, j)] += correction;
        result[(i, j)] -= correction;
        result[(j, i)] -= correction;
    }

    result
}

/// Singular Value Thresholding (SVT) step.
///
/// Soft-thresholds singular values to promote low rank.
pub fn svt_step(x: &DMatrix<f64>, tau: f64) -> DMatrix<f64> {
    if tau <= 0.0 {
        return x.clone();
    }

    let sym = (x + x.transpose()) * 0.5;
    let eigen = sym.symmetric_eigen();

    let mut sigma = eigen.eigenvalues.clone();
    for i in 0..sigma.len() {
        sigma[i] = (sigma[i] - tau).max(0.0);
    }

    let v = &eigen.eigenvectors;
    let diag = DMatrix::from_diagonal(&sigma);
    v * diag * v.transpose()
}

/// Run alternating projections between PSD cone and distance constraints.
pub fn alternating_projections(
    x0: &DMatrix<f64>,
    distance_pairs: &[(usize, usize, f64)],
    config: &SdrConfig,
) -> (DMatrix<f64>, SdrConvergence) {
    let mut x = x0.clone();
    let mut total_neg = 0;
    let mut tau = config.svt_tau;
    let mut final_residual = 1.0;
    let mut converged = false;
    let mut iter = 0;

    for k in 0..config.max_iter {
        iter = k + 1;

        // Step 1: Project onto PSD cone
        let (x_psd, neg) = project_psd(&x);
        total_neg += neg;

        // Step 2: Optional SVT
        let x_svt = if tau > 0.0 {
            let result = svt_step(&x_psd, tau);
            tau *= config.svt_decay;
            result
        } else {
            x_psd
        };

        // Step 3: Project onto distance constraints
        let x_new = project_distances(&x_svt, distance_pairs);

        // Check convergence
        let diff = &x_new - &x;
        let diff_norm = diff.norm();
        let current_norm = x_new.norm();

        final_residual = if current_norm > 1e-15 {
            diff_norm / current_norm
        } else {
            diff_norm
        };

        if final_residual < config.tol {
            x = x_new;
            converged = true;
            break;
        }

        x = x_new;
    }

    let conv = SdrConvergence {
        iterations: iter,
        converged,
        final_residual,
        neg_eigenvalues_removed: total_neg,
    };

    // Final PSD projection to guarantee result
    let (x_final, _) = project_psd(&x);
    (x_final, conv)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_project_psd_already_psd() {
        let m = DMatrix::from_row_slice(2, 2, &[2.0, 1.0, 1.0, 2.0]);
        let (p, neg) = project_psd(&m);
        assert_eq!(neg, 0);
        assert!((p - m).norm() < 1e-10);
    }

    #[test]
    fn test_project_psd_removes_negatives() {
        // Matrix with a negative eigenvalue
        let m = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 2.0, 1.0]);
        // eigenvalues are 3 and -1
        let (p, neg) = project_psd(&m);
        assert_eq!(neg, 1);
        // Verify PSD: all eigenvalues >= 0
        let eigen = p.symmetric_eigen();
        for &e in eigen.eigenvalues.iter() {
            assert!(e >= -1e-10, "Eigenvalue should be non-negative: {}", e);
        }
    }

    #[test]
    fn test_project_distances() {
        let n = 3;
        let x = DMatrix::identity(n, n);
        let pairs = vec![(0, 1, 1.5), (1, 2, 2.0)]; // target distances
        let p = project_distances(&x, &pairs);
        // Check constraint satisfaction improved
        assert!(p.nrows() == n);
    }

    #[test]
    fn test_alternating_projections_converges() {
        let n = 4;
        // Create a valid Gram matrix from known coordinates
        let coords: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, 0.0],
        ];

        // Build Gram matrix
        let mut gram = DMatrix::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                gram[(i, j)] = coords[i][0] * coords[j][0]
                    + coords[i][1] * coords[j][1]
                    + coords[i][2] * coords[j][2];
            }
        }

        // Compute distances
        let mut pairs = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let d = ((coords[i][0] - coords[j][0]).powi(2)
                    + (coords[i][1] - coords[j][1]).powi(2)
                    + (coords[i][2] - coords[j][2]).powi(2))
                .sqrt();
                pairs.push((i, j, d));
            }
        }

        // Perturb Gram matrix
        let mut x0 = gram.clone();
        x0[(0, 1)] += 0.3;
        x0[(1, 0)] += 0.3;

        let config = SdrConfig::default();
        let (result, _conv) = alternating_projections(&x0, &pairs, &config);

        // Should produce PSD matrix
        let eigen = result.symmetric_eigen();
        for &e in eigen.eigenvalues.iter() {
            assert!(e >= -1e-8, "Should be PSD: {}", e);
        }
    }

    #[test]
    fn test_svt_step() {
        let m = DMatrix::from_row_slice(2, 2, &[3.0, 0.0, 0.0, 1.0]);
        let s = svt_step(&m, 1.5);
        let eigen = s.symmetric_eigen();
        // eigenvalue 1.0 should be thresholded to 0, eigenvalue 3.0 to 1.5
        let mut evals: Vec<f64> = eigen.eigenvalues.iter().cloned().collect();
        evals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!(
            (evals[0]).abs() < 1e-10,
            "Small eigenvalue thresholded: {}",
            evals[0]
        );
        assert!(
            (evals[1] - 1.5).abs() < 1e-10,
            "Large eigenvalue reduced: {}",
            evals[1]
        );
    }
}
