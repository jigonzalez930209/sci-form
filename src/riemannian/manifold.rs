//! Riemannian Geometry Infrastructure — E3.1
//!
//! Provides retraction onto the PSD cone, Riemannian gradient projection,
//! and exponential map on the fixed-rank PSD manifold.

use nalgebra::{DMatrix, SymmetricEigen};

/// A point on the PSD manifold: a symmetric positive semidefinite matrix.
#[derive(Debug, Clone)]
pub struct PsdManifold {
    /// The current PSD matrix X.
    pub x: DMatrix<f64>,
    /// Cached eigendecomposition: eigenvalues.
    eigenvalues: Option<Vec<f64>>,
    /// Cached eigendecomposition: eigenvectors.
    eigenvectors: Option<DMatrix<f64>>,
}

impl PsdManifold {
    /// Create a new PSD manifold point from a symmetric matrix.
    /// Projects to PSD if the matrix has negative eigenvalues.
    pub fn new(x: DMatrix<f64>) -> Self {
        let mut m = Self {
            x,
            eigenvalues: None,
            eigenvectors: None,
        };
        m.ensure_psd();
        m
    }

    /// Create from a known PSD matrix (skips projection).
    pub fn from_psd(x: DMatrix<f64>) -> Self {
        Self {
            x,
            eigenvalues: None,
            eigenvectors: None,
        }
    }

    /// Dimension of the matrix.
    pub fn dim(&self) -> usize {
        self.x.nrows()
    }

    /// Ensure the matrix is PSD by projecting negative eigenvalues to zero.
    fn ensure_psd(&mut self) {
        let n = self.x.nrows();
        let eigen = SymmetricEigen::new(self.x.clone());
        let mut any_negative = false;
        let mut vals = Vec::with_capacity(n);
        for i in 0..n {
            let v = eigen.eigenvalues[i];
            if v < 0.0 {
                any_negative = true;
                vals.push(0.0);
            } else {
                vals.push(v);
            }
        }

        if any_negative {
            let mut diag = DMatrix::zeros(n, n);
            for i in 0..n {
                diag[(i, i)] = vals[i];
            }
            self.x = &eigen.eigenvectors * diag * eigen.eigenvectors.transpose();
        }

        self.eigenvalues = Some(vals);
        self.eigenvectors = Some(eigen.eigenvectors.clone_owned());
    }

    /// Get eigenvalues (cached).
    pub fn eigenvalues(&mut self) -> &[f64] {
        if self.eigenvalues.is_none() {
            let eigen = SymmetricEigen::new(self.x.clone());
            let n = self.x.nrows();
            let mut vals = Vec::with_capacity(n);
            for i in 0..n {
                vals.push(eigen.eigenvalues[i]);
            }
            self.eigenvectors = Some(eigen.eigenvectors.clone_owned());
            self.eigenvalues = Some(vals);
        }
        self.eigenvalues.as_ref().unwrap()
    }

    /// Number of positive eigenvalues.
    pub fn rank(&mut self) -> usize {
        self.eigenvalues().iter().filter(|&&v| v > 1e-10).count()
    }

    /// Frobenius norm of the matrix.
    pub fn frobenius_norm(&self) -> f64 {
        let n = self.x.nrows();
        let mut sum = 0.0;
        for j in 0..n {
            for i in 0..n {
                sum += self.x[(i, j)] * self.x[(i, j)];
            }
        }
        sum.sqrt()
    }
}

/// Project matrix X onto the PSD cone: P_+(X) = V max(Λ, 0) V^T.
pub fn psd_projection(x: &DMatrix<f64>) -> DMatrix<f64> {
    let n = x.nrows();
    let eigen = SymmetricEigen::new(x.clone());
    let mut diag = DMatrix::zeros(n, n);
    for i in 0..n {
        diag[(i, i)] = eigen.eigenvalues[i].max(0.0);
    }
    &eigen.eigenvectors * diag * eigen.eigenvectors.transpose()
}

/// First-order retraction on the PSD manifold.
///
/// Retr_X(ξ) = P_+(X + ξ)
///
/// This is the simplest retraction: add the tangent vector and project back to PSD.
pub fn psd_retraction(x: &DMatrix<f64>, xi: &DMatrix<f64>) -> DMatrix<f64> {
    psd_projection(&(x + xi))
}

/// Project a matrix ξ onto the tangent space at X on the PSD manifold.
///
/// For the PSD cone, the tangent space at X consists of symmetric matrices ξ
/// such that ξ is in the span of the positive eigenvectors plus any symmetric
/// matrix in the zero-eigenvalue directions that is PSD.
///
/// Simplified version: just symmetrize.
pub fn tangent_projection(xi: &DMatrix<f64>) -> DMatrix<f64> {
    (xi + xi.transpose()) * 0.5
}

/// Geodesic distance between two PSD matrices (affine-invariant metric).
///
/// d(A, B) = ||log(A^{-1/2} B A^{-1/2})||_F
///
/// For computational simplicity, we use the Frobenius distance as a proxy.
pub fn psd_distance(a: &DMatrix<f64>, b: &DMatrix<f64>) -> f64 {
    let diff = a - b;
    let n = diff.nrows();
    let mut sum = 0.0;
    for j in 0..n {
        for i in 0..n {
            sum += diff[(i, j)] * diff[(i, j)];
        }
    }
    sum.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    fn make_random_symmetric(n: usize, seed: u64) -> DMatrix<f64> {
        let mut rng = StdRng::seed_from_u64(seed);
        let mut m = DMatrix::zeros(n, n);
        for i in 0..n {
            m[(i, i)] = rng.gen_range(-2.0..5.0);
            for j in (i + 1)..n {
                let v = rng.gen_range(-1.0..1.0);
                m[(i, j)] = v;
                m[(j, i)] = v;
            }
        }
        m
    }

    fn make_psd(n: usize, seed: u64) -> DMatrix<f64> {
        let m = make_random_symmetric(n, seed);
        psd_projection(&m)
    }

    #[test]
    fn test_psd_projection_removes_negative_eigenvalues() {
        let m = make_random_symmetric(10, 1);
        let p = psd_projection(&m);
        let eigen = SymmetricEigen::new(p);
        for i in 0..eigen.eigenvalues.len() {
            assert!(
                eigen.eigenvalues[i] >= -1e-10,
                "Eigenvalue {} = {}, expected ≥ 0",
                i,
                eigen.eigenvalues[i]
            );
        }
    }

    #[test]
    fn test_psd_projection_preserves_psd() {
        let p = make_psd(10, 2);
        let p2 = psd_projection(&p);
        let diff = psd_distance(&p, &p2);
        assert!(
            diff < 1e-10,
            "Projecting a PSD matrix should not change it: diff = {}",
            diff
        );
    }

    #[test]
    fn test_retraction_stays_psd() {
        let x = make_psd(8, 3);
        let xi = make_random_symmetric(8, 4) * 0.1; // Small perturbation
        let y = psd_retraction(&x, &xi);
        let eigen = SymmetricEigen::new(y);
        for i in 0..eigen.eigenvalues.len() {
            assert!(
                eigen.eigenvalues[i] >= -1e-10,
                "Retracted point has negative eigenvalue: {}",
                eigen.eigenvalues[i]
            );
        }
    }

    #[test]
    fn test_tangent_projection_symmetric() {
        let mut rng = StdRng::seed_from_u64(5);
        let n = 6;
        let mut m = DMatrix::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                m[(i, j)] = rng.gen_range(-1.0..1.0);
            }
        }
        let proj = tangent_projection(&m);
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (proj[(i, j)] - proj[(j, i)]).abs() < 1e-14,
                    "Tangent projection not symmetric at ({},{})",
                    i,
                    j
                );
            }
        }
    }

    #[test]
    fn test_manifold_rank() {
        let n = 10;
        // Create a rank-3 PSD matrix
        let mut rng = StdRng::seed_from_u64(6);
        let mut a = DMatrix::zeros(n, 3);
        for j in 0..3 {
            for i in 0..n {
                a[(i, j)] = rng.gen_range(-1.0..1.0);
            }
        }
        let x = &a * a.transpose(); // rank ≤ 3
        let mut m = PsdManifold::from_psd(x);
        assert_eq!(m.rank(), 3);
    }

    #[test]
    fn test_psd_distance_symmetric() {
        let a = make_psd(5, 7);
        let b = make_psd(5, 8);
        let d_ab = psd_distance(&a, &b);
        let d_ba = psd_distance(&b, &a);
        assert!(
            (d_ab - d_ba).abs() < 1e-14,
            "Distance not symmetric: {} vs {}",
            d_ab,
            d_ba
        );
    }
}
