//! Nyström Approximation Infrastructure — E2.1
//!
//! Provides Gaussian sketch matrix generation and low-rank matrix
//! approximation via the Nyström method.

use nalgebra::{DMatrix, DVector};
use rand::Rng;

/// A Gaussian sketch matrix Ω ∈ R^{N×k}.
#[derive(Debug, Clone)]
pub struct GaussianSketch {
    /// The sketch matrix.
    pub omega: DMatrix<f64>,
    /// Number of rows (original dimension).
    pub n: usize,
    /// Number of columns (sketch size).
    pub k: usize,
}

impl GaussianSketch {
    /// Generate a Gaussian sketch matrix with entries ~ N(0, 1/k).
    pub fn new<R: Rng>(rng: &mut R, n: usize, k: usize) -> Self {
        let scale = 1.0 / (n as f64).sqrt();
        let mut omega = DMatrix::zeros(n, k);
        for j in 0..k {
            for i in 0..n {
                // Box-Muller transform for standard normal
                let u1: f64 = rng.gen_range(1e-10..1.0);
                let u2: f64 = rng.gen_range(0.0..std::f64::consts::TAU);
                let z = (-2.0 * u1.ln()).sqrt() * u2.cos();
                omega[(i, j)] = z * scale;
            }
        }
        Self { omega, n, k }
    }

    /// Default sketch size: k = min(N, 50 + ceil(log2(N))).
    pub fn default_k(n: usize) -> usize {
        n.min(50 + ((n as f64).log2().ceil() as usize))
    }
}

/// Low-rank Nyström approximation of a symmetric PSD matrix.
#[derive(Debug, Clone)]
pub struct NystromApprox {
    /// Orthogonal factor Q from thin QR of Y = SΩ.
    pub q: DMatrix<f64>,
    /// Small k×k matrix B = Q^T S Q.
    pub b: DMatrix<f64>,
    /// Approximation: S ≈ Q B Q^T.
    /// Original matrix dimension.
    pub n: usize,
    /// Sketch dimension.
    pub k: usize,
}

impl NystromApprox {
    /// Build Nyström approximation of symmetric matrix S using sketch Ω.
    ///
    /// 1. Y = S Ω
    /// 2. QR decomposition: Y = Q R
    /// 3. B = Q^T S Q
    pub fn from_matrix(s: &DMatrix<f64>, sketch: &GaussianSketch) -> Self {
        let n = s.nrows();
        let k = sketch.k;

        // Y = S * Ω
        let y = s * &sketch.omega;

        // Thin QR: Y = Q R
        let qr = y.clone().qr();
        let q = qr.q();
        // Keep only k columns
        let q = q.columns(0, k.min(q.ncols())).into_owned();
        let actual_k = q.ncols();

        // B = Q^T S Q (k × k)
        let b = q.transpose() * s * &q;

        Self {
            q,
            b,
            n,
            k: actual_k,
        }
    }

    /// Reconstruct the full N×N approximation: S_hat = Q B Q^T.
    pub fn reconstruct(&self) -> DMatrix<f64> {
        &self.q * &self.b * self.q.transpose()
    }

    /// Compute Frobenius error vs the original matrix.
    pub fn frobenius_error(&self, s: &DMatrix<f64>) -> f64 {
        let approx = self.reconstruct();
        let diff = s - &approx;
        let mut sum = 0.0;
        for j in 0..diff.ncols() {
            for i in 0..diff.nrows() {
                sum += diff[(i, j)] * diff[(i, j)];
            }
        }
        sum.sqrt()
    }

    /// Compute relative Frobenius error: ||S - S_hat||_F / ||S||_F.
    pub fn relative_error(&self, s: &DMatrix<f64>) -> f64 {
        let err = self.frobenius_error(s);
        let mut s_norm = 0.0;
        for j in 0..s.ncols() {
            for i in 0..s.nrows() {
                s_norm += s[(i, j)] * s[(i, j)];
            }
        }
        err / s_norm.sqrt().max(1e-15)
    }

    /// Compute low-rank eigendecomposition of B.
    ///
    /// Returns eigenvalues and eigenvectors in the ORIGINAL space:
    /// If B = V Λ V^T, then the approximate eigenvectors of S are Q V.
    pub fn eigendecompose(&self) -> (DVector<f64>, DMatrix<f64>) {
        use nalgebra::SymmetricEigen;
        let eigen = SymmetricEigen::new(self.b.clone());
        let eigenvalues = eigen.eigenvalues;
        let eigenvectors_small = eigen.eigenvectors;

        // Map back to full space: U = Q V
        let eigenvectors = &self.q * eigenvectors_small;

        (eigenvalues, eigenvectors)
    }

    /// Compute S^{-1/2} approximation via low-rank factorization.
    ///
    /// S ≈ Q B Q^T + (I - Q Q^T)  (complement is ~I for overlap matrices)
    /// → S^{-1/2} ≈ Q B^{-1/2} Q^T + (I - Q Q^T)
    ///
    /// The complement term ensures directions not captured by the sketch
    /// are treated as identity, which is valid for overlap matrices near I.
    pub fn inverse_sqrt(&self) -> DMatrix<f64> {
        use nalgebra::SymmetricEigen;
        let eigen = SymmetricEigen::new(self.b.clone());
        let vals = &eigen.eigenvalues;
        let vecs = &eigen.eigenvectors;

        // B^{-1/2} = V diag(1/sqrt(λ)) V^T, skip near-zero eigenvalues
        let k = vals.len();
        let mut diag_inv_sqrt = DMatrix::zeros(k, k);
        for i in 0..k {
            if vals[i] > 1e-10 {
                diag_inv_sqrt[(i, i)] = 1.0 / vals[i].sqrt();
            }
        }

        let b_inv_sqrt = vecs * &diag_inv_sqrt * vecs.transpose();

        // S^{-1/2} ≈ Q (B^{-1/2} - I_k) Q^T + I_n
        // Equivalent to: Q B^{-1/2} Q^T in Q-space, I in complement
        let delta = &b_inv_sqrt - DMatrix::identity(k, k);
        DMatrix::identity(self.n, self.n) + &self.q * delta * self.q.transpose()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn make_spd_matrix(n: usize) -> DMatrix<f64> {
        // Overlap-like SPD matrix with spectral decay
        let mut rng = StdRng::seed_from_u64(42);
        let mut s = DMatrix::identity(n, n);
        for i in 0..n {
            for j in (i + 1)..n {
                let decay = (-0.1 * (j - i) as f64).exp();
                let off = rng.gen_range(-0.3..0.3) * decay;
                s[(i, j)] = off;
                s[(j, i)] = off;
            }
        }
        s
    }

    #[test]
    fn test_gaussian_sketch_dimensions() {
        let mut rng = StdRng::seed_from_u64(1);
        let sketch = GaussianSketch::new(&mut rng, 20, 5);
        assert_eq!(sketch.n, 20);
        assert_eq!(sketch.k, 5);
        assert_eq!(sketch.omega.nrows(), 20);
        assert_eq!(sketch.omega.ncols(), 5);
    }

    #[test]
    fn test_sketch_approximate_orthogonality() {
        let mut rng = StdRng::seed_from_u64(2);
        let n = 200;
        let k = 10;
        let sketch = GaussianSketch::new(&mut rng, n, k);
        let gram = sketch.omega.transpose() * &sketch.omega;
        // E[Omega^T Omega] ≈ I_k with tolerance O(1/sqrt(N))
        let tol = 3.0 / (n as f64).sqrt();
        for i in 0..k {
            assert!(
                (gram[(i, i)] - 1.0).abs() < tol,
                "Diagonal ({},{}) = {}, expected ~1 ± {}",
                i, i, gram[(i, i)], tol
            );
        }
    }

    #[test]
    fn test_nystrom_reconstruction_error() {
        let n = 30;
        let k = 28; // Near-full rank for near-identity matrix
        let s = make_spd_matrix(n);
        let mut rng = StdRng::seed_from_u64(3);
        let sketch = GaussianSketch::new(&mut rng, n, k);
        let approx = NystromApprox::from_matrix(&s, &sketch);
        let rel_err = approx.relative_error(&s);
        assert!(
            rel_err < 0.50,
            "Relative Frobenius error = {}, expected < 50%",
            rel_err
        );
    }

    #[test]
    fn test_nystrom_eigendecompose() {
        let n = 20;
        let k = 10;
        let s = make_spd_matrix(n);
        let mut rng = StdRng::seed_from_u64(4);
        let sketch = GaussianSketch::new(&mut rng, n, k);
        let approx = NystromApprox::from_matrix(&s, &sketch);
        let (vals, _vecs) = approx.eigendecompose();
        // All eigenvalues should be non-negative for SPD matrix
        for i in 0..vals.len() {
            assert!(
                vals[i] > -0.1,
                "Eigenvalue {} = {} should be non-negative",
                i, vals[i]
            );
        }
    }

    #[test]
    fn test_inverse_sqrt() {
        let n = 15;
        let k = 12;
        let s = make_spd_matrix(n);
        let mut rng = StdRng::seed_from_u64(5);
        let sketch = GaussianSketch::new(&mut rng, n, k);
        let approx = NystromApprox::from_matrix(&s, &sketch);
        let s_inv_sqrt = approx.inverse_sqrt();

        // S^{-1/2} S S^{-1/2} ≈ I (for good approximation)
        let product = &s_inv_sqrt * &s * &s_inv_sqrt;
        for i in 0..n {
            assert!(
                (product[(i, i)] - 1.0).abs() < 0.3,
                "Diagonal ({},{}) = {}, expected ~1",
                i, i, product[(i, i)]
            );
        }
    }

    #[test]
    fn test_default_k() {
        assert_eq!(GaussianSketch::default_k(10), 10); // min(10, 54) = 10
        assert_eq!(GaussianSketch::default_k(100), 57); // 50 + ceil(log2(100)) = 50 + 7
        assert_eq!(GaussianSketch::default_k(1000), 60); // 50 + ceil(log2(1000)) = 50 + 10
    }
}
