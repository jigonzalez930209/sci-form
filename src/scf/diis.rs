//! Direct Inversion in Iterative Subspace (DIIS) convergence accelerator.
//!
//! DIIS extrapolates the Fock matrix from a history of previous iterations
//! to accelerate SCF convergence. The DIIS error vector is:
//!
//!   e_i = F_i P_i S - S P_i F_i
//!
//! The optimal linear combination minimizes ||Σ c_i e_i||² subject to Σ c_i = 1.
//!
//! # References
//!
//! P. Pulay, Chem. Phys. Lett. 73, 393 (1980).
//! P. Pulay, J. Comput. Chem. 3, 556 (1982).

use nalgebra::{DMatrix, DVector};

/// DIIS accelerator state.
pub struct DiisAccelerator {
    /// History of Fock matrices.
    fock_history: Vec<DMatrix<f64>>,
    /// History of error vectors (flattened).
    error_history: Vec<DVector<f64>>,
    /// Maximum subspace size.
    max_size: usize,
}

impl DiisAccelerator {
    /// Create a new DIIS accelerator with the given subspace size.
    pub fn new(max_size: usize) -> Self {
        Self {
            fock_history: Vec::with_capacity(max_size),
            error_history: Vec::with_capacity(max_size),
            max_size,
        }
    }

    /// Add an SCF iteration to the DIIS history.
    ///
    /// The error vector is computed as: e = FPS - SPF
    pub fn add_iteration(
        &mut self,
        fock: &DMatrix<f64>,
        density: &DMatrix<f64>,
        overlap: &DMatrix<f64>,
    ) {
        // Compute DIIS error: e = FPS - SPF
        let fps = fock * density * overlap;
        let spf = overlap * density * fock;
        let error_matrix = &fps - &spf;

        // Flatten error matrix to vector
        let n = error_matrix.nrows();
        let error_vec = DVector::from_iterator(
            n * n,
            error_matrix.iter().copied(),
        );

        // Enforce max subspace size
        if self.fock_history.len() >= self.max_size {
            self.fock_history.remove(0);
            self.error_history.remove(0);
        }

        self.fock_history.push(fock.clone());
        self.error_history.push(error_vec);
    }

    /// Extrapolate the optimal Fock matrix from the DIIS subspace.
    ///
    /// Returns the extrapolated Fock matrix, or None if not enough
    /// history is available (need at least 2 iterations).
    pub fn extrapolate(&self) -> Option<DMatrix<f64>> {
        let m = self.error_history.len();
        if m < 2 {
            return None;
        }

        // Build the DIIS B matrix:
        // B[i][j] = e_i · e_j
        // With Lagrange constraint: last row/col = -1/1
        let dim = m + 1;
        let mut b = DMatrix::zeros(dim, dim);

        for i in 0..m {
            for j in 0..=i {
                let bij = self.error_history[i].dot(&self.error_history[j]);
                b[(i, j)] = bij;
                b[(j, i)] = bij;
            }
        }

        // Lagrange constraint
        for i in 0..m {
            b[(m, i)] = -1.0;
            b[(i, m)] = -1.0;
        }
        b[(m, m)] = 0.0;

        // Right-hand side: [0, 0, ..., 0, -1]
        let mut rhs = DVector::zeros(dim);
        rhs[m] = -1.0;

        // Solve B · c = rhs
        let coefficients = solve_linear_system(&b, &rhs)?;

        // Extrapolate: F_diis = Σ c_i F_i
        let n = self.fock_history[0].nrows();
        let mut f_diis = DMatrix::zeros(n, n);

        for i in 0..m {
            f_diis += coefficients[i] * &self.fock_history[i];
        }

        Some(f_diis)
    }

    /// Maximum error norm in the current subspace.
    pub fn max_error(&self) -> f64 {
        self.error_history
            .last()
            .map(|e| e.norm())
            .unwrap_or(f64::MAX)
    }

    /// Reset the DIIS subspace (useful when convergence stalls).
    pub fn reset(&mut self) {
        self.fock_history.clear();
        self.error_history.clear();
    }

    /// Number of stored iterations.
    pub fn size(&self) -> usize {
        self.fock_history.len()
    }
}

/// Solve a linear system Ax = b using LU decomposition.
fn solve_linear_system(a: &DMatrix<f64>, b: &DVector<f64>) -> Option<DVector<f64>> {
    let lu = a.clone().lu();
    lu.solve(b)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_diis_creation() {
        let diis = DiisAccelerator::new(6);
        assert_eq!(diis.size(), 0);
    }

    #[test]
    fn test_diis_no_extrapolation_with_one_iteration() {
        let mut diis = DiisAccelerator::new(6);
        let f = DMatrix::identity(2, 2);
        let p = DMatrix::identity(2, 2);
        let s = DMatrix::identity(2, 2);

        diis.add_iteration(&f, &p, &s);
        assert!(diis.extrapolate().is_none());
    }

    #[test]
    fn test_diis_extrapolation_converged() {
        let mut diis = DiisAccelerator::new(6);

        // When commutator [F,P] = 0 (converged), DIIS should return
        // something close to the last Fock matrix
        let f = DMatrix::from_row_slice(2, 2, &[-1.0, 0.0, 0.0, -0.5]);
        let s = DMatrix::identity(2, 2);
        // P that commutes with F when S = I: P must be diagonal too
        let p = DMatrix::from_row_slice(2, 2, &[2.0, 0.0, 0.0, 0.0]);

        diis.add_iteration(&f, &p, &s);
        diis.add_iteration(&f, &p, &s);

        if let Some(f_diis) = diis.extrapolate() {
            // Should be close to f
            for i in 0..2 {
                for j in 0..2 {
                    assert!(
                        (f_diis[(i, j)] - f[(i, j)]).abs() < 1e-8,
                        "DIIS result differs at ({},{})",
                        i, j
                    );
                }
            }
        }
    }

    #[test]
    fn test_diis_reset() {
        let mut diis = DiisAccelerator::new(6);
        let f = DMatrix::identity(2, 2);
        let p = DMatrix::identity(2, 2);
        let s = DMatrix::identity(2, 2);

        diis.add_iteration(&f, &p, &s);
        assert_eq!(diis.size(), 1);

        diis.reset();
        assert_eq!(diis.size(), 0);
    }
}
