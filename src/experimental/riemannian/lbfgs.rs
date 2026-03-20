//! Riemannian L-BFGS Optimizer — E3.2
//!
//! L-BFGS on the PSD manifold for distance geometry embedding.
//! Uses retraction and tangent space operations to enforce PSD constraints.

use nalgebra::DMatrix;

use super::manifold::{psd_retraction, tangent_projection};

/// Configuration for the Riemannian L-BFGS optimizer.
#[derive(Debug, Clone)]
pub struct RiemannianConfig {
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Number of L-BFGS history pairs to store.
    pub history_size: usize,
    /// Gradient norm convergence tolerance.
    pub grad_tol: f64,
    /// Maximum per-pair distance violation (Å) for convergence.
    pub dist_tol: f64,
    /// Armijo constant for line search.
    pub armijo_c1: f64,
    /// Maximum line search steps.
    pub max_ls_steps: usize,
    /// Whether to fall back to Euclidean BFGS on failure.
    pub fallback_enabled: bool,
}

impl Default for RiemannianConfig {
    fn default() -> Self {
        Self {
            max_iter: 500,
            history_size: 5,
            grad_tol: 1e-6,
            dist_tol: 0.01,
            armijo_c1: 1e-4,
            max_ls_steps: 20,
            fallback_enabled: true,
        }
    }
}

/// Result of a Riemannian L-BFGS optimization.
#[derive(Debug, Clone)]
pub struct RiemannianResult {
    /// Final metric matrix (PSD).
    pub metric_matrix: DMatrix<f64>,
    /// Final objective function value.
    pub objective: f64,
    /// Final gradient norm.
    pub grad_norm: f64,
    /// Number of iterations used.
    pub iterations: usize,
    /// Whether the optimizer converged.
    pub converged: bool,
    /// Whether it fell back to Euclidean optimization.
    pub used_fallback: bool,
    /// Maximum per-pair distance violation.
    pub max_violation: f64,
}

/// Riemannian L-BFGS optimizer on the PSD manifold.
pub struct RiemannianLbfgs {
    config: RiemannianConfig,
}

/// A distance constraint: atoms i,j should have distance d.
#[derive(Debug, Clone, Copy)]
pub struct DistanceConstraint {
    pub i: usize,
    pub j: usize,
    pub lower: f64,
    pub upper: f64,
}

impl RiemannianLbfgs {
    pub fn new(config: RiemannianConfig) -> Self {
        Self { config }
    }

    /// Optimize a metric matrix to satisfy distance constraints.
    ///
    /// The objective is: f(X) = sum_pairs max(0, d_ij² - X_ij)² + max(0, X_ij - u_ij²)²
    /// where X_ij = X_ii + X_jj - 2*X_ij is the squared distance from the metric matrix.
    pub fn optimize(
        &self,
        initial: &DMatrix<f64>,
        constraints: &[DistanceConstraint],
    ) -> RiemannianResult {
        let mut x = initial.clone();

        // L-BFGS history
        let m = self.config.history_size;
        let mut s_hist: Vec<DMatrix<f64>> = Vec::with_capacity(m);
        let mut y_hist: Vec<DMatrix<f64>> = Vec::with_capacity(m);
        let mut rho_hist: Vec<f64> = Vec::with_capacity(m);

        let mut obj = self.objective(&x, constraints);
        let mut grad = self.gradient(&x, constraints);
        let mut grad_norm = frobenius_norm(&grad);

        let mut converged = false;
        let mut iter = 0;

        while iter < self.config.max_iter {
            if grad_norm < self.config.grad_tol {
                converged = true;
                break;
            }

            // L-BFGS two-loop recursion to compute search direction
            let direction = if s_hist.is_empty() {
                // Steepest descent on first iteration
                scale_matrix(&grad, -1.0)
            } else {
                self.lbfgs_direction(&grad, &s_hist, &y_hist, &rho_hist)
            };

            // Project direction to tangent space
            let direction = tangent_projection(&direction);

            // Armijo line search along retraction
            let (alpha, new_x, new_obj) =
                self.line_search(&x, &direction, obj, &grad, constraints);

            if alpha < 1e-16 {
                // Line search failed — try steepest descent
                let sd = scale_matrix(&grad, -1.0);
                let (alpha2, new_x2, new_obj2) =
                    self.line_search(&x, &sd, obj, &grad, constraints);
                if alpha2 < 1e-16 {
                    break; // Give up
                }
                let new_grad = self.gradient(&new_x2, constraints);
                let sk = &new_x2 - &x;
                let yk = &new_grad - &grad;
                let rho = 1.0 / inner_product(&yk, &sk).max(1e-15);

                if inner_product(&yk, &sk) > 0.0 {
                    if s_hist.len() >= m {
                        s_hist.remove(0);
                        y_hist.remove(0);
                        rho_hist.remove(0);
                    }
                    s_hist.push(sk);
                    y_hist.push(yk);
                    rho_hist.push(rho);
                }

                x = new_x2;
                obj = new_obj2;
                grad = new_grad;
                grad_norm = frobenius_norm(&grad);
                iter += 1;
                continue;
            }

            let new_grad = self.gradient(&new_x, constraints);

            // Update L-BFGS history
            let sk = &new_x - &x;
            let yk = &new_grad - &grad;
            let sy = inner_product(&yk, &sk);

            if sy > 0.0 {
                if s_hist.len() >= m {
                    s_hist.remove(0);
                    y_hist.remove(0);
                    rho_hist.remove(0);
                }
                s_hist.push(sk);
                y_hist.push(yk);
                rho_hist.push(1.0 / sy);
            }

            x = new_x;
            obj = new_obj;
            grad = new_grad;
            grad_norm = frobenius_norm(&grad);
            iter += 1;
        }

        let max_viol = self.max_violation(&x, constraints);

        RiemannianResult {
            metric_matrix: x,
            objective: obj,
            grad_norm,
            iterations: iter,
            converged,
            used_fallback: false,
            max_violation: max_viol,
        }
    }

    /// Distance geometry objective function.
    pub fn objective(&self, x: &DMatrix<f64>, constraints: &[DistanceConstraint]) -> f64 {
        let mut f = 0.0;
        for c in constraints {
            let dij_sq = x[(c.i, c.i)] + x[(c.j, c.j)] - 2.0 * x[(c.i, c.j)];
            let lower_sq = c.lower * c.lower;
            let upper_sq = c.upper * c.upper;
            if dij_sq < lower_sq {
                let v = lower_sq - dij_sq;
                f += v * v;
            } else if dij_sq > upper_sq {
                let v = dij_sq - upper_sq;
                f += v * v;
            }
        }
        f
    }

    /// Gradient of the objective with respect to the metric matrix.
    pub fn gradient(&self, x: &DMatrix<f64>, constraints: &[DistanceConstraint]) -> DMatrix<f64> {
        let n = x.nrows();
        let mut g = DMatrix::zeros(n, n);
        for c in constraints {
            let dij_sq = x[(c.i, c.i)] + x[(c.j, c.j)] - 2.0 * x[(c.i, c.j)];
            let lower_sq = c.lower * c.lower;
            let upper_sq = c.upper * c.upper;

            let factor = if dij_sq < lower_sq {
                -2.0 * (lower_sq - dij_sq)
            } else if dij_sq > upper_sq {
                2.0 * (dij_sq - upper_sq)
            } else {
                0.0
            };

            if factor.abs() > 0.0 {
                g[(c.i, c.i)] += factor;
                g[(c.j, c.j)] += factor;
                g[(c.i, c.j)] -= 2.0 * factor;
                g[(c.j, c.i)] -= 2.0 * factor;
            }
        }
        // Project to tangent space (symmetrize)
        tangent_projection(&g)
    }

    /// Maximum per-pair distance violation.
    fn max_violation(&self, x: &DMatrix<f64>, constraints: &[DistanceConstraint]) -> f64 {
        let mut max_v = 0.0f64;
        for c in constraints {
            let dij_sq = x[(c.i, c.i)] + x[(c.j, c.j)] - 2.0 * x[(c.i, c.j)];
            let dij = dij_sq.max(0.0).sqrt();
            let v = if dij < c.lower {
                c.lower - dij
            } else if dij > c.upper {
                dij - c.upper
            } else {
                0.0
            };
            max_v = max_v.max(v);
        }
        max_v
    }

    /// L-BFGS two-loop recursion.
    fn lbfgs_direction(
        &self,
        grad: &DMatrix<f64>,
        s_hist: &[DMatrix<f64>],
        y_hist: &[DMatrix<f64>],
        rho_hist: &[f64],
    ) -> DMatrix<f64> {
        let k = s_hist.len();
        let mut q = grad.clone();
        let mut alphas = vec![0.0; k];

        // Right-to-left loop
        for i in (0..k).rev() {
            alphas[i] = rho_hist[i] * inner_product(&s_hist[i], &q);
            q = &q - &scale_matrix(&y_hist[i], alphas[i]);
        }

        // Initial Hessian approximation: H_0 = γ I (projected)
        let gamma = if k > 0 {
            inner_product(&s_hist[k - 1], &y_hist[k - 1])
                / inner_product(&y_hist[k - 1], &y_hist[k - 1]).max(1e-15)
        } else {
            1.0
        };
        let mut r = scale_matrix(&q, gamma);

        // Left-to-right loop
        for i in 0..k {
            let beta = rho_hist[i] * inner_product(&y_hist[i], &r);
            r = &r + &scale_matrix(&s_hist[i], alphas[i] - beta);
        }

        // Negate for descent direction
        scale_matrix(&r, -1.0)
    }

    /// Armijo backtracking line search along manifold retraction.
    fn line_search(
        &self,
        x: &DMatrix<f64>,
        direction: &DMatrix<f64>,
        obj: f64,
        grad: &DMatrix<f64>,
        constraints: &[DistanceConstraint],
    ) -> (f64, DMatrix<f64>, f64) {
        let dg = inner_product(grad, direction);
        if dg >= 0.0 {
            // Not a descent direction
            return (0.0, x.clone(), obj);
        }

        let mut alpha = 1.0;
        for _ in 0..self.config.max_ls_steps {
            let step = scale_matrix(direction, alpha);
            let x_new = psd_retraction(x, &step);
            let obj_new = self.objective(&x_new, constraints);

            if obj_new <= obj + self.config.armijo_c1 * alpha * dg {
                return (alpha, x_new, obj_new);
            }

            alpha *= 0.5;
        }

        // Failed
        (0.0, x.clone(), obj)
    }
}

/// Frobenius inner product: Tr(A^T B).
fn inner_product(a: &DMatrix<f64>, b: &DMatrix<f64>) -> f64 {
    let n = a.nrows();
    let mut sum = 0.0;
    for j in 0..n {
        for i in 0..n {
            sum += a[(i, j)] * b[(i, j)];
        }
    }
    sum
}

/// Frobenius norm.
fn frobenius_norm(m: &DMatrix<f64>) -> f64 {
    inner_product(m, m).sqrt()
}

/// Scale a matrix by a scalar.
fn scale_matrix(m: &DMatrix<f64>, s: f64) -> DMatrix<f64> {
    m * s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_simple_constraints() -> Vec<DistanceConstraint> {
        // Triangle with side lengths 1.5, 1.5, 1.5
        vec![
            DistanceConstraint { i: 0, j: 1, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 0, j: 2, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 1, j: 2, lower: 1.4, upper: 1.6 },
        ]
    }

    #[test]
    fn test_optimizer_converges_simple() {
        let constraints = make_simple_constraints();
        let initial = DMatrix::identity(3, 3) * 2.0;
        let config = RiemannianConfig {
            max_iter: 200,
            ..Default::default()
        };
        let optimizer = RiemannianLbfgs::new(config);
        let result = optimizer.optimize(&initial, &constraints);

        assert!(
            result.converged || result.objective < 1e-4,
            "Should converge for simple triangle: obj={}, grad_norm={}, iter={}",
            result.objective, result.grad_norm, result.iterations
        );
    }

    #[test]
    fn test_optimizer_satisfies_constraints() {
        let constraints = make_simple_constraints();
        let initial = DMatrix::identity(3, 3) * 3.0;
        let config = RiemannianConfig::default();
        let optimizer = RiemannianLbfgs::new(config);
        let result = optimizer.optimize(&initial, &constraints);

        // Check that all distances are within bounds
        for c in &constraints {
            let dij_sq = result.metric_matrix[(c.i, c.i)]
                + result.metric_matrix[(c.j, c.j)]
                - 2.0 * result.metric_matrix[(c.i, c.j)];
            let dij = dij_sq.max(0.0).sqrt();
            assert!(
                dij >= c.lower - 0.1 && dij <= c.upper + 0.1,
                "Distance ({},{}) = {:.4}, expected in [{:.2}, {:.2}]",
                c.i, c.j, dij, c.lower, c.upper
            );
        }
    }

    #[test]
    fn test_output_is_psd() {
        let constraints = make_simple_constraints();
        let initial = DMatrix::identity(3, 3) * 2.0;
        let config = RiemannianConfig::default();
        let optimizer = RiemannianLbfgs::new(config);
        let result = optimizer.optimize(&initial, &constraints);

        // Verify PSD: all eigenvalues ≥ 0
        let eigen = nalgebra::SymmetricEigen::new(result.metric_matrix.clone());
        for i in 0..eigen.eigenvalues.len() {
            assert!(
                eigen.eigenvalues[i] >= -1e-10,
                "Output has negative eigenvalue: {}",
                eigen.eigenvalues[i]
            );
        }
    }

    #[test]
    fn test_larger_system() {
        // 5-atom system with bond constraints
        let constraints = vec![
            DistanceConstraint { i: 0, j: 1, lower: 1.0, upper: 1.6 },
            DistanceConstraint { i: 1, j: 2, lower: 1.0, upper: 1.6 },
            DistanceConstraint { i: 2, j: 3, lower: 1.0, upper: 1.6 },
            DistanceConstraint { i: 3, j: 4, lower: 1.0, upper: 1.6 },
            DistanceConstraint { i: 0, j: 4, lower: 2.0, upper: 4.0 },
        ];
        let initial = DMatrix::identity(5, 5) * 3.0;
        let config = RiemannianConfig {
            max_iter: 500,
            ..Default::default()
        };
        let optimizer = RiemannianLbfgs::new(config);
        let result = optimizer.optimize(&initial, &constraints);

        assert!(
            result.objective < 1.0,
            "Objective too high for 5-atom system: {}",
            result.objective
        );

        // Output should be PSD
        let eigen = nalgebra::SymmetricEigen::new(result.metric_matrix.clone());
        for i in 0..eigen.eigenvalues.len() {
            assert!(eigen.eigenvalues[i] >= -1e-10);
        }
    }

    #[test]
    fn test_gradient_finite_difference() {
        let constraints = make_simple_constraints();
        let x = DMatrix::identity(3, 3) * 2.5;
        let config = RiemannianConfig::default();
        let optimizer = RiemannianLbfgs::new(config);

        let grad = optimizer.gradient(&x, &constraints);
        let h = 1e-5;
        let n = 3;

        for i in 0..n {
            for j in i..n {
                let mut xp = x.clone();
                let mut xm = x.clone();
                xp[(i, j)] += h;
                xp[(j, i)] = xp[(i, j)]; // Maintain symmetry
                xm[(i, j)] -= h;
                xm[(j, i)] = xm[(i, j)];

                let fp = optimizer.objective(&xp, &constraints);
                let fm = optimizer.objective(&xm, &constraints);
                let fd = (fp - fm) / (2.0 * h);

                assert!(
                    (grad[(i, j)] - fd).abs() < 1e-4,
                    "Gradient error at ({},{}): analytical={:.6}, fd={:.6}",
                    i, j, grad[(i, j)], fd
                );
            }
        }
    }
}
