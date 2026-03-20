//! Integration tests for Track E3: Riemannian Optimization for ETKDG
//!
//! Validates PSD manifold operations, Riemannian L-BFGS convergence,
//! and geometric quality of optimized metric matrices.

#[cfg(feature = "experimental-riemannian")]
mod riemannian_tests {
    use sci_form::experimental::riemannian::*;
    use nalgebra::{DMatrix, SymmetricEigen};

    // ---------------------------------------------------------------
    // E3.1 — Riemannian Geometry Infrastructure
    // ---------------------------------------------------------------

    #[test]
    fn e3_1a_psd_projection_clamps_negative() {
        // Matrix with known negative eigenvalue
        let mut m = DMatrix::identity(5, 5);
        m[(0, 0)] = -2.0;
        m[(1, 1)] = -1.0;
        let p = psd_projection(&m);

        let eigen = SymmetricEigen::new(p);
        for i in 0..eigen.eigenvalues.len() {
            assert!(
                eigen.eigenvalues[i] >= -1e-10,
                "Eigenvalue {} = {}, expected ≥ 0",
                i, eigen.eigenvalues[i]
            );
        }
    }

    #[test]
    fn e3_1a_psd_projection_idempotent() {
        let m = DMatrix::identity(6, 6) * 3.0;
        let p1 = psd_projection(&m);
        let p2 = psd_projection(&p1);
        let n = p1.nrows();
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (p1[(i, j)] - p2[(i, j)]).abs() < 1e-10,
                    "PSD projection not idempotent at ({},{})",
                    i, j
                );
            }
        }
    }

    #[test]
    fn e3_1b_riemannian_gradient_consistency() {
        // The Riemannian gradient should be in the tangent space (symmetric)
        let constraints = vec![
            DistanceConstraint { i: 0, j: 1, lower: 1.0, upper: 2.0 },
            DistanceConstraint { i: 1, j: 2, lower: 1.0, upper: 2.0 },
        ];
        let x = DMatrix::identity(3, 3) * 2.0;
        let config = RiemannianConfig::default();
        let optimizer = RiemannianLbfgs::new(config);
        let grad = optimizer.gradient(&x, &constraints);

        let n = 3;
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (grad[(i, j)] - grad[(j, i)]).abs() < 1e-14,
                    "Gradient not symmetric at ({},{}): {} vs {}",
                    i, j, grad[(i, j)], grad[(j, i)]
                );
            }
        }
    }

    #[test]
    fn e3_1b_gradient_matches_finite_difference() {
        let constraints = vec![
            DistanceConstraint { i: 0, j: 1, lower: 1.3, upper: 1.7 },
            DistanceConstraint { i: 0, j: 2, lower: 2.0, upper: 2.5 },
            DistanceConstraint { i: 1, j: 2, lower: 1.3, upper: 1.7 },
        ];
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
                xp[(j, i)] = xp[(i, j)];
                xm[(i, j)] -= h;
                xm[(j, i)] = xm[(i, j)];
                let fp = optimizer.objective(&xp, &constraints);
                let fm = optimizer.objective(&xm, &constraints);
                let fd = (fp - fm) / (2.0 * h);
                assert!(
                    (grad[(i, j)] - fd).abs() < 1e-3,
                    "Gradient error at ({},{}): analytical={:.6}, fd={:.6}",
                    i, j, grad[(i, j)], fd
                );
            }
        }
    }

    #[test]
    fn e3_1c_retraction_stays_on_manifold() {
        let x = DMatrix::identity(6, 6) * 4.0;
        let mut xi = DMatrix::zeros(6, 6);
        xi[(0, 1)] = -0.5;
        xi[(1, 0)] = -0.5;
        xi[(2, 3)] = 0.3;
        xi[(3, 2)] = 0.3;

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

    // ---------------------------------------------------------------
    // E3.2 — Riemannian L-BFGS
    // ---------------------------------------------------------------

    #[test]
    fn e3_2_triangle_convergence() {
        // Equilateral triangle with sides ~1.5 Å
        let constraints = vec![
            DistanceConstraint { i: 0, j: 1, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 0, j: 2, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 1, j: 2, lower: 1.4, upper: 1.6 },
        ];
        let initial = DMatrix::identity(3, 3) * 2.0;
        let config = RiemannianConfig {
            max_iter: 300,
            ..Default::default()
        };
        let optimizer = RiemannianLbfgs::new(config);
        let result = optimizer.optimize(&initial, &constraints);

        assert!(
            result.converged || result.objective < 1e-4,
            "Triangle should converge: obj={:.6}, grad={:.6}, iter={}",
            result.objective, result.grad_norm, result.iterations
        );
    }

    #[test]
    fn e3_2_output_always_psd() {
        // Various initial conditions, output must always be PSD
        for scale in [0.5, 1.0, 2.0, 5.0, 10.0] {
            let constraints = vec![
                DistanceConstraint { i: 0, j: 1, lower: 1.0, upper: 2.0 },
                DistanceConstraint { i: 1, j: 2, lower: 1.0, upper: 2.0 },
                DistanceConstraint { i: 2, j: 3, lower: 1.0, upper: 2.0 },
            ];
            let initial = DMatrix::identity(4, 4) * scale;
            let config = RiemannianConfig {
                max_iter: 100,
                ..Default::default()
            };
            let optimizer = RiemannianLbfgs::new(config);
            let result = optimizer.optimize(&initial, &constraints);

            let eigen = SymmetricEigen::new(result.metric_matrix);
            for i in 0..eigen.eigenvalues.len() {
                assert!(
                    eigen.eigenvalues[i] >= -1e-10,
                    "Output not PSD at scale {}: eigenvalue {} = {}",
                    scale, i, eigen.eigenvalues[i]
                );
            }
        }
    }

    #[test]
    fn e3_2_chain_molecule() {
        // Linear chain of 6 atoms with bond + angle constraints
        let constraints = vec![
            // Bonds
            DistanceConstraint { i: 0, j: 1, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 1, j: 2, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 2, j: 3, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 3, j: 4, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 4, j: 5, lower: 1.4, upper: 1.6 },
            // 1-3 distances (angles)
            DistanceConstraint { i: 0, j: 2, lower: 2.2, upper: 2.6 },
            DistanceConstraint { i: 1, j: 3, lower: 2.2, upper: 2.6 },
            DistanceConstraint { i: 2, j: 4, lower: 2.2, upper: 2.6 },
            DistanceConstraint { i: 3, j: 5, lower: 2.2, upper: 2.6 },
        ];
        let initial = DMatrix::identity(6, 6) * 3.0;
        let config = RiemannianConfig {
            max_iter: 500,
            ..Default::default()
        };
        let optimizer = RiemannianLbfgs::new(config);
        let result = optimizer.optimize(&initial, &constraints);

        assert!(
            result.objective < 1.0,
            "Chain molecule objective too high: {}",
            result.objective
        );
        assert!(!result.used_fallback);
    }

    // ---------------------------------------------------------------
    // E3.3 — Validation
    // ---------------------------------------------------------------

    #[test]
    fn e3_3a_objective_decreases() {
        let constraints = vec![
            DistanceConstraint { i: 0, j: 1, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 0, j: 2, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 1, j: 2, lower: 1.4, upper: 1.6 },
        ];
        let initial = DMatrix::identity(3, 3) * 5.0;
        let config = RiemannianConfig {
            max_iter: 50,
            ..Default::default()
        };
        let optimizer = RiemannianLbfgs::new(config);

        let obj_before = optimizer.objective(&initial, &constraints);
        let result = optimizer.optimize(&initial, &constraints);

        assert!(
            result.objective <= obj_before + 1e-10,
            "Objective should decrease: before={:.6}, after={:.6}",
            obj_before, result.objective
        );
    }

    #[test]
    fn e3_3b_distance_violation_bounded() {
        let constraints = vec![
            DistanceConstraint { i: 0, j: 1, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 0, j: 2, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 1, j: 2, lower: 1.4, upper: 1.6 },
        ];
        let initial = DMatrix::identity(3, 3) * 2.0;
        let config = RiemannianConfig {
            max_iter: 500,
            ..Default::default()
        };
        let optimizer = RiemannianLbfgs::new(config);
        let result = optimizer.optimize(&initial, &constraints);

        // Max violation should be small
        assert!(
            result.max_violation < 0.5,
            "Max distance violation = {:.4} Å, expected < 0.5",
            result.max_violation
        );
    }

    #[test]
    fn e3_3c_manifold_point_rank() {
        // Verify that the solution metric matrix has appropriate rank
        let constraints = vec![
            DistanceConstraint { i: 0, j: 1, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 0, j: 2, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 1, j: 2, lower: 1.4, upper: 1.6 },
        ];
        let initial = DMatrix::identity(3, 3) * 2.0;
        let config = RiemannianConfig::default();
        let optimizer = RiemannianLbfgs::new(config);
        let result = optimizer.optimize(&initial, &constraints);

        let mut m = PsdManifold::from_psd(result.metric_matrix);
        let rank = m.rank();
        assert!(
            rank >= 2,
            "Metric matrix should have rank ≥ 2 for 3-atom triangle, got {}",
            rank
        );
    }
}
