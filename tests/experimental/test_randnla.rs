//! Integration tests for Track E2: RandNLA for EHT
//!
//! Validates Nyström approximation, randomized solver, and accuracy
//! compared to exact diagonalization.

#[cfg(feature = "experimental-randnla")]
mod randnla_tests {
    use sci_form::beta::rand_nla::*;

    use nalgebra::{DMatrix, DVector, SymmetricEigen};
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    // ---------------------------------------------------------------
    // Helpers
    // ---------------------------------------------------------------

    fn make_overlap_like(n: usize, seed: u64) -> DMatrix<f64> {
        // Build an overlap-like SPD matrix with diagonal dominance + spectral decay.
        // Real overlap matrices are close to I with off-diagonal decay,
        // so we simulate that: S = I + 0.3 * random symmetric with decaying off-diagonals.
        let mut rng = StdRng::seed_from_u64(seed);
        let mut s = DMatrix::identity(n, n);
        for i in 0..n {
            for j in (i + 1)..n {
                let decay = (-0.1 * (j - i) as f64).exp();
                let off = rng.gen_range(-0.3..0.3) * decay;
                s[(i, j)] = off;
                s[(j, i)] = off;
            }
        }
        // Ensure SPD by adding I (already have I on diagonal)
        s
    }

    fn make_hamiltonian_like(n: usize, seed: u64) -> DMatrix<f64> {
        let mut rng = StdRng::seed_from_u64(seed);
        let mut h = DMatrix::zeros(n, n);
        for i in 0..n {
            h[(i, i)] = rng.gen_range(-20.0..-5.0);
            for j in (i + 1)..n {
                let off = rng.gen_range(-3.0..0.0);
                h[(i, j)] = off;
                h[(j, i)] = off;
            }
        }
        h
    }

    fn solve_exact(h: &DMatrix<f64>, s: &DMatrix<f64>) -> (DVector<f64>, DMatrix<f64>) {
        let n = h.nrows();
        let s_eigen = SymmetricEigen::new(s.clone());
        let mut s_inv_sqrt_diag = DMatrix::zeros(n, n);
        for i in 0..n {
            let val = s_eigen.eigenvalues[i];
            if val > 1e-10 {
                s_inv_sqrt_diag[(i, i)] = 1.0 / val.sqrt();
            }
        }
        let s_inv_sqrt =
            &s_eigen.eigenvectors * &s_inv_sqrt_diag * s_eigen.eigenvectors.transpose();
        let h_prime = &s_inv_sqrt * h * &s_inv_sqrt;
        let h_eigen = SymmetricEigen::new(h_prime);
        let c = &s_inv_sqrt * &h_eigen.eigenvectors;

        let mut indices: Vec<usize> = (0..n).collect();
        let e = &h_eigen.eigenvalues;
        indices.sort_by(|&a, &b| e[a].partial_cmp(&e[b]).unwrap());

        let mut sorted_e = DVector::zeros(n);
        let mut sorted_c = DMatrix::zeros(n, n);
        for (new_i, &old_i) in indices.iter().enumerate() {
            sorted_e[new_i] = e[old_i];
            for row in 0..n {
                sorted_c[(row, new_i)] = c[(row, old_i)];
            }
        }
        (sorted_e, sorted_c)
    }

    // ---------------------------------------------------------------
    // E2.1 — Nyström Approximation
    // ---------------------------------------------------------------

    #[test]
    fn e2_1a_sketch_orthogonality() {
        let mut rng = StdRng::seed_from_u64(100);
        let n = 200;
        let k = 20;
        let sketch = GaussianSketch::new(&mut rng, n, k);
        let gram = sketch.omega.transpose() * &sketch.omega;
        let tol = 3.0 / (n as f64).sqrt();
        for i in 0..k {
            assert!(
                (gram[(i, i)] - 1.0).abs() < tol,
                "E[Ω^T Ω] diagonal ({},{}) = {:.4}, expected ~1 ± {:.4}",
                i, i, gram[(i, i)], tol
            );
        }
    }

    #[test]
    fn e2_1b_qr_projection() {
        let n = 30;
        let k = 15;
        let s = make_overlap_like(n, 200);
        let mut rng = StdRng::seed_from_u64(201);
        let sketch = GaussianSketch::new(&mut rng, n, k);
        let nys = NystromApprox::from_matrix(&s, &sketch);

        // Verify ||Y - Q(Q^T Y)||_F is small
        let y = &s * &sketch.omega;
        let q = &nys.q;
        let proj = q * (q.transpose() * &y);
        let diff = &y - &proj;
        let mut err = 0.0;
        for j in 0..diff.ncols() {
            for i in 0..diff.nrows() {
                err += diff[(i, j)] * diff[(i, j)];
            }
        }
        err = err.sqrt();
        assert!(
            err < 1e-10,
            "||Y - Q(Q^T Y)||_F = {:.2e}, expected < 1e-10",
            err
        );
    }

    #[test]
    fn e2_1c_reconstruction_error() {
        // Nyström reconstruction Q B Q^T captures the dominant subspace.
        // For a near-identity matrix, the best rank-k error is ~sqrt((n-k)/n).
        // This test verifies the approximation is not catastrophically wrong.
        let n = 30;
        let k = 28; // Near-full rank
        let s = make_overlap_like(n, 300);
        let mut rng = StdRng::seed_from_u64(301);
        let sketch = GaussianSketch::new(&mut rng, n, k);
        let nys = NystromApprox::from_matrix(&s, &sketch);
        let rel_err = nys.relative_error(&s);
        assert!(
            rel_err < 0.50,
            "Relative Frobenius error = {:.2}% with k={}, expected < 50%",
            rel_err * 100.0, k
        );
    }

    // ---------------------------------------------------------------
    // E2.2 — Randomized Solver
    // ---------------------------------------------------------------

    #[test]
    fn e2_2a_inverse_sqrt_quality() {
        let n = 20;
        let k = 20; // Full rank for this test
        let s = make_overlap_like(n, 400);
        let mut rng = StdRng::seed_from_u64(401);
        let sketch = GaussianSketch::new(&mut rng, n, k);
        let nys = NystromApprox::from_matrix(&s, &sketch);
        let s_inv_sqrt = nys.inverse_sqrt();

        // S^{-1/2} S S^{-1/2} ≈ I
        let product = &s_inv_sqrt * &s * &s_inv_sqrt;
        for i in 0..n {
            assert!(
                (product[(i, i)] - 1.0).abs() < 0.15,
                "S^(-1/2) S S^(-1/2) diagonal ({},{}) = {:.4}, expected ~1",
                i, i, product[(i, i)]
            );
        }
    }

    #[test]
    fn e2_2b_randnla_eigenvalues_vs_exact() {
        let n = 20;
        let h = make_hamiltonian_like(n, 500);
        let s = make_overlap_like(n, 501);

        let (e_exact, _) = solve_exact(&h, &s);

        // Use k = n (full) to test correctness of the solver path
        let config = RandNlaConfig {
            sketch_size: Some(n),
            seed: 502,
            fallback_enabled: false,
            max_error: 1.0,
        };
        let (e_rand, _, _info) = solve_eht_randnla(&h, &s, &config);

        // With k = n, HOMO/LUMO should match closely
        let homo_idx = n / 2 - 1;
        let homo_exact = e_exact[homo_idx];
        let homo_rand = e_rand[homo_idx];
        let gap_exact = e_exact[homo_idx + 1] - homo_exact;
        let gap_rand = e_rand[homo_idx + 1] - homo_rand;

        assert!(
            ((homo_exact - homo_rand) / homo_exact).abs() < 0.01,
            "HOMO error: exact={:.6}, rand={:.6}",
            homo_exact, homo_rand
        );
        assert!(
            (gap_exact - gap_rand).abs() < 0.5,
            "Gap error: exact={:.6}, rand={:.6}",
            gap_exact, gap_rand
        );
    }

    #[test]
    fn e2_2c_error_bound_and_fallback() {
        let n = 15;
        let h = make_hamiltonian_like(n, 600);
        let s = make_overlap_like(n, 601);

        // Tight threshold with small k → fallback expected
        let config = RandNlaConfig {
            sketch_size: Some(3),
            seed: 602,
            max_error: 1e-10,
            fallback_enabled: true,
        };
        let (_, _, info) = solve_eht_randnla(&h, &s, &config);
        assert!(
            info.used_fallback,
            "Expected fallback with k=3 and max_error=1e-10"
        );
    }

    // ---------------------------------------------------------------
    // E2.3 — Accuracy validation
    // ---------------------------------------------------------------

    #[test]
    fn e2_3a_homo_lumo_accuracy_large() {
        // 40-orbital system — randomized with k=40 (full)
        let n = 40;
        let h = make_hamiltonian_like(n, 700);
        let s = make_overlap_like(n, 701);

        let (e_exact, _) = solve_exact(&h, &s);

        let config = RandNlaConfig {
            sketch_size: Some(n),
            seed: 702,
            fallback_enabled: false,
            max_error: 1.0,
        };
        let (e_rand, _, info) = solve_eht_randnla(&h, &s, &config);

        let homo = n / 2 - 1;
        let homo_err = ((e_exact[homo] - e_rand[homo]) / e_exact[homo]).abs();
        let lumo_err = ((e_exact[homo + 1] - e_rand[homo + 1]) / e_exact[homo + 1]).abs();

        assert!(
            homo_err < 0.001,
            "HOMO relative error = {:.6}, expected < 0.1%",
            homo_err
        );
        assert!(
            lumo_err < 0.001,
            "LUMO relative error = {:.6}, expected < 0.1%",
            lumo_err
        );
        assert!(
            info.residual_error < 0.01,
            "Residual = {:.6}, expected < 1%",
            info.residual_error
        );
    }

    #[test]
    fn e2_3a_gap_accuracy() {
        let n = 30;
        let h = make_hamiltonian_like(n, 800);
        let s = make_overlap_like(n, 801);

        let (e_exact, _) = solve_exact(&h, &s);
        let config = RandNlaConfig {
            sketch_size: Some(n),
            seed: 802,
            fallback_enabled: false,
            max_error: 1.0,
        };
        let (e_rand, _, _) = solve_eht_randnla(&h, &s, &config);

        let homo = n / 2 - 1;
        let gap_exact = e_exact[homo + 1] - e_exact[homo];
        let gap_rand = e_rand[homo + 1] - e_rand[homo];

        assert!(
            (gap_exact - gap_rand).abs() < 0.01,
            "Gap error: exact={:.6}, rand={:.6}, diff={:.6}",
            gap_exact, gap_rand, (gap_exact - gap_rand).abs()
        );
    }
}
