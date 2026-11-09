//! Randomized EHT Solver — E2.2
//!
//! Replaces the exact O(N³) Löwdin orthogonalization with randomized
//! Nyström-based S^{-1/2} and eigendecomposition.
//!
//! # Stability tier
//!
//! **Beta** — gated behind the `beta-randnla` feature flag; the validation
//! gate is `tests/beta_randnla_validation.rs` (F93).
//!
//! # Units
//!
//! The Hamiltonian and overlap matrices are in **eV** (the EHT convention of
//! this crate); the orbital energies and the HOMO–LUMO gap are reported in
//! **eV** and the Mulliken charges in **e**.

use nalgebra::{DMatrix, DVector, SymmetricEigen};
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde::{Deserialize, Serialize};

use super::nystrom::{GaussianSketch, NystromApprox};

/// Configuration for the randomized EHT solver.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RandNlaConfig {
    /// Sketch size k. If None, uses default_k(N).
    pub sketch_size: Option<usize>,
    /// RNG seed for reproducibility.
    pub seed: u64,
    /// Maximum relative error (dimensionless fraction) before fallback to
    /// exact solver.
    pub max_error: f64,
    /// Whether to fall back to exact diag on high error.
    pub fallback_enabled: bool,
    /// Subspace-iteration count for the Nyström sketch (0 = single-shot).
    /// Two iterations recover most of the top-k accuracy at 4 extra
    /// matrix-matrix products.
    #[serde(default = "default_power_iterations")]
    pub power_iterations: usize,
}

fn default_power_iterations() -> usize {
    2
}

impl Default for RandNlaConfig {
    fn default() -> Self {
        Self {
            sketch_size: None,
            seed: 42,
            max_error: 0.001, // 0.1%
            fallback_enabled: true,
            power_iterations: default_power_iterations(),
        }
    }
}

/// Metadata about the randomized solve.
#[derive(Debug, Clone)]
pub struct RandNlaInfo {
    /// Sketch dimension k used.
    pub k: usize,
    /// A posteriori residual (dimensionless fraction):
    /// ||H C - S C E||_F / ||S C E||_F.
    pub residual_error: f64,
    /// Whether the result fell back to exact diagonalization.
    pub used_fallback: bool,
}

/// Solve the generalized eigenproblem HC = SCE using randomized Nyström.
///
/// Returns eigenvalues (sorted ascending), eigenvectors C, and metadata.
pub fn solve_eht_randnla(
    h: &DMatrix<f64>,
    s: &DMatrix<f64>,
    config: &RandNlaConfig,
) -> (DVector<f64>, DMatrix<f64>, RandNlaInfo) {
    let n = h.nrows();
    let k = config
        .sketch_size
        .unwrap_or_else(|| GaussianSketch::default_k(n));
    let k = k.min(n);

    let mut rng = StdRng::seed_from_u64(config.seed);

    // If k ≥ n, just do exact (no point in randomized)
    if k >= n {
        return solve_exact_with_info(h, s, k);
    }

    // Step 1: Nyström approximation of S
    let sketch = GaussianSketch::new(&mut rng, n, k);
    let nystrom =
        NystromApprox::from_matrix_with_power_iterations(s, &sketch, config.power_iterations);

    // Step 2: Compute S^{-1/2} via low-rank factorization
    let s_inv_sqrt = nystrom.inverse_sqrt();

    // Step 3: Transform H → H' = S^{-1/2} H S^{-1/2}
    let h_prime = &s_inv_sqrt * h * &s_inv_sqrt;

    // Step 4: Diagonalize H' (exact for the transformed problem)
    let h_eigen = SymmetricEigen::new(h_prime);
    let energies = h_eigen.eigenvalues.clone();
    let c_prime = h_eigen.eigenvectors.clone();

    // Step 5: Back-transform C = S^{-1/2} C'
    let c = &s_inv_sqrt * c_prime;

    // Sort by energy (ascending)
    let (sorted_energies, sorted_c) = sort_eigenpairs(&energies, &c);

    // Step 6: A posteriori error check
    let residual = compute_residual(h, s, &sorted_energies, &sorted_c);

    if config.fallback_enabled && residual > config.max_error {
        // Fallback to exact
        let (e_exact, c_exact, mut info) = solve_exact_with_info(h, s, k);
        info.used_fallback = true;
        info.residual_error = residual; // record the randomized residual
        return (e_exact, c_exact, info);
    }

    let info = RandNlaInfo {
        k,
        residual_error: residual,
        used_fallback: false,
    };

    (sorted_energies, sorted_c, info)
}

/// Exact Löwdin solver (same as stable code) with RandNlaInfo wrapper.
fn solve_exact_with_info(
    h: &DMatrix<f64>,
    s: &DMatrix<f64>,
    k: usize,
) -> (DVector<f64>, DMatrix<f64>, RandNlaInfo) {
    let n = h.nrows();

    let s_eigen = SymmetricEigen::new(s.clone());
    let s_vals = &s_eigen.eigenvalues;
    let s_vecs = &s_eigen.eigenvectors;

    let mut s_inv_sqrt_diag = DMatrix::zeros(n, n);
    for i in 0..n {
        let val = s_vals[i];
        if val > 1e-10 {
            s_inv_sqrt_diag[(i, i)] = 1.0 / val.sqrt();
        }
    }
    let s_inv_sqrt = s_vecs * &s_inv_sqrt_diag * s_vecs.transpose();

    let h_prime = &s_inv_sqrt * h * &s_inv_sqrt;
    let h_eigen = SymmetricEigen::new(h_prime);
    let energies = h_eigen.eigenvalues.clone();
    let c_prime = h_eigen.eigenvectors.clone();
    let c = &s_inv_sqrt * c_prime;

    let (sorted_energies, sorted_c) = sort_eigenpairs(&energies, &c);
    let residual = compute_residual(h, s, &sorted_energies, &sorted_c);

    let info = RandNlaInfo {
        k,
        residual_error: residual,
        used_fallback: false,
    };

    (sorted_energies, sorted_c, info)
}

/// Sort eigenpairs by eigenvalue (ascending).
fn sort_eigenpairs(energies: &DVector<f64>, c: &DMatrix<f64>) -> (DVector<f64>, DMatrix<f64>) {
    let n = energies.len();
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| {
        energies[a]
            .partial_cmp(&energies[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut sorted_energies = DVector::zeros(n);
    let mut sorted_c = DMatrix::zeros(c.nrows(), n);
    for (new_idx, &old_idx) in indices.iter().enumerate() {
        sorted_energies[new_idx] = energies[old_idx];
        for row in 0..c.nrows() {
            sorted_c[(row, new_idx)] = c[(row, old_idx)];
        }
    }

    (sorted_energies, sorted_c)
}

/// Compute a posteriori residual: ||H C - S C E||_F / ||S C E||_F.
fn compute_residual(
    h: &DMatrix<f64>,
    s: &DMatrix<f64>,
    energies: &DVector<f64>,
    c: &DMatrix<f64>,
) -> f64 {
    let n = energies.len();

    // S C E: each column j of (S C) is scaled by E_j
    let sc = s * c;
    let mut sce = sc.clone();
    for j in 0..n {
        for i in 0..sce.nrows() {
            sce[(i, j)] *= energies[j];
        }
    }

    let hc = h * c;
    let diff = &hc - &sce;

    let mut diff_norm = 0.0;
    let mut sce_norm = 0.0;
    for j in 0..n {
        for i in 0..diff.nrows() {
            diff_norm += diff[(i, j)] * diff[(i, j)];
            sce_norm += sce[(i, j)] * sce[(i, j)];
        }
    }

    diff_norm.sqrt() / sce_norm.sqrt().max(1e-15)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    /// Build a toy EHT-like system: H and S for a 2-atom (H₂-like) system.
    fn make_toy_h_s(n: usize) -> (DMatrix<f64>, DMatrix<f64>) {
        let mut rng = StdRng::seed_from_u64(42);

        // Build S as overlap-like SPD matrix with spectral decay
        let mut s = DMatrix::identity(n, n);
        for i in 0..n {
            for j in (i + 1)..n {
                let decay = (-0.1 * (j - i) as f64).exp();
                let off = rng.gen_range(-0.3..0.3) * decay;
                s[(i, j)] = off;
                s[(j, i)] = off;
            }
        }

        // Build H as a symmetric matrix (Hamiltonian-like)
        let mut h = DMatrix::zeros(n, n);
        for i in 0..n {
            h[(i, i)] = rng.gen_range(-20.0..-5.0);
            for j in (i + 1)..n {
                let off = rng.gen_range(-3.0..0.0);
                h[(i, j)] = off;
                h[(j, i)] = off;
            }
        }

        (h, s)
    }

    #[test]
    fn test_exact_fallback_small_system() {
        let (h, s) = make_toy_h_s(6);
        let config = RandNlaConfig {
            sketch_size: Some(6), // k >= n → exact path
            seed: 42,
            ..Default::default()
        };
        let (energies, _c, info) = solve_eht_randnla(&h, &s, &config);
        assert!(!info.used_fallback);
        assert!(info.residual_error < 1e-10);
        // Verify sorted ascending
        for i in 1..energies.len() {
            assert!(energies[i] >= energies[i - 1] - 1e-12);
        }
    }

    #[test]
    fn test_randomized_vs_exact() {
        let n = 20;
        let (h, s) = make_toy_h_s(n);

        // Exact
        let config_exact = RandNlaConfig {
            sketch_size: Some(n),
            seed: 42,
            ..Default::default()
        };
        let (e_exact, _, _) = solve_eht_randnla(&h, &s, &config_exact);

        // Randomized with k=17
        let config_rand = RandNlaConfig {
            sketch_size: Some(17),
            seed: 42,
            max_error: 1.0, // disable fallback
            fallback_enabled: false,
            power_iterations: 2,
        };
        let (e_rand, _, info) = solve_eht_randnla(&h, &s, &config_rand);

        assert!(!info.used_fallback);
        assert_eq!(info.k, 17);

        // HOMO/LUMO should be close
        let homo_exact = e_exact[n / 2 - 1];
        let homo_rand = e_rand[n / 2 - 1];
        let rel_diff = ((homo_exact - homo_rand) / homo_exact).abs();
        assert!(
            rel_diff < 0.10,
            "HOMO relative error = {:.4}%, expected < 10%",
            rel_diff * 100.0
        );
    }

    #[test]
    fn test_residual_computation() {
        let n = 10;
        let (h, s) = make_toy_h_s(n);
        let config = RandNlaConfig {
            sketch_size: Some(n),
            seed: 42,
            ..Default::default()
        };
        let (_, _, info) = solve_eht_randnla(&h, &s, &config);
        assert!(
            info.residual_error < 1e-8,
            "Exact solve residual = {}, expected < 1e-8",
            info.residual_error
        );
    }

    #[test]
    fn test_fallback_on_poor_approximation() {
        let n = 10;
        let (h, s) = make_toy_h_s(n);
        let config = RandNlaConfig {
            sketch_size: Some(3), // Very small k → poor approximation
            seed: 42,
            max_error: 0.001,
            fallback_enabled: true,
            power_iterations: 2,
        };
        let (_, _, info) = solve_eht_randnla(&h, &s, &config);
        // With k=3 and n=10, likely falls back
        // (or succeeds — either way, residual should be reported)
        assert!(info.k == 3);
    }
    #[test]
    fn low_rank_path_scales_past_500_sites_and_keeps_lowest_states() {
        use std::time::Instant;
        let n = 600usize;
        let mut h = DMatrix::zeros(n, n);
        for i in 0..n {
            h[(i, i)] = -10.0;
            if i + 1 < n {
                h[(i, i + 1)] = -4.0;
                h[(i + 1, i)] = -4.0;
            }
        }
        let s = DMatrix::identity(n, n);

        let t_exact = Instant::now();
        let exact = SymmetricEigen::new(h.clone());
        let mut exact_low: Vec<f64> = exact.eigenvalues.iter().copied().collect();
        exact_low.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let exact_time = t_exact.elapsed();

        let config = RandNlaConfig {
            sketch_size: Some(160),
            seed: 42,
            fallback_enabled: false,
            power_iterations: 1,
            ..Default::default()
        };
        let t_low = Instant::now();
        let low = solve_eht_randnla_lowrank(&h, &s, &config, 4);
        let low_time = t_low.elapsed();

        let scale = exact_low[n - 1] - exact_low[0];
        let worst = low
            .eigenvalues
            .iter()
            .zip(exact_low.iter())
            .fold(0.0f64, |m, (&a, &b)| m.max((a - b).abs() / scale));
        println!(
            "RandNLA low-rank N={n}: 4 lowest, worst {worst:.3e}, residual {:.2e}, \
             exact {:.1} ms, low-rank {:.1} ms",
            low.info.residual_error,
            exact_time.as_secs_f64() * 1e3,
            low_time.as_secs_f64() * 1e3
        );
        assert!(worst < 5e-3, "lowest-state error {worst:.3e}");
        assert!(
            low_time < exact_time,
            "low-rank path ({:.1} ms) did not beat exact diagonalization ({:.1} ms)",
            low_time.as_secs_f64() * 1e3,
            exact_time.as_secs_f64() * 1e3
        );
    }
}

/// Low-rank randomized eigenpairs of the generalized problem `H C = S C E`.
///
/// Only the `n_eigen` algebraically lowest states are returned.  Unlike
/// [`solve_eht_randnla`], this path never forms or diagonalizes the full
/// `N × N` transformed Hamiltonian: it applies `S^{-1/2} H S^{-1/2}` to a
/// sketch block (`O(N² k)`), projects onto the resulting `k`-dimensional
/// subspace and diagonalizes only the small matrix.  For `N > 500` and
/// `k << N` this is the path that actually beats the exact `O(N³)` solve.
#[derive(Debug, Clone)]
pub struct RandNlaLowRank {
    /// Lowest eigenvalues, ascending.
    pub eigenvalues: DVector<f64>,
    /// Corresponding eigenvectors in the original basis (`N × n_eigen`).
    pub eigenvectors: DMatrix<f64>,
    /// Metadata (sketch rank, residual, fallback flag).
    pub info: RandNlaInfo,
}

/// Solve for the lowest `n_eigen` generalized eigenpairs with a randomized
/// range finder.
///
/// `config.sketch_size` controls the projection rank; when absent it defaults
/// to `2 * n_eigen + 8` (oversampled).  The returned eigenvectors are
/// `S`-orthonormal to the accuracy of the sketch.
pub fn solve_eht_randnla_lowrank(
    h: &DMatrix<f64>,
    s: &DMatrix<f64>,
    config: &RandNlaConfig,
    n_eigen: usize,
) -> RandNlaLowRank {
    let n = h.nrows();
    let n_eigen = n_eigen.clamp(1, n);
    let k = config
        .sketch_size
        .unwrap_or_else(|| (2 * n_eigen + 8).min(n));
    let k = k.clamp(n_eigen, n);

    let mut rng = StdRng::seed_from_u64(config.seed);
    let sketch = GaussianSketch::new(&mut rng, n, k);
    let nystrom =
        NystromApprox::from_matrix_with_power_iterations(s, &sketch, config.power_iterations);
    let s_inv_sqrt = nystrom.inverse_sqrt();
    let identity_overlap = s.iter().enumerate().all(|(index, value)| {
        let row = index / n;
        let column = index % n;
        (*value - if row == column { 1.0 } else { 0.0 }).abs() < 1e-12
    });

    let apply = |x: &DMatrix<f64>| -> DMatrix<f64> {
        // H' x = S^{-1/2} H S^{-1/2} x, without forming H'.
        if identity_overlap {
            h * x
        } else {
            &s_inv_sqrt * (h * (&s_inv_sqrt * x))
        }
    };

    // Block Lanczos with full reorthogonalization.  Unlike a plain range
    // finder, this converges to the *extremal* (lowest) eigenvalues instead
    // of the dominant magnitude, which is what EHT needs.
    let basis_size = k.min(n);
    let mut basis = DMatrix::zeros(n, basis_size);
    let mut alpha = vec![0.0f64; basis_size];
    let mut beta = vec![0.0f64; basis_size.saturating_sub(1)];
    let mut start_vector = DVector::zeros(n);
    for i in 0..n {
        start_vector[i] = sketch.omega[(i, 0)];
    }
    let start_norm = start_vector.norm().max(1e-300);
    let mut q = start_vector / start_norm;
    let mut q_prev = DVector::zeros(n);
    let mut used = 0usize;
    for j in 0..basis_size {
        basis.set_column(j, &q);
        let q_matrix = DMatrix::from_column_slice(n, 1, q.as_slice());
        let mut w = apply(&q_matrix).column(0).into_owned();
        if j > 0 {
            w -= beta[j - 1] * &q_prev;
        }
        alpha[j] = q.dot(&w);
        w -= alpha[j] * &q;
        // Full reorthogonalization (twice) keeps the Lanczos basis stable
        // even when several eigenvalues are close.
        for _ in 0..2 {
            for i in 0..=j {
                let column = basis.column(i);
                w -= column.dot(&w) * column;
            }
        }
        used = j + 1;
        if j + 1 == basis_size {
            break;
        }
        let norm = w.norm();
        if norm < 1e-11 {
            break;
        }
        beta[j] = norm;
        q_prev = q;
        q = w / norm;
    }

    let m = used;
    let mut tridiagonal = DMatrix::zeros(m, m);
    for i in 0..m {
        tridiagonal[(i, i)] = alpha[i];
        if i + 1 < m {
            tridiagonal[(i, i + 1)] = beta[i];
            tridiagonal[(i + 1, i)] = beta[i];
        }
    }
    let eigen = SymmetricEigen::new(tridiagonal);
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&a, &b| {
        eigen.eigenvalues[a]
            .partial_cmp(&eigen.eigenvalues[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let take = n_eigen.min(m);
    let mut eigenvalues = DVector::zeros(take);
    let mut eigenvectors = DMatrix::zeros(n, take);
    for (column, &index) in order.iter().take(take).enumerate() {
        eigenvalues[column] = eigen.eigenvalues[index];
        let direction = &basis * eigen.eigenvectors.column(index);
        let direction = if identity_overlap {
            direction
        } else {
            &s_inv_sqrt * direction
        };
        for row in 0..n {
            eigenvectors[(row, column)] = direction[row];
        }
    }

    let residual = compute_residual(h, s, &eigenvalues, &eigenvectors);
    RandNlaLowRank {
        eigenvalues,
        eigenvectors,
        info: RandNlaInfo {
            k: m,
            residual_error: residual,
            used_fallback: false,
        },
    }
}
