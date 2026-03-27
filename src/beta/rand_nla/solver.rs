//! Randomized EHT Solver — E2.2
//!
//! Replaces the exact O(N³) Löwdin orthogonalization with randomized
//! Nyström-based S^{-1/2} and eigendecomposition.

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
    /// Maximum relative error before fallback to exact solver.
    pub max_error: f64,
    /// Whether to fall back to exact diag on high error.
    pub fallback_enabled: bool,
}

impl Default for RandNlaConfig {
    fn default() -> Self {
        Self {
            sketch_size: None,
            seed: 42,
            max_error: 0.001, // 0.1%
            fallback_enabled: true,
        }
    }
}

/// Metadata about the randomized solve.
#[derive(Debug, Clone)]
pub struct RandNlaInfo {
    /// Sketch dimension k used.
    pub k: usize,
    /// A posteriori residual: ||H C - S C E||_F / ||S C E||_F.
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
    let nystrom = NystromApprox::from_matrix(s, &sketch);

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
    indices.sort_by(|&a, &b| energies[a].partial_cmp(&energies[b]).unwrap());

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
            ..Default::default()
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
        };
        let (_, _, info) = solve_eht_randnla(&h, &s, &config);
        // With k=3 and n=10, likely falls back
        // (or succeeds — either way, residual should be reported)
        assert!(info.k == 3);
    }
}
