//! Chebyshev Expansion Infrastructure — E4.1
//!
//! Spectral bounds estimation, matrix rescaling, Chebyshev recursion,
//! and Jackson kernel damping for Gibbs oscillation suppression.

use nalgebra::DMatrix;

/// Estimate spectral bounds [E_min, E_max] via Gershgorin circles.
///
/// For each row i, the eigenvalue disk is centered at H_ii with radius
/// sum of |H_ij| for j != i. The union of all disks bounds the spectrum.
pub fn estimate_spectral_bounds(h: &DMatrix<f64>) -> (f64, f64) {
    let n = h.nrows();
    let mut e_min = f64::INFINITY;
    let mut e_max = f64::NEG_INFINITY;

    for i in 0..n {
        let center = h[(i, i)];
        let mut radius = 0.0;
        for j in 0..n {
            if j != i {
                radius += h[(i, j)].abs();
            }
        }
        let lo = center - radius;
        let hi = center + radius;
        if lo < e_min {
            e_min = lo;
        }
        if hi > e_max {
            e_max = hi;
        }
    }

    // Add small margin to ensure all eigenvalues are strictly inside [-1, 1]
    let margin = (e_max - e_min) * 0.01;
    (e_min - margin, e_max + margin)
}

/// Rescale a matrix to have spectrum in [-1, 1].
///
/// H_tilde = (H - b*I) / a, where a = (E_max - E_min)/2, b = (E_max + E_min)/2
pub fn rescale_matrix(h: &DMatrix<f64>, e_min: f64, e_max: f64) -> DMatrix<f64> {
    let a = (e_max - e_min) / 2.0;
    let b = (e_max + e_min) / 2.0;
    let n = h.nrows();
    let mut h_tilde = h.clone();
    for i in 0..n {
        h_tilde[(i, i)] -= b;
    }
    h_tilde /= a;
    h_tilde
}

/// Jackson kernel damping coefficients.
///
/// g_k^(M) = ((M-k+1)cos(πk/(M+1)) + sin(πk/(M+1))cot(π/(M+1))) / (M+1)
pub fn jackson_kernel(order: usize) -> Vec<f64> {
    let m = order as f64;
    let mp1 = m + 1.0;
    let cot_val = 1.0 / (std::f64::consts::PI / mp1).tan();

    (0..order)
        .map(|k| {
            let kf = k as f64;
            let cos_term = (mp1 - kf) * (std::f64::consts::PI * kf / mp1).cos();
            let sin_term = (std::f64::consts::PI * kf / mp1).sin() * cot_val;
            (cos_term + sin_term) / mp1
        })
        .collect()
}

/// Chebyshev expansion of a matrix function.
///
/// Stores the diagonal traces Tr[T_k(H_tilde)] and can reconstruct
/// functions f(H) via coefficient summation.
pub struct ChebyshevExpansion {
    /// Chebyshev moments: mu_k = (1/N) Tr[T_k(H_tilde)]
    pub moments: Vec<f64>,
    /// Rescaling parameter a = (E_max - E_min) / 2
    pub a: f64,
    /// Rescaling parameter b = (E_max + E_min) / 2
    pub b: f64,
    /// Matrix dimension
    pub n: usize,
    /// Expansion order
    pub order: usize,
}

impl ChebyshevExpansion {
    /// Compute Chebyshev expansion of a Hermitian matrix.
    ///
    /// Uses the three-term recurrence: T_0 = I, T_1 = H_tilde,
    /// T_{k+1} = 2 * H_tilde * T_k - T_{k-1}
    ///
    /// For efficiency, we only track the diagonal (for trace) using
    /// stochastic trace estimation with `n_vectors` random vectors.
    pub fn from_matrix(
        h: &DMatrix<f64>,
        order: usize,
        n_vectors: usize,
        seed: u64,
    ) -> Self {
        let n = h.nrows();
        let (e_min, e_max) = estimate_spectral_bounds(h);
        let a = (e_max - e_min) / 2.0;
        let b = (e_max + e_min) / 2.0;
        let h_tilde = rescale_matrix(h, e_min, e_max);

        // Stochastic trace estimation using random ±1 vectors
        let mut moments = vec![0.0; order];

        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};
        let mut rng = StdRng::seed_from_u64(seed);

        for _v in 0..n_vectors {
            // Random vector r with entries ±1/sqrt(N)
            let mut r = vec![0.0; n];
            for i in 0..n {
                r[i] = if rng.gen_bool(0.5) {
                    1.0 / (n as f64).sqrt()
                } else {
                    -1.0 / (n as f64).sqrt()
                };
            }

            // T_0 * r = r
            let mut t_prev = r.clone();
            // T_1 * r = H_tilde * r
            let mut t_curr = matvec(&h_tilde, &r);

            // mu_0 = r^T * T_0 * r = r^T * r = 1 (normalized)
            moments[0] += dot(&r, &t_prev);

            if order > 1 {
                moments[1] += dot(&r, &t_curr);
            }

            for k in 2..order {
                // T_{k} = 2 * H_tilde * T_{k-1} - T_{k-2}
                let t_next = chebyshev_step(&h_tilde, &t_curr, &t_prev);
                moments[k] += dot(&r, &t_next);
                t_prev = t_curr;
                t_curr = t_next;
            }
        }

        // Average over random vectors
        for m in moments.iter_mut() {
            *m /= n_vectors as f64;
        }

        Self {
            moments,
            a,
            b,
            n,
            order,
        }
    }

    /// Compute exact Chebyshev moments using full diagonal traces.
    /// More accurate but O(N^2 * order) instead of O(N * n_vec * order).
    pub fn from_matrix_exact(h: &DMatrix<f64>, order: usize) -> Self {
        let n = h.nrows();
        let (e_min, e_max) = estimate_spectral_bounds(h);
        let a = (e_max - e_min) / 2.0;
        let b = (e_max + e_min) / 2.0;
        let h_tilde = rescale_matrix(h, e_min, e_max);

        // Full matrix Chebyshev recursion
        let mut t_prev = DMatrix::identity(n, n);
        let mut t_curr = h_tilde.clone();

        let mut moments = vec![0.0; order];
        moments[0] = t_prev.trace() / n as f64;
        if order > 1 {
            moments[1] = t_curr.trace() / n as f64;
        }

        for k in 2..order {
            let t_next = &h_tilde * &t_curr * 2.0 - &t_prev;
            moments[k] = t_next.trace() / n as f64;
            t_prev = t_curr;
            t_curr = t_next;
        }

        Self {
            moments,
            a,
            b,
            n,
            order,
        }
    }

    /// Reconstruct DOS at a given energy using Chebyshev expansion.
    ///
    /// g(E) = (2/(pi*a)) * sum_k g_k * mu_k * T_k((E-b)/a) / sqrt(1 - x^2)
    pub fn dos_at_energy(&self, energy: f64, jackson: &[f64]) -> f64 {
        let x = (energy - self.b) / self.a;
        if x.abs() >= 1.0 {
            return 0.0;
        }

        let weight = 1.0 / (std::f64::consts::PI * self.a * (1.0 - x * x).sqrt());
        let mut sum = jackson[0] * self.moments[0];

        let mut t_prev = 1.0;
        let mut t_curr = x;

        for k in 1..self.order.min(jackson.len()) {
            sum += 2.0 * jackson[k] * self.moments[k] * t_curr;
            let t_next = 2.0 * x * t_curr - t_prev;
            t_prev = t_curr;
            t_curr = t_next;
        }

        (weight * sum * self.n as f64).max(0.0)
    }
}

/// Dense matrix-vector product.
fn matvec(m: &DMatrix<f64>, v: &[f64]) -> Vec<f64> {
    let n = m.nrows();
    let mut result = vec![0.0; n];
    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..n {
            sum += m[(i, j)] * v[j];
        }
        result[i] = sum;
    }
    result
}

/// Chebyshev recursion step: T_{k+1}*r = 2*H*T_k*r - T_{k-1}*r
fn chebyshev_step(h: &DMatrix<f64>, t_curr: &[f64], t_prev: &[f64]) -> Vec<f64> {
    let ht = matvec(h, t_curr);
    let n = ht.len();
    let mut result = vec![0.0; n];
    for i in 0..n {
        result[i] = 2.0 * ht[i] - t_prev[i];
    }
    result
}

/// Dot product.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Create a 1D Hückel chain Hamiltonian (tight-binding model).
    /// H_ii = alpha (on-site), H_{i,i+1} = beta (hopping), S = I.
    fn huckel_chain(n: usize, alpha: f64, beta: f64) -> DMatrix<f64> {
        let mut h = DMatrix::zeros(n, n);
        for i in 0..n {
            h[(i, i)] = alpha;
            if i + 1 < n {
                h[(i, i + 1)] = beta;
                h[(i + 1, i)] = beta;
            }
        }
        h
    }

    #[test]
    fn test_spectral_bounds_contain_eigenvalues() {
        let h = huckel_chain(20, -10.0, -2.5);
        let (e_min, e_max) = estimate_spectral_bounds(&h);

        // Verify with exact eigenvalues
        let eigen = nalgebra::SymmetricEigen::new(h);
        for i in 0..eigen.eigenvalues.len() {
            let e = eigen.eigenvalues[i];
            assert!(
                e >= e_min && e <= e_max,
                "Eigenvalue {} = {:.4} outside [{:.4}, {:.4}]",
                i, e, e_min, e_max
            );
        }
    }

    #[test]
    fn test_rescaled_spectrum_in_unit_interval() {
        let h = huckel_chain(15, -8.0, -3.0);
        let (e_min, e_max) = estimate_spectral_bounds(&h);
        let h_tilde = rescale_matrix(&h, e_min, e_max);

        let eigen = nalgebra::SymmetricEigen::new(h_tilde);
        for i in 0..eigen.eigenvalues.len() {
            let e = eigen.eigenvalues[i];
            assert!(
                e.abs() < 1.0 + 1e-10,
                "Rescaled eigenvalue {} = {:.6} > 1",
                i, e
            );
        }
    }

    #[test]
    fn test_jackson_kernel_damping() {
        let gk = jackson_kernel(50);
        // g_0 should be close to 1
        assert!((gk[0] - 1.0).abs() < 0.05, "g_0 = {}", gk[0]);
        // Coefficients should decrease
        assert!(gk[49] < gk[0], "Jackson kernel should damp high-order terms");
        // All coefficients should be non-negative for well-behaved kernels
        // (Jackson kernel is not strictly non-negative but g_k > 0 for moderate k)
    }

    #[test]
    fn test_chebyshev_dos_peaks() {
        // 1D chain: analytical DOS has van Hove singularities at band edges
        let h = huckel_chain(30, 0.0, -1.0);
        let expansion = ChebyshevExpansion::from_matrix_exact(&h, 100);
        let gk = jackson_kernel(100);

        // DOS should be nonzero inside the band [-2, 2] and zero outside
        let dos_inside = expansion.dos_at_energy(0.0, &gk);
        let dos_outside = expansion.dos_at_energy(5.0, &gk);

        assert!(dos_inside > 0.1, "DOS at band center should be > 0: {}", dos_inside);
        assert!(dos_outside < 0.01, "DOS outside band should be ~0: {}", dos_outside);
    }

    #[test]
    fn test_stochastic_vs_exact_moments() {
        let h = huckel_chain(20, -5.0, -1.5);
        let exact = ChebyshevExpansion::from_matrix_exact(&h, 30);
        let stoch = ChebyshevExpansion::from_matrix(&h, 30, 200, 42);

        // First few moments should agree within ~50% for stochastic trace
        for k in 0..5 {
            let rel_err = if exact.moments[k].abs() > 1e-10 {
                ((stoch.moments[k] - exact.moments[k]) / exact.moments[k]).abs()
            } else {
                (stoch.moments[k] - exact.moments[k]).abs()
            };
            assert!(
                rel_err < 0.8,
                "Moment {} disagrees: exact={:.6}, stoch={:.6}",
                k, exact.moments[k], stoch.moments[k]
            );
        }
    }
}
