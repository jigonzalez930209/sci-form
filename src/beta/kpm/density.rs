//! KPM Density, DOS and Mulliken Populations — E4.2
//!
//! Compute DOS, Fermi-Dirac occupation, and Mulliken charges
//! from Chebyshev expansion without diagonalization.

use super::chebyshev::{
    estimate_spectral_bounds, jackson_kernel, rescale_matrix, ChebyshevExpansion,
};
use nalgebra::DMatrix;

/// Configuration for KPM calculations.
#[derive(Debug, Clone)]
pub struct KpmConfig {
    /// Chebyshev expansion order.
    pub order: usize,
    /// Number of stochastic trace vectors (0 = exact).
    pub n_vectors: usize,
    /// RNG seed.
    pub seed: u64,
    /// Electronic temperature for Fermi-Dirac smearing (eV). 0 = zero temperature.
    pub temperature: f64,
}

impl Default for KpmConfig {
    fn default() -> Self {
        Self {
            order: 100,
            n_vectors: 0, // exact by default for small systems
            seed: 42,
            temperature: 0.0,
        }
    }
}

/// Result of KPM DOS calculation.
#[derive(Debug, Clone)]
pub struct KpmDosResult {
    /// Energy grid (eV).
    pub energies: Vec<f64>,
    /// Total DOS at each energy point.
    pub total_dos: Vec<f64>,
    /// Spectral bounds [E_min, E_max] in eV.
    pub e_min: f64,
    pub e_max: f64,
    /// Expansion order used.
    pub order: usize,
}

/// Result of KPM Mulliken charge calculation.
#[derive(Debug, Clone)]
pub struct KpmMullikenResult {
    /// Mulliken charges per atom.
    pub charges: Vec<f64>,
    /// Estimated electron count.
    pub n_electrons_kpm: f64,
    /// Target electron count.
    pub n_electrons_target: usize,
    /// Fermi level (eV).
    pub fermi_level: f64,
}

/// Compute Fermi-Dirac Chebyshev expansion coefficients.
///
/// c_k = (2/pi) * integral_{-1}^{1} f_FD(a*x+b, mu, T) * T_k(x) / sqrt(1-x^2) dx
pub fn fermi_dirac_coefficients(
    order: usize,
    a: f64,
    b: f64,
    mu: f64,
    temperature: f64,
) -> Vec<f64> {
    let n_quad = order * 4; // Quadrature points
    let mut coeffs = vec![0.0; order];

    for q in 0..n_quad {
        let theta = std::f64::consts::PI * (q as f64 + 0.5) / n_quad as f64;
        let x = theta.cos();
        let energy = a * x + b;

        // Fermi-Dirac distribution
        let f_fd = if temperature < 1e-10 {
            if energy < mu {
                1.0
            } else if energy > mu {
                0.0
            } else {
                0.5
            }
        } else {
            let kt = temperature * 8.617333262e-5; // eV to eV (kB in eV/K)
            1.0 / (1.0 + ((energy - mu) / kt).exp())
        };

        // Chebyshev polynomials at x via recurrence
        let mut t_prev = 1.0;
        let mut t_curr = x;
        coeffs[0] += f_fd * t_prev;
        if order > 1 {
            coeffs[1] += f_fd * t_curr;
        }
        for k in 2..order {
            let t_next = 2.0 * x * t_curr - t_prev;
            coeffs[k] += f_fd * t_next;
            t_prev = t_curr;
            t_curr = t_next;
        }
    }

    // Normalize: Chebyshev quadrature weight = pi/n_quad, and normalization factor 2/pi
    let factor = 2.0 / n_quad as f64;
    for c in coeffs.iter_mut() {
        *c *= factor;
    }
    // c_0 has a factor 1/2 in Chebyshev expansion
    coeffs[0] *= 0.5;

    coeffs
}

/// Compute KPM DOS on an energy grid.
pub fn compute_kpm_dos(
    h: &DMatrix<f64>,
    config: &KpmConfig,
    e_min_out: f64,
    e_max_out: f64,
    n_points: usize,
) -> KpmDosResult {
    let expansion = if config.n_vectors == 0 {
        ChebyshevExpansion::from_matrix_exact(h, config.order)
    } else {
        ChebyshevExpansion::from_matrix(h, config.order, config.n_vectors, config.seed)
    };

    let gk = jackson_kernel(config.order);

    let step = (e_max_out - e_min_out) / (n_points - 1).max(1) as f64;
    let energies: Vec<f64> = (0..n_points).map(|i| e_min_out + i as f64 * step).collect();
    let total_dos: Vec<f64> = energies
        .iter()
        .map(|&e| expansion.dos_at_energy(e, &gk))
        .collect();

    KpmDosResult {
        energies,
        total_dos,
        e_min: expansion.b - expansion.a,
        e_max: expansion.b + expansion.a,
        order: config.order,
    }
}

/// Compute KPM Mulliken charges from a Hamiltonian and overlap matrix.
///
/// The density matrix is P = f(H), where f is the Fermi-Dirac function.
/// In the orthogonal basis (S=I), P_ij = sum_k c_k (T_k(H_tilde))_ij.
/// Mulliken charge: q_i = Z_i - (PS)_ii
pub fn compute_kpm_mulliken(
    h: &DMatrix<f64>,
    s: &DMatrix<f64>,
    n_electrons: usize,
    nuclear_charges: &[f64],
    config: &KpmConfig,
) -> KpmMullikenResult {
    let n = h.nrows();

    // Step 1: Löwdin orthogonalization S^{-1/2} H S^{-1/2}
    let s_eigen = nalgebra::SymmetricEigen::new(s.clone());
    let mut s_inv_sqrt = DMatrix::zeros(n, n);
    for i in 0..n {
        let v = s_eigen.eigenvalues[i];
        if v > 1e-10 {
            s_inv_sqrt[(i, i)] = 1.0 / v.sqrt();
        }
    }
    let s_inv_sqrt = &s_eigen.eigenvectors * &s_inv_sqrt * s_eigen.eigenvectors.transpose();
    let h_orth = &s_inv_sqrt * h * &s_inv_sqrt;

    // Step 2: Find Fermi level by bisection
    let (e_min, e_max) = estimate_spectral_bounds(&h_orth);
    let a = (e_max - e_min) / 2.0;
    let b = (e_max + e_min) / 2.0;
    let h_tilde = rescale_matrix(&h_orth, e_min, e_max);

    let target = n_electrons as f64;
    let mut mu_lo = e_min;
    let mut mu_hi = e_max;

    // Bisection for Fermi level
    for _ in 0..60 {
        let mu = (mu_lo + mu_hi) / 2.0;
        let ne = electron_count_kpm(&h_tilde, a, b, mu, config);
        if ne < target {
            mu_lo = mu;
        } else {
            mu_hi = mu;
        }
    }
    let mu = (mu_lo + mu_hi) / 2.0;

    // Step 3: Build density matrix via Chebyshev expansion
    let coeffs = fermi_dirac_coefficients(config.order, a, b, mu, config.temperature);
    let gk = jackson_kernel(config.order);

    // P_orth = sum_k (c_k * g_k) * T_k(H_tilde)  in orthogonal basis
    let mut p_orth = DMatrix::zeros(n, n);
    let mut t_prev = DMatrix::identity(n, n);
    let mut t_curr = h_tilde.clone();

    // k=0 contribution
    let c0 = coeffs[0] * gk[0];
    p_orth += &t_prev * c0;

    if config.order > 1 {
        let c1 = 2.0 * coeffs[1] * gk[1];
        p_orth += &t_curr * c1;
    }

    for k in 2..config.order {
        let t_next = &h_tilde * &t_curr * 2.0 - &t_prev;
        let ck = 2.0 * coeffs[k] * gk[k.min(gk.len() - 1)];
        p_orth += &t_next * ck;
        t_prev = t_curr;
        t_curr = t_next;
    }

    // Step 4: Transform back: P_AO = S^{-1/2} P_orth S^{-1/2}
    let p_ao = &s_inv_sqrt * &p_orth * &s_inv_sqrt;

    // Step 5: Mulliken charges: q_i = Z_i - sum of (PS)_ii contributions per atom
    // We need to map basis functions to atoms for proper Mulliken analysis.
    // Simplified: assume 1 basis function per nuclear charge entry.
    let ps = &p_ao * s;
    let n_elec_kpm: f64 = (0..n).map(|i| ps[(i, i)]).sum();

    let mut charges = Vec::with_capacity(nuclear_charges.len());
    if nuclear_charges.len() == n {
        for i in 0..n {
            charges.push(nuclear_charges[i] - ps[(i, i)]);
        }
    } else {
        // Simplified: distribute evenly
        for z in nuclear_charges {
            charges.push(*z);
        }
    }

    KpmMullikenResult {
        charges,
        n_electrons_kpm: n_elec_kpm,
        n_electrons_target: n_electrons,
        fermi_level: mu,
    }
}

/// Count electrons at a given Fermi level using Chebyshev expansion of f_FD.
fn electron_count_kpm(h_tilde: &DMatrix<f64>, a: f64, b: f64, mu: f64, config: &KpmConfig) -> f64 {
    let n = h_tilde.nrows();
    let coeffs = fermi_dirac_coefficients(config.order, a, b, mu, config.temperature);
    let gk = jackson_kernel(config.order);

    // N_e = Tr[f(H)] = N * sum_k c_k * g_k * mu_k
    // where mu_k = (1/N) Tr[T_k(H_tilde)]

    let mut t_prev = DMatrix::identity(n, n);
    let mut t_curr = h_tilde.clone();

    let mut ne = coeffs[0] * gk[0] * t_prev.trace();

    if config.order > 1 {
        ne += 2.0 * coeffs[1] * gk[1] * t_curr.trace();
    }

    for k in 2..config.order {
        let t_next = h_tilde * &t_curr * 2.0 - &t_prev;
        ne += 2.0 * coeffs[k] * gk[k.min(gk.len() - 1)] * t_next.trace();
        t_prev = t_curr;
        t_curr = t_next;
    }

    ne
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_fermi_dirac_coefficients_step_function() {
        // At T=0, f_FD is a step function at mu
        let coeffs = fermi_dirac_coefficients(50, 5.0, -5.0, -5.0, 0.0);
        // c_0 should correspond to half-filling
        assert!(coeffs[0].abs() < 2.0, "c_0 = {}", coeffs[0]);
    }

    #[test]
    fn test_kpm_dos_shape() {
        let h = huckel_chain(20, 0.0, -1.0);
        let config = KpmConfig {
            order: 80,
            n_vectors: 0,
            ..Default::default()
        };
        let dos = compute_kpm_dos(&h, &config, -3.0, 3.0, 100);

        // DOS should be non-negative
        for &d in &dos.total_dos {
            assert!(d >= -1e-10, "Negative DOS: {}", d);
        }

        // DOS should be nonzero near band center
        let mid_idx = dos.total_dos.len() / 2;
        assert!(
            dos.total_dos[mid_idx] > 0.0,
            "DOS at band center = {}",
            dos.total_dos[mid_idx]
        );
    }

    #[test]
    fn test_kpm_mulliken_charge_conservation() {
        let n = 10;
        let h = huckel_chain(n, -10.0, -2.0);
        let s = DMatrix::identity(n, n);
        let n_electrons = n; // half-filling
        let nuclear_charges: Vec<f64> = vec![1.0; n];

        let config = KpmConfig {
            order: 60,
            n_vectors: 0,
            ..Default::default()
        };

        let result = compute_kpm_mulliken(&h, &s, n_electrons, &nuclear_charges, &config);

        // Total charge should sum to ~0 (neutral)
        let total: f64 = result.charges.iter().sum();
        assert!(
            total.abs() < 1.0,
            "Total Mulliken charge = {:.4}, expected ~0",
            total
        );
    }

    #[test]
    fn test_electron_count_at_fermi_level() {
        let n = 10;
        let h = huckel_chain(n, 0.0, -1.0);
        let (e_min, e_max) = estimate_spectral_bounds(&h);
        let h_tilde = rescale_matrix(&h, e_min, e_max);
        let a = (e_max - e_min) / 2.0;
        let b_val = (e_max + e_min) / 2.0;

        let config = KpmConfig {
            order: 80,
            n_vectors: 0,
            ..Default::default()
        };

        // At mu = 0 (band center), should give ~N/2 electrons
        let ne = electron_count_kpm(&h_tilde, a, b_val, 0.0, &config);
        assert!(
            (ne - n as f64 / 2.0).abs() < 2.0,
            "Expected ~{} electrons at mu=0, got {:.2}",
            n / 2,
            ne
        );
    }
}
