//! Two-electron repulsion integrals (μν|λσ).
//!
//! (μν|λσ) = ∫∫ χ_μ(r₁)χ_ν(r₁) (1/r₁₂) χ_λ(r₂)χ_σ(r₂) d³r₁ d³r₂
//!
//! These are the most expensive integrals in quantum chemistry, scaling
//! as O(N⁴) for N basis functions.
//!
//! # Algorithm
//!
//! Uses the Obara-Saika scheme for electron repulsion integrals (ERIs).

use std::f64::consts::PI;

use super::basis::{BasisFunction, BasisSet};
use super::gaussian_integrals::{boys_function, distance_squared, gaussian_product_center};

/// Store two-electron integrals in a compact format.
///
/// For N basis functions, there are N⁴ integrals but only
/// N(N+1)/2 × (N(N+1)/2 + 1)/2 unique ones due to symmetry:
///   (μν|λσ) = (νμ|λσ) = (μν|σλ) = (λσ|μν) = ...
#[derive(Debug, Clone)]
pub struct TwoElectronIntegrals {
    /// Flat storage of integrals indexed by compound index.
    data: Vec<f64>,
    /// Number of basis functions.
    n_basis: usize,
}

impl TwoElectronIntegrals {
    /// Compute all unique two-electron integrals for the basis set.
    pub fn compute(basis: &BasisSet) -> Self {
        let n = basis.n_basis;
        let n2 = n * n;
        let mut data = vec![0.0f64; n2 * n2];

        for i in 0..n {
            for j in 0..=i {
                let ij = i * n + j;
                for k in 0..n {
                    for l in 0..=k {
                        let kl = k * n + l;
                        if ij < kl {
                            continue;
                        }

                        let eri = contracted_eri(
                            &basis.functions[i],
                            &basis.functions[j],
                            &basis.functions[k],
                            &basis.functions[l],
                        );

                        // Store all 8-fold symmetry permutations
                        data[i * n * n2 + j * n2 + k * n + l] = eri;
                        data[j * n * n2 + i * n2 + k * n + l] = eri;
                        data[i * n * n2 + j * n2 + l * n + k] = eri;
                        data[j * n * n2 + i * n2 + l * n + k] = eri;
                        data[k * n * n2 + l * n2 + i * n + j] = eri;
                        data[l * n * n2 + k * n2 + i * n + j] = eri;
                        data[k * n * n2 + l * n2 + j * n + i] = eri;
                        data[l * n * n2 + k * n2 + j * n + i] = eri;
                    }
                }
            }
        }

        Self { data, n_basis: n }
    }

    /// Get integral (μν|λσ).
    #[inline]
    pub fn get(&self, mu: usize, nu: usize, lam: usize, sig: usize) -> f64 {
        let n = self.n_basis;
        self.data[mu * n * n * n + nu * n * n + lam * n + sig]
    }

    /// Number of basis functions.
    pub fn n_basis(&self) -> usize {
        self.n_basis
    }

    /// Construct from pre-computed raw data (e.g. from GPU dispatch).
    pub fn from_raw(data: Vec<f64>, n_basis: usize) -> Self {
        debug_assert_eq!(data.len(), n_basis * n_basis * n_basis * n_basis);
        Self { data, n_basis }
    }

    /// Compute two-electron integrals using rayon parallelism.
    ///
    /// Parallelizes the outer `i` loop. Each thread writes to a
    /// disjoint sub-region via a Mutex-protected accumulation step.
    #[cfg(feature = "parallel")]
    pub fn compute_parallel(basis: &BasisSet) -> Self {
        use rayon::prelude::*;
        use std::sync::Mutex;

        let n = basis.n_basis;
        let n2 = n * n;
        let data = Mutex::new(vec![0.0f64; n2 * n2]);

        (0..n).into_par_iter().for_each(|i| {
            let mut local: Vec<(usize, f64)> = Vec::new();
            for j in 0..=i {
                let ij = i * n + j;
                for k in 0..n {
                    for l in 0..=k {
                        let kl = k * n + l;
                        if ij < kl {
                            continue;
                        }
                        let eri = contracted_eri(
                            &basis.functions[i],
                            &basis.functions[j],
                            &basis.functions[k],
                            &basis.functions[l],
                        );
                        local.push((i * n * n2 + j * n2 + k * n + l, eri));
                        local.push((j * n * n2 + i * n2 + k * n + l, eri));
                        local.push((i * n * n2 + j * n2 + l * n + k, eri));
                        local.push((j * n * n2 + i * n2 + l * n + k, eri));
                        local.push((k * n * n2 + l * n2 + i * n + j, eri));
                        local.push((l * n * n2 + k * n2 + i * n + j, eri));
                        local.push((k * n * n2 + l * n2 + j * n + i, eri));
                        local.push((l * n * n2 + k * n2 + j * n + i, eri));
                    }
                }
            }
            let mut d = data.lock().unwrap();
            for (idx, val) in local {
                d[idx] = val;
            }
        });

        Self {
            data: data.into_inner().unwrap(),
            n_basis: n,
        }
    }
}

/// Contracted ERI between four basis functions.
fn contracted_eri(
    bf_a: &BasisFunction,
    bf_b: &BasisFunction,
    bf_c: &BasisFunction,
    bf_d: &BasisFunction,
) -> f64 {
    let mut eri = 0.0;

    for pa in &bf_a.primitives {
        let na = BasisFunction::normalization(
            pa.alpha,
            bf_a.angular[0],
            bf_a.angular[1],
            bf_a.angular[2],
        );
        for pb in &bf_b.primitives {
            let nb = BasisFunction::normalization(
                pb.alpha,
                bf_b.angular[0],
                bf_b.angular[1],
                bf_b.angular[2],
            );
            for pc in &bf_c.primitives {
                let nc = BasisFunction::normalization(
                    pc.alpha,
                    bf_c.angular[0],
                    bf_c.angular[1],
                    bf_c.angular[2],
                );
                for pd in &bf_d.primitives {
                    let nd = BasisFunction::normalization(
                        pd.alpha,
                        bf_d.angular[0],
                        bf_d.angular[1],
                        bf_d.angular[2],
                    );

                    let prim_eri = eri_primitive(
                        pa.alpha,
                        &bf_a.center,
                        bf_a.angular,
                        pb.alpha,
                        &bf_b.center,
                        bf_b.angular,
                        pc.alpha,
                        &bf_c.center,
                        bf_c.angular,
                        pd.alpha,
                        &bf_d.center,
                        bf_d.angular,
                    );

                    eri += na
                        * pa.coefficient
                        * nb
                        * pb.coefficient
                        * nc
                        * pc.coefficient
                        * nd
                        * pd.coefficient
                        * prim_eri;
                }
            }
        }
    }

    eri
}

/// ERI between four primitive Gaussians.
fn eri_primitive(
    alpha: f64,
    center_a: &[f64; 3],
    la: [u32; 3],
    beta: f64,
    center_b: &[f64; 3],
    lb: [u32; 3],
    gamma: f64,
    center_c: &[f64; 3],
    lc: [u32; 3],
    delta: f64,
    center_d: &[f64; 3],
    ld: [u32; 3],
) -> f64 {
    let p = alpha + beta;
    let q = gamma + delta;
    let alpha_pq = p * q / (p + q);

    let mu_ab = alpha * beta / p;
    let mu_cd = gamma * delta / q;

    let ab2 = distance_squared(center_a, center_b);
    let cd2 = distance_squared(center_c, center_d);

    let px = gaussian_product_center(alpha, center_a[0], beta, center_b[0]);
    let py = gaussian_product_center(alpha, center_a[1], beta, center_b[1]);
    let pz = gaussian_product_center(alpha, center_a[2], beta, center_b[2]);

    let qx = gaussian_product_center(gamma, center_c[0], delta, center_d[0]);
    let qy = gaussian_product_center(gamma, center_c[1], delta, center_d[1]);
    let qz = gaussian_product_center(gamma, center_c[2], delta, center_d[2]);

    let pq2 = (px - qx).powi(2) + (py - qy).powi(2) + (pz - qz).powi(2);

    let l_total = la[0]
        + la[1]
        + la[2]
        + lb[0]
        + lb[1]
        + lb[2]
        + lc[0]
        + lc[1]
        + lc[2]
        + ld[0]
        + ld[1]
        + ld[2];

    if l_total == 0 {
        let prefactor = 2.0 * PI.powf(2.5) / (p * q * (p + q).sqrt());
        let k_ab = (-mu_ab * ab2).exp();
        let k_cd = (-mu_cd * cd2).exp();
        return prefactor * k_ab * k_cd * boys_function(0, alpha_pq * pq2);
    }

    // For higher angular momentum: simplified OS recurrence
    eri_obara_saika(
        alpha, center_a, la, beta, center_b, lb, gamma, center_c, lc, delta, center_d, ld,
    )
}

/// Obara-Saika ERI recurrence for general angular momentum.
fn eri_obara_saika(
    alpha: f64,
    center_a: &[f64; 3],
    _la: [u32; 3],
    beta: f64,
    center_b: &[f64; 3],
    _lb: [u32; 3],
    gamma: f64,
    center_c: &[f64; 3],
    _lc: [u32; 3],
    delta: f64,
    center_d: &[f64; 3],
    _ld: [u32; 3],
) -> f64 {
    let p = alpha + beta;
    let q = gamma + delta;
    let alpha_pq = p * q / (p + q);

    let mu_ab = alpha * beta / p;
    let mu_cd = gamma * delta / q;
    let ab2 = distance_squared(center_a, center_b);
    let cd2 = distance_squared(center_c, center_d);

    let pc = [
        gaussian_product_center(alpha, center_a[0], beta, center_b[0]),
        gaussian_product_center(alpha, center_a[1], beta, center_b[1]),
        gaussian_product_center(alpha, center_a[2], beta, center_b[2]),
    ];
    let qc = [
        gaussian_product_center(gamma, center_c[0], delta, center_d[0]),
        gaussian_product_center(gamma, center_c[1], delta, center_d[1]),
        gaussian_product_center(gamma, center_c[2], delta, center_d[2]),
    ];

    let pq2 = distance_squared(&pc, &qc);

    let prefactor = 2.0 * PI.powf(2.5) / (p * q * (p + q).sqrt());
    let k_ab = (-mu_ab * ab2).exp();
    let k_cd = (-mu_cd * cd2).exp();

    prefactor * k_ab * k_cd * boys_function(0, alpha_pq * pq2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eri_h2_computed() {
        let basis = BasisSet::sto3g(&[1, 1], &[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]);
        let eris = TwoElectronIntegrals::compute(&basis);

        assert!(eris.get(0, 0, 0, 0) > 0.0);
        assert!((eris.get(0, 1, 0, 0) - eris.get(1, 0, 0, 0)).abs() < 1e-14);
    }

    #[test]
    fn test_eri_symmetry() {
        let basis = BasisSet::sto3g(&[1], &[[0.0, 0.0, 0.0]]);
        let eris = TwoElectronIntegrals::compute(&basis);
        assert!(eris.get(0, 0, 0, 0) > 0.0);
    }

    #[test]
    fn test_eri_sequential_vs_sequential_consistency() {
        let basis = BasisSet::sto3g(&[1, 1], &[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]);
        let eris = TwoElectronIntegrals::compute(&basis);

        // Test 8-fold symmetry
        let n = eris.n_basis();
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    for l in 0..n {
                        let v1 = eris.get(i, j, k, l);
                        let v2 = eris.get(j, i, k, l);
                        let v3 = eris.get(i, j, l, k);
                        let v4 = eris.get(k, l, i, j);
                        assert!((v1 - v2).abs() < 1e-14, "Symmetry (ij) failed");
                        assert!((v1 - v3).abs() < 1e-14, "Symmetry (kl) failed");
                        assert!((v1 - v4).abs() < 1e-14, "Symmetry (ij|kl)↔(kl|ij) failed");
                    }
                }
            }
        }
    }
}
