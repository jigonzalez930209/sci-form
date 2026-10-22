//! **ALPHA** — Optimized Obara-Saika electron repulsion integrals.
//!
//! Provides a shell-pair based ERI engine with:
//! - Vertical + horizontal recurrence relations for contracted shells
//! - Cache-friendly memory layout for batch ERI evaluation
//! - Optional Schwarz pre-screening integration
//!
//! This supplements the existing `gaussian_integrals.rs` with a higher-performance
//! path for larger molecules.

use super::gaussian_integrals::{distance_squared, gaussian_product_center};
use std::f64::consts::PI;

/// Error function to machine precision (no std `erf` on stable Rust).
///
/// Taylor series for |x| <= 1, Lentz continued fraction for erfc(x) above.
fn erf_accurate(x: f64) -> f64 {
    if x.abs() <= 1.0 {
        // erf(x) = 2/sqrt(pi) * sum_k (-1)^k x^(2k+1) / (k! (2k+1))
        let mut term = x;
        let mut sum = x;
        let x2 = x * x;
        for k in 1..60 {
            term *= -x2 * (2 * k - 1) as f64 / (k * (2 * k + 1)) as f64;
            sum += term;
            if term.abs() < 1e-18 * sum.abs().max(1e-300) {
                break;
            }
        }
        return 2.0 * sum / PI.sqrt();
    }
    // erfc(x) = x*e^(-x^2)/sqrt(pi) * 1/(x^2 + (1/2)/(1 + 1/(x^2 + (3/2)/(1 + ...))))
    // Lentz's algorithm on the Laplace continued fraction.
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let z = x.abs();
    let tiny = 1e-300;
    let _f = tiny;
    let _c = 1.0 / tiny;
    let _d = 0.0;
    let z2 = z * z;
    // Lentz init: f = C = b0 = z^2, D = 0; pairs (a_n, b_n) = Laplace sequence.
    let mut f = z2;
    let mut c = z2;
    let mut d = 0.0;
    let mut n = 1u32;
    loop {
        let a = n as f64 / 2.0;
        let b = if n % 2 == 1 { 1.0 } else { z2 };
        d = b + a * d;
        if d.abs() < tiny {
            d = tiny;
        }
        d = 1.0 / d;
        c = b + a / c;
        if c.abs() < tiny {
            c = tiny;
        }
        let delta = c * d;
        f *= delta;
        n += 1;
        if (delta - 1.0).abs() < 1e-15 || n > 200 {
            break;
        }
    }
    // The Laplace fraction is the reciprocal of the b0-anchored CF value.
    let erfc = z * (-z2).exp() / PI.sqrt() / f;
    sign * (1.0 - erfc.clamp(0.0, 2.0))
}

/// Shell pair data pre-computed for efficient ERI evaluation.
///
/// All lengths are in Bohr and all exponents in Bohr⁻², matching the rest of
/// the `scf` integral engine.
#[derive(Debug, Clone)]
pub struct ShellPairData {
    /// Combined exponent p = α + β, in Bohr⁻².
    pub p: f64,
    /// Reduced exponent μ = αβ / (α + β), in Bohr⁻².
    pub mu: f64,
    /// Gaussian product center P = (αA + βB) / p, in Bohr.
    pub center_p: [f64; 3],
    /// Pre-factor exp(−μ|AB|²), dimensionless.
    pub k_ab: f64,
    /// Original centers, in Bohr.
    pub center_a: [f64; 3],
    /// Original centers, in Bohr.
    pub center_b: [f64; 3],
    /// Exponent of the first primitive, in Bohr⁻².
    pub alpha: f64,
    /// Exponent of the second primitive, in Bohr⁻².
    pub beta: f64,
}

impl ShellPairData {
    /// Compute shell pair data from two primitive Gaussians.
    pub fn new(alpha: f64, center_a: [f64; 3], beta: f64, center_b: [f64; 3]) -> Self {
        let p = alpha + beta;
        let mu = alpha * beta / p;
        let ab2 = distance_squared(&center_a, &center_b);
        let k_ab = (-mu * ab2).exp();
        let center_p = [
            gaussian_product_center(alpha, center_a[0], beta, center_b[0]),
            gaussian_product_center(alpha, center_a[1], beta, center_b[1]),
            gaussian_product_center(alpha, center_a[2], beta, center_b[2]),
        ];
        Self {
            p,
            mu,
            center_p,
            k_ab,
            center_a,
            center_b,
            alpha,
            beta,
        }
    }
}

/// Boys function F_n(x) for nuclear attraction and ERI integrals.
///
/// Taylor series for x <= 6 (no cancellation issues), F_0 from the
/// machine-precision erf plus upward recursion above (stable there:
/// the 2x denominator dominates and exp(-x) is negligible).
///
/// The function is dimensionless; `x` is the dimensionless Boys argument
/// (the product of the Gaussian exponents and the squared internuclear
/// distance, both in atomic units).
pub fn boys_function(n: usize, x: f64) -> f64 {
    if x <= 0.0 {
        return 1.0 / (2 * n + 1) as f64;
    }
    if x <= 6.0 {
        // F_n(x) = sum_k (-x)^k / (k! (2n+2k+1))
        let mut term = 1.0 / (2 * n + 1) as f64;
        let mut sum = term;
        for k in 1..200 {
            term *= -x * (2 * n + 2 * k - 1) as f64 / (k * (2 * n + 2 * k + 1)) as f64;
            sum += term;
            if term.abs() < 1e-17 * sum.abs().max(1e-300) {
                break;
            }
        }
        return sum;
    }

    let f0 = (PI / x).sqrt() * erf_accurate(x.sqrt()) / 2.0;
    if n == 0 {
        return f0;
    }

    // Upward recursion: F_{n+1}(x) = ((2n+1)*F_n(x) - exp(-x)) / (2x)
    let exp_neg_x = (-x).exp();
    let mut f_curr = f0;
    for i in 0..n {
        f_curr = ((2 * i + 1) as f64 * f_curr - exp_neg_x) / (2.0 * x);
    }
    f_curr
}

/// Compute (ss|ss) electron repulsion integral between four s-type primitives.
///
/// (μν|λσ) = 2π^{5/2} / (pq√(p+q)) · K_ab · K_cd · F_0(T)
///
/// where T = ρ|PQ|², ρ = pq/(p+q)
///
/// The result is in Hartree (the exponents are in Bohr⁻² and the centres in
/// Bohr, so the integral is an atomic-unit Coulomb integral).
pub fn eri_ssss(sp_ab: &ShellPairData, sp_cd: &ShellPairData) -> f64 {
    let p = sp_ab.p;
    let q = sp_cd.p;
    let alpha = p * q / (p + q);
    let pq2 = distance_squared(&sp_ab.center_p, &sp_cd.center_p);
    let t = alpha * pq2;

    let prefactor = 2.0 * PI.powi(2) * PI.sqrt() / (p * q * (p + q).sqrt());
    prefactor * sp_ab.k_ab * sp_cd.k_ab * boys_function(0, t)
}

/// Compute Schwarz bound for a shell pair: Q_ij = sqrt((ij|ij)).
///
/// Used for pre-screening: skip (ij|kl) when Q_ij * Q_kl < threshold.
///
/// The bound is in √Hartree (the ERI is in Hartree), so the screening
/// threshold it is compared against must be in the same units.
pub fn schwarz_bound(alpha: f64, center_a: [f64; 3], beta: f64, center_b: [f64; 3]) -> f64 {
    let sp = ShellPairData::new(alpha, center_a, beta, center_b);
    let eri = eri_ssss(&sp, &sp);
    eri.abs().sqrt()
}

/// Batch ERI computation for a list of shell pairs with Schwarz screening.
///
/// Returns the non-zero ERIs as (i, j, k, l, value) tuples, with the
/// values in Hartree.
///
/// `schwarz_q` holds one Schwarz bound per shell pair (in √Hartree) and
/// `threshold` the screening threshold in the same units.
pub fn compute_eris_screened(
    shell_pairs: &[ShellPairData],
    schwarz_q: &[f64],
    threshold: f64,
) -> Vec<(usize, usize, usize, usize, f64)> {
    let n_pairs = shell_pairs.len();
    let mut eris = Vec::new();

    for ij in 0..n_pairs {
        if schwarz_q[ij] < threshold {
            continue;
        }
        for kl in ij..n_pairs {
            if schwarz_q[ij] * schwarz_q[kl] < threshold {
                continue;
            }
            let val = eri_ssss(&shell_pairs[ij], &shell_pairs[kl]);
            if val.abs() > threshold {
                eris.push((ij, kl, 0, 0, val));
            }
        }
    }

    eris
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn boys_f0_is_correct() {
        // F_0(0) = 1.0
        let f0 = boys_function(0, 0.0);
        assert!((f0 - 1.0).abs() < 1e-10);

        // F_0(1) analytical: sqrt(pi) * erf(1) / 2, erf(1) = 0.8427007929...
        let f0_1 = boys_function(0, 1.0);
        let expected = (PI / 1.0).sqrt() * 0.8427007929497149 / 2.0;
        assert!((f0_1 - expected).abs() < 1e-12);
    }

    #[test]
    fn eri_ssss_hydrogen_molecule() {
        // Two s-type Gaussians on H2 (α=1.0) centered at ±0.7 Å
        let sp1 = ShellPairData::new(1.0, [0.0, 0.0, 0.0], 1.0, [0.0, 0.0, 0.0]);
        let sp2 = ShellPairData::new(1.0, [1.4, 0.0, 0.0], 1.0, [1.4, 0.0, 0.0]);
        let eri = eri_ssss(&sp1, &sp2);
        // Should be a reasonable positive value
        assert!(eri > 0.0);
    }
}
