//! Primitive Gaussian integral evaluation using Obara-Saika recursion.
//!
//! Implements the Obara-Saika (1986) recurrence relation for computing
//! overlap, kinetic, and nuclear attraction integrals between Cartesian
//! Gaussian primitives:
//!
//!   g(α, A, l) = (x-Ax)^lx (y-Ay)^ly (z-Az)^lz exp(-α|r-A|²)
//!
//! # Algorithm
//!
//! The Gaussian product theorem states that the product of two Gaussians
//! centered at A and B is a Gaussian centered at P:
//!
//!   P = (α·A + β·B) / (α + β)
//!   μ = α·β / (α + β)
//!
//! The Obara-Saika recurrence builds higher angular momentum integrals
//! from lower ones using:
//!
//!   S(a+1,b) = (P-A)·S(a,b) + (1/2p)·[a·S(a-1,b) + b·S(a,b-1)]
//!
//! where p = α + β and S(0,0) = sqrt(π/p) · exp(-μ·|AB|²)

use std::f64::consts::PI;

/// Gaussian product center P = (α·A + β·B) / (α + β).
#[inline]
pub fn gaussian_product_center(alpha: f64, a: f64, beta: f64, b: f64) -> f64 {
    (alpha * a + beta * b) / (alpha + beta)
}

/// Distance squared between two 3D points.
#[inline]
pub fn distance_squared(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

/// Overlap integral between two s-type (l=0) Gaussian primitives.
///
/// S_00 = (π / p)^{3/2} · exp(-μ · |AB|²)
pub fn overlap_ss(alpha: f64, center_a: &[f64; 3], beta: f64, center_b: &[f64; 3]) -> f64 {
    let p = alpha + beta;
    let mu = alpha * beta / p;
    let ab2 = distance_squared(center_a, center_b);
    (PI / p).powf(1.5) * (-mu * ab2).exp()
}

/// 1D overlap integral using Obara-Saika recursion.
///
/// Returns S(la, lb)_x = ∫ (x-Ax)^la (x-Bx)^lb exp(-α(x-Ax)² - β(x-Bx)²) dx
pub fn overlap_1d(la: u32, lb: u32, pa_x: f64, pb_x: f64, p: f64, prefactor: f64) -> f64 {
    let max_a = (la + 1) as usize;
    let max_b = (lb + 1) as usize;
    let mut table = vec![vec![0.0f64; max_b]; max_a];

    table[0][0] = prefactor;

    if max_a > 1 {
        table[1][0] = pa_x * table[0][0];
    }
    for a in 2..max_a {
        table[a][0] = pa_x * table[a - 1][0] + (a - 1) as f64 / (2.0 * p) * table[a - 2][0];
    }

    for b in 1..max_b {
        for a in 0..max_a {
            table[a][b] = pb_x * table[a][b - 1];
            if b > 1 {
                table[a][b] += (b - 1) as f64 / (2.0 * p) * table[a][b - 2];
            }
            if a > 0 {
                table[a][b] += a as f64 / (2.0 * p) * table[a - 1][b - 1];
            }
        }
    }

    table[la as usize][lb as usize]
}

/// Full 3D overlap integral between two Cartesian Gaussian primitives.
pub fn overlap_primitive(
    alpha: f64,
    center_a: &[f64; 3],
    la: [u32; 3],
    beta: f64,
    center_b: &[f64; 3],
    lb: [u32; 3],
) -> f64 {
    let p = alpha + beta;
    let mu = alpha * beta / p;
    let ab2 = distance_squared(center_a, center_b);

    let prefactor_3d = (PI / p).powf(1.5) * (-mu * ab2).exp();

    let px = gaussian_product_center(alpha, center_a[0], beta, center_b[0]);
    let py = gaussian_product_center(alpha, center_a[1], beta, center_b[1]);
    let pz = gaussian_product_center(alpha, center_a[2], beta, center_b[2]);

    let pa = [px - center_a[0], py - center_a[1], pz - center_a[2]];
    let pb = [px - center_b[0], py - center_b[1], pz - center_b[2]];

    let prefactor_1d = prefactor_3d.cbrt();

    let sx = overlap_1d(la[0], lb[0], pa[0], pb[0], p, prefactor_1d);
    let sy = overlap_1d(la[1], lb[1], pa[1], pb[1], p, prefactor_1d);
    let sz = overlap_1d(la[2], lb[2], pa[2], pb[2], p, prefactor_1d);

    sx * sy * sz
}

/// Direct 1D overlap (standalone, includes all prefactors).
pub fn overlap_1d_direct(la: u32, lb: u32, alpha: f64, ax: f64, beta: f64, bx: f64) -> f64 {
    let p = alpha + beta;
    let mu = alpha * beta / p;
    let ab = ax - bx;
    let px = gaussian_product_center(alpha, ax, beta, bx);
    let pa = px - ax;
    let pb = px - bx;

    let prefactor = (PI / p).sqrt() * (-mu * ab * ab).exp();
    overlap_1d(la, lb, pa, pb, p, prefactor)
}

/// Boys function F_n(x) for nuclear attraction integrals.
///
/// F_n(x) = ∫₀¹ t^{2n} exp(-x·t²) dt
pub fn boys_function(n: u32, x: f64) -> f64 {
    if x < 1.0e-10 {
        return 1.0 / (2.0 * n as f64 + 1.0);
    }

    if x > 30.0 + n as f64 {
        // Asymptotic: F_n(x) ≈ (2n-1)!! / 2^(n+1) · √π / x^(n+1/2)
        let df = if n == 0 {
            1.0
        } else {
            double_factorial_f64(2 * n - 1)
        };
        return df * PI.sqrt() / (2.0f64.powi(n as i32 + 1) * x.powi(n as i32) * x.sqrt());
    }

    let mut sum = 0.0;
    let mut term = 1.0 / (2.0 * n as f64 + 1.0);
    sum += term;

    for k in 1..100 {
        term *= 2.0 * x / (2.0 * (n + k) as f64 + 1.0);
        sum += term;
        if term.abs() < 1.0e-15 * sum.abs() {
            break;
        }
    }

    (-x).exp() * sum
}

/// Double factorial for f64 return: (2n-1)!! = (2n)! / (2^n · n!)
fn double_factorial_f64(n: u32) -> f64 {
    if n <= 1 {
        return 1.0;
    }
    let mut result = 1.0;
    let mut k = n;
    while k > 1 {
        result *= k as f64;
        k -= 2;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_overlap_ss_same_center() {
        let s = overlap_ss(1.0, &[0.0, 0.0, 0.0], 1.0, &[0.0, 0.0, 0.0]);
        let expected = (PI / 2.0).powf(1.5);
        assert!((s - expected).abs() < 1e-10);
    }

    #[test]
    fn test_overlap_ss_decays_with_distance() {
        let a = [0.0, 0.0, 0.0];
        let s1 = overlap_ss(1.0, &a, 1.0, &[0.0, 0.0, 0.0]);
        let s2 = overlap_ss(1.0, &a, 1.0, &[1.0, 0.0, 0.0]);
        let s3 = overlap_ss(1.0, &a, 1.0, &[2.0, 0.0, 0.0]);
        assert!(s1 > s2);
        assert!(s2 > s3);
    }

    #[test]
    fn test_overlap_primitive_ss() {
        let s = overlap_primitive(
            1.0,
            &[0.0, 0.0, 0.0],
            [0, 0, 0],
            1.0,
            &[0.0, 0.0, 0.0],
            [0, 0, 0],
        );
        let expected = (PI / 2.0).powf(1.5);
        assert!((s - expected).abs() < 1e-10);
    }

    #[test]
    fn test_boys_function_zero() {
        let f = boys_function(0, 0.0);
        assert!((f - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_boys_function_large_x() {
        let x = 50.0;
        let f = boys_function(0, x);
        let expected = 0.5 * (PI / x).sqrt();
        assert!((f - expected).abs() / expected < 0.01);
    }

    #[test]
    fn test_gaussian_product_center() {
        let c = gaussian_product_center(1.0, 0.0, 1.0, 2.0);
        assert!((c - 1.0).abs() < 1e-10);
    }
}
