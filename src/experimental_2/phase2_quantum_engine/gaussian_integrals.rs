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

use crate::experimental_2::constants::PI;

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
///
/// Computes the entire (la+1) × (lb+1) table and returns element [la][lb].
pub fn overlap_1d(la: u32, lb: u32, pa_x: f64, pb_x: f64, p: f64, prefactor: f64) -> f64 {
    let max_a = (la + 1) as usize;
    let max_b = (lb + 1) as usize;
    let mut table = vec![vec![0.0f64; max_b]; max_a];

    // Base case: S(0,0) = prefactor (includes exponential)
    table[0][0] = prefactor;

    // Build upward in 'a' for b=0
    if max_a > 1 {
        table[1][0] = pa_x * table[0][0];
    }
    for a in 2..max_a {
        table[a][0] = pa_x * table[a - 1][0]
            + (a - 1) as f64 / (2.0 * p) * table[a - 2][0];
    }

    // Build upward in 'b' using transfer relation
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
///
/// S = Nx·Ny·Nz · S_x(lax, lbx) · S_y(lay, lby) · S_z(laz, lbz)
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

    // Product center
    let px = gaussian_product_center(alpha, center_a[0], beta, center_b[0]);
    let py = gaussian_product_center(alpha, center_a[1], beta, center_b[1]);
    let pz = gaussian_product_center(alpha, center_a[2], beta, center_b[2]);

    // PA and PB vectors
    let pa = [px - center_a[0], py - center_a[1], pz - center_a[2]];
    let pb = [px - center_b[0], py - center_b[1], pz - center_b[2]];

    // 1D prefactors (split from 3D)
    let prefactor_1d = prefactor_3d.cbrt();

    let sx = overlap_1d(la[0], lb[0], pa[0], pb[0], p, prefactor_1d);
    let sy = overlap_1d(la[1], lb[1], pa[1], pb[1], p, prefactor_1d);
    let sz = overlap_1d(la[2], lb[2], pa[2], pb[2], p, prefactor_1d);

    sx * sy * sz
}

/// Kinetic energy integral between two s-type primitives.
///
/// T = β(3 - 2β|PB|²) · S - 2β² · S(0,2)_corrected
///
/// For general angular momentum, uses Obara-Saika kinetic recurrence:
///   T(a,b) = -0.5·b(b-1)·S(a,b-2) + β(2b+1)·S(a,b) - 2β²·S(a,b+2)
///            applied to each Cartesian component.
pub fn kinetic_primitive(
    alpha: f64,
    center_a: &[f64; 3],
    la: [u32; 3],
    beta: f64,
    center_b: &[f64; 3],
    lb: [u32; 3],
) -> f64 {
    let mut kinetic = 0.0;

    // T = Tx·Sy·Sz + Sx·Ty·Sz + Sx·Sy·Tz
    // where Tx = kinetic integral in x, Sx = overlap in x, etc.
    for dim in 0..3 {
        let la_d = la[dim];
        let lb_d = lb[dim];

        // Kinetic component in this dimension:
        // T_d = 0.5·b(b-1)·S(a, b-2) - β(2b+1)·S(a,b) + 2β²·S(a, b+2)
        // but in the Obara-Saika formulation for the kinetic integral:
        // T(a,b) = β(2b+1)S(a,b) - 2β²S(a,b+2) at b=0:
        // T(a,0) = βS(a,0) - 2β²S(a,2) ... generalized

        let t_component = kinetic_1d(alpha, center_a, la_d, beta, center_b, lb_d, dim);

        // Overlap in the other two dimensions
        let mut s_other = 1.0;
        for dim2 in 0..3 {
            if dim2 != dim {
                s_other *= overlap_1d_direct(
                    la[dim2], lb[dim2],
                    alpha, center_a[dim2],
                    beta, center_b[dim2],
                );
            }
        }

        kinetic += t_component * s_other;
    }

    kinetic
}

/// 1D kinetic energy integral component.
fn kinetic_1d(
    alpha: f64,
    center_a: &[f64; 3],
    la: u32,
    beta: f64,
    center_b: &[f64; 3],
    lb: u32,
    dim: usize,
) -> f64 {
    let p = alpha + beta;
    let mu = alpha * beta / p;
    let ab = center_a[dim] - center_b[dim];

    // Use the relation: T = β(2lb+1)S(la,lb) - 2β²S(la,lb+2) + 0.5lb(lb-1)S(la,lb-2)
    // But we need the 1D overlap helper

    let s_direct = |a: u32, b: u32| -> f64 {
        overlap_1d_direct(a, b, alpha, center_a[dim], beta, center_b[dim])
    };

    let mut t = beta * (2.0 * lb as f64 + 1.0) * s_direct(la, lb);
    t -= 2.0 * beta * beta * s_direct(la, lb + 2);
    if lb >= 2 {
        t += 0.5 * lb as f64 * (lb - 1) as f64 * s_direct(la, lb - 2);
    }

    // Include the normalization from the 3D exponential split
    let exp_factor = (-mu * ab * ab).exp();
    t * (PI / p).sqrt() * exp_factor / (PI / p).sqrt() // normalize properly

    // Simplified: the overlap_1d_direct already handles prefactor
}

/// Direct 1D overlap (standalone, includes all prefactors).
fn overlap_1d_direct(la: u32, lb: u32, alpha: f64, ax: f64, beta: f64, bx: f64) -> f64 {
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
///
/// Uses series expansion for small x and asymptotic formula for large x.
pub fn boys_function(n: u32, x: f64) -> f64 {
    if x < 1.0e-10 {
        return 1.0 / (2.0 * n as f64 + 1.0);
    }

    if x > 30.0 + n as f64 {
        // Asymptotic: F_n(x) ≈ (2n-1)!! / 2^{n+1} · sqrt(π/x^{2n+1})
        let numerator = double_factorial_f64(2 * n);
        let denominator = 2.0f64.powi(n as i32 + 1) * x.powi(n as i32);
        return numerator * (PI / x).sqrt() / (2.0 * denominator);
    }

    // Upward recursion from F_0(x) = sqrt(π/x)·erf(sqrt(x))/2
    // OR downward recursion from high n using series expansion

    // Series expansion: F_n(x) = exp(-x) · Σ_{k=0}^{∞} (2x)^k / (2n+2k+1)!!
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
        // Two identical s-type Gaussians at same center: S = (π/2α)^{3/2}
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
            1.0, &[0.0, 0.0, 0.0], [0, 0, 0],
            1.0, &[0.0, 0.0, 0.0], [0, 0, 0],
        );
        let expected = (PI / 2.0).powf(1.5);
        assert!((s - expected).abs() < 1e-10);
    }

    #[test]
    fn test_boys_function_zero() {
        // F_0(0) = 1
        let f = boys_function(0, 0.0);
        assert!((f - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_boys_function_large_x() {
        // F_0(x) → sqrt(π/x)/2 for large x
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
