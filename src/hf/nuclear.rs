//! Nuclear attraction integral evaluator.
//!
//! Computes the matrix elements $V_{\mu\nu} = -\sum_A Z_A (\mu|\frac{1}{r_A}|\nu)$
//! where the sum runs over all nuclei A with charge $Z_A$.
//! Uses the Boys function for the incomplete gamma integral.

use super::basis::BasisSet;
use nalgebra::DMatrix;

/// Compute the nuclear attraction matrix V.
pub fn compute_nuclear_matrix(
    basis: &BasisSet,
    elements: &[u8],
    positions_bohr: &[[f64; 3]],
) -> DMatrix<f64> {
    // Use the cross-validated `scf` engine (fichas F36/F37).
    let scf_basis = super::basis::to_scf_basis(basis);
    crate::scf::nuclear_matrix::build_nuclear_matrix(&scf_basis, elements, positions_bohr)
}

pub fn boys_function(n: usize, t: f64) -> f64 {
    if t < 1e-15 {
        return 1.0 / (2.0 * n as f64 + 1.0);
    }
    if t > 30.0 {
        // Asymptotic: F_n(t) ≈ (2n-1)!! / 2^(n+1) * sqrt(π/t^(2n+1))
        let _nf = n as f64;
        return double_factorial(2 * n) as f64 / 2.0f64.powi(n as i32 + 1)
            * (std::f64::consts::PI / t.powi(2 * n as i32 + 1)).sqrt();
    }

    // Series expansion: F_n(t) = e^{-t} Σ_{k=0}^∞ t^k / (2n+2k+1)!!
    let mut sum = 0.0;
    let mut term = 1.0 / (2 * n + 1) as f64;
    sum += term;
    for k in 1..100 {
        term *= t / (n as f64 + k as f64 + 0.5);
        sum += term;
        if term.abs() < 1e-15 * sum.abs() {
            break;
        }
    }
    sum * (-t).exp()
}

fn double_factorial(n: usize) -> u64 {
    if n <= 1 {
        return 1;
    }
    let mut result = 1u64;
    let mut k = n;
    while k > 1 {
        result *= k as u64;
        k -= 2;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_boys_f0_zero() {
        let f = boys_function(0, 0.0);
        assert!((f - 1.0).abs() < 1e-10, "F_0(0) should be 1.0, got {f}");
    }

    #[test]
    fn test_boys_f0_large() {
        let f = boys_function(0, 50.0);
        let expected = (std::f64::consts::PI / 50.0).sqrt() * 0.5;
        assert!(
            (f - expected).abs() < 0.01 * expected,
            "F_0(50) = {f}, expected ≈ {expected}"
        );
    }
}
