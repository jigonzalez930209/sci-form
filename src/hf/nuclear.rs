//! Nuclear attraction integral evaluator.
//!
//! Computes the matrix elements $V_{\mu\nu} = -\sum_A Z_A (\mu|\frac{1}{r_A}|\nu)$
//! where the sum runs over all nuclei A with charge $Z_A$.
//! Uses the Boys function for the incomplete gamma integral.

use super::basis::{BasisSet, ShellType};
use nalgebra::DMatrix;

/// Compute the nuclear attraction matrix V.
pub fn compute_nuclear_matrix(
    basis: &BasisSet,
    elements: &[u8],
    positions_bohr: &[[f64; 3]],
) -> DMatrix<f64> {
    let n = basis.n_basis();
    let mut v_mat = DMatrix::zeros(n, n);
    let mut mu = 0;

    for shell_i in &basis.shells {
        let ni = shell_i.n_functions();
        let mut nu = 0;

        for shell_j in &basis.shells {
            let nj = shell_j.n_functions();
            let block = nuclear_shell_pair(shell_i, shell_j, elements, positions_bohr);

            for i in 0..ni {
                for j in 0..nj {
                    v_mat[(mu + i, nu + j)] = block[i * nj + j];
                }
            }
            nu += nj;
        }
        mu += ni;
    }
    v_mat
}

fn nuclear_shell_pair(
    a: &super::basis::Shell,
    b: &super::basis::Shell,
    elements: &[u8],
    atoms_bohr: &[[f64; 3]],
) -> Vec<f64> {
    let na = a.n_functions();
    let nb = b.n_functions();
    let mut result = vec![0.0; na * nb];
    let ab2 = dist_sq(&a.center, &b.center);

    for (&ea, &ca) in a.exponents.iter().zip(&a.coefficients) {
        for (&eb, &cb) in b.exponents.iter().zip(&b.coefficients) {
            let gamma = ea + eb;
            let p = [
                (ea * a.center[0] + eb * b.center[0]) / gamma,
                (ea * a.center[1] + eb * b.center[1]) / gamma,
                (ea * a.center[2] + eb * b.center[2]) / gamma,
            ];

            let exp_ab = (-ea * eb / gamma * ab2).exp();

            for (atom_idx, &z) in elements.iter().enumerate() {
                let c = &atoms_bohr[atom_idx];
                let pc2 = dist_sq(&p, c);
                let t = gamma * pc2;
                let f0 = boys_function(0, t);

                let prefactor = -2.0 * std::f64::consts::PI / gamma * exp_ab * z as f64 * ca * cb;

                match (a.shell_type, b.shell_type) {
                    (ShellType::S, ShellType::S) => {
                        result[0] += prefactor * f0;
                    }
                    (ShellType::S, ShellType::P) => {
                        let f1 = boys_function(1, t);
                        for d in 0..3 {
                            let pb = p[d] - b.center[d];
                            let pc = p[d] - c[d];
                            result[d] += prefactor * (pb * f0 - pc * f1);
                        }
                    }
                    (ShellType::P, ShellType::S) => {
                        let f1 = boys_function(1, t);
                        for d in 0..3 {
                            let pa = p[d] - a.center[d];
                            let pc = p[d] - c[d];
                            result[d * nb] += prefactor * (pa * f0 - pc * f1);
                        }
                    }
                    (ShellType::P, ShellType::P) => {
                        let f1 = boys_function(1, t);
                        let f2 = boys_function(2, t);
                        for i in 0..3 {
                            for j in 0..3 {
                                let pa = p[i] - a.center[i];
                                let pb = p[j] - b.center[j];
                                let pc_i = p[i] - c[i];
                                let pc_j = p[j] - c[j];
                                let mut val = pa * pb * f0 - pa * pc_j * f1 - pb * pc_i * f1
                                    + pc_i * pc_j * f2;
                                if i == j {
                                    val += (f0 - f1) / (2.0 * gamma);
                                }
                                result[i * nb + j] += prefactor * val;
                            }
                        }
                    }
                }
            }
        }
    }
    result
}

/// Boys function $F_n(t) = \int_0^1 u^{2n} e^{-tu^2} du$.
///
/// Uses series expansion for small t and asymptotic form for large t.
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

#[inline]
fn dist_sq(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
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
