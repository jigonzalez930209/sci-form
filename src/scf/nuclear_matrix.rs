//! Nuclear attraction matrix V_μν construction.
//!
//! V_μν = -Σ_A Z_A ∫ χ_μ(r) · (1/|r - R_A|) · χ_ν(r) d³r

use std::f64::consts::PI;

use nalgebra::DMatrix;

use super::basis::{BasisFunction, BasisSet};
use super::gaussian_integrals::{boys_function, distance_squared, gaussian_product_center};

/// Build the nuclear attraction matrix V for all nuclei.
pub fn build_nuclear_matrix(
    basis: &BasisSet,
    elements: &[u8],
    positions_bohr: &[[f64; 3]],
) -> DMatrix<f64> {
    let n = basis.n_basis;
    let mut v = DMatrix::zeros(n, n);

    for (&z, &rc) in elements.iter().zip(positions_bohr.iter()) {
        let z_eff = z as f64;

        for i in 0..n {
            for j in i..n {
                let v_ij = contracted_nuclear_attraction(
                    &basis.functions[i],
                    &basis.functions[j],
                    z_eff,
                    &rc,
                );
                v[(i, j)] -= v_ij;
                if i != j {
                    v[(j, i)] -= v_ij;
                }
            }
        }
    }

    v
}

/// Nuclear attraction integral between two contracted basis functions.
fn contracted_nuclear_attraction(
    bf_a: &BasisFunction,
    bf_b: &BasisFunction,
    z: f64,
    rc: &[f64; 3],
) -> f64 {
    let mut v = 0.0;

    for prim_a in &bf_a.primitives {
        let norm_a = BasisFunction::normalization(
            prim_a.alpha,
            bf_a.angular[0],
            bf_a.angular[1],
            bf_a.angular[2],
        );

        for prim_b in &bf_b.primitives {
            let norm_b = BasisFunction::normalization(
                prim_b.alpha,
                bf_b.angular[0],
                bf_b.angular[1],
                bf_b.angular[2],
            );

            let v_prim = nuclear_attraction_primitive(
                prim_a.alpha,
                &bf_a.center,
                bf_a.angular,
                prim_b.alpha,
                &bf_b.center,
                bf_b.angular,
                z,
                rc,
            );

            v += norm_a * prim_a.coefficient * norm_b * prim_b.coefficient * v_prim;
        }
    }

    v
}

/// Nuclear attraction integral between two Cartesian Gaussian primitives.
fn nuclear_attraction_primitive(
    alpha: f64,
    center_a: &[f64; 3],
    la: [u32; 3],
    beta: f64,
    center_b: &[f64; 3],
    lb: [u32; 3],
    z: f64,
    center_c: &[f64; 3],
) -> f64 {
    let p = alpha + beta;
    let mu = alpha * beta / p;
    let ab2 = distance_squared(center_a, center_b);

    let product_center = [
        gaussian_product_center(alpha, center_a[0], beta, center_b[0]),
        gaussian_product_center(alpha, center_a[1], beta, center_b[1]),
        gaussian_product_center(alpha, center_a[2], beta, center_b[2]),
    ];

    let pc2 = distance_squared(&product_center, center_c);
    let l_total = la[0] + la[1] + la[2] + lb[0] + lb[1] + lb[2];

    if l_total == 0 {
        let prefactor = 2.0 * PI / p * (-mu * ab2).exp();
        return z * prefactor * boys_function(0, p * pc2);
    }

    let prefactor = 2.0 * PI / p * (-mu * ab2).exp();

    let pa = [
        product_center[0] - center_a[0],
        product_center[1] - center_a[1],
        product_center[2] - center_a[2],
    ];
    let pb = [
        product_center[0] - center_b[0],
        product_center[1] - center_b[1],
        product_center[2] - center_b[2],
    ];
    let pc = [
        product_center[0] - center_c[0],
        product_center[1] - center_c[1],
        product_center[2] - center_c[2],
    ];

    let ex = hermite_coefficients(la[0], lb[0], pa[0], pb[0], p);
    let ey = hermite_coefficients(la[1], lb[1], pa[1], pb[1], p);
    let ez = hermite_coefficients(la[2], lb[2], pa[2], pb[2], p);

    // Boys values for R^n_000 = (−2p)^n F_n, n up to the total Hermite order.
    let max_n = l_total as usize;
    let arg = p * pc2;
    let mut boys = Vec::with_capacity(max_n + 1);
    let mut scale = 1.0;
    for n in 0..=max_n {
        boys.push(scale * boys_function(n as u32, arg));
        scale *= -2.0 * p;
    }

    let mut v = 0.0;
    for (t, et) in ex.iter().enumerate() {
        if et.abs() < 1e-18 {
            continue;
        }
        for (u, eu) in ey.iter().enumerate() {
            if eu.abs() < 1e-18 {
                continue;
            }
            for (vv, ev) in ez.iter().enumerate() {
                if ev.abs() < 1e-18 {
                    continue;
                }
                v += et * eu * ev * r_tuv(t, u, vv, 0, &pc, &boys);
            }
        }
    }

    z * prefactor * v
}

/// Cartesian Hermite Coulomb auxiliary R^n_{tuv}.
///
/// R^n_000 = (−2p)^n F_n is precomputed in `boys`.
/// R^{n}_{t+1,u,v} = t R^{n+1}_{t−1,u,v} + PC_x R^{n+1}_{t,u,v}.
fn r_tuv(t: usize, u: usize, v: usize, n: usize, pc: &[f64; 3], boys: &[f64]) -> f64 {
    if t == 0 && u == 0 && v == 0 {
        return boys[n];
    }
    if t > 0 {
        let lower = if t >= 2 {
            (t - 1) as f64 * r_tuv(t - 2, u, v, n + 1, pc, boys)
        } else {
            0.0
        };
        return lower + pc[0] * r_tuv(t - 1, u, v, n + 1, pc, boys);
    }
    if u > 0 {
        let lower = if u >= 2 {
            (u - 1) as f64 * r_tuv(t, u - 2, v, n + 1, pc, boys)
        } else {
            0.0
        };
        return lower + pc[1] * r_tuv(t, u - 1, v, n + 1, pc, boys);
    }
    let lower = if v >= 2 {
        (v - 1) as f64 * r_tuv(t, u, v - 2, n + 1, pc, boys)
    } else {
        0.0
    };
    lower + pc[2] * r_tuv(t, u, v - 1, n + 1, pc, boys)
}

/// Compute Hermite expansion coefficients E_t^{ij} for one dimension.
fn hermite_coefficients(la: u32, lb: u32, pa: f64, pb: f64, p: f64) -> Vec<f64> {
    let la = la as usize;
    let lb = lb as usize;
    let max_t = la + lb + 1;

    let mut e = vec![vec![vec![0.0f64; max_t + 1]; lb + 1]; la + 1];

    e[0][0][0] = 1.0;

    for a in 1..=la {
        for t in 0..=a {
            if t > 0 {
                e[a][0][t] += e[a - 1][0][t - 1] / (2.0 * p);
            }
            e[a][0][t] += pa * e[a - 1][0][t];
            if t + 1 < max_t {
                e[a][0][t] += (t + 1) as f64 * e[a - 1][0][t + 1];
            }
        }
    }

    for b in 1..=lb {
        for a in 0..=la {
            for t in 0..=(a + b) {
                if t > 0 {
                    e[a][b][t] += e[a][b - 1][t - 1] / (2.0 * p);
                }
                e[a][b][t] += pb * e[a][b - 1][t];
                if t + 1 < max_t {
                    e[a][b][t] += (t + 1) as f64 * e[a][b - 1][t + 1];
                }
            }
        }
    }

    (0..=la + lb).map(|t| e[la][lb][t]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nuclear_matrix_symmetric() {
        let elements = [8u8, 1, 1];
        let positions = [
            [0.0, 0.0, 0.2216],
            [0.0, 1.4313, -0.8864],
            [0.0, -1.4313, -0.8864],
        ];
        let basis = BasisSet::sto3g(&elements, &positions);
        let v = build_nuclear_matrix(&basis, &elements, &positions);

        for i in 0..basis.n_basis {
            for j in 0..basis.n_basis {
                assert!(
                    (v[(i, j)] - v[(j, i)]).abs() < 1e-12,
                    "V not symmetric at ({}, {})",
                    i,
                    j
                );
            }
        }
    }

    #[test]
    fn test_nuclear_matrix_negative() {
        let basis = BasisSet::sto3g(&[1], &[[0.0, 0.0, 0.0]]);
        let v = build_nuclear_matrix(&basis, &[1], &[[0.0, 0.0, 0.0]]);
        assert!(v[(0, 0)] < 0.0, "V_11 should be negative for attraction");
    }

    #[test]
    fn test_hermite_base_case() {
        let e = hermite_coefficients(0, 0, 0.0, 0.0, 1.0);
        assert!((e[0] - 1.0).abs() < 1e-14);
    }
}
