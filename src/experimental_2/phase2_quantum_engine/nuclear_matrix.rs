//! Nuclear attraction matrix V_μν construction.
//!
//! V_μν = -Σ_A Z_A ∫ χ_μ(r) · (1/|r - R_A|) · χ_ν(r) d³r
//!
//! Uses the McMurchie-Davidson scheme with Boys function for the
//! Coulomb potential. The nuclear attraction integral for two
//! s-type primitives is:
//!
//!   V_A = -Z_A · (2π/p) · exp(-μ|AB|²) · F_0(p|PC|²)
//!
//! where P is the Gaussian product center and C = R_A is the nucleus.

use nalgebra::DMatrix;

use super::basis_set::{BasisFunction, BasisSet};
use super::gaussian_integrals::{
    boys_function, distance_squared, gaussian_product_center,
};
use crate::experimental_2::constants::PI;

/// Build the nuclear attraction matrix V for all nuclei.
///
/// V_μν = Σ_A V_μν^A where A runs over all atoms.
pub fn build_nuclear_matrix(basis: &BasisSet, elements: &[u8], positions_bohr: &[[f64; 3]]) -> DMatrix<f64> {
    let n = basis.n_basis;
    let mut v = DMatrix::zeros(n, n);

    for (_atom_idx, (&z, &rc)) in elements.iter().zip(positions_bohr.iter()).enumerate() {
        let z_eff = z as f64; // Nuclear charge

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

/// Nuclear attraction integral between two contracted basis functions
/// for a single nucleus at position `rc` with charge `z`.
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

/// Nuclear attraction integral between two Cartesian Gaussian primitives
/// for a nucleus at position C with charge Z.
///
/// For s-type (l=0) primitives:
///   V = Z · (2π/p) · exp(-μ|AB|²) · F_0(p|PC|²)
///
/// For higher angular momentum, uses the McMurchie-Davidson expansion
/// with Hermite Gaussians and the auxiliary R integrals.
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

    // Gaussian product center P
    let product_center = [
        gaussian_product_center(alpha, center_a[0], beta, center_b[0]),
        gaussian_product_center(alpha, center_a[1], beta, center_b[1]),
        gaussian_product_center(alpha, center_a[2], beta, center_b[2]),
    ];

    let pc2 = distance_squared(&product_center, center_c);
    let l_total = la[0] + la[1] + la[2] + lb[0] + lb[1] + lb[2];

    // For s-s pairs (most common in STO-3G)
    if l_total == 0 {
        let prefactor = 2.0 * PI / p * (-mu * ab2).exp();
        return z * prefactor * boys_function(0, p * pc2);
    }

    // For higher angular momentum: McMurchie-Davidson scheme
    // Expand in Hermite Gaussians and use R auxiliary integrals
    //
    // V = Z · (2π/p) · exp(-μ|AB|²) · Σ_{t,u,v} E_t^{la,lb,x} E_u^{la,lb,y} E_v^{la,lb,z}
    //     × R_{t+u+v}(p, PC)
    //
    // where E coefficients are Hermite expansion coefficients

    let prefactor = 2.0 * PI / p * (-mu * ab2).exp();

    // PA and PB vectors
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

    // Compute Hermite expansion coefficients E for each dimension
    let ex = hermite_coefficients(la[0], lb[0], pa[0], pb[0], p);
    let ey = hermite_coefficients(la[1], lb[1], pa[1], pb[1], p);
    let ez = hermite_coefficients(la[2], lb[2], pa[2], pb[2], p);

    // Compute R auxiliary integrals
    let max_n = l_total as usize;
    let r_aux = r_auxiliary(max_n, p, &pc);

    // Contract: V = Σ_{t,u,v} E_t^x E_u^y E_v^z R_{t+u+v,0,0,0}
    let mut v = 0.0;
    for (t, et) in ex.iter().enumerate() {
        for (u, eu) in ey.iter().enumerate() {
            for (vv, ev) in ez.iter().enumerate() {
                let n = t + u + vv;
                if n <= max_n {
                    v += et * eu * ev * r_aux[n];
                }
            }
        }
    }

    z * prefactor * v
}

/// Compute Hermite expansion coefficients E_t^{ij} for one dimension.
///
/// E_t = coefficients that expand a product of two Cartesian Gaussians
/// as a sum of Hermite Gaussians:
///
///   (x-A)^i (x-B)^j exp(-p(x-P)²) = Σ_t E_t Λ_t(p, x-P)
///
/// Recursion:
///   E_0^{00} = 1
///   E_t^{i+1,j} = (1/2p)E_{t-1}^{ij} + PA·E_t^{ij} + (t+1)E_{t+1}^{ij}
fn hermite_coefficients(la: u32, lb: u32, pa: f64, pb: f64, p: f64) -> Vec<f64> {
    let la = la as usize;
    let lb = lb as usize;
    let max_t = la + lb + 1;

    // e[a][b][t]
    let mut e = vec![vec![vec![0.0f64; max_t + 1]; lb + 1]; la + 1];

    // Base case
    e[0][0][0] = 1.0;

    // Build upward in a (b=0)
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

    // Build upward in b
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

    // Extract E_t for t = 0..la+lb
    (0..=la + lb).map(|t| e[la][lb][t]).collect()
}

/// Compute R auxiliary integrals R_n(p, PC).
///
/// R_n = (-2p)^n F_n(p|PC|²)
///
/// These are the Boys-function-weighted terms needed for nuclear attraction.
fn r_auxiliary(max_n: usize, p: f64, pc: &[f64; 3]) -> Vec<f64> {
    let pc2 = pc[0] * pc[0] + pc[1] * pc[1] + pc[2] * pc[2];
    let arg = p * pc2;

    (0..=max_n)
        .map(|n| (-2.0 * p).powi(n as i32) * boys_function(n as u32, arg))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::experimental_2::phase2_quantum_engine::basis_set::BasisSet;

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
                    i, j
                );
            }
        }
    }

    #[test]
    fn test_nuclear_matrix_negative() {
        // Nuclear attraction should make diagonal elements negative
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
