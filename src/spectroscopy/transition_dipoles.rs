//! Electric dipole transition moments.
//!
//! Computes transition dipole moments (TDMs) between electronic states
//! using molecular orbital coefficients and dipole integrals:
//!
//!   μ_{0→n} = Σ_{ia} X_{ia}^n · <φ_i|r|φ_a>
//!
//! The oscillator strength is:
//!   f_n = (2/3) · ΔE_n · |μ_{0→n}|²

use nalgebra::DMatrix;

use crate::scf::basis::{BasisFunction, BasisSet};

/// Dipole integral result.
#[derive(Debug, Clone)]
pub struct DipoleIntegrals {
    /// x-component dipole matrix in AO basis.
    pub x: DMatrix<f64>,
    /// y-component dipole matrix in AO basis.
    pub y: DMatrix<f64>,
    /// z-component dipole matrix in AO basis.
    pub z: DMatrix<f64>,
}

/// Compute AO-basis dipole integrals <μ|r_k|ν> for k = x, y, z.
pub fn compute_dipole_integrals(basis: &BasisSet) -> DipoleIntegrals {
    let n = basis.n_basis;
    let mut dip_x = DMatrix::zeros(n, n);
    let mut dip_y = DMatrix::zeros(n, n);
    let mut dip_z = DMatrix::zeros(n, n);

    for mu in 0..n {
        for nu in mu..n {
            let bf_mu = &basis.functions[mu];
            let bf_nu = &basis.functions[nu];

            let mut dx = 0.0;
            let mut dy = 0.0;
            let mut dz = 0.0;

            for prim_a in &bf_mu.primitives {
                for prim_b in &bf_nu.primitives {
                    let alpha = prim_a.alpha;
                    let beta = prim_b.alpha;
                    let p = alpha + beta;
                    let mu_ab = alpha * beta / p;

                    let a = bf_mu.center;
                    let b = bf_nu.center;

                    let px = (alpha * a[0] + beta * b[0]) / p;
                    let py = (alpha * a[1] + beta * b[1]) / p;
                    let pz = (alpha * a[2] + beta * b[2]) / p;

                    let ab2 = (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2);

                    let s_prefactor = (std::f64::consts::PI / p).powf(1.5)
                        * (-mu_ab * ab2).exp()
                        * prim_a.coefficient
                        * prim_b.coefficient;

                    let norm_a = BasisFunction::normalization(
                        prim_a.alpha,
                        bf_mu.angular[0],
                        bf_mu.angular[1],
                        bf_mu.angular[2],
                    );
                    let norm_b = BasisFunction::normalization(
                        prim_b.alpha,
                        bf_nu.angular[0],
                        bf_nu.angular[1],
                        bf_nu.angular[2],
                    );
                    let s = s_prefactor * norm_a * norm_b;

                    dx += px * s;
                    dy += py * s;
                    dz += pz * s;
                }
            }

            dip_x[(mu, nu)] = dx;
            dip_x[(nu, mu)] = dx;
            dip_y[(mu, nu)] = dy;
            dip_y[(nu, mu)] = dy;
            dip_z[(mu, nu)] = dz;
            dip_z[(nu, mu)] = dz;
        }
    }

    DipoleIntegrals {
        x: dip_x,
        y: dip_y,
        z: dip_z,
    }
}

/// Transform dipole integrals from AO to MO basis.
///
///   <i|r|a> = Σ_μν C_μi · <μ|r|ν> · C_νa
pub fn transform_to_mo(
    dipole_ao: &DipoleIntegrals,
    mo_coefficients: &DMatrix<f64>,
) -> DipoleIntegrals {
    DipoleIntegrals {
        x: mo_coefficients.transpose() * &dipole_ao.x * mo_coefficients,
        y: mo_coefficients.transpose() * &dipole_ao.y * mo_coefficients,
        z: mo_coefficients.transpose() * &dipole_ao.z * mo_coefficients,
    }
}

/// Compute oscillator strength for a transition.
///
///   f = (2/3) · ΔE · (|μ_x|² + |μ_y|² + |μ_z|²)
pub fn oscillator_strength(transition_dipole: &[f64; 3], energy_hartree: f64) -> f64 {
    let mu_sq =
        transition_dipole[0].powi(2) + transition_dipole[1].powi(2) + transition_dipole[2].powi(2);

    (2.0 / 3.0) * energy_hartree * mu_sq
}

/// Compute transition dipole moment from CI vector and MO dipole integrals.
///
///   μ_{0→n} = Σ_{ia} X_{ia}^n · <i|r|a>
pub fn transition_dipole_from_ci_vector(
    ci_vector: &[f64],
    occ_indices: &[usize],
    virt_indices: &[usize],
    dipole_mo: &DipoleIntegrals,
) -> [f64; 3] {
    let n_virt = virt_indices.len();
    let mut tdm = [0.0; 3];

    for (i_l, &i) in occ_indices.iter().enumerate() {
        for (a_l, &a) in virt_indices.iter().enumerate() {
            let idx = i_l * n_virt + a_l;
            let x = ci_vector[idx];

            if x.abs() < 1e-12 {
                continue;
            }

            tdm[0] += x * dipole_mo.x[(i, a)];
            tdm[1] += x * dipole_mo.y[(i, a)];
            tdm[2] += x * dipole_mo.z[(i, a)];
        }
    }

    tdm
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_oscillator_strength_calculation() {
        let tdm = [1.0, 0.0, 0.0];
        let energy = 0.1;
        let f = oscillator_strength(&tdm, energy);
        assert!((f - 2.0 / 30.0).abs() < 1e-10);
    }

    #[test]
    fn test_oscillator_strength_zero_tdm() {
        let tdm = [0.0, 0.0, 0.0];
        let f = oscillator_strength(&tdm, 1.0);
        assert!(f.abs() < 1e-15);
    }

    #[test]
    fn test_mo_transform_identity() {
        let n = 3;
        let d = DipoleIntegrals {
            x: DMatrix::from_element(n, n, 1.0),
            y: DMatrix::from_element(n, n, 2.0),
            z: DMatrix::from_element(n, n, 3.0),
        };

        let c = DMatrix::identity(n, n);
        let d_mo = transform_to_mo(&d, &c);

        for i in 0..n {
            for j in 0..n {
                assert!((d_mo.x[(i, j)] - d.x[(i, j)]).abs() < 1e-12);
            }
        }
    }
}
