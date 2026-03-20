//! Electric dipole transition moments.
//!
//! Computes transition dipole moments (TDMs) between electronic states
//! using molecular orbital coefficients and dipole integrals:
//!
//!   μ_{0→n} = Σ_{ia} X_{ia}^n · <φ_i|r|φ_a>
//!
//! where X_{ia}^n are CI coefficients for excitation n, and <φ_i|r|φ_a>
//! are MO-basis dipole integrals.
//!
//! The oscillator strength is:
//!   f_n = (2/3) · ΔE_n · |μ_{0→n}|²

use nalgebra::DMatrix;

use crate::experimental_2::phase2_quantum_engine::basis_set::{BasisFunction, BasisSet};

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
///
/// For s-type Gaussians:
///   <s_A|x|s_B> = P_x · S_AB + (1/(2(α+β))) · ∂S/∂A_x
///
/// where P = (α·A + β·B)/(α+β) is the Gaussian product center
/// and S_AB is the overlap integral.
///
/// For higher angular momentum, use the recursion:
///   <a+1_x|r|b> = P_x · <a|r|b> + (a/(2p)) · <a-1|r|b> + (b/(2p)) · <a|r|b-1>
///   <a|x|b> = <a+1_x|b> + A_x · <a|b>  (for x-component dipole)
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

            // Contract over primitives
            for prim_a in &bf_mu.primitives {
                for prim_b in &bf_nu.primitives {
                    let alpha = prim_a.alpha;
                    let beta = prim_b.alpha;
                    let p = alpha + beta;
                    let mu_ab = alpha * beta / p;

                    let a = bf_mu.center;
                    let b = bf_nu.center;

                    // Gaussian product center
                    let px = (alpha * a[0] + beta * b[0]) / p;
                    let py = (alpha * a[1] + beta * b[1]) / p;
                    let pz = (alpha * a[2] + beta * b[2]) / p;

                    // Distance squared
                    let ab2 = (a[0] - b[0]).powi(2)
                        + (a[1] - b[1]).powi(2)
                        + (a[2] - b[2]).powi(2);

                    // Overlap prefactor for s-type
                    let s_prefactor = (std::f64::consts::PI / p).powf(1.5)
                        * (-mu_ab * ab2).exp()
                        * prim_a.coefficient
                        * prim_b.coefficient;

                    // For s-type functions (l=0):
                    // <s_A|r_k|s_B> = P_k · S_AB
                    let norm_a = BasisFunction::normalization(prim_a.alpha, bf_mu.angular[0], bf_mu.angular[1], bf_mu.angular[2]);
                    let norm_b = BasisFunction::normalization(prim_b.alpha, bf_nu.angular[0], bf_nu.angular[1], bf_nu.angular[2]);
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
    // D_mo = C^T · D_ao · C
    let x_mo = mo_coefficients.transpose() * &dipole_ao.x * mo_coefficients;
    let y_mo = mo_coefficients.transpose() * &dipole_ao.y * mo_coefficients;
    let z_mo = mo_coefficients.transpose() * &dipole_ao.z * mo_coefficients;

    DipoleIntegrals {
        x: x_mo,
        y: y_mo,
        z: z_mo,
    }
}

/// Compute oscillator strength for a transition.
///
///   f = (2/3) · ΔE · (|μ_x|² + |μ_y|² + |μ_z|²)
///
/// where ΔE is in Hartree and μ is in atomic units (e·a₀).
pub fn oscillator_strength(transition_dipole: &[f64; 3], energy_hartree: f64) -> f64 {
    let mu_sq = transition_dipole[0].powi(2)
        + transition_dipole[1].powi(2)
        + transition_dipole[2].powi(2);

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
    let n_occ = occ_indices.len();
    let n_virt = virt_indices.len();
    let mut tdm = [0.0; 3];

    for i_l in 0..n_occ {
        for a_l in 0..n_virt {
            let idx = i_l * n_virt + a_l;
            let x = ci_vector[idx];

            if x.abs() < 1e-12 {
                continue;
            }

            let i = occ_indices[i_l];
            let a = virt_indices[a_l];

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
        // Test with known values
        let tdm = [1.0, 0.0, 0.0];
        let energy = 0.1; // Hartree
        let f = oscillator_strength(&tdm, energy);

        // f = (2/3) * 0.1 * 1.0 = 0.0667
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

        // With identity MO coefficients, AO = MO
        for i in 0..n {
            for j in 0..n {
                assert!((d_mo.x[(i, j)] - d.x[(i, j)]).abs() < 1e-12);
            }
        }
    }
}
