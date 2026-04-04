//! D4 Dispersion Energy and Gradients — Core
//!
//! Two-body BJ-damped dispersion with optional three-body ATM term.

use super::params::{c8_from_c6, d4_coordination_number, dynamic_c6, get_d4_params};
use serde::{Deserialize, Serialize};

/// D4 configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct D4Config {
    /// s6 scaling factor (usually 1.0).
    pub s6: f64,
    /// s8 scaling factor.
    pub s8: f64,
    /// BJ damping a1.
    pub a1: f64,
    /// BJ damping a2 (Bohr).
    pub a2: f64,
    /// Whether to include three-body ATM term.
    pub three_body: bool,
    /// s9 scaling for ATM (usually 1.0).
    pub s9: f64,
}

impl Default for D4Config {
    fn default() -> Self {
        Self {
            s6: 1.0,
            s8: 0.95,
            a1: 0.45,
            a2: 4.0,
            three_body: false,
            s9: 1.0,
        }
    }
}

/// Result of D4 dispersion calculation.
///
/// Energy terms are in Hartree; `total_kcal_mol` provides the
/// convenience conversion. Two-body (`e2_body`) and three-body
/// (`e3_body`) terms are reported separately for transparency.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct D4Result {
    /// Two-body dispersion energy (Hartree).
    pub e2_body: f64,
    /// Three-body ATM energy (Hartree).
    pub e3_body: f64,
    /// Total dispersion energy (Hartree).
    pub total_energy: f64,
    /// Total in kcal/mol.
    pub total_kcal_mol: f64,
    /// Per-atom coordination numbers.
    pub coordination_numbers: Vec<f64>,
}

/// Compute D4 dispersion energy.
pub fn compute_d4_energy(elements: &[u8], positions: &[[f64; 3]], config: &D4Config) -> D4Result {
    let n = elements.len();
    let cn = d4_coordination_number(elements, positions);
    let ang_to_bohr = 1.0 / 0.529177;

    #[cfg(feature = "parallel")]
    let e2 = {
        use rayon::prelude::*;
        (0..n)
            .into_par_iter()
            .map(|i| {
                ((i + 1)..n)
                    .map(|j| pair_energy(i, j, elements, positions, &cn, config, ang_to_bohr))
                    .sum::<f64>()
            })
            .sum::<f64>()
    };

    #[cfg(not(feature = "parallel"))]
    let e2 = (0..n)
        .map(|i| {
            ((i + 1)..n)
                .map(|j| pair_energy(i, j, elements, positions, &cn, config, ang_to_bohr))
                .sum::<f64>()
        })
        .sum::<f64>();

    #[cfg(feature = "parallel")]
    let e3 = if config.three_body && n >= 3 {
        use rayon::prelude::*;
        (0..n)
            .into_par_iter()
            .map(|i| {
                let mut subtotal = 0.0;
                for j in (i + 1)..n {
                    for k in (j + 1)..n {
                        subtotal +=
                            triple_energy(i, j, k, elements, positions, &cn, config, ang_to_bohr);
                    }
                }
                subtotal
            })
            .sum::<f64>()
    } else {
        0.0
    };

    #[cfg(not(feature = "parallel"))]
    let e3 = if config.three_body && n >= 3 {
        let mut total = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                for k in (j + 1)..n {
                    total += triple_energy(i, j, k, elements, positions, &cn, config, ang_to_bohr);
                }
            }
        }
        total
    } else {
        0.0
    };

    let total = e2 + e3;
    let hartree_to_kcal = 627.509;

    D4Result {
        e2_body: e2,
        e3_body: e3,
        total_energy: total,
        total_kcal_mol: total * hartree_to_kcal,
        coordination_numbers: cn,
    }
}

fn pair_energy(
    i: usize,
    j: usize,
    elements: &[u8],
    positions: &[[f64; 3]],
    cn: &[f64],
    config: &D4Config,
    ang_to_bohr: f64,
) -> f64 {
    let dx = (positions[i][0] - positions[j][0]) * ang_to_bohr;
    let dy = (positions[i][1] - positions[j][1]) * ang_to_bohr;
    let dz = (positions[i][2] - positions[j][2]) * ang_to_bohr;
    let r = (dx * dx + dy * dy + dz * dz).sqrt();

    if r < 1e-10 {
        return 0.0;
    }

    let c6 = dynamic_c6(elements[i], elements[j], cn[i], cn[j]);
    let c8 = c8_from_c6(c6, elements[i], elements[j]);
    let r0 = if c6 > 1e-10 { (c8 / c6).sqrt() } else { 5.0 };
    let r_cut = config.a1 * r0 + config.a2;

    let r6 = r.powi(6);
    let damp6 = r6 / (r6 + r_cut.powi(6));
    let term6 = -config.s6 * c6 / r6 * damp6;

    let r8 = r.powi(8);
    let damp8 = r8 / (r8 + r_cut.powi(8));
    let term8 = -config.s8 * c8 / r8 * damp8;

    term6 + term8
}

fn triple_energy(
    i: usize,
    j: usize,
    k: usize,
    elements: &[u8],
    positions: &[[f64; 3]],
    cn: &[f64],
    config: &D4Config,
    ang_to_bohr: f64,
) -> f64 {
    let r_ab = distance_bohr(positions, i, j, ang_to_bohr);
    let r_bc = distance_bohr(positions, j, k, ang_to_bohr);
    let r_ca = distance_bohr(positions, k, i, ang_to_bohr);

    if r_ab < 1e-10 || r_bc < 1e-10 || r_ca < 1e-10 {
        return 0.0;
    }

    let c6_ab = dynamic_c6(elements[i], elements[j], cn[i], cn[j]);
    let c6_bc = dynamic_c6(elements[j], elements[k], cn[j], cn[k]);
    let c6_ca = dynamic_c6(elements[k], elements[i], cn[k], cn[i]);

    // C9 ≈ -sqrt(C6_AB * C6_BC * C6_CA) per Grimme D3 (JCP 132, 154104)
    let c9 = -(c6_ab * c6_bc * c6_ca).abs().sqrt();

    let cos_a = (r_ab * r_ab + r_ca * r_ca - r_bc * r_bc) / (2.0 * r_ab * r_ca);
    let cos_b = (r_ab * r_ab + r_bc * r_bc - r_ca * r_ca) / (2.0 * r_ab * r_bc);
    let cos_c = (r_bc * r_bc + r_ca * r_ca - r_ab * r_ab) / (2.0 * r_bc * r_ca);
    let angular = 3.0 * cos_a * cos_b * cos_c + 1.0;
    let r_prod = r_ab * r_bc * r_ca;

    // ATM BJ-like damping per Grimme D3 (JCP 132, 154104, Eq. 6):
    // f_damp = 1 / (1 + 6 * (R0_ABC / r_mean)^alpha)
    // where R0_ABC = (4/3)*(R_cov_A + R_cov_B)^(1/3) * ... geometric mean of covalent radii
    // and alpha = 14 (steep damping exponent)
    let r_cov_i = get_d4_params(elements[i]).r_cov;
    let r_cov_j = get_d4_params(elements[j]).r_cov;
    let r_cov_k = get_d4_params(elements[k]).r_cov;
    let r0_ab = (4.0 / 3.0) * (r_cov_i + r_cov_j);
    let r0_bc = (4.0 / 3.0) * (r_cov_j + r_cov_k);
    let r0_ca = (4.0 / 3.0) * (r_cov_k + r_cov_i);
    let r0_prod = r0_ab * r0_bc * r0_ca;
    let r9 = r_prod.powi(3);
    let r0_9 = r0_prod.powi(3);
    let fdamp = 1.0 / (1.0 + 6.0 * (r0_9 / r9));

    config.s9 * c9 * angular / r9 * fdamp
}

/// Compute numerical D4 gradient.
pub fn compute_d4_gradient(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &D4Config,
) -> Vec<[f64; 3]> {
    let n = elements.len();
    let h = 1e-5;
    let mut gradient = vec![[0.0; 3]; n];

    for i in 0..n {
        for d in 0..3 {
            let mut pos_p = positions.to_vec();
            let mut pos_m = positions.to_vec();
            pos_p[i][d] += h;
            pos_m[i][d] -= h;

            let ep = compute_d4_energy(elements, &pos_p, config).total_energy;
            let em = compute_d4_energy(elements, &pos_m, config).total_energy;

            gradient[i][d] = (ep - em) / (2.0 * h);
        }
    }

    gradient
}

fn distance_bohr(positions: &[[f64; 3]], i: usize, j: usize, ang_to_bohr: f64) -> f64 {
    let dx = (positions[i][0] - positions[j][0]) * ang_to_bohr;
    let dy = (positions[i][1] - positions[j][1]) * ang_to_bohr;
    let dz = (positions[i][2] - positions[j][2]) * ang_to_bohr;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_d4_energy_negative() {
        let elements = [6, 6];
        let pos = [[0.0, 0.0, 0.0], [3.5, 0.0, 0.0]];
        let config = D4Config::default();
        let result = compute_d4_energy(&elements, &pos, &config);
        assert!(
            result.total_energy < 0.0,
            "D4 energy should be negative: {}",
            result.total_energy
        );
    }

    #[test]
    fn test_d4_decays_with_distance() {
        let elements = [6, 6];
        let e_close = compute_d4_energy(
            &elements,
            &[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
            &D4Config::default(),
        )
        .total_energy;
        let e_far = compute_d4_energy(
            &elements,
            &[[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]],
            &D4Config::default(),
        )
        .total_energy;
        assert!(
            e_close.abs() > e_far.abs(),
            "D4 should decay: close={}, far={}",
            e_close,
            e_far
        );
    }

    #[test]
    fn test_d4_three_body() {
        let elements = [6, 6, 6];
        let pos = [[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [1.25, 2.17, 0.0]];
        let r2 = compute_d4_energy(
            &elements,
            &pos,
            &D4Config {
                three_body: false,
                ..Default::default()
            },
        );
        let r3 = compute_d4_energy(
            &elements,
            &pos,
            &D4Config {
                three_body: true,
                ..Default::default()
            },
        );
        assert!(
            (r3.total_energy - r2.total_energy).abs() > 0.0,
            "3-body should differ from 2-body"
        );
    }

    #[test]
    fn test_d4_gradient_finite() {
        let elements = [6, 8, 1, 1];
        let pos = [
            [0.0, 0.0, 0.0],
            [1.23, 0.0, 0.0],
            [-0.6, 0.9, 0.0],
            [-0.6, -0.9, 0.0],
        ];
        let grad = compute_d4_gradient(&elements, &pos, &D4Config::default());
        for g in &grad {
            for &d in g {
                assert!(d.is_finite(), "Gradient contains NaN/Inf");
            }
        }
    }
}
