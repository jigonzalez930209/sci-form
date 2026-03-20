//! D4 Dispersion Energy and Gradients — E7.2
//!
//! Two-body BJ-damped dispersion with optional three-body ATM term.

use super::params::{c8_from_c6, d4_coordination_number, dynamic_c6};

/// D4 configuration.
#[derive(Debug, Clone)]
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
        // EHT-D4 default parameters
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
#[derive(Debug, Clone)]
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
pub fn compute_d4_energy(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &D4Config,
) -> D4Result {
    let n = elements.len();
    let cn = d4_coordination_number(elements, positions);
    let ang_to_bohr = 1.0 / 0.529177;

    let mut e2 = 0.0;

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = (positions[i][0] - positions[j][0]) * ang_to_bohr;
            let dy = (positions[i][1] - positions[j][1]) * ang_to_bohr;
            let dz = (positions[i][2] - positions[j][2]) * ang_to_bohr;
            let r = (dx * dx + dy * dy + dz * dz).sqrt();

            if r < 1e-10 { continue; }

            let c6 = dynamic_c6(elements[i], elements[j], cn[i], cn[j]);
            let c8 = c8_from_c6(c6, elements[i], elements[j]);

            // BJ damping radii
            let r0 = if c6 > 1e-10 { (c8 / c6).sqrt() } else { 5.0 };
            let r_cut6 = config.a1 * r0 + config.a2;
            let r_cut8 = config.a1 * r0 + config.a2;

            // Two-body E6 term
            let r6 = r.powi(6);
            let damp6 = r6 / (r6 + r_cut6.powi(6));
            e2 -= config.s6 * c6 / r6 * damp6;

            // Two-body E8 term
            let r8 = r.powi(8);
            let damp8 = r8 / (r8 + r_cut8.powi(8));
            e2 -= config.s8 * c8 / r8 * damp8;
        }
    }

    // Three-body ATM
    let mut e3 = 0.0;
    if config.three_body && n >= 3 {
        for i in 0..n {
            for j in (i + 1)..n {
                for k in (j + 1)..n {
                    let r_ab = distance_bohr(positions, i, j, ang_to_bohr);
                    let r_bc = distance_bohr(positions, j, k, ang_to_bohr);
                    let r_ca = distance_bohr(positions, k, i, ang_to_bohr);

                    if r_ab < 1e-10 || r_bc < 1e-10 || r_ca < 1e-10 { continue; }

                    let c6_ab = dynamic_c6(elements[i], elements[j], cn[i], cn[j]);
                    let c6_bc = dynamic_c6(elements[j], elements[k], cn[j], cn[k]);
                    let c6_ca = dynamic_c6(elements[k], elements[i], cn[k], cn[i]);

                    let c9 = -(c6_ab * c6_bc * c6_ca).abs().cbrt().powi(3);

                    // Angles from law of cosines
                    let cos_a = (r_ab * r_ab + r_ca * r_ca - r_bc * r_bc) / (2.0 * r_ab * r_ca);
                    let cos_b = (r_ab * r_ab + r_bc * r_bc - r_ca * r_ca) / (2.0 * r_ab * r_bc);
                    let cos_c = (r_bc * r_bc + r_ca * r_ca - r_ab * r_ab) / (2.0 * r_bc * r_ca);

                    let angular = 3.0 * cos_a * cos_b * cos_c + 1.0;
                    let r_prod = r_ab * r_bc * r_ca;

                    e3 += config.s9 * c9 * angular / r_prod.powi(3);
                }
            }
        }
    }

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
        // Dispersion should be attractive (negative)
        let elements = [6, 6];
        let pos = [[0.0, 0.0, 0.0], [3.5, 0.0, 0.0]]; // ~3.5 Å apart
        let config = D4Config::default();
        let result = compute_d4_energy(&elements, &pos, &config);
        assert!(result.total_energy < 0.0, "D4 energy should be negative: {}", result.total_energy);
    }

    #[test]
    fn test_d4_decays_with_distance() {
        let elements = [6, 6];
        let e_close = compute_d4_energy(&elements,
            &[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]], &D4Config::default()).total_energy;
        let e_far = compute_d4_energy(&elements,
            &[[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]], &D4Config::default()).total_energy;
        assert!(e_close.abs() > e_far.abs(),
            "D4 should decay: close={}, far={}", e_close, e_far);
    }

    #[test]
    fn test_d4_three_body() {
        let elements = [6, 6, 6];
        let pos = [[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [1.25, 2.17, 0.0]];
        let r2 = compute_d4_energy(&elements, &pos, &D4Config { three_body: false, ..Default::default() });
        let r3 = compute_d4_energy(&elements, &pos, &D4Config { three_body: true, ..Default::default() });
        assert!((r3.total_energy - r2.total_energy).abs() > 0.0, "3-body should differ from 2-body");
    }

    #[test]
    fn test_d4_gradient_finite() {
        let elements = [6, 8, 1, 1];
        let pos = [
            [0.0, 0.0, 0.0], [1.23, 0.0, 0.0],
            [-0.6, 0.9, 0.0], [-0.6, -0.9, 0.0],
        ];
        let grad = compute_d4_gradient(&elements, &pos, &D4Config::default());
        for g in &grad {
            for &d in g {
                assert!(d.is_finite(), "Gradient contains NaN/Inf");
            }
        }
    }
}
