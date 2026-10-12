//! D4 Dispersion Energy and Gradients — Core
//!
//! Two-body BJ-damped dispersion with optional three-body ATM term.

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
    use crate::charges_eeq::compute_eeq2019_charges;
    use crate::xtb::d4::{D4Damping, D4Model, D4RefSet};

    // Standalone D4 routes through the cross-validated `xtb::d4::D4Model`
    // engine (validation F46): CN-dependent Gaussian reference weights with
    // charge-dependent zeta scaling, Casimir-Polder C6, BJ damping and ATM.
    // `positions` are in Angstrom; the engine works in Bohr.
    const ANG_TO_BOHR: f64 = 1.0 / 0.529177210903;
    let pos_bohr: Vec<[f64; 3]> = positions
        .iter()
        .map(|p| [p[0] * ANG_TO_BOHR, p[1] * ANG_TO_BOHR, p[2] * ANG_TO_BOHR])
        .collect();

    let damping = D4Damping {
        s6: config.s6,
        s8: config.s8,
        a1: config.a1,
        a2: config.a2,
        s9: config.s9,
    };
    let model = D4Model::with_damping_and_refs(elements, &pos_bohr, damping, D4RefSet::Dft);
    // EEQ2019 charges (exact dftd4 port, F46); positions in Bohr.
    let charges = compute_eeq2019_charges(elements, &pos_bohr, 0.0)
        .unwrap_or_else(|| vec![0.0; elements.len()]);
    let weights = model.weight_references(&charges);
    let e2 = model.get_energy(&weights);
    let e3 = if config.three_body {
        model.get_atm_energy(&pos_bohr)
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
        coordination_numbers: model.cn.clone(),
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
