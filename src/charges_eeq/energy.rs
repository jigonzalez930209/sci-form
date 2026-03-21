//! EEQ Energy and Gradients — Core
//!
//! Electrostatic energy from EEQ charges and numerical gradients.

use super::charges::{compute_eeq_charges, erf_approx, get_eeq_params, EeqConfig};

/// Result of EEQ energy calculation.
#[derive(Debug, Clone)]
pub struct EeqEnergyResult {
    /// Total electrostatic energy (kcal/mol).
    pub electrostatic_energy: f64,
    /// Per-atom charges.
    pub charges: Vec<f64>,
    /// Coordination numbers.
    pub coordination_numbers: Vec<f64>,
}

/// Compute EEQ electrostatic energy.
///
/// E_elec = Σ_i χ_i q_i + 0.5 Σ_{i,j} q_i γ(r_ij) q_j
pub fn compute_eeq_energy(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &EeqConfig,
) -> EeqEnergyResult {
    let charge_result = compute_eeq_charges(elements, positions, config);
    let q = &charge_result.charges;
    let n = elements.len();

    let mut e_elec = 0.0;

    for i in 0..n {
        let p = get_eeq_params(elements[i]);
        e_elec += p.chi * q[i];
    }

    for i in 0..n {
        let pi = get_eeq_params(elements[i]);
        for j in (i + 1)..n {
            let pj = get_eeq_params(elements[j]);
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r_ij = (dx * dx + dy * dy + dz * dz).sqrt();

            if r_ij < 1e-10 {
                continue;
            }

            let sigma_ij = (pi.r_eeq * pi.r_eeq + pj.r_eeq * pj.r_eeq).sqrt();
            let arg = std::f64::consts::SQRT_2 / sigma_ij * r_ij;
            let gij = erf_approx(arg) / r_ij;

            e_elec += q[i] * gij * q[j];
        }
    }

    let ev_to_kcal = 23.06;
    e_elec *= ev_to_kcal;

    EeqEnergyResult {
        electrostatic_energy: e_elec,
        charges: charge_result.charges,
        coordination_numbers: charge_result.coordination_numbers,
    }
}

/// Compute numerical gradient of EEQ energy.
pub fn compute_eeq_gradient(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &EeqConfig,
) -> Vec<[f64; 3]> {
    let n = elements.len();
    let h = 1e-4;
    let mut gradient = vec![[0.0; 3]; n];

    for i in 0..n {
        for d in 0..3 {
            let mut pos_p = positions.to_vec();
            let mut pos_m = positions.to_vec();
            pos_p[i][d] += h;
            pos_m[i][d] -= h;

            let ep = compute_eeq_energy(elements, &pos_p, config).electrostatic_energy;
            let em = compute_eeq_energy(elements, &pos_m, config).electrostatic_energy;

            gradient[i][d] = (ep - em) / (2.0 * h);
        }
    }

    gradient
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eeq_energy_finite() {
        let elements = [8, 1, 1];
        let pos = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let config = EeqConfig::default();
        let result = compute_eeq_energy(&elements, &pos, &config);
        assert!(result.electrostatic_energy.is_finite());
    }

    #[test]
    fn test_eeq_gradient_finite() {
        let elements = [6, 8, 1, 1];
        let pos = [
            [0.0, 0.0, 0.0],
            [1.23, 0.0, 0.0],
            [-0.6, 0.9, 0.0],
            [-0.6, -0.9, 0.0],
        ];
        let config = EeqConfig::default();
        let grad = compute_eeq_gradient(&elements, &pos, &config);
        for g in &grad {
            for &d in g {
                assert!(d.is_finite(), "Gradient contains NaN/Inf");
            }
        }
    }
}
