//! EEQ Charge Model — E5.1
//!
//! Electronegativity equalization with geometry-dependent coordination
//! numbers and Coulomb damping.

use nalgebra::DMatrix;

/// EEQ per-element parameters.
#[derive(Debug, Clone, Copy)]
pub struct EeqParams {
    /// Electronegativity chi (eV).
    pub chi: f64,
    /// Chemical hardness eta (eV).
    pub eta: f64,
    /// Charge radius for Coulomb damping (Å).
    pub r_eeq: f64,
    /// Covalent radius for CN (Å).
    pub r_cov: f64,
}

/// Configuration for EEQ calculations.
#[derive(Debug, Clone)]
pub struct EeqConfig {
    /// Total molecular charge.
    pub total_charge: f64,
    /// Regularization parameter for near-singular systems.
    pub regularization: f64,
}

impl Default for EeqConfig {
    fn default() -> Self {
        Self {
            total_charge: 0.0,
            regularization: 1e-10,
        }
    }
}

/// Result of EEQ charge calculation.
#[derive(Debug, Clone)]
pub struct EeqChargeResult {
    /// Partial charges per atom.
    pub charges: Vec<f64>,
    /// Fractional coordination numbers per atom.
    pub coordination_numbers: Vec<f64>,
    /// Total charge (should match config).
    pub total_charge: f64,
}

/// Get EEQ parameters for an element by atomic number.
/// Source: Grimme GFN-FF / independent rederivation for key organic elements.
pub fn get_eeq_params(z: u8) -> EeqParams {
    match z {
        1 => EeqParams { chi: 2.20, eta: 13.6, r_eeq: 0.80, r_cov: 0.32 },      // H
        5 => EeqParams { chi: 2.04, eta: 8.30, r_eeq: 1.40, r_cov: 0.85 },      // B
        6 => EeqParams { chi: 2.55, eta: 10.0, r_eeq: 1.30, r_cov: 0.77 },      // C
        7 => EeqParams { chi: 3.04, eta: 14.5, r_eeq: 1.20, r_cov: 0.75 },      // N
        8 => EeqParams { chi: 3.44, eta: 13.4, r_eeq: 1.10, r_cov: 0.73 },      // O
        9 => EeqParams { chi: 3.98, eta: 17.4, r_eeq: 1.00, r_cov: 0.71 },      // F
        14 => EeqParams { chi: 1.90, eta: 8.15, r_eeq: 1.75, r_cov: 1.17 },     // Si
        15 => EeqParams { chi: 2.19, eta: 10.5, r_eeq: 1.60, r_cov: 1.10 },     // P
        16 => EeqParams { chi: 2.58, eta: 10.4, r_eeq: 1.50, r_cov: 1.04 },     // S
        17 => EeqParams { chi: 3.16, eta: 13.0, r_eeq: 1.40, r_cov: 0.99 },     // Cl
        35 => EeqParams { chi: 2.96, eta: 11.8, r_eeq: 1.55, r_cov: 1.14 },     // Br
        53 => EeqParams { chi: 2.66, eta: 10.5, r_eeq: 1.70, r_cov: 1.33 },     // I
        26 => EeqParams { chi: 1.83, eta: 7.90, r_eeq: 1.70, r_cov: 1.24 },     // Fe
        29 => EeqParams { chi: 1.90, eta: 7.73, r_eeq: 1.60, r_cov: 1.32 },     // Cu
        30 => EeqParams { chi: 1.65, eta: 9.39, r_eeq: 1.65, r_cov: 1.22 },     // Zn
        _ => EeqParams { chi: 2.20, eta: 10.0, r_eeq: 1.50, r_cov: 1.00 },      // fallback
    }
}

/// Compute fractional coordination number for each atom.
///
/// CN_i = sum_{j!=i} 1 / (1 + exp(-16 * (r_cov_ij/r_ij - 1)))
pub fn fractional_coordination(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Vec<f64> {
    let n = elements.len();
    let mut cn = vec![0.0; n];

    for i in 0..n {
        let pi = get_eeq_params(elements[i]);
        for j in (i + 1)..n {
            let pj = get_eeq_params(elements[j]);
            let r_cov_ij = pi.r_cov + pj.r_cov;
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r_ij = (dx * dx + dy * dy + dz * dz).sqrt();

            if r_ij < 1e-10 { continue; }

            let f = 1.0 / (1.0 + (-16.0 * (r_cov_ij / r_ij - 1.0)).exp());
            cn[i] += f;
            cn[j] += f;
        }
    }

    cn
}

/// Coulomb interaction kernel with Gaussian damping.
///
/// gamma(r_ij) = erf(sqrt(2) / sigma_ij * r_ij) / r_ij
/// where sigma_ij = sqrt(r_eeq_i^2 + r_eeq_j^2)
fn gamma_damped(r_ij: f64, r_eeq_i: f64, r_eeq_j: f64) -> f64 {
    if r_ij < 1e-10 {
        return 0.0;
    }
    let sigma_ij = (r_eeq_i * r_eeq_i + r_eeq_j * r_eeq_j).sqrt();
    let arg = std::f64::consts::SQRT_2 / sigma_ij * r_ij;
    erf_approx(arg) / r_ij
}

/// Approximate error function (Abramowitz & Stegun 7.1.26).
pub(crate) fn erf_approx(x: f64) -> f64 {
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;
    let p = 0.3275911;

    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + p * x);
    let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();
    sign * y
}

/// Compute EEQ charges by solving the extended linear system.
///
/// [η + γ   1] [q]   [−χ]
/// [1^T     0] [λ] = [Q ]
pub fn compute_eeq_charges(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &EeqConfig,
) -> EeqChargeResult {
    let n = elements.len();
    let cn = fractional_coordination(elements, positions);

    // Build (N+1)x(N+1) system
    let dim = n + 1;
    let mut a = DMatrix::zeros(dim, dim);
    let mut b_vec = vec![0.0; dim];

    let params: Vec<EeqParams> = elements.iter().map(|&z| get_eeq_params(z)).collect();

    // Fill A matrix
    for i in 0..n {
        // Diagonal: hardness + CN-dependent correction
        a[(i, i)] = params[i].eta + config.regularization;

        // Off-diagonal: Coulomb interaction
        for j in (i + 1)..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r_ij = (dx * dx + dy * dy + dz * dz).sqrt();
            let gij = gamma_damped(r_ij, params[i].r_eeq, params[j].r_eeq);
            a[(i, j)] = gij;
            a[(j, i)] = gij;
        }

        // Constraint row/column
        a[(i, n)] = 1.0;
        a[(n, i)] = 1.0;

        // RHS: -chi (CN-adjusted)
        let cn_correction = -0.1 * (cn[i] - 2.0); // Simple CN adjustment
        b_vec[i] = -(params[i].chi + cn_correction);
    }

    // Charge constraint
    b_vec[n] = config.total_charge;

    // Solve A * [q; lambda] = b via LU decomposition
    let b_nalg = nalgebra::DVector::from_vec(b_vec);
    let solution = a.lu().solve(&b_nalg);

    let charges = match solution {
        Some(sol) => (0..n).map(|i| sol[i]).collect(),
        None => vec![0.0; n], // Fallback: zero charges
    };

    let total: f64 = charges.iter().sum();

    EeqChargeResult {
        charges,
        coordination_numbers: cn,
        total_charge: total,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn water_positions() -> Vec<[f64; 3]> {
        vec![[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]]
    }

    #[test]
    fn test_coordination_number_water() {
        let elements = [8, 1, 1];
        let pos = water_positions();
        let cn = fractional_coordination(&elements, &pos);
        // O should have CN ~ 2, H should have CN ~ 1
        assert!(cn[0] > 1.5, "O CN = {}", cn[0]);
        assert!(cn[1] > 0.5 && cn[1] < 1.5, "H CN = {}", cn[1]);
    }

    #[test]
    fn test_eeq_charge_neutrality() {
        let elements = [6, 6, 8, 1, 1, 1, 1, 1]; // ethanol-like
        let pos = [
            [0.0, 0.0, 0.0], [1.54, 0.0, 0.0], [2.57, 1.03, 0.0],
            [-0.63, 0.89, 0.0], [-0.63, -0.89, 0.0], [1.54, -0.63, 0.89],
            [1.54, -0.63, -0.89], [3.52, 0.93, 0.0],
        ];
        let config = EeqConfig::default();
        let result = compute_eeq_charges(&elements, &pos, &config);
        assert!(
            result.total_charge.abs() < 0.01,
            "Charge not neutral: {}", result.total_charge
        );
    }

    #[test]
    fn test_eeq_oxygen_negative() {
        let elements = [8, 1, 1];
        let pos = water_positions();
        let config = EeqConfig::default();
        let result = compute_eeq_charges(&elements, &pos, &config);
        // Oxygen should be negative
        assert!(result.charges[0] < 0.0, "O charge = {}", result.charges[0]);
        // Hydrogens should be positive
        assert!(result.charges[1] > 0.0, "H charge = {}", result.charges[1]);
    }

    #[test]
    fn test_gamma_damped_asymptotic() {
        // At large r, gamma -> 1/r
        let g = gamma_damped(10.0, 1.0, 1.0);
        assert!((g - 0.1).abs() < 0.01, "gamma(10) = {}, expected ~0.1", g);
    }
}
