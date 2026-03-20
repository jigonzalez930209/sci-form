//! CPM Electrochemical Surface — E10.2
//!
//! Scans the electrochemical potential μ to compute Q(μ), Ω(μ),
//! and capacitance C(μ) = ∂Q/∂μ curves.

use super::grand_potential::{compute_cpm_charges, CpmConfig};

/// Electrochemical surface from a μ scan.
#[derive(Debug, Clone)]
pub struct CpmSurface {
    /// Applied potentials (eV).
    pub mu_values: Vec<f64>,
    /// Total charge at each μ.
    pub total_charge: Vec<f64>,
    /// Grand potential at each μ (eV).
    pub free_energy: Vec<f64>,
    /// Quantum capacitance ∂Q/∂μ at each μ (e/eV).
    pub capacitance: Vec<f64>,
    /// All converged.
    pub all_converged: bool,
}

/// Compute electrochemical surface by scanning μ.
pub fn compute_cpm_surface(
    elements: &[u8],
    positions: &[[f64; 3]],
    mu_min: f64,
    mu_max: f64,
    n_points: usize,
    dielectric: f64,
) -> CpmSurface {
    let n_points = n_points.max(2);
    let step = (mu_max - mu_min) / (n_points - 1) as f64;

    let mut mu_values = Vec::with_capacity(n_points);
    let mut total_charge = Vec::with_capacity(n_points);
    let mut free_energy = Vec::with_capacity(n_points);
    let mut all_converged = true;

    for i in 0..n_points {
        let mu = mu_min + i as f64 * step;
        let config = CpmConfig {
            mu_ev: mu,
            dielectric,
            ..Default::default()
        };

        let result = compute_cpm_charges(elements, positions, &config);
        mu_values.push(mu);
        total_charge.push(result.total_charge);
        free_energy.push(result.grand_potential);
        if !result.converged {
            all_converged = false;
        }
    }

    // Compute capacitance via finite differences: C = ΔQ/Δμ
    let mut capacitance = Vec::with_capacity(n_points);
    for i in 0..n_points {
        let c = if i == 0 {
            (total_charge[1] - total_charge[0]) / (mu_values[1] - mu_values[0])
        } else if i == n_points - 1 {
            (total_charge[i] - total_charge[i - 1]) / (mu_values[i] - mu_values[i - 1])
        } else {
            (total_charge[i + 1] - total_charge[i - 1]) / (mu_values[i + 1] - mu_values[i - 1])
        };
        capacitance.push(c);
    }

    CpmSurface {
        mu_values,
        total_charge,
        free_energy,
        capacitance,
        all_converged,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_surface_dimensions() {
        let elements = vec![6u8, 6, 1, 1, 1, 1];
        let pos = vec![
            [0.0, 0.0, 0.0], [1.5, 0.0, 0.0],
            [-0.5, 0.9, 0.0], [-0.5, -0.9, 0.0],
            [2.0, 0.9, 0.0], [2.0, -0.9, 0.0],
        ];

        let surface = compute_cpm_surface(&elements, &pos, -5.5, -3.5, 10, 78.5);
        assert_eq!(surface.mu_values.len(), 10);
        assert_eq!(surface.total_charge.len(), 10);
        assert_eq!(surface.capacitance.len(), 10);
    }

    #[test]
    fn test_charge_monotonic() {
        let elements = vec![6u8, 8, 1, 1];
        let pos = vec![
            [0.0, 0.0, 0.0], [1.2, 0.0, 0.0],
            [-0.5, 0.9, 0.0], [-0.5, -0.9, 0.0],
        ];

        let surface = compute_cpm_surface(&elements, &pos, -5.5, -3.5, 20, 78.5);

        // Q should generally increase with μ (more positive potential → more electron donation)
        let q_first = surface.total_charge.first().unwrap();
        let q_last = surface.total_charge.last().unwrap();
        assert!(q_last > q_first,
            "Q should increase with μ: first={}, last={}", q_first, q_last);
    }

    #[test]
    fn test_capacitance_positive() {
        let elements = vec![6u8, 6, 8];
        let pos = vec![
            [0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0],
        ];

        let surface = compute_cpm_surface(&elements, &pos, -5.5, -3.5, 15, 78.5);

        // Capacitance should be positive (∂Q/∂μ > 0)
        for &c in &surface.capacitance {
            assert!(c > 0.0 || c.abs() < 0.1,
                "Capacitance should be ~positive: {}", c);
        }
    }
}
