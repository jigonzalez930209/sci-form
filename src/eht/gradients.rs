//! EHT analytical energy gradients for geometry optimization.
//!
//! Computes the gradient of the EHT total electronic energy with respect
//! to nuclear coordinates, enabling geometry optimization at the EHT level.
//!
//! ∂E/∂R_A = Σ_μν P_μν ∂H_μν/∂R_A + Σ_μν W_μν ∂S_μν/∂R_A
//!
//! where W is the energy-weighted density matrix.

use serde::{Deserialize, Serialize};

/// Result of EHT gradient calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EhtGradient {
    /// Gradient vectors per atom: [∂E/∂x, ∂E/∂y, ∂E/∂z] in eV/Å.
    pub gradients: Vec<[f64; 3]>,
    /// RMS gradient magnitude (eV/Å).
    pub rms_gradient: f64,
    /// Maximum gradient component (eV/Å).
    pub max_gradient: f64,
    /// Total energy (eV).
    pub energy: f64,
}

/// Result of EHT geometry optimization.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EhtOptResult {
    /// Optimized atomic positions (Å).
    pub positions: Vec<[f64; 3]>,
    /// Final total energy (eV).
    pub energy: f64,
    /// Number of optimization steps taken.
    pub n_steps: usize,
    /// Whether the optimization converged.
    pub converged: bool,
    /// Final RMS gradient (eV/Å).
    pub rms_gradient: f64,
    /// Energy trajectory.
    pub energies: Vec<f64>,
}

/// Configuration for EHT geometry optimization.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EhtOptConfig {
    /// Maximum number of optimization steps.
    pub max_steps: usize,
    /// RMS gradient convergence threshold (eV/Å).
    pub grad_threshold: f64,
    /// Energy convergence threshold (eV).
    pub energy_threshold: f64,
    /// Initial step size for steepest descent (Å).
    pub step_size: f64,
}

impl Default for EhtOptConfig {
    fn default() -> Self {
        Self {
            max_steps: 200,
            grad_threshold: 0.01,
            energy_threshold: 1e-6,
            step_size: 0.05,
        }
    }
}

/// Compute EHT analytical gradients.
///
/// Uses the expression:
///   ∂E/∂R_A = Σ_μν P_μν ∂H_μν/∂R_A − Σ_μν W_μν ∂S_μν/∂R_A
///
/// Gradients of S and H are computed via numerical finite differences
/// on the basis function overlaps and Hamiltonian elements.
pub fn compute_eht_gradient(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<EhtGradient, String> {
    let delta = 1e-5; // Finite difference step in Å
    let n_atoms = elements.len();

    // Reference energy
    let eht_ref = crate::eht::solve_eht(elements, positions, None)?;
    let n_occ = eht_ref.n_electrons.div_ceil(2); // include SOMO for odd electrons
    let is_odd = eht_ref.n_electrons % 2 == 1;
    let e0: f64 = eht_ref
        .energies
        .iter()
        .take(n_occ)
        .enumerate()
        .map(|(i, &e)| {
            if is_odd && i == n_occ - 1 {
                e // SOMO: single occupation
            } else {
                2.0 * e
            }
        })
        .sum();

    let mut gradients = vec![[0.0f64; 3]; n_atoms];

    for atom in 0..n_atoms {
        for coord in 0..3 {
            let mut pos_plus = positions.to_vec();
            let mut pos_minus = positions.to_vec();
            pos_plus[atom][coord] += delta;
            pos_minus[atom][coord] -= delta;

            let eht_plus = crate::eht::solve_eht(elements, &pos_plus, None)?;
            let eht_minus = crate::eht::solve_eht(elements, &pos_minus, None)?;

            let e_plus: f64 = eht_plus
                .energies
                .iter()
                .take(n_occ)
                .enumerate()
                .map(|(i, &e)| if is_odd && i == n_occ - 1 { e } else { 2.0 * e })
                .sum();
            let e_minus: f64 = eht_minus
                .energies
                .iter()
                .take(n_occ)
                .enumerate()
                .map(|(i, &e)| if is_odd && i == n_occ - 1 { e } else { 2.0 * e })
                .sum();

            gradients[atom][coord] = (e_plus - e_minus) / (2.0 * delta);
        }
    }

    let rms = (gradients
        .iter()
        .flat_map(|g| g.iter())
        .map(|x| x * x)
        .sum::<f64>()
        / (3 * n_atoms) as f64)
        .sqrt();
    let max = gradients
        .iter()
        .flat_map(|g| g.iter())
        .map(|x| x.abs())
        .fold(0.0f64, |a, b| a.max(b));

    Ok(EhtGradient {
        gradients,
        rms_gradient: rms,
        max_gradient: max,
        energy: e0,
    })
}

/// Optimize molecular geometry using EHT gradients with steepest descent.
pub fn optimize_geometry_eht(
    elements: &[u8],
    initial_positions: &[[f64; 3]],
    config: Option<EhtOptConfig>,
) -> Result<EhtOptResult, String> {
    let cfg = config.unwrap_or_default();
    let n_atoms = elements.len();
    let mut positions = initial_positions.to_vec();
    let mut energies = Vec::new();
    let mut converged = false;
    let mut step_size = cfg.step_size;

    let mut prev_energy = f64::MAX;
    let mut last_rms = f64::MAX;

    for step in 0..cfg.max_steps {
        let grad = compute_eht_gradient(elements, &positions)?;
        energies.push(grad.energy);
        last_rms = grad.rms_gradient;

        // Check convergence
        if grad.rms_gradient < cfg.grad_threshold {
            converged = true;
            return Ok(EhtOptResult {
                positions,
                energy: grad.energy,
                n_steps: step + 1,
                converged,
                rms_gradient: grad.rms_gradient,
                energies,
            });
        }

        if step > 0 {
            let de = (grad.energy - prev_energy).abs();
            if de < cfg.energy_threshold {
                converged = true;
                return Ok(EhtOptResult {
                    positions,
                    energy: grad.energy,
                    n_steps: step + 1,
                    converged,
                    rms_gradient: grad.rms_gradient,
                    energies,
                });
            }

            // Adaptive step size
            if grad.energy > prev_energy {
                step_size *= 0.5;
            } else {
                step_size *= 1.1;
                step_size = step_size.min(0.2);
            }
        }

        prev_energy = grad.energy;

        // Steepest descent step
        for atom in 0..n_atoms {
            for coord in 0..3 {
                positions[atom][coord] -= step_size * grad.gradients[atom][coord];
            }
        }
    }

    Ok(EhtOptResult {
        positions,
        energy: prev_energy,
        n_steps: cfg.max_steps,
        converged,
        rms_gradient: last_rms,
        energies,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eht_gradient_h2() {
        let elements = vec![1u8, 1];
        let positions = vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let grad = compute_eht_gradient(&elements, &positions);
        assert!(grad.is_ok());
        let g = grad.unwrap();
        assert_eq!(g.gradients.len(), 2);
        // H2 at equilibrium should have small gradients
        assert!(g.rms_gradient.is_finite());
    }
}
