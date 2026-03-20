//! Geometry optimizer using L-BFGS and steepest descent.
//!
//! Minimizes the total energy E(R) with respect to nuclear coordinates R
//! using gradient information. Supports two algorithms:
//!
//! - **Steepest descent**: Robust but slow. Good for far-from-minimum geometries.
//! - **L-BFGS**: Quasi-Newton with limited-memory BFGS Hessian approximation.
//!   Fast convergence near the minimum.
//!
//! Convergence criteria (following Gaussian defaults):
//! - Max gradient component < 4.5e-4 Hartree/Bohr
//! - RMS gradient < 3.0e-4 Hartree/Bohr
//! - Max displacement < 1.8e-3 Bohr
//! - RMS displacement < 1.2e-3 Bohr

use crate::experimental_2::types::{MolecularSystem, OptimizationResult};

use super::gradients::numerical_gradient;
use super::scf_loop::ScfConfig;

/// Optimization algorithm.
#[derive(Debug, Clone, Copy)]
pub enum OptAlgorithm {
    SteepestDescent,
    LBFGS { memory: usize },
}

/// Configuration for geometry optimization.
#[derive(Debug, Clone)]
pub struct OptConfig {
    /// Maximum optimization steps.
    pub max_steps: usize,
    /// Max gradient convergence threshold (Hartree/Bohr).
    pub grad_max_threshold: f64,
    /// RMS gradient convergence threshold.
    pub grad_rms_threshold: f64,
    /// Max displacement convergence threshold (Bohr).
    pub disp_max_threshold: f64,
    /// RMS displacement convergence threshold (Bohr).
    pub disp_rms_threshold: f64,
    /// Numerical gradient step size (Bohr).
    pub gradient_step: f64,
    /// Step size for steepest descent.
    pub step_size: f64,
    /// Optimization algorithm.
    pub algorithm: OptAlgorithm,
    /// SCF configuration for energy/gradient evaluation.
    pub scf_config: ScfConfig,
}

impl Default for OptConfig {
    fn default() -> Self {
        Self {
            max_steps: 100,
            grad_max_threshold: 4.5e-4,
            grad_rms_threshold: 3.0e-4,
            disp_max_threshold: 1.8e-3,
            disp_rms_threshold: 1.2e-3,
            gradient_step: 1e-3,
            step_size: 0.1,
            algorithm: OptAlgorithm::LBFGS { memory: 10 },
            scf_config: ScfConfig::default(),
        }
    }
}

/// L-BFGS storage for history vectors.
struct LbfgsHistory {
    s_history: Vec<Vec<f64>>,
    y_history: Vec<Vec<f64>>,
    rho: Vec<f64>,
    max_memory: usize,
}

impl LbfgsHistory {
    fn new(memory: usize) -> Self {
        Self {
            s_history: Vec::new(),
            y_history: Vec::new(),
            rho: Vec::new(),
            max_memory: memory,
        }
    }

    /// Add a new (s, y) pair where s = x_new - x_old, y = g_new - g_old.
    fn push(&mut self, s: Vec<f64>, y: Vec<f64>) {
        let sy: f64 = s.iter().zip(y.iter()).map(|(a, b)| a * b).sum();
        if sy.abs() < 1e-16 {
            return; // Skip degenerate updates
        }
        let rho = 1.0 / sy;

        if self.s_history.len() >= self.max_memory {
            self.s_history.remove(0);
            self.y_history.remove(0);
            self.rho.remove(0);
        }

        self.s_history.push(s);
        self.y_history.push(y);
        self.rho.push(rho);
    }

    /// L-BFGS two-loop recursion to compute search direction.
    fn compute_direction(&self, gradient: &[f64]) -> Vec<f64> {
        let n = gradient.len();
        let k = self.s_history.len();

        if k == 0 {
            // No history: return steepest descent direction
            return gradient.iter().map(|&g| -g).collect();
        }

        let mut q = gradient.to_vec();
        let mut alpha = vec![0.0; k];

        // First loop (backward)
        for i in (0..k).rev() {
            let dot: f64 = self.s_history[i].iter().zip(q.iter()).map(|(a, b)| a * b).sum();
            alpha[i] = self.rho[i] * dot;
            for j in 0..n {
                q[j] -= alpha[i] * self.y_history[i][j];
            }
        }

        // Initial Hessian approximation: H⁰ = γI where γ = s^T y / y^T y
        let last = k - 1;
        let sy: f64 = self.s_history[last]
            .iter()
            .zip(self.y_history[last].iter())
            .map(|(a, b)| a * b)
            .sum();
        let yy: f64 = self.y_history[last].iter().map(|y| y * y).sum();
        let gamma = if yy.abs() > 1e-16 { sy / yy } else { 1.0 };

        let mut r: Vec<f64> = q.iter().map(|&qi| gamma * qi).collect();

        // Second loop (forward)
        for i in 0..k {
            let dot: f64 = self.y_history[i].iter().zip(r.iter()).map(|(a, b)| a * b).sum();
            let beta = self.rho[i] * dot;
            for j in 0..n {
                r[j] += self.s_history[i][j] * (alpha[i] - beta);
            }
        }

        // Negate for descent direction
        r.iter().map(|&x| -x).collect()
    }
}

/// Flatten gradients from [[f64; 3]; N] to [f64; 3N].
fn flatten_gradients(grads: &[[f64; 3]]) -> Vec<f64> {
    grads.iter().flat_map(|g| g.iter().copied()).collect()
}

/// Apply displacement to molecular coordinates.
fn apply_displacement(system: &MolecularSystem, direction: &[f64], step: f64) -> MolecularSystem {
    let mut new_system = system.clone();
    for (atom, pos) in new_system.positions_bohr.iter_mut().enumerate() {
        for coord in 0..3 {
            pos[coord] += step * direction[atom * 3 + coord];
        }
    }
    new_system
}

/// Compute displacement statistics.
fn displacement_stats(direction: &[f64], step: f64) -> (f64, f64) {
    let displacements: Vec<f64> = direction.iter().map(|&d| (d * step).abs()).collect();
    let max_disp = displacements
        .iter()
        .copied()
        .fold(0.0f64, f64::max);
    let rms_disp = (displacements.iter().map(|d| d * d).sum::<f64>() / displacements.len() as f64).sqrt();
    (max_disp, rms_disp)
}

/// Run geometry optimization.
pub fn optimize_geometry(
    system: &MolecularSystem,
    config: &OptConfig,
) -> OptimizationResult {
    let n_atoms = system.n_atoms();
    let _n_coords = n_atoms * 3;

    let mut current = system.clone();
    let mut energies = Vec::new();
    let mut gradient_norms = Vec::new();

    // Initial gradient
    let mut grad_result = numerical_gradient(&current, &config.scf_config, config.gradient_step);
    energies.push(grad_result.energy);
    gradient_norms.push(grad_result.rms_gradient);

    let mut converged = false;
    let mut n_steps = 0;

    match config.algorithm {
        OptAlgorithm::SteepestDescent => {
            for step in 0..config.max_steps {
                n_steps = step + 1;

                let flat_grad = flatten_gradients(&grad_result.gradients);

                // Search direction = -gradient
                let direction: Vec<f64> = flat_grad.iter().map(|&g| -g).collect();

                // Apply step
                current = apply_displacement(&current, &direction, config.step_size);

                // Recompute gradient
                grad_result = numerical_gradient(&current, &config.scf_config, config.gradient_step);
                energies.push(grad_result.energy);
                gradient_norms.push(grad_result.rms_gradient);

                let (max_disp, rms_disp) = displacement_stats(&direction, config.step_size);

                // Check convergence
                if grad_result.max_gradient < config.grad_max_threshold
                    && grad_result.rms_gradient < config.grad_rms_threshold
                    && max_disp < config.disp_max_threshold
                    && rms_disp < config.disp_rms_threshold
                {
                    converged = true;
                    break;
                }
            }
        }
        OptAlgorithm::LBFGS { memory } => {
            let mut lbfgs = LbfgsHistory::new(memory);
            let mut prev_grad = flatten_gradients(&grad_result.gradients);
            let mut prev_coords: Vec<f64> = current
                .positions_bohr
                .iter()
                .flat_map(|p| p.iter().copied())
                .collect();

            for step in 0..config.max_steps {
                n_steps = step + 1;

                // Compute search direction
                let direction = lbfgs.compute_direction(&prev_grad);

                // Line search: simple backtracking with Armijo condition
                let step_size = backtracking_line_search(
                    &current,
                    &direction,
                    &prev_grad,
                    grad_result.energy,
                    &config.scf_config,
                    config.gradient_step,
                    config.step_size,
                );

                // Apply step
                current = apply_displacement(&current, &direction, step_size);

                // Recompute gradient
                grad_result = numerical_gradient(&current, &config.scf_config, config.gradient_step);
                energies.push(grad_result.energy);
                gradient_norms.push(grad_result.rms_gradient);

                let new_grad = flatten_gradients(&grad_result.gradients);
                let new_coords: Vec<f64> = current
                    .positions_bohr
                    .iter()
                    .flat_map(|p| p.iter().copied())
                    .collect();

                // Update L-BFGS history
                let s: Vec<f64> = new_coords
                    .iter()
                    .zip(prev_coords.iter())
                    .map(|(a, b)| a - b)
                    .collect();
                let y: Vec<f64> = new_grad
                    .iter()
                    .zip(prev_grad.iter())
                    .map(|(a, b)| a - b)
                    .collect();
                lbfgs.push(s, y);

                let (max_disp, rms_disp) = displacement_stats(&direction, step_size);

                // Check convergence
                if grad_result.max_gradient < config.grad_max_threshold
                    && grad_result.rms_gradient < config.grad_rms_threshold
                    && max_disp < config.disp_max_threshold
                    && rms_disp < config.disp_rms_threshold
                {
                    converged = true;
                    break;
                }

                prev_grad = new_grad;
                prev_coords = new_coords;
            }
        }
    }

    OptimizationResult {
        optimized_positions: current.positions_bohr,
        final_energy: *energies.last().unwrap_or(&0.0),
        converged,
        n_steps,
        gradient_norm: grad_result.rms_gradient,
        energy_trajectory: energies,
    }
}

/// Simple backtracking line search with Armijo condition.
///
/// Finds α such that: f(x + αd) ≤ f(x) + c₁ α ∇f·d
fn backtracking_line_search(
    system: &MolecularSystem,
    direction: &[f64],
    gradient: &[f64],
    current_energy: f64,
    scf_config: &ScfConfig,
    grad_step: f64,
    initial_step: f64,
) -> f64 {
    let c1 = 1e-4; // Armijo parameter
    let shrink = 0.5;
    let min_step = 1e-8;

    let directional_derivative: f64 = gradient
        .iter()
        .zip(direction.iter())
        .map(|(g, d)| g * d)
        .sum();

    // If direction is not a descent direction, use steepest descent
    if directional_derivative >= 0.0 {
        return initial_step * 0.01;
    }

    let mut alpha = initial_step;

    for _ in 0..10 {
        let trial = apply_displacement(system, direction, alpha);
        let trial_result = numerical_gradient(&trial, scf_config, grad_step);
        let trial_energy = trial_result.energy;

        // Armijo condition
        if trial_energy <= current_energy + c1 * alpha * directional_derivative {
            return alpha;
        }

        alpha *= shrink;
        if alpha < min_step {
            break;
        }
    }

    alpha
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lbfgs_history() {
        let mut lbfgs = LbfgsHistory::new(3);

        // Without history, direction should be steepest descent
        let grad = vec![1.0, 2.0, 3.0];
        let dir = lbfgs.compute_direction(&grad);
        assert_eq!(dir, vec![-1.0, -2.0, -3.0]);

        // Add one pair
        lbfgs.push(vec![0.1, 0.2, 0.3], vec![0.05, 0.1, 0.15]);

        // Direction should still be a descent direction
        let dir = lbfgs.compute_direction(&grad);
        let dot: f64 = dir.iter().zip(grad.iter()).map(|(d, g)| d * g).sum();
        assert!(dot < 0.0, "L-BFGS direction should be a descent direction");
    }

    #[test]
    fn test_flatten_gradients() {
        let grads = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let flat = flatten_gradients(&grads);
        assert_eq!(flat, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    }

    #[test]
    fn test_displacement_stats() {
        let dir = vec![1.0, 0.0, 0.0, 0.0, 2.0, 0.0];
        let (max_d, rms_d) = displacement_stats(&dir, 0.1);
        assert!((max_d - 0.2).abs() < 1e-10);
        assert!(rms_d > 0.0);
    }
}
