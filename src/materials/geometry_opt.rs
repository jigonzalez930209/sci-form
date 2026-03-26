//! Framework geometry optimization under periodic boundary conditions.
//!
//! Optimizes atomic positions within a fixed unit cell, useful for
//! MOF linker optimization, crystal surface relaxation, and post-assembly
//! refinement of framework structures.
//!
//! Supports:
//! - Cartesian and fractional coordinate optimization
//! - BFGS quasi-Newton with line search
//! - Steepest descent fallback
//! - Fixed-atom constraints
//! - Minimum image convention for periodic forces

use serde::{Deserialize, Serialize};

/// Optimization method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OptMethod {
    /// Steepest descent (robust, slow convergence).
    SteepestDescent,
    /// BFGS quasi-Newton (fast convergence near minimum).
    Bfgs,
}

/// Configuration for framework geometry optimization.
#[derive(Debug, Clone)]
pub struct FrameworkOptConfig {
    /// Optimization method.
    pub method: OptMethod,
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Force convergence threshold (eV/Å).
    pub force_tol: f64,
    /// Energy convergence threshold (eV).
    pub energy_tol: f64,
    /// Maximum step size (Å).
    pub max_step: f64,
    /// Indices of atoms whose positions are fixed.
    pub fixed_atoms: Vec<usize>,
    /// Unit cell lattice vectors (3×3, row-major). If None, non-periodic.
    pub lattice: Option<[[f64; 3]; 3]>,
}

impl Default for FrameworkOptConfig {
    fn default() -> Self {
        Self {
            method: OptMethod::Bfgs,
            max_iter: 200,
            force_tol: 0.05,
            energy_tol: 1e-6,
            max_step: 0.2,
            fixed_atoms: vec![],
            lattice: None,
        }
    }
}

/// Result of framework geometry optimization.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FrameworkOptResult {
    /// Optimized positions (Cartesian, Å).
    pub positions: Vec<[f64; 3]>,
    /// Final energy (eV).
    pub energy: f64,
    /// Final forces (eV/Å).
    pub forces: Vec<[f64; 3]>,
    /// Maximum force magnitude at convergence (eV/Å).
    pub max_force: f64,
    /// Number of iterations performed.
    pub n_iterations: usize,
    /// Whether optimization converged.
    pub converged: bool,
    /// Energy trajectory.
    pub energy_history: Vec<f64>,
}

/// Energy and force function type.
/// Takes (elements, positions) → (energy, forces).
pub type EnergyForceFn = dyn Fn(&[u8], &[[f64; 3]]) -> Result<(f64, Vec<[f64; 3]>), String>;

/// Optimize framework geometry using the specified method.
///
/// # Arguments
/// * `elements` - Atomic numbers
/// * `initial_positions` - Starting coordinates (Å)
/// * `energy_force_fn` - Function computing energy and forces
/// * `config` - Optimization configuration
pub fn optimize_framework(
    elements: &[u8],
    initial_positions: &[[f64; 3]],
    energy_force_fn: &EnergyForceFn,
    config: &FrameworkOptConfig,
) -> Result<FrameworkOptResult, String> {
    match config.method {
        OptMethod::SteepestDescent => {
            optimize_steepest_descent(elements, initial_positions, energy_force_fn, config)
        }
        OptMethod::Bfgs => optimize_bfgs(elements, initial_positions, energy_force_fn, config),
    }
}

/// Steepest descent optimizer.
fn optimize_steepest_descent(
    elements: &[u8],
    initial_positions: &[[f64; 3]],
    energy_force_fn: &EnergyForceFn,
    config: &FrameworkOptConfig,
) -> Result<FrameworkOptResult, String> {
    let n = elements.len();
    let mut positions: Vec<[f64; 3]> = initial_positions.to_vec();
    let mut energy_history = Vec::new();
    let mut step_size = config.max_step;

    let (mut energy, mut forces) = energy_force_fn(elements, &positions)?;
    energy_history.push(energy);

    let mut converged = false;
    let mut n_iter = 0;

    for iter in 0..config.max_iter {
        n_iter = iter + 1;

        // Zero forces on fixed atoms
        zero_fixed_forces(&mut forces, &config.fixed_atoms);

        let max_force = max_force_magnitude(&forces);
        if max_force < config.force_tol {
            converged = true;
            break;
        }

        // Take step along force direction
        let mut new_positions = positions.clone();
        for i in 0..n {
            if config.fixed_atoms.contains(&i) {
                continue;
            }
            let f_mag = (forces[i][0].powi(2) + forces[i][1].powi(2) + forces[i][2].powi(2)).sqrt();
            if f_mag < 1e-12 {
                continue;
            }
            let scale = step_size / f_mag;
            for d in 0..3 {
                new_positions[i][d] += forces[i][d] * scale;
            }
        }

        // Apply minimum image convention if periodic
        if let Some(ref lattice) = config.lattice {
            apply_pbc(&mut new_positions, lattice);
        }

        let (new_energy, new_forces) = energy_force_fn(elements, &new_positions)?;

        if new_energy < energy {
            positions = new_positions;
            energy = new_energy;
            forces = new_forces;
            step_size = (step_size * 1.2).min(config.max_step);
        } else {
            step_size *= 0.5;
            if step_size < 1e-10 {
                break;
            }
        }

        energy_history.push(energy);

        if energy_history.len() > 1 {
            let de = (energy_history[energy_history.len() - 2] - energy).abs();
            if de < config.energy_tol {
                converged = true;
                break;
            }
        }
    }

    zero_fixed_forces(&mut forces, &config.fixed_atoms);
    let max_force = max_force_magnitude(&forces);

    Ok(FrameworkOptResult {
        positions,
        energy,
        forces,
        max_force,
        n_iterations: n_iter,
        converged,
        energy_history,
    })
}

/// BFGS quasi-Newton optimizer with approximate inverse Hessian.
fn optimize_bfgs(
    elements: &[u8],
    initial_positions: &[[f64; 3]],
    energy_force_fn: &EnergyForceFn,
    config: &FrameworkOptConfig,
) -> Result<FrameworkOptResult, String> {
    let n = elements.len();
    let ndim = n * 3;
    let mut positions: Vec<[f64; 3]> = initial_positions.to_vec();
    let mut energy_history = Vec::new();

    // Initial evaluation
    let (mut energy, mut forces) = energy_force_fn(elements, &positions)?;
    zero_fixed_forces(&mut forces, &config.fixed_atoms);
    energy_history.push(energy);

    // Flatten gradient (negative force)
    let mut grad = flatten_neg_forces(&forces);

    // Initialize inverse Hessian as identity
    let mut h_inv = vec![vec![0.0f64; ndim]; ndim];
    for i in 0..ndim {
        h_inv[i][i] = 1.0;
    }

    // Zero columns/rows for fixed atoms
    for &fixed in &config.fixed_atoms {
        for d in 0..3 {
            let idx = fixed * 3 + d;
            if idx < ndim {
                for j in 0..ndim {
                    h_inv[idx][j] = 0.0;
                    h_inv[j][idx] = 0.0;
                }
            }
        }
    }

    let mut converged = false;
    let mut n_iter = 0;

    for iter in 0..config.max_iter {
        n_iter = iter + 1;

        let max_force = max_force_magnitude(&forces);
        if max_force < config.force_tol {
            converged = true;
            break;
        }

        // Search direction: p = -H_inv * grad
        let mut p = vec![0.0f64; ndim];
        for i in 0..ndim {
            for j in 0..ndim {
                p[i] -= h_inv[i][j] * grad[j];
            }
        }

        // Limit step size
        let p_norm: f64 = p.iter().map(|x| x * x).sum::<f64>().sqrt();
        if p_norm > config.max_step {
            let scale = config.max_step / p_norm;
            for x in &mut p {
                *x *= scale;
            }
        }

        // Take step with Armijo backtracking line search
        let directional_deriv: f64 = p.iter().zip(grad.iter()).map(|(a, b)| a * b).sum();
        let c_armijo = 1e-4;
        let mut alpha = 1.0;
        let mut new_positions;
        let mut new_energy;
        let mut new_forces;

        loop {
            new_positions = positions.clone();
            for i in 0..n {
                if config.fixed_atoms.contains(&i) {
                    continue;
                }
                for d in 0..3 {
                    new_positions[i][d] += alpha * p[i * 3 + d];
                }
            }

            if let Some(ref lattice) = config.lattice {
                apply_pbc(&mut new_positions, lattice);
            }

            let result = energy_force_fn(elements, &new_positions)?;
            new_energy = result.0;
            new_forces = result.1;

            // Armijo condition: f(x + α*p) <= f(x) + c * α * ∇f·p
            if new_energy <= energy + c_armijo * alpha * directional_deriv || alpha < 0.1 {
                break;
            }
            alpha *= 0.5;
        }

        zero_fixed_forces(&mut new_forces, &config.fixed_atoms);

        let new_grad = flatten_neg_forces(&new_forces);

        // BFGS update of inverse Hessian
        let s: Vec<f64> = p; // step
        let y: Vec<f64> = (0..ndim).map(|i| new_grad[i] - grad[i]).collect();

        let sy: f64 = s.iter().zip(y.iter()).map(|(a, b)| a * b).sum();

        if sy > 1e-12 {
            // H_inv update: Sherman-Morrison-Woodbury
            let mut hy = vec![0.0f64; ndim];
            for i in 0..ndim {
                for j in 0..ndim {
                    hy[i] += h_inv[i][j] * y[j];
                }
            }

            let yhy: f64 = y.iter().zip(hy.iter()).map(|(a, b)| a * b).sum();
            let rho = 1.0 / sy;

            for i in 0..ndim {
                for j in 0..ndim {
                    h_inv[i][j] +=
                        rho * ((1.0 + yhy * rho) * s[i] * s[j] - hy[i] * s[j] - s[i] * hy[j]);
                }
            }

            // Positive-definite check: if any diagonal becomes negative, reset to identity
            let has_negative_diag = (0..ndim).any(|i| h_inv[i][i] <= 0.0);
            if has_negative_diag {
                for i in 0..ndim {
                    for j in 0..ndim {
                        h_inv[i][j] = if i == j { 1.0 } else { 0.0 };
                    }
                }
            }
        }

        positions = new_positions;
        energy = new_energy;
        forces = new_forces;
        grad = new_grad;
        energy_history.push(energy);

        if energy_history.len() > 1 {
            let de = (energy_history[energy_history.len() - 2] - energy).abs();
            if de < config.energy_tol {
                converged = true;
                break;
            }
        }
    }

    let max_force = max_force_magnitude(&forces);

    Ok(FrameworkOptResult {
        positions,
        energy,
        forces,
        max_force,
        n_iterations: n_iter,
        converged,
        energy_history,
    })
}

fn zero_fixed_forces(forces: &mut [[f64; 3]], fixed: &[usize]) {
    for &idx in fixed {
        if idx < forces.len() {
            forces[idx] = [0.0, 0.0, 0.0];
        }
    }
}

fn max_force_magnitude(forces: &[[f64; 3]]) -> f64 {
    forces
        .iter()
        .map(|f| (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt())
        .fold(0.0f64, f64::max)
}

fn flatten_neg_forces(forces: &[[f64; 3]]) -> Vec<f64> {
    let mut g = Vec::with_capacity(forces.len() * 3);
    for f in forces {
        g.push(-f[0]);
        g.push(-f[1]);
        g.push(-f[2]);
    }
    g
}

/// Apply periodic boundary conditions: wrap Cartesian coordinates back into the unit cell.
fn apply_pbc(positions: &mut [[f64; 3]], lattice: &[[f64; 3]; 3]) {
    // Compute inverse lattice matrix
    let inv = invert_3x3_lattice(lattice);

    for pos in positions.iter_mut() {
        // Convert to fractional
        let frac = [
            inv[0][0] * pos[0] + inv[0][1] * pos[1] + inv[0][2] * pos[2],
            inv[1][0] * pos[0] + inv[1][1] * pos[1] + inv[1][2] * pos[2],
            inv[2][0] * pos[0] + inv[2][1] * pos[1] + inv[2][2] * pos[2],
        ];

        // Wrap to [0, 1)
        let wrapped = [
            frac[0] - frac[0].floor(),
            frac[1] - frac[1].floor(),
            frac[2] - frac[2].floor(),
        ];

        // Convert back to Cartesian
        pos[0] =
            lattice[0][0] * wrapped[0] + lattice[1][0] * wrapped[1] + lattice[2][0] * wrapped[2];
        pos[1] =
            lattice[0][1] * wrapped[0] + lattice[1][1] * wrapped[1] + lattice[2][1] * wrapped[2];
        pos[2] =
            lattice[0][2] * wrapped[0] + lattice[1][2] * wrapped[1] + lattice[2][2] * wrapped[2];
    }
}

fn invert_3x3_lattice(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

    if det.abs() < 1e-30 {
        return [[0.0; 3]; 3];
    }

    let inv_det = 1.0 / det;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ]
}

/// Convert fractional coordinates to Cartesian given lattice vectors.
pub fn frac_to_cart(frac: &[f64; 3], lattice: &[[f64; 3]; 3]) -> [f64; 3] {
    [
        lattice[0][0] * frac[0] + lattice[1][0] * frac[1] + lattice[2][0] * frac[2],
        lattice[0][1] * frac[0] + lattice[1][1] * frac[1] + lattice[2][1] * frac[2],
        lattice[0][2] * frac[0] + lattice[1][2] * frac[1] + lattice[2][2] * frac[2],
    ]
}

/// Convert Cartesian coordinates to fractional given lattice vectors.
pub fn cart_to_frac(cart: &[f64; 3], lattice: &[[f64; 3]; 3]) -> [f64; 3] {
    let inv = invert_3x3_lattice(lattice);
    [
        inv[0][0] * cart[0] + inv[0][1] * cart[1] + inv[0][2] * cart[2],
        inv[1][0] * cart[0] + inv[1][1] * cart[1] + inv[1][2] * cart[2],
        inv[2][0] * cart[0] + inv[2][1] * cart[1] + inv[2][2] * cart[2],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_harmonic_energy(
        _elements: &[u8],
        positions: &[[f64; 3]],
    ) -> Result<(f64, Vec<[f64; 3]>), String> {
        // Simple harmonic well centered at origin for each atom
        let mut energy = 0.0;
        let mut forces = Vec::with_capacity(positions.len());
        for pos in positions {
            let r2 = pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];
            energy += 0.5 * r2;
            forces.push([-pos[0], -pos[1], -pos[2]]); // F = -grad(E)
        }
        Ok((energy, forces))
    }

    #[test]
    fn test_steepest_descent() {
        let elements = vec![6u8];
        let initial = vec![[1.0, 0.5, 0.2]];
        let config = FrameworkOptConfig {
            method: OptMethod::SteepestDescent,
            max_iter: 100,
            force_tol: 0.01,
            ..Default::default()
        };

        let result =
            optimize_framework(&elements, &initial, &simple_harmonic_energy, &config).unwrap();
        assert!(result.converged);
        assert!(result.positions[0][0].abs() < 0.1);
        assert!(result.positions[0][1].abs() < 0.1);
    }

    #[test]
    fn test_bfgs() {
        let elements = vec![6u8];
        let initial = vec![[1.0, 0.5, 0.2]];
        let config = FrameworkOptConfig {
            method: OptMethod::Bfgs,
            max_iter: 50,
            force_tol: 0.01,
            ..Default::default()
        };

        let result =
            optimize_framework(&elements, &initial, &simple_harmonic_energy, &config).unwrap();
        assert!(result.converged);
        // BFGS should converge faster
        assert!(result.n_iterations < 20);
    }

    #[test]
    fn test_fixed_atoms() {
        let elements = vec![6u8, 8u8];
        let initial = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let config = FrameworkOptConfig {
            method: OptMethod::Bfgs,
            max_iter: 50,
            force_tol: 0.01,
            fixed_atoms: vec![0], // Fix first atom
            ..Default::default()
        };

        let result =
            optimize_framework(&elements, &initial, &simple_harmonic_energy, &config).unwrap();
        // First atom should remain at (1,0,0)
        assert!((result.positions[0][0] - 1.0).abs() < 1e-10);
        // Second atom should move toward origin
        assert!(result.positions[1][1].abs() < 0.2);
    }

    #[test]
    fn test_frac_cart_conversion() {
        let lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        let frac = [0.5, 0.25, 0.1];
        let cart = frac_to_cart(&frac, &lattice);
        assert!((cart[0] - 5.0).abs() < 1e-10);
        assert!((cart[1] - 2.5).abs() < 1e-10);
        assert!((cart[2] - 1.0).abs() < 1e-10);

        let back = cart_to_frac(&cart, &lattice);
        assert!((back[0] - 0.5).abs() < 1e-10);
    }
}
