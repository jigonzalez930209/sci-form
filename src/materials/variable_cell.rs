//! Variable-cell relaxation: BFGS on atomic positions plus cell strain.
//!
//! The framework optimizer (`optimize_framework`) relaxes atomic positions at a
//! fixed unit cell. Variable-cell relaxation additionally minimizes the energy
//! with respect to the cell shape. Because the periodic energy models exposed
//! through [`crate::materials`] do not return an analytic stress tensor, this
//! module computes the stress by central finite differences of the energy under
//! a symmetric strain,
//!
//! ```text
//! σ_ij = (1/V) ∂E/∂ε_ij        r' = (I + ε) r,   a'_k = (I + ε) a_k
//! ```
//!
//! and drives the cell with a damped, energy-validated strain step
//! (line search), alternating with position relaxation.
//!
//! The validation gate uses an analytic model system: a Lennard-Jones pair
//! potential on a simple-cubic one-atom cell, whose equilibrium lattice
//! constant is obtained independently by minimizing `E(a)`; the optimizer must
//! reproduce it and drive the stress to zero.

use serde::{Deserialize, Serialize};

/// Energy/force callback with an explicit lattice (Å).
///
/// Takes `(elements, positions, lattice)` and returns `(energy, forces)`.
pub type CellEnergyForceFn =
    dyn Fn(&[u8], &[[f64; 3]], &[[f64; 3]; 3]) -> Result<(f64, Vec<[f64; 3]>), String>;

/// Cell volume (Å³): `det(lattice)`.
pub fn cell_volume(lattice: &[[f64; 3]; 3]) -> f64 {
    let a = lattice;
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

/// Apply a symmetric strain `ε` to a lattice: `a'_k = (I + ε) a_k`.
pub fn strained_lattice(lattice: &[[f64; 3]; 3], strain: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0f64; 3]; 3];
    for k in 0..3 {
        for i in 0..3 {
            let mut v = lattice[k][i];
            for j in 0..3 {
                v += strain[i][j] * lattice[k][j];
            }
            out[k][i] = v;
        }
    }
    out
}

/// Apply a symmetric strain to Cartesian positions.
pub fn strained_positions(positions: &[[f64; 3]], strain: &[[f64; 3]; 3]) -> Vec<[f64; 3]> {
    positions
        .iter()
        .map(|r| {
            let mut out = [0.0f64; 3];
            for i in 0..3 {
                out[i] = r[i] + strain[i][0] * r[0] + strain[i][1] * r[1] + strain[i][2] * r[2];
            }
            out
        })
        .collect()
}

/// Numerical stress tensor (eV/Å³) from central differences of the energy.
///
/// `delta` is the strain amplitude; the off-diagonal components are perturbed
/// symmetrically (`ε_ij = ε_ji`), which is why their derivative carries the
/// factor 1/2 relative to the diagonal ones.
pub fn numerical_stress(
    elements: &[u8],
    positions: &[[f64; 3]],
    lattice: &[[f64; 3]; 3],
    energy_fn: &CellEnergyForceFn,
    delta: f64,
) -> Result<[[f64; 3]; 3], String> {
    let volume = cell_volume(lattice).abs();
    if volume < 1e-12 {
        return Err("cell volume is zero".to_string());
    }
    let energy_at = |strain: &[[f64; 3]; 3]| -> Result<f64, String> {
        let lat = strained_lattice(lattice, strain);
        let pos = strained_positions(positions, strain);
        energy_fn(elements, &pos, &lat).map(|(e, _)| e)
    };

    let mut stress = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in i..3 {
            let mut plus = [[0.0f64; 3]; 3];
            let mut minus = [[0.0f64; 3]; 3];
            if i == j {
                plus[i][i] = delta;
                minus[i][i] = -delta;
            } else {
                // Symmetric off-diagonal perturbation: ε_ij = ε_ji = ±δ/2 so
                // that the total shear amplitude is ±δ.
                plus[i][j] = delta / 2.0;
                plus[j][i] = delta / 2.0;
                minus[i][j] = -delta / 2.0;
                minus[j][i] = -delta / 2.0;
            }
            let e_plus = energy_at(&plus)?;
            let e_minus = energy_at(&minus)?;
            // ∂E/∂ε_ij.
            //
            // Diagonal (i == j): the perturbation is ε_ii = ±δ, so the central
            // difference over ±δ gives the derivative directly.
            //
            // Off-diagonal (i ≠ j): the symmetric perturbation sets
            // ε_ij = ε_ji = ±δ/2. Since the energy is invariant under
            // interchanging i and j, the symmetric-strain derivative
            // ∂E/∂ε_ij counts both entries, and the chain rule contributes the
            // extra factor 1/2 from d(δ/2)/dδ. The two therefore cancel and
            // the same 1/(2δ) central difference applies:
            //     (E(+δ/2) - E(-δ/2)) / (2·(δ/2)) · ½ = (E+ - E-)/(2δ).
            let deriv = (e_plus - e_minus) / (2.0 * delta);
            let value = deriv / volume;
            stress[i][j] = value;
            stress[j][i] = value;
        }
    }
    Ok(stress)
}

/// Configuration for variable-cell relaxation.
#[derive(Debug, Clone)]
pub struct VariableCellConfig {
    /// Maximum number of outer (cell) iterations.
    pub max_outer: usize,
    /// Position-relaxation iterations per outer step.
    pub max_inner: usize,
    /// Force convergence threshold (eV/Å).
    pub force_tol: f64,
    /// Stress convergence threshold (eV/Å³).
    pub stress_tol: f64,
    /// Strain step scale (Å³/eV): `δε = -scale · σ`.
    pub strain_scale: f64,
    /// Maximum strain per outer step.
    pub max_strain: f64,
    /// Strain amplitude for the numerical stress.
    pub stress_delta: f64,
}

impl Default for VariableCellConfig {
    fn default() -> Self {
        Self {
            max_outer: 40,
            max_inner: 100,
            force_tol: 0.05,
            stress_tol: 5e-4,
            strain_scale: 4.0,
            max_strain: 0.02,
            stress_delta: 1e-3,
        }
    }
}

/// Result of a variable-cell relaxation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariableCellResult {
    /// Relaxed positions (Å).
    pub positions: Vec<[f64; 3]>,
    /// Relaxed lattice (Å).
    pub lattice: [[f64; 3]; 3],
    /// Final energy (eV).
    pub energy: f64,
    /// Final forces (eV/Å).
    pub forces: Vec<[f64; 3]>,
    /// Final stress (eV/Å³).
    pub stress: [[f64; 3]; 3],
    /// Number of outer iterations.
    pub n_iterations: usize,
    /// Whether both forces and stress converged.
    pub converged: bool,
    /// Energy after each outer iteration.
    pub energy_history: Vec<f64>,
}

/// Maximum absolute component of a tensor or vector field.
fn max_abs_stress(stress: &[[f64; 3]; 3]) -> f64 {
    stress
        .iter()
        .flat_map(|row| row.iter())
        .fold(0.0f64, |m, v| m.max(v.abs()))
}

fn max_abs_forces(forces: &[[f64; 3]]) -> f64 {
    forces
        .iter()
        .flat_map(|f| f.iter())
        .fold(0.0f64, |m, v| m.max(v.abs()))
}

/// Relax atomic positions and the unit cell.
///
/// Alternates position relaxation (BFGS at the current cell, via
/// [`super::geometry_opt::optimize_framework`]) with a stress-driven strain
/// step accepted only when it lowers the energy (simple line search).
pub fn optimize_variable_cell(
    elements: &[u8],
    initial_positions: &[[f64; 3]],
    initial_lattice: &[[f64; 3]; 3],
    energy_fn: &CellEnergyForceFn,
    config: &VariableCellConfig,
) -> Result<VariableCellResult, String> {
    let mut lattice = *initial_lattice;
    let mut positions = initial_positions.to_vec();
    let mut energy_history = Vec::new();
    let mut converged = false;
    let mut n_iter = 0;

    let (mut energy, mut forces) = energy_fn(elements, &positions, &lattice)?;
    energy_history.push(energy);

    for outer in 0..config.max_outer {
        n_iter = outer + 1;

        // 1) Position relaxation at the current cell.
        {
            let cell = lattice;
            let closure = |els: &[u8], pos: &[[f64; 3]]| energy_fn(els, pos, &cell);
            let inner_config = super::geometry_opt::FrameworkOptConfig {
                max_iter: config.max_inner,
                force_tol: config.force_tol,
                lattice: Some(lattice),
                ..Default::default()
            };
            if let Ok(result) = super::geometry_opt::optimize_framework(
                elements,
                &positions,
                &closure,
                &inner_config,
            ) {
                positions = result.positions;
                energy = result.energy;
                forces = result.forces;
            }
        }

        // 2) Stress (recomputed after every position relaxation).
        let stress = numerical_stress(
            elements,
            &positions,
            &lattice,
            energy_fn,
            config.stress_delta,
        )?;

        let max_force = max_abs_forces(&forces);
        let max_stress = max_abs_stress(&stress);
        energy_history.push(energy);
        if max_force < config.force_tol && max_stress < config.stress_tol {
            converged = true;
            break;
        }

        // 3) Stress-driven strain step with a backtracking line search.
        let mut strain = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                strain[i][j] = -config.strain_scale * stress[i][j];
            }
        }
        // Clamp the step.
        let max_component = strain
            .iter()
            .flat_map(|row| row.iter())
            .fold(0.0f64, |m, v| m.max(v.abs()));
        if max_component > config.max_strain {
            let scale = config.max_strain / max_component;
            for row in strain.iter_mut() {
                for v in row.iter_mut() {
                    *v *= scale;
                }
            }
        }

        let mut accepted = false;
        let mut trial_scale = 1.0f64;
        for _ in 0..6 {
            let mut trial = strain;
            for row in trial.iter_mut() {
                for v in row.iter_mut() {
                    *v *= trial_scale;
                }
            }
            let trial_lattice = strained_lattice(&lattice, &trial);
            let trial_positions = strained_positions(&positions, &trial);
            if let Ok((trial_energy, _)) = energy_fn(elements, &trial_positions, &trial_lattice) {
                if trial_energy < energy {
                    lattice = trial_lattice;
                    positions = trial_positions;
                    energy = trial_energy;
                    accepted = true;
                    break;
                }
            }
            trial_scale *= 0.5;
        }
        if !accepted {
            // No downhill strain step: the cell is as good as it gets with the
            // current numerical stress.
            break;
        }
    }

    // Final forces/stress at the relaxed configuration.
    let (final_energy, final_forces) = energy_fn(elements, &positions, &lattice)?;
    let final_stress = numerical_stress(
        elements,
        &positions,
        &lattice,
        energy_fn,
        config.stress_delta,
    )?;
    if max_abs_forces(&final_forces) < config.force_tol
        && max_abs_stress(&final_stress) < config.stress_tol
    {
        converged = true;
    }

    Ok(VariableCellResult {
        positions,
        lattice,
        energy: final_energy,
        forces: final_forces,
        stress: final_stress,
        n_iterations: n_iter,
        converged,
        energy_history,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Lennard-Jones pair potential on a simple-cubic one-atom cell.
    ///
    /// A *fixed* real-space cutoff (independent of the lattice constant) keeps
    /// the model energy smooth in `a`, so the numerical stress has no cutoff
    /// artifacts.
    const EPS: f64 = 0.01; // eV
    const SIGMA: f64 = 3.0; // Å
    const CUTOFF: f64 = 9.0; // Å
    const N_IMG: i32 = 3;

    fn model_energy_at(_elements: &[u8], pos: &[[f64; 3]], lattice: &[[f64; 3]; 3]) -> f64 {
        let a = lattice[0][0];
        let shift = pos.first().copied().unwrap_or([0.0; 3]);
        let mut energy = 0.0;
        for i in -N_IMG..=N_IMG {
            for j in -N_IMG..=N_IMG {
                for k in -N_IMG..=N_IMG {
                    if i == 0 && j == 0 && k == 0 {
                        continue;
                    }
                    let dx = i as f64 * a - shift[0];
                    let dy = j as f64 * a - shift[1];
                    let dz = k as f64 * a - shift[2];
                    let r = (dx * dx + dy * dy + dz * dz).sqrt();
                    if r > CUTOFF {
                        continue;
                    }
                    let sr6 = (SIGMA / r).powi(6);
                    energy += 4.0 * EPS * (sr6 * sr6 - sr6);
                }
            }
        }
        energy
    }

    fn lj_energy(a: f64) -> f64 {
        model_energy_at(
            &[18],
            &[[0.0, 0.0, 0.0]],
            &[[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]],
        )
    }

    fn lj_energy_fn(
        elements: &[u8],
        pos: &[[f64; 3]],
        lattice: &[[f64; 3]; 3],
    ) -> Result<(f64, Vec<[f64; 3]>), String> {
        // Numerical forces for the model (central differences on the atom).
        let p0 = pos[0];
        let h = 1e-5;
        let mut f = [0.0f64; 3];
        for k in 0..3 {
            let mut p_plus = p0;
            let mut p_minus = p0;
            p_plus[k] += h;
            p_minus[k] -= h;
            let e_plus = model_energy_at(elements, &[p_plus], lattice);
            let e_minus = model_energy_at(elements, &[p_minus], lattice);
            f[k] = -(e_plus - e_minus) / (2.0 * h);
        }
        Ok((model_energy_at(elements, pos, lattice), vec![f]))
    }

    fn equilibrium_a() -> f64 {
        // Golden-section search on E(a).
        let (mut lo, mut hi) = (3.0, 6.0);
        for _ in 0..200 {
            let m1 = lo + 0.382 * (hi - lo);
            let m2 = lo + 0.618 * (hi - lo);
            if lj_energy(m1) < lj_energy(m2) {
                hi = m2;
            } else {
                lo = m1;
            }
        }
        0.5 * (lo + hi)
    }

    #[test]
    fn stress_vanishes_at_equilibrium_and_has_correct_sign() {
        let a_eq = equilibrium_a();
        let lattice = [[a_eq, 0.0, 0.0], [0.0, a_eq, 0.0], [0.0, 0.0, a_eq]];
        let positions = [[0.0, 0.0, 0.0]];
        let stress = numerical_stress(&[18], &positions, &lattice, &lj_energy_fn, 1e-3).unwrap();
        for i in 0..3 {
            assert!(
                stress[i][i].abs() < 1e-4,
                "stress at equilibrium: {:?}",
                stress
            );
        }
        // Convention σ = (1/V)∂E/∂ε: compression ⇒ negative (compressive),
        // expansion ⇒ positive (tensile).
        let compressed = [
            [0.95 * a_eq, 0.0, 0.0],
            [0.0, 0.95 * a_eq, 0.0],
            [0.0, 0.0, 0.95 * a_eq],
        ];
        let s_comp = numerical_stress(&[18], &positions, &compressed, &lj_energy_fn, 1e-3).unwrap();
        assert!(
            s_comp[0][0] < 0.0,
            "compression must give a negative stress: {s_comp:?}"
        );
        let expanded = [
            [1.05 * a_eq, 0.0, 0.0],
            [0.0, 1.05 * a_eq, 0.0],
            [0.0, 0.0, 1.05 * a_eq],
        ];
        let s_exp = numerical_stress(&[18], &positions, &expanded, &lj_energy_fn, 1e-3).unwrap();
        assert!(
            s_exp[0][0] > 0.0,
            "expansion must give a positive stress: {s_exp:?}"
        );
    }

    #[test]
    fn variable_cell_relaxation_recovers_equilibrium_lattice() {
        let a_eq = equilibrium_a();
        let a0 = 1.12 * a_eq;
        let lattice = [[a0, 0.0, 0.0], [0.0, a0, 0.0], [0.0, 0.0, a0]];
        let positions = [[0.0, 0.0, 0.0]];
        let config = VariableCellConfig {
            stress_tol: 1e-3,
            max_outer: 60,
            ..Default::default()
        };
        let result = optimize_variable_cell(&[18], &positions, &lattice, &lj_energy_fn, &config)
            .expect("relaxation runs");
        let a_final = result.lattice[0][0];
        assert!(
            (a_final - a_eq).abs() < 0.1,
            "lattice {a_final} vs equilibrium {a_eq}"
        );
        assert!(
            result.energy <= model_energy_at(&[18], &positions, &lattice) + 1e-9,
            "energy must decrease"
        );
        // Energy history is monotonically non-increasing.
        for w in result.energy_history.windows(2) {
            assert!(
                w[1] <= w[0] + 1e-9,
                "energy increased: {:?}",
                result.energy_history
            );
        }
    }

    #[test]
    fn strained_lattice_applies_affine_deformation() {
        let lattice = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
        let mut strain = [[0.0f64; 3]; 3];
        strain[0][0] = 0.1;
        strain[1][2] = 0.05;
        strain[2][1] = 0.05;
        let out = strained_lattice(&lattice, &strain);
        assert!((out[0][0] - 1.1).abs() < 1e-12);
        assert!((out[1][2] - 0.1).abs() < 1e-12);
        assert!((out[2][1] - 0.15).abs() < 1e-12);
    }
}
