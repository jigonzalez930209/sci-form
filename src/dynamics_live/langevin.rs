//! **ALPHA** — Langevin thermostat for stochastic dynamics.
//!
//! Implements a Langevin integrator that maintains a target temperature via
//! friction + random forces. Uses the BAOAB splitting scheme for accurate
//! configurational sampling.
//!
//! # Unit system
//!
//! The module shares the unit system of [`crate::dynamics_live::state`]:
//! positions in Å, velocities in Å/fs, forces in kcal/(mol·Å) and masses in
//! amu. The thermal energy therefore uses the gas constant in kcal/(mol·K)
//! and the velocity variance uses the same `amu·(Å/fs)² → kcal/mol` factor
//! as the Verlet integrator.

use crate::dynamics_live::state::{DynamicAtom, AMU_ANGFS2_TO_KCAL_MOL, R_GAS_KCAL_MOLK};

/// Langevin thermostat configuration.
#[derive(Debug, Clone)]
pub struct LangevinConfig {
    /// Target temperature (K).
    pub target_temp_k: f64,
    /// Friction coefficient γ (1/fs). Typical: 0.001–0.01 fs⁻¹.
    pub gamma: f64,
}

impl Default for LangevinConfig {
    fn default() -> Self {
        Self {
            target_temp_k: 300.0,
            gamma: 0.001,
        }
    }
}

/// Apply Langevin velocity update (O step of BAOAB):
///
///   v ← c₁ v + c₂ R
///
/// where c₁ = exp(-γ Δt), c₂ = √((1 - c₁²) kT / m), R ~ N(0,1).
pub fn langevin_o_step(
    atoms: &mut [DynamicAtom],
    config: &LangevinConfig,
    dt_fs: f64,
    rng_seed: u64,
) {
    let gamma = config.gamma;
    // Thermal energy in kcal/mol.
    let kt = R_GAS_KCAL_MOLK * config.target_temp_k;
    let c1 = (-gamma * dt_fs).exp();

    // Simple xorshift64 PRNG for reproducible Gaussian noise
    let mut state = rng_seed.wrapping_add(1);

    for atom in atoms.iter_mut() {
        // Variance of one velocity component: kT/m, converted from
        // (kcal/mol)/amu to (Å/fs)² via 1 amu·(Å/fs)² = AMU_ANGFS2_TO_KCAL_MOL.
        let c2 = ((1.0 - c1 * c1) * kt / (atom.mass * AMU_ANGFS2_TO_KCAL_MOL)).sqrt();

        for d in 0..3 {
            let r = box_muller_normal(&mut state);
            atom.velocity[d] = c1 * atom.velocity[d] + c2 * r;
        }
    }
}

/// Full BAOAB Langevin integration step.
///
/// B: half-step velocity update from forces
/// A: half-step position update
/// O: Langevin thermostat velocity randomization
/// A: half-step position update
/// B: half-step velocity update from forces (after new force evaluation)
///
/// Note: the caller must evaluate forces between the two B steps.
/// This function performs B-A-O-A (the first part), returns, then the caller
/// evaluates forces, and calls `langevin_b_step` for the final B.
pub fn langevin_baoa_step(
    atoms: &mut [DynamicAtom],
    config: &LangevinConfig,
    dt_fs: f64,
    rng_seed: u64,
) {
    let half_dt = 0.5 * dt_fs;

    // B: half-step velocity from forces
    langevin_b_step(atoms, dt_fs);

    // A: half-step position
    for atom in atoms.iter_mut() {
        for d in 0..3 {
            atom.position[d] += half_dt * atom.velocity[d];
        }
    }

    // O: Langevin thermostat
    langevin_o_step(atoms, config, dt_fs, rng_seed);

    // A: half-step position
    for atom in atoms.iter_mut() {
        for d in 0..3 {
            atom.position[d] += half_dt * atom.velocity[d];
        }
    }
}

/// Final B step: half-step velocity update from newly computed forces.
///
/// The force is given in kcal/(mol·Å); with `m` in amu the acceleration in
/// Å/fs² is `F / (m · AMU_ANGFS2_TO_KCAL_MOL)`.
pub fn langevin_b_step(atoms: &mut [DynamicAtom], dt_fs: f64) {
    let half_dt = 0.5 * dt_fs;

    for atom in atoms.iter_mut() {
        let inv_m = 1.0 / (atom.mass * AMU_ANGFS2_TO_KCAL_MOL);
        for d in 0..3 {
            atom.velocity[d] += half_dt * atom.force[d] * inv_m;
        }
    }
}

/// Box-Muller transform: generate a normal-distributed random number from a xorshift64 PRNG.
fn box_muller_normal(state: &mut u64) -> f64 {
    let u1 = xorshift64_uniform(state);
    let u2 = xorshift64_uniform(state);
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}

/// xorshift64 PRNG returning a uniform f64 in (0, 1).
fn xorshift64_uniform(state: &mut u64) -> f64 {
    let mut x = *state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *state = x;
    // Map to (0, 1)
    (x as f64) / (u64::MAX as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dynamics_live::state::DynamicAtom;

    fn hydrogen_atom() -> DynamicAtom {
        DynamicAtom {
            element: 1,
            mass: 1.008,
            charge: 0.0,
            position: [0.0, 0.0, 0.0],
            velocity: [0.01, -0.01, 0.005],
            force: [0.0, 0.0, 0.0],
        }
    }

    #[test]
    fn test_langevin_o_step_thermalizes() {
        let mut atoms = vec![hydrogen_atom(); 100];
        let config = LangevinConfig {
            target_temp_k: 300.0,
            gamma: 0.01,
        };

        // Run many O steps
        for step in 0..1000 {
            langevin_o_step(&mut atoms, &config, 1.0, step);
        }

        // Check that velocities are not all zero
        let v_sq_sum: f64 = atoms
            .iter()
            .map(|a| a.velocity[0].powi(2) + a.velocity[1].powi(2) + a.velocity[2].powi(2))
            .sum();
        assert!(v_sq_sum > 0.0);
    }

    #[test]
    fn test_langevin_reaches_target_temperature() {
        // The O step must sample the Maxwell-Boltzmann distribution at the
        // target temperature: <v²> = 3 kT / m (in Å²/fs²).
        let target_k = 300.0;
        let n = 4000;
        let mut atoms = vec![hydrogen_atom(); n];
        let config = LangevinConfig {
            target_temp_k: target_k,
            gamma: 0.05,
        };

        // Large enough step that the velocity forgets the initial value.
        for step in 0..2000 {
            langevin_o_step(&mut atoms, &config, 1.0, step);
        }

        let mean_v2: f64 = atoms
            .iter()
            .map(|a| a.velocity[0].powi(2) + a.velocity[1].powi(2) + a.velocity[2].powi(2))
            .sum::<f64>()
            / (n as f64);

        let expected = 3.0 * R_GAS_KCAL_MOLK * target_k / (atoms[0].mass * AMU_ANGFS2_TO_KCAL_MOL);
        let rel_err = (mean_v2 - expected).abs() / expected;
        assert!(
            rel_err < 0.05,
            "velocity variance {mean_v2:.6e} Å²/fs² deviates from the Maxwell-Boltzmann \
             expectation {expected:.6e} Å²/fs² by {:.2}%",
            rel_err * 100.0
        );
    }

    #[test]
    fn test_langevin_b_step_uses_kcal_unit_system() {
        // A 1 kcal/(mol·Å) force on a 1 amu atom must give the acceleration
        // implied by 1 amu·(Å/fs)² = 2390.057 kcal/mol.
        let mut atoms = vec![DynamicAtom {
            element: 1,
            mass: 1.0,
            charge: 0.0,
            position: [0.0; 3],
            velocity: [0.0; 3],
            force: [1.0, 0.0, 0.0],
        }];
        let dt = 1.0;
        langevin_b_step(&mut atoms, dt);
        let expected_v = 0.5 * dt * 1.0 / (1.0 * AMU_ANGFS2_TO_KCAL_MOL);
        assert!(
            (atoms[0].velocity[0] - expected_v).abs() < 1e-15,
            "velocity {} != expected {expected_v}",
            atoms[0].velocity[0]
        );
    }
}
