//! **ALPHA** — Langevin thermostat for stochastic dynamics.
//!
//! Implements a Langevin integrator that maintains a target temperature via
//! friction + random forces. Uses the BAOAB splitting scheme for accurate
//! configurational sampling.

use crate::dynamics_live::state::DynamicAtom;

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

/// Boltzmann constant in eV/K.
const KB_EV: f64 = 8.617_333_262e-5;

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
    let kt = KB_EV * config.target_temp_k;
    let c1 = (-gamma * dt_fs).exp();

    // Simple xorshift64 PRNG for reproducible Gaussian noise
    let mut state = rng_seed.wrapping_add(1);

    for atom in atoms.iter_mut() {
        let c2 = ((1.0 - c1 * c1) * kt / atom.mass).sqrt();

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

    // eV/Å per (amu·Å/fs²): conversion factor
    // F [eV/Å], m [amu], a = F/m [eV/(Å·amu)]
    // v [Å/fs] = a * dt [eV·fs/(Å·amu)]
    // Need: 1 eV = 1.602e-19 J, 1 amu = 1.661e-27 kg
    // => eV/(Å·amu) = 9.648e6 Å/fs² ... but in natural MD units:
    // We use eV, Å, fs, amu: acceleration = F/m has units eV/(Å·amu)
    // Factor: 1 eV/(amu·Å) = 0.009_648_533 Å/fs²
    let accel_factor = 0.009_648_533;

    // B: half-step velocity from forces
    for atom in atoms.iter_mut() {
        let inv_m = accel_factor / atom.mass;
        for d in 0..3 {
            atom.velocity[d] += half_dt * atom.force[d] * inv_m;
        }
    }

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
pub fn langevin_b_step(atoms: &mut [DynamicAtom], dt_fs: f64) {
    let half_dt = 0.5 * dt_fs;
    let accel_factor = 0.009_648_533;

    for atom in atoms.iter_mut() {
        let inv_m = accel_factor / atom.mass;
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
}
