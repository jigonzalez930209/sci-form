//! Dynamic atom and molecular system state for real-time simulation.
//!
//! Provides a contiguous flat buffer (`positions_flat`) that can be shared
//! directly with JavaScript via WASM pointer — zero serialization cost per frame.

use serde::{Deserialize, Serialize};

use crate::dynamics::{atomic_mass_amu, MdBackend};

/// Single atom carrying position, velocity, and force vectors for time-stepping.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DynamicAtom {
    /// Atomic number (1=H, 6=C, 7=N, 8=O, …).
    pub element: u8,
    /// Atomic mass in amu.
    pub mass: f64,
    /// Partial charge (e).
    pub charge: f64,
    /// Position in Å.
    pub position: [f64; 3],
    /// Velocity in Å/fs.
    pub velocity: [f64; 3],
    /// Force in kcal/(mol·Å).
    pub force: [f64; 3],
}

/// Full simulation state supporting zero-copy WASM export.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LiveMolecularSystem {
    /// Per-atom dynamic state.
    pub atoms: Vec<DynamicAtom>,
    /// Potential energy (kcal/mol).
    pub potential_energy: f64,
    /// Kinetic energy (kcal/mol).
    pub kinetic_energy: f64,
    /// Simulation time (fs).
    pub time_fs: f64,
    /// Integration step counter.
    pub step: usize,
    /// Temperature (K).
    pub temperature_k: f64,
    /// Force backend in use.
    pub backend: MdBackend,
    /// Flat position buffer: [x0,y0,z0, x1,y1,z1, …] — exposed to JS via pointer.
    pub positions_flat: Vec<f64>,
    /// Flat velocity buffer.
    pub velocities_flat: Vec<f64>,
    /// Flat force buffer.
    pub forces_flat: Vec<f64>,
    /// SMILES topology (needed for UFF).
    pub smiles: String,
    /// Bond list for force field evaluation.
    pub bonds: Vec<(usize, usize, String)>,
}

/// Conversion constants.
const AMU_ANGFS2_TO_KCAL_MOL: f64 = 2_390.057_361_533_49;
const R_GAS_KCAL_MOLK: f64 = 0.001_987_204_258_640_83;

impl LiveMolecularSystem {
    /// Initialize from a conformer result.
    pub fn from_conformer(
        smiles: &str,
        elements: &[u8],
        coords: &[f64],
        bonds: &[(usize, usize, String)],
        backend: MdBackend,
    ) -> Self {
        let n = elements.len();
        let mut atoms = Vec::with_capacity(n);
        for i in 0..n {
            atoms.push(DynamicAtom {
                element: elements[i],
                mass: atomic_mass_amu(elements[i]),
                charge: 0.0,
                position: [coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]],
                velocity: [0.0; 3],
                force: [0.0; 3],
            });
        }

        Self {
            atoms,
            potential_energy: 0.0,
            kinetic_energy: 0.0,
            time_fs: 0.0,
            step: 0,
            temperature_k: 0.0,
            backend,
            positions_flat: coords.to_vec(),
            velocities_flat: vec![0.0; n * 3],
            forces_flat: vec![0.0; n * 3],
            smiles: smiles.to_string(),
            bonds: bonds.to_vec(),
        }
    }

    /// Number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// Sync structured atoms → flat buffer for WASM pointer export.
    pub fn sync_flat_from_atoms(&mut self) {
        let n = self.atoms.len();
        self.positions_flat.resize(n * 3, 0.0);
        self.velocities_flat.resize(n * 3, 0.0);
        self.forces_flat.resize(n * 3, 0.0);
        for (i, atom) in self.atoms.iter().enumerate() {
            self.positions_flat[3 * i] = atom.position[0];
            self.positions_flat[3 * i + 1] = atom.position[1];
            self.positions_flat[3 * i + 2] = atom.position[2];
            self.velocities_flat[3 * i] = atom.velocity[0];
            self.velocities_flat[3 * i + 1] = atom.velocity[1];
            self.velocities_flat[3 * i + 2] = atom.velocity[2];
            self.forces_flat[3 * i] = atom.force[0];
            self.forces_flat[3 * i + 1] = atom.force[1];
            self.forces_flat[3 * i + 2] = atom.force[2];
        }
    }

    /// Sync flat buffers → structured atoms (after external modification).
    pub fn sync_atoms_from_flat(&mut self) {
        for (i, atom) in self.atoms.iter_mut().enumerate() {
            atom.position[0] = self.positions_flat[3 * i];
            atom.position[1] = self.positions_flat[3 * i + 1];
            atom.position[2] = self.positions_flat[3 * i + 2];
            atom.velocity[0] = self.velocities_flat[3 * i];
            atom.velocity[1] = self.velocities_flat[3 * i + 1];
            atom.velocity[2] = self.velocities_flat[3 * i + 2];
            atom.force[0] = self.forces_flat[3 * i];
            atom.force[1] = self.forces_flat[3 * i + 1];
            atom.force[2] = self.forces_flat[3 * i + 2];
        }
    }

    /// Initialize random velocities at a target temperature using Box-Muller.
    pub fn initialize_velocities(&mut self, temperature_k: f64, seed: u64) {
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};

        let mut rng = StdRng::seed_from_u64(seed);
        let n = self.atoms.len();

        for atom in &mut self.atoms {
            let sigma =
                (R_GAS_KCAL_MOLK * temperature_k / (atom.mass * AMU_ANGFS2_TO_KCAL_MOL)).sqrt();
            for d in 0..3 {
                // Box-Muller transform
                let u1 = (1.0 - rng.gen::<f64>()).max(1e-12);
                let u2 = rng.gen::<f64>();
                atom.velocity[d] =
                    sigma * (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            }
        }

        // Remove center-of-mass velocity
        let mut com_v = [0.0f64; 3];
        let mut total_mass = 0.0;
        for atom in &self.atoms {
            for d in 0..3 {
                com_v[d] += atom.mass * atom.velocity[d];
            }
            total_mass += atom.mass;
        }
        if total_mass > 0.0 {
            for d in 0..3 {
                com_v[d] /= total_mass;
            }
            for atom in &mut self.atoms {
                for d in 0..3 {
                    atom.velocity[d] -= com_v[d];
                }
            }
        }

        // Rescale to exact target temperature
        let (_ke, temp) = self.kinetic_energy_and_temperature();
        if temp > 1e-10 {
            let scale = (temperature_k / temp).sqrt();
            for atom in &mut self.atoms {
                for d in 0..3 {
                    atom.velocity[d] *= scale;
                }
            }
        }

        let (ke, temp) = self.kinetic_energy_and_temperature();
        self.kinetic_energy = ke;
        self.temperature_k = temp;
        self.sync_flat_from_atoms();

        let _ = n; // suppress unused
    }

    /// Compute kinetic energy (kcal/mol) and temperature (K).
    pub fn kinetic_energy_and_temperature(&self) -> (f64, f64) {
        let n = self.atoms.len();
        let mut ke = 0.0;
        for atom in &self.atoms {
            let v2 = atom.velocity[0].powi(2) + atom.velocity[1].powi(2) + atom.velocity[2].powi(2);
            ke += 0.5 * atom.mass * v2 * AMU_ANGFS2_TO_KCAL_MOL;
        }
        let dof = (3 * n).saturating_sub(6).max(1) as f64;
        let t = 2.0 * ke / (dof * R_GAS_KCAL_MOLK);
        (ke, t)
    }

    /// Compute forces using the configured backend.
    pub fn compute_forces(&mut self) -> Result<f64, String> {
        let elements: Vec<u8> = self.atoms.iter().map(|a| a.element).collect();
        let energy = crate::dynamics::compute_backend_energy_and_gradients(
            self.backend,
            &self.smiles,
            &elements,
            &self.positions_flat,
            &mut self.forces_flat,
        )?;

        // Negate gradients → forces
        for f in &mut self.forces_flat {
            *f = -*f;
        }

        // Sync to atoms
        for (i, atom) in self.atoms.iter_mut().enumerate() {
            atom.force[0] = self.forces_flat[3 * i];
            atom.force[1] = self.forces_flat[3 * i + 1];
            atom.force[2] = self.forces_flat[3 * i + 2];
        }

        Ok(energy)
    }

    /// Run one Velocity Verlet integration step.
    pub fn verlet_step(&mut self, dt_fs: f64) -> Result<(), String> {
        let n = self.n_atoms();

        // Half-step velocity update: v(t + dt/2) = v(t) + 0.5*dt*a(t)
        for atom in &mut self.atoms {
            let inv_m = 1.0 / (atom.mass * AMU_ANGFS2_TO_KCAL_MOL);
            for d in 0..3 {
                atom.velocity[d] += 0.5 * dt_fs * atom.force[d] * inv_m;
            }
        }

        // Position update: r(t + dt) = r(t) + dt*v(t + dt/2)
        for atom in &mut self.atoms {
            for d in 0..3 {
                atom.position[d] += dt_fs * atom.velocity[d];
            }
        }

        self.sync_flat_from_atoms();

        // Compute new forces at new positions
        let pe = self.compute_forces()?;
        self.potential_energy = pe;

        // Second half-step velocity: v(t + dt) = v(t + dt/2) + 0.5*dt*a(t+dt)
        for atom in &mut self.atoms {
            let inv_m = 1.0 / (atom.mass * AMU_ANGFS2_TO_KCAL_MOL);
            for d in 0..3 {
                atom.velocity[d] += 0.5 * dt_fs * atom.force[d] * inv_m;
            }
        }

        let (ke, temp) = self.kinetic_energy_and_temperature();
        self.kinetic_energy = ke;
        self.temperature_k = temp;
        self.time_fs += dt_fs;
        self.step += 1;

        self.sync_flat_from_atoms();

        // Divergence check
        for i in 0..n {
            if !self.positions_flat[3 * i].is_finite() {
                return Err("simulation diverged: non-finite coordinates".to_string());
            }
        }

        Ok(())
    }

    /// Run N substeps of Velocity Verlet, returning the updated state.
    pub fn integrate(&mut self, dt_fs: f64, substeps: usize) -> Result<(), String> {
        for _ in 0..substeps {
            self.verlet_step(dt_fs)?;
        }
        Ok(())
    }

    /// Apply Berendsen thermostat velocity scaling.
    pub fn berendsen_thermostat(&mut self, target_temp_k: f64, tau_fs: f64, dt_fs: f64) {
        if self.temperature_k < 1e-10 {
            return;
        }
        let lambda = (1.0 + (dt_fs / tau_fs) * (target_temp_k / self.temperature_k - 1.0))
            .sqrt()
            .clamp(0.5, 2.0);
        for atom in &mut self.atoms {
            for d in 0..3 {
                atom.velocity[d] *= lambda;
            }
        }
        let (ke, temp) = self.kinetic_energy_and_temperature();
        self.kinetic_energy = ke;
        self.temperature_k = temp;
        self.sync_flat_from_atoms();
    }

    /// Apply Nosé-Hoover chain thermostat (single chain variable for simplicity).
    pub fn nose_hoover_thermostat(
        &mut self,
        target_temp_k: f64,
        thermostat_mass: f64,
        xi: &mut f64,
        v_xi: &mut f64,
        dt_fs: f64,
    ) {
        let n = self.n_atoms();
        let dof = (3 * n).saturating_sub(6).max(1) as f64;
        let target_ke = 0.5 * dof * R_GAS_KCAL_MOLK * target_temp_k;

        // Update thermostat velocity
        *v_xi += (self.kinetic_energy - target_ke) / thermostat_mass * dt_fs;
        *xi += *v_xi * dt_fs;

        // Scale velocities
        let scale = (-(*v_xi) * dt_fs).exp();
        for atom in &mut self.atoms {
            for d in 0..3 {
                atom.velocity[d] *= scale;
            }
        }

        let (ke, temp) = self.kinetic_energy_and_temperature();
        self.kinetic_energy = ke;
        self.temperature_k = temp;
        self.sync_flat_from_atoms();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn h2_system() -> LiveMolecularSystem {
        // H2 molecule at ~0.74 Å
        let elements = [1u8, 1];
        let coords = [0.0, 0.0, 0.0, 0.74, 0.0, 0.0];
        let bonds = vec![(0, 1, "SINGLE".to_string())];
        LiveMolecularSystem::from_conformer("", &elements, &coords, &bonds, MdBackend::Uff)
    }

    #[test]
    fn from_conformer_creates_atoms() {
        let sys = h2_system();
        assert_eq!(sys.n_atoms(), 2);
        assert_eq!(sys.atoms[0].element, 1);
        assert_eq!(sys.atoms[1].element, 1);
        assert!(sys.atoms[0].mass > 0.9 && sys.atoms[0].mass < 1.1);
    }

    #[test]
    fn sync_flat_roundtrip() {
        let mut sys = h2_system();
        sys.atoms[0].position = [1.0, 2.0, 3.0];
        sys.sync_flat_from_atoms();
        assert!((sys.positions_flat[0] - 1.0).abs() < 1e-12);
        assert!((sys.positions_flat[1] - 2.0).abs() < 1e-12);
        assert!((sys.positions_flat[2] - 3.0).abs() < 1e-12);

        sys.positions_flat[3] = 5.0;
        sys.sync_atoms_from_flat();
        assert!((sys.atoms[1].position[0] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn initialize_velocities_nonzero() {
        let mut sys = h2_system();
        sys.initialize_velocities(300.0, 42);
        let v_mag: f64 = sys
            .atoms
            .iter()
            .map(|a| a.velocity.iter().map(|v| v * v).sum::<f64>())
            .sum();
        assert!(v_mag > 0.0, "velocities should be nonzero after init");
    }

    #[test]
    fn kinetic_energy_zero_for_stationary() {
        let sys = h2_system();
        let (ke, temp) = sys.kinetic_energy_and_temperature();
        assert!(ke.abs() < 1e-15, "KE should be zero for stationary atoms");
        assert!(temp.abs() < 1e-15);
    }

    #[test]
    fn berendsen_thermostat_scales_velocities() {
        let mut sys = h2_system();
        sys.initialize_velocities(300.0, 42);
        let (ke_before, _) = sys.kinetic_energy_and_temperature();
        sys.berendsen_thermostat(600.0, 100.0, 1.0);
        let (ke_after, _) = sys.kinetic_energy_and_temperature();
        // Thermostat at double temp should increase KE
        assert!(ke_after > ke_before);
    }

    #[test]
    fn n_atoms_matches_flat_buffer_length() {
        let sys = h2_system();
        assert_eq!(sys.positions_flat.len(), 3 * sys.n_atoms());
        assert_eq!(sys.velocities_flat.len(), 3 * sys.n_atoms());
        assert_eq!(sys.forces_flat.len(), 3 * sys.n_atoms());
    }
}
