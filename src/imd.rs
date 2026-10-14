//! Interactive Molecular Dynamics (IMD) steering forces.
//!
//! Provides harmonic spring potentials for interactive atom manipulation.

use serde::{Deserialize, Serialize};

/// A harmonic steering force applied to a single atom.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SteeringForce {
    /// Target atom index.
    pub atom_index: usize,
    /// Target position in Å.
    pub target_xyz: [f64; 3],
    /// Spring constant in kcal/(mol·Å²).
    pub spring_k: f64,
}

/// Collection of active IMD steering forces.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct SteeringForces {
    /// Active steering potentials.
    pub forces: Vec<SteeringForce>,
}

impl SteeringForces {
    /// Create empty steering force collection.
    pub fn new() -> Self {
        Self { forces: Vec::new() }
    }

    /// Add or update a steering force on an atom.
    pub fn apply(&mut self, atom_index: usize, target_xyz: [f64; 3], spring_k: f64) {
        // Remove existing force on this atom
        self.forces.retain(|f| f.atom_index != atom_index);
        self.forces.push(SteeringForce {
            atom_index,
            target_xyz,
            spring_k,
        });
    }

    /// Remove steering force from an atom.
    pub fn clear(&mut self, atom_index: usize) {
        self.forces.retain(|f| f.atom_index != atom_index);
    }

    /// Remove all steering forces.
    pub fn clear_all(&mut self) {
        self.forces.clear();
    }

    /// Whether any steering forces are active.
    pub fn is_active(&self) -> bool {
        !self.forces.is_empty()
    }

    /// Accumulate IMD steering energy and forces into flat buffer.
    ///
    /// Returns total steering energy (kcal/mol).
    pub fn accumulate(&self, positions_flat: &[f64], forces_flat: &mut [f64]) -> f64 {
        let mut total_energy = 0.0;

        for sf in &self.forces {
            let i = sf.atom_index;
            let dx = positions_flat[3 * i] - sf.target_xyz[0];
            let dy = positions_flat[3 * i + 1] - sf.target_xyz[1];
            let dz = positions_flat[3 * i + 2] - sf.target_xyz[2];

            // E = 0.5 * k * |r - r_target|²
            let d2 = dx * dx + dy * dy + dz * dz;
            total_energy += 0.5 * sf.spring_k * d2;

            // F = -k * (r - r_target)
            forces_flat[3 * i] -= sf.spring_k * dx;
            forces_flat[3 * i + 1] -= sf.spring_k * dy;
            forces_flat[3 * i + 2] -= sf.spring_k * dz;
        }

        total_energy
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn apply_and_clear() {
        let mut sf = SteeringForces::new();
        assert!(!sf.is_active());
        sf.apply(0, [1.0, 0.0, 0.0], 100.0);
        assert!(sf.is_active());
        sf.clear(0);
        assert!(!sf.is_active());
    }

    #[test]
    fn accumulate_force_toward_target() {
        let mut sf = SteeringForces::new();
        sf.apply(0, [1.0, 0.0, 0.0], 100.0);
        // Atom at origin, target at (1,0,0)
        let positions = [0.0, 0.0, 0.0, 5.0, 0.0, 0.0];
        let mut forces = [0.0; 6];
        let energy = sf.accumulate(&positions, &mut forces);
        // Force should pull atom 0 toward +x
        assert!(forces[0] > 0.0, "force should be in +x direction");
        assert!(energy > 0.0, "energy should be positive");
    }

    #[test]
    fn clear_all_removes_everything() {
        let mut sf = SteeringForces::new();
        sf.apply(0, [1.0, 0.0, 0.0], 100.0);
        sf.apply(1, [0.0, 1.0, 0.0], 50.0);
        sf.clear_all();
        assert!(!sf.is_active());
    }
}
