//! Berendsen barostat for NPT ensemble.

use super::pbc::SimulationBox;
use super::state::LiveMolecularSystem;

/// Configuration for Berendsen pressure coupling.
#[derive(Debug, Clone)]
pub struct BerendsenBarostat {
    /// Target pressure (bar).
    pub target_pressure_bar: f64,
    /// Pressure coupling time constant (fs).
    pub tau_p_fs: f64,
    /// Isothermal compressibility (bar⁻¹). Water: ~4.5e-5.
    pub compressibility: f64,
}

impl BerendsenBarostat {
    /// Create a new Berendsen barostat.
    pub fn new(target_pressure_bar: f64, tau_p_fs: f64, compressibility: f64) -> Self {
        Self {
            target_pressure_bar,
            tau_p_fs,
            compressibility,
        }
    }

    /// Apply isotropic Berendsen pressure coupling.
    ///
    /// This rescales coordinates and the simulation box uniformly.
    ///
    /// Requires the instantaneous pressure, which is computed from the virial:
    ///   P = (N*k_B*T + virial) / (3*V)
    pub fn apply(
        &self,
        system: &mut LiveMolecularSystem,
        sim_box: &mut SimulationBox,
        instantaneous_pressure_bar: f64,
        dt_fs: f64,
    ) {
        if !sim_box.periodic {
            return;
        }

        let dp = instantaneous_pressure_bar - self.target_pressure_bar;
        let mu = (1.0 - (self.compressibility * dt_fs / (3.0 * self.tau_p_fs)) * dp).cbrt();

        // Clamp scaling factor to prevent runaway
        let mu = mu.clamp(0.95, 1.05);

        // Scale positions
        for atom in &mut system.atoms {
            for d in 0..3 {
                atom.position[d] *= mu;
            }
        }

        // Scale box vectors
        for i in 0..3 {
            for j in 0..3 {
                sim_box.lattice[i][j] *= mu;
            }
        }

        // Recompute inverse lattice
        *sim_box =
            SimulationBox::triclinic(sim_box.lattice[0], sim_box.lattice[1], sim_box.lattice[2]);

        system.sync_flat_from_atoms();
    }
}

/// Compute instantaneous pressure from kinetic energy and virial.
///
/// P = (2*KE + virial) / (3*V)
///
/// `virial` is sum over all pairs: r_ij · f_ij (in kcal/mol·Å² units).
/// Returns pressure in bar.
pub fn compute_pressure_bar(
    kinetic_energy_kcal_mol: f64,
    virial_kcal_mol: f64,
    volume_ang3: f64,
) -> f64 {
    // 1 kcal/mol / ų = 6.9477e4 bar
    const KCAL_MOL_ANG3_TO_BAR: f64 = 6.9477e4;
    (2.0 * kinetic_energy_kcal_mol + virial_kcal_mol) / (3.0 * volume_ang3) * KCAL_MOL_ANG3_TO_BAR
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compute_pressure_bar_positive() {
        // Ideal gas: PV = nRT => P = 2*KE / (3V) in appropriate units
        let ke = 10.0; // kcal/mol
        let virial = 0.0;
        let vol = 1000.0; // ų
        let p = compute_pressure_bar(ke, virial, vol);
        assert!(p > 0.0, "pressure should be positive");
    }

    #[test]
    fn barostat_scaling() {
        let barostat = BerendsenBarostat::new(1.0, 1000.0, 4.5e-5);
        assert!((barostat.target_pressure_bar - 1.0).abs() < 1e-10);
        assert!((barostat.compressibility - 4.5e-5).abs() < 1e-15);
    }
}
