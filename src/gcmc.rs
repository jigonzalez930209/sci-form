//! Grand-canonical Monte Carlo adsorption on a periodic framework.
//!
//! The first implementation is a Lennard-Jones GCMC engine on a rigid
//! framework. It is deliberately smaller than a production MOF adsorption
//! code: there are no framework charges or Ewald sums yet, but the
//! insertion/deletion acceptance rules, periodic minimum image, guest-guest
//! interactions and isotherm averaging are the real GCMC algorithm.
//!
//! # Units
//!
//! * positions and cell lengths in **Å**;
//! * Lennard-Jones `epsilon` in **kcal/mol** and `sigma` in **Å**;
//! * pressure in **bar**;
//! * temperature in **K**;
//! * loading in **molecules per unit cell**.
//!
//! # Reference
//!
//! Frenkel & Smit, *Understanding Molecular Simulation*, 2nd ed., 2002;
//! Dubbeldam, Calero, Ellis & Snurr, *Mol. Simul.* **2016**, 42, 81 (RASPA).

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

use crate::error::SciFormError;
use crate::materials::UnitCell;

const R_GAS_KCAL_MOLK: f64 = 0.001_987_204_258_640_83;
const BAR_TO_PA: f64 = 1.0e5;
const ANGSTROM3_TO_M3: f64 = 1.0e-30;

/// Lennard-Jones parameters for one host or guest atom.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LennardJonesAtom {
    /// Position in Cartesian coordinates, Å.
    pub position: [f64; 3],
    /// Well depth in kcal/mol.
    pub epsilon_kcal_mol: f64,
    /// Zero-crossing distance in Å.
    pub sigma_angstrom: f64,
}

/// A rigid periodic framework.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GcmcFramework {
    /// Unit cell.
    pub cell: UnitCell,
    /// Framework atoms with their Lennard-Jones parameters.
    pub atoms: Vec<LennardJonesAtom>,
}

impl GcmcFramework {
    /// Build an empty cubic framework of edge `a` Å.
    pub fn empty_cubic(a: f64) -> Self {
        Self {
            cell: UnitCell::cubic(a),
            atoms: Vec::new(),
        }
    }
}

/// The adsorbate model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GuestModel {
    /// Well depth in kcal/mol.
    pub epsilon_kcal_mol: f64,
    /// Zero-crossing distance in Å.
    pub sigma_angstrom: f64,
}

/// Configuration for a GCMC isotherm.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GcmcConfig {
    /// Temperature in K.
    pub temperature_k: f64,
    /// Pressure points in bar.
    pub pressures_bar: Vec<f64>,
    /// Equilibration moves per pressure point.
    pub equilibration_steps: usize,
    /// Production moves per pressure point.
    pub production_steps: usize,
    /// RNG seed.
    pub seed: u64,
    /// Maximum number of guests in the cell.
    pub max_guests: usize,
    /// Fraction of moves that are translations rather than insert/delete.
    pub translation_probability: f64,
    /// Lennard-Jones cutoff in Å.
    pub cutoff_angstrom: f64,
    /// Sample the loading every `sample_interval` production moves.
    pub sample_interval: usize,
}

impl Default for GcmcConfig {
    fn default() -> Self {
        Self {
            temperature_k: 298.15,
            pressures_bar: vec![0.1, 1.0, 10.0],
            equilibration_steps: 20_000,
            production_steps: 100_000,
            seed: 42,
            max_guests: 100,
            translation_probability: 0.3,
            cutoff_angstrom: 12.0,
            sample_interval: 100,
        }
    }
}

impl GcmcConfig {
    fn validate(&self) -> Result<(), SciFormError> {
        if !self.temperature_k.is_finite() || self.temperature_k <= 0.0 {
            return Err(SciFormError::validation(
                "temperature_k must be finite and > 0",
            ));
        }
        if self.pressures_bar.is_empty() {
            return Err(SciFormError::validation(
                "pressures_bar must contain at least one point",
            ));
        }
        if self
            .pressures_bar
            .iter()
            .any(|pressure| !pressure.is_finite() || *pressure < 0.0)
        {
            return Err(SciFormError::validation(
                "pressures_bar must be finite and >= 0",
            ));
        }
        if self.equilibration_steps == 0 || self.production_steps == 0 {
            return Err(SciFormError::validation(
                "equilibration_steps and production_steps must be >= 1",
            ));
        }
        if !self.translation_probability.is_finite()
            || !(0.0..=1.0).contains(&self.translation_probability)
        {
            return Err(SciFormError::validation(
                "translation_probability must be in [0, 1]",
            ));
        }
        if !self.cutoff_angstrom.is_finite() || self.cutoff_angstrom <= 0.0 {
            return Err(SciFormError::validation(
                "cutoff_angstrom must be finite and > 0",
            ));
        }
        if self.max_guests == 0 {
            return Err(SciFormError::validation("max_guests must be >= 1"));
        }
        if self.sample_interval == 0 {
            return Err(SciFormError::validation("sample_interval must be >= 1"));
        }
        Ok(())
    }
}

/// One pressure point of an adsorption isotherm.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AdsorptionPoint {
    /// Pressure in bar.
    pub pressure_bar: f64,
    /// Mean loading in molecules per unit cell.
    pub mean_loading: f64,
    /// Standard error of the mean loading.
    pub standard_error: f64,
    /// Number of production samples.
    pub n_samples: usize,
}

/// A complete adsorption isotherm.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AdsorptionIsotherm {
    /// Temperature in K.
    pub temperature_k: f64,
    /// Adsorbate model.
    pub guest: GuestModel,
    /// Isotherm points.
    pub points: Vec<AdsorptionPoint>,
    /// Notes and scope caveats.
    pub notes: Vec<String>,
}

/// Compute a GCMC adsorption isotherm on a rigid Lennard-Jones framework.
///
/// The insertion acceptance is
///
/// ```text
/// P_acc = min(1, exp(-βΔU) · P V / (kBT (N+1)))
/// ```
///
/// and deletion uses the reverse factor. The framework is periodic; the
/// guest-guest and host-guest interactions use the minimum image and a
/// Lennard-Jones cutoff. The framework is rigid and neutral in this first
/// version.
pub fn compute_gcmc_isotherm(
    framework: &GcmcFramework,
    guest: &GuestModel,
    config: &GcmcConfig,
) -> Result<AdsorptionIsotherm, SciFormError> {
    config.validate()?;
    if !guest.epsilon_kcal_mol.is_finite() || guest.epsilon_kcal_mol < 0.0 {
        return Err(SciFormError::validation(
            "guest epsilon_kcal_mol must be finite and >= 0",
        ));
    }
    if !guest.sigma_angstrom.is_finite() || guest.sigma_angstrom <= 0.0 {
        return Err(SciFormError::validation(
            "guest sigma_angstrom must be finite and > 0",
        ));
    }
    if framework.atoms.iter().any(|atom| {
        !atom.position.iter().all(|value| value.is_finite())
            || !atom.epsilon_kcal_mol.is_finite()
            || atom.epsilon_kcal_mol < 0.0
            || !atom.sigma_angstrom.is_finite()
            || atom.sigma_angstrom <= 0.0
    }) {
        return Err(SciFormError::validation(
            "framework atoms must have finite positions and valid Lennard-Jones parameters",
        ));
    }
    if framework.cell.volume() <= 0.0 {
        return Err(SciFormError::validation(
            "the framework cell must have positive volume",
        ));
    }

    let mut points = Vec::with_capacity(config.pressures_bar.len());
    for &pressure in &config.pressures_bar {
        let mut rng = StdRng::seed_from_u64(config.seed);
        let mut guests: Vec<[f64; 3]> = Vec::new();
        let volume_m3 = framework.cell.volume() * ANGSTROM3_TO_M3;
        let beta = 1.0 / (R_GAS_KCAL_MOLK * config.temperature_k);

        let mut attempt = |step: usize, guests: &mut Vec<[f64; 3]>| -> Result<(), SciFormError> {
            if pressure <= 0.0 {
                guests.clear();
                return Ok(());
            }
            // P V / (kBT) with P in Pa, V in m³ and kBT in joules.
            let pv_over_kbt =
                pressure * BAR_TO_PA * volume_m3 / (1.380_649e-23 * config.temperature_k);

            let translate = rng.gen::<f64>() < config.translation_probability;
            if translate && !guests.is_empty() {
                let index = rng.gen_range(0..guests.len());
                let old = guests[index];
                let displacement = [
                    (rng.gen::<f64>() - 0.5) * 0.5,
                    (rng.gen::<f64>() - 0.5) * 0.5,
                    (rng.gen::<f64>() - 0.5) * 0.5,
                ];
                let candidate = framework.cell.wrap_cart([
                    old[0] + displacement[0],
                    old[1] + displacement[1],
                    old[2] + displacement[2],
                ]);
                let old_energy = guest_energy(framework, guest, guests, index, config);
                guests[index] = candidate;
                let new_energy = guest_energy(framework, guest, guests, index, config);
                let delta = new_energy - old_energy;
                if delta <= 0.0 || rng.gen::<f64>() < (-beta * delta).exp() {
                    // accept
                } else {
                    guests[index] = old;
                }
                return Ok(());
            }

            // Decide insertion vs deletion. At N = 0 always insert.
            let delete =
                !guests.is_empty() && (guests.len() >= config.max_guests || rng.gen::<f64>() < 0.5);
            if delete {
                let index = rng.gen_range(0..guests.len());
                let old = guests[index];
                let old_energy = guest_energy(framework, guest, guests, index, config);
                let n_before = guests.len() as f64;
                guests.swap_remove(index);
                let factor = n_before / pv_over_kbt.max(1e-30);
                let delta = -old_energy;
                if delta <= 0.0 || rng.gen::<f64>() < (factor * (-beta * delta).exp()).min(1.0) {
                    // accept
                } else {
                    guests.push(old);
                }
            } else if guests.len() < config.max_guests {
                let candidate = [rng.gen::<f64>(), rng.gen::<f64>(), rng.gen::<f64>()];
                let position = framework.cell.frac_to_cart(candidate);
                guests.push(position);
                let index = guests.len() - 1;
                let energy = guest_energy(framework, guest, guests, index, config);
                let factor = pv_over_kbt / (guests.len() as f64);
                if rng.gen::<f64>() < (factor * (-beta * energy).exp()).min(1.0) {
                    // accept
                } else {
                    guests.pop();
                }
            }
            let _ = step;
            Ok(())
        };

        for step in 0..config.equilibration_steps {
            attempt(step, &mut guests)?;
        }
        let mut samples = Vec::new();
        for step in 0..config.production_steps {
            attempt(step, &mut guests)?;
            if step % config.sample_interval == 0 {
                samples.push(guests.len() as f64);
            }
        }
        let mean = if samples.is_empty() {
            guests.len() as f64
        } else {
            samples.iter().sum::<f64>() / samples.len() as f64
        };
        let standard_error = if samples.len() > 1 {
            let variance = samples
                .iter()
                .map(|value| (value - mean).powi(2))
                .sum::<f64>()
                / (samples.len() - 1) as f64;
            (variance / samples.len() as f64).sqrt()
        } else {
            0.0
        };
        points.push(AdsorptionPoint {
            pressure_bar: pressure,
            mean_loading: mean,
            standard_error,
            n_samples: samples.len(),
        });
    }

    Ok(AdsorptionIsotherm {
        temperature_k: config.temperature_k,
        guest: guest.clone(),
        points,
        notes: vec![
            "Rigid, neutral Lennard-Jones framework; no framework charges or Ewald sums yet."
                .to_string(),
            "Guest-guest and host-guest interactions use the minimum image and a spherical cutoff."
                .to_string(),
            "Loading is molecules per unit cell; convert with the cell volume for per-volume units."
                .to_string(),
        ],
    })
}

fn guest_energy(
    framework: &GcmcFramework,
    guest: &GuestModel,
    guests: &[[f64; 3]],
    index: usize,
    config: &GcmcConfig,
) -> f64 {
    if index >= guests.len() {
        return 0.0;
    }
    let position = guests[index];
    let mut energy = 0.0;
    for atom in &framework.atoms {
        let r = framework
            .cell
            .minimum_image_distance(position, atom.position);
        if r < config.cutoff_angstrom {
            let sigma = 0.5 * (atom.sigma_angstrom + guest.sigma_angstrom);
            let epsilon = (atom.epsilon_kcal_mol * guest.epsilon_kcal_mol).sqrt();
            energy += lennard_jones(r, epsilon, sigma);
        }
    }
    for (other_index, other) in guests.iter().enumerate() {
        if other_index == index {
            continue;
        }
        let r = framework.cell.minimum_image_distance(position, *other);
        if r < config.cutoff_angstrom {
            energy += lennard_jones(r, guest.epsilon_kcal_mol, guest.sigma_angstrom);
        }
    }
    // Count guest-guest pairs twice; return the half.
    if guests.len() > 1 {
        // The caller uses this as a per-particle energy, so the pair sum is
        // already the interaction energy of particle `index` with all others.
        // No factor is applied here.
    }
    energy
}

fn lennard_jones(r: f64, epsilon: f64, sigma: f64) -> f64 {
    if epsilon == 0.0 {
        return 0.0;
    }
    let safe_r = r.max(0.25 * sigma);
    let ratio6 = (sigma / safe_r).powi(6);
    4.0 * epsilon * (ratio6 * ratio6 - ratio6)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lennard_jones_has_minimum_at_sigma_times_two_to_the_sixth() {
        let epsilon = 1.0;
        let sigma = 3.0;
        let r_min = sigma * 2.0_f64.powf(1.0 / 6.0);
        let value = lennard_jones(r_min, epsilon, sigma);
        assert!((value + epsilon).abs() < 1e-9, "value = {value}");
    }

    #[test]
    fn ideal_gas_loading_follows_pv_over_kbt() {
        let framework = GcmcFramework::empty_cubic(30.0);
        let guest = GuestModel {
            epsilon_kcal_mol: 0.0,
            sigma_angstrom: 3.0,
        };
        let config = GcmcConfig {
            temperature_k: 300.0,
            pressures_bar: vec![10.0],
            equilibration_steps: 20_000,
            production_steps: 80_000,
            seed: 7,
            max_guests: 100,
            translation_probability: 0.0,
            cutoff_angstrom: 12.0,
            sample_interval: 100,
        };
        let result = compute_gcmc_isotherm(&framework, &guest, &config).unwrap();
        let ideal =
            10.0 * BAR_TO_PA * framework.cell.volume() * ANGSTROM3_TO_M3 / (1.380_649e-23 * 300.0);
        let loading = result.points[0].mean_loading;
        assert!(
            (loading - ideal).abs() < 0.35 * ideal.max(1.0),
            "loading = {loading}, ideal = {ideal}"
        );
    }

    #[test]
    fn attractive_framework_increases_loading() {
        let guest = GuestModel {
            epsilon_kcal_mol: 0.2,
            sigma_angstrom: 3.0,
        };
        let empty = GcmcFramework::empty_cubic(20.0);
        let host = GcmcFramework {
            cell: UnitCell::cubic(20.0),
            atoms: vec![LennardJonesAtom {
                position: [10.0, 10.0, 10.0],
                epsilon_kcal_mol: 1.0,
                sigma_angstrom: 3.0,
            }],
        };
        let config = GcmcConfig {
            pressures_bar: vec![1.0],
            ..GcmcConfig::default()
        };
        let empty_result = compute_gcmc_isotherm(&empty, &guest, &config).unwrap();
        let host_result = compute_gcmc_isotherm(&host, &guest, &config).unwrap();
        assert!(host_result.points[0].mean_loading >= empty_result.points[0].mean_loading);
    }
}
