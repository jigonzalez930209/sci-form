//! Grand-canonical Monte Carlo adsorption on a periodic framework.
//!
//! The engine supports rigid frameworks with Lennard-Jones host atoms and
//! optional point charges, Ewald electrostatics, and single- or multi-site
//! guest molecules. The first implementation remains deliberately smaller than
//! a production MOF code: there is no framework flexibility and no bonded MM
//! term inside the framework.
//!
//! # Units
//!
//! * positions and cell lengths in **Å**;
//! * Lennard-Jones `epsilon` in **kcal/mol** and `sigma` in **Å**;
//! * charges in **elementary charge**;
//! * pressure in **bar**;
//! * temperature in **K**;
//! * loading in **molecules per unit cell**.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

use crate::error::SciFormError;
use crate::materials::UnitCell;

const R_GAS_KCAL_MOLK: f64 = 0.001_987_204_258_640_83;
const BAR_TO_PA: f64 = 1.0e5;
const ANGSTROM3_TO_M3: f64 = 1.0e-30;
const COULOMB_KCAL_ANGSTROM: f64 = 332.063_71;

/// Lennard-Jones parameters and optional charge for one host atom.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LennardJonesAtom {
    /// Position in Cartesian coordinates, Å.
    pub position: [f64; 3],
    /// Well depth in kcal/mol.
    pub epsilon_kcal_mol: f64,
    /// Zero-crossing distance in Å.
    pub sigma_angstrom: f64,
    /// Point charge in elementary charge.
    #[serde(default)]
    pub charge: f64,
}

/// One interaction site of a guest molecule.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GuestSite {
    /// Site offset in the molecular frame, Å.
    pub offset_angstrom: [f64; 3],
    /// Well depth in kcal/mol.
    pub epsilon_kcal_mol: f64,
    /// Zero-crossing distance in Å.
    pub sigma_angstrom: f64,
    /// Point charge in elementary charge.
    #[serde(default)]
    pub charge: f64,
}

/// A rigid periodic framework.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GcmcFramework {
    /// Unit cell.
    pub cell: UnitCell,
    /// Framework atoms with Lennard-Jones parameters and charges.
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

    fn total_charge(&self) -> f64 {
        self.atoms.iter().map(|atom| atom.charge).sum()
    }
}

/// The adsorbate model.
///
/// A single-site guest uses `epsilon_kcal_mol` and `sigma_angstrom` with an
/// empty `sites` vector. A multi-site guest sets `sites`; the legacy fields are
/// then ignored.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GuestModel {
    /// Multi-site definition. Empty means a single site at the origin.
    #[serde(default)]
    pub sites: Vec<GuestSite>,
    /// Single-site well depth in kcal/mol.
    #[serde(default)]
    pub epsilon_kcal_mol: f64,
    /// Single-site zero-crossing distance in Å.
    #[serde(default)]
    pub sigma_angstrom: f64,
    /// Single-site point charge.
    #[serde(default)]
    pub charge: f64,
}

impl GuestModel {
    /// Single-site constructor.
    pub fn single_site(epsilon_kcal_mol: f64, sigma_angstrom: f64, charge: f64) -> Self {
        Self {
            sites: Vec::new(),
            epsilon_kcal_mol,
            sigma_angstrom,
            charge,
        }
    }

    /// All interaction sites, including the legacy single-site fallback.
    pub fn interaction_sites(&self) -> Vec<GuestSite> {
        if self.sites.is_empty() {
            vec![GuestSite {
                offset_angstrom: [0.0, 0.0, 0.0],
                epsilon_kcal_mol: self.epsilon_kcal_mol,
                sigma_angstrom: self.sigma_angstrom,
                charge: self.charge,
            }]
        } else {
            self.sites.clone()
        }
    }

    /// Total molecular charge.
    pub fn total_charge(&self) -> f64 {
        self.interaction_sites()
            .iter()
            .map(|site| site.charge)
            .sum()
    }
}

/// One inserted guest molecule with its centre and orientation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GuestMolecule {
    /// Molecular centre in Å.
    pub center: [f64; 3],
    /// Row-major rotation matrix mapping molecular offsets to the cell frame.
    pub rotation: [[f64; 3]; 3],
}

impl GuestMolecule {
    /// A guest at `center` with the identity orientation.
    pub fn new(center: [f64; 3]) -> Self {
        Self {
            center,
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }

    fn site_positions(&self, sites: &[GuestSite]) -> Vec<[f64; 3]> {
        sites
            .iter()
            .map(|site| {
                let offset = site.offset_angstrom;
                [
                    self.center[0]
                        + self.rotation[0][0] * offset[0]
                        + self.rotation[0][1] * offset[1]
                        + self.rotation[0][2] * offset[2],
                    self.center[1]
                        + self.rotation[1][0] * offset[0]
                        + self.rotation[1][1] * offset[1]
                        + self.rotation[1][2] * offset[2],
                    self.center[2]
                        + self.rotation[2][0] * offset[0]
                        + self.rotation[2][1] * offset[1]
                        + self.rotation[2][2] * offset[2],
                ]
            })
            .collect()
    }
}

/// Ewald configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EwaldConfig {
    /// Ewald splitting parameter in Å⁻¹.
    pub alpha_angstrom_inv: f64,
    /// Maximum reciprocal lattice index in each direction.
    pub k_max: usize,
    /// Real-space cutoff in Å.
    pub real_cutoff_angstrom: f64,
}

impl Default for EwaldConfig {
    fn default() -> Self {
        Self {
            alpha_angstrom_inv: 0.3,
            k_max: 6,
            real_cutoff_angstrom: 8.0,
        }
    }
}

impl EwaldConfig {
    fn validate(&self) -> Result<(), SciFormError> {
        if !self.alpha_angstrom_inv.is_finite() || self.alpha_angstrom_inv <= 0.0 {
            return Err(SciFormError::validation(
                "Ewald alpha_angstrom_inv must be finite and > 0",
            ));
        }
        if self.k_max == 0 {
            return Err(SciFormError::validation("Ewald k_max must be >= 1"));
        }
        if !self.real_cutoff_angstrom.is_finite() || self.real_cutoff_angstrom <= 0.0 {
            return Err(SciFormError::validation(
                "Ewald real_cutoff_angstrom must be finite and > 0",
            ));
        }
        Ok(())
    }
}

/// Electrostatic model for the GCMC run.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum GcmcElectrostatics {
    /// No electrostatic interaction.
    None,
    /// Periodic Ewald summation with a neutralizing background.
    Ewald(EwaldConfig),
}

/// Configuration for a GCMC isotherm.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
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
    /// Fraction of translation moves that also rotate the molecule.
    #[serde(default)]
    pub rotation_probability: f64,
    /// Lennard-Jones cutoff in Å.
    pub cutoff_angstrom: f64,
    /// Sample the loading every `sample_interval` production moves.
    pub sample_interval: usize,
    /// Electrostatic model.
    #[serde(default = "default_electrostatics")]
    pub electrostatics: GcmcElectrostatics,
}

fn default_electrostatics() -> GcmcElectrostatics {
    GcmcElectrostatics::None
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
            rotation_probability: 0.0,
            cutoff_angstrom: 12.0,
            sample_interval: 100,
            electrostatics: GcmcElectrostatics::None,
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
        if !self.rotation_probability.is_finite()
            || !(0.0..=1.0).contains(&self.rotation_probability)
        {
            return Err(SciFormError::validation(
                "rotation_probability must be in [0, 1]",
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
        if let GcmcElectrostatics::Ewald(config) = &self.electrostatics {
            config.validate()?;
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

/// Compute a GCMC adsorption isotherm on a rigid periodic framework.
pub fn compute_gcmc_isotherm(
    framework: &GcmcFramework,
    guest: &GuestModel,
    config: &GcmcConfig,
) -> Result<AdsorptionIsotherm, SciFormError> {
    config.validate()?;
    validate_framework(framework)?;
    let sites = guest.interaction_sites();
    if sites.is_empty() {
        return Err(SciFormError::validation(
            "the guest model must contain at least one interaction site",
        ));
    }
    for site in &sites {
        if !site.epsilon_kcal_mol.is_finite()
            || site.epsilon_kcal_mol < 0.0
            || !site.sigma_angstrom.is_finite()
            || site.sigma_angstrom <= 0.0
            || !site.charge.is_finite()
            || !site.offset_angstrom.iter().all(|value| value.is_finite())
        {
            return Err(SciFormError::validation(
                "guest sites must have finite positions, epsilon >= 0, sigma > 0 and finite charge",
            ));
        }
    }
    if framework.cell.volume() <= 0.0 {
        return Err(SciFormError::validation(
            "the framework cell must have positive volume",
        ));
    }

    let framework_charge = framework.total_charge();
    let mut points = Vec::with_capacity(config.pressures_bar.len());
    for &pressure in &config.pressures_bar {
        let mut rng = StdRng::seed_from_u64(config.seed);
        let mut guests: Vec<GuestMolecule> = Vec::new();
        let volume_m3 = framework.cell.volume() * ANGSTROM3_TO_M3;
        let beta = 1.0 / (R_GAS_KCAL_MOLK * config.temperature_k);
        let pv_over_kbt = pressure * BAR_TO_PA * volume_m3 / (1.380_649e-23 * config.temperature_k);

        let attempt = |rng: &mut StdRng,
                       guests: &mut Vec<GuestMolecule>|
         -> Result<(), SciFormError> {
            if pressure <= 0.0 {
                guests.clear();
                return Ok(());
            }
            let old_energy = total_energy(framework, &sites, guests, config);
            let translate = rng.gen::<f64>() < config.translation_probability;
            if translate && !guests.is_empty() {
                let index = rng.gen_range(0..guests.len());
                let old = guests[index].clone();
                let displacement = [
                    (rng.gen::<f64>() - 0.5) * 0.5,
                    (rng.gen::<f64>() - 0.5) * 0.5,
                    (rng.gen::<f64>() - 0.5) * 0.5,
                ];
                guests[index].center = framework.cell.wrap_cart([
                    old.center[0] + displacement[0],
                    old.center[1] + displacement[1],
                    old.center[2] + displacement[2],
                ]);
                if rng.gen::<f64>() < config.rotation_probability {
                    guests[index].rotation = random_rotation(rng);
                }
                let new_energy = total_energy(framework, &sites, guests, config);
                let delta = new_energy - old_energy;
                if delta <= 0.0 || rng.gen::<f64>() < (-beta * delta).exp() {
                    return Ok(());
                }
                guests[index] = old;
                return Ok(());
            }

            let delete =
                !guests.is_empty() && (guests.len() >= config.max_guests || rng.gen::<f64>() < 0.5);
            if delete {
                let index = rng.gen_range(0..guests.len());
                let old = guests.remove(index);
                let new_energy = total_energy(framework, &sites, guests, config);
                let delta = new_energy - old_energy;
                let factor = (guests.len() + 1) as f64 / pv_over_kbt.max(1e-30);
                if delta <= 0.0 || rng.gen::<f64>() < (factor * (-beta * delta).exp()).min(1.0) {
                    return Ok(());
                }
                guests.insert(index, old);
            } else if guests.len() < config.max_guests {
                let center = [rng.gen::<f64>(), rng.gen::<f64>(), rng.gen::<f64>()];
                let center = framework.cell.frac_to_cart(center);
                let molecule = GuestMolecule {
                    center,
                    rotation: random_rotation(rng),
                };
                guests.push(molecule);
                let new_energy = total_energy(framework, &sites, guests, config);
                let delta = new_energy - old_energy;
                let factor = pv_over_kbt / guests.len() as f64;
                if rng.gen::<f64>() < (factor * (-beta * delta).exp()).min(1.0) {
                    return Ok(());
                }
                guests.pop();
            }
            Ok(())
        };

        for _ in 0..config.equilibration_steps {
            attempt(&mut rng, &mut guests)?;
        }
        let mut samples = Vec::new();
        for step in 0..config.production_steps {
            attempt(&mut rng, &mut guests)?;
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
            "Rigid framework; no framework flexibility or bonded MM terms.".to_string(),
            "Host/guest and guest/guest interactions use Lorentz-Berthelot mixing.".to_string(),
            "Electrostatics is optional; Ewald uses a neutralizing background.".to_string(),
            format!("Framework net charge: {framework_charge:.3} e."),
            "Loading is molecules per unit cell.".to_string(),
        ],
    })
}

fn validate_framework(framework: &GcmcFramework) -> Result<(), SciFormError> {
    if framework.atoms.iter().any(|atom| {
        !atom.position.iter().all(|value| value.is_finite())
            || !atom.epsilon_kcal_mol.is_finite()
            || atom.epsilon_kcal_mol < 0.0
            || !atom.sigma_angstrom.is_finite()
            || atom.sigma_angstrom <= 0.0
            || !atom.charge.is_finite()
    }) {
        return Err(SciFormError::validation(
            "framework atoms must have finite positions, epsilon >= 0, sigma > 0 and finite charge",
        ));
    }
    Ok(())
}

fn total_energy(
    framework: &GcmcFramework,
    sites: &[GuestSite],
    guests: &[GuestMolecule],
    config: &GcmcConfig,
) -> f64 {
    lennard_jones_energy(framework, sites, guests, config)
        + ewald_energy(framework, sites, guests, config)
}

fn lennard_jones_energy(
    framework: &GcmcFramework,
    sites: &[GuestSite],
    guests: &[GuestMolecule],
    config: &GcmcConfig,
) -> f64 {
    let mut energy = 0.0;
    for molecule in guests {
        let positions = molecule.site_positions(sites);
        for (site, position) in sites.iter().zip(positions) {
            for atom in &framework.atoms {
                let r = framework
                    .cell
                    .minimum_image_distance(position, atom.position);
                if r < config.cutoff_angstrom {
                    let sigma = 0.5 * (atom.sigma_angstrom + site.sigma_angstrom);
                    let epsilon = (atom.epsilon_kcal_mol * site.epsilon_kcal_mol).sqrt();
                    energy += lennard_jones(r, epsilon, sigma);
                }
            }
        }
    }
    for i in 0..guests.len() {
        let positions_i = guests[i].site_positions(sites);
        for j in (i + 1)..guests.len() {
            let positions_j = guests[j].site_positions(sites);
            for (site_i, position_i) in sites.iter().zip(&positions_i) {
                for (site_j, position_j) in sites.iter().zip(&positions_j) {
                    let r = framework
                        .cell
                        .minimum_image_distance(*position_i, *position_j);
                    if r < config.cutoff_angstrom {
                        let sigma = 0.5 * (site_i.sigma_angstrom + site_j.sigma_angstrom);
                        let epsilon = (site_i.epsilon_kcal_mol * site_j.epsilon_kcal_mol).sqrt();
                        energy += lennard_jones(r, epsilon, sigma);
                    }
                }
            }
        }
    }
    energy
}

fn ewald_energy(
    framework: &GcmcFramework,
    sites: &[GuestSite],
    guests: &[GuestMolecule],
    config: &GcmcConfig,
) -> f64 {
    let GcmcElectrostatics::Ewald(ewald) = &config.electrostatics else {
        return 0.0;
    };
    let mut charges = Vec::new();
    let mut positions = Vec::new();
    for atom in &framework.atoms {
        if atom.charge != 0.0 {
            charges.push(atom.charge);
            positions.push(atom.position);
        }
    }
    for molecule in guests {
        let site_positions = molecule.site_positions(sites);
        for (site, position) in sites.iter().zip(site_positions) {
            if site.charge != 0.0 {
                charges.push(site.charge);
                positions.push(position);
            }
        }
    }
    if charges.is_empty() {
        return 0.0;
    }
    let volume = framework.cell.volume();
    let alpha = ewald.alpha_angstrom_inv;
    let n = charges.len();
    let mut energy = 0.0;

    // Real space.
    for i in 0..n {
        for j in (i + 1)..n {
            let r = framework
                .cell
                .minimum_image_distance(positions[i], positions[j]);
            if r > 1e-12 && r < ewald.real_cutoff_angstrom {
                energy += charges[i] * charges[j] * libm::erfc(alpha * r) / r;
            }
        }
    }
    // Reciprocal space.
    let reciprocal = reciprocal_lattice(&framework.cell);
    let two_pi = 2.0 * std::f64::consts::PI;
    for nx in -(ewald.k_max as i32)..=ewald.k_max as i32 {
        for ny in -(ewald.k_max as i32)..=ewald.k_max as i32 {
            for nz in -(ewald.k_max as i32)..=ewald.k_max as i32 {
                if nx == 0 && ny == 0 && nz == 0 {
                    continue;
                }
                let k = [
                    two_pi
                        * (nx as f64 * reciprocal[0][0]
                            + ny as f64 * reciprocal[1][0]
                            + nz as f64 * reciprocal[2][0]),
                    two_pi
                        * (nx as f64 * reciprocal[0][1]
                            + ny as f64 * reciprocal[1][1]
                            + nz as f64 * reciprocal[2][1]),
                    two_pi
                        * (nx as f64 * reciprocal[0][2]
                            + ny as f64 * reciprocal[1][2]
                            + nz as f64 * reciprocal[2][2]),
                ];
                let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
                if k2 < 1e-12 {
                    continue;
                }
                let mut real = 0.0;
                let mut imag = 0.0;
                for (charge, position) in charges.iter().zip(&positions) {
                    let phase = k[0] * position[0] + k[1] * position[1] + k[2] * position[2];
                    real += charge * phase.cos();
                    imag += charge * phase.sin();
                }
                let structure = real * real + imag * imag;
                energy += (2.0 * std::f64::consts::PI / volume)
                    * (-k2 / (4.0 * alpha * alpha)).exp()
                    / k2
                    * structure;
            }
        }
    }
    // Self energy and neutralizing background.
    let charge_sum: f64 = charges.iter().sum();
    let self_energy = -alpha / std::f64::consts::PI.sqrt()
        * charges.iter().map(|charge| charge * charge).sum::<f64>();
    let background =
        -std::f64::consts::PI / (2.0 * volume * alpha * alpha) * charge_sum * charge_sum;
    COULOMB_KCAL_ANGSTROM * (energy + self_energy + background)
}

fn reciprocal_lattice(cell: &UnitCell) -> [[f64; 3]; 3] {
    let a = cell.lattice[0];
    let b = cell.lattice[1];
    let c = cell.lattice[2];
    let volume = cell.volume();
    let bxc = cross(b, c);
    let cxa = cross(c, a);
    let axb = cross(a, b);
    [
        [bxc[0] / volume, bxc[1] / volume, bxc[2] / volume],
        [cxa[0] / volume, cxa[1] / volume, cxa[2] / volume],
        [axb[0] / volume, axb[1] / volume, axb[2] / volume],
    ]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn random_rotation(rng: &mut StdRng) -> [[f64; 3]; 3] {
    let a = rng.gen::<f64>() * std::f64::consts::TAU;
    let b = rng.gen::<f64>() * std::f64::consts::TAU;
    let c = rng.gen::<f64>() * std::f64::consts::TAU;
    let (sa, ca) = a.sin_cos();
    let (sb, cb) = b.sin_cos();
    let (sc, cc) = c.sin_cos();
    [
        [cb * cc, -cb * sc, sb],
        [sa * sb * cc + ca * sc, -sa * sb * sc + ca * cc, -sa * cb],
        [-ca * sb * cc + sa * sc, ca * sb * sc + sa * cc, ca * cb],
    ]
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
        let guest = GuestModel::single_site(0.0, 3.0, 0.0);
        let config = GcmcConfig {
            temperature_k: 300.0,
            pressures_bar: vec![10.0],
            equilibration_steps: 20_000,
            production_steps: 80_000,
            seed: 7,
            translation_probability: 0.0,
            ..GcmcConfig::default()
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
        let guest = GuestModel::single_site(0.2, 3.0, 0.0);
        let empty = GcmcFramework::empty_cubic(20.0);
        let host = GcmcFramework {
            cell: UnitCell::cubic(20.0),
            atoms: vec![LennardJonesAtom {
                position: [10.0, 10.0, 10.0],
                epsilon_kcal_mol: 1.0,
                sigma_angstrom: 3.0,
                charge: 0.0,
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

    #[test]
    fn multi_site_guest_has_a_total_charge() {
        let guest = GuestModel {
            sites: vec![
                GuestSite {
                    offset_angstrom: [-0.5, 0.0, 0.0],
                    epsilon_kcal_mol: 0.1,
                    sigma_angstrom: 3.0,
                    charge: -0.5,
                },
                GuestSite {
                    offset_angstrom: [0.5, 0.0, 0.0],
                    epsilon_kcal_mol: 0.1,
                    sigma_angstrom: 3.0,
                    charge: 0.5,
                },
            ],
            epsilon_kcal_mol: 0.0,
            sigma_angstrom: 1.0,
            charge: 0.0,
        };
        assert!((guest.total_charge()).abs() < 1e-12);
        assert_eq!(guest.interaction_sites().len(), 2);
    }

    #[test]
    fn ewald_opposite_charges_match_direct_coulomb() {
        let cell = UnitCell::cubic(20.0);
        let framework = GcmcFramework {
            cell: cell.clone(),
            atoms: vec![
                LennardJonesAtom {
                    position: [5.0, 5.0, 5.0],
                    epsilon_kcal_mol: 0.0,
                    sigma_angstrom: 1.0,
                    charge: 1.0,
                },
                LennardJonesAtom {
                    position: [6.0, 5.0, 5.0],
                    epsilon_kcal_mol: 0.0,
                    sigma_angstrom: 1.0,
                    charge: -1.0,
                },
            ],
        };
        let guest = GuestModel::single_site(0.0, 1.0, 0.0);
        let config = GcmcConfig {
            electrostatics: GcmcElectrostatics::Ewald(EwaldConfig {
                alpha_angstrom_inv: 0.3,
                k_max: 6,
                real_cutoff_angstrom: 10.0,
            }),
            ..GcmcConfig::default()
        };
        let energy = ewald_energy(&framework, &guest.interaction_sites(), &[], &config);
        let direct = -COULOMB_KCAL_ANGSTROM;
        assert!(
            (energy - direct).abs() < 0.01 * direct.abs(),
            "Ewald={energy}, direct={direct}"
        );
    }
}
