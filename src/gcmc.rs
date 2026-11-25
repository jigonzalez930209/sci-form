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
    /// Critical temperature in K, used for the Peng-Robinson fugacity.
    #[serde(default)]
    pub critical_temperature_k: Option<f64>,
    /// Critical pressure in bar, used for the Peng-Robinson fugacity.
    #[serde(default)]
    pub critical_pressure_bar: Option<f64>,
    /// Acentric factor, used for the Peng-Robinson fugacity.
    #[serde(default)]
    pub acentric_factor: Option<f64>,
}

impl GuestModel {
    /// Single-site constructor.
    pub fn single_site(epsilon_kcal_mol: f64, sigma_angstrom: f64, charge: f64) -> Self {
        Self {
            sites: Vec::new(),
            epsilon_kcal_mol,
            sigma_angstrom,
            charge,
            critical_temperature_k: None,
            critical_pressure_bar: None,
            acentric_factor: None,
        }
    }

    /// Attach the critical constants used by the Peng-Robinson fugacity.
    pub fn with_critical_constants(
        mut self,
        critical_temperature_k: f64,
        critical_pressure_bar: f64,
        acentric_factor: f64,
    ) -> Self {
        self.critical_temperature_k = Some(critical_temperature_k);
        self.critical_pressure_bar = Some(critical_pressure_bar);
        self.acentric_factor = Some(acentric_factor);
        self
    }

    /// Peng-Robinson fugacity coefficient at `temperature_k` and `pressure_bar`.
    ///
    /// Returns `1.0` when the critical constants are absent, which keeps the
    /// ideal-gas limit used by the earlier prototype tests.
    pub fn fugacity_coefficient(
        &self,
        temperature_k: f64,
        pressure_bar: f64,
    ) -> Result<f64, SciFormError> {
        let (Some(tc), Some(pc), Some(omega)) = (
            self.critical_temperature_k,
            self.critical_pressure_bar,
            self.acentric_factor,
        ) else {
            return Ok(1.0);
        };
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(SciFormError::validation(
                "temperature_k must be finite and > 0",
            ));
        }
        if !pressure_bar.is_finite() || pressure_bar < 0.0 {
            return Err(SciFormError::validation(
                "pressure_bar must be finite and >= 0",
            ));
        }
        if pressure_bar == 0.0 {
            return Ok(1.0);
        }
        if !tc.is_finite() || tc <= 0.0 || !pc.is_finite() || pc <= 0.0 || !omega.is_finite() {
            return Err(SciFormError::validation(
                "critical constants must be finite and positive",
            ));
        }

        // Peng-Robinson in SI units.
        const R: f64 = 8.314_462_618;
        let pressure_pa = pressure_bar * BAR_TO_PA;
        let a = 0.457_24 * R * R * tc * tc / (pc * BAR_TO_PA);
        let b = 0.077_80 * R * tc / (pc * BAR_TO_PA);
        let kappa = 0.374_64 + 1.542_26 * omega - 0.269_92 * omega * omega;
        let sqrt_tr = (temperature_k / tc).sqrt();
        let alpha = (1.0 + kappa * (1.0 - sqrt_tr)).powi(2);
        let a_alpha = a * alpha;
        let big_a = a_alpha * pressure_pa / (R * R * temperature_k * temperature_k);
        let big_b = b * pressure_pa / (R * temperature_k);

        // Largest real root of Z^3 - (1-B)Z^2 + (A-3B^2-2B)Z - (AB-B^2-B^3).
        let z = largest_cubic_root(
            -(1.0 - big_b),
            big_a - 3.0 * big_b * big_b - 2.0 * big_b,
            -(big_a * big_b - big_b * big_b - big_b * big_b * big_b),
        )
        .ok_or_else(|| SciFormError::numerical("Peng-Robinson cubic has no real gas root"))?;
        if z <= big_b {
            return Err(SciFormError::numerical(
                "Peng-Robinson root is not in the gas branch",
            ));
        }
        let sqrt_two = std::f64::consts::SQRT_2;
        let log_phi = z
            - 1.0
            - (z - big_b).ln()
            - big_a / (2.0 * sqrt_two * big_b)
                * ((z + (1.0 + sqrt_two) * big_b) / (z + (1.0 - sqrt_two) * big_b)).ln();
        let phi = log_phi.exp();
        if !phi.is_finite() || phi <= 0.0 {
            return Err(SciFormError::numerical(
                "Peng-Robinson fugacity coefficient is not finite",
            ));
        }
        Ok(phi)
    }

    /// Fugacity in bar at `temperature_k` and `pressure_bar`.
    pub fn fugacity_bar(&self, temperature_k: f64, pressure_bar: f64) -> Result<f64, SciFormError> {
        Ok(pressure_bar * self.fugacity_coefficient(temperature_k, pressure_bar)?)
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
    /// Real-space Coulomb sum with the GCMC cutoff and minimum image.
    Cutoff,
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
    /// Apply the analytic Lennard-Jones long-range tail correction.
    #[serde(default = "default_tail_correction")]
    pub tail_correction: bool,
    /// Apply RASPA-style switching to real-space LJ and Coulomb terms.
    #[serde(default)]
    pub switching: bool,
    /// Shift the Lennard-Jones potential by its value at the cutoff
    /// (`U(r) - U(r_cut)`), the convention used by RASPA's "shifted"
    /// force fields.
    #[serde(default)]
    pub shift_lj: bool,
}

fn default_tail_correction() -> bool {
    true
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
            tail_correction: true,
            switching: false,
            shift_lj: false,
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
    /// Fugacity in bar used by the acceptance rule.
    #[serde(default)]
    pub fugacity_bar: f64,
    /// Mean loading in molecules per unit cell.
    pub mean_loading: f64,
    /// Standard error of the mean loading.
    pub standard_error: f64,
    /// Number of production samples.
    pub n_samples: usize,
    /// Mean guest-host plus guest-guest energy in kcal/mol per unit cell.
    #[serde(default)]
    pub mean_energy_kcal_mol: f64,
    /// Mean energy per adsorbed molecule in kcal/mol.
    #[serde(default)]
    pub mean_energy_per_molecule_kcal_mol: f64,
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

/// Uniform grid of framework atoms for the hot host-guest loop.
struct FrameworkGrid {
    lengths: [f64; 3],
    shape: [usize; 3],
    bins: Vec<Vec<usize>>,
}

impl FrameworkGrid {
    fn new(framework: &GcmcFramework, cell_size: f64) -> Option<Self> {
        let lattice = framework.cell.lattice;
        for i in 0..3 {
            for j in 0..3 {
                if i != j && lattice[i][j].abs() > 1e-12 {
                    return None;
                }
            }
        }
        let lengths = [lattice[0][0], lattice[1][1], lattice[2][2]];
        if lengths.iter().any(|length| *length <= 0.0) {
            return None;
        }
        let shape = [
            (lengths[0] / cell_size).ceil().max(1.0) as usize,
            (lengths[1] / cell_size).ceil().max(1.0) as usize,
            (lengths[2] / cell_size).ceil().max(1.0) as usize,
        ];
        let mut bins = vec![Vec::new(); shape[0] * shape[1] * shape[2]];
        for (index, atom) in framework.atoms.iter().enumerate() {
            let bin = Self::bin(&lengths, &shape, atom.position);
            bins[bin].push(index);
        }
        Some(Self {
            lengths,
            shape,
            bins,
        })
    }

    fn bin(lengths: &[f64; 3], shape: &[usize; 3], position: [f64; 3]) -> usize {
        let ix = ((position[0] / lengths[0] * shape[0] as f64).floor() as isize)
            .rem_euclid(shape[0] as isize) as usize;
        let iy = ((position[1] / lengths[1] * shape[1] as f64).floor() as isize)
            .rem_euclid(shape[1] as isize) as usize;
        let iz = ((position[2] / lengths[2] * shape[2] as f64).floor() as isize)
            .rem_euclid(shape[2] as isize) as usize;
        Self::flat(shape, ix, iy, iz)
    }

    fn flat(shape: &[usize; 3], ix: usize, iy: usize, iz: usize) -> usize {
        (ix * shape[1] + iy) * shape[2] + iz
    }

    fn neighbor_indices(base: usize, size: usize) -> ([usize; 5], usize) {
        if size == 1 {
            return ([0, 0, 0, 0, 0], 1);
        }
        let mut values = [0usize; 5];
        let mut count = 0;
        for offset in -2..=2 {
            let value = (base as isize + offset).rem_euclid(size as isize) as usize;
            if !values[..count].contains(&value) {
                values[count] = value;
                count += 1;
            }
        }
        (values, count)
    }

    fn for_each_nearby(&self, position: [f64; 3], mut visit: impl FnMut(usize)) {
        let base = [
            ((position[0] / self.lengths[0] * self.shape[0] as f64).floor() as isize)
                .rem_euclid(self.shape[0] as isize) as usize,
            ((position[1] / self.lengths[1] * self.shape[1] as f64).floor() as isize)
                .rem_euclid(self.shape[1] as isize) as usize,
            ((position[2] / self.lengths[2] * self.shape[2] as f64).floor() as isize)
                .rem_euclid(self.shape[2] as isize) as usize,
        ];
        let (xs, nx) = Self::neighbor_indices(base[0], self.shape[0]);
        let (ys, ny) = Self::neighbor_indices(base[1], self.shape[1]);
        let (zs, nz) = Self::neighbor_indices(base[2], self.shape[2]);
        for &ix in &xs[..nx] {
            for &iy in &ys[..ny] {
                for &iz in &zs[..nz] {
                    for &index in &self.bins[Self::flat(&self.shape, ix, iy, iz)] {
                        visit(index);
                    }
                }
            }
        }
    }
}

/// Precomputed host contribution to the Ewald reciprocal sum.
///
/// Host-host terms are constant during a GCMC run, so only the guest-host
/// cross terms and guest-guest terms need to be evaluated per move.
struct EwaldHost {
    alpha: f64,
    volume: f64,
    host_charges: Vec<(f64, [f64; 3])>,
    /// (k vector, 2π exp(-k²/4α²)/(V k²), Re S_host, Im S_host)
    k_vectors: Vec<([f64; 3], f64, f64, f64)>,
    host_charge: f64,
    /// Constant host-only reciprocal term, removed from reported energies.
    host_reciprocal: f64,
}

impl EwaldHost {
    fn new(framework: &GcmcFramework, config: &EwaldConfig) -> Self {
        let alpha = config.alpha_angstrom_inv;
        let volume = framework.cell.volume();
        let host_charges: Vec<(f64, [f64; 3])> = framework
            .atoms
            .iter()
            .filter(|atom| atom.charge != 0.0)
            .map(|atom| (atom.charge, atom.position))
            .collect();
        let host_charge: f64 = host_charges.iter().map(|(charge, _)| *charge).sum();
        let reciprocal = reciprocal_lattice(&framework.cell);
        let two_pi = 2.0 * std::f64::consts::PI;
        let mut k_vectors = Vec::new();
        for nx in -(config.k_max as i32)..=config.k_max as i32 {
            for ny in -(config.k_max as i32)..=config.k_max as i32 {
                for nz in -(config.k_max as i32)..=config.k_max as i32 {
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
                    let factor = (2.0 * std::f64::consts::PI / volume)
                        * (-k2 / (4.0 * alpha * alpha)).exp()
                        / k2;
                    let mut real = 0.0;
                    let mut imag = 0.0;
                    for (charge, position) in &host_charges {
                        let phase = k[0] * position[0] + k[1] * position[1] + k[2] * position[2];
                        real += charge * phase.cos();
                        imag += charge * phase.sin();
                    }
                    k_vectors.push((k, factor, real, imag));
                }
            }
        }
        let host_reciprocal = k_vectors
            .iter()
            .map(|(_, factor, real, imag)| factor * (real * real + imag * imag))
            .sum();
        Self {
            alpha,
            volume,
            host_charges,
            k_vectors,
            host_charge,
            host_reciprocal,
        }
    }
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
    let grid = FrameworkGrid::new(framework, (config.cutoff_angstrom / 2.0).max(1.0));
    let ewald_host = match &config.electrostatics {
        GcmcElectrostatics::Ewald(ewald) => Some(EwaldHost::new(framework, ewald)),
        _ => None,
    };
    let mut points = Vec::with_capacity(config.pressures_bar.len());
    for &pressure in &config.pressures_bar {
        let mut rng = StdRng::seed_from_u64(config.seed);
        let mut guests: Vec<GuestMolecule> = Vec::new();
        let volume_m3 = framework.cell.volume() * ANGSTROM3_TO_M3;
        let beta = 1.0 / (R_GAS_KCAL_MOLK * config.temperature_k);
        let fugacity_bar = guest.fugacity_bar(config.temperature_k, pressure)?;
        let fv_over_kbt =
            fugacity_bar * BAR_TO_PA * volume_m3 / (1.380_649e-23 * config.temperature_k);

        let attempt =
            |rng: &mut StdRng, guests: &mut Vec<GuestMolecule>| -> Result<(), SciFormError> {
                if pressure <= 0.0 {
                    guests.clear();
                    return Ok(());
                }
                let old_energy = total_energy(
                    framework,
                    &sites,
                    guests,
                    config,
                    grid.as_ref(),
                    ewald_host.as_ref(),
                );
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
                    let new_energy = total_energy(
                        framework,
                        &sites,
                        guests,
                        config,
                        grid.as_ref(),
                        ewald_host.as_ref(),
                    );
                    let delta = new_energy - old_energy;
                    if delta <= 0.0 || rng.gen::<f64>() < (-beta * delta).exp() {
                        return Ok(());
                    }
                    guests[index] = old;
                    return Ok(());
                }

                let delete = !guests.is_empty()
                    && (guests.len() >= config.max_guests || rng.gen::<f64>() < 0.5);
                if delete {
                    let index = rng.gen_range(0..guests.len());
                    let old = guests.remove(index);
                    let new_energy = total_energy(
                        framework,
                        &sites,
                        guests,
                        config,
                        grid.as_ref(),
                        ewald_host.as_ref(),
                    );
                    let delta = new_energy - old_energy;
                    // No `delta <= 0` shortcut: the N/(βfV) prefactor applies even
                    // when the energy does not change (ideal gas).
                    let factor = (guests.len() + 1) as f64 / fv_over_kbt.max(1e-30);
                    if rng.gen::<f64>() < (factor * (-beta * delta).exp()).min(1.0) {
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
                    let new_energy = total_energy(
                        framework,
                        &sites,
                        guests,
                        config,
                        grid.as_ref(),
                        ewald_host.as_ref(),
                    );
                    let delta = new_energy - old_energy;
                    let factor = fv_over_kbt / guests.len() as f64;
                    if rng.gen::<f64>() < (factor * (-beta * delta).exp()).min(1.0) {
                        return Ok(());
                    }
                    guests.pop();
                }
                Ok(())
            };

        let progress = std::env::var_os("SCI_FORM_GCMC_PROGRESS").is_some();
        let progress_interval = 5_000usize;
        for step in 0..config.equilibration_steps {
            attempt(&mut rng, &mut guests)?;
            if progress && step % progress_interval == 0 {
                eprintln!(
                    "[gcmc] P={pressure:.3} bar equilibration {step}/{} guests={}",
                    config.equilibration_steps,
                    guests.len()
                );
            }
        }
        let mut samples = Vec::new();
        let mut energy_samples = Vec::new();
        for step in 0..config.production_steps {
            attempt(&mut rng, &mut guests)?;
            if progress && step % progress_interval == 0 {
                let running = if samples.is_empty() {
                    guests.len() as f64
                } else {
                    samples.iter().sum::<f64>() / samples.len() as f64
                };
                eprintln!(
                    "[gcmc] P={pressure:.3} bar production {step}/{} guests={} running_mean={running:.2}",
                    config.production_steps,
                    guests.len()
                );
            }
            if step % config.sample_interval == 0 {
                samples.push(guests.len() as f64);
                energy_samples.push(total_energy(
                    framework,
                    &sites,
                    &guests,
                    config,
                    grid.as_ref(),
                    ewald_host.as_ref(),
                ));
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
        let mean_energy = if energy_samples.is_empty() {
            0.0
        } else {
            energy_samples.iter().sum::<f64>() / energy_samples.len() as f64
        };
        points.push(AdsorptionPoint {
            pressure_bar: pressure,
            fugacity_bar,
            mean_loading: mean,
            standard_error,
            n_samples: samples.len(),
            mean_energy_kcal_mol: mean_energy,
            mean_energy_per_molecule_kcal_mol: if mean > 1e-12 {
                mean_energy / mean
            } else {
                0.0
            },
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
            if guest.critical_temperature_k.is_some() {
                "Insertion/deletion uses the Peng-Robinson fugacity.".to_string()
            } else {
                "Insertion/deletion uses the ideal-gas pressure (no critical constants)."
                    .to_string()
            },
            if config.tail_correction {
                "Lennard-Jones long-range tail correction is applied.".to_string()
            } else {
                "Lennard-Jones tail correction is disabled.".to_string()
            },
            format!("Framework net charge: {framework_charge:.3} e."),
            "Loading is molecules per unit cell.".to_string(),
        ],
    })
}

/// Largest real root of the monic cubic `x^3 + a x^2 + b x + c`.
fn largest_cubic_root(a: f64, b: f64, c: f64) -> Option<f64> {
    let p = b - a * a / 3.0;
    let q = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;
    let discriminant = 0.25 * q * q + p * p * p / 27.0;
    if discriminant > 0.0 {
        let sqrt_d = discriminant.sqrt();
        let u = (-0.5 * q + sqrt_d).cbrt();
        let v = (-0.5 * q - sqrt_d).cbrt();
        let root = u + v - a / 3.0;
        return root.is_finite().then_some(root);
    }
    if p.abs() < 1e-18 {
        return Some(-a / 3.0);
    }
    let radius = (-p / 3.0).sqrt();
    let argument = (-0.5 * q / (radius * radius * radius)).clamp(-1.0, 1.0);
    let phi = argument.acos();
    let mut best = f64::NEG_INFINITY;
    for k in 0..3 {
        let angle = (phi + 2.0 * std::f64::consts::PI * k as f64) / 3.0;
        let root = 2.0 * radius * angle.cos() - a / 3.0;
        if root > best {
            best = root;
        }
    }
    best.is_finite().then_some(best)
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
    grid: Option<&FrameworkGrid>,
    ewald_host: Option<&EwaldHost>,
) -> f64 {
    lennard_jones_energy(framework, sites, guests, config, grid)
        + electrostatic_energy(framework, sites, guests, config, grid, ewald_host)
        + tail_energy(framework, sites, guests, config)
}

fn electrostatic_energy(
    framework: &GcmcFramework,
    sites: &[GuestSite],
    guests: &[GuestMolecule],
    config: &GcmcConfig,
    grid: Option<&FrameworkGrid>,
    ewald_host: Option<&EwaldHost>,
) -> f64 {
    match &config.electrostatics {
        GcmcElectrostatics::None => 0.0,
        GcmcElectrostatics::Cutoff => cutoff_coulomb_energy(framework, sites, guests, config, grid),
        GcmcElectrostatics::Ewald(_) => match ewald_host {
            Some(host) => ewald_energy_fast(framework, sites, guests, config, host),
            None => ewald_energy(framework, sites, guests, config),
        },
    }
}

fn cutoff_coulomb_energy(
    framework: &GcmcFramework,
    sites: &[GuestSite],
    guests: &[GuestMolecule],
    config: &GcmcConfig,
    grid: Option<&FrameworkGrid>,
) -> f64 {
    let mut energy = 0.0;
    // Host-host pairs are constant during a GCMC run and are intentionally
    // omitted; only guest-host and guest-guest terms change an MC move.
    for molecule in guests {
        let positions = molecule.site_positions(sites);
        for (site, position) in sites.iter().zip(positions) {
            if site.charge == 0.0 {
                continue;
            }
            if let Some(grid) = grid {
                grid.for_each_nearby(position, |index| {
                    let atom = &framework.atoms[index];
                    if atom.charge == 0.0 {
                        return;
                    }
                    let r = framework
                        .cell
                        .minimum_image_distance(position, atom.position);
                    if r > 1e-9 && r < config.cutoff_angstrom {
                        energy += switched_coulomb(
                            site.charge,
                            atom.charge,
                            r,
                            config.cutoff_angstrom,
                            config.switching,
                        );
                    }
                });
            } else {
                for atom in &framework.atoms {
                    if atom.charge == 0.0 {
                        continue;
                    }
                    let r = framework
                        .cell
                        .minimum_image_distance(position, atom.position);
                    if r > 1e-9 && r < config.cutoff_angstrom {
                        energy += switched_coulomb(
                            site.charge,
                            atom.charge,
                            r,
                            config.cutoff_angstrom,
                            config.switching,
                        );
                    }
                }
            }
        }
    }
    for i in 0..guests.len() {
        let positions_i = guests[i].site_positions(sites);
        for j in (i + 1)..guests.len() {
            let positions_j = guests[j].site_positions(sites);
            for (site_i, position_i) in sites.iter().zip(&positions_i) {
                if site_i.charge == 0.0 {
                    continue;
                }
                for (site_j, position_j) in sites.iter().zip(&positions_j) {
                    if site_j.charge == 0.0 {
                        continue;
                    }
                    let r = framework
                        .cell
                        .minimum_image_distance(*position_i, *position_j);
                    if r > 1e-9 && r < config.cutoff_angstrom {
                        energy += switched_coulomb(
                            site_i.charge,
                            site_j.charge,
                            r,
                            config.cutoff_angstrom,
                            config.switching,
                        );
                    }
                }
            }
        }
    }
    energy
}

/// Analytic Lennard-Jones energy beyond the cutoff, in kcal/mol.
///
/// For a pair with mixed `epsilon`/`sigma`, the integral of the 12-6 potential
/// over the isotropic volume from `cutoff` to infinity is
///
/// ```text
/// 16 pi epsilon [sigma^12 / (9 rc^9) - sigma^6 / (3 rc^3)].
/// ```
///
/// A host atom contributes a density `1/V`; a guest sees `N-1` other guests.
fn continuum_tail(epsilon: f64, sigma: f64, cutoff: f64) -> f64 {
    if epsilon <= 0.0 || sigma <= 0.0 || cutoff <= 0.0 {
        return 0.0;
    }
    let sigma6 = sigma.powi(6);
    let sigma12 = sigma6 * sigma6;
    let rc3 = cutoff.powi(3);
    let rc9 = rc3.powi(3);
    16.0 * std::f64::consts::PI * epsilon * (sigma12 / (9.0 * rc9) - sigma6 / (3.0 * rc3))
}

fn tail_energy(
    framework: &GcmcFramework,
    sites: &[GuestSite],
    guests: &[GuestMolecule],
    config: &GcmcConfig,
) -> f64 {
    if !config.tail_correction || guests.is_empty() {
        return 0.0;
    }
    let volume = framework.cell.volume();
    if volume <= 0.0 {
        return 0.0;
    }
    let mut energy = 0.0;
    let cutoff = config.cutoff_angstrom;

    // Guest-framework tail: each host atom has number density 1/V.
    for molecule in guests {
        let positions = molecule.site_positions(sites);
        for (site, _position) in sites.iter().zip(positions) {
            for atom in &framework.atoms {
                let sigma = 0.5 * (atom.sigma_angstrom + site.sigma_angstrom);
                let epsilon = (atom.epsilon_kcal_mol * site.epsilon_kcal_mol).sqrt();
                energy += continuum_tail(epsilon, sigma, cutoff) / volume;
            }
        }
    }

    // Guest-guest tail: each guest sees the other N-1 molecules.
    if guests.len() > 1 {
        let others = (guests.len() - 1) as f64;
        for _ in guests {
            for site in sites {
                for other in sites {
                    let sigma = 0.5 * (site.sigma_angstrom + other.sigma_angstrom);
                    let epsilon = (site.epsilon_kcal_mol * other.epsilon_kcal_mol).sqrt();
                    energy += 0.5 * others * continuum_tail(epsilon, sigma, cutoff) / volume;
                }
            }
        }
    }
    energy
}

fn lennard_jones_energy(
    framework: &GcmcFramework,
    sites: &[GuestSite],
    guests: &[GuestMolecule],
    config: &GcmcConfig,
    grid: Option<&FrameworkGrid>,
) -> f64 {
    let mut energy = 0.0;
    for molecule in guests {
        let positions = molecule.site_positions(sites);
        for (site, position) in sites.iter().zip(positions) {
            if let Some(grid) = grid {
                grid.for_each_nearby(position, |index| {
                    let atom = &framework.atoms[index];
                    let r = framework
                        .cell
                        .minimum_image_distance(position, atom.position);
                    if r < config.cutoff_angstrom {
                        let sigma = 0.5 * (atom.sigma_angstrom + site.sigma_angstrom);
                        let epsilon = (atom.epsilon_kcal_mol * site.epsilon_kcal_mol).sqrt();
                        energy += switched_lennard_jones(
                            r,
                            epsilon,
                            sigma,
                            config.cutoff_angstrom,
                            config.switching,
                            config.shift_lj,
                        );
                    }
                });
            } else {
                for atom in &framework.atoms {
                    let r = framework
                        .cell
                        .minimum_image_distance(position, atom.position);
                    if r < config.cutoff_angstrom {
                        let sigma = 0.5 * (atom.sigma_angstrom + site.sigma_angstrom);
                        let epsilon = (atom.epsilon_kcal_mol * site.epsilon_kcal_mol).sqrt();
                        energy += switched_lennard_jones(
                            r,
                            epsilon,
                            sigma,
                            config.cutoff_angstrom,
                            config.switching,
                            config.shift_lj,
                        );
                    }
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
                        energy += switched_lennard_jones(
                            r,
                            epsilon,
                            sigma,
                            config.cutoff_angstrom,
                            config.switching,
                            config.shift_lj,
                        );
                    }
                }
            }
        }
    }
    energy
}

/// Ewald energy with the host structure factor precomputed.
///
/// Host-host real-space and self terms are constant and omitted; the
/// reciprocal host factor is included so guest-host and guest-guest cross
/// terms remain exact.
fn intramolecular_coulomb(sites: &[GuestSite]) -> f64 {
    let mut energy = 0.0;
    for i in 0..sites.len() {
        for j in (i + 1)..sites.len() {
            let dx = sites[i].offset_angstrom[0] - sites[j].offset_angstrom[0];
            let dy = sites[i].offset_angstrom[1] - sites[j].offset_angstrom[1];
            let dz = sites[i].offset_angstrom[2] - sites[j].offset_angstrom[2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r > 1e-12 {
                energy += COULOMB_KCAL_ANGSTROM * sites[i].charge * sites[j].charge / r;
            }
        }
    }
    energy
}

fn ewald_energy_fast(
    framework: &GcmcFramework,
    sites: &[GuestSite],
    guests: &[GuestMolecule],
    config: &GcmcConfig,
    host: &EwaldHost,
) -> f64 {
    let GcmcElectrostatics::Ewald(ewald) = &config.electrostatics else {
        return 0.0;
    };
    let mut guest_charges = Vec::new();
    let mut guest_positions = Vec::new();
    for molecule in guests {
        let positions = molecule.site_positions(sites);
        for (site, position) in sites.iter().zip(positions) {
            if site.charge != 0.0 {
                guest_charges.push(site.charge);
                guest_positions.push(position);
            }
        }
    }
    if guest_charges.is_empty() {
        return 0.0;
    }
    let alpha = host.alpha;
    let mut energy = 0.0;
    // Real space: host-guest.
    for (guest_charge, guest_position) in guest_charges.iter().zip(&guest_positions) {
        for (host_charge, host_position) in &host.host_charges {
            let r = framework
                .cell
                .minimum_image_distance(*guest_position, *host_position);
            if r > 1e-12 && r < ewald.real_cutoff_angstrom {
                energy += guest_charge * host_charge * libm::erfc(alpha * r) / r
                    * if config.switching {
                        switch_factor(
                            r,
                            ewald.real_cutoff_angstrom,
                            0.65 * ewald.real_cutoff_angstrom,
                        )
                    } else {
                        1.0
                    };
            }
        }
    }
    // Real space: guest-guest.
    for i in 0..guest_charges.len() {
        for j in (i + 1)..guest_charges.len() {
            let r = framework
                .cell
                .minimum_image_distance(guest_positions[i], guest_positions[j]);
            if r > 1e-12 && r < ewald.real_cutoff_angstrom {
                energy += guest_charges[i] * guest_charges[j] * libm::erfc(alpha * r) / r
                    * if config.switching {
                        switch_factor(
                            r,
                            ewald.real_cutoff_angstrom,
                            0.65 * ewald.real_cutoff_angstrom,
                        )
                    } else {
                        1.0
                    };
            }
        }
    }
    // Reciprocal space: S = S_host + S_guest.
    for (k, factor, host_real, host_imag) in &host.k_vectors {
        let mut real = *host_real;
        let mut imag = *host_imag;
        for (charge, position) in guest_charges.iter().zip(&guest_positions) {
            let phase = k[0] * position[0] + k[1] * position[1] + k[2] * position[2];
            real += charge * phase.cos();
            imag += charge * phase.sin();
        }
        energy += factor * (real * real + imag * imag);
    }
    // Guest self energy and the full background. Host self is constant.
    let guest_self = -alpha / std::f64::consts::PI.sqrt()
        * guest_charges
            .iter()
            .map(|charge| charge * charge)
            .sum::<f64>();
    let guest_charge: f64 = guest_charges.iter().sum();
    let total_charge = host.host_charge + guest_charge;
    let background =
        -std::f64::consts::PI / (2.0 * host.volume * alpha * alpha) * total_charge * total_charge;
    COULOMB_KCAL_ANGSTROM * (energy - host.host_reciprocal + guest_self + background)
        - guests.len() as f64 * intramolecular_coulomb(sites)
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

fn switch_factor(r: f64, cutoff: f64, switch_start: f64) -> f64 {
    if switch_start <= 0.0 || switch_start >= cutoff {
        return 1.0;
    }
    if r <= switch_start {
        return 1.0;
    }
    if r >= cutoff {
        return 0.0;
    }
    let rc2 = cutoff * cutoff;
    let rs2 = switch_start * switch_start;
    let r2 = r * r;
    let denominator = (rc2 - rs2).powi(3);
    ((rc2 - r2).powi(2) * (rc2 + 2.0 * r2 - 3.0 * rs2) / denominator).clamp(0.0, 1.0)
}

fn switched_lennard_jones(
    r: f64,
    epsilon: f64,
    sigma: f64,
    cutoff: f64,
    switching: bool,
    shift: bool,
) -> f64 {
    let value = lennard_jones(r, epsilon, sigma);
    if shift {
        value - lennard_jones(cutoff, epsilon, sigma)
    } else if switching {
        value * switch_factor(r, cutoff, 0.9 * cutoff)
    } else {
        value
    }
}

fn switched_coulomb(
    first_charge: f64,
    second_charge: f64,
    r: f64,
    cutoff: f64,
    switching: bool,
) -> f64 {
    let value = COULOMB_KCAL_ANGSTROM * first_charge * second_charge / r;
    if switching {
        value * switch_factor(r, cutoff, 0.65 * cutoff)
    } else {
        value
    }
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
            (loading - ideal).abs() < 0.15 * ideal.max(1.0),
            "loading = {loading}, ideal = {ideal}"
        );
    }

    #[test]
    fn peng_robinson_fugacity_brackets_the_ideal_gas() {
        let guest =
            GuestModel::single_site(0.0, 3.0, 0.0).with_critical_constants(304.13, 73.77, 0.22394);
        let low = guest.fugacity_coefficient(298.15, 1.0).unwrap();
        let high = guest.fugacity_coefficient(298.15, 50.0).unwrap();
        assert!(low > 0.99 && low < 1.0, "phi(1 bar) = {low}");
        assert!(high > 0.70 && high < 0.80, "phi(50 bar) = {high}");
        assert!(high < low);
        assert_eq!(guest.fugacity_bar(298.15, 0.0).unwrap(), 0.0);
    }

    #[test]
    fn missing_critical_constants_keep_the_ideal_gas_limit() {
        let guest = GuestModel::single_site(0.0, 3.0, 0.0);
        assert_eq!(guest.fugacity_coefficient(300.0, 50.0).unwrap(), 1.0);
    }

    #[test]
    fn continuum_tail_matches_the_numerical_integral() {
        let epsilon = 0.5;
        let sigma = 3.2;
        let cutoff = 9.0;
        let analytic = continuum_tail(epsilon, sigma, cutoff);
        let steps = 200_000;
        let upper = 80.0;
        let dr = (upper - cutoff) / steps as f64;
        let mut numerical = 0.0;
        for i in 0..steps {
            let r = cutoff + (i as f64 + 0.5) * dr;
            let ratio6 = (sigma / r).powi(6);
            numerical += 4.0
                * std::f64::consts::PI
                * r
                * r
                * 4.0
                * epsilon
                * (ratio6 * ratio6 - ratio6)
                * dr;
        }
        assert!(
            (analytic - numerical).abs() < 2e-3 * analytic.abs().max(1.0),
            "analytic = {analytic}, numerical = {numerical}"
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
            equilibration_steps: 200_000,
            production_steps: 1_000_000,
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
            ..GuestModel::single_site(0.0, 1.0, 0.0)
        };
        assert!((guest.total_charge()).abs() < 1e-12);
        assert_eq!(guest.interaction_sites().len(), 2);
    }

    #[test]
    fn framework_grid_matches_the_all_atom_loop() {
        let framework = GcmcFramework {
            cell: UnitCell::cubic(20.0),
            atoms: vec![
                LennardJonesAtom {
                    position: [2.0, 3.0, 4.0],
                    epsilon_kcal_mol: 0.2,
                    sigma_angstrom: 3.1,
                    charge: 0.4,
                },
                LennardJonesAtom {
                    position: [18.0, 7.0, 5.0],
                    epsilon_kcal_mol: 0.3,
                    sigma_angstrom: 2.9,
                    charge: -0.2,
                },
                LennardJonesAtom {
                    position: [9.0, 15.0, 11.0],
                    epsilon_kcal_mol: 0.1,
                    sigma_angstrom: 3.4,
                    charge: 0.1,
                },
            ],
        };
        let guest = GuestModel::single_site(0.2, 3.0, -0.3);
        let guests = vec![
            GuestMolecule::new([3.5, 4.0, 5.5]),
            GuestMolecule::new([12.0, 6.0, 8.0]),
            GuestMolecule::new([19.0, 18.0, 2.0]),
        ];
        let config = GcmcConfig {
            cutoff_angstrom: 8.0,
            electrostatics: GcmcElectrostatics::Cutoff,
            ..GcmcConfig::default()
        };
        let grid = FrameworkGrid::new(&framework, config.cutoff_angstrom / 2.0);
        assert!(grid.is_some());
        let sites = guest.interaction_sites();
        let with_grid = total_energy(&framework, &sites, &guests, &config, grid.as_ref(), None);
        let without_grid = total_energy(&framework, &sites, &guests, &config, None, None);
        assert!(
            (with_grid - without_grid).abs() < 1e-12,
            "grid={with_grid}, all-atom={without_grid}"
        );
    }

    #[test]
    fn ewald_host_structure_factor_matches_full_energy_differences() {
        let framework = GcmcFramework {
            cell: UnitCell::cubic(18.0),
            atoms: vec![
                LennardJonesAtom {
                    position: [4.0, 4.0, 4.0],
                    epsilon_kcal_mol: 0.0,
                    sigma_angstrom: 1.0,
                    charge: 1.0,
                },
                LennardJonesAtom {
                    position: [14.0, 13.0, 6.0],
                    epsilon_kcal_mol: 0.0,
                    sigma_angstrom: 1.0,
                    charge: -1.0,
                },
            ],
        };
        let guest = GuestModel::single_site(0.0, 1.0, -0.5);
        let config = GcmcConfig {
            electrostatics: GcmcElectrostatics::Ewald(EwaldConfig {
                alpha_angstrom_inv: 0.3,
                k_max: 2,
                real_cutoff_angstrom: 8.0,
            }),
            ..GcmcConfig::default()
        };
        let GcmcElectrostatics::Ewald(ewald) = &config.electrostatics else {
            unreachable!()
        };
        let host = EwaldHost::new(&framework, ewald);
        let sites = guest.interaction_sites();
        let first = vec![GuestMolecule::new([5.0, 5.0, 5.0])];
        let second = vec![
            GuestMolecule::new([5.0, 5.0, 5.0]),
            GuestMolecule::new([12.0, 7.0, 9.0]),
        ];
        let full_delta = ewald_energy(&framework, &sites, &second, &config)
            - ewald_energy(&framework, &sites, &first, &config);
        let fast_delta = ewald_energy_fast(&framework, &sites, &second, &config, &host)
            - ewald_energy_fast(&framework, &sites, &first, &config, &host);
        assert!(
            (full_delta - fast_delta).abs() < 1e-10,
            "full={full_delta}, fast={fast_delta}"
        );
    }

    #[test]
    fn ewald_host_structure_factor_matches_full_energy_nonorthogonal() {
        let framework = GcmcFramework {
            cell: UnitCell::new([
                [18.627314, 0.0, 0.0],
                [9.313657, 16.131727, 0.0],
                [9.313657, 5.377242, 15.209138],
            ]),
            atoms: vec![
                LennardJonesAtom {
                    position: [1.248, 0.0, 8.696],
                    epsilon_kcal_mol: 0.0,
                    sigma_angstrom: 1.0,
                    charge: 1.0,
                },
                LennardJonesAtom {
                    position: [4.0, 3.0, 5.0],
                    epsilon_kcal_mol: 0.0,
                    sigma_angstrom: 1.0,
                    charge: -0.8,
                },
                LennardJonesAtom {
                    position: [10.0, 8.0, 2.0],
                    epsilon_kcal_mol: 0.0,
                    sigma_angstrom: 1.0,
                    charge: -0.2,
                },
            ],
        };
        let guest = GuestModel::single_site(0.0, 1.0, -0.5);
        let config = GcmcConfig {
            electrostatics: GcmcElectrostatics::Ewald(EwaldConfig {
                alpha_angstrom_inv: 0.3,
                k_max: 2,
                real_cutoff_angstrom: 8.0,
            }),
            ..GcmcConfig::default()
        };
        let GcmcElectrostatics::Ewald(ewald) = &config.electrostatics else {
            unreachable!()
        };
        let host = EwaldHost::new(&framework, ewald);
        let sites = guest.interaction_sites();
        let first = vec![GuestMolecule::new([5.0, 5.0, 5.0])];
        let second = vec![
            GuestMolecule::new([5.0, 5.0, 5.0]),
            GuestMolecule::new([12.0, 7.0, 9.0]),
        ];
        let full_delta = ewald_energy(&framework, &sites, &second, &config)
            - ewald_energy(&framework, &sites, &first, &config);
        let fast_delta = ewald_energy_fast(&framework, &sites, &second, &config, &host)
            - ewald_energy_fast(&framework, &sites, &first, &config, &host);
        assert!(
            (full_delta - fast_delta).abs() < 1e-9,
            "full={full_delta}, fast={fast_delta}"
        );
    }

    #[test]
    fn ewald_host_absolute_guest_energy_matches_full_minus_host() {
        let framework = GcmcFramework {
            cell: UnitCell::new([
                [18.627314, 0.0, 0.0],
                [9.313657, 16.131727, 0.0],
                [9.313657, 5.377242, 15.209138],
            ]),
            atoms: vec![
                LennardJonesAtom {
                    position: [1.248, 0.0, 8.696],
                    epsilon_kcal_mol: 0.0,
                    sigma_angstrom: 1.0,
                    charge: 1.0,
                },
                LennardJonesAtom {
                    position: [4.0, 3.0, 5.0],
                    epsilon_kcal_mol: 0.0,
                    sigma_angstrom: 1.0,
                    charge: -0.8,
                },
                LennardJonesAtom {
                    position: [10.0, 8.0, 2.0],
                    epsilon_kcal_mol: 0.0,
                    sigma_angstrom: 1.0,
                    charge: -0.2,
                },
            ],
        };
        let guest = GuestModel::single_site(0.0, 1.0, -0.5);
        let config = GcmcConfig {
            electrostatics: GcmcElectrostatics::Ewald(EwaldConfig {
                alpha_angstrom_inv: 0.3,
                k_max: 2,
                real_cutoff_angstrom: 8.0,
            }),
            ..GcmcConfig::default()
        };
        let GcmcElectrostatics::Ewald(ewald) = &config.electrostatics else {
            unreachable!()
        };
        let host = EwaldHost::new(&framework, ewald);
        let sites = guest.interaction_sites();
        let guests = vec![
            GuestMolecule::new([5.0, 5.0, 5.0]),
            GuestMolecule::new([12.0, 7.0, 9.0]),
        ];
        let host_only = ewald_energy(&framework, &sites, &[], &config);
        let full = ewald_energy(&framework, &sites, &guests, &config);
        let fast = ewald_energy_fast(&framework, &sites, &guests, &config, &host);
        assert!(
            (fast - (full - host_only)).abs() < 1e-9,
            "fast={fast}, full-host={}",
            full - host_only
        );
    }

    #[test]
    fn ewald_opposite_charges_match_direct_coulomb() {
        let cell = UnitCell::cubic(20.0);
        let framework = GcmcFramework {
            cell: cell.clone(),
            atoms: vec![LennardJonesAtom {
                position: [5.0, 5.0, 5.0],
                epsilon_kcal_mol: 0.0,
                sigma_angstrom: 1.0,
                charge: 1.0,
            }],
        };
        let guest = GuestModel::single_site(0.0, 1.0, -1.0);
        let guests = vec![GuestMolecule::new([6.0, 5.0, 5.0])];
        let config = GcmcConfig {
            electrostatics: GcmcElectrostatics::Ewald(EwaldConfig {
                alpha_angstrom_inv: 0.3,
                k_max: 6,
                real_cutoff_angstrom: 10.0,
            }),
            ..GcmcConfig::default()
        };
        let energy = ewald_energy(&framework, &guest.interaction_sites(), &guests, &config);
        let direct = -COULOMB_KCAL_ANGSTROM;
        assert!(
            (energy - direct).abs() < 0.01 * direct.abs(),
            "Ewald={energy}, direct={direct}"
        );

        let cutoff_config = GcmcConfig {
            electrostatics: GcmcElectrostatics::Cutoff,
            ..config
        };
        let cutoff = electrostatic_energy(
            &framework,
            &guest.interaction_sites(),
            &guests,
            &cutoff_config,
            None,
            None,
        );
        assert!(
            (cutoff - direct).abs() < 1e-12,
            "cutoff={cutoff}, direct={direct}"
        );
    }
}
