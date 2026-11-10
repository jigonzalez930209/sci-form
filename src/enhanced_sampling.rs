//! Well-tempered metadynamics and collective-variable sampling.
//!
//! This module adds a bias on a collective variable (CV) to a generic
//! potential. The first implementation covers distance, angle, dihedral and
//! coordinate CVs, the well-tempered Gaussian bias, a Langevin thermostat and
//! a one-dimensional free-energy surface. OPES and replica exchange are
//! documented follow-ups.
//!
//! # Units
//!
//! * coordinates in **Å**;
//! * masses in **amu**;
//! * energies in **kcal/mol**;
//! * time in **fs**;
//! * temperature in **K**;
//! * CVs are in Å for distances and radians for angles and dihedrals.
//!
//! # Example
//!
//! ```
//! use sci_form::enhanced_sampling::{
//!     CollectiveVariable, MetaDynamicsConfig, run_metadynamics, DoubleWell,
//! };
//!
//! let potential = DoubleWell::new(5.0, 1.0);
//! let config = MetaDynamicsConfig {
//!     cv: CollectiveVariable::Coordinate { index: 0 },
//!     max_steps: 2_000,
//!     ..MetaDynamicsConfig::default()
//! };
//! let result = run_metadynamics(&potential, &[-1.0], &[1.0], &config).unwrap();
//! assert_eq!(result.cv_trace.len(), 2_001);
//! ```

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

use crate::error::SciFormError;
use crate::forcefield::traits::MolecularForceField;
use crate::graph;

const AMU_ANGFS2_TO_KCAL_MOL: f64 = 2_390.057_361_533_49;
const R_GAS_KCAL_MOLK: f64 = 0.001_987_204_258_640_83;

/// A potential that returns an energy and a Cartesian gradient.
///
/// `coords` and `gradient` are flat `[x0, y0, z0, x1, ...]` arrays in Å and
/// kcal/mol/Å. The implementation must fill every gradient component.
pub trait Potential {
    /// Potential energy in kcal/mol.
    fn energy(&self, coords: &[f64]) -> f64;

    /// Gradient in kcal/mol/Å.
    fn gradient(&self, coords: &[f64], gradient: &mut [f64]);
}

/// The collective variable that is biased.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum CollectiveVariable {
    /// A single Cartesian component. Useful for one-dimensional models.
    Coordinate { index: usize },
    /// Distance between atoms `i` and `j`, in Å.
    Distance { i: usize, j: usize },
    /// Angle `i-j-k`, in radians.
    Angle { i: usize, j: usize, k: usize },
    /// Dihedral `i-j-k-l`, in radians.
    Dihedral {
        i: usize,
        j: usize,
        k: usize,
        l: usize,
    },
}

impl CollectiveVariable {
    /// Number of atom indices the CV needs.
    pub fn atom_count(&self) -> usize {
        match self {
            Self::Coordinate { .. } => 1,
            Self::Distance { .. } => 2,
            Self::Angle { .. } => 3,
            Self::Dihedral { .. } => 4,
        }
    }

    /// Evaluate the CV from flat coordinates.
    pub fn value(&self, coords: &[f64]) -> Result<f64, SciFormError> {
        match *self {
            Self::Coordinate { index } => coords
                .get(index)
                .copied()
                .ok_or_else(|| SciFormError::validation("CV coordinate index is out of range")),
            Self::Distance { i, j } => {
                let a = atom(coords, i)?;
                let b = atom(coords, j)?;
                Ok(distance(a, b))
            }
            Self::Angle { i, j, k } => {
                let a = atom(coords, i)?;
                let b = atom(coords, j)?;
                let c = atom(coords, k)?;
                Ok(angle(a, b, c))
            }
            Self::Dihedral { i, j, k, l } => {
                let a = atom(coords, i)?;
                let b = atom(coords, j)?;
                let c = atom(coords, k)?;
                let d = atom(coords, l)?;
                Ok(dihedral(a, b, c, d))
            }
        }
    }

    /// Numerical Cartesian gradient of the CV, in CV units per Å.
    ///
    /// The finite-difference step is 1e-5 Å. Distance uses the analytic
    /// gradient; the angular CVs use the central difference because their
    /// analytic forms are easy to get wrong at the linear limits.
    pub fn gradient(&self, coords: &[f64]) -> Result<Vec<f64>, SciFormError> {
        match *self {
            Self::Coordinate { index } => {
                if index >= coords.len() {
                    return Err(SciFormError::validation(
                        "CV coordinate index is out of range",
                    ));
                }
                let mut gradient = vec![0.0; coords.len()];
                gradient[index] = 1.0;
                Ok(gradient)
            }
            Self::Distance { i, j } => {
                let a = atom(coords, i)?;
                let b = atom(coords, j)?;
                let r = distance(a, b);
                if r < 1e-12 {
                    return Err(SciFormError::numerical(
                        "distance CV is undefined for coincident atoms",
                    ));
                }
                let mut gradient = vec![0.0; coords.len()];
                for (axis, (ai, bi)) in a.iter().zip(b.iter()).enumerate() {
                    let value = (ai - bi) / r;
                    gradient[3 * i + axis] = value;
                    gradient[3 * j + axis] = -value;
                }
                Ok(gradient)
            }
            Self::Angle { .. } | Self::Dihedral { .. } => {
                let mut gradient = vec![0.0; coords.len()];
                let step = 1e-5;
                for index in 0..coords.len() {
                    let mut plus = coords.to_vec();
                    let mut minus = coords.to_vec();
                    plus[index] += step;
                    minus[index] -= step;
                    let value = (self.value(&plus)? - self.value(&minus)?) / (2.0 * step);
                    gradient[index] = value;
                }
                Ok(gradient)
            }
        }
    }
}

/// A deposited Gaussian.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaussianDeposit {
    /// Integration step at which the Gaussian was deposited.
    pub step: usize,
    /// CV centre.
    pub center: f64,
    /// Gaussian height in kcal/mol.
    pub height_kcal_mol: f64,
    /// Gaussian width in CV units.
    pub sigma: f64,
    /// Importance weight used by the OPES density estimate.
    #[serde(default = "default_weight")]
    pub weight: f64,
}

fn default_weight() -> f64 {
    1.0
}

/// One point on the CV trace.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetaDynamicsPoint {
    /// Integration step.
    pub step: usize,
    /// CV value.
    pub cv: f64,
    /// Unbiased potential energy in kcal/mol.
    pub potential_kcal_mol: f64,
    /// Bias energy in kcal/mol.
    pub bias_kcal_mol: f64,
}

/// A one-dimensional free-energy surface reconstructed from the bias.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FreeEnergySurface {
    /// Grid lower bound.
    pub min: f64,
    /// Grid upper bound.
    pub max: f64,
    /// Number of bins.
    pub bins: usize,
    /// Free energy in kcal/mol, shifted so the minimum is zero.
    pub values: Vec<f64>,
}

impl FreeEnergySurface {
    /// Grid coordinate for a bin.
    pub fn coordinate(&self, index: usize) -> f64 {
        if self.bins <= 1 {
            return self.min;
        }
        self.min + (self.max - self.min) * index as f64 / (self.bins - 1) as f64
    }
}

/// Configuration for a well-tempered metadynamics run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetaDynamicsConfig {
    /// Biased collective variable.
    pub cv: CollectiveVariable,
    /// Initial Gaussian height in kcal/mol.
    pub height_kcal_mol: f64,
    /// Gaussian width in CV units.
    pub sigma: f64,
    /// Well-tempered bias factor `γ > 1`.
    pub bias_factor: f64,
    /// Deposit a Gaussian every `pace_steps` integration steps.
    pub pace_steps: usize,
    /// Temperature in K.
    pub temperature_k: f64,
    /// Integration time step in fs.
    pub dt_fs: f64,
    /// Langevin friction in ps^-1.
    pub friction_ps_inv: f64,
    /// Number of integration steps.
    pub max_steps: usize,
    /// RNG seed.
    pub seed: u64,
    /// Free-energy grid lower bound.
    pub grid_min: f64,
    /// Free-energy grid upper bound.
    pub grid_max: f64,
    /// Number of free-energy bins.
    pub grid_bins: usize,
    /// Store every `trace_stride` steps in the trace.
    pub trace_stride: usize,
}

impl Default for MetaDynamicsConfig {
    fn default() -> Self {
        Self {
            cv: CollectiveVariable::Distance { i: 0, j: 1 },
            height_kcal_mol: 0.2,
            sigma: 0.1,
            bias_factor: 10.0,
            pace_steps: 100,
            temperature_k: 300.0,
            dt_fs: 0.5,
            friction_ps_inv: 1.0,
            max_steps: 10_000,
            seed: 42,
            grid_min: 0.0,
            grid_max: 10.0,
            grid_bins: 100,
            trace_stride: 1,
        }
    }
}

impl MetaDynamicsConfig {
    fn validate(&self, coords_len: usize, masses_len: usize) -> Result<(), SciFormError> {
        if coords_len == 0 {
            return Err(SciFormError::validation("coords must not be empty"));
        }
        if masses_len != coords_len {
            return Err(SciFormError::validation(
                "masses must contain one value per Cartesian coordinate",
            ));
        }
        if !matches!(self.cv, CollectiveVariable::Coordinate { .. })
            && !coords_len.is_multiple_of(3)
        {
            return Err(SciFormError::validation(
                "atom-indexed CVs require coords to be a multiple of three",
            ));
        }
        if !self.height_kcal_mol.is_finite() || self.height_kcal_mol <= 0.0 {
            return Err(SciFormError::validation(
                "height_kcal_mol must be finite and > 0",
            ));
        }
        if !self.sigma.is_finite() || self.sigma <= 0.0 {
            return Err(SciFormError::validation("sigma must be finite and > 0"));
        }
        if !self.bias_factor.is_finite() || self.bias_factor <= 1.0 {
            return Err(SciFormError::validation(
                "bias_factor must be finite and > 1",
            ));
        }
        if self.pace_steps == 0 {
            return Err(SciFormError::validation("pace_steps must be >= 1"));
        }
        if !self.temperature_k.is_finite() || self.temperature_k <= 0.0 {
            return Err(SciFormError::validation(
                "temperature_k must be finite and > 0",
            ));
        }
        if !self.dt_fs.is_finite() || self.dt_fs <= 0.0 {
            return Err(SciFormError::validation("dt_fs must be finite and > 0"));
        }
        if !self.friction_ps_inv.is_finite() || self.friction_ps_inv < 0.0 {
            return Err(SciFormError::validation(
                "friction_ps_inv must be finite and >= 0",
            ));
        }
        if self.max_steps == 0 {
            return Err(SciFormError::validation("max_steps must be >= 1"));
        }
        if !self.grid_min.is_finite()
            || !self.grid_max.is_finite()
            || self.grid_max <= self.grid_min
        {
            return Err(SciFormError::validation(
                "grid_max must be finite and greater than grid_min",
            ));
        }
        if self.grid_bins < 2 {
            return Err(SciFormError::validation("grid_bins must be >= 2"));
        }
        if self.trace_stride == 0 {
            return Err(SciFormError::validation("trace_stride must be >= 1"));
        }
        let n_atoms = coords_len / 3;
        let check = |index: usize| {
            if index >= n_atoms {
                return Err(SciFormError::validation(format!(
                    "CV atom index {index} is out of range for {n_atoms} atoms"
                )));
            }
            Ok(())
        };
        match self.cv {
            CollectiveVariable::Coordinate { .. } => {}
            CollectiveVariable::Distance { i, j } => {
                check(i)?;
                check(j)?;
            }
            CollectiveVariable::Angle { i, j, k } => {
                check(i)?;
                check(j)?;
                check(k)?;
            }
            CollectiveVariable::Dihedral { i, j, k, l } => {
                check(i)?;
                check(j)?;
                check(k)?;
                check(l)?;
            }
        }
        Ok(())
    }
}

/// Result of a metadynamics run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetaDynamicsResult {
    /// Biased collective variable.
    pub cv: CollectiveVariable,
    /// Number of integration steps.
    pub steps: usize,
    /// Integration time step in fs.
    pub dt_fs: f64,
    /// Sampled CV points.
    pub cv_trace: Vec<MetaDynamicsPoint>,
    /// Deposited Gaussians.
    pub gaussians: Vec<GaussianDeposit>,
    /// Reconstructed free-energy surface.
    pub free_energy: FreeEnergySurface,
    /// Histogram-reweighted free-energy surface.
    pub reweighted_free_energy: FreeEnergySurface,
    /// Final flat coordinates in Å.
    pub final_coords: Vec<f64>,
    /// Final CV value.
    pub final_cv: f64,
    /// Human-readable notes and caveats.
    pub notes: Vec<String>,
}

/// Run well-tempered metadynamics on a generic potential.
///
/// The integrator is velocity Verlet with an Ornstein-Uhlenbeck Langevin
/// thermostat. The bias is a sum of Gaussians; its height decays as
/// `h0 exp(-V(s)/(R ΔT))` with `ΔT = (γ - 1) T`.
pub fn run_metadynamics<P: Potential>(
    potential: &P,
    initial_coords: &[f64],
    masses_amu: &[f64],
    config: &MetaDynamicsConfig,
) -> Result<MetaDynamicsResult, SciFormError> {
    config.validate(initial_coords.len(), masses_amu.len())?;
    let n_coords = initial_coords.len();
    let mut x = initial_coords.to_vec();
    if !x.iter().all(|value| value.is_finite()) {
        return Err(SciFormError::validation(
            "initial_coords contain a non-finite value",
        ));
    }
    if !masses_amu
        .iter()
        .all(|mass| mass.is_finite() && *mass > 0.0)
    {
        return Err(SciFormError::validation("masses must be finite and > 0"));
    }

    let mut rng = StdRng::seed_from_u64(config.seed);
    let mut v = vec![0.0; n_coords];
    for index in 0..n_coords {
        let sigma = (R_GAS_KCAL_MOLK * config.temperature_k
            / (masses_amu[index] * AMU_ANGFS2_TO_KCAL_MOL))
            .sqrt();
        v[index] = sigma * sample_standard_normal(&mut rng);
    }

    let mut gradient = vec![0.0; n_coords];
    let mut gaussians = Vec::new();
    let mut trace = Vec::with_capacity(config.max_steps / config.trace_stride + 1);

    let mut cv_value = config.cv.value(&x)?;
    let mut bias = 0.0;
    let mut potential_energy = potential.energy(&x);
    trace.push(MetaDynamicsPoint {
        step: 0,
        cv: cv_value,
        potential_kcal_mol: potential_energy,
        bias_kcal_mol: bias,
    });

    deposit_gaussian(&mut gaussians, 0, cv_value, config, bias)?;

    let friction_fs_inv = config.friction_ps_inv * 1e-3;
    let kbt = R_GAS_KCAL_MOLK * config.temperature_k;

    for step in 1..=config.max_steps {
        potential.gradient(&x, &mut gradient);
        let cv_gradient = config.cv.gradient(&x)?;
        let (_, bias_derivative) = bias_and_derivative(&gaussians, cv_value, config.sigma);

        let mut force = vec![0.0; n_coords];
        for index in 0..n_coords {
            force[index] = -gradient[index] - bias_derivative * cv_gradient[index];
        }

        for index in 0..n_coords {
            let acceleration = force[index] / (masses_amu[index] * AMU_ANGFS2_TO_KCAL_MOL);
            v[index] += 0.5 * acceleration * config.dt_fs;
            x[index] += v[index] * config.dt_fs;
        }

        potential.gradient(&x, &mut gradient);
        cv_value = config.cv.value(&x)?;
        let (bias_energy, bias_derivative) =
            bias_and_derivative(&gaussians, cv_value, config.sigma);
        bias = bias_energy;
        let cv_gradient = config.cv.gradient(&x)?;
        for index in 0..n_coords {
            force[index] = -gradient[index] - bias_derivative * cv_gradient[index];
        }
        for index in 0..n_coords {
            let acceleration = force[index] / (masses_amu[index] * AMU_ANGFS2_TO_KCAL_MOL);
            v[index] += 0.5 * acceleration * config.dt_fs;
        }

        // Ornstein-Uhlenbeck thermostat: exact velocity update for the
        // harmonic bath, applied per Cartesian component.
        if friction_fs_inv > 0.0 {
            let decay = (-friction_fs_inv * config.dt_fs).exp();
            let noise_scale = (kbt * (1.0 - decay * decay)).sqrt();
            for index in 0..n_coords {
                let scale = noise_scale / (masses_amu[index] * AMU_ANGFS2_TO_KCAL_MOL).sqrt();
                v[index] = decay * v[index] + scale * sample_standard_normal(&mut rng);
            }
        }

        if step % config.pace_steps == 0 {
            deposit_gaussian(&mut gaussians, step, cv_value, config, bias)?;
        }
        if step % config.trace_stride == 0 {
            potential_energy = potential.energy(&x);
            trace.push(MetaDynamicsPoint {
                step,
                cv: cv_value,
                potential_kcal_mol: potential_energy,
                bias_kcal_mol: bias,
            });
        }
        if !x.iter().all(|value| value.is_finite())
            || !v.iter().all(|value| value.is_finite())
            || !potential_energy.is_finite()
        {
            return Err(SciFormError::numerical(format!(
                "metadynamics became non-finite at step {step}"
            )));
        }
    }

    let free_energy = build_free_energy(&gaussians, config);
    let reweighted_free_energy = reweight_free_energy(
        &trace,
        &ReweightConfig {
            grid_min: config.grid_min,
            grid_max: config.grid_max,
            grid_bins: config.grid_bins,
            temperature_k: config.temperature_k,
            block_size: 1,
        },
    )?;
    let mut notes = vec![
        "Well-tempered metadynamics with a Gaussian bias and a Langevin thermostat.".to_string(),
        "The free-energy surface is the negative of the converged bias, shifted to its minimum."
            .to_string(),
        "OPES and replica exchange are follow-up implementations.".to_string(),
    ];
    if trace.len() < 10 {
        notes.push("Very short run: the free-energy surface is not converged.".to_string());
    }

    Ok(MetaDynamicsResult {
        cv: config.cv.clone(),
        steps: config.max_steps,
        dt_fs: config.dt_fs,
        cv_trace: trace,
        gaussians,
        free_energy,
        reweighted_free_energy,
        final_coords: x,
        final_cv: cv_value,
        notes,
    })
}

/// Configuration for histogram reweighting.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReweightConfig {
    /// Grid lower bound.
    pub grid_min: f64,
    /// Grid upper bound.
    pub grid_max: f64,
    /// Number of bins.
    pub grid_bins: usize,
    /// Temperature in K.
    pub temperature_k: f64,
    /// Number of consecutive samples per block for correlated reweighting.
    #[serde(default = "default_block_size")]
    pub block_size: usize,
}

fn default_block_size() -> usize {
    1
}

impl ReweightConfig {
    fn validate(&self) -> Result<(), SciFormError> {
        if !self.grid_min.is_finite()
            || !self.grid_max.is_finite()
            || self.grid_max <= self.grid_min
        {
            return Err(SciFormError::validation(
                "grid_max must be finite and greater than grid_min",
            ));
        }
        if self.grid_bins < 2 {
            return Err(SciFormError::validation("grid_bins must be >= 2"));
        }
        if !self.temperature_k.is_finite() || self.temperature_k <= 0.0 {
            return Err(SciFormError::validation(
                "temperature_k must be finite and > 0",
            ));
        }
        if self.block_size == 0 {
            return Err(SciFormError::validation("block_size must be >= 1"));
        }
        Ok(())
    }
}

/// Reweight a biased CV trace into an unbiased free-energy histogram.
///
/// Each sample is weighted by `exp(+V_bias/kBT)`, the standard importance
/// weight for a biased trajectory. The result is shifted so its minimum is zero.
pub fn reweight_free_energy(
    trace: &[MetaDynamicsPoint],
    config: &ReweightConfig,
) -> Result<FreeEnergySurface, SciFormError> {
    config.validate()?;
    if trace.is_empty() {
        return Err(SciFormError::validation(
            "cannot reweight an empty CV trace",
        ));
    }
    let bins = config.grid_bins;
    let width = (config.grid_max - config.grid_min) / (bins - 1) as f64;
    let mut weights = vec![0.0f64; bins];
    let kbt = R_GAS_KCAL_MOLK * config.temperature_k;
    for point in trace {
        if !point.cv.is_finite() || !point.bias_kcal_mol.is_finite() {
            return Err(SciFormError::validation(
                "CV trace contains a non-finite value",
            ));
        }
        let index = ((point.cv - config.grid_min) / width)
            .round()
            .clamp(0.0, (bins - 1) as f64) as usize;
        weights[index] += (point.bias_kcal_mol / kbt).exp();
    }
    if weights.iter().all(|weight| *weight <= 0.0) {
        return Err(SciFormError::numerical(
            "no CV sample fell inside the reweighting grid",
        ));
    }
    let mut values: Vec<f64> = weights
        .iter()
        .map(|weight| {
            if *weight > 0.0 {
                -kbt * weight.ln()
            } else {
                f64::INFINITY
            }
        })
        .collect();
    let minimum = values.iter().copied().fold(f64::INFINITY, f64::min);
    if minimum.is_finite() {
        for value in &mut values {
            *value -= minimum;
        }
    }
    Ok(FreeEnergySurface {
        min: config.grid_min,
        max: config.grid_max,
        bins,
        values,
    })
}

/// Result of a correlation-aware reweighting.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReweightResult {
    /// Reweighted free-energy surface.
    pub surface: FreeEnergySurface,
    /// Per-bin standard error in kcal/mol.
    pub standard_error: Vec<f64>,
    /// Number of blocks.
    pub n_blocks: usize,
    /// Samples per block.
    pub block_size: usize,
}

/// Reweight a biased trace with block-averaged uncertainty.
///
/// The surface is the same histogram estimate as [`reweight_free_energy`]; the
/// standard error comes from the block-to-block variance of the per-bin
/// weights. This is a simple correlation-aware diagnostic, not a full
/// autocorrelation analysis.
pub fn reweight_free_energy_correlated(
    trace: &[MetaDynamicsPoint],
    config: &ReweightConfig,
) -> Result<ReweightResult, SciFormError> {
    let surface = reweight_free_energy(trace, config)?;
    config.validate()?;
    let bins = config.grid_bins;
    let width = (config.grid_max - config.grid_min) / (bins - 1) as f64;
    let kbt = R_GAS_KCAL_MOLK * config.temperature_k;
    let block_size = config.block_size;
    let n_blocks = trace.len().div_ceil(block_size);
    let mut block_weights = vec![vec![0.0f64; bins]; n_blocks];
    for (index, point) in trace.iter().enumerate() {
        let block = index / block_size;
        let bin = ((point.cv - config.grid_min) / width)
            .round()
            .clamp(0.0, (bins - 1) as f64) as usize;
        block_weights[block][bin] += (point.bias_kcal_mol / kbt).exp();
    }
    let mut standard_error = vec![0.0f64; bins];
    for (bin, error) in standard_error.iter_mut().enumerate() {
        let values: Vec<f64> = block_weights.iter().map(|block| block[bin]).collect();
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        if mean > 0.0 && values.len() > 1 {
            let variance = values
                .iter()
                .map(|value| (value - mean).powi(2))
                .sum::<f64>()
                / (values.len() - 1) as f64;
            *error = kbt * (variance / values.len() as f64).sqrt() / mean;
        }
    }
    Ok(ReweightResult {
        surface,
        standard_error,
        n_blocks,
        block_size,
    })
}

/// Configuration for the OPES-style kernel-density bias.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OpesConfig {
    /// Biased collective variable.
    pub cv: CollectiveVariable,
    /// Kernel width in CV units.
    pub sigma: f64,
    /// Bias factor `γ > 1`.
    pub bias_factor: f64,
    /// Deposit a kernel every `pace_steps` steps.
    pub pace_steps: usize,
    /// Temperature in K.
    pub temperature_k: f64,
    /// Integration time step in fs.
    pub dt_fs: f64,
    /// Langevin friction in ps^-1.
    pub friction_ps_inv: f64,
    /// Number of integration steps.
    pub max_steps: usize,
    /// RNG seed.
    pub seed: u64,
    /// Grid lower bound.
    pub grid_min: f64,
    /// Grid upper bound.
    pub grid_max: f64,
    /// Number of grid bins.
    pub grid_bins: usize,
    /// Trace stride.
    pub trace_stride: usize,
    /// Regularization added to the density ratio.
    pub regularization: f64,
}

impl Default for OpesConfig {
    fn default() -> Self {
        Self {
            cv: CollectiveVariable::Distance { i: 0, j: 1 },
            sigma: 0.1,
            bias_factor: 10.0,
            pace_steps: 100,
            temperature_k: 300.0,
            dt_fs: 0.5,
            friction_ps_inv: 1.0,
            max_steps: 10_000,
            seed: 42,
            grid_min: 0.0,
            grid_max: 10.0,
            grid_bins: 100,
            trace_stride: 1,
            regularization: 1e-8,
        }
    }
}

impl OpesConfig {
    fn validate(&self, coords_len: usize, masses_len: usize) -> Result<(), SciFormError> {
        if coords_len == 0 || masses_len != coords_len {
            return Err(SciFormError::validation(
                "coords and masses must be non-empty and of equal length",
            ));
        }
        if !matches!(self.cv, CollectiveVariable::Coordinate { .. })
            && !coords_len.is_multiple_of(3)
        {
            return Err(SciFormError::validation(
                "atom-indexed CVs require coords to be a multiple of three",
            ));
        }
        if !self.sigma.is_finite() || self.sigma <= 0.0 {
            return Err(SciFormError::validation("sigma must be finite and > 0"));
        }
        if !self.bias_factor.is_finite() || self.bias_factor <= 1.0 {
            return Err(SciFormError::validation(
                "bias_factor must be finite and > 1",
            ));
        }
        if !self.regularization.is_finite() || self.regularization <= 0.0 {
            return Err(SciFormError::validation(
                "regularization must be finite and > 0",
            ));
        }
        if self.pace_steps == 0 || self.max_steps == 0 || self.trace_stride == 0 {
            return Err(SciFormError::validation(
                "pace_steps, max_steps and trace_stride must be >= 1",
            ));
        }
        if !self.temperature_k.is_finite() || self.temperature_k <= 0.0 {
            return Err(SciFormError::validation(
                "temperature_k must be finite and > 0",
            ));
        }
        if !self.dt_fs.is_finite() || self.dt_fs <= 0.0 {
            return Err(SciFormError::validation("dt_fs must be finite and > 0"));
        }
        if !self.friction_ps_inv.is_finite() || self.friction_ps_inv < 0.0 {
            return Err(SciFormError::validation(
                "friction_ps_inv must be finite and >= 0",
            ));
        }
        if !self.grid_min.is_finite()
            || !self.grid_max.is_finite()
            || self.grid_max <= self.grid_min
        {
            return Err(SciFormError::validation(
                "grid_max must be finite and greater than grid_min",
            ));
        }
        if self.grid_bins < 2 {
            return Err(SciFormError::validation("grid_bins must be >= 2"));
        }
        Ok(())
    }
}

/// Result of an OPES-style run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OpesResult {
    /// Biased CV.
    pub cv: CollectiveVariable,
    /// Integration steps.
    pub steps: usize,
    /// Time step in fs.
    pub dt_fs: f64,
    /// CV trace.
    pub cv_trace: Vec<MetaDynamicsPoint>,
    /// Deposited kernels.
    pub kernels: Vec<GaussianDeposit>,
    /// Bias-derived free-energy surface.
    pub free_energy: FreeEnergySurface,
    /// Histogram-reweighted free-energy surface.
    pub reweighted_free_energy: FreeEnergySurface,
    /// Final coordinates.
    pub final_coords: Vec<f64>,
    /// Final CV.
    pub final_cv: f64,
    /// Notes and caveats.
    pub notes: Vec<String>,
}

/// Run an OPES-style kernel-density bias.
///
/// The bias is
///
/// ```text
/// V(s) = -(1 - 1/γ) kBT ln(P(s)/P_max + ε),
/// ```
///
/// where `P` is a Gaussian kernel density estimate of the visited CV values.
/// This is the OPES probability-enhancement idea with a bounded-grid
/// normalisation; it is not a full implementation of the OPES flow equations.
pub fn run_opes<P: Potential>(
    potential: &P,
    initial_coords: &[f64],
    masses_amu: &[f64],
    config: &OpesConfig,
) -> Result<OpesResult, SciFormError> {
    config.validate(initial_coords.len(), masses_amu.len())?;
    let n_coords = initial_coords.len();
    let mut x = initial_coords.to_vec();
    let mut rng = StdRng::seed_from_u64(config.seed);
    let mut v = vec![0.0; n_coords];
    for index in 0..n_coords {
        let sigma = (R_GAS_KCAL_MOLK * config.temperature_k
            / (masses_amu[index] * AMU_ANGFS2_TO_KCAL_MOL))
            .sqrt();
        v[index] = sigma * sample_standard_normal(&mut rng);
    }
    let mut gradient = vec![0.0; n_coords];
    let mut kernels = Vec::new();
    let mut trace = Vec::with_capacity(config.max_steps / config.trace_stride + 1);
    let mut cv_value = config.cv.value(&x)?;
    let mut potential_energy = potential.energy(&x);
    trace.push(MetaDynamicsPoint {
        step: 0,
        cv: cv_value,
        potential_kcal_mol: potential_energy,
        bias_kcal_mol: 0.0,
    });
    let kbt = R_GAS_KCAL_MOLK * config.temperature_k;
    let friction_fs_inv = config.friction_ps_inv * 1e-3;
    let mut bias;

    for step in 1..=config.max_steps {
        potential.gradient(&x, &mut gradient);
        let cv_gradient = config.cv.gradient(&x)?;
        let (_, bias_derivative) = opes_bias_and_derivative(&kernels, cv_value, config, kbt);
        let mut force = vec![0.0; n_coords];
        for index in 0..n_coords {
            force[index] = -gradient[index] - bias_derivative * cv_gradient[index];
        }
        for index in 0..n_coords {
            let acceleration = force[index] / (masses_amu[index] * AMU_ANGFS2_TO_KCAL_MOL);
            v[index] += 0.5 * acceleration * config.dt_fs;
            x[index] += v[index] * config.dt_fs;
        }
        potential.gradient(&x, &mut gradient);
        cv_value = config.cv.value(&x)?;
        let (bias_energy, bias_derivative) =
            opes_bias_and_derivative(&kernels, cv_value, config, kbt);
        bias = bias_energy;
        let cv_gradient = config.cv.gradient(&x)?;
        for index in 0..n_coords {
            force[index] = -gradient[index] - bias_derivative * cv_gradient[index];
        }
        for index in 0..n_coords {
            let acceleration = force[index] / (masses_amu[index] * AMU_ANGFS2_TO_KCAL_MOL);
            v[index] += 0.5 * acceleration * config.dt_fs;
        }
        if friction_fs_inv > 0.0 {
            let decay = (-friction_fs_inv * config.dt_fs).exp();
            let noise_scale = (kbt * (1.0 - decay * decay)).sqrt();
            for index in 0..n_coords {
                let scale = noise_scale / (masses_amu[index] * AMU_ANGFS2_TO_KCAL_MOL).sqrt();
                v[index] = decay * v[index] + scale * sample_standard_normal(&mut rng);
            }
        }
        if step % config.pace_steps == 0 {
            kernels.push(GaussianDeposit {
                step,
                center: cv_value,
                height_kcal_mol: 1.0,
                sigma: config.sigma,
                weight: (bias / kbt).exp(),
            });
        }
        if step % config.trace_stride == 0 {
            potential_energy = potential.energy(&x);
            trace.push(MetaDynamicsPoint {
                step,
                cv: cv_value,
                potential_kcal_mol: potential_energy,
                bias_kcal_mol: bias,
            });
        }
        if !x.iter().all(|value| value.is_finite())
            || !v.iter().all(|value| value.is_finite())
            || !potential_energy.is_finite()
        {
            return Err(SciFormError::numerical(format!(
                "OPES became non-finite at step {step}"
            )));
        }
    }

    let free_energy = opes_free_energy(&kernels, config, kbt);
    let reweighted_free_energy = reweight_free_energy(
        &trace,
        &ReweightConfig {
            grid_min: config.grid_min,
            grid_max: config.grid_max,
            grid_bins: config.grid_bins,
            temperature_k: config.temperature_k,
            block_size: 1,
        },
    )?;
    let mut notes = vec![
        "OPES-style kernel-density bias with a Langevin thermostat.".to_string(),
        "Bias-derived and histogram-reweighted surfaces are both reported.".to_string(),
        "This is a bounded-grid OPES approximation, not the full OPES flow equations.".to_string(),
    ];
    if trace.len() < 10 {
        notes.push("Very short run: the free-energy surface is not converged.".to_string());
    }
    Ok(OpesResult {
        cv: config.cv.clone(),
        steps: config.max_steps,
        dt_fs: config.dt_fs,
        cv_trace: trace,
        kernels,
        free_energy,
        reweighted_free_energy,
        final_coords: x,
        final_cv: cv_value,
        notes,
    })
}

fn opes_kernel_density(kernels: &[GaussianDeposit], cv: f64, sigma: f64) -> (f64, f64) {
    if kernels.is_empty() {
        return (0.0, 0.0);
    }
    let mut density = 0.0;
    let mut derivative = 0.0;
    let mut total_weight = 0.0;
    for kernel in kernels {
        let weight = kernel.weight.max(1e-300);
        let delta = cv - kernel.center;
        let value = (-(delta * delta) / (2.0 * sigma * sigma)).exp();
        density += weight * value;
        derivative += weight * -value * delta / (sigma * sigma);
        total_weight += weight;
    }
    let norm = total_weight.max(1e-300);
    (density / norm, derivative / norm)
}

fn opes_bias_and_derivative(
    kernels: &[GaussianDeposit],
    cv: f64,
    config: &OpesConfig,
    kbt: f64,
) -> (f64, f64) {
    if kernels.is_empty() {
        return (0.0, 0.0);
    }
    let (density, derivative) = opes_kernel_density(kernels, cv, config.sigma);
    let max_density = kernels
        .iter()
        .map(|kernel| opes_kernel_density(kernels, kernel.center, config.sigma).0)
        .fold(0.0f64, f64::max)
        .max(config.regularization);
    let ratio = density / max_density + config.regularization;
    let factor = 1.0 - 1.0 / config.bias_factor;
    let bias = factor * kbt * ratio.ln();
    let derivative = factor * kbt * derivative / (density + config.regularization * max_density);
    (bias, derivative)
}

fn opes_free_energy(
    kernels: &[GaussianDeposit],
    config: &OpesConfig,
    kbt: f64,
) -> FreeEnergySurface {
    let mut values = Vec::with_capacity(config.grid_bins);
    for index in 0..config.grid_bins {
        let cv = config.grid_min
            + (config.grid_max - config.grid_min) * index as f64 / (config.grid_bins - 1) as f64;
        let (bias, _) = opes_bias_and_derivative(kernels, cv, config, kbt);
        values.push(-bias);
    }
    let minimum = values.iter().copied().fold(f64::INFINITY, f64::min);
    if minimum.is_finite() {
        for value in &mut values {
            *value -= minimum;
        }
    }
    FreeEnergySurface {
        min: config.grid_min,
        max: config.grid_max,
        bins: config.grid_bins,
        values,
    }
}

/// UFF-backed potential for molecular metadynamics.
pub struct UffPotential {
    force_field: MolecularForceField,
    n_atoms: usize,
}

impl UffPotential {
    /// Build a UFF potential from a SMILES string.
    pub fn from_smiles(smiles: &str) -> Result<Self, SciFormError> {
        let molecule = graph::Molecule::from_smiles(smiles)?;
        let n_atoms = molecule.graph.node_count();
        let force_field = crate::forcefield::builder::build_uff_force_field(&molecule);
        Ok(Self {
            force_field,
            n_atoms,
        })
    }

    /// Number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.n_atoms
    }
}

impl Potential for UffPotential {
    fn energy(&self, coords: &[f64]) -> f64 {
        let mut gradient = vec![0.0; coords.len()];
        self.force_field
            .compute_system_energy_and_gradients(coords, &mut gradient)
    }

    fn gradient(&self, coords: &[f64], gradient: &mut [f64]) {
        self.force_field
            .compute_system_energy_and_gradients(coords, gradient);
    }
}

/// Atomic masses repeated for the three Cartesian components.
pub fn cartesian_masses(smiles: &str) -> Result<Vec<f64>, SciFormError> {
    let molecule = graph::Molecule::from_smiles(smiles)?;
    let mut masses = Vec::with_capacity(molecule.graph.node_count() * 3);
    for atom in molecule.graph.node_weights() {
        let mass = crate::dynamics::atomic_mass_amu(atom.element);
        if mass <= 0.0 {
            return Err(SciFormError::unsupported(format!(
                "no atomic mass for element Z={}",
                atom.element
            )));
        }
        masses.extend_from_slice(&[mass, mass, mass]);
    }
    Ok(masses)
}

/// Run metadynamics on a UFF molecule.
pub fn run_uff_metadynamics(
    smiles: &str,
    initial_coords: &[f64],
    config: &MetaDynamicsConfig,
) -> Result<MetaDynamicsResult, SciFormError> {
    let potential = UffPotential::from_smiles(smiles)?;
    if initial_coords.len() != potential.n_atoms() * 3 {
        return Err(SciFormError::geometry(format!(
            "coords length {} != 3 * atoms {}",
            initial_coords.len(),
            potential.n_atoms()
        )));
    }
    let masses = cartesian_masses(smiles)?;
    run_metadynamics(&potential, initial_coords, &masses, config)
}

/// Run OPES on a UFF molecule.
pub fn run_uff_opes(
    smiles: &str,
    initial_coords: &[f64],
    config: &OpesConfig,
) -> Result<OpesResult, SciFormError> {
    let potential = UffPotential::from_smiles(smiles)?;
    if initial_coords.len() != potential.n_atoms() * 3 {
        return Err(SciFormError::geometry(format!(
            "coords length {} != 3 * atoms {}",
            initial_coords.len(),
            potential.n_atoms()
        )));
    }
    let masses = cartesian_masses(smiles)?;
    run_opes(&potential, initial_coords, &masses, config)
}

/// A one-dimensional double well for tests and demonstrations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DoubleWell {
    /// Barrier height in kcal/mol.
    pub barrier_kcal_mol: f64,
    /// Minimum position in Å.
    pub minimum: f64,
}

impl DoubleWell {
    /// Create a symmetric double well.
    pub fn new(barrier_kcal_mol: f64, minimum: f64) -> Self {
        Self {
            barrier_kcal_mol,
            minimum,
        }
    }

    fn coefficient(&self) -> f64 {
        self.barrier_kcal_mol / self.minimum.powi(4)
    }
}

impl Potential for DoubleWell {
    fn energy(&self, coords: &[f64]) -> f64 {
        let x = coords[0];
        self.coefficient() * (x * x - self.minimum * self.minimum).powi(2)
    }

    fn gradient(&self, coords: &[f64], gradient: &mut [f64]) {
        let x = coords[0];
        gradient[0] = 4.0 * self.coefficient() * x * (x * x - self.minimum * self.minimum);
    }
}

fn deposit_gaussian(
    gaussians: &mut Vec<GaussianDeposit>,
    step: usize,
    center: f64,
    config: &MetaDynamicsConfig,
    current_bias: f64,
) -> Result<(), SciFormError> {
    let delta_t = (config.bias_factor - 1.0) * config.temperature_k;
    if delta_t <= 0.0 {
        return Err(SciFormError::validation(
            "bias_factor must produce a positive Delta T",
        ));
    }
    let height = config.height_kcal_mol * (-current_bias / (R_GAS_KCAL_MOLK * delta_t)).exp();
    if !height.is_finite() || height <= 0.0 {
        return Err(SciFormError::numerical(
            "metadynamics Gaussian height became non-finite",
        ));
    }
    gaussians.push(GaussianDeposit {
        step,
        center,
        height_kcal_mol: height,
        sigma: config.sigma,
        weight: 1.0,
    });
    Ok(())
}

fn bias_and_derivative(gaussians: &[GaussianDeposit], cv: f64, sigma: f64) -> (f64, f64) {
    let mut energy = 0.0;
    let mut derivative = 0.0;
    for gaussian in gaussians {
        let delta = cv - gaussian.center;
        let exponent = -(delta * delta) / (2.0 * sigma * sigma);
        let value = gaussian.height_kcal_mol * exponent.exp();
        energy += value;
        derivative += -value * delta / (sigma * sigma);
    }
    (energy, derivative)
}

fn build_free_energy(
    gaussians: &[GaussianDeposit],
    config: &MetaDynamicsConfig,
) -> FreeEnergySurface {
    let mut values = Vec::with_capacity(config.grid_bins);
    for index in 0..config.grid_bins {
        let cv = config.grid_min
            + (config.grid_max - config.grid_min) * index as f64 / (config.grid_bins - 1) as f64;
        let (bias, _) = bias_and_derivative(gaussians, cv, config.sigma);
        values.push(-bias);
    }
    let minimum = values.iter().copied().fold(f64::INFINITY, f64::min);
    if minimum.is_finite() {
        for value in &mut values {
            *value -= minimum;
        }
    }
    FreeEnergySurface {
        min: config.grid_min,
        max: config.grid_max,
        bins: config.grid_bins,
        values,
    }
}

fn atom(coords: &[f64], index: usize) -> Result<[f64; 3], SciFormError> {
    let start = index * 3;
    if start + 2 >= coords.len() {
        return Err(SciFormError::validation(format!(
            "atom index {index} is out of range"
        )));
    }
    Ok([coords[start], coords[start + 1], coords[start + 2]])
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn angle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let u = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let v = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
    let dot = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    let nu = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
    let nv = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    (dot / (nu * nv).max(1e-12)).clamp(-1.0, 1.0).acos()
}

fn dihedral(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let b1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let b2 = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
    let b3 = [d[0] - c[0], d[1] - c[1], d[2] - c[2]];
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let m1 = cross(n1, normalize(b2));
    let x = dot(n1, n2);
    let y = dot(m1, n2);
    y.atan2(x)
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn normalize(a: [f64; 3]) -> [f64; 3] {
    let norm = dot(a, a).sqrt().max(1e-12);
    [a[0] / norm, a[1] / norm, a[2] / norm]
}

fn sample_standard_normal(rng: &mut StdRng) -> f64 {
    let u1 = (1.0 - rng.gen::<f64>()).max(1e-12);
    let u2 = rng.gen::<f64>();
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}

/// One accepted or rejected replica-exchange proposal.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReplicaSwap {
    /// Integration step.
    pub step: usize,
    /// First replica index.
    pub replica_i: usize,
    /// Second replica index.
    pub replica_j: usize,
    /// Whether the configurations were exchanged.
    pub accepted: bool,
    /// Metropolis exponent `(β_i-β_j)(U_i-U_j)`.
    pub delta: f64,
}

/// One potential-energy sample from a replica.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReplicaPoint {
    /// Integration step.
    pub step: usize,
    /// Potential energy in kcal/mol.
    pub potential_kcal_mol: f64,
}

/// Configuration for Hamiltonian/temperature replica exchange.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReplicaExchangeConfig {
    /// Temperatures in K, one per replica.
    pub temperatures_k: Vec<f64>,
    /// Integration time step in fs.
    pub dt_fs: f64,
    /// Langevin friction in ps^-1.
    pub friction_ps_inv: f64,
    /// Number of integration steps.
    pub max_steps: usize,
    /// Attempt adjacent swaps every `swap_interval` steps.
    pub swap_interval: usize,
    /// RNG seed.
    pub seed: u64,
    /// Store a sample every `trace_stride` steps.
    pub trace_stride: usize,
}

impl Default for ReplicaExchangeConfig {
    fn default() -> Self {
        Self {
            temperatures_k: vec![300.0, 400.0, 500.0, 600.0],
            dt_fs: 0.5,
            friction_ps_inv: 1.0,
            max_steps: 10_000,
            swap_interval: 50,
            seed: 42,
            trace_stride: 10,
        }
    }
}

impl ReplicaExchangeConfig {
    fn validate(&self, coords_len: usize, masses_len: usize) -> Result<(), SciFormError> {
        if coords_len == 0 || masses_len != coords_len {
            return Err(SciFormError::validation(
                "coords and masses must be non-empty and of equal length",
            ));
        }
        if self.temperatures_k.len() < 2 {
            return Err(SciFormError::validation(
                "at least two replica temperatures are required",
            ));
        }
        if self
            .temperatures_k
            .iter()
            .any(|temperature| !temperature.is_finite() || *temperature <= 0.0)
        {
            return Err(SciFormError::validation(
                "replica temperatures must be finite and > 0",
            ));
        }
        if !self.dt_fs.is_finite() || self.dt_fs <= 0.0 {
            return Err(SciFormError::validation("dt_fs must be finite and > 0"));
        }
        if !self.friction_ps_inv.is_finite() || self.friction_ps_inv < 0.0 {
            return Err(SciFormError::validation(
                "friction_ps_inv must be finite and >= 0",
            ));
        }
        if self.max_steps == 0 || self.swap_interval == 0 || self.trace_stride == 0 {
            return Err(SciFormError::validation(
                "max_steps, swap_interval and trace_stride must be >= 1",
            ));
        }
        Ok(())
    }
}

/// Result of a replica-exchange run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReplicaExchangeResult {
    /// Temperatures in K.
    pub temperatures_k: Vec<f64>,
    /// Potential-energy traces per replica.
    pub traces: Vec<Vec<ReplicaPoint>>,
    /// Swap proposals.
    pub swaps: Vec<ReplicaSwap>,
    /// Final coordinates per replica.
    pub final_coords: Vec<Vec<f64>>,
    /// Notes and caveats.
    pub notes: Vec<String>,
}

/// Run temperature replica exchange with a Langevin thermostat.
///
/// Each replica follows velocity Verlet at its own temperature. Every
/// `swap_interval` steps, adjacent replicas attempt a configuration exchange
/// with probability
///
/// ```text
/// P_acc = min(1, exp((β_i - β_j)(U_i - U_j))).
/// ```
pub fn run_replica_exchange<P: Potential>(
    potential: &P,
    initial_coords: &[f64],
    masses_amu: &[f64],
    config: &ReplicaExchangeConfig,
) -> Result<ReplicaExchangeResult, SciFormError> {
    config.validate(initial_coords.len(), masses_amu.len())?;
    let n_replicas = config.temperatures_k.len();
    let n_coords = initial_coords.len();
    let mut rng = StdRng::seed_from_u64(config.seed);
    let mut coords = vec![initial_coords.to_vec(); n_replicas];
    let mut velocities = vec![vec![0.0; n_coords]; n_replicas];
    for (replica, temperature) in config.temperatures_k.iter().enumerate() {
        for index in 0..n_coords {
            let sigma = (R_GAS_KCAL_MOLK * temperature
                / (masses_amu[index] * AMU_ANGFS2_TO_KCAL_MOL))
                .sqrt();
            velocities[replica][index] = sigma * sample_standard_normal(&mut rng);
        }
    }
    let mut traces: Vec<Vec<ReplicaPoint>> = (0..n_replicas)
        .map(|replica| {
            vec![ReplicaPoint {
                step: 0,
                potential_kcal_mol: potential.energy(&coords[replica]),
            }]
        })
        .collect();
    let mut swaps = Vec::new();
    let friction_fs_inv = config.friction_ps_inv * 1e-3;

    for step in 1..=config.max_steps {
        for replica in 0..n_replicas {
            let temperature = config.temperatures_k[replica];
            let kbt = R_GAS_KCAL_MOLK * temperature;
            let mut gradient = vec![0.0; n_coords];
            potential.gradient(&coords[replica], &mut gradient);
            for index in 0..n_coords {
                let acceleration = -gradient[index] / (masses_amu[index] * AMU_ANGFS2_TO_KCAL_MOL);
                velocities[replica][index] += 0.5 * acceleration * config.dt_fs;
                coords[replica][index] += velocities[replica][index] * config.dt_fs;
            }
            potential.gradient(&coords[replica], &mut gradient);
            for index in 0..n_coords {
                let acceleration = -gradient[index] / (masses_amu[index] * AMU_ANGFS2_TO_KCAL_MOL);
                velocities[replica][index] += 0.5 * acceleration * config.dt_fs;
            }
            if friction_fs_inv > 0.0 {
                let decay = (-friction_fs_inv * config.dt_fs).exp();
                let noise_scale = (kbt * (1.0 - decay * decay)).sqrt();
                for index in 0..n_coords {
                    let scale = noise_scale / (masses_amu[index] * AMU_ANGFS2_TO_KCAL_MOL).sqrt();
                    velocities[replica][index] = decay * velocities[replica][index]
                        + scale * sample_standard_normal(&mut rng);
                }
            }
        }

        if step % config.swap_interval == 0 {
            for pair in 0..(n_replicas / 2) {
                let i = 2 * pair;
                let j = i + 1;
                let energy_i = potential.energy(&coords[i]);
                let energy_j = potential.energy(&coords[j]);
                let beta_i = 1.0 / (R_GAS_KCAL_MOLK * config.temperatures_k[i]);
                let beta_j = 1.0 / (R_GAS_KCAL_MOLK * config.temperatures_k[j]);
                let delta = (beta_i - beta_j) * (energy_i - energy_j);
                let accepted = delta >= 0.0 || rng.gen::<f64>() < delta.exp();
                swaps.push(ReplicaSwap {
                    step,
                    replica_i: i,
                    replica_j: j,
                    accepted,
                    delta,
                });
                if accepted {
                    coords.swap(i, j);
                    velocities.swap(i, j);
                }
            }
        }

        if step % config.trace_stride == 0 {
            for replica in 0..n_replicas {
                traces[replica].push(ReplicaPoint {
                    step,
                    potential_kcal_mol: potential.energy(&coords[replica]),
                });
            }
        }
    }

    Ok(ReplicaExchangeResult {
        temperatures_k: config.temperatures_k.clone(),
        traces,
        swaps,
        final_coords: coords,
        notes: vec![
            "Temperature replica exchange with a Langevin thermostat.".to_string(),
            "Adjacent pairs are exchanged with the Metropolis criterion.".to_string(),
            "The output keeps one coordinate set per temperature index.".to_string(),
        ],
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn double_well_config() -> MetaDynamicsConfig {
        MetaDynamicsConfig {
            cv: CollectiveVariable::Coordinate { index: 0 },
            height_kcal_mol: 0.2,
            sigma: 0.15,
            bias_factor: 10.0,
            pace_steps: 50,
            temperature_k: 300.0,
            dt_fs: 0.5,
            friction_ps_inv: 1.0,
            max_steps: 4_000,
            seed: 7,
            grid_min: -1.5,
            grid_max: 1.5,
            grid_bins: 61,
            trace_stride: 10,
        }
    }

    #[test]
    fn cv_distance_and_gradient_are_consistent() {
        let coords = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let cv = CollectiveVariable::Distance { i: 0, j: 1 };
        assert!((cv.value(&coords).unwrap() - 1.0).abs() < 1e-12);
        let gradient = cv.gradient(&coords).unwrap();
        assert!((gradient[0] + 1.0).abs() < 1e-12);
        assert!((gradient[3] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn double_well_metadynamics_samples_both_minima() {
        let potential = DoubleWell::new(2.0, 1.0);
        let result = run_metadynamics(&potential, &[-1.0], &[1.0], &double_well_config()).unwrap();
        let min = result
            .cv_trace
            .iter()
            .map(|point| point.cv)
            .fold(f64::INFINITY, f64::min);
        let max = result
            .cv_trace
            .iter()
            .map(|point| point.cv)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(min < -0.5, "min = {min}");
        assert!(max > 0.5, "max = {max}");
        assert!(!result.gaussians.is_empty());
        assert!(result
            .free_energy
            .values
            .iter()
            .all(|value| value.is_finite()));
    }

    #[test]
    fn uff_gradient_matches_finite_difference() {
        let potential = UffPotential::from_smiles("CCO").unwrap();
        let conformer = crate::embed("CCO", 42);
        assert!(conformer.error.is_none());
        let coords = conformer.coords;
        let mut gradient = vec![0.0; coords.len()];
        potential.gradient(&coords, &mut gradient);
        let step = 1e-5;
        for index in 0..coords.len() {
            let mut plus = coords.clone();
            let mut minus = coords.clone();
            plus[index] += step;
            minus[index] -= step;
            let numerical = (potential.energy(&plus) - potential.energy(&minus)) / (2.0 * step);
            assert!(
                (gradient[index] - numerical).abs() < 1e-4,
                "index {index}: analytic {} numerical {}",
                gradient[index],
                numerical
            );
        }
    }

    #[test]
    fn uff_metadynamics_smoke() {
        let conformer = crate::embed("CCO", 42);
        assert!(conformer.error.is_none());
        let config = MetaDynamicsConfig {
            cv: CollectiveVariable::Distance { i: 0, j: 1 },
            max_steps: 200,
            pace_steps: 25,
            trace_stride: 10,
            grid_min: 1.0,
            grid_max: 2.0,
            grid_bins: 21,
            ..MetaDynamicsConfig::default()
        };
        let result = run_uff_metadynamics("CCO", &conformer.coords, &config).unwrap();
        assert_eq!(result.cv_trace.len(), 21);
        assert!(result.final_coords.iter().all(|value| value.is_finite()));
        assert!(!result.gaussians.is_empty());
    }

    #[test]
    fn replica_exchange_runs_and_records_swaps() {
        let potential = DoubleWell::new(2.0, 1.0);
        let config = ReplicaExchangeConfig {
            temperatures_k: vec![200.0, 300.0, 400.0, 500.0],
            max_steps: 500,
            swap_interval: 20,
            ..ReplicaExchangeConfig::default()
        };
        let result = run_replica_exchange(&potential, &[-1.0], &[1.0], &config).unwrap();
        assert_eq!(result.traces.len(), 4);
        assert!(!result.swaps.is_empty());
        assert!(result
            .final_coords
            .iter()
            .all(|coords| coords.iter().all(|value| value.is_finite())));
    }

    #[test]
    fn reweighting_a_constant_bias_is_flat() {
        let trace = vec![
            MetaDynamicsPoint {
                step: 0,
                cv: 0.0,
                potential_kcal_mol: 0.0,
                bias_kcal_mol: 0.0,
            },
            MetaDynamicsPoint {
                step: 1,
                cv: 1.0,
                potential_kcal_mol: 0.0,
                bias_kcal_mol: 0.0,
            },
            MetaDynamicsPoint {
                step: 2,
                cv: 2.0,
                potential_kcal_mol: 0.0,
                bias_kcal_mol: 0.0,
            },
        ];
        let surface = reweight_free_energy(
            &trace,
            &ReweightConfig {
                grid_min: 0.0,
                grid_max: 2.0,
                grid_bins: 3,
                temperature_k: 300.0,
                block_size: 1,
            },
        )
        .unwrap();
        assert!(surface.values.iter().all(|value| value.abs() < 1e-12));
    }

    #[test]
    fn correlated_reweighting_reports_block_error() {
        let trace = vec![
            MetaDynamicsPoint {
                step: 0,
                cv: 0.0,
                potential_kcal_mol: 0.0,
                bias_kcal_mol: 0.0,
            },
            MetaDynamicsPoint {
                step: 1,
                cv: 0.5,
                potential_kcal_mol: 0.0,
                bias_kcal_mol: 0.0,
            },
            MetaDynamicsPoint {
                step: 2,
                cv: 1.0,
                potential_kcal_mol: 0.0,
                bias_kcal_mol: 0.0,
            },
        ];
        let result = reweight_free_energy_correlated(
            &trace,
            &ReweightConfig {
                grid_min: 0.0,
                grid_max: 1.0,
                grid_bins: 3,
                temperature_k: 300.0,
                block_size: 2,
            },
        )
        .unwrap();
        assert_eq!(result.n_blocks, 2);
        assert!(result
            .standard_error
            .iter()
            .all(|value| value.is_finite() && *value >= 0.0));
    }

    #[test]
    fn opes_double_well_samples_both_minima() {
        let potential = DoubleWell::new(2.0, 1.0);
        let config = OpesConfig {
            cv: CollectiveVariable::Coordinate { index: 0 },
            sigma: 0.15,
            bias_factor: 10.0,
            pace_steps: 50,
            max_steps: 4_000,
            grid_min: -1.5,
            grid_max: 1.5,
            grid_bins: 61,
            trace_stride: 10,
            ..OpesConfig::default()
        };
        let result = run_opes(&potential, &[-1.0], &[1.0], &config).unwrap();
        let min = result
            .cv_trace
            .iter()
            .map(|p| p.cv)
            .fold(f64::INFINITY, f64::min);
        let max = result
            .cv_trace
            .iter()
            .map(|p| p.cv)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(min < -0.5, "min = {min}");
        assert!(max > 0.5, "max = {max}");
        assert!(result.free_energy.values.iter().all(|v| v.is_finite()));
        assert!(result
            .reweighted_free_energy
            .values
            .iter()
            .any(|v| v.is_finite()));
    }

    #[test]
    fn invalid_config_is_rejected() {
        let potential = DoubleWell::new(2.0, 1.0);
        let mut config = double_well_config();
        config.bias_factor = 1.0;
        assert!(run_metadynamics(&potential, &[0.0], &[1.0], &config).is_err());
    }
}
