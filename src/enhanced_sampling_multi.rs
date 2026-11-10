//! Two-dimensional OPES and multi-CV reweighting.
//!
//! This module extends the 1D OPES prototype to two collective variables. It
//! uses the same self-consistent importance weights, a 2D Gaussian kernel
//! density, a Langevin integrator and a 2D histogram reweighting surface.
//!
//! # Units
//!
//! Coordinates are Å, masses amu, energies kcal/mol, time fs, temperatures K.
//! CV values are Å for distances/coordinates or radians for angles.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

use crate::enhanced_sampling::{CollectiveVariable, Potential};
use crate::error::SciFormError;

const AMU_ANGFS2_TO_KCAL_MOL: f64 = 2_390.057_361_533_49;
const R_GAS_KCAL_MOLK: f64 = 0.001_987_204_258_640_83;

/// Configuration for 2D OPES.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiOpesConfig {
    /// First collective variable.
    pub cv_x: CollectiveVariable,
    /// Second collective variable.
    pub cv_y: CollectiveVariable,
    /// Kernel width along x.
    pub sigma_x: f64,
    /// Kernel width along y.
    pub sigma_y: f64,
    /// Bias factor `γ > 1`.
    pub bias_factor: f64,
    /// Deposit a kernel every `pace_steps` steps.
    pub pace_steps: usize,
    /// Temperature in K.
    pub temperature_k: f64,
    /// Time step in fs.
    pub dt_fs: f64,
    /// Langevin friction in ps⁻¹.
    pub friction_ps_inv: f64,
    /// Number of integration steps.
    pub max_steps: usize,
    /// RNG seed.
    pub seed: u64,
    /// x grid lower bound.
    pub x_min: f64,
    /// x grid upper bound.
    pub x_max: f64,
    /// Number of x bins.
    pub x_bins: usize,
    /// y grid lower bound.
    pub y_min: f64,
    /// y grid upper bound.
    pub y_max: f64,
    /// Number of y bins.
    pub y_bins: usize,
    /// Trace stride.
    pub trace_stride: usize,
    /// Density regularization.
    pub regularization: f64,
    /// Fixed-point flow iterations.
    pub flow_iterations: usize,
}

impl Default for MultiOpesConfig {
    fn default() -> Self {
        Self {
            cv_x: CollectiveVariable::Coordinate { index: 0 },
            cv_y: CollectiveVariable::Coordinate { index: 1 },
            sigma_x: 0.2,
            sigma_y: 0.2,
            bias_factor: 10.0,
            pace_steps: 100,
            temperature_k: 300.0,
            dt_fs: 0.5,
            friction_ps_inv: 1.0,
            max_steps: 10_000,
            seed: 42,
            x_min: -2.0,
            x_max: 2.0,
            x_bins: 41,
            y_min: -2.0,
            y_max: 2.0,
            y_bins: 41,
            trace_stride: 1,
            regularization: 1e-8,
            flow_iterations: 5,
        }
    }
}

impl MultiOpesConfig {
    fn validate(&self) -> Result<(), SciFormError> {
        if !self.sigma_x.is_finite() || self.sigma_x <= 0.0 {
            return Err(SciFormError::validation("sigma_x must be finite and > 0"));
        }
        if !self.sigma_y.is_finite() || self.sigma_y <= 0.0 {
            return Err(SciFormError::validation("sigma_y must be finite and > 0"));
        }
        if !self.bias_factor.is_finite() || self.bias_factor <= 1.0 {
            return Err(SciFormError::validation(
                "bias_factor must be finite and > 1",
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
        if self.pace_steps == 0 || self.max_steps == 0 || self.trace_stride == 0 {
            return Err(SciFormError::validation(
                "pace_steps, max_steps and trace_stride must be >= 1",
            ));
        }
        if self.x_bins < 2 || self.y_bins < 2 {
            return Err(SciFormError::validation("grid bins must be >= 2"));
        }
        if !self.x_min.is_finite()
            || !self.x_max.is_finite()
            || self.x_max <= self.x_min
            || !self.y_min.is_finite()
            || !self.y_max.is_finite()
            || self.y_max <= self.y_min
        {
            return Err(SciFormError::validation(
                "grid maxima must be finite and greater than minima",
            ));
        }
        if !self.regularization.is_finite() || self.regularization <= 0.0 {
            return Err(SciFormError::validation(
                "regularization must be finite and > 0",
            ));
        }
        if self.flow_iterations == 0 {
            return Err(SciFormError::validation("flow_iterations must be >= 1"));
        }
        Ok(())
    }
}

/// One 2D kernel.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiKernel {
    /// Deposit step.
    pub step: usize,
    /// CV centre.
    pub center: [f64; 2],
    /// Importance weight.
    pub weight: f64,
}

/// One trace point.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiTracePoint {
    /// Step.
    pub step: usize,
    /// CV values.
    pub cv: [f64; 2],
    /// Potential energy.
    pub potential_kcal_mol: f64,
    /// Bias energy.
    pub bias_kcal_mol: f64,
}

/// A 2D free-energy grid.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Surface2D {
    /// x lower bound.
    pub x_min: f64,
    /// x upper bound.
    pub x_max: f64,
    /// y lower bound.
    pub y_min: f64,
    /// y upper bound.
    pub y_max: f64,
    /// x bins.
    pub nx: usize,
    /// y bins.
    pub ny: usize,
    /// Row-major values, shifted so the minimum is zero.
    pub values: Vec<f64>,
}

impl Surface2D {
    /// Value at a grid index.
    pub fn at(&self, ix: usize, iy: usize) -> f64 {
        self.values[iy * self.nx + ix]
    }

    /// x coordinate of a bin.
    pub fn x(&self, ix: usize) -> f64 {
        self.x_min + (self.x_max - self.x_min) * ix as f64 / (self.nx - 1) as f64
    }

    /// y coordinate of a bin.
    pub fn y(&self, iy: usize) -> f64 {
        self.y_min + (self.y_max - self.y_min) * iy as f64 / (self.ny - 1) as f64
    }
}

/// Result of a 2D OPES run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiOpesResult {
    /// First CV.
    pub cv_x: CollectiveVariable,
    /// Second CV.
    pub cv_y: CollectiveVariable,
    /// Trace.
    pub trace: Vec<MultiTracePoint>,
    /// Deposited kernels.
    pub kernels: Vec<MultiKernel>,
    /// Bias-derived surface.
    pub free_energy: Surface2D,
    /// Reweighted surface.
    pub reweighted_free_energy: Surface2D,
    /// Final coordinates.
    pub final_coords: Vec<f64>,
    /// Final CV values.
    pub final_cv: [f64; 2],
    /// Notes.
    pub notes: Vec<String>,
}

/// Configuration for 2D reweighting.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiReweightConfig {
    /// x lower bound.
    pub x_min: f64,
    /// x upper bound.
    pub x_max: f64,
    /// x bins.
    pub nx: usize,
    /// y lower bound.
    pub y_min: f64,
    /// y upper bound.
    pub y_max: f64,
    /// y bins.
    pub ny: usize,
    /// Temperature in K.
    pub temperature_k: f64,
}

/// Run 2D OPES with a self-consistent kernel-weight flow.
pub fn run_multi_opes<P: Potential>(
    potential: &P,
    initial_coords: &[f64],
    masses_amu: &[f64],
    config: &MultiOpesConfig,
) -> Result<MultiOpesResult, SciFormError> {
    config.validate()?;
    let n_coords = initial_coords.len();
    if n_coords == 0 || masses_amu.len() != n_coords {
        return Err(SciFormError::validation(
            "coords and masses must be non-empty and of equal length",
        ));
    }
    let mut x = initial_coords.to_vec();
    let mut rng = StdRng::seed_from_u64(config.seed);
    let mut v = vec![0.0; n_coords];
    for index in 0..n_coords {
        let sigma = (R_GAS_KCAL_MOLK * config.temperature_k
            / (masses_amu[index] * AMU_ANGFS2_TO_KCAL_MOL))
            .sqrt();
        v[index] = sigma * normal(&mut rng);
    }
    let mut gradient = vec![0.0; n_coords];
    let mut kernels: Vec<MultiKernel> = Vec::new();
    let mut trace = Vec::new();
    let kbt = R_GAS_KCAL_MOLK * config.temperature_k;
    let friction_fs_inv = config.friction_ps_inv * 1e-3;
    let mut cv = [config.cv_x.value(&x)?, config.cv_y.value(&x)?];
    let mut potential_energy = potential.energy(&x);
    trace.push(MultiTracePoint {
        step: 0,
        cv,
        potential_kcal_mol: potential_energy,
        bias_kcal_mol: 0.0,
    });

    for step in 1..=config.max_steps {
        potential.gradient(&x, &mut gradient);
        let gx = config.cv_x.gradient(&x)?;
        let gy = config.cv_y.gradient(&x)?;
        let (_, dx, dy) = multi_bias(&kernels, cv, config, kbt);
        let mut force = vec![0.0; n_coords];
        for index in 0..n_coords {
            force[index] = -gradient[index] - dx * gx[index] - dy * gy[index];
        }
        for index in 0..n_coords {
            let acceleration = force[index] / (masses_amu[index] * AMU_ANGFS2_TO_KCAL_MOL);
            v[index] += 0.5 * acceleration * config.dt_fs;
            x[index] += v[index] * config.dt_fs;
        }
        potential.gradient(&x, &mut gradient);
        cv = [config.cv_x.value(&x)?, config.cv_y.value(&x)?];
        let (bias, dx, dy) = multi_bias(&kernels, cv, config, kbt);
        let gx = config.cv_x.gradient(&x)?;
        let gy = config.cv_y.gradient(&x)?;
        for index in 0..n_coords {
            force[index] = -gradient[index] - dx * gx[index] - dy * gy[index];
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
                v[index] = decay * v[index] + scale * normal(&mut rng);
            }
        }
        if step % config.pace_steps == 0 {
            kernels.push(MultiKernel {
                step,
                center: cv,
                weight: 1.0,
            });
            update_multi_weights(&mut kernels, config, kbt);
        }
        if step % config.trace_stride == 0 {
            potential_energy = potential.energy(&x);
            trace.push(MultiTracePoint {
                step,
                cv,
                potential_kcal_mol: potential_energy,
                bias_kcal_mol: bias,
            });
        }
        if !x.iter().all(|value| value.is_finite())
            || !v.iter().all(|value| value.is_finite())
            || !potential_energy.is_finite()
        {
            return Err(SciFormError::numerical(format!(
                "2D OPES became non-finite at step {step}"
            )));
        }
    }

    let free_energy = multi_free_energy(&kernels, config, kbt);
    let reweighted_free_energy = reweight_free_energy_2d(
        &trace,
        &MultiReweightConfig {
            x_min: config.x_min,
            x_max: config.x_max,
            nx: config.x_bins,
            y_min: config.y_min,
            y_max: config.y_max,
            ny: config.y_bins,
            temperature_k: config.temperature_k,
        },
    )?;
    Ok(MultiOpesResult {
        cv_x: config.cv_x.clone(),
        cv_y: config.cv_y.clone(),
        trace,
        kernels,
        free_energy,
        reweighted_free_energy,
        final_coords: x,
        final_cv: cv,
        notes: vec![
            "2D OPES with self-consistent importance weights.".to_string(),
            "Bias-derived and histogram-reweighted surfaces are both reported.".to_string(),
        ],
    })
}

/// Reweight a 2D OPES trace into an unbiased surface.
pub fn reweight_free_energy_2d(
    trace: &[MultiTracePoint],
    config: &MultiReweightConfig,
) -> Result<Surface2D, SciFormError> {
    if trace.is_empty() {
        return Err(SciFormError::validation("cannot reweight an empty trace"));
    }
    if config.nx < 2
        || config.ny < 2
        || config.x_max <= config.x_min
        || config.y_max <= config.y_min
    {
        return Err(SciFormError::validation("invalid 2D reweighting grid"));
    }
    if !config.temperature_k.is_finite() || config.temperature_k <= 0.0 {
        return Err(SciFormError::validation("temperature must be positive"));
    }
    let kbt = R_GAS_KCAL_MOLK * config.temperature_k;
    let dx = (config.x_max - config.x_min) / (config.nx - 1) as f64;
    let dy = (config.y_max - config.y_min) / (config.ny - 1) as f64;
    let mut weights = vec![0.0f64; config.nx * config.ny];
    for point in trace {
        let ix = ((point.cv[0] - config.x_min) / dx)
            .round()
            .clamp(0.0, (config.nx - 1) as f64) as usize;
        let iy = ((point.cv[1] - config.y_min) / dy)
            .round()
            .clamp(0.0, (config.ny - 1) as f64) as usize;
        weights[iy * config.nx + ix] += (point.bias_kcal_mol / kbt).exp();
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
    Ok(Surface2D {
        x_min: config.x_min,
        x_max: config.x_max,
        y_min: config.y_min,
        y_max: config.y_max,
        nx: config.nx,
        ny: config.ny,
        values,
    })
}

fn multi_density(
    kernels: &[MultiKernel],
    cv: [f64; 2],
    config: &MultiOpesConfig,
) -> (f64, f64, f64) {
    if kernels.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let mut density = 0.0;
    let mut dx = 0.0;
    let mut dy = 0.0;
    let mut total = 0.0;
    for kernel in kernels {
        let weight = kernel.weight.max(1e-300);
        let mx = cv[0] - kernel.center[0];
        let my = cv[1] - kernel.center[1];
        let value = (-mx * mx / (2.0 * config.sigma_x * config.sigma_x)
            - my * my / (2.0 * config.sigma_y * config.sigma_y))
            .exp();
        density += weight * value;
        dx += weight * -value * mx / (config.sigma_x * config.sigma_x);
        dy += weight * -value * my / (config.sigma_y * config.sigma_y);
        total += weight;
    }
    let norm = total.max(1e-300);
    (density / norm, dx / norm, dy / norm)
}

fn multi_bias(
    kernels: &[MultiKernel],
    cv: [f64; 2],
    config: &MultiOpesConfig,
    kbt: f64,
) -> (f64, f64, f64) {
    if kernels.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let (density, dx, dy) = multi_density(kernels, cv, config);
    let max_density = kernels
        .iter()
        .map(|kernel| multi_density(kernels, kernel.center, config).0)
        .fold(0.0f64, f64::max)
        .max(config.regularization);
    let ratio = density / max_density + config.regularization;
    let factor = 1.0 - 1.0 / config.bias_factor;
    let denominator = density + config.regularization * max_density;
    (
        factor * kbt * ratio.ln(),
        factor * kbt * dx / denominator,
        factor * kbt * dy / denominator,
    )
}

fn update_multi_weights(kernels: &mut [MultiKernel], config: &MultiOpesConfig, kbt: f64) {
    let mut previous: Vec<f64> = kernels.iter().map(|kernel| kernel.weight).collect();
    for _ in 0..config.flow_iterations {
        let mut next = Vec::with_capacity(kernels.len());
        for kernel in kernels.iter() {
            let (bias, _, _) = multi_bias(kernels, kernel.center, config, kbt);
            next.push((bias / kbt).exp().max(1e-12));
        }
        let maximum = next.iter().copied().fold(0.0f64, f64::max).max(1e-300);
        for value in &mut next {
            *value /= maximum;
        }
        for (kernel, weight) in kernels.iter_mut().zip(&next) {
            kernel.weight = *weight;
        }
        let change = next
            .iter()
            .zip(&previous)
            .map(|(new, old)| (new - old).abs())
            .fold(0.0f64, f64::max);
        previous = next;
        if change < 1e-10 {
            break;
        }
    }
}

fn multi_free_energy(kernels: &[MultiKernel], config: &MultiOpesConfig, kbt: f64) -> Surface2D {
    let mut values = Vec::with_capacity(config.x_bins * config.y_bins);
    for iy in 0..config.y_bins {
        let y =
            config.y_min + (config.y_max - config.y_min) * iy as f64 / (config.y_bins - 1) as f64;
        for ix in 0..config.x_bins {
            let x = config.x_min
                + (config.x_max - config.x_min) * ix as f64 / (config.x_bins - 1) as f64;
            let (bias, _, _) = multi_bias(kernels, [x, y], config, kbt);
            values.push(-bias);
        }
    }
    let minimum = values.iter().copied().fold(f64::INFINITY, f64::min);
    if minimum.is_finite() {
        for value in &mut values {
            *value -= minimum;
        }
    }
    Surface2D {
        x_min: config.x_min,
        x_max: config.x_max,
        y_min: config.y_min,
        y_max: config.y_max,
        nx: config.x_bins,
        ny: config.y_bins,
        values,
    }
}

fn normal(rng: &mut StdRng) -> f64 {
    let u1 = (1.0 - rng.gen::<f64>()).max(1e-12);
    let u2 = rng.gen::<f64>();
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Harmonic2D;
    impl Potential for Harmonic2D {
        fn energy(&self, coords: &[f64]) -> f64 {
            0.5 * (coords[0] * coords[0] + coords[1] * coords[1])
        }
        fn gradient(&self, coords: &[f64], gradient: &mut [f64]) {
            gradient[0] = coords[0];
            gradient[1] = coords[1];
        }
    }

    #[test]
    fn two_dimensional_opes_runs_and_reweights() {
        let config = MultiOpesConfig {
            max_steps: 2_000,
            pace_steps: 50,
            trace_stride: 10,
            x_min: -1.0,
            x_max: 1.0,
            y_min: -1.0,
            y_max: 1.0,
            flow_iterations: 3,
            ..MultiOpesConfig::default()
        };
        let result = run_multi_opes(&Harmonic2D, &[0.0, 0.0], &[1.0, 1.0], &config).unwrap();
        assert!(!result.kernels.is_empty());
        assert!(result
            .trace
            .iter()
            .all(|point| point.cv.iter().all(|v| v.is_finite())));
        assert!(result.free_energy.values.iter().any(|v| v.is_finite()));
        assert!(result
            .reweighted_free_energy
            .values
            .iter()
            .any(|v| v.is_finite()));
    }

    #[test]
    fn constant_bias_2d_reweights_flat() {
        let trace = vec![
            MultiTracePoint {
                step: 0,
                cv: [0.0, 0.0],
                potential_kcal_mol: 0.0,
                bias_kcal_mol: 0.0,
            },
            MultiTracePoint {
                step: 1,
                cv: [1.0, 1.0],
                potential_kcal_mol: 0.0,
                bias_kcal_mol: 0.0,
            },
        ];
        let surface = reweight_free_energy_2d(
            &trace,
            &MultiReweightConfig {
                x_min: 0.0,
                x_max: 1.0,
                nx: 3,
                y_min: 0.0,
                y_max: 1.0,
                ny: 3,
                temperature_k: 300.0,
            },
        )
        .unwrap();
        assert!(surface
            .values
            .iter()
            .filter(|value| value.is_finite())
            .all(|value| value.abs() < 1e-12));
    }

    #[test]
    fn invalid_2d_config_is_rejected() {
        let config = MultiOpesConfig {
            x_min: 1.0,
            x_max: -1.0,
            ..MultiOpesConfig::default()
        };
        assert!(run_multi_opes(&Harmonic2D, &[0.0, 0.0], &[1.0, 1.0], &config).is_err());
    }
}
