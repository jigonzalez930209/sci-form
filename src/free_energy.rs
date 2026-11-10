//! Alchemical free-energy estimators: FEP, BAR and MBAR.
//!
//! This module implements the statistical estimator layer, not the alchemical
//! molecular dynamics itself. A caller supplies samples with their reduced
//! potentials `u_k = β U_k` at every alchemical state; the estimator returns
//! the dimensionless free-energy differences and, when a temperature is given,
//! the corresponding kcal/mol values.
//!
//! # Units
//!
//! Reduced potentials and `Δf` are in units of `kBT` and are dimensionless.
//! `ΔG` is in kcal/mol when `temperature_k` is set. Overlap and effective
//! sample size are dimensionless.
//!
//! # Reference
//!
//! Zwanzig, *J. Chem. Phys.* **1954**, 22, 1420 (FEP);
//! Bennett, *J. Comput. Phys.* **1976**, 22, 245 (BAR);
//! Shirts & Chodera, *J. Chem. Phys.* **2008**, 129, 124105 (MBAR).

use serde::{Deserialize, Serialize};

use crate::error::SciFormError;

const R_KCAL_PER_MOLK: f64 = 0.001_987_204_258_640_83;

/// Free-energy estimator.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FreeEnergyEstimator {
    /// One-directional exponential average from state 0.
    Fep,
    /// Two-state Bennett acceptance ratio, implemented through the MBAR
    /// two-state solution.
    Bar,
    /// Multistate Bennett acceptance ratio.
    Mbar,
}

/// One sample and its reduced potentials at every state.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlchemicalSample {
    /// State that generated the sample (0-based).
    pub origin_state: usize,
    /// Reduced potential `βU_k` at every state, in state order.
    pub reduced_potentials: Vec<f64>,
}

/// Configuration for a free-energy estimate.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FreeEnergyConfig {
    /// Number of alchemical states.
    pub n_states: usize,
    /// Temperature in K. When set, `delta_g_kcal_mol` is filled.
    pub temperature_k: Option<f64>,
    /// MBAR convergence tolerance in `kBT`.
    pub tolerance: f64,
    /// Maximum MBAR iterations.
    pub max_iterations: usize,
}

impl Default for FreeEnergyConfig {
    fn default() -> Self {
        Self {
            n_states: 2,
            temperature_k: Some(298.15),
            tolerance: 1e-10,
            max_iterations: 10_000,
        }
    }
}

/// Result of a free-energy estimate.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FreeEnergyResult {
    /// Estimator used.
    pub estimator: FreeEnergyEstimator,
    /// Dimensionless free energies relative to state 0, in `kBT`.
    pub delta_f_kbt: Vec<f64>,
    /// Free energies in kcal/mol when a temperature was supplied.
    pub delta_g_kcal_mol: Option<Vec<f64>>,
    /// Symmetric overlap matrix `O_kl = Σ_n W_nk W_nl`.
    pub overlap: Vec<Vec<f64>>,
    /// Effective sample size per state, `1 / Σ_n W_nk²`.
    pub effective_sample_size: Vec<f64>,
    /// Whether MBAR converged.
    pub converged: bool,
    /// MBAR iterations used.
    pub iterations: usize,
    /// Notes and caveats.
    pub notes: Vec<String>,
}

impl FreeEnergyConfig {
    fn validate(&self) -> Result<(), SciFormError> {
        if self.n_states < 2 {
            return Err(SciFormError::validation("n_states must be >= 2"));
        }
        if !self.tolerance.is_finite() || self.tolerance <= 0.0 {
            return Err(SciFormError::validation("tolerance must be finite and > 0"));
        }
        if self.max_iterations == 0 {
            return Err(SciFormError::validation("max_iterations must be >= 1"));
        }
        if let Some(temperature) = self.temperature_k {
            if !temperature.is_finite() || temperature <= 0.0 {
                return Err(SciFormError::validation(
                    "temperature_k must be finite and > 0",
                ));
            }
        }
        Ok(())
    }
}

/// Estimate a free-energy difference from alchemical samples.
pub fn estimate_free_energy(
    samples: &[AlchemicalSample],
    estimator: FreeEnergyEstimator,
    config: &FreeEnergyConfig,
) -> Result<FreeEnergyResult, SciFormError> {
    config.validate()?;
    if samples.is_empty() {
        return Err(SciFormError::validation(
            "at least one alchemical sample is required",
        ));
    }
    for (index, sample) in samples.iter().enumerate() {
        if sample.origin_state >= config.n_states {
            return Err(SciFormError::validation(format!(
                "sample {index} has origin_state {} outside 0..{}",
                sample.origin_state, config.n_states
            )));
        }
        if sample.reduced_potentials.len() != config.n_states {
            return Err(SciFormError::validation(format!(
                "sample {index} has {} reduced potentials, expected {}",
                sample.reduced_potentials.len(),
                config.n_states
            )));
        }
        if sample
            .reduced_potentials
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(SciFormError::validation(format!(
                "sample {index} has a non-finite reduced potential"
            )));
        }
    }

    let mut counts = vec![0usize; config.n_states];
    for sample in samples {
        counts[sample.origin_state] += 1;
    }
    if counts.contains(&0) {
        return Err(SciFormError::validation(
            "every alchemical state must contribute at least one sample",
        ));
    }

    let delta_f = match estimator {
        FreeEnergyEstimator::Fep => fep(samples, config.n_states)?,
        FreeEnergyEstimator::Bar | FreeEnergyEstimator::Mbar => mbar(samples, &counts, config)?,
    };

    let delta_g = config.temperature_k.map(|temperature| {
        delta_f
            .iter()
            .map(|value| value * R_KCAL_PER_MOLK * temperature)
            .collect()
    });

    let (overlap, effective_sample_size, converged, iterations) = match estimator {
        FreeEnergyEstimator::Fep => (
            vec![vec![1.0; config.n_states]; config.n_states],
            counts.iter().map(|count| *count as f64).collect(),
            true,
            1,
        ),
        FreeEnergyEstimator::Bar | FreeEnergyEstimator::Mbar => {
            let (weights, converged, iterations) = mbar_weights(samples, &counts, config)?;
            let mut overlap = vec![vec![0.0; config.n_states]; config.n_states];
            let mut ess = vec![0.0; config.n_states];
            for k in 0..config.n_states {
                let mut norm = 0.0;
                for row in &weights {
                    norm += row[k] * row[k];
                }
                ess[k] = if norm > 0.0 { 1.0 / norm } else { 0.0 };
                for l in 0..config.n_states {
                    overlap[k][l] = weights.iter().map(|row| row[k] * row[l]).sum();
                }
            }
            (overlap, ess, converged, iterations)
        }
    };

    let mut notes = vec![
        format!("Estimator: {:?}.", estimator),
        "Reduced potentials and Δf are in kBT.".to_string(),
    ];
    if estimator == FreeEnergyEstimator::Fep {
        notes.push(
            "FEP uses only the samples from state 0 and is biased when the forward and reverse distributions do not overlap."
                .to_string(),
        );
    }
    if !converged {
        notes.push("MBAR did not converge within max_iterations.".to_string());
    }

    Ok(FreeEnergyResult {
        estimator,
        delta_f_kbt: delta_f,
        delta_g_kcal_mol: delta_g,
        overlap,
        effective_sample_size,
        converged,
        iterations,
        notes,
    })
}

fn fep(samples: &[AlchemicalSample], n_states: usize) -> Result<Vec<f64>, SciFormError> {
    if n_states != 2 {
        return Err(SciFormError::validation(
            "FEP requires exactly two alchemical states",
        ));
    }
    let work: Vec<f64> = samples
        .iter()
        .filter(|sample| sample.origin_state == 0)
        .map(|sample| sample.reduced_potentials[1] - sample.reduced_potentials[0])
        .collect();
    if work.is_empty() {
        return Err(SciFormError::validation(
            "FEP requires at least one sample from state 0",
        ));
    }
    let negative_work: Vec<f64> = work.iter().map(|value| -*value).collect();
    let delta = -log_sum_exp(&negative_work) + (work.len() as f64).ln();
    Ok(vec![0.0, delta])
}

fn mbar(
    samples: &[AlchemicalSample],
    counts: &[usize],
    config: &FreeEnergyConfig,
) -> Result<Vec<f64>, SciFormError> {
    let (_, _, _) = mbar_weights(samples, counts, config)?;
    let (f, _, _) = solve_mbar(samples, counts, config)?;
    let reference = f[0];
    Ok(f.iter().map(|value| value - reference).collect())
}

fn mbar_weights(
    samples: &[AlchemicalSample],
    counts: &[usize],
    config: &FreeEnergyConfig,
) -> Result<(Vec<Vec<f64>>, bool, usize), SciFormError> {
    let (f, converged, iterations) = solve_mbar(samples, counts, config)?;
    let n_states = config.n_states;
    let mut weights = vec![vec![0.0; n_states]; samples.len()];
    for (sample_index, sample) in samples.iter().enumerate() {
        let ln_denom = log_sum_exp(
            &(0..n_states)
                .map(|k| (counts[k] as f64).ln() + f[k] - sample.reduced_potentials[k])
                .collect::<Vec<_>>(),
        );
        for k in 0..n_states {
            weights[sample_index][k] = (f[k] - sample.reduced_potentials[k] - ln_denom).exp();
        }
    }
    Ok((weights, converged, iterations))
}

fn solve_mbar(
    samples: &[AlchemicalSample],
    counts: &[usize],
    config: &FreeEnergyConfig,
) -> Result<(Vec<f64>, bool, usize), SciFormError> {
    let n_states = config.n_states;
    let mut f = vec![0.0; n_states];
    let mut converged = false;
    let mut iterations = 0;
    for iteration in 1..=config.max_iterations {
        iterations = iteration;
        let mut log_denom = Vec::with_capacity(samples.len());
        for sample in samples {
            let terms: Vec<f64> = (0..n_states)
                .map(|k| (counts[k] as f64).ln() + f[k] - sample.reduced_potentials[k])
                .collect();
            log_denom.push(log_sum_exp(&terms));
        }
        let mut next = vec![0.0; n_states];
        for k in 0..n_states {
            let terms: Vec<f64> = samples
                .iter()
                .enumerate()
                .map(|(index, sample)| -sample.reduced_potentials[k] - log_denom[index])
                .collect();
            next[k] = -log_sum_exp(&terms);
        }
        let reference = next[0];
        for value in &mut next {
            *value -= reference;
        }
        let change = next
            .iter()
            .zip(&f)
            .map(|(new, old)| (new - old).abs())
            .fold(0.0, f64::max);
        f = next;
        if change < config.tolerance {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err(SciFormError::numerical(
            "MBAR did not converge within max_iterations",
        ));
    }
    Ok((f, converged, iterations))
}

fn log_sum_exp(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NEG_INFINITY;
    }
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !max.is_finite() {
        return max;
    }
    max + values
        .iter()
        .map(|value| (value - max).exp())
        .sum::<f64>()
        .ln()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn constant_offset_samples(offset: f64, n_per_state: usize) -> Vec<AlchemicalSample> {
        let mut samples = Vec::new();
        for _ in 0..n_per_state {
            samples.push(AlchemicalSample {
                origin_state: 0,
                reduced_potentials: vec![0.0, offset],
            });
            samples.push(AlchemicalSample {
                origin_state: 1,
                reduced_potentials: vec![0.0, offset],
            });
        }
        samples
    }

    #[test]
    fn fep_recovers_a_constant_offset() {
        let samples = constant_offset_samples(1.5, 10);
        let result = estimate_free_energy(
            &samples,
            FreeEnergyEstimator::Fep,
            &FreeEnergyConfig::default(),
        )
        .unwrap();
        assert!((result.delta_f_kbt[1] - 1.5).abs() < 1e-12);
    }

    #[test]
    fn mbar_recovers_a_constant_offset() {
        let samples = constant_offset_samples(1.5, 10);
        let result = estimate_free_energy(
            &samples,
            FreeEnergyEstimator::Mbar,
            &FreeEnergyConfig::default(),
        )
        .unwrap();
        assert!((result.delta_f_kbt[1] - 1.5).abs() < 1e-10);
        assert!(result.converged);
    }

    #[test]
    fn mbar_recovers_three_state_offsets() {
        let mut samples = Vec::new();
        for _ in 0..10 {
            samples.push(AlchemicalSample {
                origin_state: 0,
                reduced_potentials: vec![0.0, 1.0, 2.0],
            });
            samples.push(AlchemicalSample {
                origin_state: 1,
                reduced_potentials: vec![0.0, 1.0, 2.0],
            });
            samples.push(AlchemicalSample {
                origin_state: 2,
                reduced_potentials: vec![0.0, 1.0, 2.0],
            });
        }
        let config = FreeEnergyConfig {
            n_states: 3,
            ..FreeEnergyConfig::default()
        };
        let result = estimate_free_energy(&samples, FreeEnergyEstimator::Mbar, &config).unwrap();
        assert!((result.delta_f_kbt[1] - 1.0).abs() < 1e-10);
        assert!((result.delta_f_kbt[2] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn invalid_sample_shape_is_rejected() {
        let samples = vec![AlchemicalSample {
            origin_state: 0,
            reduced_potentials: vec![0.0],
        }];
        assert!(estimate_free_energy(
            &samples,
            FreeEnergyEstimator::Mbar,
            &FreeEnergyConfig::default()
        )
        .is_err());
    }
}
