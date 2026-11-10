//! Kinetic Monte Carlo for elementary-step networks.
//!
//! This is the stochastic companion to the deterministic microkinetic solver in
//! [`crate::kinetics`]. It consumes the same [`ElementaryStep`] and
//! [`ThermodynamicState`] contracts, evaluates HTST rates, and propagates
//! integer species counts with the Gillespie direct method.
//!
//! # Units
//!
//! Rates are in s⁻¹, time in seconds and species counts are dimensionless
//! integers. The network's stoichiometry is integer-valued; fractional
//! coefficients are rejected rather than silently rounded.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

use crate::error::SciFormError;
use crate::kinetics::{evaluate_htst_rate, ElementaryStep, ThermodynamicState};

/// One species in a kMC network.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmcSpecies {
    /// Stable species identifier.
    pub id: String,
    /// Initial integer population.
    pub initial_count: u64,
}

/// One irreversible elementary step.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmcStep {
    /// HTST step definition.
    pub step: ElementaryStep,
    /// Reactant index and integer stoichiometric coefficient.
    pub reactants: Vec<(usize, u32)>,
    /// Product index and integer stoichiometric coefficient.
    pub products: Vec<(usize, u32)>,
}

/// Configuration for a kMC run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmcConfig {
    /// Species and initial populations.
    pub species: Vec<KmcSpecies>,
    /// Irreversible elementary steps.
    pub steps: Vec<KmcStep>,
    /// Thermodynamic state.
    pub state: ThermodynamicState,
    /// Maximum number of events.
    pub max_events: usize,
    /// Maximum simulated time in seconds.
    pub max_time_s: f64,
    /// RNG seed.
    pub seed: u64,
    /// Store a frame every `frame_stride` accepted events.
    pub frame_stride: usize,
}

impl Default for KmcConfig {
    fn default() -> Self {
        Self {
            species: Vec::new(),
            steps: Vec::new(),
            state: ThermodynamicState::default(),
            max_events: 100_000,
            max_time_s: 1.0,
            seed: 42,
            frame_stride: 100,
        }
    }
}

/// One accepted kMC event.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmcEvent {
    /// Event index (1-based).
    pub event: usize,
    /// Time in seconds.
    pub time_s: f64,
    /// Step index that fired.
    pub step_index: usize,
    /// Propensity of the fired step before the event, s⁻¹.
    pub propensity_s_inv: f64,
    /// Species counts after the event.
    pub counts: Vec<u64>,
}

/// A species-count snapshot.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmcFrame {
    /// Simulated time in seconds.
    pub time_s: f64,
    /// Event index.
    pub event: usize,
    /// Species identifiers.
    pub species: Vec<String>,
    /// Species counts.
    pub counts: Vec<u64>,
}

/// Result of a kMC run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmcResult {
    /// Accepted events.
    pub events: Vec<KmcEvent>,
    /// Sampled frames.
    pub frames: Vec<KmcFrame>,
    /// Final counts in the same order as `species`.
    pub final_counts: Vec<u64>,
    /// Total simulated time in seconds.
    pub final_time_s: f64,
    /// Whether the run stopped on `max_events` or `max_time_s`.
    pub stopped_on_limit: bool,
    /// Notes and caveats.
    pub notes: Vec<String>,
}

impl KmcConfig {
    fn validate(&self) -> Result<(), SciFormError> {
        if self.species.is_empty() {
            return Err(SciFormError::validation("at least one species is required"));
        }
        if self.steps.is_empty() {
            return Err(SciFormError::validation("at least one step is required"));
        }
        if self.max_events == 0 {
            return Err(SciFormError::validation("max_events must be >= 1"));
        }
        if !self.max_time_s.is_finite() || self.max_time_s <= 0.0 {
            return Err(SciFormError::validation(
                "max_time_s must be finite and > 0",
            ));
        }
        if self.frame_stride == 0 {
            return Err(SciFormError::validation("frame_stride must be >= 1"));
        }
        for (step_index, step) in self.steps.iter().enumerate() {
            if step.reactants.is_empty() {
                return Err(SciFormError::validation(format!(
                    "step {step_index} has no reactants"
                )));
            }
            for &(species, coefficient) in step.reactants.iter().chain(&step.products) {
                if species >= self.species.len() {
                    return Err(SciFormError::validation(format!(
                        "step {step_index} refers to unknown species index {species}"
                    )));
                }
                if coefficient == 0 {
                    return Err(SciFormError::validation(format!(
                        "step {step_index} has a zero stoichiometric coefficient"
                    )));
                }
            }
        }
        Ok(())
    }
}

/// Run a Gillespie direct-method kMC simulation.
///
/// The propensity of an irreversible step is
///
/// ```text
/// a_j = k_j · Π_r C(n_r, ν_r)
/// ```
///
/// where `k_j` is the HTST forward rate and `C` is the binomial coefficient for
/// consuming `ν_r` copies of a species. Reversible chemistry is represented by
/// two irreducible steps.
pub fn run_kinetic_monte_carlo(config: &KmcConfig) -> Result<KmcResult, SciFormError> {
    config.validate()?;
    let rates = config
        .steps
        .iter()
        .map(|step| evaluate_htst_rate(&step.step, config.state).map_err(SciFormError::from))
        .collect::<Result<Vec<_>, _>>()?;

    let mut counts: Vec<u64> = config
        .species
        .iter()
        .map(|species| species.initial_count)
        .collect();
    let species_ids: Vec<String> = config
        .species
        .iter()
        .map(|species| species.id.clone())
        .collect();
    let mut rng = StdRng::seed_from_u64(config.seed);
    let mut events = Vec::new();
    let mut frames = vec![KmcFrame {
        time_s: 0.0,
        event: 0,
        species: species_ids.clone(),
        counts: counts.clone(),
    }];
    let mut time = 0.0_f64;
    let mut stopped_on_limit = false;

    for event_index in 1..=config.max_events {
        let propensities: Vec<f64> = config
            .steps
            .iter()
            .zip(&rates)
            .map(|(step, rate)| {
                let mut propensity = rate.forward_rate_s_inv;
                for &(species, coefficient) in &step.reactants {
                    propensity *= binomial(counts[species], coefficient);
                }
                propensity
            })
            .collect();
        let total: f64 = propensities.iter().sum();
        if !total.is_finite() || total <= 0.0 {
            break;
        }
        let tau = -rng.gen::<f64>().max(1e-300).ln() / total;
        if time + tau > config.max_time_s {
            stopped_on_limit = true;
            break;
        }
        time += tau;
        let mut threshold = rng.gen::<f64>() * total;
        let mut selected = propensities.len() - 1;
        for (index, propensity) in propensities.iter().enumerate() {
            if threshold <= *propensity {
                selected = index;
                break;
            }
            threshold -= *propensity;
        }
        let step = &config.steps[selected];
        for &(species, coefficient) in &step.reactants {
            counts[species] = counts[species].saturating_sub(coefficient as u64);
        }
        for &(species, coefficient) in &step.products {
            counts[species] = counts[species].saturating_add(coefficient as u64);
        }
        events.push(KmcEvent {
            event: event_index,
            time_s: time,
            step_index: selected,
            propensity_s_inv: propensities[selected],
            counts: counts.clone(),
        });
        if event_index % config.frame_stride == 0 {
            frames.push(KmcFrame {
                time_s: time,
                event: event_index,
                species: species_ids.clone(),
                counts: counts.clone(),
            });
        }
    }

    if events.len() == config.max_events {
        stopped_on_limit = true;
    }
    let mut notes = vec![
        "Gillespie direct method with irreversible elementary steps.".to_string(),
        "Rates come from the existing HTST evaluator.".to_string(),
        "Reversible chemistry must be encoded as two irreversible steps.".to_string(),
    ];
    if events.is_empty() {
        notes.push("No event fired; check the initial populations and rates.".to_string());
    }

    Ok(KmcResult {
        events,
        frames,
        final_counts: counts,
        final_time_s: time,
        stopped_on_limit,
        notes,
    })
}

fn binomial(n: u64, k: u32) -> f64 {
    if k == 0 {
        return 1.0;
    }
    if n < k as u64 {
        return 0.0;
    }
    let mut value = 1.0;
    for index in 0..k {
        value *= (n - index as u64) as f64 / (index + 1) as f64;
    }
    value
}

#[cfg(test)]
mod tests {
    use super::*;

    fn step(activation_ev: f64) -> ElementaryStep {
        ElementaryStep {
            step_id: "A->B".into(),
            activation_free_energy_ev: activation_ev,
            reaction_free_energy_ev: -0.1,
            prefactor_s_inv: Some(1.0e12),
        }
    }

    fn config() -> KmcConfig {
        KmcConfig {
            species: vec![
                KmcSpecies {
                    id: "A".into(),
                    initial_count: 100,
                },
                KmcSpecies {
                    id: "B".into(),
                    initial_count: 0,
                },
            ],
            steps: vec![KmcStep {
                step: step(0.05),
                reactants: vec![(0, 1)],
                products: vec![(1, 1)],
            }],
            state: ThermodynamicState {
                temperature_k: 500.0,
                pressure_bar: 1.0,
            },
            max_events: 10_000,
            max_time_s: 10.0,
            seed: 7,
            frame_stride: 100,
        }
    }

    #[test]
    fn irreversible_conversion_conserves_mass() {
        let result = run_kinetic_monte_carlo(&config()).unwrap();
        assert!(!result.events.is_empty());
        assert_eq!(result.final_counts[0] + result.final_counts[1], 100);
        assert!(result.final_time_s <= 10.0);
    }

    #[test]
    fn high_barrier_network_does_not_fire_quickly() {
        let mut config = config();
        config.steps[0].step = step(2.0);
        config.max_time_s = 1e-12;
        let result = run_kinetic_monte_carlo(&config).unwrap();
        assert!(result.events.is_empty());
        assert!(result.stopped_on_limit);
    }

    #[test]
    fn binomial_coefficients_are_correct() {
        assert_eq!(binomial(5, 0), 1.0);
        assert_eq!(binomial(5, 1), 5.0);
        assert_eq!(binomial(5, 2), 10.0);
        assert_eq!(binomial(2, 3), 0.0);
    }

    #[test]
    fn invalid_species_index_is_rejected() {
        let mut config = config();
        config.steps[0].reactants = vec![(9, 1)];
        assert!(run_kinetic_monte_carlo(&config).is_err());
    }
}
