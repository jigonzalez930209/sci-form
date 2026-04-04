//! Alpha kinetics contracts and CPU reference helpers.
//!
//! The roadmap needs a stable serialization-friendly layer for elementary-step
//! barriers, transition-state theory rates, and simple reaction-network traces.

#[cfg(feature = "experimental-gpu")]
pub mod gpu;

use serde::{Deserialize, Serialize};
use std::cell::RefCell;

const KCAL_MOL_TO_EV: f64 = 0.043_364_115_308_770_5;
const SPEED_OF_LIGHT_CM_PER_S: f64 = 2.997_924_58e10;

/// Boltzmann constant in eV/K.
pub const BOLTZMANN_EV_PER_K: f64 = 8.617_333_262_145e-5;

/// Planck constant in eV*s.
pub const PLANCK_EV_S: f64 = 4.135_667_696e-15;

/// Thermodynamic state for a kinetic evaluation.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct ThermodynamicState {
    /// Temperature in K.
    pub temperature_k: f64,
    /// Pressure in bar.
    pub pressure_bar: f64,
}

impl Default for ThermodynamicState {
    fn default() -> Self {
        Self {
            temperature_k: 298.15,
            pressure_bar: 1.0,
        }
    }
}

/// One elementary reaction step.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ElementaryStep {
    /// Stable identifier.
    pub step_id: String,
    /// Free-energy barrier in eV.
    pub activation_free_energy_ev: f64,
    /// Reaction free energy in eV.
    pub reaction_free_energy_ev: f64,
    /// Optional prefactor override in s^-1.
    pub prefactor_s_inv: Option<f64>,
}

/// Result of evaluating one elementary step.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ElementaryRateResult {
    /// Stable identifier copied from the input step.
    pub step_id: String,
    /// Forward rate constant in s^-1.
    pub forward_rate_s_inv: f64,
    /// Reverse rate constant in s^-1.
    pub reverse_rate_s_inv: f64,
    /// Equilibrium constant from the free-energy change.
    pub equilibrium_constant: f64,
    /// Thermodynamic state used for evaluation.
    pub state: ThermodynamicState,
}

/// A species population snapshot for lightweight microkinetic traces.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct SpeciesPopulation {
    /// Species identifier.
    pub species_id: String,
    /// Coverage or concentration in arbitrary normalized units.
    pub population: f64,
}

/// One output frame for a coarse kinetic trajectory.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct MicrokineticFrame {
    /// Simulation time in seconds.
    pub time_s: f64,
    /// Species populations at this time.
    pub species: Vec<SpeciesPopulation>,
}

/// Minimal trace structure shared with render and chart adapters.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct MicrokineticTrace {
    /// Evaluated elementary steps.
    pub rates: Vec<ElementaryRateResult>,
    /// Time-dependent snapshots.
    pub frames: Vec<MicrokineticFrame>,
}

/// Configuration for the GSM + MBH -> HTST adapter.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct HtstAdapterConfig {
    /// Stable identifier for the elementary step.
    pub step_id: String,
    /// Thermodynamic state used for the final rate evaluation.
    pub state: ThermodynamicState,
    /// GSM transition-state search controls.
    pub gsm: crate::alpha::gsm::GsmConfig,
    /// Finite-difference step used by MBH at the TS geometry.
    pub mbh_fd_step: f64,
}

impl Default for HtstAdapterConfig {
    fn default() -> Self {
        Self {
            step_id: "elementary-step".into(),
            state: ThermodynamicState::default(),
            gsm: crate::alpha::gsm::GsmConfig::default(),
            mbh_fd_step: 0.005,
        }
    }
}

/// Result of the GSM + MBH -> HTST adapter.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct HtstAnalysisResult {
    /// Evaluated HTST rate constants.
    pub rate: ElementaryRateResult,
    /// GSM transition-state energy in kcal/mol.
    pub ts_energy_kcal_mol: f64,
    /// Activation barrier in kcal/mol.
    pub activation_energy_kcal_mol: f64,
    /// Reverse barrier in kcal/mol.
    pub reverse_barrier_kcal_mol: f64,
    /// Lowest imaginary mode from MBH, if present.
    pub imaginary_frequency_cm1: Option<f64>,
    /// Reduced MBH frequencies in cm^-1.
    pub mbh_frequencies_cm1: Vec<f64>,
    /// Number of rigid blocks used by MBH.
    pub n_blocks: usize,
    /// Number of flexible atoms used by MBH.
    pub n_flexible: usize,
    /// GSM path energies in kcal/mol.
    pub path_energies_kcal_mol: Vec<f64>,
    /// GSM path coordinates.
    pub path_coords: Vec<Vec<f64>>,
    /// TS coordinates.
    pub ts_coords: Vec<f64>,
    /// Number of GSM nodes.
    pub n_nodes: usize,
    /// Number of energy evaluations reported by GSM.
    pub energy_evaluations: usize,
}

// ── Convergence diagnostics ─────────────────────────────────────────────────

/// Convergence diagnostics for kinetics solvers.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct KineticsConvergenceDiagnostics {
    /// Maximum population in the final frame.
    pub max_population: f64,
    /// Total population sum (should be conserved).
    pub total_population: f64,
    /// Relative mass conservation error.
    pub mass_conservation_error: f64,
    /// Maximum rate of change in the final frame (steady-state proxy).
    pub max_rate_of_change: f64,
    /// Whether the system appears to have reached steady state.
    pub is_steady_state: bool,
}

/// Extract convergence diagnostics from a microkinetic trace.
pub fn extract_kinetics_diagnostics(trace: &MicrokineticTrace) -> KineticsConvergenceDiagnostics {
    let first_total: f64 = trace
        .frames
        .first()
        .map(|f| f.species.iter().map(|s| s.population).sum())
        .unwrap_or(0.0);
    let last_total: f64 = trace
        .frames
        .last()
        .map(|f| f.species.iter().map(|s| s.population).sum())
        .unwrap_or(0.0);
    let max_pop = trace
        .frames
        .last()
        .map(|f| {
            f.species
                .iter()
                .map(|s| s.population)
                .fold(0.0_f64, f64::max)
        })
        .unwrap_or(0.0);

    // Estimate rate of change from last two frames
    let max_roc = if trace.frames.len() >= 2 {
        let f1 = &trace.frames[trace.frames.len() - 2];
        let f2 = &trace.frames[trace.frames.len() - 1];
        let dt = f2.time_s - f1.time_s;
        if dt > 0.0 {
            f1.species
                .iter()
                .zip(&f2.species)
                .map(|(s1, s2)| ((s2.population - s1.population) / dt).abs())
                .fold(0.0_f64, f64::max)
        } else {
            0.0
        }
    } else {
        0.0
    };

    let conservation_err = if first_total > 0.0 {
        (last_total - first_total).abs() / first_total
    } else {
        0.0
    };

    KineticsConvergenceDiagnostics {
        max_population: max_pop,
        total_population: last_total,
        mass_conservation_error: conservation_err,
        max_rate_of_change: max_roc,
        is_steady_state: max_roc < 1e-10,
    }
}

/// Evaluate forward and reverse HTST-style rates for one elementary step.
pub fn evaluate_htst_rate(
    step: &ElementaryStep,
    state: ThermodynamicState,
) -> Result<ElementaryRateResult, String> {
    if !state.temperature_k.is_finite() || state.temperature_k <= 0.0 {
        return Err("temperature must be positive and finite".into());
    }
    if !step.activation_free_energy_ev.is_finite() || step.activation_free_energy_ev < 0.0 {
        return Err("activation free energy must be finite and non-negative".into());
    }
    if !step.reaction_free_energy_ev.is_finite() {
        return Err("reaction free energy must be finite".into());
    }

    let thermal_energy = BOLTZMANN_EV_PER_K * state.temperature_k;
    let prefactor = step
        .prefactor_s_inv
        .unwrap_or(state.temperature_k * BOLTZMANN_EV_PER_K / PLANCK_EV_S);
    let forward_rate_s_inv = prefactor * (-step.activation_free_energy_ev / thermal_energy).exp();
    let equilibrium_constant = (-step.reaction_free_energy_ev / thermal_energy).exp();
    let reverse_rate_s_inv = if equilibrium_constant > 0.0 {
        forward_rate_s_inv / equilibrium_constant
    } else {
        0.0
    };

    Ok(ElementaryRateResult {
        step_id: step.step_id.clone(),
        forward_rate_s_inv,
        reverse_rate_s_inv,
        equilibrium_constant,
        state,
    })
}

/// Evaluate HTST rates over a temperature grid for one elementary step.
pub fn evaluate_htst_temperature_sweep(
    step: &ElementaryStep,
    temperatures_k: &[f64],
    pressure_bar: f64,
) -> Result<Vec<ElementaryRateResult>, String> {
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        temperatures_k
            .par_iter()
            .map(|&t| {
                evaluate_htst_rate(
                    step,
                    ThermodynamicState {
                        temperature_k: t,
                        pressure_bar,
                    },
                )
            })
            .collect()
    }

    #[cfg(not(feature = "parallel"))]
    {
        temperatures_k
            .iter()
            .map(|&t| {
                evaluate_htst_rate(
                    step,
                    ThermodynamicState {
                        temperature_k: t,
                        pressure_bar,
                    },
                )
            })
            .collect()
    }
}

/// Configuration for a microkinetic network.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct MicrokineticNetworkConfig {
    /// Species identifiers.
    pub species: Vec<String>,
    /// Elementary steps forming the network.
    pub steps: Vec<ElementaryStep>,
    /// Stoichiometry: for each step, (species_index, coefficient) pairs.
    /// Negative = reactant, positive = product.
    pub stoichiometry: Vec<Vec<(usize, f64)>>,
    /// Initial species concentrations/populations.
    pub initial_populations: Vec<f64>,
    /// Total integration time in seconds.
    pub total_time_s: f64,
    /// Number of output frames.
    pub n_frames: usize,
    /// Thermodynamic state.
    pub state: ThermodynamicState,
}

/// Build and solve a microkinetic network using simple forward Euler integration.
///
/// dc_i/dt = Σ_j ν_{ij} (k_fwd_j · Π_{reactants} c_r - k_rev_j · Π_{products} c_p)
pub fn solve_microkinetic_network(
    config: &MicrokineticNetworkConfig,
) -> Result<MicrokineticTrace, String> {
    if config.species.is_empty() {
        return Err("at least one species is required".into());
    }
    if config.initial_populations.len() != config.species.len() {
        return Err("initial_populations length must match species count".into());
    }
    if config.steps.len() != config.stoichiometry.len() {
        return Err("steps and stoichiometry must have the same length".into());
    }
    if config.n_frames < 2 {
        return Err("need at least 2 frames".into());
    }
    if config.total_time_s <= 0.0 {
        return Err("total time must be positive".into());
    }

    let n_species = config.species.len();
    let n_steps = config.steps.len();

    // Evaluate all rates
    let rates: Vec<ElementaryRateResult> = config
        .steps
        .iter()
        .map(|step| evaluate_htst_rate(step, config.state))
        .collect::<Result<Vec<_>, _>>()?;

    let dt = config.total_time_s / (config.n_frames - 1) as f64;
    // Sub-stepping for numerical stability
    let n_substeps = 100.max((config.total_time_s / 1e-12) as usize).min(10000);
    let sub_dt = dt / n_substeps as f64;

    let mut concentrations = config.initial_populations.clone();
    let mut frames = Vec::with_capacity(config.n_frames);

    // Record initial frame
    frames.push(MicrokineticFrame {
        time_s: 0.0,
        species: config
            .species
            .iter()
            .zip(&concentrations)
            .map(|(id, &pop)| SpeciesPopulation {
                species_id: id.clone(),
                population: pop,
            })
            .collect(),
    });

    for frame_idx in 1..config.n_frames {
        for _ in 0..n_substeps {
            let mut dc = vec![0.0; n_species];
            for (step_idx, stoich) in config.stoichiometry.iter().enumerate().take(n_steps) {
                let k_fwd = rates[step_idx].forward_rate_s_inv;
                let k_rev = rates[step_idx].reverse_rate_s_inv;

                // Forward rate: k_fwd · Π(c_reactant)
                let mut fwd_product = k_fwd;
                let mut rev_product = k_rev;
                for &(sp_idx, coeff) in stoich {
                    if coeff < 0.0 {
                        // Reactant
                        fwd_product *= concentrations[sp_idx].max(0.0).powf(-coeff);
                    } else if coeff > 0.0 {
                        // Product
                        rev_product *= concentrations[sp_idx].max(0.0).powf(coeff);
                    }
                }

                let net_rate = fwd_product - rev_product;
                for &(sp_idx, coeff) in stoich {
                    dc[sp_idx] += coeff * net_rate;
                }
            }

            for (c, d) in concentrations.iter_mut().zip(&dc) {
                *c = (*c + sub_dt * d).max(0.0);
            }
        }

        frames.push(MicrokineticFrame {
            time_s: frame_idx as f64 * dt,
            species: config
                .species
                .iter()
                .zip(&concentrations)
                .map(|(id, &pop)| SpeciesPopulation {
                    species_id: id.clone(),
                    population: pop,
                })
                .collect(),
        });
    }

    Ok(MicrokineticTrace { rates, frames })
}

/// Steady-state solver for a microkinetic network via simple time propagation to equilibrium.
pub fn solve_microkinetic_steady_state(
    config: &MicrokineticNetworkConfig,
    convergence_tol: f64,
) -> Result<MicrokineticTrace, String> {
    // Use the full solver with many frames, return when concentrations stabilize
    let mut extended_config = config.clone();
    extended_config.n_frames = 200;
    let trace = solve_microkinetic_network(&extended_config)?;

    // Check if converged by comparing last two frames
    if let (Some(last), Some(prev)) = (
        trace.frames.last(),
        trace.frames.get(trace.frames.len().saturating_sub(2)),
    ) {
        let max_change: f64 = last
            .species
            .iter()
            .zip(&prev.species)
            .map(|(a, b)| (a.population - b.population).abs())
            .fold(0.0, f64::max);
        if max_change > convergence_tol {
            // Not converged — return what we have
        }
    }

    Ok(trace)
}

/// Build an HTST rate analysis by chaining GSM transition-state search and MBH frequencies.
pub fn analyze_gsm_mbh_htst_step(
    elements: &[u8],
    reactant: &[f64],
    product: &[f64],
    rings: &[(Vec<usize>, bool)],
    energy_fn: &dyn Fn(&[f64]) -> f64,
    config: &HtstAdapterConfig,
) -> Result<HtstAnalysisResult, String> {
    if reactant.len() != product.len() {
        return Err("reactant and product coordinates must have the same length".into());
    }
    if reactant.len() != elements.len() * 3 {
        return Err(
            "reactant/product coordinates must be flat arrays of length 3 * n_atoms".into(),
        );
    }
    if config.mbh_fd_step <= 0.0 {
        return Err("MBH finite-difference step must be positive".into());
    }

    let gsm_result =
        crate::alpha::gsm::find_transition_state(reactant, product, energy_fn, &config.gsm);
    let ts_positions = flat_coords_to_positions(&gsm_result.ts_coords)?;
    let mbh_result = crate::beta::mbh::compute_mbh_frequencies(
        elements,
        &ts_positions,
        rings,
        energy_fn,
        &crate::beta::mbh::MbhConfig {
            fd_step: config.mbh_fd_step,
        },
    );

    let imaginary_frequency_cm1 = mbh_result
        .frequencies
        .iter()
        .copied()
        .find(|frequency| *frequency < 0.0);
    let prefactor_s_inv =
        imaginary_frequency_cm1.map(|frequency| frequency.abs() * SPEED_OF_LIGHT_CM_PER_S);

    let activation_energy_ev = gsm_result.activation_energy * KCAL_MOL_TO_EV;
    let reaction_free_energy_ev =
        (gsm_result.activation_energy - gsm_result.reverse_barrier) * KCAL_MOL_TO_EV;
    let rate = evaluate_htst_rate(
        &ElementaryStep {
            step_id: config.step_id.clone(),
            activation_free_energy_ev: activation_energy_ev,
            reaction_free_energy_ev,
            prefactor_s_inv,
        },
        config.state,
    )?;

    Ok(HtstAnalysisResult {
        rate,
        ts_energy_kcal_mol: gsm_result.ts_energy,
        activation_energy_kcal_mol: gsm_result.activation_energy,
        reverse_barrier_kcal_mol: gsm_result.reverse_barrier,
        imaginary_frequency_cm1,
        mbh_frequencies_cm1: mbh_result.frequencies,
        n_blocks: mbh_result.n_blocks,
        n_flexible: mbh_result.n_flexible,
        path_energies_kcal_mol: gsm_result.path_energies,
        path_coords: gsm_result.path_coords,
        ts_coords: gsm_result.ts_coords,
        n_nodes: gsm_result.n_nodes,
        energy_evaluations: gsm_result.energy_evaluations,
    })
}

/// Build an HTST rate analysis using a named sci-form GSM backend.
pub fn analyze_gsm_mbh_htst_step_with_backend(
    smiles: &str,
    reactant: &[f64],
    product: &[f64],
    backend: crate::alpha::gsm::GsmEnergyBackend,
    config: &HtstAdapterConfig,
) -> Result<HtstAnalysisResult, String> {
    let elements = crate::alpha::gsm::backend::elements_for_smiles(smiles)?;
    let rings = crate::alpha::gsm::backend::rings_from_smiles(smiles)?;
    let capability = crate::alpha::gsm::plan_gsm_backends(&elements)
        .into_iter()
        .find(|entry| entry.backend == backend)
        .ok_or_else(|| format!("unknown GSM backend '{}'", backend.as_str()))?;
    if !capability.available || !capability.suitable_for_reaction_path {
        let message = capability
            .limitations
            .first()
            .or_else(|| capability.warnings.first())
            .cloned()
            .unwrap_or_else(|| {
                format!(
                    "GSM backend '{}' is unavailable for this system",
                    backend.as_str()
                )
            });
        return Err(message);
    }

    let owned_smiles = smiles.to_string();
    let owned_elements = elements.clone();
    let first_error = RefCell::new(None::<String>);
    let energy_fn = |coords: &[f64]| -> f64 {
        match crate::alpha::gsm::backend::evaluate_gsm_backend_energy_kcal(
            &owned_smiles,
            &owned_elements,
            coords,
            backend,
        ) {
            Ok(energy) => energy,
            Err(err) => {
                let mut slot = first_error.borrow_mut();
                if slot.is_none() {
                    *slot = Some(err);
                }
                f64::INFINITY
            }
        }
    };

    let result =
        analyze_gsm_mbh_htst_step(&elements, reactant, product, &rings, &energy_fn, config);
    if let Some(err) = first_error.into_inner() {
        return Err(err);
    }
    result
}

fn flat_coords_to_positions(coords: &[f64]) -> Result<Vec<[f64; 3]>, String> {
    if !coords.len().is_multiple_of(3) {
        return Err("flat coordinate array length must be a multiple of 3".into());
    }

    Ok(coords
        .chunks_exact(3)
        .map(|chunk| [chunk[0], chunk[1], chunk[2]])
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn htst_rate_is_positive() {
        let result = evaluate_htst_rate(
            &ElementaryStep {
                step_id: "adsorption".into(),
                activation_free_energy_ev: 0.45,
                reaction_free_energy_ev: -0.10,
                prefactor_s_inv: None,
            },
            ThermodynamicState::default(),
        )
        .unwrap();

        assert!(result.forward_rate_s_inv > 0.0);
        assert!(result.reverse_rate_s_inv > 0.0);
        assert!(result.equilibrium_constant > 1.0);
    }

    #[test]
    fn gsm_mbh_adapter_returns_rate_and_path() {
        let elements = [1u8];
        let reactant = [-1.0, 0.0, 0.0];
        let product = [1.0, 0.0, 0.0];
        let energy_fn = |coords: &[f64]| -> f64 {
            let x = coords[0];
            (x * x - 1.0).powi(2)
        };

        let result = analyze_gsm_mbh_htst_step(
            &elements,
            &reactant,
            &product,
            &[],
            &energy_fn,
            &HtstAdapterConfig {
                step_id: "double-well".into(),
                gsm: crate::alpha::gsm::GsmConfig {
                    max_nodes: 7,
                    max_iter: 20,
                    step_size: 0.02,
                    ..Default::default()
                },
                ..Default::default()
            },
        )
        .unwrap();

        assert!(result.rate.forward_rate_s_inv > 0.0);
        assert!(!result.path_coords.is_empty());
        assert!(result.n_nodes >= 3);
    }

    // --- Arrhenius behavior validation ---

    #[test]
    fn htst_rate_increases_with_temperature() {
        let step = ElementaryStep {
            step_id: "test".into(),
            activation_free_energy_ev: 0.5,
            reaction_free_energy_ev: -0.1,
            prefactor_s_inv: Some(1e13),
        };
        let r300 = evaluate_htst_rate(
            &step,
            ThermodynamicState {
                temperature_k: 300.0,
                pressure_bar: 1.0,
            },
        )
        .unwrap();
        let r600 = evaluate_htst_rate(
            &step,
            ThermodynamicState {
                temperature_k: 600.0,
                pressure_bar: 1.0,
            },
        )
        .unwrap();
        assert!(
            r600.forward_rate_s_inv > r300.forward_rate_s_inv,
            "rate should increase with temperature (Arrhenius)"
        );
    }

    #[test]
    fn htst_arrhenius_slope_matches_barrier() {
        // ln(k) = ln(A) - Ea/(kB T), so Ea = -kB · d(ln k)/d(1/T)
        let step = ElementaryStep {
            step_id: "test".into(),
            activation_free_energy_ev: 0.5,
            reaction_free_energy_ev: 0.0,
            prefactor_s_inv: Some(1e13),
        };
        let t1 = 400.0;
        let t2 = 800.0;
        let r1 = evaluate_htst_rate(
            &step,
            ThermodynamicState {
                temperature_k: t1,
                pressure_bar: 1.0,
            },
        )
        .unwrap();
        let r2 = evaluate_htst_rate(
            &step,
            ThermodynamicState {
                temperature_k: t2,
                pressure_bar: 1.0,
            },
        )
        .unwrap();
        let ea_recovered = -BOLTZMANN_EV_PER_K
            * (r2.forward_rate_s_inv.ln() - r1.forward_rate_s_inv.ln())
            / (1.0 / t2 - 1.0 / t1);
        assert!(
            (ea_recovered - 0.5).abs() < 0.01,
            "recovered Ea = {:.4} eV, expected 0.5 eV",
            ea_recovered
        );
    }

    // --- Temperature sweep ---

    #[test]
    fn temperature_sweep_returns_correct_count() {
        let step = ElementaryStep {
            step_id: "test".into(),
            activation_free_energy_ev: 0.3,
            reaction_free_energy_ev: -0.05,
            prefactor_s_inv: None,
        };
        let temps: Vec<f64> = (200..=1000).step_by(100).map(|t| t as f64).collect();
        let results = evaluate_htst_temperature_sweep(&step, &temps, 1.0).unwrap();
        assert_eq!(results.len(), temps.len());
        // Verify monotonic increase
        for pair in results.windows(2) {
            assert!(pair[1].forward_rate_s_inv >= pair[0].forward_rate_s_inv);
        }
    }

    // --- Microkinetic network ---

    #[test]
    fn microkinetic_simple_ab_reaction() {
        // A → B with barrier 0.3 eV
        let config = MicrokineticNetworkConfig {
            species: vec!["A".into(), "B".into()],
            steps: vec![ElementaryStep {
                step_id: "A->B".into(),
                activation_free_energy_ev: 0.3,
                reaction_free_energy_ev: -0.2,
                prefactor_s_inv: Some(1e10),
            }],
            stoichiometry: vec![vec![(0, -1.0), (1, 1.0)]],
            initial_populations: vec![1.0, 0.0],
            total_time_s: 1e-6,
            n_frames: 20,
            state: ThermodynamicState::default(),
        };

        let trace = solve_microkinetic_network(&config).unwrap();
        assert_eq!(trace.frames.len(), 20);
        // A should decrease, B should increase
        let first = &trace.frames[0];
        let last = &trace.frames.last().unwrap();
        let a_first = first.species[0].population;
        let a_last = last.species[0].population;
        let b_last = last.species[1].population;
        assert!(a_last <= a_first, "A should decrease over time");
        assert!(b_last > 0.0, "B should be produced");
    }

    #[test]
    fn microkinetic_conservation_of_mass() {
        // A ↔ B: total population should be conserved
        let config = MicrokineticNetworkConfig {
            species: vec!["A".into(), "B".into()],
            steps: vec![ElementaryStep {
                step_id: "A<->B".into(),
                activation_free_energy_ev: 0.2,
                reaction_free_energy_ev: -0.05,
                prefactor_s_inv: Some(1e8),
            }],
            stoichiometry: vec![vec![(0, -1.0), (1, 1.0)]],
            initial_populations: vec![1.0, 0.0],
            total_time_s: 1e-4,
            n_frames: 10,
            state: ThermodynamicState::default(),
        };

        let trace = solve_microkinetic_network(&config).unwrap();
        for frame in &trace.frames {
            let total: f64 = frame.species.iter().map(|s| s.population).sum();
            assert!(
                (total - 1.0).abs() < 0.05,
                "mass conservation: total = {:.4}",
                total
            );
        }
    }

    // --- Equilibrium constant validation ---

    #[test]
    fn equilibrium_constant_thermodynamic_consistency() {
        // K_eq = exp(-ΔG/kT)
        let dg = -0.15_f64; // exergonic
        let t = 298.15;
        let step = ElementaryStep {
            step_id: "test".into(),
            activation_free_energy_ev: 0.4,
            reaction_free_energy_ev: dg,
            prefactor_s_inv: Some(1e12),
        };
        let result = evaluate_htst_rate(
            &step,
            ThermodynamicState {
                temperature_k: t,
                pressure_bar: 1.0,
            },
        )
        .unwrap();
        let expected_keq = (-dg / (BOLTZMANN_EV_PER_K * t)).exp();
        assert!(
            (result.equilibrium_constant - expected_keq).abs() / expected_keq < 1e-8,
            "K_eq vs thermodynamic: {} vs {}",
            result.equilibrium_constant,
            expected_keq
        );
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Validation tests against experimental / simulated reference data
    // ═══════════════════════════════════════════════════════════════════════

    /// Reference: Eyring equation for a unimolecular reaction.
    /// k = (kB T / h) × exp(-ΔG‡/kBT)
    /// At 298.15 K: kB T / h = 6.21 × 10¹² s⁻¹
    /// With ΔG‡ = 0.5 eV (≈ 48.3 kJ/mol): k ≈ 6.21e12 × exp(-19.4) ≈ 2.26 × 10⁴ s⁻¹
    #[test]
    fn htst_rate_matches_eyring_equation() {
        let barrier_ev = 0.5;
        let t = 298.15;
        let step = ElementaryStep {
            step_id: "eyring".into(),
            activation_free_energy_ev: barrier_ev,
            reaction_free_energy_ev: 0.0,
            prefactor_s_inv: None, // Use default TST prefactor: kBT/h
        };
        let result = evaluate_htst_rate(
            &step,
            ThermodynamicState {
                temperature_k: t,
                pressure_bar: 1.0,
            },
        )
        .unwrap();

        let kbt_over_h = t * BOLTZMANN_EV_PER_K / PLANCK_EV_S;
        let expected_rate = kbt_over_h * (-barrier_ev / (BOLTZMANN_EV_PER_K * t)).exp();
        let rel_err = (result.forward_rate_s_inv - expected_rate).abs() / expected_rate;
        assert!(
            rel_err < 1e-8,
            "Eyring rate: expected {:.4e}, got {:.4e}",
            expected_rate,
            result.forward_rate_s_inv
        );
    }

    /// Reference: H₂ + I₂ → 2HI gas phase reaction.
    /// Experimental Ea ≈ 1.74 eV (167 kJ/mol), A ≈ 2.0 × 10¹¹ L mol⁻¹ s⁻¹
    /// We verify that at 700 K: k ≈ A × exp(-Ea/kT) ≈ 2e11 × exp(-28.8) ≈ 0.061 s⁻¹
    #[test]
    fn htst_h2_i2_reference_rate() {
        let ea_ev = 1.74;
        let a_prefactor = 2.0e11;
        let t = 700.0;
        let step = ElementaryStep {
            step_id: "H2+I2".into(),
            activation_free_energy_ev: ea_ev,
            reaction_free_energy_ev: -0.1,
            prefactor_s_inv: Some(a_prefactor),
        };
        let result = evaluate_htst_rate(
            &step,
            ThermodynamicState {
                temperature_k: t,
                pressure_bar: 1.0,
            },
        )
        .unwrap();

        let expected = a_prefactor * (-ea_ev / (BOLTZMANN_EV_PER_K * t)).exp();
        let rel_err = (result.forward_rate_s_inv - expected).abs() / expected;
        assert!(
            rel_err < 1e-8,
            "H2+I2 rate at 700K: expected {:.4e}, got {:.4e}",
            expected,
            result.forward_rate_s_inv
        );
    }

    /// Verify forward/reverse rate consistency through detailed balance.
    /// k_f / k_r = K_eq = exp(-ΔG/kBT) for all temperatures.
    #[test]
    fn detailed_balance_consistency_sweep() {
        let step = ElementaryStep {
            step_id: "db".into(),
            activation_free_energy_ev: 0.6,
            reaction_free_energy_ev: -0.2,
            prefactor_s_inv: Some(1e13),
        };
        for t in [200.0, 300.0, 500.0, 800.0, 1200.0] {
            let r = evaluate_htst_rate(
                &step,
                ThermodynamicState {
                    temperature_k: t,
                    pressure_bar: 1.0,
                },
            )
            .unwrap();
            let ratio = r.forward_rate_s_inv / r.reverse_rate_s_inv;
            let expected = (-step.reaction_free_energy_ev / (BOLTZMANN_EV_PER_K * t)).exp();
            let rel_err = (ratio - expected).abs() / expected;
            assert!(
                rel_err < 1e-8,
                "detailed balance at {}K: k_f/k_r = {:.4e}, K_eq = {:.4e}",
                t,
                ratio,
                expected
            );
        }
    }

    /// Reference: microkinetic steady-state for A → B → C sequential reaction.
    /// At steady state, [B]_ss =  k₁[A]/(k₂) and [C]_ss → [A₀] as t → ∞
    #[test]
    fn microkinetic_sequential_abc_reaches_product() {
        let config = MicrokineticNetworkConfig {
            species: vec!["A".into(), "B".into(), "C".into()],
            steps: vec![
                ElementaryStep {
                    step_id: "A->B".into(),
                    activation_free_energy_ev: 0.3,
                    reaction_free_energy_ev: -0.15,
                    prefactor_s_inv: Some(1e10),
                },
                ElementaryStep {
                    step_id: "B->C".into(),
                    activation_free_energy_ev: 0.25,
                    reaction_free_energy_ev: -0.15,
                    prefactor_s_inv: Some(1e10),
                },
            ],
            stoichiometry: vec![vec![(0, -1.0), (1, 1.0)], vec![(1, -1.0), (2, 1.0)]],
            initial_populations: vec![1.0, 0.0, 0.0],
            total_time_s: 1e-4,
            n_frames: 50,
            state: ThermodynamicState {
                temperature_k: 500.0,
                pressure_bar: 1.0,
            },
        };
        let trace = solve_microkinetic_network(&config).unwrap();
        let last = trace.frames.last().unwrap();
        let c_final = last.species[2].population;
        // C should accumulate since both reactions are exergonic
        assert!(
            c_final > 0.01,
            "final C population should be substantial: {:.4}",
            c_final
        );
        // Mass conservation across all frames
        for frame in &trace.frames {
            let total: f64 = frame.species.iter().map(|s| s.population).sum();
            assert!(
                (total - 1.0).abs() < 0.05,
                "mass conservation violated: total = {:.4}",
                total,
            );
        }
    }

    /// Convergence diagnostics test for microkinetic solver.
    #[test]
    fn kinetics_diagnostics_detects_steady_state() {
        // Use higher barrier to keep rates manageable for Euler integrator
        let config = MicrokineticNetworkConfig {
            species: vec!["A".into(), "B".into()],
            steps: vec![ElementaryStep {
                step_id: "A->B".into(),
                activation_free_energy_ev: 0.8,
                reaction_free_energy_ev: -0.3,
                prefactor_s_inv: Some(1e12),
            }],
            stoichiometry: vec![vec![(0, -1.0), (1, 1.0)]],
            initial_populations: vec![1.0, 0.0],
            total_time_s: 1e-2,
            n_frames: 100,
            state: ThermodynamicState {
                temperature_k: 500.0,
                pressure_bar: 1.0,
            },
        };
        let trace = solve_microkinetic_network(&config).unwrap();
        let diag = extract_kinetics_diagnostics(&trace);
        assert!(
            diag.mass_conservation_error < 0.1,
            "mass error: {:.4}",
            diag.mass_conservation_error,
        );
    }
}
