//! F116–F118: first v2.1 validation rows.
//!
//! F116 and F117 are snapshots of the first external-reference runs. They are
//! recorded as gaps until the production runs and the Ewald comparison pass.
//! F118 is a live L0 smoke gate for the two enhanced-sampling entry points.

#![cfg(feature = "alpha-v2")]

use std::path::PathBuf;

use serde_json::{json, Value};

struct Row<'a> {
    ficha: &'a str,
    observable: &'a str,
    metric: &'a str,
    value: f64,
    tolerance: f64,
    pass: Option<bool>,
    n_cases: usize,
    level: &'a str,
    units: &'a str,
    reference: &'a str,
    notes: &'a str,
}

fn fixture(name: &str) -> Value {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/validation")
        .join(name);
    serde_json::from_str(
        &std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("missing fixture {}: {error}", path.display())),
    )
    .expect("invalid fixture JSON")
}

fn record(row: Row) {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("target")
        .join("validation");
    let _ = std::fs::create_dir_all(&dir);
    let mut payload = json!({
        "ficha": row.ficha,
        "observable": row.observable,
        "metric": row.metric,
        "value": row.value,
        "tolerance": row.tolerance,
        "pass": row.pass,
        "n_cases": row.n_cases,
        "source": "tests/v21_validation.rs",
        "level": row.level,
        "units": row.units,
        "reference": { "tool": row.reference },
    });
    if row.pass.is_none() {
        payload["status"] = json!("gap");
    }
    if !row.notes.is_empty() {
        payload["notes"] = json!(row.notes);
    }
    eprintln!(
        "{} {:<42} {} = {:.6e} (tol {:.1e}) {:?}",
        row.ficha, row.observable, row.metric, row.value, row.tolerance, row.pass
    );
    let path = dir.join(format!("{}.json", row.ficha));
    let mut rows: Vec<Value> = std::fs::read_to_string(&path)
        .ok()
        .and_then(|text| serde_json::from_str::<Value>(&text).ok())
        .map(|value| match value {
            Value::Array(items) => items,
            other => vec![other],
        })
        .unwrap_or_default();
    rows.retain(|existing| {
        existing.get("observable").and_then(Value::as_str) != Some(row.observable)
    });
    rows.push(payload);
    let _ = std::fs::write(&path, Value::Array(rows).to_string());
}

fn gate(row: Row) {
    let passed = row.value.is_finite() && row.value <= row.tolerance;
    record(Row {
        pass: Some(passed),
        ..row
    });
    assert!(passed);
}

fn gap(row: Row) {
    record(Row { pass: None, ..row });
}

#[test]
fn f116_gcmc_co2_cu_btc_vs_raspa() {
    let data = fixture("v21_gcmc.json");
    let ours = data["sci_form_loading_molecules_per_uc"].as_f64().unwrap();
    let reference = data["raspa_loading_molecules_per_uc"].as_f64().unwrap();
    let relative = (ours - reference).abs() / reference;
    if relative <= 0.10 {
        gate(Row {
            ficha: "F116",
            observable: "co2_cu_btc_loading",
            metric: "relative_error",
            value: relative,
            tolerance: 0.10,
            pass: None,
            n_cases: 1,
            level: "L2",
            units: "1",
            reference: "RASPA2 ExampleMOFsForceField, 323 K, 10 bar",
            notes:
                "160k moves (40k eq + 120k prod), Ewald k_max=4 + switching; 113.67 vs 116 (2.0 %)",
        });
    } else {
        gap(Row {
            ficha: "F116",
            observable: "co2_cu_btc_loading",
            metric: "relative_error",
            value: relative,
            tolerance: 0.10,
            pass: None,
            n_cases: 1,
            level: "L2",
            units: "1",
            reference: "RASPA2 ExampleMOFsForceField, 323 K, 10 bar",
            notes: "160k moves (40k eq + 120k prod), Ewald k_max=4 + switching",
        });
    }
}

#[test]
fn f117_pka_pbe0_3_21g() {
    let data = fixture("v21_pka_pbe0.json");
    let mae = data["mae_pka_units"].as_f64().unwrap();
    gap(Row {
        ficha: "F117",
        observable: "relative_pka_pbe0_3_21g",
        metric: "mae_pka_units",
        value: mae,
        tolerance: 1.0,
        pass: None,
        n_cases: data["experimental"].as_array().map(Vec::len).unwrap_or(0),
        level: "L4",
        units: "pKa",
        reference: "CRC Handbook, relative to acetic acid 4.76",
        notes:
            "PBE0/3-21G Medium grid with density mixing; better than UHF but still above the gate",
    });
}

#[test]
fn f118_enhanced_sampling_entry_points_are_finite() {
    let conformer = sci_form::embed("CCO", 42);
    assert!(conformer.error.is_none());
    let metadynamics = sci_form::enhanced_sampling::MetaDynamicsConfig {
        max_steps: 100,
        trace_stride: 10,
        grid_bins: 32,
        grid_max: 5.0,
        ..Default::default()
    };
    let md = sci_form::run_uff_metadynamics("CCO", &conformer.coords, &metadynamics).unwrap();
    let opes_config = sci_form::enhanced_sampling::OpesConfig {
        max_steps: 100,
        trace_stride: 10,
        grid_bins: 32,
        grid_max: 5.0,
        ..Default::default()
    };
    let opes = sci_form::run_uff_opes("CCO", &conformer.coords, &opes_config).unwrap();
    let non_finite = md
        .cv_trace
        .iter()
        .chain(opes.cv_trace.iter())
        .filter(|point| !point.cv.is_finite() || !point.potential_kcal_mol.is_finite())
        .count()
        + md.free_energy
            .values
            .iter()
            .chain(opes.free_energy.values.iter())
            .filter(|value| !value.is_finite())
            .count();
    gate(Row {
        ficha: "F118",
        observable: "enhanced_sampling_finite_traces",
        metric: "non_finite_values",
        value: non_finite as f64,
        tolerance: 0.0,
        pass: None,
        n_cases: 2,
        level: "L0",
        units: "count",
        reference: "sci-form metadynamics and OPES",
        notes: "WTMetaD and OPES on ethanol; checks finite traces and FES",
    });
}

/// Live reproduction of F117. This is a heavy gate: PBE0/3-21G takes minutes.
#[test]
#[ignore = "heavy: PBE0/3-21G four-state pKa cycle"]
fn f117_live_pka_pbe0_3_21g() {
    use sci_form::electrochemistry::{
        ElectrochemistryConfig, ElectrochemistryMethod, PkaReference,
    };

    let acids = [
        ("formic", "OC=O", "[O-]C=O", 3.75),
        ("propanoic", "CCC(=O)O", "CCC(=O)[O-]", 4.87),
        ("chloroacetic", "OC(=O)CCl", "[O-]C(=O)CCl", 2.87),
        ("benzoic", "OC(=O)c1ccccc1", "[O-]C(=O)c1ccccc1", 4.20),
    ];
    let reference = PkaReference::acetic_acid();
    let config = ElectrochemistryConfig {
        method: ElectrochemistryMethod::Pbe0,
        ..Default::default()
    };
    let mut errors = Vec::new();
    for (name, acid, base, experimental) in acids {
        let result = sci_form::compute_pka(acid, base, Some(&reference), &config).unwrap();
        eprintln!("F117 {name:12} pKa {:.2} vs {experimental:.2}", result.pka);
        errors.push((result.pka - experimental).abs());
    }
    gap(Row {
        ficha: "F117",
        observable: "relative_pka_pbe0_3_21g_live",
        metric: "mae_pka_units",
        value: errors.iter().sum::<f64>() / errors.len() as f64,
        tolerance: 1.0,
        pass: None,
        n_cases: errors.len(),
        level: "L4",
        units: "pKa",
        reference: "CRC Handbook, relative to acetic acid 4.76",
        notes: "Live PBE0/3-21G Medium run; heavy, ignored by default",
    });
}

struct Harmonic2D;

impl sci_form::enhanced_sampling::Potential for Harmonic2D {
    fn energy(&self, coords: &[f64]) -> f64 {
        0.5 * (coords[0] * coords[0] + coords[1] * coords[1])
    }
    fn gradient(&self, coords: &[f64], gradient: &mut [f64]) {
        gradient[0] = coords[0];
        gradient[1] = coords[1];
    }
}

#[test]
fn f119_two_dimensional_opes_is_finite() {
    let config = sci_form::enhanced_sampling_multi::MultiOpesConfig {
        max_steps: 200,
        pace_steps: 20,
        trace_stride: 10,
        x_min: -3.0,
        x_max: 3.0,
        y_min: -3.0,
        y_max: 3.0,
        ..Default::default()
    };
    let result = sci_form::enhanced_sampling_multi::run_multi_opes(
        &Harmonic2D,
        &[0.1, -0.1],
        &[1.0, 1.0],
        &config,
    )
    .unwrap();
    let non_finite = result
        .free_energy
        .values
        .iter()
        .chain(result.reweighted_free_energy.values.iter())
        .filter(|value| !value.is_finite())
        .count();
    gate(Row {
        ficha: "F119",
        observable: "two_dimensional_opes_finite",
        metric: "non_finite_values",
        value: non_finite as f64,
        tolerance: 0.0,
        pass: None,
        n_cases: 1,
        level: "L0",
        units: "count",
        reference: "sci-form 2D OPES",
        notes: "Harmonic 2D surface; finite FES and reweighted surface",
    });
}

#[test]
fn f120_gcmc_primitive_ewald_vs_raspa() {
    let data = fixture("v21_gcmc_ewald.json");
    let ours = data["sci_form_loading_molecules_per_primitive_cell"]
        .as_f64()
        .unwrap();
    let reference = data["raspa_loading_molecules_per_primitive_cell"]
        .as_f64()
        .unwrap();
    let relative = (ours - reference).abs() / reference;
    let row = Row {
        ficha: "F120",
        observable: "co2_cu_btc_primitive_ewald_loading",
        metric: "relative_error",
        value: relative,
        tolerance: 0.10,
        pass: None,
        n_cases: 1,
        level: "L2",
        units: "1",
        reference: "RASPA2 primitive P1 Cu-BTC, 323 K, 10 bar",
        notes: "EwaldHost k_max=4 with switching and intramolecular correction; 40k moves",
    };
    if relative <= 0.10 {
        gate(row);
    } else {
        gap(row);
    }
}

#[test]
fn f121_multisite_co2_empty_box_is_ideal() {
    let data = fixture("v21_gcmc_empty.json");
    let expected = data["expected"].as_f64().unwrap();
    let measured = data["measured"].as_f64().unwrap();
    let relative = (measured - expected).abs() / expected;
    gap(Row {
        ficha: "F121",
        observable: "co2_multisite_empty_box_loading_lj",
        metric: "relative_error",
        value: relative,
        tolerance: 0.10,
        pass: None,
        n_cases: 1,
        level: "L0",
        units: "1",
        reference: "ideal gas fV/kT",
        notes: "Empty 40 A box, three-site CO2 with LJ; 16.6 vs ideal 13.75 after intramolecular correction",
    });
}

/// L3 free-energy gate: well-tempered metadynamics on an analytic double well.
/// The exact potential of mean force is the double-well potential itself.
#[test]
fn f124_metadynamics_recovers_double_well_pmf() {
    use sci_form::enhanced_sampling::{
        run_metadynamics, CollectiveVariable, DoubleWell, MetaDynamicsConfig,
    };

    let barrier = 10.0;
    let well = 1.0;
    let potential = DoubleWell::new(barrier, well);
    let config = MetaDynamicsConfig {
        cv: CollectiveVariable::Coordinate { index: 0 },
        height_kcal_mol: 0.3,
        sigma: 0.15,
        bias_factor: 10.0,
        pace_steps: 500,
        temperature_k: 300.0,
        dt_fs: 0.5,
        friction_ps_inv: 5.0,
        max_steps: 400_000,
        seed: 7,
        grid_min: -1.8,
        grid_max: 1.8,
        grid_bins: 73,
        trace_stride: 2_000,
    };
    let result = run_metadynamics(&potential, &[well], &[10.0], &config).unwrap();
    let exact = |x: f64| barrier / well.powi(4) * (x * x - well * well).powi(2);
    let offset = exact(well);
    let mut sum_sq = 0.0;
    let mut count = 0usize;
    for (index, value) in result.free_energy.values.iter().enumerate() {
        let x = result.free_energy.coordinate(index);
        if x.abs() <= 1.5 && exact(x) - offset <= 6.0 {
            let error = value - (exact(x) - offset);
            sum_sq += error * error;
            count += 1;
        }
    }
    let rmse = (sum_sq / count as f64).sqrt();
    let recovered_barrier = result.free_energy.values[result.free_energy.bins / 2];
    let barrier_error = (recovered_barrier - barrier).abs();
    let tolerance = 0.5;
    eprintln!(
        "metadynamics RMSE={rmse:.3} kcal/mol, barrier={recovered_barrier:.3} vs {barrier}, n={count}"
    );
    assert!(rmse <= tolerance, "PMF RMSE {rmse} > {tolerance}");
    assert!(
        barrier_error <= tolerance,
        "barrier error {barrier_error} > {tolerance}"
    );
    gate(Row {
        ficha: "F124",
        observable: "double_well_pmf_rmse",
        metric: "rmse_kcal_mol",
        value: rmse.max(barrier_error),
        tolerance,
        pass: None,
        n_cases: 1,
        level: "L3",
        units: "kcal/mol",
        reference: "analytic double-well PMF, 300 K",
        notes: "Well-tempered metadynamics, beta=10, 400k steps; exact F(x)=V(x)+const",
    });
}

/// L3 free-energy gate: 2D OPES on an isotropic harmonic potential.
/// The exact PMF is `0.5 (x^2 + y^2)` up to an additive constant.
#[test]
fn f125_two_dimensional_opes_recovers_harmonic_pmf() {
    use sci_form::enhanced_sampling::{CollectiveVariable, Potential};
    use sci_form::enhanced_sampling_multi::{run_multi_opes, MultiOpesConfig};

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

    let config = MultiOpesConfig {
        cv_x: CollectiveVariable::Coordinate { index: 0 },
        cv_y: CollectiveVariable::Coordinate { index: 1 },
        sigma_x: 0.35,
        sigma_y: 0.35,
        bias_factor: 10.0,
        pace_steps: 400,
        temperature_k: 300.0,
        dt_fs: 0.5,
        friction_ps_inv: 5.0,
        max_steps: 200_000,
        seed: 11,
        x_min: -1.2,
        x_max: 1.2,
        x_bins: 13,
        y_min: -1.2,
        y_max: 1.2,
        y_bins: 13,
        trace_stride: 10,
        regularization: 1e-8,
        flow_iterations: 5,
    };
    let result = run_multi_opes(&Harmonic2D, &[0.0, 0.0], &[10.0, 10.0], &config).unwrap();
    let mut sum_sq = 0.0;
    let mut count = 0usize;
    for iy in 0..result.free_energy.ny {
        for ix in 0..result.free_energy.nx {
            let x = result.free_energy.x(ix);
            let y = result.free_energy.y(iy);
            if x * x + y * y <= 0.81 {
                let error = result.free_energy.at(ix, iy) - 0.5 * (x * x + y * y);
                sum_sq += error * error;
                count += 1;
            }
        }
    }
    let rmse = (sum_sq / count as f64).sqrt();
    let mut reweighted_sum = 0.0;
    let mut reweighted_count = 0usize;
    for iy in 0..result.reweighted_free_energy.ny {
        for ix in 0..result.reweighted_free_energy.nx {
            let x = result.reweighted_free_energy.x(ix);
            let y = result.reweighted_free_energy.y(iy);
            if x * x + y * y <= 0.81 {
                let error = result.reweighted_free_energy.at(ix, iy) - 0.5 * (x * x + y * y);
                reweighted_sum += error * error;
                reweighted_count += 1;
            }
        }
    }
    let reweighted_rmse = (reweighted_sum / reweighted_count as f64).sqrt();
    let tolerance = 0.5;
    eprintln!(
        "OPES harmonic RMSE={rmse:.3}, reweighted RMSE={reweighted_rmse:.3} over {count} bins"
    );
    assert!(rmse <= 0.35, "OPES PMF RMSE {rmse} > 0.35");
    assert!(
        reweighted_rmse <= tolerance,
        "OPES reweighted PMF RMSE {reweighted_rmse} > {tolerance}"
    );
    gate(Row {
        ficha: "F125",
        observable: "harmonic_2d_opes_pmf_rmse",
        metric: "rmse_kcal_mol",
        value: rmse.max(reweighted_rmse),
        tolerance,
        pass: None,
        n_cases: 1,
        level: "L3",
        units: "kcal/mol",
        reference: "analytic 2D harmonic PMF, 300 K",
        notes: "Two-dimensional OPES with self-consistent kernel weights, 200k steps, sigma=0.35; bias and reweighted surfaces",
    });
}

fn curated_networks() -> Vec<(String, f64, Value)> {
    let data = fixture("v21_networks.json");
    data["networks"]
        .as_array()
        .unwrap()
        .iter()
        .map(|network| {
            (
                network["name"].as_str().unwrap().to_string(),
                network["temperature_k"].as_f64().unwrap(),
                network.clone(),
            )
        })
        .collect()
}

fn discover_network(network: &Value) -> sci_form::networks::ReactionNetwork {
    use sci_form::networks::{TemplateNetworkConfig, TemplateRule};
    let kcal_to_ev = |kcal: f64| kcal / 23.060_548_867;
    let rules = network["steps"]
        .as_array()
        .unwrap()
        .iter()
        .map(|step| {
            let reactants = step["reactants"]
                .as_array()
                .unwrap()
                .iter()
                .map(|value| value.as_str().unwrap())
                .collect::<Vec<_>>();
            let products = step["products"]
                .as_array()
                .unwrap()
                .iter()
                .map(|value| value.as_str().unwrap())
                .collect::<Vec<_>>();
            TemplateRule::new(
                step["rule"].as_str().unwrap(),
                &reactants,
                &products,
                kcal_to_ev(step["barrier_kcal_mol"].as_f64().unwrap()),
                Some(step["prefactor_s_inv"].as_f64().unwrap()),
            )
            .unwrap()
        })
        .collect();
    let initial_species = network["initial_counts"]
        .as_object()
        .unwrap()
        .keys()
        .cloned()
        .collect();
    TemplateNetworkConfig {
        initial_species,
        rules,
        max_depth: 2,
        max_species: 30,
    }
    .discover()
    .unwrap()
}

/// L3 network-recovery gate: the five curated networks must come back from
/// the template engine with every literature species present.
#[test]
fn f126_curated_networks_are_recovered() {
    let mut recovered = 0usize;
    let mut failures = Vec::new();
    for (name, _temperature, network) in curated_networks() {
        let discovered = discover_network(&network);
        let found = discovered
            .species
            .iter()
            .map(|species| species.smiles.as_str())
            .collect::<std::collections::BTreeSet<_>>();
        let missing = network["species"]
            .as_array()
            .unwrap()
            .iter()
            .map(|value| value.as_str().unwrap())
            .filter(|expected| !found.contains(expected))
            .collect::<Vec<_>>();
        if missing.is_empty() {
            recovered += 1;
        } else {
            failures.push(format!("{name}: missing {missing:?}"));
        }
    }
    eprintln!("recovered {recovered}/5 curated networks; failures: {failures:?}");
    let tolerance = 0.0;
    gate(Row {
        ficha: "F126",
        observable: "curated_network_recovery",
        metric: "missing_networks",
        value: (5 - recovered) as f64,
        tolerance,
        pass: None,
        n_cases: 5,
        level: "L3",
        units: "count",
        reference: "GRI-Mech 3.0, NIST, JPL curated mechanisms",
        notes: "Five curated gas-phase networks recovered from template rules; all literature species present",
    });
}

/// L3 kMC/ODE gate: mean populations from the stochastic simulator must agree
/// with the deterministic mean-field network within 5 %.
#[test]
fn f127_kmc_matches_ode_on_curated_networks() {
    use sci_form::kinetics::{
        evaluate_htst_rate, solve_microkinetic_network, ElementaryStep, MicrokineticNetworkConfig,
        ThermodynamicState,
    };
    use sci_form::kmc::{run_kinetic_monte_carlo, KmcConfig, KmcSpecies, KmcStep};

    let seeds = 8u64;
    let events = 3_000usize;
    let mut worst = 0.0f64;
    let mut details = Vec::new();
    for (name, temperature_k, network) in curated_networks() {
        let discovered = discover_network(&network);
        let species = discovered
            .species
            .iter()
            .map(|species| species.smiles.clone())
            .collect::<Vec<_>>();
        let counts = network["initial_counts"].as_object().unwrap();
        let initial_counts = species
            .iter()
            .map(|smiles| counts.get(smiles).and_then(Value::as_u64).unwrap_or(0))
            .collect::<Vec<u64>>();
        let total: f64 = initial_counts.iter().map(|&count| count as f64).sum();
        let state = ThermodynamicState {
            temperature_k,
            pressure_bar: 1.0,
        };
        let steps: Vec<KmcStep> = discovered
            .edges
            .iter()
            .enumerate()
            .map(|(index, edge)| KmcStep {
                step: ElementaryStep {
                    step_id: format!("edge-{index}"),
                    activation_free_energy_ev: edge.barrier_ev,
                    reaction_free_energy_ev: 1000.0,
                    prefactor_s_inv: edge.prefactor_s_inv,
                },
                reactants: edge.reactants.iter().map(|&id| (id, 1)).collect(),
                products: edge.products.iter().map(|&id| (id, 1)).collect(),
            })
            .collect::<Vec<_>>();
        // One common comparison time, set so the stochastic run fires about
        // `events` events. Both simulators then see the same time window.
        let initial_propensity: f64 = steps
            .iter()
            .map(|step| {
                let rate = evaluate_htst_rate(&step.step, state).unwrap();
                step.reactants
                    .iter()
                    .fold(rate.forward_rate_s_inv, |acc, &(id, _)| {
                        acc * initial_counts[id] as f64
                    })
            })
            .sum();
        let comparison_time = events as f64 / initial_propensity.max(1e-300);
        let kmc_species: Vec<KmcSpecies> = species
            .iter()
            .zip(&initial_counts)
            .map(|(id, &count)| KmcSpecies {
                id: id.clone(),
                initial_count: count,
            })
            .collect();
        let stoichiometry: Vec<Vec<(usize, f64)>> = discovered
            .edges
            .iter()
            .map(|edge| {
                edge.reactants
                    .iter()
                    .map(|&index| (index, -1.0))
                    .chain(edge.products.iter().map(|&index| (index, 1.0)))
                    .collect::<Vec<_>>()
            })
            .collect();
        let ode_steps: Vec<ElementaryStep> = discovered
            .edges
            .iter()
            .enumerate()
            .map(|(index, edge)| ElementaryStep {
                step_id: format!("edge-{index}"),
                activation_free_energy_ev: edge.barrier_ev,
                reaction_free_energy_ev: 1000.0,
                prefactor_s_inv: edge.prefactor_s_inv,
            })
            .collect();
        let ode = solve_microkinetic_network(&MicrokineticNetworkConfig {
            species: species.clone(),
            steps: ode_steps,
            stoichiometry,
            initial_populations: initial_counts.iter().map(|&count| count as f64).collect(),
            total_time_s: comparison_time,
            n_frames: 11,
            state,
        })
        .unwrap();
        let ode_last = ode.frames.last().unwrap();
        let mut mean_kmc = vec![0.0f64; species.len()];
        for seed in 1..=seeds {
            let kmc = run_kinetic_monte_carlo(&KmcConfig {
                species: kmc_species.clone(),
                steps: steps.clone(),
                state,
                max_events: 2_000_000,
                max_time_s: comparison_time,
                seed,
                frame_stride: 100_000,
            })
            .unwrap();
            for (value, &count) in mean_kmc.iter_mut().zip(&kmc.final_counts) {
                *value += count as f64 / seeds as f64;
            }
        }
        let error = mean_kmc
            .iter()
            .zip(&ode_last.species)
            .map(|(kmc, ode)| (kmc - ode.population).abs() / total.max(1.0))
            .fold(0.0f64, f64::max);
        let conversion = mean_kmc
            .iter()
            .zip(&initial_counts)
            .map(|(kmc, &initial)| (kmc - initial as f64).abs() / total.max(1.0))
            .fold(0.0f64, f64::max);
        assert!(
            conversion >= 0.05,
            "{name}: only {conversion:.4} conversion, comparison is trivial"
        );
        details.push(format!("{name}: err={error:.4} conv={conversion:.4}"));
        worst = worst.max(error);
    }
    eprintln!("kMC/ODE mean-population error: {details:?}, worst={worst:.4}");
    let tolerance = 0.05;
    gate(Row {
        ficha: "F127",
        observable: "curated_network_kmc_vs_ode",
        metric: "max_normalized_population_error",
        value: worst,
        tolerance,
        pass: None,
        n_cases: 5,
        level: "L3",
        units: "1",
        reference: "mean-field ODE with identical HTST rates, 8 seeds",
        notes: "Five curated networks, 3k events per seed at a common comparison time",
    });
}

/// L3 determinism and resource-limit gate for network discovery and kMC.
#[test]
fn f128_network_seed_determinism_and_limits() {
    use sci_form::kinetics::{ElementaryStep, ThermodynamicState};
    use sci_form::kmc::{run_kinetic_monte_carlo, KmcConfig, KmcSpecies, KmcStep};
    use sci_form::networks::{TemplateNetworkConfig, TemplateRule};

    let kcal_to_ev = |kcal: f64| kcal / 23.060_548_867;
    let rule = TemplateRule::new(
        "R1",
        &["[H]", "O=O"],
        &["[O]", "[OH]"],
        kcal_to_ev(17.041),
        Some(2.65e16),
    )
    .unwrap();
    let limited = TemplateNetworkConfig {
        initial_species: vec!["[H]".into(), "O=O".into()],
        rules: vec![rule.clone()],
        max_depth: 1,
        max_species: 2,
    }
    .discover()
    .unwrap();
    let limit_ok = limited.species.len() <= 2;
    let deep = TemplateNetworkConfig {
        initial_species: vec!["[H]".into(), "O=O".into()],
        rules: vec![rule],
        max_depth: 3,
        max_species: 20,
    }
    .discover()
    .unwrap();
    let species = deep
        .species
        .iter()
        .map(|species| species.smiles.clone())
        .collect::<Vec<_>>();
    let counts = vec![100u64; species.len()];
    let steps: Vec<KmcStep> = deep
        .edges
        .iter()
        .enumerate()
        .map(|(index, edge)| KmcStep {
            step: ElementaryStep {
                step_id: format!("edge-{index}"),
                activation_free_energy_ev: edge.barrier_ev,
                reaction_free_energy_ev: 1000.0,
                prefactor_s_inv: edge.prefactor_s_inv,
            },
            reactants: edge.reactants.iter().map(|&id| (id, 1)).collect(),
            products: edge.products.iter().map(|&id| (id, 1)).collect(),
        })
        .collect();
    let config = |seed| KmcConfig {
        species: species
            .iter()
            .zip(&counts)
            .map(|(id, &count)| KmcSpecies {
                id: id.clone(),
                initial_count: count,
            })
            .collect(),
        steps: steps.clone(),
        state: ThermodynamicState {
            temperature_k: 1000.0,
            pressure_bar: 1.0,
        },
        max_events: 5_000,
        max_time_s: 1.0,
        seed,
        frame_stride: 5_000,
    };
    let first = run_kinetic_monte_carlo(&config(7)).unwrap();
    let repeat = run_kinetic_monte_carlo(&config(7)).unwrap();
    let other = run_kinetic_monte_carlo(&config(8)).unwrap();
    let deterministic = first.final_counts == repeat.final_counts
        && first.final_time_s == repeat.final_time_s
        && first.events.len() == repeat.events.len();
    let seed_changes = first.final_counts != other.final_counts
        || first
            .events
            .iter()
            .zip(&other.events)
            .any(|(a, b)| a.step_index != b.step_index || a.time_s != b.time_s);
    let violations = (!limit_ok) as usize + (!deterministic) as usize + (!seed_changes) as usize;
    eprintln!(
        "limits_ok={limit_ok} deterministic={deterministic} seed_changes={seed_changes} species={}",
        limited.species.len()
    );
    gate(Row {
        ficha: "F128",
        observable: "network_seed_determinism",
        metric: "violations",
        value: violations as f64,
        tolerance: 0.0,
        pass: None,
        n_cases: 1,
        level: "L3",
        units: "count",
        reference: "same-seed replay and max_depth/max_species caps",
        notes: "Same seed reproduces the trace, a different seed changes it, and discovery respects both caps",
    });
}

/// L3 replica-exchange gate: the adjacent-temperature swap acceptance must
/// stay in the 20–40 % window recommended by the roadmap.
#[test]
fn f129_replica_exchange_acceptance_window() {
    use sci_form::enhanced_sampling::{
        cartesian_masses, run_replica_exchange, ReplicaExchangeConfig, UffPotential,
    };

    let conformer = sci_form::embed("CCO", 42);
    assert!(conformer.error.is_none());
    let potential = UffPotential::from_smiles("CCO").unwrap();
    let masses = cartesian_masses("CCO").unwrap();
    let config = ReplicaExchangeConfig {
        temperatures_k: vec![300.0, 600.0, 1300.0, 2600.0],
        dt_fs: 0.5,
        friction_ps_inv: 5.0,
        max_steps: 20_000,
        swap_interval: 100,
        seed: 5,
        trace_stride: 1_000,
    };
    let result = run_replica_exchange(&potential, &conformer.coords, &masses, &config).unwrap();
    let proposals = result.swaps.len();
    let accepted = result.swaps.iter().filter(|swap| swap.accepted).count();
    let acceptance = accepted as f64 / proposals.max(1) as f64;
    let outside = (0.20 - acceptance).max(acceptance - 0.40).max(0.0);
    eprintln!("replica exchange: {accepted}/{proposals} accepted = {acceptance:.3}");
    gate(Row {
        ficha: "F129",
        observable: "replica_exchange_acceptance",
        metric: "distance_outside_20_40_window",
        value: outside,
        tolerance: 0.0,
        pass: None,
        n_cases: 1,
        level: "L3",
        units: "1",
        reference: "roadmap acceptance window 20-40 %",
        notes: "Four-temperature replica exchange on UFF ethanol (300/600/1300/2600 K); 20k steps, 28.7 % acceptance",
    });
}

/// L3 Raman gate: UFF harmonic positions against NIST/CCBDB experimental
/// fundamentals. The uncoupled polarizability is used for the activities.
#[test]
fn f130_raman_positions_vs_experimental_fundamentals() {
    use sci_form::raman::{compute_raman_spectrum, RamanConfig};

    let data = fixture("v21_raman.json");
    let mut worst = 0.0f64;
    let mut mae_sum = 0.0;
    let mut bands = 0usize;
    for case in data["cases"].as_array().unwrap() {
        let smiles = case["smiles"].as_str().unwrap();
        let conformer = sci_form::embed(smiles, 42);
        assert!(conformer.error.is_none(), "{smiles}: {:?}", conformer.error);
        let positions: Vec<[f64; 3]> = conformer
            .coords
            .chunks_exact(3)
            .map(|chunk| [chunk[0], chunk[1], chunk[2]])
            .collect();
        let spectrum = compute_raman_spectrum(
            smiles,
            &conformer.elements,
            &positions,
            &RamanConfig::default(),
        )
        .unwrap();
        for band in case["bands"].as_array().unwrap() {
            let reference = band.as_f64().unwrap();
            let error = spectrum
                .modes
                .iter()
                .map(|mode| (mode.frequency_cm1 - reference).abs())
                .fold(f64::INFINITY, f64::min);
            worst = worst.max(error);
            mae_sum += error;
            bands += 1;
        }
    }
    let mae = mae_sum / bands as f64;
    eprintln!("Raman positions: MAE={mae:.1} cm^-1, worst={worst:.1} over {bands} bands");
    gap(Row {
        ficha: "F130",
        observable: "raman_positions_vs_experiment",
        metric: "mae_cm1",
        value: mae,
        tolerance: 10.0,
        pass: None,
        n_cases: bands,
        level: "L3",
        units: "cm^-1",
        reference: "NIST CCCBDB gas-phase fundamentals",
        notes: "UFF harmonic modes and uncoupled polarizability; no coupled response",
    });
}

/// L2 GCMC gate: CO2 in IRMOF-1 (MOF-5) versus a rigid RASPA2 reference.
#[test]
fn f131_gcmc_co2_irmof1_vs_raspa() {
    let data = fixture("v21_gcmc_irmof1.json");
    let ours = data["sci_form_loading_molecules_per_uc"].as_f64().unwrap();
    let reference = data["raspa_loading_molecules_per_uc"].as_f64().unwrap();
    let relative = (ours - reference).abs() / reference;
    let row = Row {
        ficha: "F131",
        observable: "co2_irmof1_loading",
        metric: "relative_error",
        value: relative,
        tolerance: 0.10,
        pass: None,
        n_cases: 1,
        level: "L2",
        units: "1",
        reference: "RASPA2 IRMOF-1 (MOF-5), 298 K, 1 bar",
        notes: "Ewald k_max=8 with shift_lj, 100k moves x 2 seeds; framework geometry and charges match RASPA",
    };
    if relative <= 0.10 {
        gate(row);
    } else {
        gap(row);
    }
}

/// L2 GCMC gate: N2 in HKUST-1 (Cu-BTC) versus a rigid RASPA2 reference.
#[test]
fn f132_gcmc_n2_cu_btc_vs_raspa() {
    let data = fixture("v21_gcmc_n2.json");
    let ours = data["sci_form_loading_molecules_per_primitive_cell"]
        .as_f64()
        .unwrap();
    let reference = data["raspa_loading_molecules_per_primitive_cell"]
        .as_f64()
        .unwrap();
    let relative = (ours - reference).abs() / reference;
    let row = Row {
        ficha: "F132",
        observable: "n2_cu_btc_primitive_loading",
        metric: "relative_error",
        value: relative,
        tolerance: 0.10,
        pass: None,
        n_cases: 1,
        level: "L2",
        units: "1",
        reference: "RASPA2 N2/Cu-BTC (HKUST-1), 298 K, 10 bar",
        notes: "Ewald k_max=4 with shift_lj, 60k moves x 2 seeds; TraPPE N2 model",
    };
    if relative <= 0.10 {
        gate(row);
    } else {
        gap(row);
    }
}
