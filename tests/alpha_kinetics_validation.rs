//! F87: HTST/Eyring + microkinetics vs closed form + scipy.
//!
//! Fixture `tests/fixtures/validation/kinetics_ref.json` (oracle
//! `scripts/validation/oracle_kinetics.py`): Eyring closed form with
//! CODATA constants + scipy solve_ivp (RK45, rtol=1e-10) references.
//! Cantera cross-check is the unmet promotion criterion (needs full
//! mechanism/thermo I/O). Gates: rates (1e-12), traces (1e-4),
//! steady-state equilibrium ratio == K.

#![cfg(feature = "alpha-kinetics")]

use sci_form::alpha::kinetics::{
    evaluate_htst_rate, solve_microkinetic_network, solve_microkinetic_steady_state,
    ElementaryStep, MicrokineticNetworkConfig, ThermodynamicState,
};

#[derive(serde::Deserialize)]
struct HtstCase {
    #[serde(rename = "T_K")]
    t: f64,
    #[serde(rename = "dg_ts_eV")]
    dg_ts_ev: f64,
    #[serde(rename = "dg_rxn_eV")]
    dg_rxn_ev: f64,
    pref_override: Option<f64>,
    kf: f64,
    kr: f64,
    #[serde(rename = "K")]
    keq: f64,
    /// Mirrors the fixture schema; the reference prefactor is not compared.
    #[allow(dead_code)]
    pref: f64,
}

#[derive(serde::Deserialize)]
struct StepSpec {
    step_id: String,
    #[serde(rename = "dg_ts_eV")]
    dg_ts_ev: f64,
    #[serde(rename = "dg_rxn_eV")]
    dg_rxn_ev: f64,
    pref: Option<f64>,
}

#[derive(serde::Deserialize)]
struct NetCase {
    id: String,
    species: Vec<String>,
    steps: Vec<StepSpec>,
    stoichiometry: Vec<Vec<(usize, f64)>>,
    initial: Vec<f64>,
    #[serde(rename = "T_K")]
    t: f64,
    t_end_s: f64,
    times_s: Vec<f64>,
    populations: Vec<Vec<f64>>,
    #[serde(rename = "K_ref")]
    k_ref: Option<f64>,
}

#[derive(serde::Deserialize)]
struct Fixture {
    htst_tolerance_rel: f64,
    trace_tolerance_abs: f64,
    htst: Vec<HtstCase>,
    networks: Vec<NetCase>,
}

fn load() -> Fixture {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/validation/kinetics_ref.json"
    );
    serde_json::from_str(&std::fs::read_to_string(path).expect("kinetics_ref.json"))
        .expect("parse fixture")
}

fn to_step(s: &StepSpec) -> ElementaryStep {
    ElementaryStep {
        step_id: s.step_id.clone(),
        activation_free_energy_ev: s.dg_ts_ev,
        reaction_free_energy_ev: s.dg_rxn_ev,
        prefactor_s_inv: s.pref,
    }
}

/// F87a: Eyring rates + equilibrium constant vs CODATA closed form.
#[test]
fn f87_htst_matches_eyring() {
    let fix = load();
    assert!(fix.htst.len() >= 60);
    let mut worst = 0.0f64;
    for c in &fix.htst {
        let step = ElementaryStep {
            step_id: "t".into(),
            activation_free_energy_ev: c.dg_ts_ev,
            reaction_free_energy_ev: c.dg_rxn_ev,
            prefactor_s_inv: c.pref_override,
        };
        let r = evaluate_htst_rate(
            &step,
            ThermodynamicState {
                temperature_k: c.t,
                pressure_bar: 1.0,
            },
        )
        .expect("htst rate");
        for (got, want, name) in [
            (r.forward_rate_s_inv, c.kf, "kf"),
            (r.reverse_rate_s_inv, c.kr, "kr"),
            (r.equilibrium_constant, c.keq, "K"),
        ] {
            let d = (got - want).abs() / want.abs().max(1e-300);
            worst = worst.max(d);
            assert!(
                d <= fix.htst_tolerance_rel,
                "{name} T={} dg#={} dgr={}: got {got:.6e} want {want:.6e} d={d:.3e}",
                c.t,
                c.dg_ts_ev,
                c.dg_rxn_ev
            );
        }
    }
    println!(
        "F87a Eyring vs closed form: n={} worst {:.3e}",
        fix.htst.len(),
        worst
    );
}

/// F87b: microkinetic traces vs scipy solve_ivp.
#[test]
fn f87_traces_match_solve_ivp() {
    let fix = load();
    assert!(fix.networks.len() >= 2);
    for n in &fix.networks {
        let cfg = MicrokineticNetworkConfig {
            species: n.species.clone(),
            steps: n.steps.iter().map(to_step).collect(),
            stoichiometry: n.stoichiometry.clone(),
            initial_populations: n.initial.clone(),
            total_time_s: n.t_end_s,
            n_frames: n.times_s.len(),
            state: ThermodynamicState {
                temperature_k: n.t,
                pressure_bar: 1.0,
            },
        };
        let trace = solve_microkinetic_network(&cfg).expect("network");
        assert_eq!(trace.frames.len(), n.times_s.len());
        let mut worst = 0.0f64;
        for (f, (t, want)) in trace
            .frames
            .iter()
            .zip(n.times_s.iter().zip(&n.populations))
        {
            assert!((f.time_s - t).abs() <= 1e-9 * t.max(1e-300));
            for (sp, w) in f.species.iter().zip(want) {
                worst = worst.max((sp.population - w).abs());
            }
        }
        // Mass conservation on every frame.
        let total0: f64 = n.initial.iter().sum();
        for f in &trace.frames {
            let tot: f64 = f.species.iter().map(|s| s.population).sum();
            assert!((tot - total0).abs() <= 1e-9, "mass not conserved");
        }
        println!("F87b {}: worst {:.3e}", n.id, worst);
        assert!(worst <= fix.trace_tolerance_abs, "{}: {worst:.3e}", n.id);
    }
}

/// F87c: steady state of A<->B honours the equilibrium constant.
#[test]
fn f87_steady_state_matches_equilibrium() {
    let fix = load();
    let n = fix
        .networks
        .iter()
        .find(|x| x.id == "AB-reversible")
        .expect("AB case");
    let k_ref = n.k_ref.expect("K_ref");
    let cfg = MicrokineticNetworkConfig {
        species: n.species.clone(),
        steps: n.steps.iter().map(to_step).collect(),
        stoichiometry: n.stoichiometry.clone(),
        initial_populations: n.initial.clone(),
        total_time_s: 50.0 * n.t_end_s,
        n_frames: 11,
        state: ThermodynamicState {
            temperature_k: n.t,
            pressure_bar: 1.0,
        },
    };
    let trace = solve_microkinetic_steady_state(&cfg, 1e-6).expect("steady state");
    let last = trace.frames.last().expect("frames");
    let (ca, cb) = (last.species[0].population, last.species[1].population);
    let ratio = cb / ca.max(1e-300);
    let d = (ratio - k_ref).abs() / k_ref;
    println!("F87c steady state: cB/cA={ratio:.6} K={k_ref:.6} d={d:.2e}");
    assert!(d <= 1e-3, "equilibrium ratio {d:.3e}");
}
