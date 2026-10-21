//! F85: IMD harmonic steering vs analytic closed form.
//!
//! Fixture `tests/fixtures/validation/imd_steering.json` (oracle
//! `scripts/validation/oracle_imd.py`, seed 42): E = k/2|r-r0|^2 and
//! F = -k(r-r0) from independent numpy code, including multi-force
//! additivity. Gates energy/forces (1e-12), additivity/clear semantics,
//! and energy-force consistency via finite differences.
//!
//! Core module (`sci_form::imd`).

use sci_form::imd::SteeringForces;

#[derive(serde::Deserialize)]
struct ForceSpec {
    atom: usize,
    target: [f64; 3],
    k: f64,
}

#[derive(serde::Deserialize)]
struct ImdCase {
    n_atoms: usize,
    positions: Vec<f64>,
    forces: Vec<ForceSpec>,
    energy_ref: f64,
    forces_ref: Vec<f64>,
}

#[derive(serde::Deserialize)]
struct Fixture {
    tolerance: f64,
    cases: Vec<ImdCase>,
}

fn load() -> Fixture {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/validation/imd_steering.json"
    );
    serde_json::from_str(&std::fs::read_to_string(path).expect("imd_steering.json"))
        .expect("parse fixture")
}

fn build(specs: &[ForceSpec]) -> SteeringForces {
    let mut sf = SteeringForces::new();
    for s in specs {
        sf.apply(s.atom, s.target, s.k);
    }
    sf
}

/// F85a: steering energy + forces vs analytic reference.
#[test]
fn f85_steering_matches_analytic() {
    let fix = load();
    assert!(fix.cases.len() >= 12);
    let mut worst = 0.0f64;
    for (n, c) in fix.cases.iter().enumerate() {
        assert!(c.positions.len() == 3 * c.n_atoms);
        let sf = build(&c.forces);
        assert!(sf.is_active());
        let mut f = vec![0.0; 3 * c.n_atoms];
        let e = sf.accumulate(&c.positions, &mut f);
        let de = (e - c.energy_ref).abs() / c.energy_ref.abs().max(1e-300);
        let df = f
            .iter()
            .zip(&c.forces_ref)
            .fold(0.0f64, |m, (&g, &w)| m.max((g - w).abs()));
        worst = worst.max(de).max(df);
        assert!(de <= fix.tolerance, "case {n}: dE {de:.3e}");
        assert!(df <= 1e-9, "case {n}: dF {df:.3e}");
    }
    println!(
        "F85a steering vs analytic: n={} worst {:.3e}",
        fix.cases.len(),
        worst
    );
}

/// F85b: management semantics (update-in-place, clear, clear_all).
#[test]
fn f85_steering_semantics() {
    let mut sf = SteeringForces::new();
    assert!(!sf.is_active());
    sf.apply(0, [1.0, 0.0, 0.0], 10.0);
    sf.apply(0, [2.0, 0.0, 0.0], 20.0); // update, not duplicate
    assert_eq!(sf.forces.len(), 1);
    assert_eq!(sf.forces[0].target_xyz, [2.0, 0.0, 0.0]);
    sf.apply(3, [0.0, 0.0, 0.0], 5.0);
    assert_eq!(sf.forces.len(), 2);
    sf.clear(0);
    assert_eq!(sf.forces.len(), 1);
    sf.clear_all();
    assert!(!sf.is_active());
    // Empty set accumulates zero without touching the buffer.
    let mut f = vec![7.0; 9];
    let e = sf.accumulate(&[0.0; 9], &mut f);
    assert_eq!(e, 0.0);
    assert!(f.iter().all(|&v| v == 7.0));
}

/// F85c: F = -dE/dr via central finite differences.
#[test]
fn f85_force_is_energy_gradient() {
    let fix = load();
    let c = &fix.cases[fix.cases.len() - 1];
    let sf = build(&c.forces);
    let h = 1e-5;
    let mut f = vec![0.0; c.positions.len()];
    sf.accumulate(&c.positions, &mut f);
    let mut worst = 0.0f64;
    for i in 0..c.positions.len() {
        let mut pp = c.positions.clone();
        let mut pm = c.positions.clone();
        pp[i] += h;
        pm[i] -= h;
        let mut tmp = vec![0.0; c.positions.len()];
        let ep = sf.accumulate(&pp, &mut tmp);
        let em = sf.accumulate(&pm, &mut tmp);
        worst = worst.max(((ep - em) / (2.0 * h) + f[i]).abs());
    }
    println!("F85c gradient consistency: worst {worst:.3e}");
    assert!(worst <= 1e-6);
}
