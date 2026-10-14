//! F82: growing-string TS on Muller-Brown vs ASE CI-NEB.
//!
//! Fixture `tests/fixtures/validation/gsm_mullerbrown.json` (oracle
//! `scripts/validation/oracle_gsm.py`): minima from scipy, MEP + TS from
//! ASE CI-NEB/FIRE (external path optimizer) on the textbook MB surface.
//! The gate runs OUR `find_transition_state` on a Rust MB closure over
//! the same endpoints and must reproduce barrier/TS location.
//! (MB units treated as kcal/mol; the energy closure is the adapter,
//! the path search is what is gated.)

#![cfg(feature = "alpha-gsm")]

use sci_form::alpha::gsm::{find_transition_state, grow_string, interpolate_node, GsmConfig};

#[derive(serde::Deserialize)]
struct Minimum {
    coords: [f64; 2],
    energy: f64,
}

#[derive(serde::Deserialize)]
struct MepPoint {
    coords: [f64; 2],
    energy: f64,
}

#[derive(serde::Deserialize)]
struct Fixture {
    minima: Vec<Minimum>,
    ts_energy: f64,
    ts_coords: [f64; 2],
    barrier: f64,
    mep: Vec<MepPoint>,
    barrier_tolerance_rel: f64,
    ts_location_tolerance: f64,
}

fn load() -> Fixture {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/validation/gsm_mullerbrown.json"
    );
    serde_json::from_str(&std::fs::read_to_string(path).expect("gsm_mullerbrown.json"))
        .expect("parse fixture")
}

// Muller-Brown surface, textbook constants.
fn mb(p: &[f64]) -> f64 {
    const A: [f64; 4] = [-200.0, -100.0, -170.0, 15.0];
    const AA: [f64; 4] = [-1.0, -1.0, -6.5, 0.7];
    const BB: [f64; 4] = [0.0, 0.0, 11.0, 0.6];
    const X0: [f64; 4] = [1.0, 0.0, -0.5, -1.0];
    const Y0: [f64; 4] = [0.0, 0.5, 1.5, 1.0];
    const CC: [f64; 4] = [-10.0, -10.0, -6.5, 0.7];
    let (x, y) = (p[0], p[1]);
    let mut v = 0.0;
    for k in 0..4 {
        v += A[k]
            * (AA[k] * (x - X0[k]).powi(2)
                + BB[k] * (x - X0[k]) * (y - Y0[k])
                + CC[k] * (y - Y0[k]).powi(2))
            .exp();
    }
    v
}

/// F82a: interpolation exactness (endpoints + midpoint).
#[test]
fn f82_interpolate_exact() {
    let fix = load();
    let r = [fix.minima[0].coords[0], fix.minima[0].coords[1]];
    let p = [fix.minima[1].coords[0], fix.minima[1].coords[1]];
    for (t, want) in [(0.0, &r), (1.0, &p)] {
        let got = interpolate_node(&r, &p, t);
        assert!((got[0] - want[0]).abs() <= 1e-12, "t={t}");
        assert!((got[1] - want[1]).abs() <= 1e-12, "t={t}");
    }
    let mid = interpolate_node(&r, &p, 0.5);
    assert!((mid[0] - (r[0] + p[0]) / 2.0).abs() <= 1e-15);
    assert!((mid[1] - (r[1] + p[1]) / 2.0).abs() <= 1e-15);
    // String grows from both ends and joins.
    let g = grow_string(&r, &p, &mb, &GsmConfig::default());
    assert!(g.nodes.len() >= 2);
    println!("F82a: string nodes={} joined={} iters={}", g.nodes.len(), g.joined, g.iterations);
}

/// MB-tuned config: gradients ~1e2/mol per unit length need small steps.
fn mb_config() -> GsmConfig {
    GsmConfig { max_nodes: 21, step_size: 0.001, max_iter: 2000, ..GsmConfig::default() }
}

/// F82b: GSM transition state vs ASE CI-NEB (barrier + location).
#[test]
fn f82_ts_matches_ase_neb() {
    let fix = load();
    let reactant = [fix.minima[0].coords[0], fix.minima[0].coords[1]];
    let product = [fix.minima[1].coords[0], fix.minima[1].coords[1]];
    let cfg = mb_config();
    let res = find_transition_state(&reactant, &product, &mb, &cfg);
    let d_bar = (res.activation_energy - fix.barrier).abs() / fix.barrier.abs();
    let d_ts = ((res.ts_coords[0] - fix.ts_coords[0]).powi(2)
        + (res.ts_coords[1] - fix.ts_coords[1]).powi(2))
    .sqrt();
    let d_e = (res.ts_energy - fix.ts_energy).abs() / fix.ts_energy.abs();
    println!(
        "F82b GSM vs ASE: barrier {:.4} (ref {:.4}, d={:.3e}) TS [{:.4},{:.4}] (d={:.4}) nodes={}",
        res.activation_energy, fix.barrier, d_bar,
        res.ts_coords[0], res.ts_coords[1], d_ts, res.n_nodes
    );
    assert!(d_bar <= fix.barrier_tolerance_rel, "barrier {d_bar:.3e}");
    assert!(d_ts <= fix.ts_location_tolerance, "location {d_ts:.4}");
    assert!(d_e <= fix.barrier_tolerance_rel, "ts energy {d_e:.3e}");
}
