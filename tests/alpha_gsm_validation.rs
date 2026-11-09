//! F82: growing-string TS on Muller-Brown vs ASE CI-NEB.
//!
//! Fixture `tests/fixtures/validation/gsm_mullerbrown.json` (oracle
//! `scripts/validation/oracle_gsm.py`): minima from scipy, MEP + TS from
//! ASE CI-NEB/FIRE (external path optimizer) on the textbook MB surface.
//! The gate runs OUR `find_transition_state` on a Rust MB closure over
//! the same endpoints and must reproduce barrier/TS location.
//! (MB units treated as kcal/mol; the energy closure is the adapter,
//! the path search is what is gated.)

use sci_form::gsm::{find_transition_state, grow_string, interpolate_node, GsmConfig};

#[derive(serde::Deserialize)]
struct Minimum {
    /// Mirrors the fixture schema; the reference minimum is only used for
    /// the endpoint checks, not the sampled energy.
    #[allow(dead_code)]
    coords: [f64; 2],
    #[allow(dead_code)]
    energy: f64,
}

#[derive(serde::Deserialize)]
struct MepPoint {
    /// Mirrors the fixture schema; the reference MEP is only used for the
    /// shape checks, not the sampled values.
    #[allow(dead_code)]
    coords: [f64; 2],
    #[allow(dead_code)]
    energy: f64,
}

#[derive(serde::Deserialize)]
struct Fixture {
    minima: Vec<Minimum>,
    ts_energy: f64,
    ts_coords: [f64; 2],
    barrier: f64,
    /// Mirrors the fixture schema; the reference MEP is only used for the
    /// shape checks, not the sampled values.
    #[allow(dead_code)]
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
    println!(
        "F82a: string nodes={} joined={} iters={}",
        g.nodes.len(),
        g.joined,
        g.iterations
    );
}

/// MB-tuned config: gradients ~1e2/mol per unit length need small steps.
fn mb_config() -> GsmConfig {
    GsmConfig {
        max_nodes: 21,
        step_size: 0.001,
        max_iter: 2000,
        ..GsmConfig::default()
    }
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
        res.activation_energy,
        fix.barrier,
        d_bar,
        res.ts_coords[0],
        res.ts_coords[1],
        d_ts,
        res.n_nodes
    );
    assert!(d_bar <= fix.barrier_tolerance_rel, "barrier {d_bar:.3e}");
    assert!(d_ts <= fix.ts_location_tolerance, "location {d_ts:.4}");
    assert!(d_e <= fix.barrier_tolerance_rel, "ts energy {d_e:.3e}");
}

// ── Five curved double wells: GSM vs a discrete NEB reference ──────────────
//
// V(x,y) = (x² − a²)² + ½·(y − c·x²)² has minima at (±a, c·a²), a saddle at
// the origin and an exact barrier a⁴. The MEP is the parabola y = c·x², so a
// straight-line interpolation is intentionally wrong and the path optimizer
// has to bend towards y=0 at the saddle.
mod gsm_vs_neb_five_reactions {
    use super::*;

    fn curved_well(a: f64, c: f64, p: &[f64]) -> f64 {
        let (x, y) = (p[0], p[1]);
        (x * x - a * a).powi(2) + 0.5 * (y - c * x * x).powi(2)
    }

    fn gradient(a: f64, c: f64, p: &[f64]) -> [f64; 2] {
        let h = 1e-6;
        let ex = curved_well(a, c, &[p[0] + h, p[1]]);
        let emx = curved_well(a, c, &[p[0] - h, p[1]]);
        let ey = curved_well(a, c, &[p[0], p[1] + h]);
        let emy = curved_well(a, c, &[p[0], p[1] - h]);
        [(ex - emx) / (2.0 * h), (ey - emy) / (2.0 * h)]
    }

    /// Minimal discrete NEB: fixed endpoints, spring coupling along the
    /// tangent, gradient descent on the perpendicular force.
    fn neb_barrier(a: f64, c: f64, n_images: usize, step: f64, iters: usize) -> f64 {
        let start = [-a, c * a * a];
        let end = [a, c * a * a];
        let mut images: Vec<[f64; 2]> = (0..n_images)
            .map(|i| {
                let t = i as f64 / (n_images - 1) as f64;
                [
                    start[0] + t * (end[0] - start[0]),
                    start[1] + t * (end[1] - start[1]),
                ]
            })
            .collect();
        let spring = 0.1;
        for _ in 0..iters {
            let old = images.clone();
            for i in 1..n_images - 1 {
                let prev = old[i - 1];
                let next = old[i + 1];
                let mut tangent = [next[0] - prev[0], next[1] - prev[1]];
                let norm = (tangent[0] * tangent[0] + tangent[1] * tangent[1]).sqrt();
                if norm < 1e-12 {
                    continue;
                }
                tangent = [tangent[0] / norm, tangent[1] / norm];
                let g = gradient(a, c, &old[i]);
                let g_perp = [
                    g[0] - (g[0] * tangent[0] + g[1] * tangent[1]) * tangent[0],
                    g[1] - (g[0] * tangent[0] + g[1] * tangent[1]) * tangent[1],
                ];
                let d_prev = ((old[i][0] - prev[0]).powi(2) + (old[i][1] - prev[1]).powi(2)).sqrt();
                let d_next = ((next[0] - old[i][0]).powi(2) + (next[1] - old[i][1]).powi(2)).sqrt();
                let spring_force = spring * (d_next - d_prev);
                images[i][0] -= step * (g_perp[0] - spring_force * tangent[0]);
                images[i][1] -= step * (g_perp[1] - spring_force * tangent[1]);
            }
        }
        let e0 = curved_well(a, c, &start);
        images
            .iter()
            .map(|p| curved_well(a, c, p))
            .fold(f64::NEG_INFINITY, f64::max)
            - e0
    }

    #[test]
    fn gsm_matches_neb_on_five_curved_wells() {
        let config = GsmConfig {
            max_nodes: 41,
            step_size: 0.01,
            max_iter: 5000,
            ..GsmConfig::default()
        };
        let mut worst_barrier = 0.0f64;
        let mut worst_gsm_vs_neb = 0.0f64;
        for (idx, &a) in [0.8f64, 1.0, 1.2, 1.5, 2.0].iter().enumerate() {
            let c = 0.5;
            let start = [-a, c * a * a];
            let end = [a, c * a * a];
            let potential = |p: &[f64]| curved_well(a, c, p);
            let gsm = find_transition_state(&start, &end, &potential, &config);
            let neb = neb_barrier(a, c, 21, 0.002, 4000);
            let exact = a.powi(4);
            let gsm_err = (gsm.activation_energy - exact).abs() / exact;
            let neb_err = (neb - exact).abs() / exact;
            let cross = (gsm.activation_energy - neb).abs() / exact;
            eprintln!(
                "GSM/NEB #{idx} a={a}: GSM {:.5} NEB {:.5} exact {:.5} (gsm err {:.2e}, neb err {:.2e}, cross {:.2e})",
                gsm.activation_energy, neb, exact, gsm_err, neb_err, cross
            );
            worst_barrier = worst_barrier.max(gsm_err).max(neb_err);
            worst_gsm_vs_neb = worst_gsm_vs_neb.max(cross);
            assert!(
                (gsm.ts_coords[0]).abs() < 0.15 * a,
                "a={a}: GSM TS x={}",
                gsm.ts_coords[0]
            );
        }
        assert!(
            worst_barrier < 0.05,
            "worst analytic barrier error {worst_barrier:.3e}"
        );
        assert!(
            worst_gsm_vs_neb < 0.05,
            "worst GSM/NEB cross-check {worst_gsm_vs_neb:.3e}"
        );
    }
}
