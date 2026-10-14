//! F94: PSD cone ops + Gram-matrix optimizer vs numpy/scipy.
//!
//! Fixture `tests/fixtures/validation/riemannian_numpy.json` (oracle
//! `scripts/validation/oracle_riemannian.py`, seed 42).
//! Documented gap: `psd_distance` is the Frobenius proxy, NOT the
//! affine-invariant geodesic distance (geomstats cross-check is the
//! unmet promotion criterion). Gated: projection/retraction (1e-12),
//! distance proxy (1e-12), objective/gradient pointwise (1e-9),
//! L-BFGS optimizer vs scipy L-BFGS-B on feasible distance geometry.
//!
//! Core module (`sci_form::riemannian`).

use nalgebra::DMatrix;
use sci_form::riemannian::{
    psd_distance, psd_projection, psd_retraction, DistanceConstraint, RiemannianLbfgs,
};

#[derive(serde::Deserialize)]
struct ProjCase {
    n: usize,
    x: Vec<Vec<f64>>,
    xi: Vec<Vec<f64>>,
    projection_ref: Vec<Vec<f64>>,
    retraction_ref: Vec<Vec<f64>>,
}

#[derive(serde::Deserialize)]
struct DistCase {
    n: usize,
    a: Vec<Vec<f64>>,
    b: Vec<Vec<f64>>,
    frobenius_ref: f64,
}

#[derive(serde::Deserialize)]
struct EvalCase {
    n: usize,
    x: Vec<Vec<f64>>,
    constraints: Vec<[f64; 4]>,
    objective_ref: f64,
    gradient_ref: Vec<Vec<f64>>,
}

#[derive(serde::Deserialize)]
struct OptCase {
    n: usize,
    initial: Vec<Vec<f64>>,
    constraints: Vec<[f64; 4]>,
    scipy_objective: f64,
    scipy_max_violation_A: f64,
}

#[derive(serde::Deserialize)]
struct Fixture {
    projection_tolerance: f64,
    distance_tolerance: f64,
    objective_tolerance: f64,
    optimizer_violation_tolerance_A: f64,
    projections: Vec<ProjCase>,
    distances: Vec<DistCase>,
    #[serde(rename = "eval")]
    evals: Vec<EvalCase>,
    optimizer: Vec<OptCase>,
}

fn load() -> Fixture {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/validation/riemannian_numpy.json"
    );
    serde_json::from_str(&std::fs::read_to_string(path).expect("riemannian_numpy.json"))
        .expect("parse fixture")
}

fn mat(n: usize, rows: &[Vec<f64>]) -> DMatrix<f64> {
    DMatrix::from_row_slice(n, n, &rows.concat())
}

fn frob_diff(got: &DMatrix<f64>, want: &DMatrix<f64>) -> f64 {
    assert_eq!((got.nrows(), got.ncols()), (want.nrows(), want.ncols()));
    let num: f64 = (got - want).iter().map(|&v| v * v).sum();
    let den: f64 = want.iter().map(|&v| v * v).sum::<f64>().max(1e-300);
    (num / den).sqrt()
}

fn cons(raw: &[[f64; 4]]) -> Vec<DistanceConstraint> {
    raw.iter()
        .map(|c| DistanceConstraint { i: c[0] as usize, j: c[1] as usize, lower: c[2], upper: c[3] })
        .collect()
}

/// F94a: PSD projection + retraction vs numpy eigh-clip.
#[test]
fn f94_projection_matches_numpy() {
    let fix = load();
    assert!(fix.projections.len() >= 6);
    let mut worst = 0.0f64;
    for c in &fix.projections {
        let x = mat(c.n, &c.x);
        let xi = mat(c.n, &c.xi);
        let p = psd_projection(&x);
        let r = psd_retraction(&x, &xi);
        // Retraction from zero tangent = projection.
        let r0 = psd_retraction(&x, &DMatrix::zeros(c.n, c.n));
        worst = worst
            .max(frob_diff(&p, &mat(c.n, &c.projection_ref)))
            .max(frob_diff(&r, &mat(c.n, &c.retraction_ref)))
            .max(frob_diff(&r0, &p));
        // Projection output is PSD: P_+(P) = P.
        let pp = psd_projection(&p);
        assert!(frob_diff(&pp, &p) <= 1e-12, "projection not idempotent");
    }
    println!("F94a projection vs numpy: worst {worst:.3e}");
    assert!(worst <= fix.projection_tolerance);
}

/// F94b: Frobenius-proxy distance vs numpy norm.
#[test]
fn f94_distance_matches_numpy() {
    let fix = load();
    let mut worst = 0.0f64;
    for c in &fix.distances {
        let a = mat(c.n, &c.a);
        let b = mat(c.n, &c.b);
        let got = psd_distance(&a, &b);
        worst = worst.max((got - c.frobenius_ref).abs() / c.frobenius_ref.max(1e-300));
        // Symmetry + identity of indiscernibles.
        assert!((psd_distance(&b, &a) - got).abs() <= 1e-12);
        assert!(psd_distance(&a, &a) == 0.0);
    }
    println!("F94b distance vs numpy: worst {worst:.3e}");
    assert!(worst <= fix.distance_tolerance);
}

/// F94c: objective + analytic gradient vs independent numpy evaluation.
#[test]
fn f94_objective_gradient_match() {
    let fix = load();
    assert!(fix.evals.len() >= 2);
    let opt = RiemannianLbfgs::new(Default::default());
    for c in &fix.evals {
        let x = mat(c.n, &c.x);
        let cs = cons(&c.constraints);
        let f = opt.objective(&x, &cs);
        let g = opt.gradient(&x, &cs);
        let df = (f - c.objective_ref).abs() / c.objective_ref.abs().max(1e-300);
        let dg = frob_diff(&g, &mat(c.n, &c.gradient_ref));
        println!("F94c N={}: dobj {df:.3e} dgrad {dg:.3e}", c.n);
        assert!(df <= fix.objective_tolerance, "N={}: {df:.3e}", c.n);
        assert!(dg <= fix.objective_tolerance, "N={}: {dg:.3e}", c.n);
    }
}

/// Max bound violation (Å) of a Gram matrix.
fn max_violation(x: &DMatrix<f64>, cs: &[DistanceConstraint]) -> f64 {
    cs.iter()
        .map(|c| {
            let d2 = x[(c.i, c.i)] + x[(c.j, c.j)] - 2.0 * x[(c.i, c.j)];
            let d = d2.max(0.0).sqrt();
            if d < c.lower {
                c.lower - d
            } else if d > c.upper {
                d - c.upper
            } else {
                0.0
            }
        })
        .fold(0.0f64, |m, v| m.max(v))
}

/// F94d: L-BFGS optimizer reaches the feasible set like scipy L-BFGS-B.
#[test]
fn f94_optimizer_matches_scipy() {
    let fix = load();
    assert!(fix.optimizer.len() >= 2);
    let opt = RiemannianLbfgs::new(Default::default());
    for c in &fix.optimizer {
        let x0 = mat(c.n, &c.initial);
        let cs = cons(&c.constraints);
        let r = opt.optimize(&x0, &cs);
        let viol = max_violation(&r.metric_matrix, &cs);
        println!(
            "F94d N={}: obj {:.3e} (scipy {:.3e}) viol {:.3e} (scipy {:.3e}) iters={} conv={} fb={}",
            c.n,
            r.objective,
            c.scipy_objective,
            viol,
            c.scipy_max_violation_A,
            r.iterations,
            r.converged,
            r.used_fallback
        );
        assert!(
            viol <= fix.optimizer_violation_tolerance_A,
            "N={}: violation {viol:.3e}",
            c.n
        );
        assert!(r.objective <= 1e-6, "N={}: objective {}", c.n, r.objective);
    }
}
