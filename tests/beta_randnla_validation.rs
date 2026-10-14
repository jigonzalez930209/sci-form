//! F93: randomized Nyström vs numpy/scipy deterministic references.
//!
//! Fixture `tests/fixtures/validation/randnla_numpy.json` (oracle
//! `scripts/validation/oracle_randnla.py`, seed 42).
//! Finding (F95-style quantified gap): single-shot rank-N/4 sketches
//! resolve clustered overlap-like spectra to ~1e-1, not 1e-3 — the 1e-3
//! promotion criterion needs power iterations/oversampling. Gated:
//! characterization (tol 0.5, locks measured accuracy), randomized path
//! in the near-identity regime with fallback OFF (5e-3), default config
//! with fallback ON (exact, 1e-6).

#![cfg(feature = "beta-randnla")]

use nalgebra::DMatrix;
use rand::{rngs::StdRng, SeedableRng};
use sci_form::beta::rand_nla::{solve_eht_randnla, GaussianSketch, NystromApprox, RandNlaConfig};

#[derive(serde::Deserialize)]
struct CharCase {
    n: usize,
    k: usize,
    matrix: Vec<Vec<f64>>,
    top_eigenvalues: Vec<f64>,
}

#[derive(serde::Deserialize)]
struct ValidCase {
    n: usize,
    k: usize,
    matrix: Vec<Vec<f64>>,
    top_eigenvalues: Vec<f64>,
    inv_sqrt: Vec<Vec<f64>>,
    h: Vec<Vec<f64>>,
    energies: Vec<f64>,
}

#[derive(serde::Deserialize)]
struct ProdCase {
    n: usize,
    h: Vec<Vec<f64>>,
    s: Vec<Vec<f64>>,
    energies: Vec<f64>,
}

#[derive(serde::Deserialize)]
struct Fixture {
    characterization_tolerance_rel: f64,
    valid_regime_tolerance: f64,
    inv_sqrt_tolerance_rel: f64,
    production_tolerance: f64,
    characterization: Vec<CharCase>,
    valid_regime: Vec<ValidCase>,
    production: Vec<ProdCase>,
}

fn load() -> Fixture {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/validation/randnla_numpy.json"
    );
    serde_json::from_str(&std::fs::read_to_string(path).expect("randnla_numpy.json"))
        .expect("parse fixture")
}

fn mat(n: usize, rows: &[Vec<f64>]) -> DMatrix<f64> {
    DMatrix::from_row_slice(n, n, &rows.concat())
}

fn sketch_approx(s: &DMatrix<f64>, n: usize, k: usize) -> NystromApprox {
    let mut rng = StdRng::seed_from_u64(42);
    let sketch = GaussianSketch::new(&mut rng, n, k);
    NystromApprox::from_matrix(s, &sketch)
}

/// F93a: characterization — top-k on realistic overlaps at rank N/4.
/// Locks measured accuracy (~1e-1); 1e-3 needs power iterations.
#[test]
fn f93_characterization_rank_quarter() {
    let fix = load();
    assert!(fix.characterization.len() >= 2);
    for c in &fix.characterization {
        let s = mat(c.n, &c.matrix);
        let approx = sketch_approx(&s, c.n, c.k);
        let (evals, _) = approx.eigendecompose();
        let mut got: Vec<f64> = evals.iter().copied().collect();
        got.sort_by(|a, b| b.partial_cmp(a).unwrap());
        let worst = got[..c.k]
            .iter()
            .zip(&c.top_eigenvalues)
            .fold(0.0f64, |m, (&g, &w)| m.max((g - w).abs() / w.abs().max(1e-300)));
        println!("F93a N={} k={}: top-k worst {:.3e}", c.n, c.k, worst);
        assert!(worst <= fix.characterization_tolerance_rel, "N={}: {worst:.3e}", c.n);
    }
}

/// F93b: randomized path (fallback OFF) in the near-identity regime.
#[test]
fn f93_valid_regime_randomized() {
    let fix = load();
    assert!(fix.valid_regime.len() >= 3);
    for c in &fix.valid_regime {
        let s = mat(c.n, &c.matrix);
        let h = mat(c.n, &c.h);
        // Top-k sanity in the valid regime.
        let approx = sketch_approx(&s, c.n, c.k);
        let (evals, _) = approx.eigendecompose();
        let mut got: Vec<f64> = evals.iter().copied().collect();
        got.sort_by(|a, b| b.partial_cmp(a).unwrap());
        let topk = got[..c.k]
            .iter()
            .zip(&c.top_eigenvalues)
            .fold(0.0f64, |m, (&g, &w)| m.max((g - w).abs() / w.abs().max(1e-300)));
        // Generalized solve, randomized path only.
        let cfg = RandNlaConfig {
            sketch_size: Some(c.k),
            seed: 42,
            max_error: 1.0,
            fallback_enabled: false,
        };
        let (energies, _, info) = solve_eht_randnla(&h, &s, &cfg);
        assert!(!info.used_fallback);
        let width = c.energies.last().unwrap() - c.energies.first().unwrap();
        let worst = energies
            .iter()
            .zip(&c.energies)
            .fold(0.0f64, |m, (&g, &w)| m.max((g - w).abs() / width));
        // Inverse square root vs scipy sqrtm.
        let x = approx.inverse_sqrt();
        let want = mat(c.n, &c.inv_sqrt);
        let isr = (&x - &want).norm() / want.norm().max(1e-300);
        println!("F93b N={}: topk {topk:.3e} gen {worst:.3e} (resid {:.2e}) isqrt {isr:.3e}",
            c.n, info.residual_error);
        assert!(worst <= fix.valid_regime_tolerance, "N={} gen: {worst:.3e}", c.n);
        assert!(isr <= fix.inv_sqrt_tolerance_rel, "N={} isqrt: {isr:.3e}", c.n);
    }
}

/// F93c: production default config (fallback ON) matches scipy exactly.
#[test]
fn f93_production_default_exact() {
    let fix = load();
    assert!(fix.production.len() >= 3);
    for c in &fix.production {
        let h = mat(c.n, &c.h);
        let s = mat(c.n, &c.s);
        let cfg = RandNlaConfig { seed: 42, ..RandNlaConfig::default() };
        let (energies, _, info) = solve_eht_randnla(&h, &s, &cfg);
        let width = c.energies.last().unwrap() - c.energies.first().unwrap();
        let worst = energies
            .iter()
            .zip(&c.energies)
            .fold(0.0f64, |m, (&g, &w)| m.max((g - w).abs() / width));
        println!("F93c N={}: worst {worst:.3e} fallback={} k={}", c.n, info.used_fallback, info.k);
        assert!(worst <= fix.production_tolerance, "N={}: {worst:.3e}", c.n);
    }
}
