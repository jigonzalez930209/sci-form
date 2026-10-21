//! F80: Obara-Saika Boys/ERI/Schwarz vs scipy + PySCF/libcint.
//!
//! Fixture `tests/fixtures/validation/obara_saika.json` (oracle
//! `scripts/validation/oracle_obara_saika.py`, seed 42): Boys references
//! from the closed form 1/2·x^-(n+1/2)·γ(n+1/2,x) (scipy, no erf
//! approximation) and primitive (ss|ss) references from PySCF `int2e`
//! over single-primitive custom bases (libcint, independent of our
//! recurrence). Units: bohr / bohr^-2 / Hartree.
//!
//! Core module (`sci_form::scf::obara_saika`).

use sci_form::scf::obara_saika::{boys_function, eri_ssss, schwarz_bound, ShellPairData};

#[derive(serde::Deserialize)]
struct BoysCase {
    n: usize,
    x: f64,
    #[serde(rename = "ref")]
    reference: f64,
}

#[derive(serde::Deserialize)]
struct EriCase {
    alpha: f64,
    ca: [f64; 3],
    beta: f64,
    cb: [f64; 3],
    gamma: f64,
    cc: [f64; 3],
    delta: f64,
    cd: [f64; 3],
    #[serde(rename = "ref")]
    reference: f64,
}

#[derive(serde::Deserialize)]
struct SchwarzCase {
    alpha: f64,
    ca: [f64; 3],
    beta: f64,
    cb: [f64; 3],
    bound_ref: f64,
}

#[derive(serde::Deserialize)]
struct Fixture {
    boys_tolerance_rel: f64,
    eri_tolerance_rel: f64,
    boys: Vec<BoysCase>,
    eris: Vec<EriCase>,
    schwarz: Vec<SchwarzCase>,
}

fn load() -> Fixture {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/validation/obara_saika.json"
    );
    serde_json::from_str(&std::fs::read_to_string(path).expect("obara_saika.json"))
        .expect("parse fixture")
}

/// F80a: Boys function vs scipy closed form (n = 0..6, x = 0..100,
/// including the x = 30 asymptotic-branch boundary on both sides).
#[test]
fn f80_boys_matches_scipy() {
    let fix = load();
    assert!(fix.boys.len() >= 100);
    let mut worst = (0.0f64, String::new());
    for c in &fix.boys {
        let got = boys_function(c.n, c.x);
        let denom = c.reference.abs().max(1e-300);
        let rel = (got - c.reference).abs() / denom;
        if rel > worst.0 {
            worst = (rel, format!("n={} x={}", c.n, c.x));
        }
    }
    println!(
        "F80a Boys vs scipy: n={} worst rel err {:.3e} ({})",
        fix.boys.len(),
        worst.0,
        worst.1
    );
    assert!(
        worst.0 <= fix.boys_tolerance_rel,
        "Boys worst rel err {:.3e} at {} (tol {})",
        worst.0,
        worst.1,
        fix.boys_tolerance_rel
    );
}

/// F80b: primitive (ss|ss) vs PySCF/libcint on random exponent/center
/// quartets (seed 42).
#[test]
fn f80_eri_matches_pyscf() {
    let fix = load();
    assert!(fix.eris.len() >= 10);
    let mut worst = (0.0f64, String::new());
    for c in &fix.eris {
        let ab = ShellPairData::new(c.alpha, c.ca, c.beta, c.cb);
        let cd = ShellPairData::new(c.gamma, c.cc, c.delta, c.cd);
        let got = eri_ssss(&ab, &cd);
        let rel = (got - c.reference).abs() / c.reference.abs().max(1e-300);
        if rel > worst.0 {
            worst = (
                rel,
                format!(
                    "a={:.3} b={:.3} g={:.3} d={:.3}",
                    c.alpha, c.beta, c.gamma, c.delta
                ),
            );
        }
    }
    println!(
        "F80b ERI vs PySCF: n={} worst rel err {:.3e} ({})",
        fix.eris.len(),
        worst.0,
        worst.1
    );
    assert!(
        worst.0 <= fix.eri_tolerance_rel,
        "ERI worst rel err {:.3e} at {} (tol {})",
        worst.0,
        worst.1,
        fix.eri_tolerance_rel
    );
}

/// F80c: Schwarz bounds match sqrt((ij|ij)) from libcint, and the bound
/// is never violated: Q_ij·Q_kl >= |(ij|kl)| on a deterministic sweep.
#[test]
fn f80_schwarz_never_violated() {
    let fix = load();
    assert!(fix.schwarz.len() >= 8);
    let mut worst = 0.0f64;
    for c in &fix.schwarz {
        let got = schwarz_bound(c.alpha, c.ca, c.beta, c.cb);
        let rel = (got - c.bound_ref).abs() / c.bound_ref.abs().max(1e-300);
        worst = worst.max(rel);
    }
    println!(
        "F80c Schwarz vs libcint: n={} worst rel err {:.3e}",
        fix.schwarz.len(),
        worst
    );
    assert!(
        worst <= 1e-10,
        "Schwarz worst rel err {worst:.3e} (tol 1e-10)"
    );

    // Deterministic pseudo-random sweep (xorshift64, seed 42): the
    // screening inequality must hold everywhere, with margin.
    let mut state: u64 = 42;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state as f64) / (u64::MAX as f64)
    };
    let mut pairs = Vec::new();
    for _ in 0..12 {
        let a = 10f64.powf(next() * (5.0f64.log10() - 0.2f64.log10()) + 0.2f64.log10());
        let b = 10f64.powf(next() * (5.0f64.log10() - 0.2f64.log10()) + 0.2f64.log10());
        let ca = [next() * 4.0 - 2.0, next() * 4.0 - 2.0, next() * 4.0 - 2.0];
        let cb = [next() * 4.0 - 2.0, next() * 4.0 - 2.0, next() * 4.0 - 2.0];
        pairs.push(ShellPairData::new(a, ca, b, cb));
    }
    let q: Vec<f64> = pairs
        .iter()
        .map(|sp| schwarz_bound(sp.alpha, sp.center_a, sp.beta, sp.center_b))
        .collect();
    let mut min_margin = f64::INFINITY;
    for (i, pi) in pairs.iter().enumerate() {
        for (j, pj) in pairs.iter().enumerate() {
            let direct = eri_ssss(pi, pj).abs();
            let bound = q[i] * q[j];
            // 1e-12 slack: for i == j, sqrt(e)^2 can round 1 ulp below e.
            assert!(
                bound >= direct * (1.0 - 1e-12),
                "Schwarz violated: Q{i}·Q{j} = {bound:.6e} < |(ij|kl)| = {direct:.6e}"
            );
            if direct > 1e-300 {
                min_margin = min_margin.min(bound / direct);
            }
        }
    }
    println!("F80c Schwarz sweep: 12 pairs, 144 bounds, min margin {min_margin:.3}");
    assert!(min_margin >= 1.0 - 1e-12);
}
