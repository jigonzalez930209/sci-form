//! F83: SDR projections + embedding vs numpy + CVXPY.
//!
//! Fixture `tests/fixtures/validation/sdr_ref.json` (oracle
//! `scripts/validation/oracle_sdr.py`, seed 42): numpy double-centering,
//! eigh-clip PSD, eigenvalue soft-threshold, top-3 Gram reconstruction,
//! and CVXPY/SCS convex-SDP optima (max distance error) on exact and
//! noisy distance sets. Coordinates are compared via rotation-invariant
//! Gram reconstruction (eigenvector signs are arbitrary).
//!
//! Core module (`sci_form::sdr`).

use nalgebra::DMatrix;
use sci_form::sdr::{alternating_projections, extract_coordinates, project_distances, project_psd, svt_step, warm_start_gram, SdrConfig};

#[derive(serde::Deserialize)]
struct ProjCase {
    n: usize,
    x: Vec<Vec<f64>>,
    tau: f64,
    psd_ref: Vec<Vec<f64>>,
    neg_ref: usize,
    svt_ref: Vec<Vec<f64>>,
}

#[derive(serde::Deserialize)]
struct EmbedCase {
    n: usize,
    noise: f64,
    pairs: Vec<(usize, usize, f64)>,
    gram_ref: Vec<Vec<f64>>,
    gram_recon_ref: Vec<Vec<f64>>,
    cvxpy_max_err: f64,
}

#[derive(serde::Deserialize)]
struct Fixture {
    projection_tolerance: f64,
    gram_tolerance: f64,
    embed_tolerance_A: f64,
    projections: Vec<ProjCase>,
    embeddings: Vec<EmbedCase>,
}

fn load() -> Fixture {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/validation/sdr_ref.json"
    );
    serde_json::from_str(&std::fs::read_to_string(path).expect("sdr_ref.json"))
        .expect("parse fixture")
}

fn mat(n: usize, rows: &[Vec<f64>]) -> DMatrix<f64> {
    DMatrix::from_row_slice(n, n, &rows.concat())
}

fn frob_rel(got: &DMatrix<f64>, want: &DMatrix<f64>) -> f64 {
    let num: f64 = (got - want).iter().map(|&v| v * v).sum();
    let den: f64 = want.iter().map(|&v| v * v).sum::<f64>().max(1e-300);
    (num / den).sqrt()
}

/// F83a: PSD clip (+ removal count) and SVT vs numpy.
#[test]
fn f83_projections_match_numpy() {
    let fix = load();
    assert!(fix.projections.len() >= 2);
    for c in &fix.projections {
        let x = mat(c.n, &c.x);
        let (p, neg) = project_psd(&x);
        assert_eq!(neg, c.neg_ref, "removal count");
        let e = frob_rel(&p, &mat(c.n, &c.psd_ref));
        assert!(e <= fix.projection_tolerance, "psd {e:.3e}");
        // Idempotence + output PSD (re-project == itself).
        let (pp, neg2) = project_psd(&p);
        assert_eq!(neg2, 0);
        assert!(frob_rel(&pp, &p) <= 1e-12);
        let s = svt_step(&x, c.tau);
        let es = frob_rel(&s, &mat(c.n, &c.svt_ref));
        assert!(es <= fix.projection_tolerance, "svt {es:.3e}");
        // tau <= 0 is the identity.
        assert_eq!(svt_step(&x, 0.0), x);
    }
    println!("F83a projections: {} cases", fix.projections.len());
}

/// F83b: distance projection satisfies the constraints it is given.
#[test]
fn f83_distance_projection_satisfies() {
    let fix = load();
    let c = &fix.embeddings[0];
    let x = mat(c.n, &c.gram_ref);
    let p = project_distances(&x, &c.pairs);
    let mut worst = 0.0f64;
    for &(i, j, d) in &c.pairs {
        let got = (p[(i, i)] + p[(j, j)] - 2.0 * p[(i, j)]).max(0.0).sqrt();
        worst = worst.max((got - d).abs());
    }
    println!("F83b distance projection residual: {worst:.3e}");
    assert!(worst <= 1e-9);
}

/// F83c: full pipeline — warm start, Gram reconstruction, and distance
/// errors vs numpy + CVXPY optimum.
#[test]
fn f83_embed_matches_cvxpy() {
    let fix = load();
    assert!(fix.embeddings.len() >= 3);
    let cfg = SdrConfig { max_iter: 500, tol: 1e-9, ..SdrConfig::default() };
    for c in &fix.embeddings {
        let g0 = warm_start_gram(c.n, &c.pairs);
        let eg = frob_rel(&g0, &mat(c.n, &c.gram_ref));
        assert!(eg <= fix.gram_tolerance, "warm start {eg:.3e}");
        let (gram, conv) = alternating_projections(&g0, &c.pairs, &cfg);
        let coords = extract_coordinates(&gram);
        assert_eq!(coords.len(), 3 * c.n);
        // Rotation-invariant check: reconstructed Gram vs reference.
        let n = c.n;
        let cc: DMatrix<f64> = DMatrix::from_fn(n, 3, |i, j| coords[3 * i + j]);
        let recon = &cc * cc.transpose();
        let er = frob_rel(&recon, &mat(n, &c.gram_recon_ref));
        // Distance errors from our coords.
        let mut max_err = 0.0f64;
        for &(i, j, d) in &c.pairs {
            let dx = coords[3 * i] - coords[3 * j];
            let dy = coords[3 * i + 1] - coords[3 * j + 1];
            let dz = coords[3 * i + 2] - coords[3 * j + 2];
            max_err = max_err.max(((dx * dx + dy * dy + dz * dz).sqrt() - d).abs());
        }
        println!(
            "F83c n={} noise={}: recon {er:.3e} max_err {max_err:.3e} (cvxpy {:.3e}) iters={} conv={}",
            c.n, c.noise, c.cvxpy_max_err, conv.iterations, conv.converged
        );
        assert!(max_err <= fix.embed_tolerance_A.max(c.cvxpy_max_err * 3.0),
            "n={}: {max_err:.3e} vs cvxpy {:.3e}", c.n, c.cvxpy_max_err);
    }
}
