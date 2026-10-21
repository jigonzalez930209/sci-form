//! F88: MP meshes (ASE), Bloch phases, 1D chain bands (analytic).
//!
//! Fixture `tests/fixtures/validation/periodic_ref.json` (oracle
//! `scripts/validation/oracle_periodic.py`): ASE `monkhorst_pack` grids,
//! textbook Gamma grids, e^{i2pik.R} phases, E(k) = e + 2t cos(2 pi k)
//! chain bands. Documented gap: `build_bloch_hamiltonian` keeps the cos
//! part only (complex H(k) unimplemented).
//!
//! Core module (`sci_form::kmesh`).

use nalgebra::DMatrix;
use sci_form::kmesh::{
    assemble_periodic_operators, bloch_phase, build_bloch_hamiltonian, bz_integrate_scalar,
    bz_integrate_vector, monkhorst_pack_mesh, validate_kmesh_weights, KMeshCentering, KMeshConfig,
};
use std::collections::HashMap;

#[derive(serde::Deserialize)]
struct MeshCase {
    grid: [usize; 3],
    centering: String,
    ase_ref: Option<Vec<Vec<f64>>>,
    closed_ref_1d: HashMap<String, Vec<f64>>,
}

#[derive(serde::Deserialize)]
struct PhaseCase {
    k: [f64; 3],
    #[serde(rename = "R")]
    r: [i32; 3],
    cos_ref: f64,
    sin_ref: f64,
}

#[derive(serde::Deserialize)]
struct ChainCase {
    eps: f64,
    t: f64,
    k_path: Vec<f64>,
    bands_ref: Vec<f64>,
}

#[derive(serde::Deserialize)]
struct Chain2Case {
    es: f64,
    ts: f64,
    ep: f64,
    tp: f64,
    k_path: Vec<f64>,
    band1_ref: Vec<f64>,
    band2_ref: Vec<f64>,
}

#[derive(serde::Deserialize)]
struct Fixture {
    mesh_tolerance: f64,
    phase_tolerance: f64,
    band_tolerance: f64,
    meshes: Vec<MeshCase>,
    phases: Vec<PhaseCase>,
    chains: Vec<ChainCase>,
    chains2: Vec<Chain2Case>,
}

fn load() -> Fixture {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/validation/periodic_ref.json"
    );
    serde_json::from_str(&std::fs::read_to_string(path).expect("periodic_ref.json"))
        .expect("parse fixture")
}

/// F88a: MP meshes vs ASE (+ closed form); Gamma vs closed form; weights sum to 1.
#[test]
fn f88_meshes_match_ase() {
    let fix = load();
    assert!(fix.meshes.len() >= 8);
    for m in &fix.meshes {
        let centering = match m.centering.as_str() {
            "MonkhorstPack" => KMeshCentering::MonkhorstPack,
            "GammaCentered" => KMeshCentering::GammaCentered,
            other => panic!("centering {other}"),
        };
        let mesh = monkhorst_pack_mesh(&KMeshConfig {
            grid: m.grid,
            centering,
        })
        .expect("mesh");
        assert_eq!(mesh.points.len(), m.grid.iter().product::<usize>());
        assert!(validate_kmesh_weights(&mesh, 1e-12));
        assert!((mesh.total_weight - 1.0).abs() <= 1e-12);
        // Uniform weights.
        let w0 = mesh.points[0].weight;
        assert!(mesh.points.iter().all(|p| (p.weight - w0).abs() <= 1e-15));
        // Compare as sorted sets (ordering convention-independent).
        let mut got: Vec<Vec<f64>> = mesh.points.iter().map(|p| p.fractional.to_vec()).collect();
        for v in got.iter_mut() {
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        }
        got.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if let Some(ase) = &m.ase_ref {
            assert_eq!(got.len(), ase.len());
            let e = got.iter().zip(ase).fold(0.0f64, |mm, (g, w)| {
                mm.max(
                    g.iter()
                        .zip(w)
                        .fold(0.0, |m2, (&a, &b)| m2.max((a - b).abs())),
                )
            });
            assert!(e <= fix.mesh_tolerance, "MP {:?}: {e:.3e}", m.grid);
        }
        // Per-axis closed form.
        for (ax, n) in m.grid.iter().enumerate() {
            let mut axis: Vec<f64> = mesh.points.iter().map(|p| p.fractional[ax]).collect();
            axis.sort_by(|a, b| a.partial_cmp(b).unwrap());
            axis.dedup();
            let want = &m.closed_ref_1d[&n.to_string()];
            assert_eq!(axis.len(), want.len(), "axis count {n}");
            for (g, w) in axis.iter().zip(want) {
                assert!((g - w).abs() <= fix.mesh_tolerance, "axis value {g} vs {w}");
            }
        }
    }
    println!("F88a meshes: {} cases", fix.meshes.len());
}

/// F88b: Bloch phases vs closed form.
#[test]
fn f88_bloch_phases_match() {
    let fix = load();
    assert!(fix.phases.len() >= 10);
    let mut worst = 0.0f64;
    for c in &fix.phases {
        let (cos_p, sin_p) = bloch_phase(c.k, c.r);
        worst = worst
            .max((cos_p - c.cos_ref).abs())
            .max((sin_p - c.sin_ref).abs());
        // Unit modulus.
        assert!((cos_p * cos_p + sin_p * sin_p - 1.0).abs() <= 1e-12);
    }
    println!("F88b phases: worst {worst:.3e}");
    assert!(worst <= fix.phase_tolerance);
}

/// F88c: 1D chain bands vs analytic cosine dispersion (1 + 2 orbitals).
#[test]
fn f88_chain_bands_match_analytic() {
    let fix = load();
    for c in &fix.chains {
        let h0 = DMatrix::from_row_slice(1, 1, &[c.eps]);
        let ht = DMatrix::from_row_slice(1, 1, &[c.t]);
        let space = vec![(&h0, [0, 0, 0]), (&ht, [1, 0, 0]), (&ht, [-1, 0, 0])];
        let mut worst = 0.0f64;
        for (k, want) in c.k_path.iter().zip(&c.bands_ref) {
            let hk = build_bloch_hamiltonian(&space, [*k, 0.0, 0.0]);
            worst = worst.max((hk[(0, 0)] - want).abs());
        }
        println!("F88c chain e={} t={}: worst {worst:.3e}", c.eps, c.t);
        assert!(worst <= fix.band_tolerance);
    }
    for c in &fix.chains2 {
        let h0 = DMatrix::from_row_slice(2, 2, &[c.es, 0.0, 0.0, c.ep]);
        let ht = DMatrix::from_row_slice(2, 2, &[c.ts, 0.0, 0.0, c.tp]);
        let space = vec![(&h0, [0, 0, 0]), (&ht, [1, 0, 0]), (&ht, [-1, 0, 0])];
        for ((k, w1), w2) in c.k_path.iter().zip(&c.band1_ref).zip(&c.band2_ref) {
            let hk = build_bloch_hamiltonian(&space, [*k, 0.0, 0.0]);
            let ev = hk.symmetric_eigen().eigenvalues;
            let (lo, hi) = if ev[0] <= ev[1] {
                (ev[0], ev[1])
            } else {
                (ev[1], ev[0])
            };
            assert!((lo - w1).abs() <= fix.band_tolerance, "band1 {lo} vs {w1}");
            assert!((hi - w2).abs() <= fix.band_tolerance, "band2 {hi} vs {w2}");
        }
    }
}

/// F88d: BZ integration + operator assembly consistency.
#[test]
fn f88_bz_integration_and_assembly() {
    let fix = load();
    let mesh = monkhorst_pack_mesh(&KMeshConfig {
        grid: [2, 2, 2],
        centering: KMeshCentering::MonkhorstPack,
    })
    .expect("mesh");
    let w: Vec<f64> = mesh.points.iter().map(|p| p.weight).collect();
    // Constant integrates to itself; linear k.x integrates to ~0 by symmetry.
    let ones = vec![3.0; 8];
    assert!((bz_integrate_scalar(&ones, &w) - 3.0).abs() <= 1e-12);
    let kx: Vec<f64> = mesh.points.iter().map(|p| p.fractional[0]).collect();
    assert!(bz_integrate_scalar(&kx, &w).abs() <= 1e-12);
    let vecs = vec![vec![1.0, 2.0]; 8];
    assert_eq!(bz_integrate_vector(&vecs, &w), vec![1.0, 2.0]);
    assert!(bz_integrate_vector(&[], &w).is_empty());
    // Assembly: one operator per k-point; Gamma H = sum_R H(R).
    let c = &fix.chains[0];
    let h0 = DMatrix::from_row_slice(1, 1, &[c.eps]);
    let ht = DMatrix::from_row_slice(1, 1, &[c.t]);
    let s0 = DMatrix::from_row_slice(1, 1, &[1.0]);
    let space = vec![(&h0, [0, 0, 0]), (&ht, [1, 0, 0]), (&ht, [-1, 0, 0])];
    let sspace = vec![(&s0, [0, 0, 0])];
    let (hs, ss) = assemble_periodic_operators(&space, &sspace, &mesh);
    assert_eq!(hs.len(), 8);
    assert_eq!(ss.len(), 8);
    let g0 = build_bloch_hamiltonian(&space, [0.0, 0.0, 0.0]);
    assert!((g0[(0, 0)] - (c.eps + 2.0 * c.t)).abs() <= 1e-12);
}
