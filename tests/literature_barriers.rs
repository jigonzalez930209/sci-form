//! Literature barrier validation (TODO item 179).
//!
//! Compares computed barriers against experimental gas-phase references
//! (`tests/fixtures/validation/literature_barriers.json`) using rigid scans
//! with exact-symmetry endpoints through the validated single-point
//! backends (`neb_backend_energy_kcal` with `uff`/`gfn2`).
//!
//! Design notes (all verified by measurement, 2026-10-13):
//! - Torsions rotate one top rigidly (0..120 deg); the 60-deg eclipsed TS is
//!   always sampled. Ammonia mirrors N through the H3 plane (planar TS
//!   sampled). Rigid profiles are upper bounds to the relaxed MEP.
//! - The simplified Cartesian NEB is NOT used: linear interpolation cuts
//!   chords through rotating tops (C-H compressed to ~cos(60 deg) mid-path,
//!   energies diverge to thousands of kcal/mol) and over-relaxation collapses
//!   the path (ammonia 9x100 -> 0.44 kcal/mol). The experimental CI-NEB
//!   diverged on these paths (UFF energies ~1e78, GFN2 ethane collapsed to
//!   0.24). All documented in the fixture; numbers come from scans.
//! - The `pm3` NEB backend is excluded: it runs the unvalidated legacy PM3
//!   engine (`solve_pm3_with_state`), not the faithful F42 port (e.g.
//!   methanol rigid gives 29.4, ammonia inverts).

use serde::Deserialize;

#[derive(Deserialize)]
struct BarrierCase {
    name: String,
    smiles: String,
    literature_barrier_kcal_mol: f64,
    backends: Vec<String>,
}

#[derive(Deserialize)]
struct BarrierFixture {
    tolerance_kcal_mol: f64,
    cases: Vec<BarrierCase>,
}

fn load_fixture() -> BarrierFixture {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/validation/literature_barriers.json"
    );
    let text = std::fs::read_to_string(path).expect("literature_barriers.json");
    serde_json::from_str(&text).expect("parse fixture")
}

fn get_pos(c: &[f64], i: usize) -> [f64; 3] {
    [c[3 * i], c[3 * i + 1], c[3 * i + 2]]
}

fn set_pos(c: &mut [f64], i: usize, p: [f64; 3]) {
    c[3 * i] = p[0];
    c[3 * i + 1] = p[1];
    c[3 * i + 2] = p[2];
}

/// Rigid rotation of point `p` about the axis through `pt` by `theta`.
fn rot_about_axis(p: [f64; 3], ax: [f64; 3], pt: [f64; 3], theta: f64) -> [f64; 3] {
    let n = (ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2]).sqrt();
    let u = [ax[0] / n, ax[1] / n, ax[2] / n];
    let v = [p[0] - pt[0], p[1] - pt[1], p[2] - pt[2]];
    let (c, s) = (theta.cos(), theta.sin());
    let dot = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    let cr = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    [
        pt[0] + v[0] * c + cr[0] * s + u[0] * dot * (1.0 - c),
        pt[1] + v[1] * c + cr[1] * s + u[1] * dot * (1.0 - c),
        pt[2] + v[2] * c + cr[2] * s + u[2] * dot * (1.0 - c),
    ]
}

/// Barrier from a rigid torsion scan (methyl top `rot_atoms` about the
/// bond `axis_a`-`axis_b`, 0..120 deg in 20-deg steps).
fn torsion_barrier(
    smiles: &str,
    method: &str,
    start: &[f64],
    axis_a: usize,
    axis_b: usize,
    rot_atoms: std::ops::Range<usize>,
    pivot: usize,
) -> f64 {
    let pa = get_pos(start, axis_a);
    let pb = get_pos(start, axis_b);
    let ax = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
    let pp = get_pos(start, pivot);
    let mut es = Vec::new();
    for d in 0..=6 {
        let mut g = start.to_vec();
        for i in rot_atoms.clone() {
            set_pos(
                &mut g,
                i,
                rot_about_axis(
                    get_pos(start, i),
                    ax,
                    pp,
                    d as f64 * std::f64::consts::PI / 9.0,
                ),
            );
        }
        es.push(sci_form::neb_backend_energy_kcal(method, smiles, &g).expect("single point"));
    }
    let mx = es.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    mx - 0.5 * (es[0] + es[6])
}

/// Barrier from the ammonia umbrella scan (N between exact mirrors).
fn ammonia_barrier(method: &str, start: &[f64]) -> f64 {
    let (h0, h1, h2) = (get_pos(start, 1), get_pos(start, 2), get_pos(start, 3));
    let u = [h1[0] - h0[0], h1[1] - h0[1], h1[2] - h0[2]];
    let v = [h2[0] - h0[0], h2[1] - h0[1], h2[2] - h0[2]];
    let mut n = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    let nl = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    n = [n[0] / nl, n[1] / nl, n[2] / nl];
    let nat = get_pos(start, 0);
    let d0 = (nat[0] - h0[0]) * n[0] + (nat[1] - h0[1]) * n[1] + (nat[2] - h0[2]) * n[2];
    let mut es = Vec::new();
    for k in 0..=6 {
        let dd = d0 * (1.0 - 2.0 * k as f64 / 6.0);
        let mut g = start.to_vec();
        set_pos(
            &mut g,
            0,
            [
                nat[0] + n[0] * (dd - d0),
                nat[1] + n[1] * (dd - d0),
                nat[2] + n[2] * (dd - d0),
            ],
        );
        es.push(sci_form::neb_backend_energy_kcal(method, "N", &g).expect("single point"));
    }
    let mx = es.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    mx - 0.5 * (es[0] + es[6])
}

#[test]
fn barriers_match_literature() {
    let fix = load_fixture();
    // Seed-42 embedding atom order (verified): CC -> [C0, C1, H2..4, H5..7],
    // CO -> [C0, O1, H2..4 (methyl), H5 (hydroxyl)], N -> [N0, H1..3].
    let mut worst = (0.0f64, String::new());
    for case in &fix.cases {
        let conf = sci_form::embed(&case.smiles, 42);
        assert!(conf.error.is_none(), "{}: embed failed", case.name);
        for backend in &case.backends {
            let b = match case.name.as_str() {
                "ethane_torsion" => {
                    torsion_barrier(&case.smiles, backend, &conf.coords, 0, 1, 5..8, 1)
                }
                "methanol_torsion" => {
                    torsion_barrier(&case.smiles, backend, &conf.coords, 0, 1, 2..5, 0)
                }
                "ammonia_inversion" => ammonia_barrier(backend, &conf.coords),
                other => panic!("unknown case {other}"),
            };
            let d = (b - case.literature_barrier_kcal_mol).abs();
            eprintln!(
                "LIT {:17} {:4}: barrier={:7.3} lit={:6.2} |d|={:.3}",
                case.name, backend, b, case.literature_barrier_kcal_mol, d
            );
            if d > worst.0 {
                worst = (d, format!("{}:{}", case.name, backend));
            }
        }
    }
    assert!(
        worst.0 <= fix.tolerance_kcal_mol,
        "worst literature deviation {:.3} kcal/mol ({}) exceeds {:.1}",
        worst.0,
        worst.1,
        fix.tolerance_kcal_mol
    );
}
