//! F86: EDL Gouy-Chapman / GCS / Helmholtz vs analytic electrochemistry.
//!
//! Fixture `tests/fixtures/validation/edl_analytic.json` (oracle
//! `scripts/validation/oracle_edl.py`): independent numpy implementation
//! of the textbook closed forms — GC potential/C_diff, Grahame-inverted
//! Stern split, Helmholtz compact capacitor, Stern series capacitance.
//! Same equations, independently coded constants and unit handling, so
//! the gate catches wrong constants, Å↔m slips, and missing factors.

#![cfg(feature = "alpha-edl")]

use sci_form::alpha::edl::{
    compute_edl_profile, compute_gcs_profile, compute_gouy_chapman_profile,
    compute_helmholtz_profile, scan_edl_capacitance, EdlConfig, EdlModel,
};

#[derive(serde::Deserialize)]
struct GcCase {
    #[serde(rename = "T_K")]
    t: f64,
    #[serde(rename = "I_M")]
    ionic: f64,
    eps_b: f64,
    #[serde(rename = "phi0_V")]
    phi0: f64,
    n: usize,
    #[serde(rename = "extent_A")]
    extent: f64,
    #[serde(rename = "x_A")]
    x: Vec<f64>,
    #[serde(rename = "phi_V")]
    phi: Vec<f64>,
    #[serde(rename = "field_V_per_m")]
    field: Vec<f64>,
    #[serde(rename = "rho_C_per_m3")]
    rho: Vec<f64>,
    #[serde(rename = "c_diffuse_F_per_m2")]
    c_diff: f64,
}

#[derive(serde::Deserialize)]
struct GcsCase {
    #[serde(rename = "T_K")]
    t: f64,
    #[serde(rename = "I_M")]
    ionic: f64,
    eps_b: f64,
    eps_c: f64,
    #[serde(rename = "stern_A")]
    stern: f64,
    #[serde(rename = "sigma_C_per_m2")]
    sigma: f64,
    n: usize,
    #[serde(rename = "extent_A")]
    extent: f64,
    #[serde(rename = "x_A")]
    x: Vec<f64>,
    #[serde(rename = "phi_V")]
    phi: Vec<f64>,
    #[serde(rename = "field_V_per_m")]
    field: Vec<f64>,
    #[serde(rename = "rho_C_per_m3")]
    rho: Vec<f64>,
    #[serde(rename = "phi_compact_V")]
    phi_c: f64,
    #[serde(rename = "phi_diffuse_V")]
    phi_d: f64,
    #[serde(rename = "c_compact_F_per_m2")]
    c_c: f64,
    #[serde(rename = "c_diffuse_F_per_m2")]
    c_d: f64,
    #[serde(rename = "c_total_F_per_m2")]
    c_t: f64,
}

#[derive(serde::Deserialize)]
struct HCase {
    #[serde(rename = "sigma_C_per_m2")]
    sigma: f64,
    eps_c: f64,
    #[serde(rename = "stern_A")]
    stern: f64,
    #[serde(rename = "c_compact_F_per_m2")]
    c_c: f64,
    #[serde(rename = "drop_V")]
    drop: f64,
}

#[derive(serde::Deserialize)]
struct ScanCase {
    #[serde(rename = "T_K")]
    t: f64,
    #[serde(rename = "I_M")]
    ionic: f64,
    eps_b: f64,
    eps_c: f64,
    #[serde(rename = "stern_A")]
    stern: f64,
    sigma_grid: Vec<f64>,
    #[serde(rename = "c_total_F_per_m2")]
    c_tot: Vec<f64>,
}

#[derive(serde::Deserialize)]
struct Fixture {
    tolerance_rel: f64,
    gc: Vec<GcCase>,
    gcs: Vec<GcsCase>,
    helmholtz: Vec<HCase>,
    scans: Vec<ScanCase>,
}

fn load() -> Fixture {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/validation/edl_analytic.json"
    );
    serde_json::from_str(&std::fs::read_to_string(path).expect("edl_analytic.json"))
        .expect("parse fixture")
}

fn base_config(t: f64, ionic: f64, eps_b: f64, n: usize, extent: f64) -> EdlConfig {
    EdlConfig {
        model: EdlModel::GouyChapman,
        temperature_k: t,
        ionic_strength_m: ionic,
        stern_thickness_angstrom: 3.0,
        compact_dielectric: 6.0,
        bulk_dielectric: eps_b,
        ..EdlConfig::default()
    }
}

/// Normalized max error: max|a-b| / max|b| (robust for decaying tails).
fn nerr(got: &[f64], want: &[f64]) -> f64 {
    assert_eq!(got.len(), want.len());
    let scale = want.iter().fold(0.0f64, |m, &v| m.max(v.abs())).max(1e-300);
    got.iter()
        .zip(want)
        .fold(0.0f64, |m, (&g, &w)| m.max((g - w).abs() / scale))
}

/// F86a: Gouy-Chapman profiles (φ, E, ρ) + diffuse capacitance, 3 T × 3 I × 6 φ₀.
#[test]
fn f86_gc_matches_analytic() {
    let fix = load();
    assert!(fix.gc.len() >= 50);
    let mut worst = (0.0f64, String::new());
    for c in &fix.gc {
        let mut cfg = base_config(c.t, c.ionic, c.eps_b, c.n, c.extent);
        cfg.numerics.n_points = c.n;
        cfg.numerics.extent_angstrom = c.extent;
        let r = compute_gouy_chapman_profile(c.phi0, &cfg).expect("gc profile");
        let e = nerr(&r.distance_axis_angstrom, &c.x)
            .max(nerr(&r.electrostatic_potential_v, &c.phi))
            .max(nerr(&r.field_strength_v_per_m, &c.field))
            .max(nerr(&r.charge_density_c_per_m3, &c.rho));
        if e > worst.0 {
            worst = (e, format!("T={} I={} phi0={}", c.t, c.ionic, c.phi0));
        }
        let dc = (r.differential_capacitance.diffuse_f_per_m2 - c.c_diff).abs()
            / c.c_diff.abs().max(1e-300);
        assert!(dc <= fix.tolerance_rel, "GC C_diff rel err {dc:.3e}");
        // Surface value exact: φ(0) = φ₀.
        assert!((r.electrostatic_potential_v[0] - c.phi0).abs() <= 1e-9);
    }
    println!("F86a GC vs analytic: n={} worst {:.3e}", fix.gc.len(), worst.0);
    assert!(worst.0 <= fix.tolerance_rel, "GC worst {:.3e} at {}", worst.0, worst.1);
}

/// F86b: GCS profiles + Stern split + series capacitance.
#[test]
fn f86_gcs_matches_analytic() {
    let fix = load();
    assert!(fix.gcs.len() >= 50);
    let mut worst = (0.0f64, String::new());
    for c in &fix.gcs {
        let mut cfg = base_config(c.t, c.ionic, c.eps_b, c.n, c.extent);
        cfg.model = EdlModel::GouyChapmanStern;
        cfg.compact_dielectric = c.eps_c;
        cfg.stern_thickness_angstrom = c.stern;
        cfg.numerics.n_points = c.n;
        cfg.numerics.extent_angstrom = c.extent;
        let r = compute_gcs_profile(c.sigma, &cfg).expect("gcs profile");
        let e = nerr(&r.electrostatic_potential_v, &c.phi)
            .max(nerr(&r.field_strength_v_per_m, &c.field))
            .max(nerr(&r.charge_density_c_per_m3, &c.rho));
        if e > worst.0 {
            worst = (e, format!("T={} I={} sigma={}", c.t, c.ionic, c.sigma));
        }
        for (got, want, name) in [
            (r.compact_layer_drop_v, c.phi_c, "phi_c"),
            (r.diffuse_layer_drop_v, c.phi_d, "phi_d"),
            (r.differential_capacitance.compact_f_per_m2, c.c_c, "C_H"),
            (r.differential_capacitance.diffuse_f_per_m2, c.c_d, "C_diff"),
            (r.differential_capacitance.total_f_per_m2, c.c_t, "C_tot"),
        ] {
            let d = (got - want).abs() / want.abs().max(1e-300);
            assert!(d <= fix.tolerance_rel, "GCS {name} rel err {d:.3e}");
        }
    }
    println!("F86b GCS vs analytic: n={} worst {:.3e}", fix.gcs.len(), worst.0);
    assert!(worst.0 <= fix.tolerance_rel, "GCS worst {:.3e} at {}", worst.0, worst.1);
}

/// F86c: Helmholtz compact capacitor (C_H + drop linear in σ).
#[test]
fn f86_helmholtz_matches_analytic() {
    let fix = load();
    assert!(!fix.helmholtz.is_empty());
    let mut worst = 0.0f64;
    for c in &fix.helmholtz {
        let mut cfg = EdlConfig::default();
        cfg.compact_dielectric = c.eps_c;
        cfg.stern_thickness_angstrom = c.stern;
        let r = compute_helmholtz_profile(c.sigma, &cfg).expect("helmholtz");
        let dc = (r.differential_capacitance.compact_f_per_m2 - c.c_c).abs() / c.c_c;
        let dd = (r.compact_layer_drop_v - c.drop).abs() / c.drop.abs().max(1e-300);
        worst = worst.max(dc).max(dd);
    }
    println!("F86c Helmholtz vs analytic: n={} worst {:.3e}", fix.helmholtz.len(), worst);
    assert!(worst <= fix.tolerance_rel);
}

/// F86d: capacitance scan C(σ) matches reference curve; even with minimum at σ = 0.
#[test]
fn f86_capacitance_scan_shape() {
    let fix = load();
    assert!(!fix.scans.is_empty());
    for c in &fix.scans {
        let mut cfg = base_config(c.t, c.ionic, c.eps_b, 65, 50.0);
        cfg.model = EdlModel::GouyChapmanStern;
        cfg.compact_dielectric = c.eps_c;
        cfg.stern_thickness_angstrom = c.stern;
        let lo = c.sigma_grid[0];
        let hi = *c.sigma_grid.last().unwrap();
        let scan = scan_edl_capacitance(lo, hi, c.sigma_grid.len(), &cfg).expect("scan");
        assert_eq!(scan.len(), c.c_tot.len());
        let scale = c.c_tot.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
        let e = scan
            .iter()
            .zip(&c.c_tot)
            .fold(0.0f64, |m, ((_, g), &w)| m.max((g - w).abs() / scale));
        println!("F86d scan T={} I={}: worst {:.3e}", c.t, c.ionic, e);
        assert!(e <= fix.tolerance_rel);
        // Physics shape: even curve, minimum at zero charge.
        let mid = scan.len() / 2;
        for i in 0..mid {
            let (sl, cl) = scan[mid - 1 - i];
            let (sr, cr) = scan[mid + 1 + i];
            assert!((sl + sr).abs() <= 1e-12, "scan grid not symmetric");
            assert!((cl - cr).abs() / cl <= 1e-9, "C(σ) not even");
        }
        let c0 = scan[mid].1;
        for (_, ci) in &scan {
            assert!(*ci >= c0 * (1.0 - 1e-9), "C minimum not at σ=0");
        }
    }
    // compute_edl_profile dispatches to the same backends.
    let cfg = EdlConfig::default();
    let _ = compute_edl_profile(0.05, &cfg).expect("dispatch");
}
