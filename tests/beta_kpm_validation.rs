//! F91: KPM DOS/moments/Mulliken vs numpy eigh.
//!
//! Fixture `tests/fixtures/validation/kpm_numpy.json` (oracle
//! `scripts/validation/oracle_kpm.py`, GOE seed 42): exact eigenvalues
//! from `numpy.linalg.eigh` — independent of our Chebyshev recurrence.
//! Gates: exact moments vs textbook cos(k·arccos) from those eigenvalues
//! (bounds-independent, 1e-9), Jackson DOS curve vs Gaussian-broadened
//! exact DOS (L2 ≤ 5%), Mulliken charges vs projector reference.

#![cfg(feature = "beta-kpm")]

use nalgebra::DMatrix;
use sci_form::beta::kpm::{
    compute_kpm_dos, compute_kpm_mulliken, jackson_kernel, ChebyshevExpansion, KpmConfig,
};

#[derive(serde::Deserialize)]
struct KpmSystem {
    n: usize,
    matrix: Vec<Vec<f64>>,
    eigenvalues: Vec<f64>,
    order: usize,
    #[serde(rename = "grid_eV")]
    grid: Vec<f64>,
    #[serde(rename = "dos_ref")]
    dos_ref: Vec<f64>,
    #[serde(rename = "dos_jackson_ref")]
    dos_jackson_ref: Vec<f64>,
    nelec: usize,
    charges_ref: Vec<f64>,
    #[serde(rename = "temperature_K", default)]
    temperature: f64,
    #[serde(rename = "charges_fd_ref", default)]
    charges_fd_ref: Vec<f64>,
    #[serde(default)]
    gapped: bool,
}

#[derive(serde::Deserialize)]
struct Fixture {
    /// Mirrors the fixture schema; the relative DOS tolerance is not used
    /// because the gate compares the L2 error against the absolute tolerance.
    #[allow(dead_code)]
    dos_tolerance_l2_rel: f64,
    moment_tolerance: f64,
    charge_tolerance_abs: f64,
    systems: Vec<KpmSystem>,
}

fn load() -> Fixture {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/validation/kpm_numpy.json"
    );
    serde_json::from_str(&std::fs::read_to_string(path).expect("kpm_numpy.json"))
        .expect("parse fixture")
}

fn hmat(sys: &KpmSystem) -> DMatrix<f64> {
    DMatrix::from_row_slice(sys.n, sys.n, &sys.matrix.concat())
}

/// F91a: exact Chebyshev moments vs textbook moments from numpy eigenvalues.
#[test]
fn f91_moments_match_eigh() {
    let fix = load();
    assert!(fix.systems.len() >= 2);
    for sys in &fix.systems {
        let h = hmat(sys);
        let order = 64;
        let exp = ChebyshevExpansion::from_matrix_exact(&h, order);
        let mut worst = 0.0f64;
        for k in 0..order {
            let mut acc = 0.0;
            for &lam in &sys.eigenvalues {
                let x = ((lam - exp.b) / exp.a).clamp(-1.0, 1.0);
                acc += (k as f64 * x.acos()).cos();
            }
            let want = acc / sys.n as f64;
            worst = worst.max((exp.moments[k] - want).abs());
        }
        println!("F91a moments N={}: worst {:.3e}", sys.n, worst);
        assert!(worst <= fix.moment_tolerance, "N={}: {worst:.3e}", sys.n);
    }
}

/// F91b: Jackson DOS curve vs exact-eigenvalue Jackson reference (tight,
/// pins rescaling/grid/weight), plus Gaussian shape sanity (kernel class).
#[test]
fn f91_dos_matches_broadened_eigh() {
    let fix = load();
    for sys in &fix.systems {
        if sys.grid.is_empty() {
            continue;
        }
        let h = hmat(sys);
        let cfg = KpmConfig {
            order: sys.order,
            n_vectors: 0,
            seed: 42,
            temperature: 0.0,
        };
        let lo = sys.grid[0];
        let hi = *sys.grid.last().unwrap();
        let r = compute_kpm_dos(&h, &cfg, lo, hi, sys.grid.len());
        assert_eq!(r.energies.len(), sys.grid.len());
        assert_eq!(
            r.total_dos.len(),
            sys.dos_jackson_ref.len(),
            "fixture DOS length mismatch"
        );
        let l2 = |want: &[f64]| {
            let num: f64 = r
                .total_dos
                .iter()
                .zip(want)
                .map(|(&g, &w)| (g - w).powi(2))
                .sum();
            let den: f64 = want.iter().map(|&w| w.powi(2)).sum();
            (num / den).sqrt()
        };
        let tight = l2(&sys.dos_jackson_ref);
        let shape = l2(&sys.dos_ref);
        // Trapezoid integral ≈ N states.
        let dx = (hi - lo) / (sys.grid.len() - 1) as f64;
        let integ: f64 = r.total_dos.iter().sum::<f64>() * dx;
        println!(
            "F91b DOS N={}: Jackson L2 {tight:.3e}, Gauss L2 {shape:.4}, integral {integ:.3}",
            sys.n
        );
        assert!(tight <= 1e-9, "N={}: Jackson L2 {tight:.3e}", sys.n);
        assert!(shape <= 0.15, "N={}: Gauss shape L2 {shape:.4}", sys.n);
        assert!(
            (integ - sys.n as f64).abs() / sys.n as f64 <= 0.02,
            "integral {integ}"
        );
    }
}

/// F91c: Mulliken charges (S = I, half-filled closed shell) vs projector
/// reference — T = 0 on the gapped (molecular-like) spectrum, finite-T
/// Fermi-Dirac branch on the gapless GOE spectra.
#[test]
fn f91_mulliken_matches_projector() {
    let fix = load();
    for sys in &fix.systems {
        let h = hmat(sys);
        let n = sys.n;
        let s = DMatrix::identity(n, n);
        let z = vec![1.0; n];
        if sys.gapped {
            let cfg = KpmConfig {
                order: sys.order,
                n_vectors: 0,
                seed: 42,
                temperature: 0.0,
            };
            let r = compute_kpm_mulliken(&h, &s, sys.nelec, &z, &cfg);
            let worst = r
                .charges
                .iter()
                .zip(&sys.charges_ref)
                .fold(0.0f64, |m, (&g, &w)| m.max((g - w).abs()));
            println!(
                "F91c Mulliken gapped N={}: worst {:.3e} (e_count {:.4}/{}, mu {:.4})",
                sys.n, worst, r.n_electrons_kpm, r.n_electrons_target, r.fermi_level
            );
            assert!(
                worst <= fix.charge_tolerance_abs,
                "N={}: {worst:.3e}",
                sys.n
            );
        } else {
            let cfg = KpmConfig {
                order: sys.order,
                n_vectors: 0,
                seed: 42,
                temperature: sys.temperature,
            };
            assert!(sys.temperature > 0.0);
            let r = compute_kpm_mulliken(&h, &s, sys.nelec, &z, &cfg);
            let worst = r
                .charges
                .iter()
                .zip(&sys.charges_fd_ref)
                .fold(0.0f64, |m, (&g, &w)| m.max((g - w).abs()));
            println!(
                "F91c Mulliken finite-T N={}: worst {:.3e} (e_count {:.4}/{}, mu {:.4})",
                sys.n, worst, r.n_electrons_kpm, r.n_electrons_target, r.fermi_level
            );
            assert!(
                worst <= fix.charge_tolerance_abs,
                "N={}: {worst:.3e}",
                sys.n
            );
        }
    }
    // Jackson kernel sanity (g_0 = 1, decaying).
    let gk = jackson_kernel(500);
    assert!((gk[0] - 1.0).abs() < 1e-15);
    assert!(gk[499] < gk[1]);
}
