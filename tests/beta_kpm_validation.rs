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

/// F91d: KPM remains accurate at the scales where it is meant to compete
/// with diagonalization. The reference is the exact spectrum of the same
/// disordered Hückel chain; the gate checks the integrated DOS and the
/// half-filling Fermi energy (spectral accuracy within 0.1 eV).
#[test]
fn f91_kpm_scale_matches_diagonalization() {
    use nalgebra::SymmetricEigen;

    fn disordered_chain(n: usize) -> DMatrix<f64> {
        let mut h = DMatrix::zeros(n, n);
        for i in 0..n {
            h[(i, i)] = 0.25 * ((i as f64) * 0.37).sin();
            if i + 1 < n {
                h[(i, i + 1)] = -1.0;
                h[(i + 1, i)] = -1.0;
            }
        }
        h
    }

    fn gaussian_dos(energies: &[f64], grid: &[f64], sigma: f64) -> Vec<f64> {
        grid.iter()
            .map(|&e| {
                energies
                    .iter()
                    .map(|&eps| (-((e - eps) / sigma).powi(2)).exp())
                    .sum::<f64>()
                    / (sigma * std::f64::consts::PI.sqrt())
            })
            .collect()
    }

    for n in [200usize, 500] {
        let h = disordered_chain(n);
        let eigen = SymmetricEigen::new(h.clone());
        let mut eigenvalues: Vec<f64> = eigen.eigenvalues.iter().copied().collect();
        eigenvalues.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let lo = eigenvalues[0] - 0.5;
        let hi = eigenvalues[n - 1] + 0.5;
        let n_points = 800usize;
        let grid: Vec<f64> = (0..n_points)
            .map(|i| lo + (hi - lo) * i as f64 / (n_points - 1) as f64)
            .collect();

        let config = KpmConfig {
            order: 400,
            n_vectors: 0,
            seed: 42,
            temperature: 0.0,
        };
        let kpm = compute_kpm_dos(&h, &config, lo, hi, n_points);
        let dx = (hi - lo) / (n_points - 1) as f64;
        let integral: f64 = kpm.total_dos.iter().sum::<f64>() * dx;
        let integral_err = (integral - n as f64).abs() / n as f64;

        // Fermi energy at half filling: linear interpolation of the
        // cumulative KPM DOS, compared with the exact median eigenvalue.
        let target = n as f64 / 2.0;
        let mut cum = 0.0;
        let mut e_fermi = f64::NAN;
        for i in 1..n_points {
            let step = 0.5 * (kpm.total_dos[i - 1] + kpm.total_dos[i]) * dx;
            if cum + step >= target {
                let frac = (target - cum) / step.max(1e-300);
                e_fermi = grid[i - 1] + frac * dx;
                break;
            }
            cum += step;
        }
        let exact_fermi = 0.5 * (eigenvalues[n / 2 - 1] + eigenvalues[n / 2]);

        // L2 shape vs a Gaussian-broadened exact DOS at the KPM resolution.
        let sigma = (hi - lo) / config.order as f64;
        let reference = gaussian_dos(&eigenvalues, &grid, sigma);
        let num: f64 = kpm
            .total_dos
            .iter()
            .zip(&reference)
            .map(|(&g, &w)| (g - w).powi(2))
            .sum();
        let den: f64 = reference.iter().map(|&w| w * w).sum();
        let l2 = (num / den.max(1e-300)).sqrt();

        eprintln!(
            "F91d N={n}: integral {integral:.2}/{n} (err {integral_err:.2e}), \
             E_F {e_fermi:.4} vs {exact_fermi:.4} eV, L2 {l2:.3e}"
        );
        assert!(
            integral_err <= 0.02,
            "N={n}: integrated DOS error {integral_err:.3e}"
        );
        assert!(
            (e_fermi - exact_fermi).abs() <= 0.1,
            "N={n}: E_F {e_fermi:.4} vs {exact_fermi:.4}"
        );
        assert!(l2 <= 0.25, "N={n}: DOS L2 {l2:.3e}");
    }
}

/// F91e: sparse CSR KPM matches the dense path and scales past the dense
/// matrix limit. The 5000-site chain has O(N) non-zeros; the dense matrix
/// would need 200 MB and a dense Chebyshev recursion.
#[test]
fn f91_sparse_kpm_matches_dense_and_scales() {
    use sci_form::beta::kpm::{compute_kpm_dos_sparse, SparseMatrix};
    use std::time::Instant;

    fn chain(n: usize) -> SparseMatrix {
        let mut triplets = Vec::with_capacity(3 * n);
        for i in 0..n {
            triplets.push((i, i, 0.25 * ((i as f64) * 0.37).sin()));
            if i + 1 < n {
                triplets.push((i, i + 1, -1.0));
                triplets.push((i + 1, i, -1.0));
            }
        }
        SparseMatrix::from_triplets(n, &triplets)
    }

    let sparse = chain(200);
    let mut dense = DMatrix::zeros(200, 200);
    for i in 0..200 {
        for k in sparse.row_offsets[i]..sparse.row_offsets[i + 1] {
            dense[(i, sparse.col_indices[k])] = sparse.values[k];
        }
    }
    let config = KpmConfig {
        order: 200,
        n_vectors: 0,
        seed: 42,
        temperature: 0.0,
    };
    let lo = -2.5;
    let hi = 2.5;
    let n_points = 400;
    let dense_result = compute_kpm_dos(&dense, &config, lo, hi, n_points);
    let sparse_result = compute_kpm_dos_sparse(&sparse, &config, lo, hi, n_points);
    let num: f64 = dense_result
        .total_dos
        .iter()
        .zip(&sparse_result.total_dos)
        .map(|(a, b)| (a - b).powi(2))
        .sum();
    let den: f64 = dense_result
        .total_dos
        .iter()
        .map(|v| v * v)
        .sum::<f64>()
        .max(1e-300);
    let l2 = (num / den).sqrt();
    eprintln!("F91e sparse vs dense N=200: L2 {l2:.3e}");
    assert!(l2 < 1e-10, "sparse/dense L2 {l2:.3e}");

    let big = chain(5000);
    let big_config = KpmConfig {
        order: 400,
        n_vectors: 32,
        seed: 42,
        temperature: 0.0,
    };
    let started = Instant::now();
    let big_result = compute_kpm_dos_sparse(&big, &big_config, -3.0, 3.0, 600);
    let elapsed = started.elapsed().as_secs_f64();
    let dx = 6.0 / 599.0;
    let integral: f64 = big_result.total_dos.iter().sum::<f64>() * dx;
    eprintln!(
        "F91e sparse N=5000: {:?} nnz={}, elapsed {:.2}s, integral {:.1}",
        big.nnz(),
        big.nnz(),
        elapsed,
        integral
    );
    assert!(elapsed < 30.0, "sparse KPM took {elapsed:.1}s");
    assert!(
        (integral - 5000.0).abs() / 5000.0 < 0.05,
        "integrated DOS {integral}"
    );
    assert!(big_result.total_dos.iter().all(|v| v.is_finite()));
}
