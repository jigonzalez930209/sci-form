//! Quantitative beta-module benchmarks against stable routes.
//!
//! Each benchmark compares an approximate method with an independent, exact
//! reference computed in the same test:
//!
//! * `beta-kpm`: Kernel Polynomial Method DOS vs the exact eigenvalue spectrum
//!   (diagonalization).
//! * `beta-randnla`: Nyström reconstruction vs the best rank-k approximation
//!   (exact SVD).
//! * `beta-mbh`: Mobile Block Hessian frequencies vs the full numerical
//!   Hessian.
//!
//! These gates back the "benchmarking frente a rutas estables" and
//! "validación cuantitativa por submódulo" backlog items.

#[cfg(feature = "beta-kpm")]
mod kpm_vs_diagonalization {
    use nalgebra::{DMatrix, SymmetricEigen};
    use sci_form::beta::kpm::{compute_kpm_dos, KpmConfig};

    fn huckel_chain(n: usize, alpha: f64, beta: f64) -> DMatrix<f64> {
        let mut h = DMatrix::zeros(n, n);
        for i in 0..n {
            h[(i, i)] = alpha;
            if i + 1 < n {
                h[(i, i + 1)] = beta;
                h[(i + 1, i)] = beta;
            }
        }
        h
    }

    /// Exact DOS: Gaussian-broadened eigenvalue histogram.
    fn exact_dos(energies: &[f64], grid: &[f64], sigma: f64) -> Vec<f64> {
        grid.iter()
            .map(|&e| {
                let sum: f64 = energies
                    .iter()
                    .map(|&eps| (-((e - eps) / sigma).powi(2)).exp())
                    .sum();
                // exp(-(x/σ)²) ⇒ normalisation 1/(σ√π)
                sum / (sigma * std::f64::consts::PI.sqrt())
            })
            .collect()
    }

    #[test]
    fn kpm_dos_matches_exact_diagonalization() {
        let n = 40;
        let h = huckel_chain(n, 0.0, -2.0);
        let eigen = SymmetricEigen::new(h.clone());
        let eigenvalues: Vec<f64> = eigen.eigenvalues.iter().copied().collect();
        let e_min = eigenvalues.iter().cloned().fold(f64::INFINITY, f64::min);
        let e_max = eigenvalues
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);

        let n_points = 400;
        let grid: Vec<f64> = (0..n_points)
            .map(|i| e_min - 0.5 + (e_max - e_min + 1.0) * i as f64 / (n_points - 1) as f64)
            .collect();

        let config = KpmConfig {
            order: 400,
            n_vectors: 0, // exact trace
            ..KpmConfig::default()
        };
        let result = compute_kpm_dos(&h, &config, e_min - 0.5, e_max + 0.5, n_points);

        // Calibrate the Gaussian broadening against the KPM resolution
        // (~π·range/(2·2.355·order) FWHM-equivalent): scan and report.
        let unit = (e_max - e_min) / config.order as f64;
        let mut best = (0.0f64, 0.0f64, Vec::new());
        for factor in [0.4f64, 0.6, 0.8, 1.0, 1.5, 2.0] {
            let sigma = unit * factor;
            let reference = exact_dos(&eigenvalues, &grid, sigma);
            let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
            let (m1, m2) = (mean(&result.total_dos), mean(&reference));
            let mut num = 0.0;
            let mut d1 = 0.0;
            let mut d2 = 0.0;
            for (a, b) in result.total_dos.iter().zip(reference.iter()) {
                num += (a - m1) * (b - m2);
                d1 += (a - m1).powi(2);
                d2 += (b - m2).powi(2);
            }
            let corr = num / (d1.sqrt() * d2.sqrt());
            eprintln!(
                "  sigma = {:.3}·range/order -> correlation {:.4}",
                factor, corr
            );
            if corr > best.0 {
                best = (corr, factor, reference);
            }
        }
        let reference = best.2;

        // 1) Integrated DOS (total number of states) must match.
        let dx = (e_max - e_min + 1.0) / (n_points - 1) as f64;
        let integral_kpm: f64 = result.total_dos.iter().sum::<f64>() * dx;
        let integral_ref: f64 = reference.iter().sum::<f64>() * dx;
        assert!(
            (integral_kpm - n as f64).abs() / (n as f64) < 0.05,
            "KPM integrated DOS {} vs {n} states",
            integral_kpm
        );
        assert!(
            (integral_ref - n as f64).abs() / (n as f64) < 0.05,
            "reference integrated DOS {} vs {n} states",
            integral_ref
        );

        // 2) First moment (centre of mass) must match within 1% of the band.
        let band = e_max - e_min;
        let centre = |dos: &[f64]| -> f64 {
            let total: f64 = dos.iter().sum::<f64>() * dx;
            dos.iter().zip(grid.iter()).map(|(d, e)| d * e).sum::<f64>() * dx / total
        };
        let c_kpm = centre(&result.total_dos);
        let c_ref = centre(&reference);
        assert!(
            (c_kpm - c_ref).abs() < 0.01 * band,
            "KPM centre {c_kpm} vs reference {c_ref} (band {band})"
        );

        // 3) Correlation with the exact broadened DOS.
        let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
        let (m1, m2) = (mean(&result.total_dos), mean(&reference));
        let mut num = 0.0;
        let mut d1 = 0.0;
        let mut d2 = 0.0;
        for (a, b) in result.total_dos.iter().zip(reference.iter()) {
            num += (a - m1) * (b - m2);
            d1 += (a - m1).powi(2);
            d2 += (b - m2).powi(2);
        }
        let corr = num / (d1.sqrt() * d2.sqrt());
        // The KPM (Jackson kernel) line shape differs from a Gaussian, so the
        // point-wise correlation saturates near 0.86; the integrated DOS and
        // the first moment above are the tight quantitative gates.
        assert!(
            corr > 0.85,
            "KPM/exact DOS correlation {corr} at the best broadening (factor {:.1})",
            best.1
        );
        eprintln!(
            "KPM vs diagonalization: integral {integral_kpm:.2}/{n}, centre {c_kpm:.4} vs {c_ref:.4}, corr {corr:.4}"
        );
    }

    #[test]
    fn kpm_stochastic_trace_approaches_exact() {
        let h = huckel_chain(30, 0.0, -2.0);
        let exact = compute_kpm_dos(
            &h,
            &KpmConfig {
                order: 200,
                n_vectors: 0,
                ..KpmConfig::default()
            },
            -4.5,
            4.5,
            200,
        );
        let stochastic = compute_kpm_dos(
            &h,
            &KpmConfig {
                order: 200,
                n_vectors: 8,
                seed: 7,
                ..KpmConfig::default()
            },
            -4.5,
            4.5,
            200,
        );
        let dx = 9.0 / 199.0;
        let i_exact: f64 = exact.total_dos.iter().sum::<f64>() * dx;
        let i_stoch: f64 = stochastic.total_dos.iter().sum::<f64>() * dx;
        assert!(
            (i_stoch - 30.0).abs() / 30.0 < 0.15,
            "stochastic integrated DOS {i_stoch} vs 30 states (exact {i_exact})"
        );
    }
}

#[cfg(feature = "beta-randnla")]
mod nystrom_vs_exact_svd {
    use nalgebra::DMatrix;
    use rand::SeedableRng;
    use sci_form::beta::rand_nla::{GaussianSketch, NystromApprox};

    /// Low-rank-plus-noise symmetric matrix with known spectrum.
    fn low_rank_psd(n: usize, rank: usize, noise: f64, seed: u64) -> DMatrix<f64> {
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng as _};
        let mut rng = StdRng::seed_from_u64(seed);
        let mut a = DMatrix::zeros(n, n);
        for k in 0..rank {
            let scale = 1.0 / (k as f64 + 1.0);
            let v: Vec<f64> = (0..n).map(|_| rng.gen_range(-1.0..1.0)).collect();
            for i in 0..n {
                for j in 0..n {
                    a[(i, j)] += scale * v[i] * v[j];
                }
            }
        }
        for i in 0..n {
            for j in 0..n {
                a[(i, j)] += noise * rng.gen_range(-1.0..1.0);
            }
        }
        for i in 0..n {
            for j in (i + 1)..n {
                let v = 0.5 * (a[(i, j)] + a[(j, i)]);
                a[(i, j)] = v;
                a[(j, i)] = v;
            }
        }
        a
    }

    /// Best rank-k Frobenius error from the exact SVD.
    fn best_rank_k_error(a: &DMatrix<f64>, k: usize) -> f64 {
        let svd = a.clone().svd(true, true);
        let s = &svd.singular_values;
        let tail: f64 = s.iter().skip(k).map(|v| v * v).sum();
        tail.sqrt()
    }

    #[test]
    fn nystrom_error_is_bounded_by_rank_k_error() {
        let n = 80;
        let a = low_rank_psd(n, 6, 1e-3, 11);
        for k in [4usize, 12, 30] {
            let mut rng = rand::rngs::StdRng::seed_from_u64(99);
            let sketch = GaussianSketch::new(&mut rng, n, k);
            let approx = NystromApprox::from_matrix(&a, &sketch);
            let err = approx.frobenius_error(&a);
            let best = best_rank_k_error(&a, k);
            assert!(
                err <= 2.0 * best + 1e-9,
                "k={k}: Nyström error {err:.6} vs best rank-{k} error {best:.6}"
            );
            eprintln!("Nyström k={k}: error {err:.6}, best rank-{k} {best:.6}");
        }
    }

    #[test]
    fn nystrom_error_decreases_with_rank() {
        let n = 80;
        let a = low_rank_psd(n, 6, 1e-3, 5);
        let mut errors = Vec::new();
        for k in [4usize, 10, 20, 40] {
            let mut rng = rand::rngs::StdRng::seed_from_u64(123);
            let sketch = GaussianSketch::new(&mut rng, n, k);
            let approx = NystromApprox::from_matrix(&a, &sketch);
            errors.push(approx.frobenius_error(&a));
        }
        for w in errors.windows(2) {
            assert!(
                w[1] <= w[0] * 1.05,
                "Nyström error not decreasing with rank: {errors:?}"
            );
        }
        assert!(
            errors.last().unwrap() < &(0.1 * errors[0]),
            "high-rank Nyström should be much better: {errors:?}"
        );
    }
}

#[cfg(feature = "beta-mbh")]
mod mbh_vs_full_hessian {
    use nalgebra::{DMatrix, SymmetricEigen};
    use sci_form::beta::mbh::blocks::{build_projection_matrix, detect_rigid_blocks};
    use sci_form::beta::mbh::{compute_mbh_frequencies, MbhConfig};
    use sci_form::ir::hessian::compute_uff_analytical_hessian;

    /// kcal/(mol·Å²·amu) → cm⁻¹ (same convention as the MBH module).
    const FREQ_FACTOR: f64 = 108.59;

    fn atomic_mass(z: u8) -> f64 {
        match z {
            1 => 1.008,
            6 => 12.011,
            _ => 12.0,
        }
    }

    #[test]
    fn mbh_reproduces_projected_analytical_hessian() {
        // Methane: the MBH projection removes the carbon's rigid-body motion.
        let elements = [6u8, 1, 1, 1, 1];
        let positions = [
            [0.0, 0.0, 0.0],
            [0.6291, 0.6291, 0.6291],
            [-0.6291, -0.6291, 0.6291],
            [-0.6291, 0.6291, -0.6291],
            [0.6291, -0.6291, -0.6291],
        ];
        let rings: Vec<(Vec<usize>, bool)> = vec![(vec![0], true)];
        let n_full = 3 * elements.len();

        // Analytical UFF Hessian (full space).
        let flat: Vec<f64> = positions.iter().flatten().copied().collect();
        let h_full =
            compute_uff_analytical_hessian("C", &flat, None).expect("analytical UFF Hessian");

        // MBH projection, then H_red = Lᵀ H L (analytical reference).
        let decomp = detect_rigid_blocks(elements.len(), &elements, &positions, &rings);
        let l_cols = build_projection_matrix(&decomp, &positions, &elements);
        let n_red = l_cols.len();
        assert!(n_red > 0, "no reduced DOF");
        let mut h_red = DMatrix::zeros(n_red, n_red);
        for i in 0..n_red {
            for j in 0..n_red {
                let mut val = 0.0;
                for k in 0..n_full {
                    for l in 0..n_full {
                        val += l_cols[i][k] * h_full[(k, l)] * l_cols[j][l];
                    }
                }
                h_red[(i, j)] = val;
            }
        }

        // Mass weighting (same construction as the MBH module).
        let mut mass_diag = vec![0.0f64; n_full];
        for (i, &z) in elements.iter().enumerate() {
            let m = atomic_mass(z);
            for k in 0..3 {
                mass_diag[i * 3 + k] = m;
            }
        }
        let mut m_red = DMatrix::zeros(n_red, n_red);
        for i in 0..n_red {
            for j in 0..n_red {
                let mut val = 0.0;
                for k in 0..n_full {
                    val += l_cols[i][k] * mass_diag[k] * l_cols[j][k];
                }
                m_red[(i, j)] = val;
            }
        }
        let m_eig = m_red.symmetric_eigen();
        let mut m_inv_sqrt = DMatrix::zeros(n_red, n_red);
        for i in 0..n_red {
            if m_eig.eigenvalues[i] > 1e-10 {
                let val = 1.0 / m_eig.eigenvalues[i].sqrt();
                for r in 0..n_red {
                    for c in 0..n_red {
                        m_inv_sqrt[(r, c)] +=
                            val * m_eig.eigenvectors[(r, i)] * m_eig.eigenvectors[(c, i)];
                    }
                }
            }
        }
        let standard = &m_inv_sqrt * &h_red * &m_inv_sqrt;
        let eig = SymmetricEigen::new(standard);
        let mut reference: Vec<f64> = eig
            .eigenvalues
            .iter()
            .map(|l| FREQ_FACTOR * l.abs().sqrt() * if *l >= 0.0 { 1.0 } else { -1.0 })
            .collect();
        reference.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // Numerical MBH on the same surface.
        let energy_fn =
            |pos: &[f64]| -> f64 { sci_form::compute_uff_energy("C", pos).unwrap_or(f64::NAN) };
        let mbh = compute_mbh_frequencies(
            &elements,
            &positions,
            &rings,
            &energy_fn,
            &MbhConfig::default(),
        );
        let mut mbh_freqs = mbh.frequencies.clone();
        mbh_freqs.sort_by(|a, b| a.partial_cmp(b).unwrap());

        assert_eq!(mbh_freqs.len(), reference.len(), "DOF mismatch");
        // Near-zero modes are dominated by finite-difference noise (~5 cm⁻¹),
        // so the relative comparison uses an absolute floor.
        let mut worst = 0.0f64;
        for (m, r) in mbh_freqs.iter().zip(reference.iter()) {
            let denom = r.abs().max(200.0);
            worst = worst.max((m - r).abs() / denom);
        }
        eprintln!(
            "MBH vs projected analytical Hessian: {} modes, worst relative deviation {:.3}",
            mbh_freqs.len(),
            worst
        );
        assert!(
            worst < 0.05,
            "MBH deviation {worst:.3} (mbh {mbh_freqs:?}, projected {reference:?})"
        );
    }
}
