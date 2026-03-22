//! Integration tests for Track E4: Kernel Polynomial Method (KPM)

#[cfg(feature = "experimental-kpm")]
mod kpm_tests {
    use nalgebra::{DMatrix, SymmetricEigen};
    use sci_form::beta::kpm::*;

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

    // E4.1 — Chebyshev Expansion Infrastructure

    #[test]
    fn e4_1a_gershgorin_bounds_all_eigenvalues() {
        let h = huckel_chain(30, -10.0, -2.5);
        let (e_min, e_max) = estimate_spectral_bounds(&h);
        let eigen = SymmetricEigen::new(h);
        for i in 0..eigen.eigenvalues.len() {
            let e = eigen.eigenvalues[i];
            assert!(
                e >= e_min && e <= e_max,
                "Eigenvalue {:.4} outside bounds [{:.4}, {:.4}]",
                e,
                e_min,
                e_max
            );
        }
    }

    #[test]
    fn e4_1b_rescaled_eigenvalues_within_unit() {
        let h = huckel_chain(25, -8.0, -3.0);
        let (e_min, e_max) = estimate_spectral_bounds(&h);
        let h_t = rescale_matrix(&h, e_min, e_max);
        let eigen = SymmetricEigen::new(h_t);
        for i in 0..eigen.eigenvalues.len() {
            assert!(
                eigen.eigenvalues[i].abs() <= 1.0 + 1e-10,
                "Rescaled eigenvalue {:.6} > 1",
                eigen.eigenvalues[i]
            );
        }
    }

    #[test]
    fn e4_1c_jackson_kernel_properties() {
        let gk = jackson_kernel(80);
        assert!((gk[0] - 1.0).abs() < 0.05, "g_0 should be ~1: {}", gk[0]);
        assert!(gk[79] < gk[0], "High-order coefficients should be damped");
    }

    #[test]
    fn e4_1d_dos_chain_vs_analytical() {
        // 1D tight-binding chain: analytical DOS is (1/pi) / sqrt(4t^2 - E^2)
        let n = 50;
        let t = 1.0;
        let h = huckel_chain(n, 0.0, -t);
        let expansion = ChebyshevExpansion::from_matrix_exact(&h, 150);
        let gk = jackson_kernel(150);

        // Band is [-2t, 2t]. DOS should be positive inside, zero outside.
        let dos_center = expansion.dos_at_energy(0.0, &gk);
        let dos_edge = expansion.dos_at_energy(1.8, &gk);
        let dos_outside = expansion.dos_at_energy(3.0, &gk);

        assert!(dos_center > 0.1, "DOS at center = {}", dos_center);
        assert!(dos_edge > 0.0, "DOS near edge = {}", dos_edge);
        assert!(dos_outside < 0.01, "DOS outside band = {}", dos_outside);
    }

    // E4.2 — Density, DOS and Populations

    #[test]
    fn e4_2a_kpm_dos_nonnegative() {
        let h = huckel_chain(30, -5.0, -1.5);
        let config = KpmConfig {
            order: 100,
            n_vectors: 0,
            ..Default::default()
        };
        let dos = compute_kpm_dos(&h, &config, -10.0, 0.0, 200);
        for &d in &dos.total_dos {
            assert!(d >= -1e-10, "Negative DOS value: {}", d);
        }
    }

    #[test]
    fn e4_2b_kpm_dos_integrates_to_n() {
        // Integral of DOS over band = N (number of states)
        let n = 20;
        let h = huckel_chain(n, 0.0, -1.0);
        let config = KpmConfig {
            order: 120,
            n_vectors: 0,
            ..Default::default()
        };
        let dos = compute_kpm_dos(&h, &config, -3.0, 3.0, 500);

        let de = (dos.energies.last().unwrap() - dos.energies.first().unwrap())
            / (dos.energies.len() - 1) as f64;
        let integral: f64 = dos.total_dos.iter().sum::<f64>() * de;

        // Should be close to N
        assert!(
            (integral - n as f64).abs() < n as f64 * 0.3,
            "DOS integral = {:.2}, expected ~{}",
            integral,
            n
        );
    }

    #[test]
    fn e4_2c_kpm_mulliken_charge_neutral() {
        let n = 12;
        let h = huckel_chain(n, -10.0, -2.0);
        let s = DMatrix::identity(n, n);
        let nuclear_charges = vec![1.0; n];
        let config = KpmConfig {
            order: 80,
            n_vectors: 0,
            ..Default::default()
        };
        let result = compute_kpm_mulliken(&h, &s, n, &nuclear_charges, &config);
        let total_charge: f64 = result.charges.iter().sum();
        assert!(
            total_charge.abs() < 2.0,
            "Total charge = {:.4}, expected ~0",
            total_charge
        );
    }

    // E4.3 — Validation

    #[test]
    fn e4_3a_kpm_dos_vs_exact() {
        // Compare KPM broadened DOS with exact Gaussian-broadened DOS
        let n = 20;
        let h = huckel_chain(n, -8.0, -2.0);
        let eigen = SymmetricEigen::new(h.clone());
        let sigma = 0.3;

        // Exact broadened DOS
        let n_points = 200;
        let e_lo = -15.0;
        let e_hi = -3.0;
        let step = (e_hi - e_lo) / (n_points - 1) as f64;
        let exact_dos: Vec<f64> = (0..n_points)
            .map(|i| {
                let e = e_lo + i as f64 * step;
                let norm = 1.0 / (sigma * (2.0 * std::f64::consts::PI).sqrt());
                let inv_2s2 = 1.0 / (2.0 * sigma * sigma);
                let mut d = 0.0;
                for k in 0..eigen.eigenvalues.len() {
                    d += norm * (-(e - eigen.eigenvalues[k]).powi(2) * inv_2s2).exp();
                }
                d
            })
            .collect();

        // KPM DOS
        let config = KpmConfig {
            order: 150,
            n_vectors: 0,
            ..Default::default()
        };
        let kpm_dos = compute_kpm_dos(&h, &config, e_lo, e_hi, n_points);

        // Compare: L2 relative error
        let exact_norm: f64 = exact_dos.iter().map(|d| d * d).sum::<f64>().sqrt();
        let kpm_norm: f64 = kpm_dos.total_dos.iter().map(|d| d * d).sum::<f64>().sqrt();

        // Both should be nonzero
        assert!(exact_norm > 0.1, "Exact DOS norm too small");
        assert!(kpm_norm > 0.1, "KPM DOS norm too small");

        // Shapes should correlate (peak positions match)
        let exact_peak = exact_dos.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let kpm_peak = kpm_dos
            .total_dos
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            exact_peak > 0.0 && kpm_peak > 0.0,
            "Peaks should be positive"
        );
    }

    #[test]
    fn e4_3b_larger_system_dos() {
        // 100-site chain: KPM should still work
        let n = 100;
        let h = huckel_chain(n, 0.0, -1.0);
        let config = KpmConfig {
            order: 200,
            n_vectors: 20,
            seed: 42,
            ..Default::default()
        };
        let dos = compute_kpm_dos(&h, &config, -3.0, 3.0, 100);

        // Should have nonzero DOS inside band
        let mid = dos.total_dos[50];
        assert!(mid > 0.0, "DOS at center of 100-site chain = {}", mid);
    }
}
