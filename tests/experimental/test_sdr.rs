//! Integration tests for E8: SDR Embedding

#[cfg(feature = "experimental-sdr")]
mod sdr_tests {
    use sci_form::alpha::sdr::*;

    #[test]
    fn test_psd_projection_identity() {
        let m = nalgebra::DMatrix::identity(3, 3);
        let (p, neg) = project_psd(&m);
        assert_eq!(neg, 0);
        assert!((p - m).norm() < 1e-10);
    }

    #[test]
    fn test_psd_projection_removes_negatives() {
        let m = nalgebra::DMatrix::from_row_slice(
            3,
            3,
            &[2.0, 1.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 1.0],
        );
        let (p, neg) = project_psd(&m);
        assert!(neg > 0, "Should have removed negative eigenvalues");
        let eigen = p.symmetric_eigen();
        for &e in eigen.eigenvalues.iter() {
            assert!(e >= -1e-10, "All eigenvalues should be >= 0: {}", e);
        }
    }

    #[test]
    fn test_warm_start_gram_preserves_structure() {
        // Known coordinates
        let coords: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.866, 0.0]];
        let mut pairs = Vec::new();
        for i in 0..3 {
            for j in (i + 1)..3 {
                let d = ((coords[i][0] - coords[j][0]).powi(2)
                    + (coords[i][1] - coords[j][1]).powi(2)
                    + (coords[i][2] - coords[j][2]).powi(2))
                .sqrt();
                pairs.push((i, j, d));
            }
        }
        let gram = warm_start_gram(3, &pairs);
        assert_eq!(gram.nrows(), 3);
        // Gram matrix from valid distances should be approximately PSD
        let eigen = gram.symmetric_eigen();
        let min_eval = eigen.eigenvalues.min();
        assert!(
            min_eval > -0.1,
            "Should be approximately PSD: min={}",
            min_eval
        );
    }

    #[test]
    fn test_extract_coordinates_preserves_distances() {
        // Build Gram from 4 known 3D points
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [0.0, 1.5, 0.0],
            [0.0, 0.0, 1.5],
        ];
        let n = 4;
        let mut gram = nalgebra::DMatrix::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                gram[(i, j)] = pts[i].iter().zip(&pts[j]).map(|(a, b)| a * b).sum::<f64>();
            }
        }
        let coords = extract_coordinates(&gram);

        // Verify pairwise distances match
        for i in 0..n {
            for j in (i + 1)..n {
                let d_orig = ((pts[i][0] - pts[j][0]).powi(2)
                    + (pts[i][1] - pts[j][1]).powi(2)
                    + (pts[i][2] - pts[j][2]).powi(2))
                .sqrt();
                let d_new = ((coords[i * 3] - coords[j * 3]).powi(2)
                    + (coords[i * 3 + 1] - coords[j * 3 + 1]).powi(2)
                    + (coords[i * 3 + 2] - coords[j * 3 + 2]).powi(2))
                .sqrt();
                assert!(
                    (d_orig - d_new).abs() < 0.2,
                    "Distance ({},{}) orig={:.3} new={:.3}",
                    i,
                    j,
                    d_orig,
                    d_new
                );
            }
        }
    }

    #[test]
    fn test_sdr_embed_equilateral_triangle() {
        let d = 1.5;
        let pairs = vec![(0, 1, d), (0, 2, d), (1, 2, d)];
        let result = sdr_embed(3, &pairs, &SdrConfig::default());
        assert_eq!(result.num_atoms, 3);
        assert_eq!(result.coords.len(), 9);
    }

    #[test]
    fn test_sdr_embed_tetrahedron_distances() {
        let d = 2.0;
        let pairs = vec![
            (0, 1, d),
            (0, 2, d),
            (0, 3, d),
            (1, 2, d),
            (1, 3, d),
            (2, 3, d),
        ];
        let result = sdr_embed(4, &pairs, &SdrConfig::default());
        // All pairwise distances should be approximately d
        for &(i, j, dt) in &pairs {
            let dx = result.coords[i * 3] - result.coords[j * 3];
            let dy = result.coords[i * 3 + 1] - result.coords[j * 3 + 1];
            let dz = result.coords[i * 3 + 2] - result.coords[j * 3 + 2];
            let da = (dx * dx + dy * dy + dz * dz).sqrt();
            assert!(
                (da - dt).abs() < 1.0,
                "Distance ({},{}) target={:.2} actual={:.2}",
                i,
                j,
                dt,
                da
            );
        }
    }

    #[test]
    fn test_sdr_convergence_info() {
        let pairs = vec![(0, 1, 1.0), (0, 2, 1.0), (1, 2, 1.0)];
        let result = sdr_embed(
            3,
            &pairs,
            &SdrConfig {
                max_iter: 50,
                ..Default::default()
            },
        );
        assert!(result.convergence.iterations > 0);
        assert!(result.convergence.iterations <= 50);
    }

    #[test]
    fn test_alternating_projections_improves() {
        // Start with a non-PSD matrix and check that AP fixes it
        let n = 3;
        let x =
            nalgebra::DMatrix::from_row_slice(n, n, &[1.0, 2.0, 0.5, 2.0, 1.0, 0.5, 0.5, 0.5, 1.0]);
        let pairs = vec![(0, 1, 1.0), (0, 2, 1.5), (1, 2, 1.2)];
        let config = SdrConfig {
            max_iter: 100,
            tol: 1e-4,
            ..Default::default()
        };
        let (result, _conv) = alternating_projections(&x, &pairs, &config);

        // Result should be PSD
        let eigen = result.symmetric_eigen();
        for &e in eigen.eigenvalues.iter() {
            assert!(e >= -1e-8, "Should be PSD after AP: {}", e);
        }
    }
}
