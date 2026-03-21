//! Integration tests for E9: Mobile Block Hessian

#[cfg(feature = "experimental-mbh")]
mod mbh_tests {
    use sci_form::beta::mbh::*;

    fn harmonic_energy(coords: &[f64]) -> f64 {
        let k = 100.0;
        let r0 = 1.4;
        let n = coords.len() / 3;
        let mut e = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = coords[i * 3] - coords[j * 3];
                let dy = coords[i * 3 + 1] - coords[j * 3 + 1];
                let dz = coords[i * 3 + 2] - coords[j * 3 + 2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                if r < 3.0 {
                    e += 0.5 * k * (r - r0).powi(2);
                }
            }
        }
        e
    }

    fn benzene_like() -> (Vec<u8>, Vec<[f64; 3]>) {
        let elements = vec![6u8; 12]; // 6 ring + 6 external
        let positions: Vec<[f64; 3]> = (0..12).map(|i| {
            if i < 6 {
                let a = i as f64 * std::f64::consts::PI / 3.0;
                [1.4 * a.cos(), 1.4 * a.sin(), 0.0]
            } else {
                let a = (i - 6) as f64 * std::f64::consts::PI / 3.0;
                [2.5 * a.cos(), 2.5 * a.sin(), 0.0]
            }
        }).collect();
        (elements, positions)
    }

    #[test]
    fn test_block_detection_benzene() {
        let (el, pos) = benzene_like();
        let rings = vec![(vec![0, 1, 2, 3, 4, 5], true)];
        let decomp = detect_rigid_blocks(12, &el, &pos, &rings);
        assert_eq!(decomp.blocks.len(), 1);
        assert_eq!(decomp.flexible_atoms.len(), 6);
        assert!(decomp.blocks[0].is_aromatic);
    }

    #[test]
    fn test_block_detection_no_rings() {
        let elements = vec![6, 1, 1, 1, 1];
        let pos = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0], [0.0, -1.0, 0.0],
        ];
        let decomp = detect_rigid_blocks(5, &elements, &pos, &[]);
        assert_eq!(decomp.blocks.len(), 0);
        assert_eq!(decomp.flexible_atoms.len(), 5);
    }

    #[test]
    fn test_projection_matrix_shape() {
        let (el, pos) = benzene_like();
        let rings = vec![(vec![0, 1, 2, 3, 4, 5], true)];
        let decomp = detect_rigid_blocks(12, &el, &pos, &rings);
        let l = build_projection_matrix(&decomp, &pos, &el);
        // 6 DOF for block + 18 DOF for 6 flexible atoms = 24 columns
        assert_eq!(l.len(), 24);
        for col in &l {
            assert_eq!(col.len(), 36); // 3 * 12 atoms
        }
    }

    #[test]
    fn test_mbh_speedup() {
        let (el, pos) = benzene_like();
        let rings = vec![(vec![0, 1, 2, 3, 4, 5], true)];
        let result = compute_mbh_frequencies(&el, &pos, &rings, &harmonic_energy, &MbhConfig::default());
        // With 1 block of 6, 6 flexible: reduced=24, full=36 → speedup=1.5
        assert!(result.speedup > 1.0, "Speedup: {}", result.speedup);
    }

    #[test]
    fn test_mbh_fewer_dof() {
        let (el, pos) = benzene_like();
        let rings = vec![(vec![0, 1, 2, 3, 4, 5], true)];
        let result = compute_mbh_frequencies(&el, &pos, &rings, &harmonic_energy, &MbhConfig::default());
        assert!(result.n_dof_reduced < result.n_dof_full);
    }

    #[test]
    fn test_mbh_frequencies_sorted() {
        let (el, pos) = benzene_like();
        let rings = vec![(vec![0, 1, 2, 3, 4, 5], true)];
        let result = compute_mbh_frequencies(&el, &pos, &rings, &harmonic_energy, &MbhConfig::default());
        for w in result.frequencies.windows(2) {
            assert!(w[0] <= w[1], "Frequencies should be sorted: {} > {}", w[0], w[1]);
        }
    }

    #[test]
    fn test_mbh_all_flexible_matches_full() {
        let elements = vec![6, 1, 1];
        let pos = vec![[0.0, 0.0, 0.0], [1.4, 0.0, 0.0], [-1.4, 0.0, 0.0]];
        let result = compute_mbh_frequencies(&elements, &pos, &[], &harmonic_energy, &MbhConfig::default());
        assert_eq!(result.n_dof_reduced, result.n_dof_full);
        assert!((result.speedup - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_mbh_frequencies_finite() {
        let (el, pos) = benzene_like();
        let rings = vec![(vec![0, 1, 2, 3, 4, 5], true)];
        let result = compute_mbh_frequencies(&el, &pos, &rings, &harmonic_energy, &MbhConfig::default());
        for &f in &result.frequencies {
            assert!(f.is_finite(), "All frequencies should be finite");
        }
    }
}
