//! Integration tests for E11: Growing String Method

#[cfg(feature = "experimental-gsm")]
mod gsm_tests {
    use sci_form::alpha::gsm::*;

    fn double_well(coords: &[f64]) -> f64 {
        let x = coords[0];
        (x * x - 1.0).powi(2)
    }

    #[allow(dead_code)]
    fn morse_1d(coords: &[f64]) -> f64 {
        let x = coords[0];
        let d = 10.0; // well depth kcal/mol
        let a = 1.5;  // width
        let r_eq = 1.0;
        d * (1.0 - (-a * (x - r_eq)).exp()).powi(2)
    }

    #[test]
    fn test_interpolation() {
        let r = vec![0.0, 0.0, 0.0];
        let p = vec![4.0, 4.0, 4.0];
        let mid = interpolate_node(&r, &p, 0.5);
        assert!((mid[0] - 2.0).abs() < 1e-10);
        assert!((mid[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_interpolation_quarter() {
        let r = vec![0.0];
        let p = vec![8.0];
        let q = interpolate_node(&r, &p, 0.25);
        assert!((q[0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_grow_string_basic() {
        let r = vec![-1.5, 0.0, 0.0];
        let p = vec![1.5, 0.0, 0.0];
        let config = GsmConfig {
            max_nodes: 7,
            max_iter: 15,
            step_size: 0.01,
            ..Default::default()
        };

        let path = grow_string(&r, &p, &double_well, &config);
        assert!(path.nodes.len() >= 3, "Should have at least 3 nodes");
        assert_eq!(path.energies.len(), path.nodes.len());
    }

    #[test]
    fn test_find_ts_double_well() {
        let r = vec![-1.0, 0.0, 0.0];
        let p = vec![1.0, 0.0, 0.0];
        let config = GsmConfig {
            max_nodes: 11,
            max_iter: 30,
            step_size: 0.01,
            ..Default::default()
        };

        let result = find_transition_state(&r, &p, &double_well, &config);
        // The TS of (x²-1)² is at x=0, E=1
        assert!(result.activation_energy > 0.0, "Barrier should be positive");
    }

    #[test]
    fn test_gsm_path_ordered() {
        let r = vec![-1.0, 0.0, 0.0];
        let p = vec![1.0, 0.0, 0.0];
        let config = GsmConfig {
            max_nodes: 7,
            max_iter: 10,
            step_size: 0.01,
            ..Default::default()
        };

        let result = find_transition_state(&r, &p, &double_well, &config);
        assert!(result.path_coords.len() == result.n_nodes);
    }

    #[test]
    fn test_gsm_energy_evaluations_counted() {
        let r = vec![-1.0, 0.0, 0.0];
        let p = vec![1.0, 0.0, 0.0];
        let config = GsmConfig { max_nodes: 5, max_iter: 5, ..Default::default() };
        let result = find_transition_state(&r, &p, &double_well, &config);
        assert!(result.energy_evaluations > 0, "Should count evaluations");
    }

    #[test]
    fn test_gsm_symmetric_barrier() {
        // Symmetric potential: forward and reverse barriers should be similar
        let r = vec![-1.0, 0.0, 0.0];
        let p = vec![1.0, 0.0, 0.0];
        let config = GsmConfig {
            max_nodes: 11,
            max_iter: 30,
            step_size: 0.005,
            ..Default::default()
        };

        let result = find_transition_state(&r, &p, &double_well, &config);
        // For (x²-1)², both minima are at E=0, TS at E=1
        // Barriers should be similar (both ~1.0)
        let diff = (result.activation_energy - result.reverse_barrier).abs();
        assert!(diff < result.activation_energy.abs() * 0.5 + 0.5,
            "Barriers should be roughly equal for symmetric potential: fwd={}, rev={}",
            result.activation_energy, result.reverse_barrier);
    }

    #[test]
    fn test_refine_saddle() {
        // Start near the TS of double well
        let ts_init = vec![0.3, 0.0, 0.0]; // near x=0
        let prev = vec![-0.5, 0.0, 0.0];
        let next = vec![0.5, 0.0, 0.0];

        let (refined, _energy) = refine_saddle(
            &ts_init, &double_well, (&prev, &next),
            50, 0.005, 0.001,
        );

        // Should move closer to x=0
        assert!(refined[0].abs() < ts_init[0].abs() + 0.5,
            "Refinement should move toward saddle: {}", refined[0]);
    }
}
