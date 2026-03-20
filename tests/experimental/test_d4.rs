//! Integration tests for E7: D4 Dispersion Correction

#[cfg(feature = "experimental-d4")]
mod d4_tests {
    use sci_form::experimental::d4::*;

    fn methane() -> (Vec<u8>, Vec<[f64; 3]>) {
        let elements = vec![6, 1, 1, 1, 1];
        let pos = vec![
            [0.0, 0.0, 0.0],
            [0.63, 0.63, 0.63],
            [-0.63, -0.63, 0.63],
            [-0.63, 0.63, -0.63],
            [0.63, -0.63, -0.63],
        ];
        (elements, pos)
    }

    fn ethane() -> (Vec<u8>, Vec<[f64; 3]>) {
        let elements = vec![6, 6, 1, 1, 1, 1, 1, 1];
        let pos = vec![
            [0.0, 0.0, 0.0],
            [1.54, 0.0, 0.0],
            [-0.36, 1.03, 0.0],
            [-0.36, -0.52, 0.89],
            [-0.36, -0.52, -0.89],
            [1.90, 1.03, 0.0],
            [1.90, -0.52, 0.89],
            [1.90, -0.52, -0.89],
        ];
        (elements, pos)
    }

    #[test]
    fn test_d4_energy_attractive() {
        let (el, pos) = methane();
        let result = compute_d4_energy(&el, &pos, &D4Config::default());
        assert!(result.total_energy < 0.0,
            "Dispersion should be attractive: {}", result.total_energy);
    }

    #[test]
    fn test_d4_coordination_numbers() {
        let (el, pos) = methane();
        let result = compute_d4_energy(&el, &pos, &D4Config::default());
        // Carbon should have CN > 1 (bonded to 4H), hydrogens > 0.3
        assert!(result.coordination_numbers[0] > 1.0,
            "C should have CN>1: {}", result.coordination_numbers[0]);
        for i in 1..5 {
            assert!(result.coordination_numbers[i] > 0.3,
                "H should have CN>0.3: {}", result.coordination_numbers[i]);
        }
    }

    #[test]
    fn test_d4_larger_molecule_more_dispersion() {
        let (el1, pos1) = methane();
        let (el2, pos2) = ethane();
        let e1 = compute_d4_energy(&el1, &pos1, &D4Config::default()).total_energy.abs();
        let e2 = compute_d4_energy(&el2, &pos2, &D4Config::default()).total_energy.abs();
        assert!(e2 > e1, "Ethane should have more dispersion than methane");
    }

    #[test]
    fn test_d4_distance_decay() {
        let el = vec![6u8, 6];
        let near = compute_d4_energy(&el,
            &[[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]], &D4Config::default()).total_energy;
        let far = compute_d4_energy(&el,
            &[[0.0, 0.0, 0.0], [12.0, 0.0, 0.0]], &D4Config::default()).total_energy;
        assert!(near.abs() > far.abs(), "D4 should decay with distance");
    }

    #[test]
    fn test_d4_three_body_contribution() {
        let (el, pos) = methane();
        let r2 = compute_d4_energy(&el, &pos, &D4Config { three_body: false, ..Default::default() });
        let r3 = compute_d4_energy(&el, &pos, &D4Config { three_body: true, ..Default::default() });
        // 3-body should add a (small) correction
        assert!((r3.total_energy - r2.total_energy).abs() > 1e-12,
            "3-body contribution should be nonzero");
        // 3-body is typically much smaller than 2-body
        assert!(r3.e3_body.abs() < r3.e2_body.abs(),
            "3-body should be smaller than 2-body");
    }

    #[test]
    fn test_d4_gradient_consistency() {
        let (el, pos) = methane();
        let config = D4Config::default();
        let grad = compute_d4_gradient(&el, &pos, &config);
        assert_eq!(grad.len(), el.len());
        // Check gradients are finite and not all zero
        let total_mag: f64 = grad.iter()
            .flat_map(|g| g.iter())
            .map(|v| v * v)
            .sum();
        assert!(total_mag > 0.0, "Gradients should not be all zero");
    }

    #[test]
    fn test_d4_s8_scaling() {
        let (el, pos) = ethane();
        let e1 = compute_d4_energy(&el, &pos, &D4Config { s8: 0.5, ..Default::default() });
        let e2 = compute_d4_energy(&el, &pos, &D4Config { s8: 1.5, ..Default::default() });
        // Higher s8 → more C8 contribution → more negative
        assert!(e2.total_energy < e1.total_energy,
            "Higher s8 should give more dispersion");
    }

    #[test]
    fn test_d4_single_atom_zero() {
        let result = compute_d4_energy(&[6], &[[0.0, 0.0, 0.0]], &D4Config::default());
        assert!((result.total_energy).abs() < 1e-15, "Single atom should have zero dispersion");
    }
}
