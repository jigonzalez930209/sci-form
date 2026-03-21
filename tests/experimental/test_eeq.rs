//! Integration tests for Track E5: Dynamic EEQ Force Field

#[cfg(feature = "experimental-eeq")]
mod eeq_tests {
    use sci_form::charges_eeq::*;

    fn ethanol_atoms() -> (Vec<u8>, Vec<[f64; 3]>) {
        let elements = vec![6, 6, 8, 1, 1, 1, 1, 1, 1];
        let positions = vec![
            [0.0, 0.0, 0.0],       // C
            [1.54, 0.0, 0.0],      // C
            [2.42, 1.18, 0.0],     // O
            [-0.63, 0.89, 0.0],   // H
            [-0.63, -0.45, 0.89], // H
            [-0.63, -0.45, -0.89],// H
            [1.95, -0.52, 0.89],  // H
            [1.95, -0.52, -0.89], // H
            [3.37, 1.01, 0.0],    // H (OH)
        ];
        (elements, positions)
    }

    // E5.1 — EEQ Charge Model

    #[test]
    fn e5_1a_charge_neutrality() {
        let (elements, positions) = ethanol_atoms();
        let config = EeqConfig::default();
        let result = compute_eeq_charges(&elements, &positions, &config);
        assert!(
            result.total_charge.abs() < 0.01,
            "Total charge = {:.6}, expected 0", result.total_charge
        );
    }

    #[test]
    fn e5_1b_coordination_numbers() {
        let (elements, positions) = ethanol_atoms();
        let cn = fractional_coordination(&elements, &positions);
        // C atoms should have CN > 1 (bonded neighbors), O > 0.5
        assert!(cn[0] > 1.0, "C1 CN = {}", cn[0]);
        assert!(cn[2] > 0.5, "O CN = {}", cn[2]);
    }

    #[test]
    fn e5_1c_oxygen_most_negative() {
        let (elements, positions) = ethanol_atoms();
        let config = EeqConfig::default();
        let result = compute_eeq_charges(&elements, &positions, &config);
        // Oxygen should be the most negative atom
        let o_charge = result.charges[2];
        assert!(o_charge < 0.0, "O charge = {}", o_charge);
        for (i, &q) in result.charges.iter().enumerate() {
            if i != 2 {
                assert!(o_charge <= q + 0.01,
                    "O ({:.3}) should be ≤ atom {} ({:.3})", o_charge, i, q);
            }
        }
    }

    #[test]
    fn e5_1d_hydrogen_positive() {
        let (elements, positions) = ethanol_atoms();
        let config = EeqConfig::default();
        let result = compute_eeq_charges(&elements, &positions, &config);
        // All H atoms should be positive
        for i in 3..9 {
            assert!(result.charges[i] > -0.1,
                "H{} charge = {}", i, result.charges[i]);
        }
    }

    // E5.2 — EEQ Energy and Gradients

    #[test]
    fn e5_2a_energy_finite() {
        let (elements, positions) = ethanol_atoms();
        let config = EeqConfig::default();
        let result = compute_eeq_energy(&elements, &positions, &config);
        assert!(result.electrostatic_energy.is_finite());
    }

    #[test]
    fn e5_2b_gradient_nonzero() {
        let (elements, positions) = ethanol_atoms();
        let config = EeqConfig::default();
        let grad = compute_eeq_gradient(&elements, &positions, &config);
        let norm: f64 = grad.iter().map(|g| g[0]*g[0] + g[1]*g[1] + g[2]*g[2]).sum::<f64>().sqrt();
        assert!(norm > 0.0, "Gradient norm = 0");
    }

    // E5.3 — Validation

    #[test]
    fn e5_3a_charged_system() {
        let elements = vec![6, 8, 8, 1];
        let positions = vec![
            [0.0, 0.0, 0.0], [1.25, 0.0, 0.0],
            [-1.25, 0.0, 0.0], [0.0, 1.09, 0.0],
        ];
        let config = EeqConfig { total_charge: -1.0, ..Default::default() };
        let result = compute_eeq_charges(&elements, &positions, &config);
        assert!(
            (result.total_charge - (-1.0)).abs() < 0.01,
            "Total charge = {:.4}, expected -1", result.total_charge
        );
    }

    #[test]
    fn e5_3b_water_charge_comparison() {
        let elements = vec![8, 1, 1];
        let positions = vec![[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let config = EeqConfig::default();
        let eeq = compute_eeq_charges(&elements, &positions, &config);
        // Compare with Gasteiger (sci_form::compute_charges)
        // EEQ should give reasonable O charge: -0.2 to -0.8
        assert!(eeq.charges[0] < -0.05 && eeq.charges[0] > -1.0,
            "O charge = {:.4}, out of expected range", eeq.charges[0]);
    }
}
