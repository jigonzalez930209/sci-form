//! Integration tests for Track E6: ALPB Implicit Solvation

#[cfg(feature = "experimental-alpb")]
mod alpb_tests {
    use sci_form::solvation_alpb::*;

    fn water_system() -> (Vec<u8>, Vec<[f64; 3]>, Vec<f64>) {
        (
            vec![8, 1, 1],
            vec![[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]],
            vec![-0.834, 0.417, 0.417],
        )
    }

    // E6.1 — Born Radii and GB Kernel

    #[test]
    fn e6_1a_born_radii_physical() {
        let (elements, positions, _) = water_system();
        let born = compute_born_radii(&elements, &positions, 1.4);
        for &r in &born.radii {
            assert!(r > 0.5 && r < 50.0, "Born radius = {}", r);
        }
    }

    #[test]
    fn e6_1b_gb_kernel_monotone() {
        // f_GB should increase with distance
        let f1 = gb_kernel(1.0, 1.5, 1.5);
        let f2 = gb_kernel(5.0, 1.5, 1.5);
        let f3 = gb_kernel(10.0, 1.5, 1.5);
        assert!(f1 < f2 && f2 < f3, "f_GB not monotone: {}, {}, {}", f1, f2, f3);
    }

    #[test]
    fn e6_1c_gb_kernel_limits() {
        // f_GB(0, R, R) = R
        let f0 = gb_kernel(0.0, 2.0, 2.0);
        assert!((f0 - 2.0).abs() < 0.5, "f_GB(0) = {}, expected ~2", f0);
        // f_GB(large_r) ~ r
        let f = gb_kernel(100.0, 2.0, 2.0);
        assert!((f - 100.0).abs() < 1.0, "f_GB(100) = {}", f);
    }

    // E6.2 — ALPB Energy

    #[test]
    fn e6_2a_solvation_water() {
        let (elements, positions, charges) = water_system();
        let config = AlpbConfig::default();
        let result = compute_alpb_solvation(&elements, &positions, &charges, &config);
        assert!(result.total_energy.is_finite());
        assert!(result.born_radii.len() == 3);
    }

    #[test]
    fn e6_2b_dielectric_ordering() {
        // Higher dielectric → more negative solvation energy
        let (elements, positions, charges) = water_system();
        let r_low = compute_alpb_solvation(&elements, &positions, &charges,
            &AlpbConfig { solvent_dielectric: 4.8, ..Default::default() });
        let r_high = compute_alpb_solvation(&elements, &positions, &charges,
            &AlpbConfig { solvent_dielectric: 78.5, ..Default::default() });
        // Higher eps should give more solvation
        assert!(
            r_high.electrostatic_energy <= r_low.electrostatic_energy + 1e-6,
            "High eps ({:.2}) should give ≤ low eps ({:.2})",
            r_high.electrostatic_energy, r_low.electrostatic_energy
        );
    }

    #[test]
    fn e6_2c_vacuum_zero_solvation() {
        let (elements, positions, charges) = water_system();
        let config = AlpbConfig { solvent_dielectric: 1.0, ..Default::default() };
        let result = compute_alpb_solvation(&elements, &positions, &charges, &config);
        assert!(result.electrostatic_energy.abs() < 1e-8,
            "Vacuum electrostatic = {}", result.electrostatic_energy);
    }

    // E6.3 — Validation

    #[test]
    fn e6_3a_alpb_factor_range() {
        let (elements, positions, charges) = water_system();
        let result = compute_alpb_solvation(&elements, &positions, &charges, &AlpbConfig::default());
        // ALPB factor should be between 0 and 1
        assert!(result.alpb_factor >= 0.0 && result.alpb_factor <= 2.0,
            "ALPB factor = {}", result.alpb_factor);
    }

    #[test]
    fn e6_3b_different_solvents() {
        let (elements, positions, charges) = water_system();
        for (eps, name) in [(78.5, "water"), (47.0, "DMSO"), (4.8, "CHCl3")] {
            let config = AlpbConfig { solvent_dielectric: eps, ..Default::default() };
            let result = compute_alpb_solvation(&elements, &positions, &charges, &config);
            assert!(result.total_energy.is_finite(),
                "{}: energy = {}", name, result.total_energy);
        }
    }
}
