//! Integration tests for E10: Constant Potential Method

#[cfg(feature = "experimental-cpm")]
mod cpm_tests {
    use sci_form::experimental::cpm::*;

    fn water() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![8, 1, 1],
            vec![[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]],
        )
    }

    fn ethanol() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![6, 6, 8, 1, 1, 1, 1, 1, 1],
            vec![
                [0.0, 0.0, 0.0], [1.52, 0.0, 0.0], [2.14, 1.21, 0.0],
                [-0.36, 1.03, 0.0], [-0.36, -0.52, 0.89], [-0.36, -0.52, -0.89],
                [1.88, -0.52, 0.89], [1.88, -0.52, -0.89], [3.10, 1.21, 0.0],
            ],
        )
    }

    #[test]
    fn test_cpm_converges() {
        let (el, pos) = water();
        let result = compute_cpm_charges(&el, &pos, &CpmConfig::default());
        assert!(result.converged);
    }

    #[test]
    fn test_cpm_nonzero_charge() {
        let (el, pos) = water();
        // At non-neutral potential, total charge should deviate from zero
        let r1 = compute_cpm_charges(&el, &pos, &CpmConfig { mu_ev: -3.0, ..Default::default() });
        let r2 = compute_cpm_charges(&el, &pos, &CpmConfig { mu_ev: -6.0, ..Default::default() });
        assert!((r1.total_charge - r2.total_charge).abs() > 0.01,
            "Different μ should give different Q");
    }

    #[test]
    fn test_cpm_charge_increases_with_mu() {
        let (el, pos) = ethanol();
        let r_low = compute_cpm_charges(&el, &pos, &CpmConfig { mu_ev: -5.5, ..Default::default() });
        let r_high = compute_cpm_charges(&el, &pos, &CpmConfig { mu_ev: -3.5, ..Default::default() });
        assert!(r_high.total_charge > r_low.total_charge);
    }

    #[test]
    fn test_cpm_grand_potential_finite() {
        let (el, pos) = water();
        let result = compute_cpm_charges(&el, &pos, &CpmConfig::default());
        assert!(result.grand_potential.is_finite());
        assert!(result.electrostatic_energy.is_finite());
    }

    #[test]
    fn test_cpm_surface_scan() {
        let (el, pos) = water();
        let surface = compute_cpm_surface(&el, &pos, -5.5, -3.5, 10, 78.5);
        assert_eq!(surface.mu_values.len(), 10);
        assert_eq!(surface.total_charge.len(), 10);
        assert_eq!(surface.capacitance.len(), 10);
    }

    #[test]
    fn test_cpm_surface_charge_monotonic() {
        let (el, pos) = ethanol();
        let surface = compute_cpm_surface(&el, &pos, -5.5, -3.5, 20, 78.5);
        // Q should generally increase with μ
        let first = surface.total_charge[0];
        let last = *surface.total_charge.last().unwrap();
        assert!(last > first, "Q should increase: first={}, last={}", first, last);
    }

    #[test]
    fn test_cpm_surface_capacitance_positive() {
        let (el, pos) = water();
        let surface = compute_cpm_surface(&el, &pos, -5.0, -4.0, 15, 78.5);
        let positive_count = surface.capacitance.iter().filter(|&&c| c > 0.0).count();
        assert!(positive_count > surface.capacitance.len() / 2,
            "Most capacitance values should be positive");
    }

    #[test]
    fn test_cpm_dielectric_effect() {
        let (el, pos) = water();
        let r_water = compute_cpm_charges(&el, &pos, &CpmConfig { dielectric: 78.5, mu_ev: -4.0, ..Default::default() });
        let r_vacuum = compute_cpm_charges(&el, &pos, &CpmConfig { dielectric: 1.0, mu_ev: -4.0, ..Default::default() });
        // In vacuum (ε=1), Coulomb coupling is stronger → different charge distribution
        assert!((r_water.total_charge - r_vacuum.total_charge).abs() > 0.001,
            "Dielectric should affect charges");
    }
}
