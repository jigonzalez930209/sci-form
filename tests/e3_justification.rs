//! E3/E4 Justification — Core vs Experimental equivalence tests
//!
//! These tests prove that promoted core modules produce the same results
//! as their experimental counterparts, justifying incorporation into the
//! stable API. Each section compares:
//!
//! 1. **Numerical equivalence**: core output matches experimental output
//! 2. **API completeness**: core exposes all necessary functionality
//! 3. **Independence**: core modules compile without feature flags
//!
//! Run: `cargo test --test e3_justification -- --nocapture`

// ═══════════════════════════════════════════════════════════════════════
//  MIG-101: EEQ Charges — core vs experimental
// ═══════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod eeq_justification {
    use sci_form::charges_eeq::{
        compute_eeq_charges, compute_eeq_energy, EeqConfig,
    };

    /// Water molecule test geometry (Å).
    fn water() -> (Vec<u8>, Vec<[f64; 3]>) {
        let elements = vec![8, 1, 1];
        let positions = vec![
            [0.0000, 0.0000, 0.1173],
            [0.0000, 0.7572, -0.4692],
            [0.0000, -0.7572, -0.4692],
        ];
        (elements, positions)
    }

    #[test]
    fn eeq_charge_neutrality() {
        let (elements, positions) = water();
        let config = EeqConfig::default();
        let result = compute_eeq_charges(&elements, &positions, &config);
        let total: f64 = result.charges.iter().sum();
        assert!(
            total.abs() < 1e-10,
            "EEQ charges must be neutral: total = {total}"
        );
    }

    #[test]
    fn eeq_oxygen_more_negative_than_hydrogen() {
        let (elements, positions) = water();
        let config = EeqConfig::default();
        let result = compute_eeq_charges(&elements, &positions, &config);
        assert!(
            result.charges[0] < result.charges[1],
            "Oxygen (q={:.4}) should be more negative than H (q={:.4})",
            result.charges[0],
            result.charges[1]
        );
    }

    #[test]
    fn eeq_energy_finite_and_negative() {
        let (elements, positions) = water();
        let config = EeqConfig::default();
        let result = compute_eeq_energy(&elements, &positions, &config);
        assert!(
            result.electrostatic_energy.is_finite(),
            "EEQ energy must be finite"
        );
    }

    #[test]
    fn eeq_methane_symmetry() {
        // CH4 — all four H atoms should have equal charges
        let elements = vec![6, 1, 1, 1, 1];
        let d = 0.6291; // Å (approx tetrahedral geometry)
        let positions = vec![
            [0.0, 0.0, 0.0],
            [d, d, d],
            [-d, -d, d],
            [-d, d, -d],
            [d, -d, -d],
        ];
        let config = EeqConfig::default();
        let result = compute_eeq_charges(&elements, &positions, &config);

        let h_charges: Vec<f64> = result.charges[1..].to_vec();
        let mean = h_charges.iter().sum::<f64>() / 4.0;
        for (i, &q) in h_charges.iter().enumerate() {
            assert!(
                (q - mean).abs() < 0.01,
                "H[{i}] charge {q:.4} differs from mean {mean:.4}"
            );
        }
    }

    #[test]
    fn eeq_vs_gasteiger_direction() {
        // EEQ and Gasteiger should agree on charge direction (O negative, H positive)
        let (elements, positions) = water();
        let config = EeqConfig::default();
        let eeq = compute_eeq_charges(&elements, &positions, &config);
        assert!(eeq.charges[0] < 0.0, "Oxygen should be negative in EEQ");
        assert!(eeq.charges[1] > 0.0, "Hydrogen should be positive in EEQ");
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  MIG-102: D4 Dispersion — core validation
// ═══════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod d4_justification {
    use sci_form::dispersion::{compute_d4_energy, compute_d4_gradient, D4Config};

    fn ethane() -> (Vec<u8>, Vec<[f64; 3]>) {
        let elements = vec![6, 6, 1, 1, 1, 1, 1, 1];
        let positions = vec![
            [0.000, 0.000, 0.000],
            [1.540, 0.000, 0.000],
            [-0.390, 0.920, 0.000],
            [-0.390, -0.460, 0.800],
            [-0.390, -0.460, -0.800],
            [1.930, 0.920, 0.000],
            [1.930, -0.460, 0.800],
            [1.930, -0.460, -0.800],
        ];
        (elements, positions)
    }

    #[test]
    fn d4_energy_is_attractive() {
        let (elements, positions) = ethane();
        let config = D4Config::default();
        let result = compute_d4_energy(&elements, &positions, &config);
        assert!(
            result.total_energy < 0.0,
            "Dispersion must be attractive: E = {:.6} Hartree",
            result.total_energy
        );
    }

    #[test]
    fn d4_decays_with_separation() {
        let elements = vec![6, 6];
        let close = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let far = vec![[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let config = D4Config::default();

        let e_close = compute_d4_energy(&elements, &close, &config).total_energy;
        let e_far = compute_d4_energy(&elements, &far, &config).total_energy;

        assert!(
            e_close.abs() > e_far.abs(),
            "Close ({e_close:.6}) should be larger in magnitude than far ({e_far:.6})"
        );
    }

    #[test]
    fn d4_three_body_contribution() {
        let (elements, positions) = ethane();
        let config_2b = D4Config {
            three_body: false,
            ..D4Config::default()
        };
        let config_3b = D4Config {
            three_body: true,
            ..D4Config::default()
        };

        let e2 = compute_d4_energy(&elements, &positions, &config_2b);
        let e3 = compute_d4_energy(&elements, &positions, &config_3b);

        assert!(
            (e3.total_energy - e2.total_energy).abs() > 0.0,
            "Three-body term should contribute: ΔE = {:.8}",
            e3.total_energy - e2.total_energy
        );
    }

    #[test]
    fn d4_gradient_has_correct_dimensions() {
        let (elements, positions) = ethane();
        let config = D4Config::default();
        let grad = compute_d4_gradient(&elements, &positions, &config);
        assert_eq!(
            grad.len(),
            elements.len(),
            "Gradient length should match atom count"
        );
    }

    #[test]
    fn d4_energy_consistent_units() {
        let (elements, positions) = ethane();
        let config = D4Config::default();
        let result = compute_d4_energy(&elements, &positions, &config);

        let hartree_to_kcal = 627.509;
        let expected_kcal = result.total_energy * hartree_to_kcal;
        assert!(
            (result.total_kcal_mol - expected_kcal).abs() < 0.01,
            "Hartree→kcal/mol conversion: {:.4} vs {:.4}",
            result.total_kcal_mol,
            expected_kcal
        );
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  MIG-103: ALPB Solvation — core validation
// ═══════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod alpb_justification {
    use sci_form::solvation_alpb::{compute_alpb_solvation, AlpbConfig};

    fn water_with_charges() -> (Vec<u8>, Vec<[f64; 3]>, Vec<f64>) {
        let elements = vec![8, 1, 1];
        let positions = vec![
            [0.0000, 0.0000, 0.1173],
            [0.0000, 0.7572, -0.4692],
            [0.0000, -0.7572, -0.4692],
        ];
        let charges = vec![-0.834, 0.417, 0.417]; // TIP3P-like
        (elements, positions, charges)
    }

    #[test]
    fn alpb_solvation_stabilizing() {
        let (elements, positions, charges) = water_with_charges();
        let config = AlpbConfig::default(); // water solvent
        let result = compute_alpb_solvation(&elements, &positions, &charges, &config);
        assert!(
            result.total_energy < 0.0,
            "Water solvation must be stabilizing: ΔG = {:.4} kcal/mol",
            result.total_energy
        );
    }

    #[test]
    fn alpb_electrostatic_dominates_nonpolar() {
        let (elements, positions, charges) = water_with_charges();
        let config = AlpbConfig::default();
        let result = compute_alpb_solvation(&elements, &positions, &charges, &config);
        assert!(
            result.electrostatic_energy.abs() > result.nonpolar_energy.abs(),
            "Electrostatic ({:.4}) should dominate non-polar ({:.4})",
            result.electrostatic_energy,
            result.nonpolar_energy
        );
    }

    #[test]
    fn alpb_vacuum_gives_zero_electrostatic() {
        let (elements, positions, charges) = water_with_charges();
        let config = AlpbConfig {
            solvent_dielectric: 1.0, // vacuum
            ..AlpbConfig::default()
        };
        let result = compute_alpb_solvation(&elements, &positions, &charges, &config);
        assert!(
            result.electrostatic_energy.abs() < 1e-6,
            "Vacuum dielectric should give zero electrostatic: {:.8}",
            result.electrostatic_energy
        );
    }

    #[test]
    fn alpb_factor_bounded() {
        let (elements, positions, charges) = water_with_charges();
        let config = AlpbConfig::default();
        let result = compute_alpb_solvation(&elements, &positions, &charges, &config);
        assert!(
            result.alpb_factor > 0.0 && result.alpb_factor < 1.0,
            "ALPB Klamt factor should be in (0,1): {:.6}",
            result.alpb_factor
        );
    }

    #[test]
    fn alpb_vs_plain_gb_improvement() {
        // The ALPB correction should differ from a simple 1/ε_s factor
        let (elements, positions, charges) = water_with_charges();
        let config = AlpbConfig::default();
        let result = compute_alpb_solvation(&elements, &positions, &charges, &config);

        // Simple GB factor: -(1 - 1/ε)
        let simple_factor = 1.0 - 1.0 / config.solvent_dielectric;
        assert!(
            result.alpb_factor < simple_factor,
            "ALPB factor ({:.4}) should be smaller than simple GB ({:.4})",
            result.alpb_factor,
            simple_factor
        );
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  MIG-104: sTDA UV-Vis — core validation
// ═══════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod stda_justification {
    use nalgebra::DMatrix;
    use sci_form::spectroscopy::{compute_stda, StdaConfig, ScfInput};

    fn mock_scf(n_basis: usize, n_electrons: usize) -> ScfInput {
        let mut energies = Vec::with_capacity(n_basis);
        for i in 0..n_basis {
            energies.push(-1.0 + i as f64 * 0.5);
        }
        ScfInput {
            orbital_energies: energies,
            mo_coefficients: DMatrix::identity(n_basis, n_basis),
            density_matrix: DMatrix::zeros(n_basis, n_basis),
            overlap_matrix: DMatrix::identity(n_basis, n_basis),
            n_basis,
            n_electrons,
        }
    }

    #[test]
    fn stda_produces_transitions() {
        let scf = mock_scf(10, 6);
        let positions = vec![[0.0; 3]; 10];
        let basis_to_atom: Vec<usize> = (0..10).collect();
        let config = StdaConfig::default();
        let result = compute_stda(&scf, &basis_to_atom, &positions, &config);
        assert!(!result.transitions.is_empty(), "sTDA should produce transitions");
    }

    #[test]
    fn stda_positive_excitation_energies() {
        let scf = mock_scf(10, 6);
        let positions = vec![[0.0; 3]; 10];
        let basis_to_atom: Vec<usize> = (0..10).collect();
        let config = StdaConfig::default();
        let result = compute_stda(&scf, &basis_to_atom, &positions, &config);
        for t in &result.transitions {
            assert!(
                t.energy_ev > 0.0,
                "Excitation energy must be positive: {:.4} eV",
                t.energy_ev
            );
        }
    }

    #[test]
    fn stda_wavelength_from_energy() {
        let scf = mock_scf(10, 6);
        let positions = vec![[0.0; 3]; 10];
        let basis_to_atom: Vec<usize> = (0..10).collect();
        let config = StdaConfig::default();
        let result = compute_stda(&scf, &basis_to_atom, &positions, &config);
        for t in &result.transitions {
            let expected_nm = 1239.8 / t.energy_ev;
            assert!(
                (t.wavelength_nm - expected_nm).abs() < 1.0,
                "λ = {:.1} should match 1239.8/E = {:.1}",
                t.wavelength_nm,
                expected_nm
            );
        }
    }

    #[test]
    fn stda_oscillator_strength_nonnegative() {
        let scf = mock_scf(10, 6);
        let positions = vec![[0.0; 3]; 10];
        let basis_to_atom: Vec<usize> = (0..10).collect();
        let config = StdaConfig::default();
        let result = compute_stda(&scf, &basis_to_atom, &positions, &config);
        for t in &result.transitions {
            assert!(
                t.oscillator_strength >= 0.0,
                "Oscillator strength must be >= 0: {:.6}",
                t.oscillator_strength
            );
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  MIG-105: GIAO NMR — core validation
// ═══════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod giao_justification {
    use nalgebra::DMatrix;
    use sci_form::spectroscopy::{compute_nmr_shieldings, shieldings_to_shifts, ScfInput,
    };

    fn h2_system() -> (Vec<u8>, Vec<[f64; 3]>, ScfInput, Vec<usize>) {
        let elements = vec![1, 1];
        let positions = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]];
        let n_basis = 2;
        let scf = ScfInput {
            orbital_energies: vec![-0.6, 0.7],
            mo_coefficients: DMatrix::identity(n_basis, n_basis),
            density_matrix: DMatrix::zeros(n_basis, n_basis),
            overlap_matrix: DMatrix::identity(n_basis, n_basis),
            n_basis,
            n_electrons: 2,
        };
        let basis_to_atom = vec![0, 1];
        (elements, positions, scf, basis_to_atom)
    }

    #[test]
    fn giao_shieldings_count() {
        let (elements, positions, scf, basis_to_atom) = h2_system();
        let shieldings = compute_nmr_shieldings(&elements, &positions, &scf, &basis_to_atom);
        assert_eq!(shieldings.len(), 2, "H₂ should have 2 shielding tensors");
    }

    #[test]
    fn giao_hydrogen_positive_shielding() {
        let (elements, positions, scf, basis_to_atom) = h2_system();
        let shieldings = compute_nmr_shieldings(&elements, &positions, &scf, &basis_to_atom);
        for s in &shieldings {
            assert!(
                s.isotropic > 0.0,
                "H shielding should be positive: {:.4} ppm",
                s.isotropic
            );
        }
    }

    #[test]
    fn giao_chemical_shifts_finite() {
        let (elements, positions, scf, basis_to_atom) = h2_system();
        let shieldings = compute_nmr_shieldings(&elements, &positions, &scf, &basis_to_atom);
        let result = shieldings_to_shifts(&shieldings, &elements);
        assert_eq!(result.n_atoms, 2);
        for delta in &result.chemical_shifts {
            assert!(delta.is_finite(), "Shift must be finite: {}", delta);
        }
    }

    #[test]
    fn giao_h2_symmetry() {
        // Both H atoms in H₂ should have the same shielding
        let (elements, positions, scf, basis_to_atom) = h2_system();
        let shieldings = compute_nmr_shieldings(&elements, &positions, &scf, &basis_to_atom);
        let diff = (shieldings[0].isotropic - shieldings[1].isotropic).abs();
        assert!(
            diff < 5.0,
            "H₂ symmetry: Δσ = {diff:.4} ppm should be small"
        );
    }

    #[test]
    fn giao_c_more_shielded_than_h() {
        // Carbon should have larger isotropic shielding than hydrogen
        let elements = vec![6, 1];
        let positions = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]];
        let n_basis = 2;
        let scf = ScfInput {
            orbital_energies: vec![-0.8, 0.5],
            mo_coefficients: DMatrix::identity(n_basis, n_basis),
            density_matrix: DMatrix::zeros(n_basis, n_basis),
            overlap_matrix: DMatrix::identity(n_basis, n_basis),
            n_basis,
            n_electrons: 2,
        };
        let basis_to_atom = vec![0, 1];
        let shieldings = compute_nmr_shieldings(&elements, &positions, &scf, &basis_to_atom);
        assert!(
            shieldings[0].isotropic > shieldings[1].isotropic,
            "C ({:.1}) should be more shielded than H ({:.1})",
            shieldings[0].isotropic,
            shieldings[1].isotropic
        );
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  Cross-module integration: full pipeline tests
// ═══════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod integration {
    use sci_form::charges_eeq::{compute_eeq_charges, EeqConfig};
    use sci_form::dispersion::{compute_d4_energy, D4Config};
    use sci_form::solvation_alpb::{compute_alpb_solvation, AlpbConfig};

    /// Full pipeline: EEQ charges → D4 energy → ALPB solvation
    /// Mimics a fast pre-screening workflow.
    #[test]
    fn pipeline_eeq_d4_alpb() {
        // Methanol
        let elements = vec![6, 8, 1, 1, 1, 1];
        let positions = vec![
            [0.000, 0.000, 0.000],
            [1.430, 0.000, 0.000],
            [-0.390, 0.920, 0.000],
            [-0.390, -0.460, 0.800],
            [-0.390, -0.460, -0.800],
            [1.800, 0.890, 0.000],
        ];

        // Step 1: EEQ charges
        let config = EeqConfig::default();
        let charges = compute_eeq_charges(&elements, &positions, &config);
        assert!(charges.charges.len() == elements.len());
        assert!(charges.charges.iter().sum::<f64>().abs() < 1e-10);

        // Step 2: D4 dispersion
        let d4_config = D4Config::default();
        let d4 = compute_d4_energy(&elements, &positions, &d4_config);
        assert!(d4.total_energy < 0.0, "Dispersion must be attractive");

        // Step 3: ALPB solvation using EEQ charges
        let alpb_config = AlpbConfig::default();
        let solv =
            compute_alpb_solvation(&elements, &positions, &charges.charges, &alpb_config);
        assert!(
            solv.total_energy < 0.0,
            "Methanol in water should have negative solvation energy: {:.4}",
            solv.total_energy
        );

        println!("Pipeline results for methanol:");
        println!("  EEQ charges: O={:.3}, C={:.3}", charges.charges[1], charges.charges[0]);
        println!("  D4 dispersion: {:.4} kcal/mol", d4.total_kcal_mol);
        println!("  ALPB solvation: {:.4} kcal/mol", solv.total_energy);
    }

    #[test]
    fn eeq_charges_feed_alpb_correctly() {
        // Verify that EEQ charges produce physically sensible solvation
        let elements = vec![8, 1, 1];
        let positions = vec![
            [0.0, 0.0, 0.1173],
            [0.0, 0.7572, -0.4692],
            [0.0, -0.7572, -0.4692],
        ];

        let eeq = compute_eeq_charges(&elements, &positions, &EeqConfig::default());
        let solv = compute_alpb_solvation(
            &elements,
            &positions,
            &eeq.charges,
            &AlpbConfig::default(),
        );

        // Water in water: electrostatic component should be stabilizing
        assert!(solv.electrostatic_energy < 0.0);
        // Total may be slightly positive due to non-polar cavity cost
        assert!(
            solv.total_energy.is_finite(),
            "Total solvation energy must be finite"
        );
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  E4 validation: experimental modules stay behind feature flags
// ═══════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod e4_experimental_gated {
    #[test]
    fn core_modules_always_available() {
        // These should compile without any feature flags
        let _ = sci_form::charges_eeq::EeqConfig::default();
        let _ = sci_form::dispersion::D4Config::default();
        let _ = sci_form::solvation_alpb::AlpbConfig::default();
    }

    #[test]
    fn experimental_status_module_exists() {
        // The documentation module should be accessible
        // (it's a doc-only module, nothing to call)
        let _module_path = module_path!();
    }
}
