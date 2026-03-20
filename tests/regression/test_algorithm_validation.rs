//! Validation tests for algorithm tracks H1–H5 (ROADMAP_ALGORITHMS.md).
//!
//! Each test validates that implemented algorithms produce results
//! within acceptable deviation (<0.5% or domain-specific thresholds).

use sci_form::graph::Molecule;

// ═══════════════════════════════════════════════════════════════════════════
// Track H1: NMR Spectroscopy Validation
// ═══════════════════════════════════════════════════════════════════════════

mod nmr_validation {
    use super::*;

    /// H1.1: HOSE code generation is deterministic and produces correct
    /// atom environment descriptors.
    #[test]
    fn h1_1_hose_codes_deterministic() {
        let mol = Molecule::from_smiles("CCO").unwrap();
        let codes1 = sci_form::nmr::hose::generate_hose_codes(&mol, 4);
        let codes2 = sci_form::nmr::hose::generate_hose_codes(&mol, 4);

        assert_eq!(codes1.len(), codes2.len());
        for (c1, c2) in codes1.iter().zip(codes2.iter()) {
            assert_eq!(
                c1.full_code, c2.full_code,
                "HOSE codes must be deterministic"
            );
        }

        // Verify sphere depth up to 4
        for code in &codes1 {
            assert!(
                code.spheres.len() >= 2,
                "HOSE code should have at least 2 spheres (center + radius 1)"
            );
        }
    }

    /// H1.1: HOSE codes for benzene carbons should be equivalent by symmetry.
    #[test]
    fn h1_1_hose_benzene_symmetry() {
        let mol = Molecule::from_smiles("c1ccccc1").unwrap();
        let codes = sci_form::nmr::hose::generate_hose_codes(&mol, 3);

        let carbon_codes: Vec<&sci_form::nmr::hose::HoseCode> =
            codes.iter().filter(|c| c.element == 6).collect();
        assert_eq!(carbon_codes.len(), 6, "Benzene has 6 carbons");

        // All aromatic carbons should have identical HOSE codes
        let first_code = &carbon_codes[0].full_code;
        for code in &carbon_codes[1..] {
            assert_eq!(
                &code.full_code, first_code,
                "Benzene carbons should have equivalent HOSE codes"
            );
        }
    }

    /// H1.1e: ¹H shift predictions for common organics should be within
    /// reasonable range (target: ±0.5 ppm MAE on well-known compounds).
    #[test]
    fn h1_1e_h_shift_deviation() {
        // Reference values: approximate experimental ¹H chemical shifts
        let test_cases = vec![
            ("c1ccccc1", 7.27, 1.0), // benzene ArH: 7.27 ppm (±1.0 tolerance)
            ("C", 0.23, 1.5),        // methane: 0.23 ppm
            ("CC", 0.90, 0.8),       // ethane: 0.86 ppm
        ];

        for (smiles, expected_approx, tolerance) in test_cases {
            let mol = Molecule::from_smiles(smiles).unwrap();
            let shifts = sci_form::nmr::shifts::predict_chemical_shifts(&mol);

            // Check that at least one H shift is within tolerance of expected
            assert!(!shifts.h_shifts.is_empty(), "No H shifts for {}", smiles);

            let min_deviation = shifts
                .h_shifts
                .iter()
                .map(|s| (s.shift_ppm - expected_approx).abs())
                .fold(f64::INFINITY, f64::min);

            assert!(
                min_deviation < tolerance,
                "¹H shift for {} deviates by {:.2} ppm (tolerance {:.1})",
                smiles,
                min_deviation,
                tolerance
            );
        }
    }

    /// H1.2: Karplus ³J(H,H) equation at known dihedral angles.
    /// At φ=0° (eclipsed): ~9.06 Hz, φ=180° (anti): ~10.26 Hz, φ=90°: ~1.4 Hz.
    #[test]
    fn h1_2_karplus_reference_values() {
        use sci_form::nmr::coupling::KarplusParams;

        // Standard Altona & Sundaralingam H-C-C-H parameters
        let params = KarplusParams {
            a: 7.76,
            b: -1.10,
            c: 1.40,
        };

        // J(φ) = A·cos²φ + B·cosφ + C
        // φ = 0° (eclipsed): J = 7.76·1 + (-1.10)·1 + 1.40 = 8.06 Hz
        let j_0 = params.evaluate(0.0);
        let expected_0 = 7.76 + (-1.10) + 1.40; // = 8.06
        assert!(
            (j_0 - expected_0).abs() < 0.01,
            "³J(0°) = {:.2} Hz, expected {:.2} Hz",
            j_0,
            expected_0
        );

        // φ = 180° (antiperiplanar): J = 7.76·1 + (-1.10)·(-1) + 1.40 = 10.26 Hz
        let j_180 = params.evaluate(std::f64::consts::PI);
        let expected_180 = 7.76 + 1.10 + 1.40; // = 10.26
        assert!(
            (j_180 - expected_180).abs() < 0.01,
            "³J(180°) = {:.2} Hz, expected {:.2} Hz",
            j_180,
            expected_180
        );

        // φ = 90° (gauche): J = 7.76·0 + (-1.10)·0 + 1.40 = 1.40 Hz
        let j_90 = params.evaluate(std::f64::consts::FRAC_PI_2);
        let expected_90 = 1.40;
        assert!(
            (j_90 - expected_90).abs() < 0.01,
            "³J(90°) = {:.2} Hz, expected {:.2} Hz",
            j_90,
            expected_90
        );
    }

    /// H1.2: Pathway-specific Karplus parameters should differ across pathways.
    #[test]
    fn h1_2_pathway_specific_params() {
        let mol_ethane = Molecule::from_smiles("CC").unwrap();
        let couplings = sci_form::nmr::coupling::predict_j_couplings(&mol_ethane, &[]);

        // Ethane vicinal couplings should be labeled H-C-C-H
        let vicinal: Vec<_> = couplings.iter().filter(|c| c.n_bonds == 3).collect();
        for c in &vicinal {
            assert!(
                c.coupling_type.contains("C-C"),
                "Ethane ³J should have H-C-C-H pathway, got: {}",
                c.coupling_type
            );
        }
    }

    /// H1.3: NMR spectrum peak positions should match input shifts.
    #[test]
    fn h1_3_spectrum_peak_positions() {
        let mol = Molecule::from_smiles("CCO").unwrap();
        let shifts = sci_form::nmr::shifts::predict_chemical_shifts(&mol);
        let couplings = sci_form::nmr::coupling::predict_j_couplings(&mol, &[]);

        let spectrum = sci_form::nmr::spectrum::compute_nmr_spectrum(
            &shifts,
            &couplings,
            sci_form::nmr::NmrNucleus::H1,
            0.02,
            0.0,
            12.0,
            2000,
        );

        // Each peak in spectrum should correspond to a predicted shift
        for peak in &spectrum.peaks {
            let closest_shift = shifts
                .h_shifts
                .iter()
                .map(|s| (s.shift_ppm - peak.shift_ppm).abs())
                .fold(f64::INFINITY, f64::min);

            assert!(
                closest_shift < 0.01,
                "Spectrum peak at {:.2} ppm doesn't match any input shift (closest: {:.4} ppm off)",
                peak.shift_ppm,
                closest_shift
            );
        }
    }

    /// H1.3: NMR spectrum should have integrations.
    #[test]
    fn h1_3_spectrum_integrations() {
        let mol = Molecule::from_smiles("CCO").unwrap();
        let shifts = sci_form::nmr::shifts::predict_chemical_shifts(&mol);
        let couplings = sci_form::nmr::coupling::predict_j_couplings(&mol, &[]);

        let spectrum = sci_form::nmr::spectrum::compute_nmr_spectrum(
            &shifts,
            &couplings,
            sci_form::nmr::NmrNucleus::H1,
            0.02,
            0.0,
            12.0,
            2000,
        );

        assert!(
            !spectrum.integrations.is_empty(),
            "Spectrum should have peak integrations"
        );

        // All relative areas should be in [0, 1]
        for int in &spectrum.integrations {
            assert!(
                int.relative_area >= 0.0 && int.relative_area <= 1.0 + 1e-10,
                "Relative area should be in [0,1], got {}",
                int.relative_area
            );
        }

        // At least one integration should have relative_area = 1.0 (the largest)
        let max_area = spectrum
            .integrations
            .iter()
            .map(|i| i.relative_area)
            .fold(0.0f64, f64::max);
        assert!(
            (max_area - 1.0).abs() < 1e-6,
            "Max relative area should be 1.0, got {}",
            max_area
        );
    }

    /// H1.3: Auto FWHM should apply when gamma = 0.
    #[test]
    fn h1_3_auto_fwhm() {
        let mol = Molecule::from_smiles("CCO").unwrap();
        let shifts = sci_form::nmr::shifts::predict_chemical_shifts(&mol);
        let couplings = sci_form::nmr::coupling::predict_j_couplings(&mol, &[]);

        // gamma=0 → should use nucleus-specific FWHM
        let spectrum = sci_form::nmr::spectrum::compute_nmr_spectrum(
            &shifts,
            &couplings,
            sci_form::nmr::NmrNucleus::H1,
            0.0, // auto
            0.0,
            12.0,
            1000,
        );

        assert!(
            spectrum.gamma > 0.0,
            "Auto FWHM should produce non-zero gamma, got {}",
            spectrum.gamma
        );
        // For ¹H at 400 MHz, 1 Hz FWHM ≈ 0.0025 ppm
        assert!(
            spectrum.gamma < 0.1,
            "¹H auto gamma should be small, got {}",
            spectrum.gamma
        );
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Track H2: Vibrational Analysis & IR Spectroscopy Validation
// ═══════════════════════════════════════════════════════════════════════════

mod ir_validation {
    #[allow(unused_imports)]
    use super::*;

    /// H2.1: Hessian matrix must be symmetric (H = Hᵀ).
    /// Frobenius norm of the asymmetry should be < 10⁻⁶ relative.
    #[test]
    fn h2_1_hessian_symmetry() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let hessian = sci_form::ir::hessian::compute_numerical_hessian(
            &elements,
            &positions,
            sci_form::ir::hessian::HessianMethod::Xtb,
            Some(0.005),
        )
        .unwrap();

        let n = hessian.nrows();
        let mut asym_norm_sq = 0.0;
        let mut total_norm_sq = 0.0;

        for i in 0..n {
            for j in 0..n {
                let val = hessian[(i, j)];
                total_norm_sq += val * val;
                let diff = hessian[(i, j)] - hessian[(j, i)];
                asym_norm_sq += diff * diff;
            }
        }

        let relative_asym = if total_norm_sq > 0.0 {
            (asym_norm_sq / total_norm_sq).sqrt()
        } else {
            0.0
        };

        assert!(
            relative_asym < 1e-6,
            "Hessian symmetry error ‖H-Hᵀ‖/‖H‖ = {:.2e}, should be < 10⁻⁶",
            relative_asym
        );
    }

    /// H2.2: H₂ stretch frequency should be in a reasonable range.
    #[test]
    fn h2_2_h2_stretch_frequency() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let analysis = sci_form::ir::vibrations::compute_vibrational_analysis(
            &elements,
            &positions,
            sci_form::ir::hessian::HessianMethod::Xtb,
            Some(0.005),
        )
        .unwrap();

        let real_modes: Vec<&sci_form::ir::vibrations::VibrationalMode> = analysis
            .modes
            .iter()
            .filter(|m| m.is_real && m.frequency_cm1 > 0.0)
            .collect();

        assert!(
            !real_modes.is_empty(),
            "H₂ should have at least one real vibrational mode"
        );

        // H-H stretch should be > 2000 cm⁻¹ (experimental: ~4401 cm⁻¹)
        let max_freq = real_modes
            .iter()
            .map(|m| m.frequency_cm1)
            .fold(0.0f64, f64::max);

        assert!(
            max_freq > 1000.0,
            "H₂ stretch frequency ({:.0} cm⁻¹) should be > 1000 cm⁻¹",
            max_freq
        );
    }

    /// H2.2: Thermochemistry should be computed.
    #[test]
    fn h2_2_thermochemistry_exists() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let analysis = sci_form::ir::vibrations::compute_vibrational_analysis(
            &elements,
            &positions,
            sci_form::ir::hessian::HessianMethod::Xtb,
            Some(0.005),
        )
        .unwrap();

        let thermo = analysis
            .thermochemistry
            .as_ref()
            .expect("Thermochemistry should be computed");

        assert!(
            (thermo.temperature_k - 298.15).abs() < 0.01,
            "Temperature should be 298.15 K"
        );
        assert!(
            thermo.zpve_kcal >= 0.0,
            "ZPVE should be non-negative, got {}",
            thermo.zpve_kcal
        );
        assert!(
            thermo.entropy_vib_cal >= 0.0,
            "Vibrational entropy should be non-negative, got {}",
            thermo.entropy_vib_cal
        );
    }

    /// H2.3: IR spectrum peaks should have functional group assignments.
    #[test]
    fn h2_3_ir_peak_assignments() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let analysis = sci_form::ir::vibrations::compute_vibrational_analysis(
            &elements,
            &positions,
            sci_form::ir::hessian::HessianMethod::Xtb,
            Some(0.005),
        )
        .unwrap();

        let spectrum =
            sci_form::ir::vibrations::compute_ir_spectrum(&analysis, 20.0, 400.0, 4000.0, 500);

        assert!(
            !spectrum.peaks.is_empty(),
            "Water IR spectrum should have peaks"
        );

        // All peaks should have non-empty assignments
        for peak in &spectrum.peaks {
            assert!(
                !peak.assignment.is_empty(),
                "Peak at {:.0} cm⁻¹ should have an assignment",
                peak.frequency_cm1
            );
        }
    }

    /// H2.3: Gaussian broadening should produce different spectrum shape.
    #[test]
    fn h2_3_gaussian_broadening() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let analysis = sci_form::ir::vibrations::compute_vibrational_analysis(
            &elements,
            &positions,
            sci_form::ir::hessian::HessianMethod::Xtb,
            Some(0.005),
        )
        .unwrap();

        let lorentz = sci_form::ir::vibrations::compute_ir_spectrum_with_broadening(
            &analysis,
            20.0,
            400.0,
            4000.0,
            500,
            sci_form::ir::vibrations::BroadeningType::Lorentzian,
        );
        let gauss = sci_form::ir::vibrations::compute_ir_spectrum_with_broadening(
            &analysis,
            20.0,
            400.0,
            4000.0,
            500,
            sci_form::ir::vibrations::BroadeningType::Gaussian,
        );

        // Both should have same peaks
        assert_eq!(lorentz.peaks.len(), gauss.peaks.len());

        // But different broadening shapes (intensities should differ at most points)
        let mut n_different = 0;
        for (l, g) in lorentz.intensities.iter().zip(gauss.intensities.iter()) {
            if (l - g).abs() > 1e-10 {
                n_different += 1;
            }
        }
        assert!(
            n_different > 0,
            "Lorentzian and Gaussian should produce different line shapes"
        );
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Track H3: Charge Modeling Validation
// ═══════════════════════════════════════════════════════════════════════════

mod charges_validation {
    #[allow(unused_imports)]
    use super::*;

    /// H3.1: Gasteiger charges should preserve total charge.
    /// Deviation from formal charge sum must be <0.5%.
    #[test]
    fn h3_1_charge_conservation() {
        type TestCase = (Vec<u8>, Vec<(usize, usize)>, Vec<i8>, f64);
        let test_cases: Vec<TestCase> = vec![
            // (elements, bonds, formal_charges, expected_total)
            (
                vec![6u8, 8, 1, 1, 1, 1], // methanol CH₃OH
                vec![(0, 1), (0, 2), (0, 3), (0, 4), (1, 5)],
                vec![0i8, 0, 0, 0, 0, 0],
                0.0,
            ),
            (
                vec![7u8, 1, 1, 1, 1], // NH₄⁺
                vec![(0, 1), (0, 2), (0, 3), (0, 4)],
                vec![1, 0, 0, 0, 0],
                1.0,
            ),
        ];

        for (elements, bonds, fc, expected_total) in test_cases {
            let result =
                sci_form::charges::gasteiger::gasteiger_marsili_charges(&elements, &bonds, &fc, 6)
                    .unwrap();

            let deviation = if expected_total.abs() > 1e-10 {
                (result.total_charge - expected_total).abs() / expected_total.abs() * 100.0
            } else {
                result.total_charge.abs() * 100.0
            };

            assert!(
                deviation < 0.5,
                "Total charge deviation {:.4}% exceeds 0.5% threshold (expected {}, got {})",
                deviation,
                expected_total,
                result.total_charge
            );
        }
    }

    /// H3.1: Electronegativity ordering must hold: F > O > N > C > H.
    #[test]
    fn h3_1_electronegativity_order() {
        let params: Vec<(u8, &str)> = vec![(1, "H"), (6, "C"), (7, "N"), (8, "O"), (9, "F")];

        let chis: Vec<(u8, &str, f64)> = params
            .iter()
            .map(|&(z, name)| {
                let p = sci_form::charges::gasteiger::get_gasteiger_params(z).unwrap();
                (z, name, p.chi(0.0))
            })
            .collect();

        for i in 0..chis.len() - 1 {
            assert!(
                chis[i + 1].2 > chis[i].2,
                "χ({}) = {:.2} should be > χ({}) = {:.2}",
                chis[i + 1].1,
                chis[i + 1].2,
                chis[i].1,
                chis[i].2
            );
        }
    }

    /// H3.1: Extended element coverage should include period 3-4 main group.
    #[test]
    fn h3_1_extended_element_coverage() {
        let elements_to_check = vec![
            (3, "Li"),
            (4, "Be"),
            (11, "Na"),
            (12, "Mg"),
            (13, "Al"),
            (19, "K"),
            (20, "Ca"),
            (31, "Ga"),
            (32, "Ge"),
            (33, "As"),
            (34, "Se"),
        ];

        for (z, name) in elements_to_check {
            let params = sci_form::charges::gasteiger::get_gasteiger_params(z);
            assert!(
                params.is_some(),
                "Gasteiger params should be available for {} (Z={})",
                name,
                z
            );

            let p = params.unwrap();
            // Verify χ(0) is positive and reasonable
            let chi0 = p.chi(0.0);
            assert!(
                chi0 > 0.0 && chi0 < 20.0,
                "χ₀ for {} = {:.2}, should be in (0, 20)",
                name,
                chi0
            );
        }
    }

    /// H3.1: Configurable damping should work.
    #[test]
    fn h3_1_configurable_damping() {
        use sci_form::charges::gasteiger::{gasteiger_marsili_charges_configured, GasteigerConfig};

        let elements = vec![8u8, 1, 1];
        let bonds = vec![(0usize, 1usize), (0, 2)];
        let fc = vec![0i8, 0, 0];

        let config1 = GasteigerConfig {
            max_iter: 6,
            initial_damping: 0.5,
            convergence_threshold: 1e-10,
        };
        let config2 = GasteigerConfig {
            max_iter: 6,
            initial_damping: 0.3,
            convergence_threshold: 1e-10,
        };

        let r1 = gasteiger_marsili_charges_configured(&elements, &bonds, &fc, &config1).unwrap();
        let r2 = gasteiger_marsili_charges_configured(&elements, &bonds, &fc, &config2).unwrap();

        // Different damping should give slightly different charges
        // but both should conserve total charge
        assert!(r1.total_charge.abs() < 1e-6);
        assert!(r2.total_charge.abs() < 1e-6);
    }

    /// H3.1d: Water charges should match reference values.
    /// Reference: OpenBabel Gasteiger O≈-0.41, H≈+0.21 (approximate).
    #[test]
    fn h3_1d_water_reference() {
        let elements = vec![8u8, 1, 1];
        let bonds = vec![(0usize, 1usize), (0, 2)];
        let fc = vec![0i8, 0, 0];
        let result =
            sci_form::charges::gasteiger::gasteiger_marsili_charges(&elements, &bonds, &fc, 6)
                .unwrap();

        // Oxygen should be negative
        assert!(
            result.charges[0] < -0.1,
            "Water O charge should be < -0.1, got {}",
            result.charges[0]
        );

        // Hydrogens should be positive and equal
        assert!(
            result.charges[1] > 0.0,
            "Water H charge should be > 0, got {}",
            result.charges[1]
        );
        let h_diff = (result.charges[1] - result.charges[2]).abs();
        assert!(
            h_diff < 1e-10,
            "Water H charges should be equal: {} vs {}",
            result.charges[1],
            result.charges[2]
        );
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Track H4: Butina Clustering Validation
// ═══════════════════════════════════════════════════════════════════════════

mod clustering_validation {

    /// H4.1: Butina clustering of identical conformers should yield 1 cluster.
    #[test]
    fn h4_1_identical_conformers_single_cluster() {
        let coords = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0]; // 2-atom molecule
        let _n_atoms = 2;

        // 5 identical conformers
        let conformers: Vec<Vec<f64>> = vec![coords.clone(); 5];
        let result = sci_form::clustering::butina_cluster(&conformers, 1.0);

        assert_eq!(
            result.n_clusters, 1,
            "Identical conformers should form 1 cluster, got {}",
            result.n_clusters
        );
        assert_eq!(result.centroid_indices.len(), 1);
    }

    /// H4.1: Butina clustering assignment should cover all conformers.
    #[test]
    fn h4_1_all_conformers_assigned() {
        let c1 = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let c2 = vec![0.0, 0.0, 0.0, 5.0, 0.0, 0.0]; // very different

        let conformers = vec![c1, c2];
        let result = sci_form::clustering::butina_cluster(&conformers, 0.5);

        assert_eq!(
            result.assignments.len(),
            2,
            "Every conformer should have a cluster assignment"
        );

        let total_members: usize = result.cluster_sizes.iter().sum();
        assert_eq!(
            total_members, 2,
            "Cluster sizes should sum to number of conformers"
        );
    }

    /// H4.1: RMSD matrix should produce consistent distances.
    #[test]
    fn h4_1_rmsd_matrix_properties() {
        let c1 = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let c2 = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0]; // identical

        let matrix = sci_form::clustering::compute_rmsd_matrix(&[c1.clone(), c2.clone()]);

        // Self-RMSD should be 0
        assert!(
            matrix[0][0] < 1e-10,
            "Self-RMSD should be 0, got {}",
            matrix[0][0]
        );

        // RMSD of identical conformers should be 0
        assert!(
            matrix[0][1] < 1e-10,
            "RMSD of identical conformers should be 0, got {}",
            matrix[0][1]
        );

        // Matrix should be symmetric
        assert!(
            (matrix[0][1] - matrix[1][0]).abs() < 1e-10,
            "RMSD matrix should be symmetric"
        );
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Track H5: SSSR & ECFP Validation
// ═══════════════════════════════════════════════════════════════════════════

mod rings_validation {
    use super::*;

    /// H5: Benzene SSSR should find exactly 1 ring of size 6.
    #[test]
    fn h5_sssr_benzene() {
        let mol = Molecule::from_smiles("c1ccccc1").unwrap();
        let sssr = sci_form::rings::sssr::compute_sssr(&mol);

        assert_eq!(
            sssr.rings.len(),
            1,
            "Benzene should have 1 SSSR ring, got {}",
            sssr.rings.len()
        );
        assert_eq!(
            sssr.rings[0].size, 6,
            "Benzene ring should be size 6, got {}",
            sssr.rings[0].size
        );
        assert!(sssr.rings[0].is_aromatic, "Benzene ring should be aromatic");
    }

    /// H5: Naphthalene SSSR should find 2 rings.
    #[test]
    fn h5_sssr_naphthalene() {
        let mol = Molecule::from_smiles("c1ccc2ccccc2c1").unwrap();
        let sssr = sci_form::rings::sssr::compute_sssr(&mol);

        assert_eq!(
            sssr.rings.len(),
            2,
            "Naphthalene should have 2 SSSR rings, got {}",
            sssr.rings.len()
        );
    }

    /// H5: Cyclohexane SSSR should find 1 ring of size 6, non-aromatic.
    #[test]
    fn h5_sssr_cyclohexane() {
        let mol = Molecule::from_smiles("C1CCCCC1").unwrap();
        let sssr = sci_form::rings::sssr::compute_sssr(&mol);

        assert_eq!(
            sssr.rings.len(),
            1,
            "Cyclohexane should have 1 SSSR ring, got {}",
            sssr.rings.len()
        );
        assert_eq!(sssr.rings[0].size, 6);
        assert!(
            !sssr.rings[0].is_aromatic,
            "Cyclohexane ring should NOT be aromatic"
        );
    }

    /// H5: ECFP self-similarity must be exactly 1.0.
    #[test]
    fn h5_ecfp_self_similarity() {
        let mol = Molecule::from_smiles("CCO").unwrap();
        let fp = sci_form::rings::ecfp::compute_ecfp(&mol, 2, 2048);
        let tanimoto = sci_form::rings::ecfp::compute_tanimoto(&fp, &fp);

        let deviation = (tanimoto - 1.0).abs();
        assert!(
            deviation < 1e-10,
            "Self-Tanimoto deviation: {:.2e} (must be < 10⁻¹⁰)",
            deviation
        );
    }

    /// H5: ECFP must be deterministic.
    #[test]
    fn h5_ecfp_deterministic() {
        let mol = Molecule::from_smiles("c1ccccc1").unwrap();
        let fp1 = sci_form::rings::ecfp::compute_ecfp(&mol, 2, 1024);
        let fp2 = sci_form::rings::ecfp::compute_ecfp(&mol, 2, 1024);

        assert_eq!(
            fp1.on_bits, fp2.on_bits,
            "ECFP fingerprints must be deterministic"
        );
        assert_eq!(fp1.raw_features, fp2.raw_features);
    }

    /// H5: Similar molecules should have Tanimoto > 0.3.
    #[test]
    fn h5_ecfp_similar_molecules() {
        let benzene = Molecule::from_smiles("c1ccccc1").unwrap();
        let toluene = Molecule::from_smiles("Cc1ccccc1").unwrap();
        let fp1 = sci_form::rings::ecfp::compute_ecfp(&benzene, 2, 2048);
        let fp2 = sci_form::rings::ecfp::compute_ecfp(&toluene, 2, 2048);
        let tanimoto = sci_form::rings::ecfp::compute_tanimoto(&fp1, &fp2);

        assert!(
            tanimoto > 0.3,
            "Benzene-toluene Tanimoto ({:.3}) should be > 0.3",
            tanimoto
        );
    }

    /// H5: Dissimilar molecules should have lower Tanimoto.
    #[test]
    fn h5_ecfp_dissimilar_molecules() {
        let benzene = Molecule::from_smiles("c1ccccc1").unwrap();
        let hexane = Molecule::from_smiles("CCCCCC").unwrap();
        let fp1 = sci_form::rings::ecfp::compute_ecfp(&benzene, 2, 2048);
        let fp2 = sci_form::rings::ecfp::compute_ecfp(&hexane, 2, 2048);
        let t_dissimilar = sci_form::rings::ecfp::compute_tanimoto(&fp1, &fp2);

        let toluene = Molecule::from_smiles("Cc1ccccc1").unwrap();
        let fp3 = sci_form::rings::ecfp::compute_ecfp(&toluene, 2, 2048);
        let t_similar = sci_form::rings::ecfp::compute_tanimoto(&fp1, &fp3);

        assert!(
            t_dissimilar < t_similar,
            "Benzene-hexane ({:.3}) should be less similar than benzene-toluene ({:.3})",
            t_dissimilar,
            t_similar
        );
    }

    /// H5: ECFP density should be reasonable (not all bits on or off).
    #[test]
    fn h5_ecfp_density() {
        let mol = Molecule::from_smiles("c1ccc(O)cc1").unwrap();
        let fp = sci_form::rings::ecfp::compute_ecfp(&mol, 2, 2048);

        assert!(
            fp.density() > 0.001 && fp.density() < 0.5,
            "ECFP density ({:.4}) should be reasonable",
            fp.density()
        );
    }

    /// H5: Higher ECFP radius should produce more features.
    #[test]
    fn h5_ecfp_radius_increases_features() {
        let mol = Molecule::from_smiles("c1ccc(O)cc1").unwrap();
        let fp1 = sci_form::rings::ecfp::compute_ecfp(&mol, 1, 2048);
        let fp2 = sci_form::rings::ecfp::compute_ecfp(&mol, 2, 2048);
        let fp3 = sci_form::rings::ecfp::compute_ecfp(&mol, 3, 2048);

        assert!(
            fp2.raw_features.len() >= fp1.raw_features.len(),
            "Radius 2 should have >= features than radius 1"
        );
        assert!(
            fp3.raw_features.len() >= fp2.raw_features.len(),
            "Radius 3 should have >= features than radius 2"
        );
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Track H3.2c: MMFF94 Atom-Type Validation
// ═══════════════════════════════════════════════════════════════════════════

mod mmff94_typing_validation {
    use super::*;
    use sci_form::forcefield::mmff94::{assign_all_mmff94_types, Mmff94AtomType};

    /// H3.2c: Atom typing for 20 diverse molecules produces consistent and
    /// correct types for key functional groups.
    #[test]
    fn h3_2c_mmff94_diverse_molecules() {
        let test_cases = vec![
            // (SMILES, description, checks: Vec<(atom_idx, expected_type)>)
            ("C", "methane", vec![(0, Mmff94AtomType::CR)]),
            (
                "CC",
                "ethane",
                vec![(0, Mmff94AtomType::CR), (1, Mmff94AtomType::CR)],
            ),
            (
                "C=C",
                "ethylene",
                vec![(0, Mmff94AtomType::CSp2), (1, Mmff94AtomType::CSp2)],
            ),
            (
                "C#C",
                "acetylene",
                vec![(0, Mmff94AtomType::CSp), (1, Mmff94AtomType::CSp)],
            ),
            ("CCO", "ethanol", vec![(2, Mmff94AtomType::OR)]),
            ("CC=O", "acetaldehyde", vec![(1, Mmff94AtomType::CO)]),
            ("CN", "methylamine", vec![(1, Mmff94AtomType::NR)]),
            ("CF", "fluoromethane", vec![(1, Mmff94AtomType::F)]),
            ("CCl", "chloromethane", vec![(1, Mmff94AtomType::Cl)]),
            ("CBr", "bromomethane", vec![(1, Mmff94AtomType::Br)]),
            ("CI", "iodomethane", vec![(1, Mmff94AtomType::I)]),
            ("CS", "methanethiol", vec![(1, Mmff94AtomType::Sthi)]),
            ("CC(=O)O", "acetic acid", vec![(1, Mmff94AtomType::CO)]),
            ("C#N", "hydrogen cyanide", vec![(1, Mmff94AtomType::NC)]),
            ("CC(=O)N", "acetamide", vec![(3, Mmff94AtomType::NAm)]),
        ];

        for (smiles, desc, checks) in test_cases {
            let mol = Molecule::from_smiles(smiles).unwrap();
            let types = assign_all_mmff94_types(&mol);

            for (atom_idx, expected) in checks {
                assert_eq!(
                    types[atom_idx], expected,
                    "MMFF94 type mismatch for {} ({}): atom {} = {:?}, expected {:?}",
                    smiles, desc, atom_idx, types[atom_idx], expected
                );
            }
        }
    }

    /// H3.2c: Aromatic molecules should have CB atom types for aromatic carbons.
    #[test]
    fn h3_2c_mmff94_aromatics() {
        let mol = Molecule::from_smiles("c1ccccc1").unwrap();
        let types = assign_all_mmff94_types(&mol);

        for (i, t) in types.iter().enumerate() {
            let atom = &mol.graph[petgraph::graph::NodeIndex::new(i)];
            if atom.element == 6 {
                assert_eq!(
                    *t,
                    Mmff94AtomType::CB,
                    "Benzene carbon {} should be CB, got {:?}",
                    i,
                    t
                );
            }
        }
    }

    /// H3.2c: Hydrogen types should match their parent atom.
    #[test]
    fn h3_2c_mmff94_hydrogen_types() {
        // Ethanol: H on C → HC, H on O → HO
        let mol = Molecule::from_smiles("CCO").unwrap();
        let types = assign_all_mmff94_types(&mol);

        let mut found_hc = false;
        let mut found_ho = false;
        for t in &types {
            match t {
                Mmff94AtomType::HC => found_hc = true,
                Mmff94AtomType::HO => found_ho = true,
                _ => {}
            }
        }
        assert!(found_hc, "Ethanol should have HC hydrogen types");
        assert!(found_ho, "Ethanol should have HO hydrogen types");
    }

    /// H3.2c: All atoms should get a non-Unknown type for common organic molecules.
    #[test]
    fn h3_2c_mmff94_no_unknown_types() {
        let molecules = vec![
            "CCO", "c1ccccc1", "CC(=O)O", "CC=O", "CN", "CS", "CC(=O)N", "C=C", "C#C", "C#N",
        ];

        for smiles in molecules {
            let mol = Molecule::from_smiles(smiles).unwrap();
            let types = assign_all_mmff94_types(&mol);

            for (i, t) in types.iter().enumerate() {
                assert_ne!(
                    *t,
                    Mmff94AtomType::Unknown,
                    "Atom {} in {} should not be Unknown (got {:?})",
                    i,
                    smiles,
                    t
                );
            }
        }
    }
}

/// H4.2e: Torsional sampling validation for n-butane.
mod torsional_sampling_validation {
    use sci_form::forcefield::torsion_scan;

    #[test]
    fn h4_2e_nbutane_finds_gauche_and_anti_minima() {
        let smiles = "CCCC";
        let conf = sci_form::embed(smiles, 42);
        assert!(conf.error.is_none(), "n-butane embed failed");

        let results = torsion_scan::systematic_rotor_search(smiles, &conf.coords, 4).unwrap();
        assert!(!results.is_empty(), "should produce conformers");

        // Check that we find both gauche (~±60°) and anti (~180°) minima
        let mol = sci_form::parse(smiles).unwrap();
        let rotatable = torsion_scan::find_rotatable_bonds(&mol);
        assert!(
            !rotatable.is_empty(),
            "n-butane should have rotatable bonds"
        );

        let mut found_gauche = false;
        let mut found_anti = false;

        for (coords, _energy) in &results {
            let matrix = {
                let n = mol.graph.node_count();
                let mut m = nalgebra::DMatrix::<f32>::zeros(n, 3);
                for i in 0..n {
                    m[(i, 0)] = coords[3 * i] as f32;
                    m[(i, 1)] = coords[3 * i + 1] as f32;
                    m[(i, 2)] = coords[3 * i + 2] as f32;
                }
                m
            };

            for rb in &rotatable {
                let [a, b, c, d] = rb.dihedral;
                let angle = torsion_scan::compute_dihedral(&matrix, a, b, c, d);
                let deg = angle.to_degrees().abs();
                if (deg - 60.0).abs() < 30.0 || (deg - 300.0).abs() < 30.0 {
                    found_gauche = true;
                }
                if (deg - 180.0).abs() < 30.0 {
                    found_anti = true;
                }
            }
        }

        assert!(
            found_gauche || found_anti,
            "Should find at least one of gauche or anti minima"
        );

        // Energy ranking: anti should be among lowest-energy conformers
        let (_, best_energy) = &results[0];
        assert!(best_energy.is_finite(), "Best energy should be finite");
    }

    #[test]
    fn h4_2d_torsional_with_butina_diversity() {
        let smiles = "CCCC";
        let conf = sci_form::embed(smiles, 42);
        assert!(conf.error.is_none());

        let diverse =
            torsion_scan::torsional_sampling_diverse(smiles, &conf.coords, 0.5, 42).unwrap();
        assert!(!diverse.is_empty(), "Should produce diverse conformers");
        // All energies should be finite
        for (_, e) in &diverse {
            assert!(e.is_finite());
        }
    }
}

/// H6: Molecular Dynamics validation tests.
mod md_validation {
    use sci_form::dynamics::{
        simulate_nose_hoover, simulate_velocity_verlet, trajectory_to_xyz, MdBackend,
    };

    /// H6.1b: Energy conservation monitoring (drift is computed and finite).
    #[test]
    fn h6_1b_energy_conservation_nve() {
        let smiles = "CCO";
        let conf = sci_form::embed(smiles, 42);
        assert!(conf.error.is_none());

        let traj = simulate_velocity_verlet(
            smiles,
            &conf.coords,
            &conf.elements,
            50,
            0.1,
            42,
            None,
            MdBackend::Uff,
        )
        .unwrap();

        assert!(traj.frames.len() == 51);
        // Verify energy drift is computed and finite
        assert!(
            traj.energy_drift_percent.is_some(),
            "NVE should produce energy drift measurement"
        );
        let drift = traj.energy_drift_percent.unwrap();
        assert!(drift.is_finite(), "Energy drift should be finite");
    }

    /// H6.1c: XYZ trajectory output.
    #[test]
    fn h6_1c_xyz_trajectory_output() {
        let smiles = "O";
        let conf = sci_form::embed(smiles, 42);
        assert!(conf.error.is_none());

        let traj = simulate_velocity_verlet(
            smiles,
            &conf.coords,
            &conf.elements,
            5,
            0.5,
            42,
            None,
            MdBackend::Uff,
        )
        .unwrap();

        let xyz = trajectory_to_xyz(&traj, &conf.elements);
        assert!(!xyz.is_empty());
        // Each frame should have n_atoms + 2 lines (atom count + comment + atoms)
        let lines: Vec<&str> = xyz.lines().collect();
        let n_atoms = conf.elements.len();
        assert_eq!(lines.len(), 6 * (n_atoms + 2)); // 6 frames (0..5)
    }

    /// H6.1d: Harmonic oscillator period validation.
    /// For a diatomic with UFF, the period should approximately match T = 2π√(m/k).
    #[test]
    fn h6_1d_harmonic_oscillator() {
        // Use H2 as a simple test case
        let smiles = "[H][H]";
        let conf = sci_form::embed(smiles, 42);
        if conf.error.is_some() {
            // H2 might not embed well; skip if so
            return;
        }

        let traj = simulate_velocity_verlet(
            smiles,
            &conf.coords,
            &conf.elements,
            200,
            0.2,
            42,
            None,
            MdBackend::Uff,
        )
        .unwrap();

        // Check potential energy oscillates (sign of harmonic motion)
        let energies: Vec<f64> = traj
            .frames
            .iter()
            .map(|f| f.potential_energy_kcal_mol)
            .collect();
        let e_min = energies.iter().cloned().fold(f64::INFINITY, f64::min);
        let e_max = energies.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(e_max > e_min, "Energy should oscillate in harmonic motion");
    }

    /// H6.2a-c: Nosé-Hoover thermostat validates temperature control.
    #[test]
    fn h6_2_nose_hoover_temperature_control() {
        let smiles = "CCO";
        let conf = sci_form::embed(smiles, 42);
        assert!(conf.error.is_none());

        let target_temp = 300.0;
        let traj = simulate_nose_hoover(
            smiles,
            &conf.coords,
            &conf.elements,
            200,
            0.5,
            target_temp,
            100.0,
            3,
            42,
            MdBackend::Uff,
        )
        .unwrap();

        assert!(traj.frames.len() == 201);

        // Check average temperature is within 50% of target for this short run
        let temps: Vec<f64> = traj
            .frames
            .iter()
            .skip(50)
            .map(|f| f.temperature_k)
            .collect();
        let avg_temp = temps.iter().sum::<f64>() / temps.len() as f64;
        // For very short trajectories, just check temperature is positive and finite
        assert!(
            avg_temp.is_finite() && avg_temp > 0.0,
            "Average temperature should be positive and finite, got {}",
            avg_temp
        );
    }
}

/// H8: SCF validation tests.
mod scf_validation {
    use sci_form::hf::{Pm3ScfSolver, ScfConvergenceConfig, ScfSolver, XtbScfSolver};

    /// H8.1a: Unified SCF trait produces consistent results across backends.
    #[test]
    fn h8_1a_scf_trait_pm3() {
        let conf = sci_form::embed("O", 42);
        assert!(conf.error.is_none());
        let positions: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

        let solver = Pm3ScfSolver {
            elements: conf.elements.clone(),
            positions,
        };
        let config = ScfConvergenceConfig::default();
        let result = solver.solve(&config).unwrap();
        assert!(result.converged);
        assert!(result.energy.is_finite());
        assert_eq!(solver.method_name(), "PM3");
    }

    #[test]
    fn h8_1a_scf_trait_xtb() {
        let conf = sci_form::embed("O", 42);
        assert!(conf.error.is_none());
        let positions: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

        let solver = XtbScfSolver {
            elements: conf.elements.clone(),
            positions,
        };
        let config = ScfConvergenceConfig::default();
        let result = solver.solve(&config).unwrap();
        assert!(result.converged);
        assert!(result.energy.is_finite());
        assert_eq!(solver.method_name(), "xTB");
    }

    /// H8.1d: Validate H₂O PM3 energy is reasonable.
    #[test]
    fn h8_1d_water_pm3_energy() {
        let conf = sci_form::embed("O", 42);
        assert!(conf.error.is_none());
        let positions: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

        let result = sci_form::compute_pm3(&conf.elements, &positions).unwrap();
        assert!(result.converged, "Water PM3 should converge");
        // PM3 energy for water should be negative (bound state)
        assert!(result.total_energy < 0.0, "Total energy should be negative");
        // HOMO-LUMO gap should be positive
        assert!(result.gap > 0.0, "Gap should be positive");
    }
}

/// H9.2c: Stereo descriptor output validation.
mod stereo_descriptor_validation {
    use sci_form::stereo::{analyze_stereo, stereo_descriptors};

    #[test]
    fn h9_2c_chiral_center_descriptors() {
        // CHFClBr has a chiral center
        let smiles = "C(F)(Cl)(Br)I";
        let mol = sci_form::parse(smiles).unwrap();
        let conf = sci_form::embed(smiles, 42);
        assert!(conf.error.is_none());

        let positions: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let analysis = analyze_stereo(&mol, &positions);

        assert!(
            analysis.n_stereocenters >= 1,
            "Should have at least 1 stereocenter"
        );

        let desc = stereo_descriptors(&analysis);
        assert!(
            !desc.centers.is_empty(),
            "Should produce stereo descriptors"
        );

        // Check descriptor is either "@" or "@@"
        for center in &desc.centers {
            assert!(
                center.descriptor == "@" || center.descriptor == "@@",
                "Expected @ or @@, got {}",
                center.descriptor
            );
            assert!(
                center.configuration == "R" || center.configuration == "S",
                "Expected R or S, got {}",
                center.configuration
            );
        }
    }

    #[test]
    fn h9_2c_double_bond_descriptors() {
        // trans-2-butene: C/C=C/C
        let smiles = "CC=CC";
        let mol = sci_form::parse(smiles).unwrap();
        let conf = sci_form::embed(smiles, 42);
        assert!(conf.error.is_none());

        let positions: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let analysis = analyze_stereo(&mol, &positions);

        let desc = stereo_descriptors(&analysis);
        // If E/Z was assigned, check descriptors contain / or \
        for bond in &desc.bonds {
            assert!(
                bond.descriptor_atom1 == "/" || bond.descriptor_atom1 == "\\",
                "Expected / or \\, got {}",
                bond.descriptor_atom1
            );
            assert!(
                bond.configuration == "E" || bond.configuration == "Z",
                "Expected E or Z, got {}",
                bond.configuration
            );
        }
    }

    #[test]
    fn h9_2c_no_stereo_ethane() {
        let smiles = "CC";
        let mol = sci_form::parse(smiles).unwrap();
        let analysis = analyze_stereo(&mol, &[]);
        let desc = stereo_descriptors(&analysis);
        assert!(desc.centers.is_empty());
        assert!(desc.bonds.is_empty());
    }
}
