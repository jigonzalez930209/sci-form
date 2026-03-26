//! Integration tests for Track D spectroscopy features:
//! UV-Vis (sTDA), IR (numerical Hessian), and NMR (empirical shifts + Karplus).

// ─── Helpers ─────────────────────────────────────────────────────────────────

fn water_molecule() -> (Vec<u8>, Vec<[f64; 3]>) {
    let elements = vec![8, 1, 1];
    let positions = vec![
        [0.0, 0.0, 0.1173],
        [0.0, 0.7572, -0.4692],
        [0.0, -0.7572, -0.4692],
    ];
    (elements, positions)
}

fn embed_smiles(smiles: &str) -> (Vec<u8>, Vec<[f64; 3]>) {
    let conf = sci_form::embed(smiles, 42);
    assert!(conf.error.is_none(), "embed failed: {:?}", conf.error);
    let pos: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    (conf.elements, pos)
}

// ═══════════════════════════════════════════════════════════════════════════════
// UV-Vis sTDA tests
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn test_stda_uvvis_benzene_gaussian() {
    let (elems, pos) = embed_smiles("c1ccccc1");
    let spectrum = sci_form::compute_stda_uvvis(
        &elems,
        &pos,
        0.3,
        1.0,
        8.0,
        500,
        sci_form::reactivity::BroadeningType::Gaussian,
    )
    .expect("sTDA UV-Vis should succeed for benzene");

    assert_eq!(spectrum.energies_ev.len(), 500);
    assert_eq!(spectrum.absorptivity.len(), 500);
    assert!(
        !spectrum.excitations.is_empty(),
        "benzene should have excitations"
    );

    // Benzene π→π* transitions are in the 4–7 eV range
    for exc in &spectrum.excitations {
        assert!(exc.energy_ev > 0.0, "excitation energy must be positive");
        assert!(exc.wavelength_nm > 0.0, "wavelength must be positive");
        assert!(
            exc.oscillator_strength >= 0.0,
            "oscillator strength must be non-negative"
        );
    }

    // Spectrum should have some non-zero absorption
    let max_abs = spectrum
        .absorptivity
        .iter()
        .cloned()
        .fold(0.0_f64, f64::max);
    assert!(max_abs > 0.0, "spectrum should have non-zero absorption");
}

#[test]
fn test_stda_uvvis_lorentzian_broadening() {
    let (elems, pos) = embed_smiles("c1ccccc1");
    let gauss = sci_form::compute_stda_uvvis(
        &elems,
        &pos,
        0.3,
        1.0,
        8.0,
        500,
        sci_form::reactivity::BroadeningType::Gaussian,
    )
    .unwrap();
    let lorentz = sci_form::compute_stda_uvvis(
        &elems,
        &pos,
        0.3,
        1.0,
        8.0,
        500,
        sci_form::reactivity::BroadeningType::Lorentzian,
    )
    .unwrap();

    // Both spectra should have the same grid size and same excitations
    assert_eq!(gauss.energies_ev.len(), lorentz.energies_ev.len());
    assert_eq!(gauss.excitations.len(), lorentz.excitations.len());

    // But different broadening shape → typically different intensity profiles
    // (at least one point should differ unless intensity is all zero)
    let any_differ = gauss
        .absorptivity
        .iter()
        .zip(&lorentz.absorptivity)
        .any(|(g, l)| (g - l).abs() > 1e-10);
    // If there's absorption at all, shapes must differ
    let max_gauss = gauss.absorptivity.iter().cloned().fold(0.0_f64, f64::max);
    if max_gauss > 1e-6 {
        assert!(any_differ, "Gaussian and Lorentzian profiles should differ");
    }
}

#[test]
fn test_stda_uvvis_ethanol() {
    let (elems, pos) = embed_smiles("CCO");
    let spectrum = sci_form::compute_stda_uvvis(
        &elems,
        &pos,
        0.3,
        1.0,
        10.0,
        300,
        sci_form::reactivity::BroadeningType::Gaussian,
    )
    .expect("sTDA should succeed for ethanol");

    assert_eq!(spectrum.energies_ev.len(), 300);
    assert_eq!(spectrum.absorptivity.len(), 300);
}

#[test]
fn test_stda_uvvis_excitation_properties() {
    let (elems, pos) = embed_smiles("c1ccccc1");
    let spectrum = sci_form::compute_stda_uvvis(
        &elems,
        &pos,
        0.3,
        1.0,
        8.0,
        500,
        sci_form::reactivity::BroadeningType::Gaussian,
    )
    .unwrap();

    for exc in &spectrum.excitations {
        // Energy-wavelength consistency: E(eV) = 1239.84 / λ(nm)
        let expected_nm = 1239.84 / exc.energy_ev;
        assert!(
            (exc.wavelength_nm - expected_nm).abs() < 1.0,
            "wavelength {:.1} nm inconsistent with energy {:.3} eV (expected {:.1} nm)",
            exc.wavelength_nm,
            exc.energy_ev,
            expected_nm
        );
    }
}

#[test]
fn test_stda_uvvis_water_is_deep_uv_only() {
    let (elems, pos) = water_molecule();
    let spectrum = sci_form::compute_stda_uvvis(
        &elems,
        &pos,
        0.3,
        1.0,
        10.0,
        500,
        sci_form::reactivity::BroadeningType::Gaussian,
    )
    .expect("sTDA UV-Vis should succeed for water");

    let strongest_long_wavelength_excitation = spectrum
        .excitations
        .iter()
        .filter(|excitation| excitation.wavelength_nm >= 220.0)
        .map(|excitation| excitation.oscillator_strength)
        .fold(0.0_f64, f64::max);

    assert!(
        strongest_long_wavelength_excitation < 1e-3,
        "water should be effectively transparent above ~220 nm in this low-cost model, got f_max={}",
        strongest_long_wavelength_excitation
    );

    let strongest_excitation = spectrum
        .excitations
        .iter()
        .map(|excitation| excitation.wavelength_nm)
        .fold(0.0_f64, f64::max);
    assert!(
        strongest_excitation < 220.0,
        "water UV-Vis bands should remain in the deep UV, got strongest wavelength {} nm",
        strongest_excitation
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// IR Spectroscopy tests
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn test_ir_vibrational_analysis_water_eht() {
    let (elems, pos) = water_molecule();
    let analysis = sci_form::compute_vibrational_analysis(&elems, &pos, "eht", None)
        .expect("vibrational analysis should succeed for water");

    assert_eq!(analysis.n_atoms, 3);
    assert_eq!(analysis.method, "EHT");

    // Water: 3N – 6 = 3 real vibrational modes (non-linear)
    assert!(
        analysis.n_real_modes > 0,
        "water should have real vibrational modes"
    );

    // ZPVE should be positive for real modes
    assert!(analysis.zpve_ev >= 0.0, "ZPVE must be non-negative");

    // Check that some modes have significant frequency magnitude
    let significant_modes: Vec<_> = analysis.modes.iter().filter(|m| m.is_real).collect();
    assert!(
        !significant_modes.is_empty(),
        "should have at least one significant mode"
    );

    // Positive-frequency modes are real vibrations;
    // negative frequencies (imaginary) can occur when the geometry is
    // not a true minimum on the semi-empirical PES.
    let positive_modes: Vec<_> = significant_modes
        .iter()
        .filter(|m| m.frequency_cm1 > 0.0)
        .collect();
    assert!(
        !positive_modes.is_empty(),
        "should have at least one positive-frequency mode"
    );
}

#[test]
fn test_ir_vibrational_analysis_ethanol() {
    let (elems, pos) = embed_smiles("CCO");
    let analysis = sci_form::compute_vibrational_analysis(&elems, &pos, "eht", None)
        .expect("vibrational analysis should succeed for ethanol");

    // Ethanol has 9 atoms → 3*9 - 6 = 21 vibrational modes
    assert_eq!(analysis.n_atoms, 9);
    assert!(
        !analysis.modes.is_empty(),
        "ethanol should produce vibrational modes"
    );

    // At least some modes should have IR intensity > 0
    let active_modes: Vec<_> = analysis
        .modes
        .iter()
        .filter(|m| m.is_real && m.ir_intensity > 0.001)
        .collect();
    assert!(
        !active_modes.is_empty(),
        "ethanol should have IR-active modes"
    );
}

#[test]
fn test_ir_spectrum_generation_water() {
    let (elems, pos) = water_molecule();
    let analysis = sci_form::compute_vibrational_analysis(&elems, &pos, "eht", None).unwrap();

    let spectrum = sci_form::compute_ir_spectrum(&analysis, 15.0, 0.0, 4000.0, 1000);

    assert_eq!(spectrum.wavenumbers.len(), 1000);
    assert_eq!(spectrum.intensities.len(), 1000);
    assert_eq!(spectrum.gamma, 15.0);

    // Wavenumber grid first and last values should be finite
    assert!(spectrum.wavenumbers[0].is_finite());
    assert!(spectrum.wavenumbers[999].is_finite());
    // Grid spacing should be approximately (wn_max - wn_min) / n_points
    let step = (spectrum.wavenumbers[1] - spectrum.wavenumbers[0]).abs();
    assert!(step > 0.0 && step < 100.0, "grid step should be reasonable");

    // Peaks should have valid data
    for peak in &spectrum.peaks {
        assert!(peak.frequency_cm1.is_finite());
        assert!(peak.ir_intensity >= 0.0);
    }
}

#[test]
fn test_ir_spectrum_different_resolutions() {
    let (elems, pos) = water_molecule();
    let analysis = sci_form::compute_vibrational_analysis(&elems, &pos, "eht", None).unwrap();

    let sp_100 = sci_form::compute_ir_spectrum(&analysis, 15.0, 0.0, 4000.0, 100);
    let sp_500 = sci_form::compute_ir_spectrum(&analysis, 15.0, 0.0, 4000.0, 500);

    assert_eq!(sp_100.wavenumbers.len(), 100);
    assert_eq!(sp_500.wavenumbers.len(), 500);

    // Same peaks regardless of resolution
    assert_eq!(sp_100.peaks.len(), sp_500.peaks.len());
}

#[test]
fn test_ir_method_selection() {
    let (elems, pos) = water_molecule();

    // EHT should work
    let res_eht = sci_form::compute_vibrational_analysis(&elems, &pos, "eht", None);
    assert!(res_eht.is_ok());

    // PM3 should work (water is small enough)
    let res_pm3 = sci_form::compute_vibrational_analysis(&elems, &pos, "pm3", None);
    assert!(res_pm3.is_ok());

    // xTB should work
    let res_xtb = sci_form::compute_vibrational_analysis(&elems, &pos, "xtb", None);
    assert!(res_xtb.is_ok());

    // Unknown method should fail
    let res_bad = sci_form::compute_vibrational_analysis(&elems, &pos, "garbage", None);
    assert!(res_bad.is_err());
}

#[test]
fn test_ir_custom_step_size() {
    let (elems, pos) = water_molecule();

    let default_step = sci_form::compute_vibrational_analysis(&elems, &pos, "eht", None).unwrap();
    let small_step =
        sci_form::compute_vibrational_analysis(&elems, &pos, "eht", Some(0.001)).unwrap();

    // Both should produce modes
    assert!(default_step.n_real_modes > 0);
    assert!(small_step.n_real_modes > 0);
}

#[test]
fn test_ir_vibrational_analysis_uff_water() {
    let (elems, pos) = water_molecule();
    let analysis = sci_form::compute_vibrational_analysis_uff("O", &elems, &pos, None)
        .expect("UFF vibrational analysis should succeed for water");

    assert_eq!(analysis.method, "UFF");
    assert_eq!(analysis.elements, elems);
    assert_eq!(
        analysis.n_real_modes, 3,
        "water should expose 3 vibrational modes after rigid-body projection"
    );

    let positive_modes: Vec<f64> = analysis
        .modes
        .iter()
        .filter(|mode| mode.is_real && mode.frequency_cm1 > 0.0)
        .map(|mode| mode.frequency_cm1)
        .collect();

    assert_eq!(
        positive_modes.len(),
        3,
        "water should retain 3 positive vibrational bands with UFF"
    );
    assert!(
        positive_modes
            .iter()
            .any(|&frequency| (1200.0..1900.0).contains(&frequency)),
        "water should show a qualitative bending mode in the mid-IR, got {:?}",
        positive_modes
    );
    assert!(
        positive_modes
            .iter()
            .filter(|&&frequency| (3000.0..4200.0).contains(&frequency))
            .count()
            >= 2,
        "water should show two O-H stretching modes near 3600 cm^-1, got {:?}",
        positive_modes
    );
    assert!(
        analysis
            .notes
            .iter()
            .any(|note| note.contains("Gasteiger") || note.contains("neutral-charge")),
        "UFF analysis should explain the IR-intensity approximation"
    );
    assert!(
        analysis
            .notes
            .iter()
            .any(|note| note.contains("projected out before diagonalization")),
        "UFF analysis should mention rigid-body projection"
    );
}

#[test]
fn test_ir_peak_assignment_identifies_water_bend() {
    let result = sci_form::ir::assign_peaks(&[1640.0, 3600.0], &[12.0, 25.0], &[8, 1, 1], None);

    assert!(
        result
            .functional_groups
            .iter()
            .any(|group| group == "H-O-H bend"),
        "water bending band should be identified explicitly"
    );
    assert!(
        result
            .assignments
            .iter()
            .flat_map(|assignment| assignment.assignments.iter())
            .any(|group| group.group.contains("O-H")),
        "water stretching region should still match an O-H assignment"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// NMR Spectroscopy tests
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn test_nmr_shifts_ethanol() {
    let result = sci_form::predict_nmr_shifts("CCO").expect("NMR shifts should work for ethanol");

    // Ethanol should have ¹H shifts
    assert!(
        !result.h_shifts.is_empty(),
        "ethanol has protons → non-empty ¹H shifts"
    );

    // Ethanol has 2 carbons
    assert!(
        !result.c_shifts.is_empty(),
        "ethanol has carbons → non-empty ¹³C shifts"
    );

    // Check that ¹H shifts are in a reasonable range (0–15 ppm)
    for shift in &result.h_shifts {
        assert!(
            shift.shift_ppm >= -1.0 && shift.shift_ppm <= 15.0,
            "unreasonable ¹H shift: {} ppm",
            shift.shift_ppm
        );
        assert_eq!(shift.element, 1, "¹H shifts should be for hydrogen");
    }

    // Check that ¹³C shifts are in a reasonable range (0–220 ppm)
    for shift in &result.c_shifts {
        assert!(
            shift.shift_ppm >= -10.0 && shift.shift_ppm <= 230.0,
            "unreasonable ¹³C shift: {} ppm",
            shift.shift_ppm
        );
        assert_eq!(shift.element, 6, "¹³C shifts should be for carbon");
    }
}

#[test]
fn test_nmr_shifts_benzene() {
    let result =
        sci_form::predict_nmr_shifts("c1ccccc1").expect("NMR shifts should work for benzene");

    // All ¹H in benzene should be around 7.27 ppm (aromatic)
    for shift in &result.h_shifts {
        assert!(
            shift.shift_ppm > 5.0 && shift.shift_ppm < 10.0,
            "benzene ¹H should be aromatic (~7 ppm), got {}",
            shift.shift_ppm
        );
    }

    // All ¹³C in benzene should be around 128 ppm (aromatic)
    for shift in &result.c_shifts {
        assert!(
            shift.shift_ppm > 100.0 && shift.shift_ppm < 160.0,
            "benzene ¹³C should be aromatic (~128 ppm), got {}",
            shift.shift_ppm
        );
    }
}

#[test]
fn test_nmr_shifts_acetic_acid() {
    let result =
        sci_form::predict_nmr_shifts("CC(=O)O").expect("NMR shifts should work for acetic acid");

    assert!(!result.h_shifts.is_empty());
    assert!(!result.c_shifts.is_empty());

    // Should have notes about the prediction method
    // (notes may or may not be populated depending on implementation)
}

#[test]
fn test_nmr_couplings_ethanol() {
    let (_, pos) = embed_smiles("CCO");
    let couplings =
        sci_form::predict_nmr_couplings("CCO", &pos).expect("J-couplings should work for ethanol");

    // Ethanol has H-C-C-H vicinal couplings (³J) and possibly geminal (²J)
    for coupling in &couplings {
        assert!(
            coupling.j_hz.abs() < 20.0,
            "unreasonable J-coupling: {} Hz",
            coupling.j_hz
        );
        assert!(
            coupling.n_bonds == 2 || coupling.n_bonds == 3,
            "expected geminal or vicinal coupling, got {} bonds",
            coupling.n_bonds
        );
    }
}

#[test]
fn test_nmr_couplings_without_3d() {
    // Should still work with empty positions (topological estimate)
    let couplings = sci_form::predict_nmr_couplings("CCO", &[])
        .expect("J-couplings should work without 3D coords");

    // Still should find vicinal H-H couplings topologically
    for coupling in &couplings {
        assert!(coupling.j_hz.abs() < 20.0);
    }
}

#[test]
fn test_nmr_spectrum_h1_ethanol() {
    let spectrum = sci_form::compute_nmr_spectrum("CCO", "1H", 0.02, 0.0, 12.0, 1000)
        .expect("¹H NMR spectrum should succeed for ethanol");

    assert_eq!(spectrum.ppm_axis.len(), 1000);
    assert_eq!(spectrum.intensities.len(), 1000);

    // NMR convention: ppm axis runs high to low
    assert!(
        spectrum.ppm_axis[0] > spectrum.ppm_axis[999],
        "ppm axis should decrease (NMR convention)"
    );

    // Should have peaks
    assert!(
        !spectrum.peaks.is_empty(),
        "ethanol ¹H spectrum should have peaks"
    );

    // Peaks should be in range
    for peak in &spectrum.peaks {
        assert!(
            peak.shift_ppm >= 0.0 && peak.shift_ppm <= 12.0,
            "peak at {} ppm outside spectral window",
            peak.shift_ppm
        );
    }
}

#[test]
fn test_nmr_spectrum_with_coords_uses_public_geometry_path() {
    let (_, pos) = embed_smiles("CCO");

    let couplings_2d = sci_form::predict_nmr_couplings("CCO", &[]).unwrap();

    let couplings_3d = sci_form::predict_nmr_couplings("CCO", &pos).unwrap();

    assert!(
        couplings_2d
            .iter()
            .zip(couplings_3d.iter())
            .any(|(topo, geom)| (topo.j_hz - geom.j_hz).abs() > 0.05),
        "embedded 3D geometry should change at least one vicinal coupling"
    );

    let spectrum =
        sci_form::compute_nmr_spectrum_with_coords("CCO", &pos, "1H", 0.02, 0.0, 12.0, 1000)
            .expect("coordinate-aware ¹H NMR spectrum should succeed for ethanol");
    assert_eq!(spectrum.ppm_axis.len(), 1000);
    assert!(!spectrum.peaks.is_empty());
}

#[test]
fn test_nmr_spectrum_c13_benzene() {
    let spectrum = sci_form::compute_nmr_spectrum("c1ccccc1", "13C", 0.5, 0.0, 220.0, 2000)
        .expect("¹³C NMR spectrum should succeed for benzene");

    assert_eq!(spectrum.ppm_axis.len(), 2000);
    assert_eq!(spectrum.intensities.len(), 2000);
    assert!(!spectrum.peaks.is_empty());

    // All carbons in benzene are equivalent → one peak region around 128 ppm
    let peak_ppms: Vec<f64> = spectrum.peaks.iter().map(|p| p.shift_ppm).collect();
    assert!(
        peak_ppms.iter().all(|&p| p > 100.0 && p < 160.0),
        "benzene ¹³C peaks should be in aromatic range: {:?}",
        peak_ppms
    );
}

#[test]
fn test_nmr_spectrum_nucleus_aliases() {
    // All these should work for ¹H
    for alias in &["1H", "H1", "h1", "1h", "proton"] {
        let res = sci_form::compute_nmr_spectrum("C", alias, 0.02, 0.0, 12.0, 100);
        assert!(res.is_ok(), "nucleus alias '{}' should work", alias);
    }
    // All these should work for ¹³C
    for alias in &["13C", "C13", "c13", "13c", "carbon"] {
        let res = sci_form::compute_nmr_spectrum("C", alias, 0.5, 0.0, 220.0, 100);
        assert!(res.is_ok(), "nucleus alias '{}' should work", alias);
    }
    // All these should also work for ¹⁹F
    for alias in &["19F", "F19", "f19", "19f", "fluorine"] {
        let res = sci_form::compute_nmr_spectrum("FC", alias, 0.5, -250.0, 0.0, 100);
        assert!(res.is_ok(), "nucleus alias '{}' should work", alias);
    }
    // Expanded nuclei should also work
    for (smiles, alias, ppm_min, ppm_max) in [
        ("[H]", "2H", -2.0, 14.0),
        ("[Cl]", "35Cl", -400.0, 1600.0),
        ("[Br]", "79Br", -600.0, 2200.0),
        ("[Pt]", "195Pt", -3000.0, 3000.0),
    ] {
        let res = sci_form::compute_nmr_spectrum(smiles, alias, 0.0, ppm_min, ppm_max, 100);
        assert!(res.is_ok(), "expanded nucleus alias '{}' should work", alias);
    }
    // Unknown nucleus should fail
    let res = sci_form::compute_nmr_spectrum("C", "999X", 0.02, 0.0, 12.0, 100);
    assert!(res.is_err());
}

#[test]
fn test_nmr_shifts_for_expanded_nuclei() {
    let cases = [
        ("1H", "[H]"),
        ("2H", "[H]"),
        ("3H", "[H]"),
        ("3He", "[He]"),
        ("6Li", "[Li]"),
        ("7Li", "[Li]"),
        ("9Be", "[Be]"),
        ("10B", "[B]"),
        ("11B", "[B]"),
        ("13C", "[C]"),
        ("14N", "[N]"),
        ("15N", "[N]"),
        ("17O", "[O]"),
        ("19F", "[F]"),
        ("23Na", "[Na]"),
        ("25Mg", "[Mg]"),
        ("27Al", "[Al]"),
        ("29Si", "[Si]"),
        ("31P", "[P]"),
        ("33S", "[S]"),
        ("35Cl", "[Cl]"),
        ("37Cl", "[Cl]"),
        ("39K", "[K]"),
        ("40K", "[K]"),
        ("41K", "[K]"),
        ("43Ca", "[Ca]"),
        ("45Sc", "[Sc]"),
        ("47Ti", "[Ti]"),
        ("49Ti", "[Ti]"),
        ("50V", "[V]"),
        ("51V", "[V]"),
        ("53Cr", "[Cr]"),
        ("55Mn", "[Mn]"),
        ("57Fe", "[Fe]"),
        ("59Co", "[Co]"),
        ("61Ni", "[Ni]"),
        ("63Cu", "[Cu]"),
        ("65Cu", "[Cu]"),
        ("67Zn", "[Zn]"),
        ("69Ga", "[Ga]"),
        ("71Ga", "[Ga]"),
        ("73Ge", "[Ge]"),
        ("75As", "[As]"),
        ("77Se", "[Se]"),
        ("79Br", "[Br]"),
        ("81Br", "[Br]"),
        ("85Rb", "[Rb]"),
        ("87Rb", "[Rb]"),
        ("87Sr", "[Sr]"),
        ("91Zr", "[Zr]"),
        ("93Nb", "[Nb]"),
        ("95Mo", "[Mo]"),
        ("97Mo", "[Mo]"),
        ("99Ru", "[Ru]"),
        ("101Ru", "[Ru]"),
        ("103Rh", "[Rh]"),
        ("105Pd", "[Pd]"),
        ("107Ag", "[Ag]"),
        ("109Ag", "[Ag]"),
        ("111Cd", "[Cd]"),
        ("113Cd", "[Cd]"),
        ("113In", "[In]"),
        ("115In", "[In]"),
        ("115Sn", "[Sn]"),
        ("117Sn", "[Sn]"),
        ("119Sn", "[Sn]"),
        ("121Sb", "[Sb]"),
        ("123Sb", "[Sb]"),
        ("123Te", "[Te]"),
        ("125Te", "[Te]"),
        ("127I", "[I]"),
        ("129Xe", "[Xe]"),
        ("131Xe", "[Xe]"),
        ("133Cs", "[Cs]"),
        ("135Ba", "[Ba]"),
        ("137Ba", "[Ba]"),
        ("183W", "[W]"),
        ("195Pt", "[Pt]"),
        ("197Au", "[Au]"),
        ("199Hg", "[Hg]"),
        ("201Hg", "[Hg]"),
        ("203Tl", "[Tl]"),
        ("205Tl", "[Tl]"),
        ("207Pb", "[Pb]"),
        ("209Bi", "[Bi]"),
    ];

    for (nucleus, smiles) in cases {
        let result = sci_form::predict_nmr_shifts_for_nucleus(smiles, nucleus)
            .unwrap_or_else(|err| panic!("{} on {} failed: {}", nucleus, smiles, err));
        assert!(
            !result.is_empty(),
            "{} on {} should yield at least one shift",
            nucleus,
            smiles
        );
        assert!(
            result.iter().all(|shift| shift.shift_ppm.is_finite() && shift.confidence.is_finite()),
            "{} should produce finite relative shifts",
            nucleus
        );
    }
}

#[test]
fn test_public_giao_nmr_water_h1() {
    let (elements, positions) = water_molecule();
    let result = sci_form::compute_giao_nmr(&elements, &positions, "1H")
        .expect("public GIAO NMR should succeed for water proton shifts");

    assert_eq!(result.nucleus, "1H");
    assert_eq!(result.target_atomic_number, 1);
    assert_eq!(result.n_target_atoms, 2);
    assert_eq!(result.shieldings.len(), 2);
    assert_eq!(result.chemical_shifts.len(), 2);
    assert!(result.scf_iterations > 0);
    assert!(result
        .chemical_shifts
        .iter()
        .all(|shift| shift.is_finite()));
    assert!(result
        .notes
        .iter()
        .any(|note| note.contains("RHF/STO-3G")));
}

#[test]
fn test_public_giao_nmr_rejects_fallback_basis_by_default() {
    let error = sci_form::compute_giao_nmr(&[78], &[[0.0, 0.0, 0.0]], "195Pt")
        .expect_err("unsupported heavy elements should fail before SCF by default");

    assert!(error.contains("fallback-only elements"), "{error}");
}

#[test]
fn test_public_giao_nmr_reports_inadequate_forced_fallback_basis() {
    let config = sci_form::GiaoNmrConfig {
        allow_basis_fallback: true,
        ..Default::default()
    };
    let error = sci_form::compute_giao_nmr_configured(
        &[78],
        &[[0.0, 0.0, 0.0]],
        "195Pt",
        &config,
    )
    .expect_err("forcing fallback should still fail when the basis cannot host the electrons");

    assert!(error.contains("basis functions"), "{error}");
}

#[test]
fn test_hose_codes_ethanol() {
    let codes = sci_form::compute_hose_codes("CCO", 2).expect("HOSE codes should work");

    // Ethanol has 9 atoms → 9 HOSE codes
    assert!(!codes.is_empty());

    for code in &codes {
        assert!(!code.full_code.is_empty(), "HOSE code should be non-empty");
        assert!(!code.spheres.is_empty(), "should have at least one sphere");
    }
}

#[test]
fn test_hose_codes_benzene() {
    let codes = sci_form::compute_hose_codes("c1ccccc1", 3).expect("HOSE codes should work");

    assert!(!codes.is_empty());

    // All carbons in benzene should produce similar HOSE environments
    let carbon_codes: Vec<_> = codes.iter().filter(|c| c.element == 6).collect();
    assert!(
        carbon_codes.len() == 6,
        "benzene has 6 carbons, got {} carbon HOSE codes",
        carbon_codes.len()
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// End-to-end pipeline tests
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn test_full_spectroscopy_pipeline_ethanol() {
    // Step 1: Generate 3D conformer
    let conf = sci_form::embed("CCO", 42);
    assert!(conf.error.is_none());
    let pos: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

    // Step 2: UV-Vis (requires 3D)
    let uvvis = sci_form::compute_stda_uvvis(
        &conf.elements,
        &pos,
        0.3,
        1.0,
        10.0,
        300,
        sci_form::reactivity::BroadeningType::Gaussian,
    )
    .expect("UV-Vis should succeed");
    assert_eq!(uvvis.energies_ev.len(), 300);

    // Step 3: IR (requires 3D)
    let ir_analysis = sci_form::compute_vibrational_analysis(&conf.elements, &pos, "eht", None)
        .expect("IR analysis should succeed");
    let ir_spectrum = sci_form::compute_ir_spectrum(&ir_analysis, 15.0, 400.0, 4000.0, 500);
    assert_eq!(ir_spectrum.wavenumbers.len(), 500);

    // Step 4: NMR (topology only — no 3D needed for shifts)
    let nmr_shifts = sci_form::predict_nmr_shifts("CCO").expect("NMR shifts should succeed");
    assert!(!nmr_shifts.h_shifts.is_empty());

    // Step 5: NMR J-couplings (benefits from 3D for Karplus)
    let _nmr_couplings =
        sci_form::predict_nmr_couplings("CCO", &pos).expect("J-couplings should succeed");

    // Step 6: NMR spectrum
    let nmr_h = sci_form::compute_nmr_spectrum("CCO", "1H", 0.02, 0.0, 12.0, 1000)
        .expect("¹H NMR spectrum should succeed");
    assert_eq!(nmr_h.ppm_axis.len(), 1000);

    let nmr_c = sci_form::compute_nmr_spectrum("CCO", "13C", 0.5, 0.0, 220.0, 1000)
        .expect("¹³C NMR spectrum should succeed");
    assert_eq!(nmr_c.ppm_axis.len(), 1000);

    // All spectra should have valid data
    assert!(uvvis.absorptivity.iter().all(|v| v.is_finite()));
    assert!(ir_spectrum.intensities.iter().all(|v| v.is_finite()));
    assert!(nmr_h.intensities.iter().all(|v| v.is_finite()));
    assert!(nmr_c.intensities.iter().all(|v| v.is_finite()));
}

#[test]
fn test_full_spectroscopy_pipeline_benzene() {
    let (elems, pos) = embed_smiles("c1ccccc1");

    // UV-Vis
    let uvvis = sci_form::compute_stda_uvvis(
        &elems,
        &pos,
        0.3,
        1.0,
        8.0,
        500,
        sci_form::reactivity::BroadeningType::Lorentzian,
    )
    .expect("UV-Vis should succeed for benzene");
    assert!(!uvvis.excitations.is_empty());

    // IR
    let ir = sci_form::compute_vibrational_analysis(&elems, &pos, "eht", None)
        .expect("IR should succeed for benzene");
    assert!(ir.n_real_modes > 0);

    // NMR — ¹H
    let nmr_h = sci_form::compute_nmr_spectrum("c1ccccc1", "1H", 0.02, 5.0, 10.0, 500)
        .expect("¹H NMR should succeed for benzene");
    assert!(!nmr_h.peaks.is_empty());

    // NMR — ¹³C
    let nmr_c = sci_form::compute_nmr_spectrum("c1ccccc1", "13C", 0.5, 100.0, 160.0, 500)
        .expect("¹³C NMR should succeed for benzene");
    assert!(!nmr_c.peaks.is_empty());
}

#[test]
fn test_nmr_serialization() {
    // Verify that all NMR types are serializable (they derive Serialize)
    let result = sci_form::predict_nmr_shifts("CCO").unwrap();
    let json = serde_json::to_string(&result).expect("NmrShiftResult should serialize");
    assert!(json.contains("shift_ppm"));

    let spectrum = sci_form::compute_nmr_spectrum("CCO", "1H", 0.02, 0.0, 12.0, 100).unwrap();
    let json = serde_json::to_string(&spectrum).expect("NmrSpectrum should serialize");
    assert!(json.contains("ppm_axis"));
}

#[test]
fn test_ir_serialization() {
    let (elems, pos) = water_molecule();
    let analysis = sci_form::compute_vibrational_analysis(&elems, &pos, "eht", None).unwrap();
    let json = serde_json::to_string(&analysis).expect("VibrationalAnalysis should serialize");
    assert!(json.contains("modes"));

    let spectrum = sci_form::compute_ir_spectrum(&analysis, 15.0, 0.0, 4000.0, 100);
    let json = serde_json::to_string(&spectrum).expect("IrSpectrum should serialize");
    assert!(json.contains("wavenumbers"));
}

#[test]
fn test_uvvis_serialization() {
    let (elems, pos) = embed_smiles("c1ccccc1");
    let spectrum = sci_form::compute_stda_uvvis(
        &elems,
        &pos,
        0.3,
        1.0,
        8.0,
        50,
        sci_form::reactivity::BroadeningType::Gaussian,
    )
    .unwrap();
    let json = serde_json::to_string(&spectrum).expect("StdaUvVisSpectrum should serialize");
    assert!(json.contains("excitations"));
    assert!(json.contains("absorptivity"));
}
