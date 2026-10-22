//! Regression tests for the v1.0 crash blockers (see
//! `reports/v1.0-blockers.md`, section "2. A. Robustez y numérico").
//!
//! Every test in this file reproduces one of the blockers that used to
//! panic, hang, or silently return a wrong/empty result.  They assert the
//! *new* behaviour: an actionable `Err` instead of a panic.

use sci_form::charges::gasteiger::{gasteiger_marsili_charges_configured, GasteigerConfig};
use sci_form::graph::Molecule;

// ─── Helpers ─────────────────────────────────────────────────────────────────

/// Water geometry in Ångström.
fn water() -> (Vec<u8>, Vec<[f64; 3]>) {
    (
        vec![8, 1, 1],
        vec![
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469],
        ],
    )
}

/// Same geometry with one `NaN` coordinate.
fn water_with_nan() -> (Vec<u8>, Vec<[f64; 3]>) {
    let (elements, mut positions) = water();
    positions[1][0] = f64::NAN;
    (elements, positions)
}

/// Same geometry with one `Inf` coordinate.
fn water_with_inf() -> (Vec<u8>, Vec<[f64; 3]>) {
    let (elements, mut positions) = water();
    positions[2][2] = f64::INFINITY;
    (elements, positions)
}

// ─── A-C1 · `Molecule::from_smiles("[")` must not panic ──────────────────────

#[test]
fn a_c1_from_smiles_empty_bracket_atom_returns_error() {
    let err = match Molecule::from_smiles("[") {
        Ok(_) => panic!("an empty bracket atom must be a parse error, not a success"),
        Err(e) => e,
    };
    assert!(
        err.contains("empty bracket atom"),
        "unexpected error message: {}",
        err
    );
}

#[test]
fn a_c1_from_smiles_isotope_only_bracket_returns_error() {
    let err = match Molecule::from_smiles("[13]") {
        Ok(_) => panic!("a bracket without an element symbol must be a parse error"),
        Err(e) => e,
    };
    assert!(
        err.contains("empty bracket atom"),
        "unexpected error message: {}",
        err
    );
}

#[test]
fn a_c1_public_wrappers_reject_the_empty_bracket_atom() {
    // `parse`, `compute_charges` and friends all funnel into the SMILES
    // parser: they must report the error instead of panicking.
    assert!(sci_form::parse("[").is_err());
    assert!(sci_form::compute_charges("[").is_err());
    assert!(sci_form::compute_ecfp("[", 2, 1024).is_err());
    assert!(sci_form::compute_hose_codes("[", 3).is_err());
}

// ─── A-C2 · a NaN/Inf coordinate must be rejected, never panic ──────────────

#[test]
fn a_c2_eht_rejects_nan_coordinates() {
    let (elements, positions) = water_with_nan();
    let err = sci_form::compute_population(&elements, &positions)
        .expect_err("a NaN coordinate must be rejected");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

#[test]
fn a_c2_eht_rejects_inf_coordinates() {
    let (elements, positions) = water_with_inf();
    let err = sci_form::compute_population(&elements, &positions)
        .expect_err("an Inf coordinate must be rejected");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

#[test]
fn a_c2_eht_rejects_length_mismatch() {
    let (elements, positions) = water();
    let err = sci_form::compute_population(&elements, &positions[..2])
        .expect_err("a length mismatch must be rejected");
    assert!(
        err.contains("positions has 2 entries"),
        "unexpected error: {}",
        err
    );
}

#[test]
fn a_c2_uhf_rejects_nan_coordinates() {
    let (elements, positions) = water_with_nan();
    let err = sci_form::compute_uhf(&elements, &positions, 0, 1)
        .expect_err("a NaN coordinate must be rejected by UHF");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

#[test]
fn a_c2_rohf_rejects_nan_coordinates() {
    let (elements, positions) = water_with_nan();
    let err = sci_form::compute_rohf(&elements, &positions, 0, 1)
        .expect_err("a NaN coordinate must be rejected by ROHF");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

#[test]
fn a_c2_hf3c_rejects_nan_coordinates() {
    let (elements, positions) = water_with_nan();
    let err = sci_form::compute_hf3c(&elements, &positions, &sci_form::hf::HfConfig::default())
        .expect_err("a NaN coordinate must be rejected by HF-3c");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

#[test]
fn a_c2_pm3_rejects_nan_coordinates() {
    let (elements, positions) = water_with_nan();
    let err = sci_form::compute_pm3(&elements, &positions)
        .expect_err("a NaN coordinate must be rejected by PM3");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

#[test]
fn a_c2_xtb_rejects_nan_coordinates() {
    let (elements, positions) = water_with_nan();
    let err = sci_form::compute_xtb(&elements, &positions)
        .expect_err("a NaN coordinate must be rejected by xTB");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

#[test]
fn a_c2_ani_rejects_nan_coordinates() {
    let (elements, positions) = water_with_nan();
    let err = sci_form::compute_ani(&elements, &positions)
        .expect_err("a NaN coordinate must be rejected by ANI");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

#[test]
fn a_c2_ani_rejects_inf_coordinates() {
    let (elements, positions) = water_with_inf();
    let err = sci_form::compute_ani(&elements, &positions)
        .expect_err("an Inf coordinate must be rejected by ANI");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

#[test]
fn a_c2_orbital_mesh_rejects_nan_coordinates() {
    let (elements, positions) = water_with_nan();
    let err = sci_form::compute_orbital_mesh(&elements, &positions, "eht", 0, 0.4, 3.0, 0.02)
        .expect_err("a NaN coordinate must be rejected");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

#[test]
fn a_c2_dos_rejects_nan_coordinates() {
    let (elements, positions) = water_with_nan();
    let err = sci_form::compute_dos(&elements, &positions, 0.3, -30.0, 5.0, 100)
        .expect_err("a NaN coordinate must be rejected");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

#[test]
fn a_c2_stda_uvvis_rejects_nan_coordinates() {
    let (elements, positions) = water_with_nan();
    let err = sci_form::compute_stda_uvvis(
        &elements,
        &positions,
        0.3,
        1.0,
        20.0,
        50,
        sci_form::reactivity::BroadeningType::Gaussian,
    )
    .expect_err("a NaN coordinate must be rejected");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

#[test]
fn a_c2_fluorescence_rejects_nan_coordinates() {
    let (elements, positions) = water_with_nan();
    let err = sci_form::compute_fluorescence(
        &elements,
        &positions,
        sci_form::spectroscopy::FluorescenceConfig::default(),
    )
    .expect_err("a NaN coordinate must be rejected");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

#[test]
fn a_c2_vibrational_analysis_rejects_nan_coordinates() {
    let (elements, positions) = water_with_nan();
    let err = sci_form::compute_vibrational_analysis(&elements, &positions, "eht", None)
        .expect_err("a NaN coordinate must be rejected");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

#[test]
fn a_c2_band_structure_rejects_nan_coordinates() {
    let (elements, positions) = water_with_nan();
    let lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
    let config = sci_form::eht::band_structure::BandStructureConfig::default();
    let err = sci_form::eht::band_structure::compute_band_structure(
        &elements, &positions, &lattice, &config, 8,
    )
    .expect_err("a NaN coordinate must be rejected");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

#[test]
fn a_c2_esp_rejects_nan_coordinates() {
    let (elements, positions) = water_with_nan();
    let err = sci_form::compute_esp(&elements, &positions, 0.5, 3.0)
        .expect_err("a NaN coordinate must be rejected");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

#[test]
fn a_c2_flat_coordinate_helper_reports_the_documented_message() {
    // The flat-coordinate helpers use the documented message style.
    let err = sci_form::compute_sasa(&[1, 1, 1], &[0.0; 8], None)
        .expect_err("a length mismatch must be rejected");
    assert_eq!(
        err,
        "coords has 8 values but 3 atoms require 9 (3 atoms x 3 coordinates)"
    );
}

// ─── A-C4 / A-H7 · UHF charge and multiplicity validation ───────────────────
#[test]
fn a_c2_scf_reports_a_nan_geometry_instead_of_panicking() {
    use sci_form::scf::scf_loop::{run_scf, try_run_scf, ScfConfig};
    use sci_form::scf::types::MolecularSystem;

    let system = MolecularSystem::from_angstrom(
        &[1, 1],
        &[[0.0, 0.0, 0.0], [f64::NAN, 0.0, 0.0]],
        0,
        1,
    );

    let err = try_run_scf(&system, &ScfConfig::default())
        .err()
        .expect("a NaN coordinate must be rejected by the SCF driver");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);

    // The infallible driver must not panic either: it reports a
    // non-converged result carrying the reason.
    let result = run_scf(&system, &ScfConfig::default());
    assert!(!result.converged, "an invalid system must not converge");
    let error = result.error.expect("the reason must be reported");
    assert!(error.contains("non-finite"), "unexpected error: {}", error);
}

#[test]
fn a_c4_uhf_rejects_a_charge_larger_than_the_nuclear_charge() {
    let (elements, positions) = water();
    // Water has 10 electrons; a charge of +99 is impossible.
    let err = sci_form::compute_uhf(&elements, &positions, 99, 1)
        .expect_err("an impossible charge must be rejected");
    assert!(err.contains("charge=99"), "unexpected error: {}", err);
    assert!(
        err.contains("incompatible with 10 nuclear charges"),
        "unexpected error: {}",
        err
    );
}

#[test]
fn a_c4_uhf_rejects_an_incompatible_multiplicity() {
    let (elements, positions) = water();
    // 10 electrons, multiplicity 4 → 3 unpaired electrons, 10 ± 3 is odd.
    let err = sci_form::compute_uhf(&elements, &positions, 0, 4)
        .expect_err("an incompatible multiplicity must be rejected");
    assert!(err.contains("multiplicity=4"), "unexpected error: {}", err);
    assert!(
        err.contains("incompatible with 10 electrons"),
        "unexpected error: {}",
        err
    );
}

#[test]
fn a_c4_uhf_rejects_a_zero_multiplicity() {
    let (elements, positions) = water();
    let err = sci_form::compute_uhf(&elements, &positions, 0, 0)
        .expect_err("a zero multiplicity must be rejected");
    assert!(
        err.contains("multiplicity must be >= 1"),
        "unexpected error: {}",
        err
    );
}

#[test]
fn a_c4_molecular_system_electron_count_never_wraps() {
    use sci_form::scf::types::MolecularSystem;
    let system =
        MolecularSystem::from_angstrom(&[1, 1], &[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]], 99, 1);
    assert_eq!(
        system.n_electrons(),
        0,
        "a negative electron count must not wrap around to a huge usize"
    );
}

// ─── A-C5 / A-C6 · GFN0 helpers report errors instead of panicking ──────────

#[test]
fn a_c5_gfn0_repulsion_rejects_an_unsupported_element() {
    let err = sci_form::xtb::gfn0::gfn0_repulsion(&[0], &[[0.0, 0.0, 0.0]])
        .expect_err("Z=0 has no GFN0 parameters");
    assert!(
        err.contains("GFN0 parameters not available for Z=0"),
        "unexpected error: {}",
        err
    );
}

#[test]
fn a_c5_gfn0_eeq_rejects_an_unsupported_element() {
    let err = sci_form::xtb::gfn0::gfn0_eeq(&[0], &[[0.0, 0.0, 0.0]], &[0.0], 0.0)
        .expect_err("Z=0 has no GFN0 parameters");
    assert!(
        err.contains("GFN0 parameters not available for Z=0"),
        "unexpected error: {}",
        err
    );
}

#[test]
fn a_c5_gfn0_h0_rejects_an_unsupported_element() {
    let err = match sci_form::xtb::gfn0::gfn0_h0(&[0], &[[0.0, 0.0, 0.0]], &[0.0], &[0.0]) {
        Ok(_) => panic!("Z=0 has no GFN0 parameters, the call must fail"),
        Err(e) => e,
    };
    assert!(
        err.contains("GFN0 parameters not available for Z=0"),
        "unexpected error: {}",
        err
    );
}

#[test]
fn a_c6_solve_gfn0_rejects_an_empty_molecule() {
    let err = sci_form::xtb::gfn0::solve_gfn0(&[], &[], 0.0)
        .expect_err("an empty molecule must be rejected, not panic");
    assert!(
        err.contains("empty geometry") || err.contains("empty molecule"),
        "unexpected error: {}",
        err
    );
}

// ─── A-C7 · PM3 length mismatch and empty molecule ──────────────────────────

#[test]
fn a_c7_pm3_rejects_a_length_mismatch() {
    let (elements, positions) = water();
    let err = sci_form::compute_pm3(&elements, &positions[..2])
        .expect_err("a length mismatch must be rejected");
    assert!(
        err.contains("positions has 2 entries"),
        "unexpected error: {}",
        err
    );
}

#[test]
fn a_c7_pm3_rejects_an_empty_molecule() {
    let err = sci_form::compute_pm3(&[], &[]).expect_err("an empty molecule must be rejected");
    assert!(err.contains("empty geometry"), "unexpected error: {}", err);
}

// ─── A-C8 · empty results are replaced by explicit errors ───────────────────

#[test]
fn a_c8_dos_rejects_zero_points() {
    let (elements, positions) = water();
    let err = sci_form::compute_dos(&elements, &positions, 0.3, -30.0, 5.0, 0)
        .expect_err("n_points = 0 must be rejected");
    assert!(
        err.contains("n_points must be >= 2"),
        "unexpected error: {}",
        err
    );
}

#[test]
fn a_c8_low_level_dos_rejects_zero_points() {
    let err = sci_form::dos::compute_dos(&[0.0], 0.1, -1.0, 1.0, 0)
        .expect_err("n_points = 0 must be rejected");
    assert!(
        err.contains("n_points must be >= 2"),
        "unexpected error: {}",
        err
    );
}

#[test]
fn a_c8_low_level_dos_rejects_a_zero_sigma() {
    let err = sci_form::dos::compute_dos(&[0.0], 0.0, -1.0, 1.0, 201)
        .expect_err("sigma = 0 must be rejected");
    assert!(
        err.contains("sigma must be a finite positive number"),
        "unexpected error: {}",
        err
    );
}

#[test]
fn a_c8_orbital_mesh_rejects_a_zero_spacing() {
    let (elements, positions) = water();
    let err = sci_form::compute_orbital_mesh(&elements, &positions, "eht", 0, 0.0, 3.0, 0.02)
        .expect_err("a zero spacing must be rejected instead of returning an empty grid");
    assert!(
        err.contains("spacing must be a finite positive number"),
        "unexpected error: {}",
        err
    );
}

#[test]
fn a_c8_orbital_mesh_rejects_an_out_of_range_mo_index() {
    let (elements, positions) = water();
    let err = sci_form::compute_orbital_mesh(&elements, &positions, "eht", 9999, 0.4, 3.0, 0.02)
        .expect_err("an out-of-range mo_index must be rejected");
    assert!(
        err.contains("mo_index 9999 is out of range"),
        "unexpected error: {}",
        err
    );
}

#[test]
fn a_c8_low_level_dos_pdos_rejects_a_length_mismatch() {
    let err = sci_form::dos::compute_pdos(
        &[1, 1],
        &[0.0; 3],
        &[0.0],
        &[vec![0.0]],
        2,
        0.3,
        -1.0,
        1.0,
        10,
    )
    .expect_err("a length mismatch must be rejected");
    assert!(
        err.contains("coords has 3 values"),
        "unexpected error: {}",
        err
    );
}

// ─── A-C9 · ANI must not abort on a non-finite coordinate ───────────────────

#[test]
fn a_c9_ani_cell_list_rejects_a_non_finite_coordinate() {
    let err = match sci_form::ani::neighbor::CellList::try_new(
        &[[0.0, 0.0, 0.0], [f64::INFINITY, 0.0, 0.0]],
        5.2,
    ) {
        Ok(_) => panic!("a non-finite coordinate must be rejected"),
        Err(e) => e,
    };
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

#[test]
fn a_c9_ani_cell_list_never_panics_on_a_non_finite_coordinate() {
    // The infallible constructor must not panic either.
    let cl =
        sci_form::ani::neighbor::CellList::new(&[[0.0, 0.0, 0.0], [f64::INFINITY, 0.0, 0.0]], 5.2);
    let pairs = cl.find_neighbors(&[[0.0, 0.0, 0.0], [f64::INFINITY, 0.0, 0.0]]);
    assert!(pairs.is_empty());
}

// ─── A-C10 · `embed_batch` must not hang on a huge thread count ─────────────

#[test]
fn a_c10_embed_batch_clamps_the_thread_count() {
    let config = sci_form::ConformerConfig {
        seed: 42,
        num_threads: 1 << 20,
    };
    let results = sci_form::embed_batch(&["CCO"], &config)
        .expect("embed_batch must not fail because of an absurd thread count");
    assert_eq!(results.len(), 1);
}

// ─── A-C11 · Gasteiger configuration validation ─────────────────────────────

#[test]
fn a_c11_gasteiger_rejects_zero_max_iter() {
    let config = GasteigerConfig {
        max_iter: 0,
        ..GasteigerConfig::default()
    };
    let err = gasteiger_marsili_charges_configured(&[1, 1], &[(0, 1)], &[0, 0], &config)
        .expect_err("max_iter = 0 must be rejected");
    assert_eq!(err, "max_iter must be >= 1 (got 0)");
}

#[test]
fn a_c11_gasteiger_rejects_a_zero_initial_damping() {
    let config = GasteigerConfig {
        initial_damping: 0.0,
        ..GasteigerConfig::default()
    };
    let err = gasteiger_marsili_charges_configured(&[1, 1], &[(0, 1)], &[0, 0], &config)
        .expect_err("initial_damping = 0 must be rejected");
    assert_eq!(err, "initial_damping must be in (0, 1] (got 0)");
}

#[test]
fn a_c11_gasteiger_rejects_an_out_of_range_initial_damping() {
    let config = GasteigerConfig {
        initial_damping: 1.5,
        ..GasteigerConfig::default()
    };
    let err = gasteiger_marsili_charges_configured(&[1, 1], &[(0, 1)], &[0, 0], &config)
        .expect_err("initial_damping > 1 must be rejected");
    assert_eq!(err, "initial_damping must be in (0, 1] (got 1.5)");
}

#[test]
fn a_c11_gasteiger_still_computes_with_a_valid_config() {
    let config = GasteigerConfig::default();
    let result = gasteiger_marsili_charges_configured(&[1, 1], &[(0, 1)], &[0, 0], &config)
        .expect("the default configuration must keep working");
    assert_eq!(result.charges.len(), 2);
}

// ─── A-M4 · the SMARTS parser must report parse errors, never panic ─────────

#[test]
fn a_m4_smarts_parser_reports_errors_instead_of_panicking() {
    let malformed = [
        "[", "]", "[C", "[C:", "C:", "(C", "C)", "[$(]", "[C;]", "[#6", "[!]", "[]",
    ];
    for pattern in malformed {
        // The only requirement is that the parser returns instead of
        // panicking; both `Ok` and `Err` are acceptable outcomes.
        let _ = sci_form::smarts::parse_smarts(pattern);
    }
}

#[test]
fn a_m4_smarts_parser_reports_an_error_for_an_unterminated_bracket() {
    let err = sci_form::smarts::parse_smarts("[C")
        .expect_err("an unterminated bracket must be a parse error");
    assert!(!err.is_empty());
}

// ─── Powder XRD and ML descriptors ──────────────────────────────────────────

#[test]
fn powder_xrd_rejects_a_length_mismatch() {
    let cell = sci_form::materials::UnitCell::cubic(4.0);
    let err = sci_form::materials::simulate_powder_xrd(&cell, &[6, 6], &[[0.0, 0.0, 0.0]], 90.0)
        .expect_err("a length mismatch must be rejected");
    assert!(
        err.contains("positions has 1 entries"),
        "unexpected error: {}",
        err
    );
}

#[test]
fn powder_xrd_rejects_a_non_finite_fractional_coordinate() {
    let cell = sci_form::materials::UnitCell::cubic(4.0);
    let err = sci_form::materials::simulate_powder_xrd(&cell, &[6], &[[f64::NAN, 0.0, 0.0]], 90.0)
        .expect_err("a non-finite coordinate must be rejected");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

#[test]
fn ml_3d_descriptors_reject_a_coordinate_length_mismatch() {
    let err = sci_form::ml::compute_3d_descriptors(&[6, 1, 1], &[0.0; 6])
        .expect_err("a length mismatch must be rejected, not silently zeroed");
    assert_eq!(
        err,
        "coords has 6 values but 3 atoms require 9 (3 atoms x 3 coordinates)"
    );
}

#[test]
fn ml_3d_descriptors_reject_a_non_finite_coordinate() {
    let err = sci_form::ml::compute_3d_descriptors(&[6], &[f64::NAN, 0.0, 0.0])
        .expect_err("a non-finite coordinate must be rejected");
    assert!(err.contains("non-finite"), "unexpected error: {}", err);
}

// ─── NPA input validation ───────────────────────────────────────────────────

#[test]
fn npa_rejects_an_out_of_range_basis_atom_map() {
    use nalgebra::DMatrix;
    let overlap = DMatrix::<f64>::identity(2, 2);
    let density = DMatrix::<f64>::identity(2, 2);
    let err = sci_form::population::compute_npa(&[6], &overlap, &density, &[0, 7])
        .expect_err("an out-of-range basis_atom_map entry must be rejected");
    assert!(
        err.contains("references an atom outside"),
        "unexpected error: {}",
        err
    );
}
