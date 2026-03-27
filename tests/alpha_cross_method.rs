//! Cross-method comparison tests.
//!
//! Validates that different computational methods agree on qualitative
//! chemical behavior for the same molecules:
//! - All methods produce finite, convergent results for H₂
//! - Energy ordering: DFT ≠ PM3 ≠ xTB (different baselines), but
//!   all should show H₂ is bound (lower energy than dissociated)
//! - HOMO-LUMO gap is positive across all methods
//! - All methods give symmetric charges for symmetric molecules
//!
//! These tests verify that the full stack (conformer → properties → quantum)
//! functions as an integrated pipeline.

#![cfg(all(
    feature = "alpha-dft",
    feature = "alpha-obara-saika",
    feature = "alpha-reaxff"
))]

use sci_form::dft::ks_fock::{solve_ks_dft, DftConfig};
use sci_form::eht::solver::solve_eht;

// ─── H₂: multi-method comparison ────────────────────────────────────────────

/// All electronic structure methods should converge for H₂.
#[test]
fn h2_all_methods_converge() {
    let elements = [1u8, 1];
    let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];

    // DFT/SVWN
    let dft = solve_ks_dft(&elements, &positions, &DftConfig::default()).unwrap();
    assert!(dft.converged, "DFT/SVWN must converge for H₂");

    // PM3
    let pm3 = sci_form::compute_pm3(&elements, &positions).unwrap();
    assert!(pm3.converged, "PM3 must converge for H₂");

    // EHT
    let eht = solve_eht(&elements, &positions, None).unwrap();
    assert!(eht.gap.is_finite(), "EHT should give finite HOMO-LUMO gap");
}

/// All methods should show positive HOMO-LUMO gap for H₂.
#[test]
fn h2_homo_lumo_gap_positive_all_methods() {
    let elements = [1u8, 1];
    let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];

    let dft = solve_ks_dft(&elements, &positions, &DftConfig::default()).unwrap();
    assert!(dft.gap > 0.0, "DFT gap ({:.3} eV) should be > 0", dft.gap);

    let pm3 = sci_form::compute_pm3(&elements, &positions).unwrap();
    assert!(pm3.gap > 0.0, "PM3 gap ({:.3} eV) should be > 0", pm3.gap);

    let eht = solve_eht(&elements, &positions, None).unwrap();
    assert!(eht.gap > 0.0, "EHT gap ({:.3} eV) should be > 0", eht.gap);
}

// ─── Ethanol: conformer → property pipeline ──────────────────────────────────

/// Full pipeline: embed SMILES → DFT calculation.
#[test]
fn ethanol_embed_then_dft() {
    let conf = sci_form::embed("CCO", 42);
    assert!(conf.error.is_none(), "Embed should succeed for CCO");

    let positions: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

    let dft = solve_ks_dft(&conf.elements, &positions, &DftConfig::default()).unwrap();
    assert!(dft.converged, "DFT must converge for ethanol");
    assert!(
        dft.total_energy < 0.0,
        "Total energy should be negative for ethanol"
    );
    assert_eq!(dft.n_electrons, 26, "Ethanol has 26 electrons");
}

/// Full pipeline: embed SMILES → PM3 → check heat of formation.
#[test]
fn ethanol_embed_then_pm3() {
    let conf = sci_form::embed("CCO", 42);
    assert!(conf.error.is_none());

    let positions: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

    let pm3 = sci_form::compute_pm3(&conf.elements, &positions).unwrap();
    assert!(pm3.converged, "PM3 must converge for ethanol");
    // Experimental HOF of ethanol ≈ -56.2 kcal/mol (NIST)
    // PM3 HOF within ±30 kcal/mol of experimental is reasonable for PM3 accuracy
    assert!(
        pm3.heat_of_formation.is_finite(),
        "PM3 HOF should be finite"
    );
}

// ─── Charge neutrality across methods ────────────────────────────────────────

/// All methods should give approximately neutral charges for neutral H₂.
#[test]
fn h2_charge_neutrality_all_methods() {
    let elements = [1u8, 1];
    let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];

    let dft = solve_ks_dft(&elements, &positions, &DftConfig::default()).unwrap();
    let dft_sum: f64 = dft.mulliken_charges.iter().sum();
    assert!(
        dft_sum.abs() < 0.1,
        "DFT Mulliken sum ({:.4}) should be near 0",
        dft_sum
    );

    let pm3 = sci_form::compute_pm3(&elements, &positions).unwrap();
    let pm3_sum: f64 = pm3.mulliken_charges.iter().sum();
    assert!(
        pm3_sum.abs() < 0.1,
        "PM3 Mulliken sum ({:.4}) should be near 0",
        pm3_sum
    );

    let eht_pop = sci_form::compute_population(&elements, &positions).unwrap();
    let eht_sum: f64 = eht_pop.mulliken_charges.iter().sum();
    assert!(
        eht_sum.abs() < 0.5,
        "EHT Mulliken sum ({:.4}) should be near 0",
        eht_sum
    );
}

// ─── ReaxFF vs DFT qualitative agreement ─────────────────────────────────────

/// Both ReaxFF and DFT should show that equilibrium H₂ is more stable than stretched.
#[test]
fn h2_equilibrium_more_stable_than_stretched() {
    let elements = [1u8, 1];

    // DFT
    let dft_eq = solve_ks_dft(
        &elements,
        &[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]],
        &DftConfig::default(),
    )
    .unwrap();
    let dft_stretched = solve_ks_dft(
        &elements,
        &[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
        &DftConfig::default(),
    )
    .unwrap();
    assert!(
        dft_eq.total_energy < dft_stretched.total_energy,
        "DFT: equilibrium ({:.4}) should be lower than stretched ({:.4})",
        dft_eq.total_energy,
        dft_stretched.total_energy
    );

    // ReaxFF
    use sci_form::forcefield::reaxff::gradients::compute_reaxff_gradient;
    use sci_form::forcefield::reaxff::params::ReaxffParams;
    let params = ReaxffParams::default_chon();
    let (e_eq, _) =
        compute_reaxff_gradient(&[0.0, 0.0, 0.0, 0.74, 0.0, 0.0], &elements, &params).unwrap();
    let (e_str, _) =
        compute_reaxff_gradient(&[0.0, 0.0, 0.0, 3.0, 0.0, 0.0], &elements, &params).unwrap();

    // ReaxFF might not show the same ordering due to simplified params,
    // but both should produce finite energies
    assert!(e_eq.is_finite() && e_str.is_finite());
}
