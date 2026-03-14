//! Integration tests for Gasteiger-Marsili partial charges.
//!
//! Tests the full SMILES → parse → charges pipeline via the public API.

#[test]
fn test_methane_charges_neutral() {
    let result = sci_form::compute_charges("C").expect("methane should succeed");
    // CH4: 5 atoms (C + 4 implicit H → explicit after parsing)
    assert!(!result.charges.is_empty(), "should have at least 1 atom");
    // Total charge should be near zero
    assert!(
        result.total_charge.abs() < 1e-6,
        "total charge should be ~0, got {}",
        result.total_charge
    );
}

#[test]
fn test_ethanol_charge_polarisation() {
    let result = sci_form::compute_charges("CCO").expect("ethanol should succeed");
    // Oxygen (element 8) should be the most negative atom
    // Parse CCO → C(0), C(1), O(2), then Hs
    // O at index 2 should have negative charge
    assert!(
        result.charges[2] < 0.0,
        "oxygen should have negative charge, got {}",
        result.charges[2]
    );
    // Total charge near zero
    assert!(
        result.total_charge.abs() < 1e-6,
        "total charge should be ~0, got {}",
        result.total_charge
    );
}

#[test]
fn test_acetic_acid_charges() {
    let result = sci_form::compute_charges("CC(=O)O").expect("acetic acid should succeed");
    assert!(result.charges.len() >= 4, "need at least 4 heavy atoms");
    assert!(result.total_charge.abs() < 1e-6);
    assert_eq!(result.iterations, 6);
}

#[test]
fn test_benzene_symmetric_charges() {
    let result = sci_form::compute_charges("c1ccccc1").expect("benzene should succeed");
    // All carbons should have similar charges (symmetric molecule)
    let carbon_charges: Vec<f64> = result.charges[..6].to_vec();
    let mean = carbon_charges.iter().sum::<f64>() / 6.0;
    for (i, &q) in carbon_charges.iter().enumerate() {
        assert!(
            (q - mean).abs() < 0.01,
            "carbon {} charge {} deviates from mean {}",
            i,
            q,
            mean
        );
    }
}

#[test]
fn test_charged_molecule() {
    // Ammonium: [NH4+]
    let result = sci_form::compute_charges("[NH4+]").expect("ammonium should succeed");
    // Total charge should be close to +1
    assert!(
        (result.total_charge - 1.0).abs() < 0.1,
        "ammonium total charge should be ~1, got {}",
        result.total_charge
    );
}

#[test]
fn test_invalid_smiles_returns_error() {
    let result = sci_form::compute_charges("NOT_SMILES");
    assert!(result.is_err(), "invalid SMILES should return error");
}
