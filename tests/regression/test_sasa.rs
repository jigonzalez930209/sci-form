//! Integration tests for solvent-accessible surface area (SASA).
//!
//! Tests the public API `compute_sasa()` with various molecular geometries.

#[test]
fn test_single_carbon_sasa() {
    // Single carbon atom: SASA = 4π(r_vdw + r_probe)²
    let elements = vec![6u8];
    let coords = vec![0.0, 0.0, 0.0];
    let result = sci_form::compute_sasa(&elements, &coords, None).expect("should succeed");
    // C vdW = 1.7 Å, probe = 1.4 Å → r = 3.1 → SASA ≈ 4π(3.1)² ≈ 120.8 Ų
    assert!(
        result.total_sasa > 100.0 && result.total_sasa < 140.0,
        "single C SASA should be ~121 Ų, got {}",
        result.total_sasa
    );
    assert_eq!(result.atom_sasa.len(), 1);
}

#[test]
fn test_two_bonded_atoms_less_than_two_isolated() {
    // Two close atoms should have less SASA than two isolated atoms
    let elements = vec![6u8, 6];
    let close_coords = vec![0.0, 0.0, 0.0, 1.5, 0.0, 0.0];
    let far_coords = vec![0.0, 0.0, 0.0, 100.0, 0.0, 0.0];

    let close = sci_form::compute_sasa(&elements, &close_coords, None).unwrap();
    let far = sci_form::compute_sasa(&elements, &far_coords, None).unwrap();

    assert!(
        close.total_sasa < far.total_sasa,
        "bonded SASA ({}) should be less than isolated ({})",
        close.total_sasa,
        far.total_sasa
    );
}

#[test]
fn test_water_molecule_sasa() {
    let elements = vec![8u8, 1, 1];
    let coords = vec![0.0, 0.0, 0.0, 0.96, 0.0, 0.0, -0.24, 0.93, 0.0];
    let result = sci_form::compute_sasa(&elements, &coords, None).unwrap();
    assert!(
        result.total_sasa > 30.0 && result.total_sasa < 200.0,
        "water SASA should be reasonable, got {}",
        result.total_sasa
    );
    assert_eq!(result.atom_sasa.len(), 3);
    assert!((result.probe_radius - 1.4).abs() < 1e-10);
}

#[test]
fn test_custom_probe_radius() {
    let elements = vec![6u8];
    let coords = vec![0.0, 0.0, 0.0];
    let small = sci_form::compute_sasa(&elements, &coords, Some(0.5)).unwrap();
    let large = sci_form::compute_sasa(&elements, &coords, Some(2.0)).unwrap();
    assert!(
        small.total_sasa < large.total_sasa,
        "larger probe → larger SASA sphere"
    );
}

#[test]
fn test_coords_length_mismatch() {
    let elements = vec![6u8, 7];
    let coords = vec![0.0, 0.0, 0.0]; // only 3 values for 2 atoms
    let result = sci_form::compute_sasa(&elements, &coords, None);
    assert!(result.is_err(), "mismatched lengths should error");
}

#[test]
fn test_sasa_per_atom_sums_to_total() {
    let elements = vec![6u8, 7, 8];
    let coords = vec![0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 3.0, 0.0, 0.0];
    let result = sci_form::compute_sasa(&elements, &coords, None).unwrap();
    let sum: f64 = result.atom_sasa.iter().sum();
    assert!(
        (sum - result.total_sasa).abs() < 1e-6,
        "atom SASA sum ({}) should equal total ({})",
        sum,
        result.total_sasa
    );
}
