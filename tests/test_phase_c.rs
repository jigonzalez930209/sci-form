use sci_form;

#[test]
fn test_population_water() {
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let result = sci_form::compute_population(&elements, &positions).unwrap();
    assert_eq!(result.num_atoms, 3);
    assert_eq!(result.mulliken_charges.len(), 3);
    assert_eq!(result.lowdin_charges.len(), 3);
    // Oxygen should have negative charge
    assert!(result.mulliken_charges[0] < 0.0, "O should be negative");
    // Hydrogens should be positive
    assert!(result.mulliken_charges[1] > 0.0, "H should be positive");
    assert!(result.mulliken_charges[2] > 0.0, "H should be positive");
}

#[test]
fn test_population_h2_symmetric() {
    let elements = vec![1u8, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
    let result = sci_form::compute_population(&elements, &positions).unwrap();
    let diff = (result.mulliken_charges[0] - result.mulliken_charges[1]).abs();
    assert!(diff < 0.01, "H₂ charges should be symmetric, diff={diff}");
}

#[test]
fn test_population_charge_conservation() {
    let elements = vec![6u8, 1, 1, 1, 1]; // CH4
    let positions = vec![
        [0.0, 0.0, 0.0],
        [0.63, 0.63, 0.63],
        [-0.63, -0.63, 0.63],
        [-0.63, 0.63, -0.63],
        [0.63, -0.63, -0.63],
    ];
    let result = sci_form::compute_population(&elements, &positions).unwrap();
    assert!(
        result.total_charge_mulliken.abs() < 0.5,
        "Neutral molecule charge should be near 0, got {}",
        result.total_charge_mulliken
    );
}

#[test]
fn test_dipole_water_nonzero() {
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let result = sci_form::compute_dipole(&elements, &positions).unwrap();
    assert!(result.magnitude > 0.1, "Water should have nonzero dipole");
    assert_eq!(result.unit, "Debye");
}

#[test]
fn test_dipole_h2_near_zero() {
    let elements = vec![1u8, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
    let result = sci_form::compute_dipole(&elements, &positions).unwrap();
    assert!(
        result.magnitude < 0.5,
        "H₂ dipole should be ~0 D, got {:.3}",
        result.magnitude
    );
}

#[test]
fn test_dos_water() {
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let result = sci_form::compute_dos(&elements, &positions, 0.3, -30.0, 5.0, 200).unwrap();
    assert_eq!(result.energies.len(), 200);
    assert_eq!(result.total_dos.len(), 200);
    assert_eq!(result.pdos.len(), 3); // 3 atoms
    // Total DOS should have some nonzero values
    let max_dos = result.total_dos.iter().cloned().fold(0.0f64, f64::max);
    assert!(max_dos > 0.1, "DOS should have peaks");
}

#[test]
fn test_rmsd_identical() {
    let coords = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let rmsd = sci_form::compute_rmsd(&coords, &coords);
    assert!(rmsd < 1e-10);
}

#[test]
fn test_rmsd_translated() {
    let reference = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let shifted: Vec<f64> = reference.iter().map(|x| x + 3.0).collect();
    let rmsd = sci_form::compute_rmsd(&shifted, &reference);
    assert!(rmsd < 1e-8, "Pure translation RMSD should be 0, got {rmsd}");
}

#[test]
fn test_uff_energy_methane() {
    // Methane with known tetrahedral geometry
    let smiles = "C";
    let coords = vec![
        0.0, 0.0, 0.0,         // C
        0.63, 0.63, 0.63,      // H
        -0.63, -0.63, 0.63,    // H
        -0.63, 0.63, -0.63,    // H
        0.63, -0.63, -0.63,    // H
    ];
    let energy = sci_form::compute_uff_energy(smiles, &coords).unwrap();
    assert!(energy.is_finite(), "UFF energy should be finite");
}

#[test]
fn test_esp_water() {
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let grid = sci_form::compute_esp(&elements, &positions, 0.5, 2.0).unwrap();
    assert!(grid.values.len() > 0);
    assert_eq!(grid.values.len(), grid.dims[0] * grid.dims[1] * grid.dims[2]);
    // No NaN values
    assert!(grid.values.iter().all(|v| v.is_finite()));
}
