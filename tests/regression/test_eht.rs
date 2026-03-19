//! Integration tests for the EHT (Extended Hückel Theory) pipeline.
//!
//! Tests the full SMILES → EHT → orbitals → volume → mesh pipeline.

use sci_form::eht::{
    self, evaluate_orbital_on_grid, marching_cubes, solve_eht, EhtResult, IsosurfaceMesh,
};

// ─── Helper: build molecule from elements + positions ────────────────────────

fn h2_molecule() -> (Vec<u8>, Vec<[f64; 3]>) {
    // H₂ at ~0.74 Å bond length
    let elements = vec![1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
    (elements, positions)
}

fn water_molecule() -> (Vec<u8>, Vec<[f64; 3]>) {
    // H₂O with O at origin, H-O-H angle ~104.5°
    let elements = vec![8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
    (elements, positions)
}

fn methane_molecule() -> (Vec<u8>, Vec<[f64; 3]>) {
    // CH₄ tetrahedral geometry
    let elements = vec![6, 1, 1, 1, 1];
    let positions = vec![
        [0.0, 0.0, 0.0],
        [0.629, 0.629, 0.629],
        [-0.629, -0.629, 0.629],
        [-0.629, 0.629, -0.629],
        [0.629, -0.629, -0.629],
    ];
    (elements, positions)
}

// ─── Full pipeline tests ─────────────────────────────────────────────────────

#[test]
fn test_h2_full_pipeline() {
    let (elems, pos) = h2_molecule();
    let result = solve_eht(&elems, &pos, None).expect("H₂ EHT should succeed");

    // H₂ has 2 basis functions → 2 MOs
    assert_eq!(result.energies.len(), 2, "H₂ should have 2 MOs");

    // HOMO < LUMO
    assert!(
        result.homo_energy < result.lumo_energy,
        "HOMO ({}) should be below LUMO ({})",
        result.homo_energy,
        result.lumo_energy
    );

    // Gap should be positive
    assert!(
        result.gap > 0.0,
        "H₂ gap should be positive: {}",
        result.gap
    );

    // 2 electrons → HOMO index = 0
    assert_eq!(result.homo_index, 0);
    assert_eq!(result.lumo_index, 1);
}

#[test]
fn test_water_full_pipeline() {
    let (elems, pos) = water_molecule();
    let result = solve_eht(&elems, &pos, None).expect("H₂O EHT should succeed");

    // O has 4 basis functions (2s, 2px, 2py, 2pz), each H has 1 → 6 MOs
    assert_eq!(result.energies.len(), 6, "H₂O should have 6 MOs");

    // 8 valence electrons → 4 occupied MOs
    assert_eq!(result.n_electrons, 8);
    assert_eq!(result.homo_index, 3); // 0,1,2,3 occupied
    assert_eq!(result.lumo_index, 4);

    // Energies sorted ascending
    for i in 1..result.energies.len() {
        assert!(
            result.energies[i] >= result.energies[i - 1],
            "Energies not sorted: {} > {} at index {}",
            result.energies[i - 1],
            result.energies[i],
            i
        );
    }
}

#[test]
fn test_methane_full_pipeline() {
    let (elems, pos) = methane_molecule();
    let result = solve_eht(&elems, &pos, None).expect("CH₄ EHT should succeed");

    // C: 4 basis (2s,2px,2py,2pz) + 4H: 4×1s = 8 MOs
    assert_eq!(result.energies.len(), 8, "CH₄ should have 8 MOs");

    // 4 + 4×1 = 8 valence electrons → 4 occupied MOs
    assert_eq!(result.n_electrons, 8);
    assert_eq!(result.homo_index, 3);
    assert_eq!(result.lumo_index, 4);

    // First MO should be bonding (lower energy than atomic)
    assert!(
        result.energies[0] < -10.0,
        "Lowest MO should be deep bonding: {}",
        result.energies[0]
    );
}

// ─── Orbital visualization pipeline ─────────────────────────────────────────

#[test]
fn test_h2_orbital_to_mesh() {
    let (elems, pos) = h2_molecule();
    let result = solve_eht(&elems, &pos, None).unwrap();

    let basis = eht::basis::build_basis(&elems, &pos);

    // Generate volumetric grid for bonding orbital (MO 0)
    let grid = evaluate_orbital_on_grid(&basis, &result.coefficients, 0, &pos, 0.2, 3.0);

    assert!(!grid.values.is_empty(), "Grid should have values");
    assert_eq!(
        grid.values.len(),
        grid.dims[0] * grid.dims[1] * grid.dims[2],
        "Grid size mismatch"
    );

    // Bonding orbital should have nonzero values between the atoms
    let max_val = grid
        .values
        .iter()
        .cloned()
        .fold(f32::NEG_INFINITY, f32::max);
    let min_val = grid.values.iter().cloned().fold(f32::INFINITY, f32::min);
    let max_abs = max_val.abs().max(min_val.abs());
    assert!(
        max_abs > 1e-6,
        "Bonding orbital should have nonzero density, max_abs={}",
        max_abs
    );

    // Extract isosurface at a fraction of the peak value (accounting for sign)
    let iso = if max_val.abs() > min_val.abs() {
        max_val * 0.5
    } else {
        min_val * 0.5
    };
    let mesh = marching_cubes(&grid, iso);
    assert!(
        mesh.num_triangles > 0,
        "Bonding orbital isosurface should produce triangles"
    );
    assert_eq!(mesh.vertices.len(), mesh.normals.len());
    assert_eq!(mesh.indices.len(), mesh.num_triangles * 3);
}

#[test]
fn test_water_homo_visualization() {
    let (elems, pos) = water_molecule();
    let result = solve_eht(&elems, &pos, None).unwrap();
    let basis = eht::basis::build_basis(&elems, &pos);

    let homo = result.homo_index;
    let grid = evaluate_orbital_on_grid(&basis, &result.coefficients, homo, &pos, 0.25, 3.0);

    // HOMO should have both positive and negative lobes
    let has_positive = grid.values.iter().any(|&v| v > 0.01);
    let has_negative = grid.values.iter().any(|&v| v < -0.01);

    assert!(
        has_positive || has_negative,
        "HOMO should have significant orbital density"
    );

    // Extract both positive and negative isosurfaces
    let pos_mesh = marching_cubes(&grid, 0.02);
    let neg_mesh = marching_cubes(&grid, -0.02);

    // At least one lobe should produce triangles
    assert!(
        pos_mesh.num_triangles > 0 || neg_mesh.num_triangles > 0,
        "HOMO should produce at least one isosurface lobe"
    );
}

// ─── Serialization round-trip ────────────────────────────────────────────────

#[test]
fn test_eht_result_serialization() {
    let (elems, pos) = h2_molecule();
    let result = solve_eht(&elems, &pos, None).unwrap();

    let json = serde_json::to_string(&result).expect("EhtResult should serialize");
    let deser: EhtResult = serde_json::from_str(&json).expect("EhtResult should deserialize");

    assert_eq!(result.energies.len(), deser.energies.len());
    assert_eq!(result.n_electrons, deser.n_electrons);
    assert_eq!(result.homo_index, deser.homo_index);
    assert_eq!(result.lumo_index, deser.lumo_index);
}

#[test]
fn test_mesh_serialization() {
    let (elems, pos) = h2_molecule();
    let result = solve_eht(&elems, &pos, None).unwrap();
    let basis = eht::basis::build_basis(&elems, &pos);
    let grid = evaluate_orbital_on_grid(&basis, &result.coefficients, 0, &pos, 0.3, 3.0);
    let mesh = marching_cubes(&grid, 0.02);

    let json = serde_json::to_string(&mesh).expect("IsosurfaceMesh should serialize");
    let deser: IsosurfaceMesh =
        serde_json::from_str(&json).expect("IsosurfaceMesh should deserialize");
    assert_eq!(mesh.num_triangles, deser.num_triangles);
    assert_eq!(mesh.vertices.len(), deser.vertices.len());
}

// ─── Edge cases ──────────────────────────────────────────────────────────────

#[test]
fn test_single_atom_eht() {
    // Single carbon atom
    let elems = vec![6];
    let pos = vec![[0.0, 0.0, 0.0]];
    let result = solve_eht(&elems, &pos, None).expect("Single atom should succeed");

    // C has 4 basis functions
    assert_eq!(result.energies.len(), 4);
    assert_eq!(result.n_electrons, 4); // C has 4 valence electrons
}

#[test]
fn test_unsupported_element_error() {
    // Helium (Z=2) is not in the parameter table
    let elems = vec![2];
    let pos = vec![[0.0, 0.0, 0.0]];
    let result = solve_eht(&elems, &pos, None);
    assert!(result.is_err(), "Unsupported element should return error");
}
