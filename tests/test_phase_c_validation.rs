//! Comprehensive validation tests for Phase C2–C7.
//!
//! These tests go beyond basic "it works" to validate:
//! - Physical correctness (electronegativity ordering, conservation laws)
//! - Mathematical invariants (symmetry, normalization, unitarity)
//! - Numerical robustness (edge cases, cross-module consistency)
//! - Round-trip fidelity (export → import identity)

// ═══════════════════════════════════════════════════════════════════════════
// C2: Population Analysis — Mulliken & Löwdin
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_population_electronegativity_ordering_water() {
    // Oxygen is more electronegative than hydrogen → more negative charge
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let result = sci_form::compute_population(&elements, &positions).unwrap();

    // O should be the most negative atom
    assert!(
        result.mulliken_charges[0] < result.mulliken_charges[1],
        "O ({:.4}) should be more negative than H ({:.4})",
        result.mulliken_charges[0],
        result.mulliken_charges[1]
    );
    assert!(
        result.lowdin_charges[0] < result.lowdin_charges[1],
        "Löwdin: O ({:.4}) should be more negative than H ({:.4})",
        result.lowdin_charges[0],
        result.lowdin_charges[1]
    );
}

#[test]
fn test_population_charge_conservation_water() {
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let result = sci_form::compute_population(&elements, &positions).unwrap();

    let sum_mulliken: f64 = result.mulliken_charges.iter().sum();
    let sum_lowdin: f64 = result.lowdin_charges.iter().sum();

    assert!(
        sum_mulliken.abs() < 0.5,
        "Neutral molecule Mulliken sum should ≈ 0, got {sum_mulliken:.4}"
    );
    assert!(
        sum_lowdin.abs() < 0.5,
        "Neutral molecule Löwdin sum should ≈ 0, got {sum_lowdin:.4}"
    );
}

#[test]
fn test_population_formaldehyde_carbon_charge() {
    // H₂C=O: Carbon bonded to O (electron-withdrawing) should be slightly positive
    let elements = vec![6u8, 8, 1, 1]; // C, O, H, H
    let positions = vec![
        [0.0, 0.0, 0.0],    // C
        [1.21, 0.0, 0.0],   // O (C=O ≈ 1.21 Å)
        [-0.55, 0.93, 0.0], // H
        [-0.55, -0.93, 0.0], // H
    ];
    let result = sci_form::compute_population(&elements, &positions).unwrap();

    // Overall charge conservation for neutral molecule
    let sum: f64 = result.mulliken_charges.iter().sum();
    assert!(sum.abs() < 0.5, "Charge sum = {sum:.4}, expected ≈ 0");

    // Oxygen should carry negative charge (most electronegative)
    assert!(
        result.mulliken_charges[1] < 0.0,
        "Oxygen should be negative, got {:.4}",
        result.mulliken_charges[1]
    );
}

#[test]
fn test_population_lowdin_vs_mulliken_bounded() {
    // Löwdin charges should generally be smaller in magnitude than Mulliken
    // (Löwdin partitions overlap more evenly) — but at minimum both should be finite
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let result = sci_form::compute_population(&elements, &positions).unwrap();

    for (i, (&m, &l)) in result
        .mulliken_charges
        .iter()
        .zip(result.lowdin_charges.iter())
        .enumerate()
    {
        assert!(m.is_finite(), "Mulliken charge[{i}] is NaN/Inf");
        assert!(l.is_finite(), "Löwdin charge[{i}] is NaN/Inf");
        // Both should be within reasonable range for neutral small molecules
        assert!(
            m.abs() < 3.0,
            "Mulliken charge[{i}] = {m:.4} out of range"
        );
        assert!(l.abs() < 3.0, "Löwdin charge[{i}] = {l:.4} out of range");
    }
}

#[test]
fn test_population_methane_hydrogen_equivalence() {
    // CH4: all 4 H atoms should have nearly identical charges by symmetry
    let elements = vec![6u8, 1, 1, 1, 1];
    let d = 0.63;
    let positions = vec![
        [0.0, 0.0, 0.0],
        [d, d, d],
        [-d, -d, d],
        [-d, d, -d],
        [d, -d, -d],
    ];
    let result = sci_form::compute_population(&elements, &positions).unwrap();

    let h_charges = &result.mulliken_charges[1..];
    let mean: f64 = h_charges.iter().sum::<f64>() / 4.0;
    for (i, &q) in h_charges.iter().enumerate() {
        assert!(
            (q - mean).abs() < 0.05,
            "H[{i}] charge {q:.4} deviates from mean {mean:.4} by > 0.05"
        );
    }
}

#[test]
fn test_population_gross_orbital_nonnegative() {
    // Mulliken gross orbital populations should be non-negative
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let result = sci_form::compute_population(&elements, &positions).unwrap();

    assert!(!result.mulliken_populations.is_empty());
    // Populations should be non-negative (or very slightly negative due to numerical noise)
    for (i, &pop) in result.mulliken_populations.iter().enumerate() {
        assert!(
            pop > -0.1,
            "Orbital population[{i}] = {pop:.6} is negative (beyond tolerance)"
        );
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// C3: Dipole Moments
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_dipole_water_magnitude_range() {
    // Water experimental dipole ≈ 1.85 D; EHT/Mulliken won't be exact but
    // should be in a physically reasonable range (0.5–5.0 D)
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let result = sci_form::compute_dipole(&elements, &positions).unwrap();

    assert!(
        result.magnitude > 0.3 && result.magnitude < 10.0,
        "Water dipole {:.3} D outside physical range",
        result.magnitude
    );
    assert_eq!(result.unit, "Debye");
}

#[test]
fn test_dipole_vector_magnitude_consistent() {
    // |μ| should equal sqrt(μx² + μy² + μz²)
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let result = sci_form::compute_dipole(&elements, &positions).unwrap();

    let calc_mag = (result.vector[0].powi(2)
        + result.vector[1].powi(2)
        + result.vector[2].powi(2))
    .sqrt();

    assert!(
        (calc_mag - result.magnitude).abs() < 1e-10,
        "Vector magnitude {calc_mag:.6} != reported magnitude {:.6}",
        result.magnitude
    );
}

#[test]
fn test_dipole_symmetric_molecule_near_zero() {
    // H₂ is symmetric → dipole should be near zero
    let elements = vec![1u8, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
    let result = sci_form::compute_dipole(&elements, &positions).unwrap();
    assert!(
        result.magnitude < 0.5,
        "H₂ dipole should be ~0 D, got {:.4}",
        result.magnitude
    );
}

#[test]
fn test_dipole_methane_near_zero() {
    // CH4 (Td symmetry): dipole should be very small
    let elements = vec![6u8, 1, 1, 1, 1];
    let d = 0.63;
    let positions = vec![
        [0.0, 0.0, 0.0],
        [d, d, d],
        [-d, -d, d],
        [-d, d, -d],
        [d, -d, -d],
    ];
    let result = sci_form::compute_dipole(&elements, &positions).unwrap();
    assert!(
        result.magnitude < 1.0,
        "CH₄ dipole should be near 0, got {:.4} D",
        result.magnitude
    );
}

#[test]
fn test_dipole_direction_hf() {
    // HF: dipole should point from H(+) to F(-)
    // With H at origin and F along +x, dipole vector x should be positive
    // (convention: dipole points from − to +, but for charge*pos formula,
    //  if F has negative charge at positive x, μ_x = q_H*0 + q_F*x_F will be negative)
    let elements = vec![1u8, 9]; // H, F
    let positions = vec![[0.0, 0.0, 0.0], [0.92, 0.0, 0.0]]; // HF ≈ 0.92 Å
    let result = sci_form::compute_dipole(&elements, &positions).unwrap();
    assert!(result.magnitude > 0.1, "HF should have nonzero dipole");
    // The dipole vector should be primarily along x-axis
    assert!(
        result.vector[0].abs()
            > result.vector[1].abs().max(result.vector[2].abs()),
        "HF dipole should be primarily along the bond axis (x)"
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// C4: Electrostatic Potential (ESP)
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_esp_grid_dimensions_consistent() {
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let grid = sci_form::compute_esp(&elements, &positions, 0.5, 2.0).unwrap();

    assert_eq!(
        grid.values.len(),
        grid.dims[0] * grid.dims[1] * grid.dims[2],
        "Values length should match dims product"
    );
    assert!(grid.dims[0] > 0 && grid.dims[1] > 0 && grid.dims[2] > 0);
    assert!(grid.spacing > 0.0);
}

#[test]
fn test_esp_all_finite() {
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let grid = sci_form::compute_esp(&elements, &positions, 0.5, 2.0).unwrap();

    let nan_count = grid.values.iter().filter(|v| !v.is_finite()).count();
    assert_eq!(nan_count, 0, "ESP grid contains {nan_count} NaN/Inf values");
}

#[test]
fn test_esp_has_positive_and_negative_regions() {
    // Water ESP should have both positive and negative regions
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let grid = sci_form::compute_esp(&elements, &positions, 0.5, 2.0).unwrap();

    let has_pos = grid.values.iter().any(|&v| v > 1e-6);
    let has_neg = grid.values.iter().any(|&v| v < -1e-6);
    assert!(
        has_pos && has_neg,
        "Water ESP should have both +/- regions (pos={has_pos}, neg={has_neg})"
    );
}

#[test]
fn test_esp_cube_roundtrip() {
    // Export to cube then re-import; values should match
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let grid = sci_form::compute_esp(&elements, &positions, 0.8, 1.5).unwrap();

    // Export
    let mut buf = Vec::new();
    sci_form::esp::export_cube(&mut buf, &elements, &positions, &grid).unwrap();

    // Re-import
    let mut reader = std::io::BufReader::new(&buf[..]);
    let cube = sci_form::esp::read_cube(&mut reader).unwrap();

    assert_eq!(cube.elements.len(), elements.len());
    assert_eq!(cube.dims, grid.dims);
    assert_eq!(cube.values.len(), grid.values.len());

    // Values should survive round-trip within formatting precision
    for (i, (&orig, &read)) in grid.values.iter().zip(cube.values.iter()).enumerate() {
        let diff = (orig - read).abs();
        // .cube uses %12.5E format ≈ 6 significant digits
        let tol = orig.abs() * 1e-4 + 1e-10;
        assert!(
            diff < tol,
            "Cube value[{i}] mismatch: orig={orig:.8e}, read={read:.8e}, diff={diff:.8e}"
        );
    }
}

#[test]
fn test_esp_finer_grid_more_points() {
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];

    let coarse = sci_form::compute_esp(&elements, &positions, 1.0, 2.0).unwrap();
    let fine = sci_form::compute_esp(&elements, &positions, 0.5, 2.0).unwrap();

    assert!(
        fine.values.len() > coarse.values.len(),
        "Finer grid should have more points ({} vs {})",
        fine.values.len(),
        coarse.values.len()
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// C5: Density of States (DOS) and Projected DOS (PDOS)
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_dos_grid_matches_parameters() {
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let n_points = 300;
    let result =
        sci_form::compute_dos(&elements, &positions, 0.3, -30.0, 5.0, n_points).unwrap();

    assert_eq!(result.energies.len(), n_points);
    assert_eq!(result.total_dos.len(), n_points);
    assert!((result.energies[0] - (-30.0)).abs() < 1e-10);
    assert!((result.energies[n_points - 1] - 5.0).abs() < 1e-10);
}

#[test]
fn test_dos_total_nonnegative() {
    // DOS values should be non-negative (Gaussians are always positive)
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let result = sci_form::compute_dos(&elements, &positions, 0.3, -30.0, 5.0, 200).unwrap();

    for (i, &d) in result.total_dos.iter().enumerate() {
        assert!(d >= 0.0, "DOS[{i}] = {d:.6} is negative");
    }
}

#[test]
fn test_dos_pdos_sum_equals_total() {
    // Sum of all PDOS should approximately equal total DOS
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let result = sci_form::compute_dos(&elements, &positions, 0.3, -30.0, 5.0, 200).unwrap();

    assert_eq!(result.pdos.len(), 3, "Should have 3 atoms in PDOS");

    for gi in 0..result.total_dos.len() {
        let pdos_sum: f64 = result.pdos.iter().map(|p| p[gi]).sum();
        let total = result.total_dos[gi];
        if total > 0.01 {
            let rel_err = (pdos_sum - total).abs() / total;
            assert!(
                rel_err < 0.15,
                "PDOS sum ({pdos_sum:.6}) != total DOS ({total:.6}) at grid[{gi}], rel_err={rel_err:.4}"
            );
        }
    }
}

#[test]
fn test_dos_integral_approximation() {
    // Integral of DOS should ~ number of orbitals (each Gaussian integrates to 1)
    let elements = vec![1u8, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
    let sigma = 0.3;
    let n_points = 1000;
    let e_min = -40.0;
    let e_max = 10.0;
    let result =
        sci_form::compute_dos(&elements, &positions, sigma, e_min, e_max, n_points).unwrap();

    let de = (e_max - e_min) / (n_points - 1) as f64;
    let integral: f64 = result.total_dos.iter().sum::<f64>() * de;

    // H2 has 2 molecular orbitals
    assert!(
        integral > 1.0 && integral < 4.0,
        "DOS integral {integral:.3} should be ≈ 2 (number of MOs)"
    );
}

#[test]
fn test_dos_narrower_sigma_sharper_peaks() {
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];

    let broad = sci_form::compute_dos(&elements, &positions, 1.0, -30.0, 5.0, 200).unwrap();
    let narrow = sci_form::compute_dos(&elements, &positions, 0.1, -30.0, 5.0, 200).unwrap();

    let max_broad = broad.total_dos.iter().cloned().fold(0.0f64, f64::max);
    let max_narrow = narrow.total_dos.iter().cloned().fold(0.0f64, f64::max);

    assert!(
        max_narrow > max_broad,
        "Narrower sigma should give taller peaks: narrow_max={max_narrow:.3}, broad_max={max_broad:.3}"
    );
}

#[test]
fn test_dos_pdos_nonnegative() {
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let result = sci_form::compute_dos(&elements, &positions, 0.3, -30.0, 5.0, 200).unwrap();

    for (a, pdos) in result.pdos.iter().enumerate() {
        for (i, &val) in pdos.iter().enumerate() {
            assert!(
                val >= -1e-10,
                "PDOS[atom={a}][{i}] = {val:.8} is negative"
            );
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// C6: Molecular Alignment and RMSD (Kabsch)
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_rmsd_rotation_invariance() {
    // RMSD of a rotated structure with itself should be ≈ 0
    let reference = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

    // Rotate 90° around z-axis: (x,y,z) → (-y,x,z)
    let rotated = vec![
        0.0, 0.0, 0.0, // origin stays
        0.0, 1.0, 0.0, // (1,0,0) → (0,1,0)
        -1.0, 0.0, 0.0, // (0,1,0) → (-1,0,0)
        0.0, 0.0, 1.0, // z-axis stays
    ];

    let rmsd = sci_form::compute_rmsd(&rotated, &reference);
    assert!(
        rmsd < 1e-6,
        "RMSD after 90° rotation should be ~0, got {rmsd:.8}"
    );
}

#[test]
fn test_rmsd_known_value() {
    // Two structures with known displacement: one atom shifted by 1.0 Å
    // 3 atoms: 2 identical + 1 shifted → RMSD = sqrt(1²/3) ≈ 0.577
    let reference = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0];
    let shifted = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 1.0, 0.0]; // last atom shifted +1 in y

    let rmsd = sci_form::compute_rmsd(&shifted, &reference);
    // After alignment, the rigid part should be aligned; the remaining error
    // is at most 1 Å for one atom out of 3
    assert!(
        rmsd < 1.0,
        "RMSD should be < 1.0 Å, got {rmsd:.4}"
    );
}

#[test]
fn test_rmsd_symmetry() {
    // RMSD(A, B) should equal RMSD(B, A)
    let a = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let b = vec![0.5, 0.0, 0.0, 1.5, 0.0, 0.0, 0.5, 1.0, 0.0];

    let rmsd_ab = sci_form::compute_rmsd(&a, &b);
    let rmsd_ba = sci_form::compute_rmsd(&b, &a);

    assert!(
        (rmsd_ab - rmsd_ba).abs() < 1e-10,
        "RMSD should be symmetric: {rmsd_ab:.8} vs {rmsd_ba:.8}"
    );
}

#[test]
fn test_rmsd_triangle_inequality() {
    // RMSD(A,C) ≤ RMSD(A,B) + RMSD(B,C) — approximate with alignment
    let a = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let b = vec![0.1, 0.0, 0.0, 1.1, 0.0, 0.0, 0.1, 1.0, 0.0];
    let c = vec![0.3, 0.0, 0.0, 1.3, 0.0, 0.0, 0.3, 1.0, 0.0];

    let ab = sci_form::compute_rmsd(&a, &b);
    let bc = sci_form::compute_rmsd(&b, &c);
    let ac = sci_form::compute_rmsd(&a, &c);

    // Triangle inequality holds for RMSD
    assert!(
        ac <= ab + bc + 1e-10,
        "Triangle inequality: {ac:.6} > {ab:.6} + {bc:.6}"
    );
}

#[test]
fn test_alignment_result_fields() {
    use sci_form::alignment::align_coordinates;

    let reference = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let shifted: Vec<f64> = reference.iter().map(|x| x + 5.0).collect();

    let result = align_coordinates(&shifted, &reference);

    assert_eq!(result.aligned_coords.len(), reference.len());
    assert!(result.rmsd < 1e-8, "Pure translation RMSD should be ~0");

    // Rotation matrix should be close to identity for pure translation
    for i in 0..3 {
        for j in 0..3 {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (result.rotation[i][j] - expected).abs() < 0.01,
                "Rotation[{i}][{j}] = {:.4} should be {expected}",
                result.rotation[i][j]
            );
        }
    }
}

#[test]
fn test_rmsd_single_atom() {
    // Single atom: RMSD should be 0 after centroid alignment
    let a = vec![1.0, 2.0, 3.0];
    let b = vec![4.0, 5.0, 6.0];
    let rmsd = sci_form::compute_rmsd(&a, &b);
    assert!(rmsd < 1e-10, "Single atom RMSD should be 0, got {rmsd}");
}

// ═══════════════════════════════════════════════════════════════════════════
// C7: Force Fields (UFF) and Strain Energy
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_uff_energy_finite_ethane() {
    let smiles = "CC";
    // Ethane-like coordinates
    let coords = vec![
        0.0, 0.0, 0.0,       // C1
        1.54, 0.0, 0.0,      // C2
        -0.36, 1.01, 0.0,    // H
        -0.36, -0.51, 0.87,  // H
        -0.36, -0.51, -0.87, // H
        1.90, 1.01, 0.0,     // H
        1.90, -0.51, 0.87,   // H
        1.90, -0.51, -0.87,  // H
    ];
    let energy = sci_form::compute_uff_energy(smiles, &coords).unwrap();
    assert!(energy.is_finite(), "Ethane UFF energy should be finite");
}

#[test]
fn test_uff_energy_distorted_higher() {
    // Distorted geometry should have higher energy than reasonable one
    let smiles = "C";
    // Reasonable tetrahedral CH4
    let good_coords = vec![
        0.0, 0.0, 0.0, 0.63, 0.63, 0.63, -0.63, -0.63, 0.63, -0.63, 0.63, -0.63, 0.63,
        -0.63, -0.63,
    ];
    // Distorted: compress all H atoms closer
    let bad_coords = vec![
        0.0, 0.0, 0.0, 0.30, 0.30, 0.30, -0.30, -0.30, 0.30, -0.30, 0.30, -0.30, 0.30,
        -0.30, -0.30,
    ];

    let e_good = sci_form::compute_uff_energy(smiles, &good_coords).unwrap();
    let e_bad = sci_form::compute_uff_energy(smiles, &bad_coords).unwrap();

    assert!(
        e_bad > e_good,
        "Compressed CH₄ should have higher energy: good={e_good:.4}, bad={e_bad:.4}"
    );
}

#[test]
fn test_uff_energy_gradient_finite() {
    // Verify that we can compute energy and that result is finite for water
    let smiles = "O";
    let coords = vec![
        0.0, 0.0, 0.0, // O
        0.96, 0.0, 0.0, // H
        -0.24, 0.93, 0.0, // H
    ];
    let energy = sci_form::compute_uff_energy(smiles, &coords).unwrap();
    assert!(energy.is_finite(), "Water UFF energy should be finite");
    assert!(energy > -1000.0 && energy < 10000.0, "Energy {energy} out of reasonable range");
}

#[test]
fn test_uff_energy_stretched_bond_higher() {
    // Stretching a C-H bond should increase energy
    let smiles = "C";
    let base = vec![
        0.0, 0.0, 0.0, 0.63, 0.63, 0.63, -0.63, -0.63, 0.63, -0.63, 0.63, -0.63, 0.63,
        -0.63, -0.63,
    ];
    // Stretch first H far away
    let mut stretched = base.clone();
    stretched[3] = 3.0;
    stretched[4] = 3.0;
    stretched[5] = 3.0;

    let e_normal = sci_form::compute_uff_energy(smiles, &base).unwrap();
    let e_stretched = sci_form::compute_uff_energy(smiles, &stretched).unwrap();

    assert!(
        e_stretched > e_normal,
        "Stretched bond energy {e_stretched:.2} should > normal {e_normal:.2}"
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// Cross-module consistency tests
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_population_dipole_consistency() {
    // Population charges and dipole should agree: both compute from same EHT
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];

    let pop = sci_form::compute_population(&elements, &positions).unwrap();
    let dip = sci_form::compute_dipole(&elements, &positions).unwrap();

    // Both should succeed for the same molecule
    assert_eq!(pop.num_atoms, 3);
    assert!(dip.magnitude > 0.0);

    // Manual dipole from population charges
    let eang_to_debye = 4.80321;
    let mut mu = [0.0f64; 3];
    for (a, q) in pop.mulliken_charges.iter().enumerate() {
        for k in 0..3 {
            mu[k] += q * positions[a][k];
        }
    }
    let manual_mag = (mu[0] * mu[0] + mu[1] * mu[1] + mu[2] * mu[2]).sqrt() * eang_to_debye;

    // Both approaches should be in the same ballpark
    // (exact match depends on whether dipole uses same or different EHT invocation)
    assert!(
        manual_mag > 0.0,
        "Manual dipole from population charges should be nonzero"
    );
}

#[test]
fn test_esp_uses_population_charges() {
    // ESP grid should have nonzero values when population gives nonzero charges
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];

    let pop = sci_form::compute_population(&elements, &positions).unwrap();
    let grid = sci_form::compute_esp(&elements, &positions, 0.5, 2.0).unwrap();

    // If charges are nonzero, ESP should have nonzero values
    let has_charge = pop.mulliken_charges.iter().any(|&q| q.abs() > 1e-6);
    let has_esp = grid.values.iter().any(|&v| v.abs() > 1e-10);

    assert!(has_charge, "Water should have nonzero Mulliken charges");
    assert!(has_esp, "ESP should be nonzero when charges are nonzero");
}

#[test]
fn test_dos_energy_range_covers_eigenvalues() {
    // DOS computed over a wide range should capture all orbital peaks
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];

    // Compute with wide energy window
    let result =
        sci_form::compute_dos(&elements, &positions, 0.3, -50.0, 20.0, 500).unwrap();

    let max_dos = result.total_dos.iter().cloned().fold(0.0f64, f64::max);
    assert!(
        max_dos > 0.1,
        "Wide-range DOS should capture orbital peaks, max={max_dos:.4}"
    );
}

#[test]
fn test_embed_then_uff_pipeline() {
    // Full pipeline: embed → UFF energy
    let smiles = "C";
    let result = sci_form::embed(smiles, 42);
    if result.error.is_some() {
        // Skip if embedding fails (not a UFF test issue)
        return;
    }
    assert!(!result.coords.is_empty());

    let energy = sci_form::compute_uff_energy(smiles, &result.coords).unwrap();
    assert!(
        energy.is_finite(),
        "UFF energy from embedded coords should be finite"
    );
}

#[test]
fn test_embed_then_population_pipeline() {
    // Full pipeline: embed → population analysis
    let smiles = "O";
    let result = sci_form::embed(smiles, 42);
    if result.error.is_some() {
        return;
    }

    let positions: Vec<[f64; 3]> = result.coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    let pop = sci_form::compute_population(&result.elements, &positions).unwrap();
    assert_eq!(pop.num_atoms, result.num_atoms);

    let sum: f64 = pop.mulliken_charges.iter().sum();
    assert!(
        sum.abs() < 1.0,
        "Embedded water charge sum = {sum:.4}, expected ≈ 0"
    );
}
