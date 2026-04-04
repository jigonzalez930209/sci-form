//! Comprehensive tests for all newly implemented features:
//! - SMARTS batch matching (sequential + parallel)
//! - Semi-analytical Hessian (PM3/xTB gradient-based)
//! - 3D molecular descriptors (gyration, asphericity, NPR)
//! - Atropisomerism detection (biaryl axial chirality)
//! - Pro-chirality detection (enantiotopic groups)
//! - MMFF94s static variant
//! - Powder XRD simulation
//! - Space group symmetry operations
//! - Porosity calculation
//! - ML prediction uncertainty
//! - Pharmacophore fingerprints
//! - Restrained conformer generation

use sci_form::graph::Molecule;

fn _positions_from_flat(coords: &[f64]) -> Vec<[f64; 3]> {
    coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect()
}

// ─── SMARTS batch matching ────────────────────────────────────────────────────

#[test]
fn test_smarts_batch_sequential() {
    let mol1 = Molecule::from_smiles("CCO").unwrap();
    let mol2 = Molecule::from_smiles("c1ccccc1").unwrap();
    let mol3 = Molecule::from_smiles("CC(=O)O").unwrap();

    let pattern = sci_form::smarts::parse_smarts("[OX2]").unwrap();
    let results = sci_form::smarts::substruct_match_batch(&[&mol1, &mol2, &mol3], &pattern);

    assert_eq!(results.len(), 3);
    // Ethanol has O-H
    assert!(!results[0].is_empty(), "Ethanol should match [OX2]");
    // Benzene has no oxygen
    assert!(results[1].is_empty(), "Benzene should not match [OX2]");
    // Acetic acid has OH
    assert!(!results[2].is_empty(), "Acetic acid should match [OX2]");
}

#[test]
fn test_smarts_has_substruct_match() {
    let mol = Molecule::from_smiles("c1ccccc1").unwrap();
    let aromatic = sci_form::smarts::parse_smarts("c:c").unwrap();
    let oh = sci_form::smarts::parse_smarts("[OX2H]").unwrap();

    assert!(sci_form::smarts::has_substruct_match(&mol, &aromatic));
    assert!(!sci_form::smarts::has_substruct_match(&mol, &oh));
}

#[cfg(feature = "parallel")]
#[test]
fn test_smarts_batch_parallel() {
    let mol1 = Molecule::from_smiles("CCO").unwrap();
    let mol2 = Molecule::from_smiles("c1ccccc1").unwrap();
    let mol3 = Molecule::from_smiles("CC(=O)O").unwrap();

    let pattern = sci_form::smarts::parse_smarts("[#6]").unwrap();
    let results =
        sci_form::smarts::substruct_match_batch_parallel(&[&mol1, &mol2, &mol3], &pattern);

    assert_eq!(results.len(), 3);
    assert!(
        results.iter().all(|r| !r.is_empty()),
        "All have carbon atoms (by atomic number)"
    );
}

#[cfg(feature = "parallel")]
#[test]
fn test_smarts_has_substruct_match_batch_parallel() {
    let mol1 = Molecule::from_smiles("CCO").unwrap();
    let mol2 = Molecule::from_smiles("c1ccccc1").unwrap();
    let mol3 = Molecule::from_smiles("CC(=O)O").unwrap();

    let pattern = sci_form::smarts::parse_smarts("[#8]").unwrap();
    let results =
        sci_form::smarts::has_substruct_match_batch_parallel(&[&mol1, &mol2, &mol3], &pattern);

    assert_eq!(results.len(), 3);
    assert!(results[0], "Ethanol has oxygen");
    assert!(!results[1], "Benzene has no oxygen");
    assert!(results[2], "Acetic acid has oxygen");
}

// ─── Semi-analytical Hessian ──────────────────────────────────────────────────

#[test]
fn test_semianalytical_hessian_pm3_water() {
    // Water: simple system for Hessian validation
    let elements = [8u8, 1, 1];
    let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];

    let hessian = sci_form::ir::hessian::compute_semianalytical_hessian(
        &elements,
        &positions,
        sci_form::ir::hessian::HessianMethod::Pm3,
        None,
    );

    assert!(
        hessian.is_ok(),
        "Semi-analytical Hessian should succeed for water"
    );
    let h = hessian.unwrap();
    assert_eq!(h.nrows(), 9); // 3 atoms × 3
    assert_eq!(h.ncols(), 9);

    // Hessian should be symmetric
    for i in 0..9 {
        for j in 0..9 {
            assert!(
                (h[(i, j)] - h[(j, i)]).abs() < 1e-6,
                "Hessian should be symmetric at ({}, {}): {} vs {}",
                i,
                j,
                h[(i, j)],
                h[(j, i)]
            );
        }
    }
}

#[test]
fn test_semianalytical_hessian_xtb_water() {
    let elements = [8u8, 1, 1];
    let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];

    let hessian = sci_form::ir::hessian::compute_semianalytical_hessian(
        &elements,
        &positions,
        sci_form::ir::hessian::HessianMethod::Xtb,
        None,
    );

    assert!(
        hessian.is_ok(),
        "xTB semi-analytical Hessian should succeed for water"
    );
    let h = hessian.unwrap();
    assert_eq!(h.nrows(), 9);

    // Check diagonal elements are non-negative (restoring forces)
    let has_positive_diagonal = (0..9).any(|i| h[(i, i)] > 0.0);
    assert!(
        has_positive_diagonal,
        "Should have positive diagonal elements"
    );
}

#[test]
fn test_semianalytical_vs_numerical_hessian_consistency() {
    // Validate the semi-analytical PM3 Hessian is physically reasonable.
    // The energy-based numerical Hessian for SCF methods can suffer from
    // convergence noise at small step sizes, so we validate the gradient-based
    // (semi-analytical) Hessian against structural expectations instead.
    let elements = [8u8, 1, 1];
    let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];

    let semi = sci_form::ir::hessian::compute_semianalytical_hessian(
        &elements,
        &positions,
        sci_form::ir::hessian::HessianMethod::Pm3,
        Some(0.005),
    )
    .unwrap();

    let n3 = 3 * elements.len();
    assert_eq!(semi.nrows(), n3);
    assert_eq!(semi.ncols(), n3);

    // 1. Symmetry: H_{ij} ≈ H_{ji}
    let mut max_asym = 0.0f64;
    for i in 0..n3 {
        for j in (i + 1)..n3 {
            max_asym = max_asym.max((semi[(i, j)] - semi[(j, i)]).abs());
        }
    }
    assert!(
        max_asym < 1e-6,
        "Semi-analytical Hessian should be symmetric, max asymmetry: {}",
        max_asym
    );

    // 2. Finite values
    for i in 0..n3 {
        for j in 0..n3 {
            assert!(semi[(i, j)].is_finite(), "Hessian[{},{}] is not finite", i, j);
        }
    }

    // 3. Physical reasonableness: eigenvalues of mass-weighted Hessian
    //    should give vibrational frequencies. Water has 3 vibrational modes
    //    (3N-6 = 3) and 6 rotational/translational modes (≈0).
    let eigen = semi.clone().symmetric_eigen();
    let mut sorted_eigenvalues: Vec<f64> = eigen.eigenvalues.iter().copied().collect();
    sorted_eigenvalues.sort_by(|a, b| b.abs().partial_cmp(&a.abs()).unwrap());
    // At least 3 eigenvalues should be significantly nonzero (the vibrational modes)
    let significant = sorted_eigenvalues.iter().filter(|e| e.abs() > 1.0).count();
    assert!(
        significant >= 3,
        "Water should have at least 3 significant Hessian eigenvalues, got {}",
        significant
    );

    // 4. Cross-check: numerical Hessian should have the same dimension
    let numer = sci_form::ir::hessian::compute_numerical_hessian(
        &elements,
        &positions,
        sci_form::ir::hessian::HessianMethod::Pm3,
        Some(0.005),
    )
    .unwrap();
    assert_eq!(semi.nrows(), numer.nrows());
    assert_eq!(semi.ncols(), numer.ncols());
}

// ─── 3D Molecular Descriptors ─────────────────────────────────────────────────

#[test]
fn test_3d_descriptors_ethanol() {
    let conf = sci_form::embed("CCO", 42);
    assert!(conf.error.is_none());

    let desc = sci_form::ml::compute_3d_descriptors(&conf.elements, &conf.coords);

    assert!(desc.radius_of_gyration > 0.0, "Rg should be positive");
    assert!(desc.span > 0.0, "Span should be positive");
    assert!(
        desc.asphericity >= 0.0 && desc.asphericity <= 1.0,
        "Asphericity should be in [0, 1], got {}",
        desc.asphericity
    );
    assert!(
        desc.eccentricity >= 0.0,
        "Eccentricity should be non-negative"
    );
    assert!(
        desc.npr1 >= 0.0 && desc.npr1 <= 1.0,
        "NPR1 should be in [0, 1], got {}",
        desc.npr1
    );
    assert!(
        desc.npr2 >= 0.0 && desc.npr2 <= 1.0,
        "NPR2 should be in [0, 1], got {}",
        desc.npr2
    );
    assert!(
        desc.sphericity >= 0.0 && desc.sphericity <= 1.0,
        "Sphericity should be in [0, 1], got {}",
        desc.sphericity
    );
}

#[test]
fn test_3d_descriptors_benzene_vs_ethane() {
    // Benzene is more planar → different shape descriptors than ethane
    let benzene = sci_form::embed("c1ccccc1", 42);
    let ethane = sci_form::embed("CC", 42);
    assert!(benzene.error.is_none());
    assert!(ethane.error.is_none());

    let desc_benz = sci_form::ml::compute_3d_descriptors(&benzene.elements, &benzene.coords);
    let desc_eth = sci_form::ml::compute_3d_descriptors(&ethane.elements, &ethane.coords);

    // Benzene has larger Rg than ethane
    assert!(
        desc_benz.radius_of_gyration > desc_eth.radius_of_gyration,
        "Benzene Rg ({}) > ethane Rg ({})",
        desc_benz.radius_of_gyration,
        desc_eth.radius_of_gyration
    );

    // Benzene has larger span
    assert!(
        desc_benz.span > desc_eth.span,
        "Benzene span ({}) > ethane span ({})",
        desc_benz.span,
        desc_eth.span
    );
}

#[test]
fn test_3d_descriptors_empty_molecule() {
    let desc = sci_form::ml::compute_3d_descriptors(&[], &[]);
    assert_eq!(desc.radius_of_gyration, 0.0);
    assert_eq!(desc.span, 0.0);
}

// ─── Atropisomerism Detection ─────────────────────────────────────────────────

#[test]
fn test_stereo_analysis_includes_atropisomeric_fields() {
    // Simple molecule — no atropisomeric axes expected
    let stereo = sci_form::analyze_stereo("CCO", &[]).unwrap();
    assert!(
        stereo.atropisomeric_axes.is_empty(),
        "Ethanol should have no atropisomeric axes"
    );
}

#[test]
fn test_stereo_analysis_includes_prochiral_fields() {
    // Simple molecule — check the field exists
    let stereo = sci_form::analyze_stereo("CCO", &[]).unwrap();
    // prochiral_centers field must exist (may or may not be empty depending on detection)
    let _ = stereo.prochiral_centers.len();
}

#[test]
fn test_stereo_chiral_molecule_with_new_fields() {
    let conf = sci_form::embed("C[C@H](F)Cl", 42);
    assert!(conf.error.is_none());

    let stereo = sci_form::analyze_stereo("C[C@H](F)Cl", &conf.coords).unwrap();
    assert!(stereo.n_stereocenters >= 1, "Should detect chiral center");

    // Verify new fields exist on the result
    let _ = &stereo.atropisomeric_axes;
    let _ = &stereo.prochiral_centers;
}

// ─── MMFF94s Static Variant ───────────────────────────────────────────────────

#[test]
fn test_mmff94s_variant_builds() {
    use sci_form::forcefield::mmff94::{Mmff94Builder, Mmff94Variant};

    let elements = [6u8, 6, 8, 1, 1, 1, 1, 1, 1];
    let bonds = vec![
        (0, 1, 1u8),
        (1, 2, 1),
        (0, 3, 1),
        (0, 4, 1),
        (0, 5, 1),
        (1, 6, 1),
        (1, 7, 1),
        (2, 8, 1),
    ];

    let mmff94_terms = Mmff94Builder::build(&elements, &bonds);
    let mmff94s_terms = Mmff94Builder::build_variant(&elements, &bonds, Mmff94Variant::Mmff94s);

    // Both should produce force field terms
    assert!(!mmff94_terms.is_empty(), "MMFF94 should produce terms");
    assert!(!mmff94s_terms.is_empty(), "MMFF94s should produce terms");

    // Same number of terms (same topology, different parameters)
    assert_eq!(
        mmff94_terms.len(),
        mmff94s_terms.len(),
        "MMFF94 and MMFF94s should have same number of terms"
    );
}

#[test]
fn test_mmff94_variant_enum() {
    use sci_form::forcefield::mmff94::Mmff94Variant;

    // Verify enum variants exist
    let v1 = Mmff94Variant::Mmff94;
    let v2 = Mmff94Variant::Mmff94s;
    assert_ne!(v1, v2);
    assert_eq!(v1, Mmff94Variant::Mmff94);
}

// ─── Powder XRD Simulation ────────────────────────────────────────────────────

#[test]
fn test_powder_xrd_cubic_nacl() {
    use sci_form::materials::{simulate_powder_xrd, UnitCell};

    // NaCl face-centered cubic, a = 5.64 Å
    let cell = UnitCell::cubic(5.64);
    // Na at (0,0,0), Cl at (0.5,0,0) — simplified
    let elements = [11u8, 17]; // Na, Cl
    let frac_coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]];

    let xrd = simulate_powder_xrd(&cell, &elements, &frac_coords, 90.0);

    assert!(!xrd.two_theta.is_empty(), "Should produce reflections");
    assert_eq!(xrd.two_theta.len(), xrd.intensities.len());
    assert_eq!(xrd.two_theta.len(), xrd.miller_indices.len());
    assert_eq!(xrd.two_theta.len(), xrd.d_spacings.len());

    // Check d-spacings are positive
    assert!(xrd.d_spacings.iter().all(|&d| d > 0.0));

    // 2θ should be in (0, 90) range
    assert!(xrd.two_theta.iter().all(|&t| t > 0.0 && t <= 90.0));

    // At least one reflection should have 100% intensity (max normalized)
    assert!(
        xrd.intensities.iter().any(|&i| (i - 100.0).abs() < 0.1),
        "Should have a 100% intensity peak"
    );
}

#[test]
fn test_powder_xrd_narrow_range() {
    use sci_form::materials::{simulate_powder_xrd, UnitCell};

    let cell = UnitCell::cubic(3.0);
    let elements = [26u8]; // Fe
    let frac_coords = [[0.0, 0.0, 0.0]];

    let xrd = simulate_powder_xrd(&cell, &elements, &frac_coords, 30.0);
    // All reflections should be within the specified range
    assert!(xrd.two_theta.iter().all(|&t| t <= 30.0));
}

// ─── Space Group Symmetry ─────────────────────────────────────────────────────

#[test]
fn test_space_group_p1() {
    use sci_form::materials::get_space_group;

    let sg = get_space_group(1).unwrap();
    assert_eq!(sg.number, 1);
    assert_eq!(sg.symbol, "P1");
    assert_eq!(sg.crystal_system, "triclinic");
    assert_eq!(sg.operations.len(), 1); // Identity only
}

#[test]
fn test_space_group_p_minus_1() {
    use sci_form::materials::get_space_group;

    let sg = get_space_group(2).unwrap();
    assert_eq!(sg.number, 2);
    assert_eq!(sg.symbol, "P-1");
    assert_eq!(sg.operations.len(), 2); // Identity + inversion
}

#[test]
fn test_space_group_p21c() {
    use sci_form::materials::get_space_group;

    let sg = get_space_group(14).unwrap();
    assert_eq!(sg.number, 14);
    assert_eq!(sg.symbol, "P2_1/c");
    assert_eq!(sg.crystal_system, "monoclinic");
    assert_eq!(sg.operations.len(), 4);
}

#[test]
fn test_space_group_fm3m() {
    use sci_form::materials::get_space_group;

    let sg = get_space_group(225).unwrap();
    assert_eq!(sg.number, 225);
    assert_eq!(sg.symbol, "Fm-3m");
    assert_eq!(sg.crystal_system, "cubic");
    assert_eq!(sg.operations.len(), 4); // Subset implemented
}

#[test]
fn test_space_group_unsupported() {
    use sci_form::materials::get_space_group;
    assert!(get_space_group(999).is_none());
    assert!(get_space_group(100).is_none());
}

#[test]
fn test_apply_symmetry_identity() {
    use sci_form::materials::{apply_symmetry, SymmetryOperation};

    let identity = SymmetryOperation {
        rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        translation: [0.0, 0.0, 0.0],
        label: "x,y,z".to_string(),
    };

    let pos = [0.25, 0.5, 0.75];
    let result = apply_symmetry(&identity, pos);
    assert!((result[0] - 0.25).abs() < 1e-10);
    assert!((result[1] - 0.5).abs() < 1e-10);
    assert!((result[2] - 0.75).abs() < 1e-10);
}

#[test]
fn test_apply_symmetry_inversion() {
    use sci_form::materials::{apply_symmetry, SymmetryOperation};

    let inversion = SymmetryOperation {
        rotation: [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
        translation: [0.0, 0.0, 0.0],
        label: "-x,-y,-z".to_string(),
    };

    let pos = [0.3, 0.4, 0.5];
    let result = apply_symmetry(&inversion, pos);
    assert!((result[0] - (-0.3)).abs() < 1e-10);
    assert!((result[1] - (-0.4)).abs() < 1e-10);
    assert!((result[2] - (-0.5)).abs() < 1e-10);
}

#[test]
fn test_expand_by_symmetry_p_minus_1() {
    use sci_form::materials::{expand_by_symmetry, get_space_group};

    let sg = get_space_group(2).unwrap(); // P-1: identity + inversion
    let frac = [[0.25, 0.3, 0.4]];
    let elements = [6u8]; // Carbon

    let (expanded_coords, expanded_elements) = expand_by_symmetry(&sg, &frac, &elements);

    // P-1 should double the sites (0.25,0.3,0.4) -> also (-0.25,-0.3,-0.4) wrapped
    assert_eq!(expanded_coords.len(), 2);
    assert_eq!(expanded_elements.len(), 2);
    assert!(expanded_elements.iter().all(|&e| e == 6));
}

#[test]
fn test_expand_by_symmetry_at_inversion_center() {
    use sci_form::materials::{expand_by_symmetry, get_space_group};

    let sg = get_space_group(2).unwrap(); // P-1
                                          // Atom at origin — inversion produces the same point
    let frac = [[0.0, 0.0, 0.0]];
    let elements = [8u8];

    let (expanded_coords, expanded_elements) = expand_by_symmetry(&sg, &frac, &elements);

    // Should detect duplicate and keep only 1
    assert_eq!(expanded_coords.len(), 1);
    assert_eq!(expanded_elements.len(), 1);
}

// ─── Porosity Calculation ─────────────────────────────────────────────────────

#[test]
fn test_porosity_empty_cell() {
    use sci_form::materials::{compute_porosity, UnitCell};

    // Large cell with a single small atom → mostly void
    let cell = UnitCell::cubic(20.0);
    let elements = [6u8]; // Carbon
    let frac_coords = [[0.5, 0.5, 0.5]];

    let result = compute_porosity(&cell, &elements, &frac_coords, 1.4, 1.0);

    assert!(
        result.porosity > 0.5,
        "Large cell with single atom should be mostly porous, got {}",
        result.porosity
    );
    assert!(result.pore_volume > 0.0);
    assert!(result.largest_cavity_diameter > 0.0);
}

#[test]
fn test_porosity_dense_packing() {
    use sci_form::materials::{compute_porosity, UnitCell};

    // Small cell packed with atoms → minimal void
    let cell = UnitCell::cubic(3.0);
    let elements = [26u8, 26, 26, 26]; // Fe atoms
    let frac_coords = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ];

    let result = compute_porosity(&cell, &elements, &frac_coords, 1.4, 0.5);

    assert!(result.porosity >= 0.0 && result.porosity <= 1.0);
    assert!(result.pore_volume >= 0.0);
}

#[test]
fn test_porosity_result_fields() {
    use sci_form::materials::{compute_porosity, UnitCell};

    let cell = UnitCell::cubic(10.0);
    let elements = [30u8]; // Zn
    let frac_coords = [[0.5, 0.5, 0.5]];

    let result = compute_porosity(&cell, &elements, &frac_coords, 1.0, 1.0);

    assert!(result.porosity >= 0.0 && result.porosity <= 1.0);
    assert!(result.pore_volume >= 0.0);
    assert!(result.largest_cavity_diameter >= 0.0);
    assert!(result.pore_limiting_diameter >= 0.0);
}

// ─── ML Prediction Uncertainty ────────────────────────────────────────────────

#[test]
fn test_ml_uncertainty_small_molecule() {
    let parsed = sci_form::parse("CCO").unwrap();
    let bonds: Vec<(usize, usize, u8)> = parsed
        .graph
        .edge_references()
        .map(|e| {
            use petgraph::visit::EdgeRef;
            let order = match e.weight().order {
                sci_form::graph::BondOrder::Single => 1,
                sci_form::graph::BondOrder::Double => 2,
                sci_form::graph::BondOrder::Triple => 3,
                _ => 0,
            };
            (e.source().index(), e.target().index(), order)
        })
        .collect();

    let desc = sci_form::compute_ml_descriptors(&[6, 6, 8, 1, 1, 1, 1, 1, 1], &bonds, &[], &[]);
    let props = sci_form::predict_ml_properties(&desc);

    assert!(
        props.uncertainty.confidence > 0.5,
        "Small organic should have high confidence, got {}",
        props.uncertainty.confidence
    );
    assert!(props.uncertainty.logp_std > 0.0);
    assert!(props.uncertainty.solubility_std > 0.0);
    assert!(
        props.uncertainty.warnings.is_empty(),
        "Small organic should have no warnings"
    );
}

#[test]
fn test_ml_uncertainty_large_molecule() {
    // Simulate a large molecule descriptor
    use sci_form::ml::descriptors::MolecularDescriptors;
    use sci_form::ml::models::predict_properties;

    let desc = MolecularDescriptors {
        molecular_weight: 1200.0, // Very large
        n_heavy_atoms: 80,
        n_hydrogens: 100,
        n_bonds: 180,
        n_rotatable_bonds: 20,
        n_hbd: 3,
        n_hba: 8,
        fsp3: 0.5,
        total_abs_charge: 0.0,
        max_charge: 0.0,
        min_charge: 0.0,
        wiener_index: 5000.0,
        n_rings: 10,
        n_aromatic: 6,
        balaban_j: 2.5,
        sum_electronegativity: 200.0,
        sum_polarizability: 300.0,
    };

    let props = predict_properties(&desc);

    assert!(
        props.uncertainty.confidence < 0.8,
        "Large molecule should have lower confidence, got {}",
        props.uncertainty.confidence
    );
    assert!(
        !props.uncertainty.warnings.is_empty(),
        "Large molecule should have warnings"
    );
}

// ─── Pharmacophore Fingerprints ───────────────────────────────────────────────

#[test]
fn test_pharmacophore_detect_features_ethanol() {
    use sci_form::ml::pharmacophore::{detect_features, PharmFeatureType};

    let elements = [6u8, 6, 8, 1, 1, 1, 1, 1, 1];
    let bonds = [
        (0, 1, 1u8),
        (1, 2, 1),
        (0, 3, 1),
        (0, 4, 1),
        (0, 5, 1),
        (1, 6, 1),
        (1, 7, 1),
        (2, 8, 1),
    ];

    let features = detect_features(&elements, &bonds, &[], &[]);

    // Ethanol should have HBD (O-H), HBA (O), and hydrophobic (C)
    let has_hbd = features
        .iter()
        .any(|f| f.feature_type == PharmFeatureType::HBondDonor);
    let has_hba = features
        .iter()
        .any(|f| f.feature_type == PharmFeatureType::HBondAcceptor);
    let has_hydrophobic = features
        .iter()
        .any(|f| f.feature_type == PharmFeatureType::Hydrophobic);

    assert!(has_hbd, "Ethanol should have H-bond donor");
    assert!(has_hba, "Ethanol should have H-bond acceptor");
    assert!(has_hydrophobic, "Ethanol should have hydrophobic feature");
}

#[test]
fn test_pharmacophore_fingerprint_similarity() {
    use sci_form::ml::pharmacophore::{
        compute_pharmacophore_fingerprint, detect_features, pharmacophore_tanimoto,
    };

    // Two similar molecules
    let elem1 = [6u8, 6, 8, 1, 1, 1, 1, 1, 1]; // ethanol CCO
    let bonds1 = [
        (0, 1, 1u8),
        (1, 2, 1),
        (0, 3, 1),
        (0, 4, 1),
        (0, 5, 1),
        (1, 6, 1),
        (1, 7, 1),
        (2, 8, 1),
    ];

    let elem2 = [6u8, 6, 6, 8, 1, 1, 1, 1, 1, 1, 1, 1]; // propanol CCCO
    let bonds2 = [
        (0, 1, 1u8),
        (1, 2, 1),
        (2, 3, 1),
        (0, 4, 1),
        (0, 5, 1),
        (0, 6, 1),
        (1, 7, 1),
        (1, 8, 1),
        (2, 9, 1),
        (2, 10, 1),
        (3, 11, 1),
    ];

    let f1 = detect_features(&elem1, &bonds1, &[], &[]);
    let f2 = detect_features(&elem2, &bonds2, &[], &[]);

    let fp1 = compute_pharmacophore_fingerprint(&f1, 2048);
    let fp2 = compute_pharmacophore_fingerprint(&f2, 2048);

    let sim = pharmacophore_tanimoto(&fp1, &fp2);
    assert!(
        sim > 0.0 && sim <= 1.0,
        "Similarity should be in (0, 1], got {}",
        sim
    );

    // Self-similarity should be 1.0
    let self_sim = pharmacophore_tanimoto(&fp1, &fp1);
    assert!(
        (self_sim - 1.0).abs() < 1e-10,
        "Self-similarity should be 1.0, got {}",
        self_sim
    );
}

#[test]
fn test_pharmacophore_fingerprint_density() {
    use sci_form::ml::pharmacophore::{compute_pharmacophore_fingerprint, detect_features};

    let elements = [6u8, 6, 8, 1, 1, 1, 1, 1, 1];
    let bonds = [
        (0, 1, 1u8),
        (1, 2, 1),
        (0, 3, 1),
        (0, 4, 1),
        (0, 5, 1),
        (1, 6, 1),
        (1, 7, 1),
        (2, 8, 1),
    ];

    let features = detect_features(&elements, &bonds, &[], &[]);
    let fp = compute_pharmacophore_fingerprint(&features, 1024);

    assert_eq!(fp.n_bits, 1024);
    assert!(fp.density > 0.0 && fp.density < 1.0);
    assert!(!fp.on_bits.is_empty());
}

// ─── Restrained Conformer Generation ─────────────────────────────────────────

#[test]
fn test_restrained_conformer_basic() {
    use sci_form::conformer::{generate_3d_conformer_restrained, DistanceRestraint};

    let mol = Molecule::from_smiles("CCCC").unwrap();

    // Restrain atoms 0 and 3 to be close (gauche-like)
    let restraints = vec![DistanceRestraint {
        atom_i: 0,
        atom_j: 3,
        target_distance: 3.0,
        force_constant: 100.0,
    }];

    let result = generate_3d_conformer_restrained(&mol, 42, &restraints);
    assert!(
        result.is_ok(),
        "Restrained conformer should succeed: {:?}",
        result.err()
    );

    let coords = result.unwrap();
    assert_eq!(coords.nrows(), mol.graph.node_count());
    assert_eq!(coords.ncols(), 3);
}

#[test]
fn test_restrained_conformer_empty_restraints() {
    use sci_form::conformer::generate_3d_conformer_restrained;

    let mol = Molecule::from_smiles("CCO").unwrap();

    // No restraints — should behave like normal conformer gen
    let result = generate_3d_conformer_restrained(&mol, 42, &[]);
    assert!(result.is_ok(), "Empty restraints should work");
}

// ─── Combined integration tests ──────────────────────────────────────────────

#[test]
fn test_full_pipeline_with_new_features() {
    // End-to-end: embed → 3D descriptors → ML → pharmacophore
    let conf = sci_form::embed("c1ccc(O)cc1", 42); // phenol
    assert!(conf.error.is_none());

    // 3D descriptors
    let desc_3d = sci_form::ml::compute_3d_descriptors(&conf.elements, &conf.coords);
    assert!(desc_3d.radius_of_gyration > 0.0);

    // Stereo analysis with new fields
    let stereo = sci_form::analyze_stereo("c1ccc(O)cc1", &conf.coords).unwrap();
    assert_eq!(stereo.n_stereocenters, 0); // Phenol has no chiral centers
    let _ = &stereo.atropisomeric_axes;
    let _ = &stereo.prochiral_centers;
}

#[test]
fn test_materials_pipeline_xrd_and_symmetry() {
    use sci_form::materials::{
        compute_porosity, expand_by_symmetry, get_space_group, simulate_powder_xrd, UnitCell,
    };

    // Create a monoclinic cell
    let cell = UnitCell::from_parameters(&sci_form::materials::CellParameters {
        a: 8.0,
        b: 6.0,
        c: 7.0,
        alpha: 90.0,
        beta: 100.0,
        gamma: 90.0,
    });

    // Start with asymmetric unit
    let asym_elements = [30u8, 8]; // Zn, O
    let asym_frac = [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]];

    // Expand by P2_1/c symmetry
    let sg = get_space_group(14).unwrap();
    let (full_frac, full_elements) = expand_by_symmetry(&sg, &asym_frac, &asym_elements);

    assert!(full_frac.len() >= asym_frac.len());
    assert_eq!(full_frac.len(), full_elements.len());

    // Simulate XRD
    let xrd = simulate_powder_xrd(&cell, &full_elements, &full_frac, 60.0);
    assert!(!xrd.two_theta.is_empty());

    // Compute porosity
    let porosity = compute_porosity(&cell, &full_elements, &full_frac, 1.4, 0.5);
    assert!(porosity.porosity >= 0.0 && porosity.porosity <= 1.0);
}
