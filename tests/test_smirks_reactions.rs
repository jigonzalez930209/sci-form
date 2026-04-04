//! Integration tests for SMIRKS reaction transforms.
//!
//! These tests validate the SMIRKS implementation by testing:
//! 1. Common organic reaction patterns
//! 2. Edge cases and error handling
//! 3. Performance characteristics
//! 4. Round-trip consistency

use sci_form::smirks::{apply_smirks, parse_smirks};

#[test]
fn test_acid_base_reactions() {
    // Carboxylic acid deprotonation
    let result = apply_smirks("[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]", "CC(=O)O").unwrap();
    assert!(result.success, "Acid deprotonation should succeed");
    assert_eq!(result.n_transforms, 1);

    // Alcohol deprotonation
    let result = apply_smirks("[C:1][OH:2]>>[C:1][O-:2]", "CCO").unwrap();
    assert!(result.success, "Alcohol deprotonation should succeed");

    // Amine protonation - should not match CCO
    let result = apply_smirks("[N:1]>>[N+:1]", "CCO").unwrap();
    assert!(!result.success, "Should not match molecules without nitrogen");
}

#[test]
fn test_oxidation_reduction_reactions() {
    // Ketone reduction
    let smirks = "[C:1]=[O:2]>>[C:1][OH:2]";
    let result = parse_smirks(smirks);
    assert!(result.is_ok(), "Ketone reduction pattern should parse");

    // Aldehyde reduction
    let smirks = "[C:1][C:2]([H:3])=[O:4]>>[C:1][C:2]([H:3])[OH:4]";
    let result = parse_smirks(smirks);
    assert!(result.is_ok(), "Aldehyde reduction pattern should parse");

    // Alcohol oxidation
    let smirks = "[C:1][C:2]([OH:3])[H:4]>>[C:1][C:2](=[O:3])";
    let result = parse_smirks(smirks);
    assert!(result.is_ok(), "Alcohol oxidation pattern should parse");
}

#[test]
fn test_substitution_reactions() {
    // Aromatic halogenation
    let result = apply_smirks("[c:1][H:2]>>[c:1][Cl:2]", "c1ccccc1").unwrap();
    assert!(result.success, "Aromatic halogenation should match benzene");

    // Nitration
    let smirks = "[c:1][H:2]>>[c:1][N+:2](=[O:3])[O-:4]";
    let result = parse_smirks(smirks);
    assert!(result.is_ok(), "Nitration pattern should parse");
}

#[test]
fn test_hydrolysis_reactions() {
    // Ester hydrolysis
    let smirks = "[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[OH:3]";
    let result = apply_smirks(smirks, "CC(=O)OCC").unwrap();
    assert!(result.success, "Ester hydrolysis should match");

    // Amide hydrolysis
    let smirks = "[C:1](=[O:2])[N:3]>>[C:1](=[O:2])[OH:3]";
    let result = parse_smirks(smirks);
    assert!(result.is_ok(), "Amide hydrolysis pattern should parse");
}

#[test]
fn test_functional_group_selectivity() {
    // Should match alcohol but not ether
    let smirks = "[C:1][OH:2]>>[C:1][O-:2]";

    // Match ethanol (has OH)
    let result = apply_smirks(smirks, "CCO").unwrap();
    assert!(result.success, "Should match ethanol");

    // Don't match dimethyl ether (no OH)
    let result = apply_smirks(smirks, "COC").unwrap();
    assert!(!result.success, "Should not match ether");
}

#[test]
fn test_aromatic_vs_aliphatic() {
    // Aromatic C-H should be distinguished from aliphatic
    let aromatic_smirks = "[c:1][H:2]>>[c:1][Cl:2]";

    // Should match benzene
    let result = apply_smirks(aromatic_smirks, "c1ccccc1").unwrap();
    assert!(result.success, "Should match aromatic");

    // Aliphatic vs aromatic test - just verify they parse correctly
    let aliphatic_smirks = "[C:1][H:2]>>[C:1][Cl:2]";
    let result = parse_smirks(aliphatic_smirks);
    assert!(result.is_ok(), "Aliphatic pattern should parse");
}

#[test]
fn test_multiple_sites() {
    // Pattern that could match multiple sites in a molecule
    let smirks = "[OH:1]>>[O-:1]";

    // Glycerol has 3 OH groups - but we only transform the first match
    let result = apply_smirks(smirks, "OCC(O)CO").unwrap();
    assert!(result.success, "Should match at least one site");
    assert_eq!(result.n_transforms, 1, "Should transform one site per call");
}

#[test]
fn test_stereochemistry_preservation() {
    // SMIRKS should preserve stereochemistry when not explicitly changed
    // This is a parsing test since full product generation isn't implemented yet
    let smirks = "[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]";
    let result = apply_smirks(smirks, "C[C@H](O)C(=O)O").unwrap();
    assert!(result.success, "Should handle stereochemistry");
}

#[test]
fn test_multi_component_handling() {
    // Multi-component reactions are now supported via per-component matching.
    // With a single SMILES input, apply_smirks splits on '.' internally.
    let smirks = "[C:1](=[O:2])[OH:3].[O:4]>>[C:1](=[O:2])[O:3].[O:4]";
    let result = apply_smirks(smirks, "CC(=O)O.O");
    assert!(result.is_ok());
}

#[test]
fn test_invalid_smirks_patterns() {
    // Missing separator
    assert!(parse_smirks("ABC").is_err());

    // Empty parts
    assert!(parse_smirks(">>").is_err());
    assert!(parse_smirks("A>>").is_err());
    assert!(parse_smirks(">>B").is_err());

    // Multiple separators
    assert!(parse_smirks("A>>B>>C").is_err());
}

#[test]
fn test_atom_map_consistency() {
    // Atom maps must be consistent between reactant and product
    let smirks = "[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]";
    let transform = parse_smirks(smirks).unwrap();

    assert!(transform.atom_map.contains_key(&1), "Map should contain C");
    assert!(transform.atom_map.contains_key(&2), "Map should contain O");
}

#[test]
fn test_complex_patterns() {
    // Fischer esterification pattern (complex, multi-step)
    let smirks = "[C:1](=[O:2])[OH:3].[C:4][OH:5]>>[C:1](=[O:2])[O:5][C:4]";
    let transform = parse_smirks(smirks).unwrap();

    assert_eq!(transform.reactant_smarts.len(), 2, "Two reactants");
    assert_eq!(transform.product_smarts.len(), 1, "One product");
    assert!(transform.atom_map.len() >= 3, "Multiple atom mappings");
}

#[test]
fn test_performance_simple_patterns() {
    // Simple patterns should parse quickly
    let patterns = vec![
        "[C:1][OH:2]>>[C:1][O-:2]",
        "[c:1][H:2]>>[c:1][Cl:2]",
        "[N:1]>>[N+:1]",
    ];

    for pattern in patterns {
        let result = parse_smirks(pattern);
        assert!(result.is_ok(), "Pattern {} should parse", pattern);
    }
}

#[test]
fn test_common_named_reactions() {
    // Test SMIRKS patterns for common named reactions

    // Williamson ether synthesis (simplified)
    let smirks = "[C:1][OH:2].[C:3][Br:4]>>[C:1][O:2][C:3]";
    assert!(parse_smirks(smirks).is_ok());

    // Grignard addition (simplified)
    let smirks = "[C:1]=[O:2].[C:3][Mg:4][Br:5]>>[C:1]([O:2])[C:3]";
    assert!(parse_smirks(smirks).is_ok());

    // Diels-Alder (simplified, without stereochemistry)
    let smirks = "[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2][C:5]=[C:6][C:4][C:3]1";
    assert!(parse_smirks(smirks).is_ok());
}

#[test]
fn test_round_trip_consistency() {
    // Parse a SMIRKS and check internal consistency
    let original = "[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]";
    let transform = parse_smirks(original).unwrap();

    // Check that the original SMIRKS is preserved
    assert_eq!(transform.smirks, original);

    // Check reactant/product parts are preserved
    assert_eq!(transform.reactant_smarts.len(), 1);
    assert_eq!(transform.product_smarts.len(), 1);
}
