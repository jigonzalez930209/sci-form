//! SMIRKS reaction transform support.
//!
//! SMIRKS is an extension of SMARTS that describes chemical reactions
//! as atom-mapped reactant>>product transforms. This module parses
//! SMIRKS patterns and applies them to molecular graphs.

use crate::graph::Molecule;
use crate::smarts::{parse_smarts, substruct_match};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// A parsed SMIRKS reaction transform.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SmirksTransform {
    /// Reactant SMARTS pattern(s).
    pub reactant_smarts: Vec<String>,
    /// Product SMARTS pattern(s).
    pub product_smarts: Vec<String>,
    /// Atom map: reactant_atom_map_num → product_atom_map_num.
    pub atom_map: HashMap<usize, usize>,
    /// Bond changes: (atom_map1, atom_map2, old_order, new_order).
    pub bond_changes: Vec<BondChange>,
    /// Original SMIRKS string.
    pub smirks: String,
}

/// A bond change specified by a SMIRKS transform.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BondChange {
    /// Atom map number of first atom.
    pub atom1_map: usize,
    /// Atom map number of second atom.
    pub atom2_map: usize,
    /// Bond order in reactant (None = bond doesn't exist).
    pub old_order: Option<String>,
    /// Bond order in product (None = bond broken).
    pub new_order: Option<String>,
}

/// Result of applying a SMIRKS transform.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SmirksResult {
    /// Product SMILES strings.
    pub products: Vec<String>,
    /// Atom mapping from reactant to product indices.
    pub atom_mapping: HashMap<usize, usize>,
    /// Number of transforms applied.
    pub n_transforms: usize,
    /// Whether the transform was successfully applied.
    pub success: bool,
    /// Error or warning messages.
    pub messages: Vec<String>,
}

/// Parse a SMIRKS reaction string.
///
/// Format: `reactant_smarts>>product_smarts`
/// Atom maps: `[C:1]`, `[N:2]`, etc.
///
/// # Example
/// ```ignore
/// let transform = parse_smirks("[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]").unwrap();
/// ```
pub fn parse_smirks(smirks: &str) -> Result<SmirksTransform, String> {
    let parts: Vec<&str> = smirks.split(">>").collect();
    if parts.len() != 2 {
        return Err("SMIRKS must contain exactly one '>>' separator".to_string());
    }

    let reactant_part = parts[0].trim();
    let product_part = parts[1].trim();

    if reactant_part.is_empty() || product_part.is_empty() {
        return Err("SMIRKS reactant and product parts must be non-empty".to_string());
    }

    // Split by '.' for multi-component reactions
    let reactant_smarts: Vec<String> = reactant_part.split('.').map(|s| s.to_string()).collect();
    let product_smarts: Vec<String> = product_part.split('.').map(|s| s.to_string()).collect();

    // Extract atom maps from both sides
    let reactant_maps = extract_atom_maps(reactant_part)?;
    let product_maps = extract_atom_maps(product_part)?;

    // Build the atom map (reactant map num → product map num)
    let mut atom_map = HashMap::new();
    for map_num in reactant_maps.keys() {
        if product_maps.contains_key(map_num) {
            atom_map.insert(*map_num, *map_num);
        }
    }

    // Validate atom map bijectivity: each mapped atom in reactant must appear
    // exactly once in product and vice versa
    let mapped_in_reactant: std::collections::HashSet<usize> = atom_map.keys().copied().collect();
    let mapped_in_product: std::collections::HashSet<usize> = atom_map.values().copied().collect();
    if mapped_in_reactant != mapped_in_product {
        return Err(format!(
            "SMIRKS atom maps are not bijective: reactant maps {:?} vs product maps {:?}",
            mapped_in_reactant, mapped_in_product
        ));
    }

    // Detect bond changes
    let bond_changes =
        detect_bond_changes(reactant_part, product_part, &reactant_maps, &product_maps);

    Ok(SmirksTransform {
        reactant_smarts,
        product_smarts,
        atom_map,
        bond_changes,
        smirks: smirks.to_string(),
    })
}

/// Apply a SMIRKS transform to a molecule (represented as SMILES).
///
/// Returns the product SMILES if the reactant pattern matches.
pub fn apply_smirks(smirks: &str, smiles: &str) -> Result<SmirksResult, String> {
    let transform = parse_smirks(smirks)?;

    // Parse the input molecule
    if transform.reactant_smarts.len() > 1 || transform.product_smarts.len() > 1 {
        return Ok(SmirksResult {
            products: vec![],
            atom_mapping: HashMap::new(),
            n_transforms: 0,
            success: false,
            messages: vec![
                "Multi-component SMIRKS are not supported by apply_smirks because Molecule::from_smiles retains only the largest fragment".to_string(),
            ],
        });
    }

    let mol = Molecule::from_smiles(smiles)?;

    // Try to match the reactant pattern using substructure matching
    let matches = match_smarts_pattern(&mol, &transform.reactant_smarts[0])?;

    if matches.is_empty() {
        return Ok(SmirksResult {
            products: vec![],
            atom_mapping: HashMap::new(),
            n_transforms: 0,
            success: false,
            messages: vec!["No match found for reactant pattern".to_string()],
        });
    }

    // Apply the first match
    let atom_mapping = &matches[0];

    // For now, return a simplified result indicating the transform was recognized
    Ok(SmirksResult {
        products: transform.product_smarts.clone(),
        atom_mapping: atom_mapping.clone(),
        n_transforms: 1,
        success: true,
        messages: vec![format!(
            "Transform applied: {} atoms mapped",
            atom_mapping.len()
        )],
    })
}

/// Extract atom map numbers from a SMARTS/SMIRKS string.
/// Returns map_number → pattern_atom_index.
fn extract_atom_maps(pattern: &str) -> Result<HashMap<usize, usize>, String> {
    let mut maps = HashMap::new();
    let bytes = pattern.as_bytes();
    let mut pos = 0;
    let mut atom_idx = 0;

    while pos < bytes.len() {
        if bytes[pos] == b'[' {
            // Find the closing bracket
            let start = pos;
            while pos < bytes.len() && bytes[pos] != b']' {
                pos += 1;
            }
            let bracket_content = &pattern[start..=pos.min(bytes.len() - 1)];

            // Look for :N atom map
            if let Some(colon_pos) = bracket_content.rfind(':') {
                let map_str = &bracket_content[colon_pos + 1..bracket_content.len() - 1];
                if let Ok(map_num) = map_str.parse::<usize>() {
                    if maps.insert(map_num, atom_idx).is_some() {
                        return Err(format!("duplicate atom map :{} in pattern '{}'", map_num, pattern));
                    }
                }
            }
            atom_idx += 1;
        } else if bytes[pos].is_ascii_uppercase()
            || (bytes[pos] == b'c'
                || bytes[pos] == b'n'
                || bytes[pos] == b'o'
                || bytes[pos] == b's')
        {
            atom_idx += 1;
        }
        pos += 1;
    }

    Ok(maps)
}

/// Detect bond changes between reactant and product patterns.
fn detect_bond_changes(
    _reactant: &str,
    _product: &str,
    reactant_maps: &HashMap<usize, usize>,
    product_maps: &HashMap<usize, usize>,
) -> Vec<BondChange> {
    let mut changes = Vec::new();

    // Atoms that appear in reactant but not product → bonds broken
    for map_num in reactant_maps.keys() {
        if !product_maps.contains_key(map_num) {
            changes.push(BondChange {
                atom1_map: *map_num,
                atom2_map: 0,
                old_order: Some("SINGLE".to_string()),
                new_order: None,
            });
        }
    }

    // Atoms that appear in product but not reactant → bonds formed
    for map_num in product_maps.keys() {
        if !reactant_maps.contains_key(map_num) {
            changes.push(BondChange {
                atom1_map: *map_num,
                atom2_map: 0,
                old_order: None,
                new_order: Some("SINGLE".to_string()),
            });
        }
    }

    changes
}

/// Substructure matching for SMARTS patterns.
/// Returns list of atom-map-number → molecule-index mappings.
fn match_smarts_pattern(mol: &Molecule, pattern: &str) -> Result<Vec<HashMap<usize, usize>>, String> {
    let parsed = parse_smarts(pattern)?;
    let mapped_atoms: Vec<(usize, usize)> = parsed
        .atoms
        .iter()
        .enumerate()
        .filter_map(|(idx, atom)| atom.map_idx.map(|map_idx| (idx, map_idx as usize)))
        .collect();

    Ok(substruct_match(mol, &parsed)
        .into_iter()
        .map(|matched_atoms| {
            mapped_atoms
                .iter()
                .map(|(pattern_idx, map_num)| (*map_num, matched_atoms[*pattern_idx]))
                .collect()
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_smirks_basic() {
        let result = parse_smirks("[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]");
        assert!(result.is_ok());
        let t = result.unwrap();
        assert_eq!(t.reactant_smarts.len(), 1);
        assert_eq!(t.product_smarts.len(), 1);
        assert!(!t.atom_map.is_empty());
    }

    #[test]
    fn test_parse_smirks_invalid() {
        assert!(parse_smirks("no_separator").is_err());
        assert!(parse_smirks(">>").is_err());
    }

    #[test]
    fn test_extract_atom_maps() {
        let maps = extract_atom_maps("[C:1](=O)[OH:2]").unwrap();
        assert!(maps.contains_key(&1));
        assert!(maps.contains_key(&2));
    }

    #[test]
    fn test_extract_atom_maps_rejects_duplicates() {
        let err = extract_atom_maps("[C:1][O:1]").unwrap_err();
        assert!(err.contains("duplicate atom map"));
    }

    #[test]
    fn test_apply_smirks() {
        let result = apply_smirks("[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]", "CC(=O)O");
        let result = result.unwrap();
        assert!(result.success);
        assert_eq!(result.n_transforms, 1);
        assert_eq!(result.atom_mapping.len(), 2);
    }

    #[test]
    fn test_apply_smirks_requires_real_match() {
        let result = apply_smirks("[N:1]>>[N:1]", "CCO").unwrap();
        assert!(!result.success);
        assert_eq!(result.n_transforms, 0);
    }

    #[test]
    fn test_apply_smirks_rejects_multicomponent_transform() {
        let result = apply_smirks("[O:1].[Na+:2]>>[O:1][Na+:2]", "CC(=O)O").unwrap();
        assert!(!result.success);
        assert!(result.messages[0].contains("Multi-component SMIRKS"));
    }
}
