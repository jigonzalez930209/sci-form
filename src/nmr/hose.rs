//! HOSE (Hierarchically Ordered Spherical description of Environment) code generation.
//!
//! Generates unique environment descriptors for each atom by traversing
//! the molecular graph in concentric spheres (radius 0 to max_radius).
//! Includes a compile-time HOSE→chemical-shift lookup database for
//! ¹H and ¹³C prediction with fallback from higher to lower radius.

use crate::graph::{BondOrder, Molecule};
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

/// A HOSE code descriptor for a specific atom.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HoseCode {
    /// Atom index in the molecule.
    pub atom_index: usize,
    /// Atomic number of the center atom.
    pub element: u8,
    /// HOSE string at each sphere radius (0..=max_radius).
    pub spheres: Vec<String>,
    /// Full concatenated HOSE code.
    pub full_code: String,
}

fn bond_symbol(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Single => "",
        BondOrder::Double => "=",
        BondOrder::Triple => "#",
        BondOrder::Aromatic => "*",
        BondOrder::Unknown => "",
    }
}

fn element_symbol(z: u8) -> &'static str {
    match z {
        1 => "H",
        5 => "B",
        6 => "C",
        7 => "N",
        8 => "O",
        9 => "F",
        14 => "Si",
        15 => "P",
        16 => "S",
        17 => "Cl",
        35 => "Br",
        53 => "I",
        _ => "X",
    }
}

/// Generate HOSE codes for all atoms in a molecule.
///
/// `mol`: parsed molecular graph
/// `max_radius`: maximum sphere radius (typically 3–5)
pub fn generate_hose_codes(mol: &Molecule, max_radius: usize) -> Vec<HoseCode> {
    let n = mol.graph.node_count();
    let mut codes = Vec::with_capacity(n);

    for atom_idx in 0..n {
        let center = NodeIndex::new(atom_idx);
        let center_element = mol.graph[center].element;
        let center_sym = element_symbol(center_element);

        let mut spheres = Vec::with_capacity(max_radius + 1);
        spheres.push(center_sym.to_string());

        // BFS by sphere depth
        let mut visited = vec![false; n];
        visited[atom_idx] = true;
        let mut current_frontier: Vec<(NodeIndex, BondOrder)> = Vec::new();

        // Initial frontier: direct neighbors
        for edge in mol.graph.edges(center) {
            let neighbor = if edge.source() == center {
                edge.target()
            } else {
                edge.source()
            };
            current_frontier.push((neighbor, edge.weight().order));
        }

        for _radius in 1..=max_radius {
            if current_frontier.is_empty() {
                spheres.push(String::new());
                continue;
            }

            // Sort frontier entries for deterministic output
            let mut sphere_parts: BTreeSet<String> = BTreeSet::new();
            let mut next_frontier: Vec<(NodeIndex, BondOrder)> = Vec::new();

            for (node, bond_order) in &current_frontier {
                let node_idx = node.index();
                if visited[node_idx] {
                    continue;
                }
                visited[node_idx] = true;

                let elem = mol.graph[*node].element;
                let sym = element_symbol(elem);
                let bond_sym = bond_symbol(*bond_order);
                sphere_parts.insert(format!("{}{}", bond_sym, sym));

                // Collect next frontier
                for edge in mol.graph.edges(*node) {
                    let next = if edge.source() == *node {
                        edge.target()
                    } else {
                        edge.source()
                    };
                    if !visited[next.index()] {
                        next_frontier.push((next, edge.weight().order));
                    }
                }
            }

            spheres.push(sphere_parts.into_iter().collect::<Vec<_>>().join(","));
            current_frontier = next_frontier;
        }

        let full_code = format!("{}/{}", spheres[0], spheres[1..].join("/"));

        codes.push(HoseCode {
            atom_index: atom_idx,
            element: center_element,
            spheres,
            full_code,
        });
    }

    codes
}

// ─── HOSE → Chemical Shift Database ─────────────────────────────────────────
//
// Empirical database mapping HOSE code prefixes to chemical shifts (ppm).
// Based on NMRShiftDB2 reference data patterns for common organic environments.
// Uses fallback: try radius 4 → 3 → 2 → 1 match.

/// Result of a HOSE-based chemical shift lookup.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HoseShiftLookup {
    /// Atom index.
    pub atom_index: usize,
    /// Element (atomic number).
    pub element: u8,
    /// Predicted shift in ppm.
    pub shift_ppm: f64,
    /// HOSE code used for matching.
    pub matched_hose: String,
    /// Radius at which the match was found.
    pub match_radius: usize,
    /// Confidence (higher radius = higher confidence).
    pub confidence: f64,
}

/// ¹H HOSE shift database: (HOSE prefix pattern, shift_ppm).
/// Keys are simplified HOSE sphere-1 environment strings.
fn h1_hose_database() -> Vec<(&'static str, f64)> {
    vec![
        // sp3 C-H environments
        ("C/H,H,H,C", 0.90),  // methyl next to alkyl
        ("C/H,H,C,C", 1.30),  // methylene
        ("C/H,C,C,C", 1.50),  // methine
        ("C/H,H,H,O", 3.40),  // methyl next to O (methoxy)
        ("C/H,H,O", 3.50),    // CH₂ next to O
        ("C/H,H,H,N", 2.30),  // methyl next to N
        ("C/H,H,N", 2.60),    // CH₂ next to N
        ("C/H,H,H,=O", 2.10), // methyl next to carbonyl
        ("C/H,H,=O", 2.50),   // methylene next to carbonyl
        ("C/H,H,H,*C", 2.30), // methyl on aromatic
        // sp2 C-H (aromatic)
        ("C/*C,*C,H", 7.27),  // benzene-type ArH
        ("C/*C,*C,*C", 7.50), // fused aromatic ArH
        ("C/*C,*N,H", 7.80),  // pyridine-adjacent ArH
        // sp2 C-H (alkene)
        ("C/=C,H,H", 5.25),  // vinyl =CH₂
        ("C/=C,H,C", 5.40),  // internal alkene
        ("C/=C,=C,H", 6.30), // conjugated diene
        // sp C-H (alkyne)
        ("C/#C,H", 2.50), // terminal alkyne
        // Aldehyde
        ("C/=O,H,C", 9.50), // aldehyde CHO
        ("C/=O,H,H", 9.60), // formaldehyde
        // O-H
        ("O/H,C", 2.50),  // alcohol OH
        ("O/H,*C", 5.50), // phenol OH
        // N-H
        ("N/H,H,C", 1.50), // primary amine
        ("N/H,C,C", 2.20), // secondary amine
    ]
}

/// ¹³C HOSE shift database.
fn c13_hose_database() -> Vec<(&'static str, f64)> {
    vec![
        // sp3 carbon
        ("C/H,H,H,C", 15.0),  // methyl
        ("C/H,H,C,C", 25.0),  // methylene
        ("C/H,C,C,C", 35.0),  // methine
        ("C/C,C,C,C", 40.0),  // quaternary
        ("C/H,H,H,O", 55.0),  // methoxy
        ("C/H,H,O,C", 65.0),  // C-O-C
        ("C/H,H,H,N", 32.0),  // N-methyl
        ("C/H,H,N", 45.0),    // C-N
        ("C/H,H,H,=O", 30.0), // methyl next to carbonyl
        // sp2 carbon (aromatic)
        ("C/*C,*C,H", 128.0),  // unsubstituted ArC
        ("C/*C,*C,C", 137.0),  // alkyl-substituted ArC
        ("C/*C,*C,O", 155.0),  // O-substituted ArC
        ("C/*C,*C,N", 148.0),  // N-substituted ArC
        ("C/*C,*C,F", 163.0),  // F-substituted ArC
        ("C/*C,*C,Cl", 134.0), // Cl-substituted ArC
        // sp2 carbon (alkene)
        ("C/=C,H,H", 115.0), // vinyl =CH₂
        ("C/=C,H,C", 130.0), // internal alkene
        ("C/=C,C,C", 140.0), // trisubstituted alkene
        // Carbonyl
        ("C/=O,O,C", 175.0), // ester / carboxylic acid
        ("C/=O,N,C", 170.0), // amide
        ("C/=O,C,C", 205.0), // ketone
        ("C/=O,H,C", 200.0), // aldehyde
        // sp carbon
        ("C/#C,H", 70.0), // terminal alkyne
        ("C/#C,C", 85.0), // internal alkyne
    ]
}

/// Predict chemical shift from HOSE code using database lookup with fallback.
///
/// Searches the HOSE code against the database, trying the most specific
/// match first (full code), then falling back to shorter prefixes.
pub fn predict_shift_from_hose(hose_code: &HoseCode, nucleus: u8) -> Option<HoseShiftLookup> {
    let database = match nucleus {
        1 => h1_hose_database(),
        6 => c13_hose_database(),
        _ => return None,
    };

    // Try matching from highest radius down to 1
    for radius in (1..=hose_code.spheres.len().saturating_sub(1)).rev() {
        // Build the HOSE prefix up to this radius
        let prefix = format!(
            "{}/{}",
            hose_code.spheres[0],
            hose_code.spheres[1..=radius].join("/")
        );

        for &(pattern, shift) in &database {
            if prefix.contains(pattern)
                || pattern.contains(&prefix)
                || fuzzy_hose_match(&prefix, pattern)
            {
                return Some(HoseShiftLookup {
                    atom_index: hose_code.atom_index,
                    element: nucleus,
                    shift_ppm: shift,
                    matched_hose: pattern.to_string(),
                    match_radius: radius,
                    confidence: 0.5 + 0.1 * radius as f64,
                });
            }
        }
    }

    None
}

/// Fuzzy matching: check if the sphere-1 fragment of the HOSE code
/// matches the sphere-1 fragment of the database pattern.
fn fuzzy_hose_match(hose: &str, pattern: &str) -> bool {
    let hose_parts: Vec<&str> = hose.split('/').collect();
    let pat_parts: Vec<&str> = pattern.split('/').collect();

    if hose_parts.len() < 2 || pat_parts.len() < 2 {
        return false;
    }

    // Match center element
    if hose_parts[0] != pat_parts[0] {
        return false;
    }

    // Fuzzy match sphere-1: check if sorted neighbor sets overlap
    let hose_neighbors: BTreeSet<&str> = hose_parts[1].split(',').collect();
    let pat_neighbors: BTreeSet<&str> = pat_parts[1].split(',').collect();

    let intersection = hose_neighbors.intersection(&pat_neighbors).count();
    let union = hose_neighbors.union(&pat_neighbors).count();

    if union == 0 {
        return false;
    }

    // Require at least 50% overlap
    (intersection as f64 / union as f64) > 0.5
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hose_codes_ethanol() {
        let mol = Molecule::from_smiles("CCO").unwrap();
        let codes = generate_hose_codes(&mol, 3);

        // Should have codes for all atoms (including H)
        assert_eq!(codes.len(), mol.graph.node_count());

        // All codes should have non-empty first sphere
        for code in &codes {
            assert!(!code.spheres[0].is_empty());
            assert!(!code.full_code.is_empty());
        }
    }

    #[test]
    fn test_hose_codes_benzene() {
        let mol = Molecule::from_smiles("c1ccccc1").unwrap();
        let codes = generate_hose_codes(&mol, 3);

        assert_eq!(codes.len(), mol.graph.node_count());

        // All carbon atoms should have similar HOSE codes due to symmetry
        let carbon_codes: Vec<&HoseCode> = codes.iter().filter(|c| c.element == 6).collect();
        assert!(!carbon_codes.is_empty());
    }

    #[test]
    fn test_hose_codes_deterministic() {
        let mol = Molecule::from_smiles("CC(=O)O").unwrap();
        let codes1 = generate_hose_codes(&mol, 3);
        let codes2 = generate_hose_codes(&mol, 3);

        for (c1, c2) in codes1.iter().zip(codes2.iter()) {
            assert_eq!(c1.full_code, c2.full_code);
        }
    }
}
