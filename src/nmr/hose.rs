//! HOSE (Hierarchically Ordered Spherical description of Environment) code generation.
//!
//! Generates unique environment descriptors for each atom by traversing
//! the molecular graph in concentric spheres (radius 0 to max_radius).

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

        let full_code = format!(
            "{}/{}",
            spheres[0],
            spheres[1..].join("/")
        );

        codes.push(HoseCode {
            atom_index: atom_idx,
            element: center_element,
            spheres,
            full_code,
        });
    }

    codes
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
