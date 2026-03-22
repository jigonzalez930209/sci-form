//! Canonical SMILES generation using the Morgan algorithm.
//!
//! Implements deterministic SMILES output via:
//! 1. Morgan-like atom ranking for canonical numbering
//! 2. DFS traversal of the molecular graph
//! 3. Ring closure handling with minimal digit assignment

use crate::graph::{BondOrder, Molecule};
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use std::collections::{HashMap, HashSet};

/// Generate a canonical SMILES string for a molecule.
///
/// Uses Morgan-like invariants to produce a deterministic atom ordering,
/// then performs DFS traversal to build the SMILES string.
pub fn to_canonical_smiles(mol: &Molecule) -> String {
    let n = mol.graph.node_count();
    if n == 0 {
        return String::new();
    }

    // 1. Compute Morgan-like canonical ranking
    let ranks = compute_morgan_ranks(mol);

    // 2. Find the starting atom (lowest rank, break ties by atomic number)
    let start = (0..n)
        .filter(|&i| mol.graph[NodeIndex::new(i)].element != 1)
        .min_by_key(|&i| {
            let atom = &mol.graph[NodeIndex::new(i)];
            (ranks[i], atom.element)
        })
        .unwrap_or(0);

    // 3. DFS traversal to build SMILES
    let mut smiles = String::new();
    let mut visited = vec![false; n];
    let mut ring_closures: HashMap<(usize, usize), u32> = HashMap::new();
    let mut next_ring_num = 1u32;

    // Pre-compute ring closure assignments by detecting back-edges in DFS
    {
        let mut in_stack = vec![false; n];
        let mut dfs_visited = vec![false; n];
        let mut dfs_parent = vec![usize::MAX; n];

        // Find back-edges using iterative DFS
        dfs_visited[start] = true;
        in_stack[start] = true;
        let mut actual_stack = vec![(start, 0usize)]; // (node, neighbor_idx)
        while let Some(&mut (node, ref mut ni)) = actual_stack.last_mut() {
            let neighbors: Vec<usize> = sorted_neighbors(mol, node, &ranks)
                .into_iter()
                .filter(|&nb| mol.graph[NodeIndex::new(nb)].element != 1)
                .collect();
            if *ni < neighbors.len() {
                let nb = neighbors[*ni];
                *ni += 1;
                if !dfs_visited[nb] {
                    dfs_visited[nb] = true;
                    in_stack[nb] = true;
                    dfs_parent[nb] = node;
                    actual_stack.push((nb, 0));
                } else if nb != dfs_parent[node] && in_stack[nb] {
                    // Back-edge: need a ring closure between nb and node
                    let key = if nb < node { (nb, node) } else { (node, nb) };
                    if let std::collections::hash_map::Entry::Vacant(entry) =
                        ring_closures.entry(key)
                    {
                        entry.insert(next_ring_num);
                        next_ring_num += 1;
                    }
                }
            } else {
                in_stack[node] = false;
                actual_stack.pop();
            }
        }
    }

    // 4. Build SMILES string via DFS
    dfs_write_smiles(
        mol,
        start,
        &mut visited,
        &ranks,
        &ring_closures,
        &mut smiles,
    );

    smiles
}

/// Compute Morgan-like canonical ranks for atoms.
///
/// Initial invariant: (atomic_number, degree, h_count, formal_charge, is_aromatic).
/// Iterate until ranks stabilize (or max iterations reached).
fn compute_morgan_ranks(mol: &Molecule) -> Vec<u64> {
    let n = mol.graph.node_count();

    // Initial invariants
    let mut ranks: Vec<u64> = (0..n)
        .map(|i| {
            let idx = NodeIndex::new(i);
            let atom = &mol.graph[idx];
            let degree = mol.graph.edges(idx).count() as u64;
            let h_count = mol
                .graph
                .neighbors(idx)
                .filter(|nb| mol.graph[*nb].element == 1)
                .count() as u64;
            let is_aromatic = mol
                .graph
                .edges(idx)
                .any(|e| matches!(e.weight().order, BondOrder::Aromatic))
                as u64;
            // Combine into initial rank
            (atom.element as u64) * 1000000
                + degree * 10000
                + h_count * 100
                + (atom.formal_charge as i64 + 10) as u64 * 10
                + is_aromatic
        })
        .collect();

    // Iterate Morgan refinement
    for _ in 0..n.min(20) {
        let mut new_ranks = vec![0u64; n];
        for i in 0..n {
            let idx = NodeIndex::new(i);
            let mut neighbor_ranks: Vec<u64> = mol
                .graph
                .neighbors(idx)
                .map(|nb| ranks[nb.index()])
                .collect();
            neighbor_ranks.sort();

            // Hash current rank with sorted neighbor ranks
            let mut hash: u64 = 0xcbf29ce484222325;
            for &val in std::iter::once(&ranks[i]).chain(neighbor_ranks.iter()) {
                for byte in val.to_le_bytes() {
                    hash ^= byte as u64;
                    hash = hash.wrapping_mul(0x100000001b3);
                }
            }
            new_ranks[i] = hash;
        }

        // Check if number of distinct ranks changed
        let old_distinct: HashSet<u64> = ranks.iter().copied().collect();
        let new_distinct: HashSet<u64> = new_ranks.iter().copied().collect();
        ranks = new_ranks;
        if old_distinct.len() == new_distinct.len() {
            break;
        }
    }

    ranks
}

/// Get sorted neighbor indices based on canonical ranks.
fn sorted_neighbors(mol: &Molecule, node: usize, ranks: &[u64]) -> Vec<usize> {
    let idx = NodeIndex::new(node);
    let mut neighbors: Vec<usize> = mol.graph.neighbors(idx).map(|nb| nb.index()).collect();
    neighbors.sort_by_key(|&nb| ranks[nb]);
    neighbors
}

/// Element symbol for common organic atoms.
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
        _ => "?",
    }
}

/// Check if an atom can be written in the aromatic organic subset (lowercase).
fn is_aromatic_organic(mol: &Molecule, idx: NodeIndex) -> bool {
    let atom = &mol.graph[idx];
    let has_aromatic = mol
        .graph
        .edges(idx)
        .any(|e| matches!(e.weight().order, BondOrder::Aromatic));
    has_aromatic && matches!(atom.element, 6 | 7 | 8 | 16)
}

/// DFS-based SMILES writing.
fn dfs_write_smiles(
    mol: &Molecule,
    node: usize,
    visited: &mut [bool],
    ranks: &[u64],
    ring_closures: &HashMap<(usize, usize), u32>,
    smiles: &mut String,
) {
    visited[node] = true;
    let idx = NodeIndex::new(node);
    let atom = &mol.graph[idx];

    // Skip explicit hydrogens (they're implicit in SMILES)
    if atom.element == 1 {
        return;
    }

    // Write atom symbol
    if is_aromatic_organic(mol, idx) {
        let sym = element_symbol(atom.element).to_lowercase();
        smiles.push_str(&sym);
    } else {
        smiles.push_str(element_symbol(atom.element));
    }

    // Write ring closures for this atom (both opening and closing ends)
    for (&(a, b), &ring_num) in ring_closures {
        if a == node || b == node {
            let other = if a == node { b } else { a };
            if visited[other] {
                // Closing end: write bond symbol for the ring closure if needed
                write_bond_symbol(mol, node, other, smiles);
            }
            // Write ring digit at both the opening and closing ends
            if ring_num < 10 {
                smiles.push(char::from(b'0' + ring_num as u8));
            } else {
                smiles.push('%');
                smiles.push_str(&ring_num.to_string());
            }
        }
    }

    // Get neighbors (excluding visited ones, H atoms, and ring closure edges)
    let neighbors = sorted_neighbors(mol, node, ranks);
    let unvisited: Vec<usize> = neighbors
        .iter()
        .filter(|&&nb| {
            if visited[nb] || mol.graph[NodeIndex::new(nb)].element == 1 {
                return false;
            }
            // Don't follow back-edges (ring closure edges) as tree edges
            let key = if node < nb { (node, nb) } else { (nb, node) };
            !ring_closures.contains_key(&key)
        })
        .copied()
        .collect();

    // Write branches and continue DFS
    if unvisited.len() <= 1 {
        // No branching needed
        for &nb in &unvisited {
            // Write bond symbol if not single/aromatic
            write_bond_symbol(mol, node, nb, smiles);
            dfs_write_smiles(mol, nb, visited, ranks, ring_closures, smiles);
        }
    } else {
        // Multiple branches: main chain = last, branches = parenthesized
        for (i, &nb) in unvisited.iter().enumerate() {
            if i < unvisited.len() - 1 {
                smiles.push('(');
                write_bond_symbol(mol, node, nb, smiles);
                dfs_write_smiles(mol, nb, visited, ranks, ring_closures, smiles);
                smiles.push(')');
            } else {
                write_bond_symbol(mol, node, nb, smiles);
                dfs_write_smiles(mol, nb, visited, ranks, ring_closures, smiles);
            }
        }
    }
}

/// Write bond symbol between two atoms (omit for single/aromatic).
fn write_bond_symbol(mol: &Molecule, from: usize, to: usize, smiles: &mut String) {
    let from_idx = NodeIndex::new(from);
    let to_idx = NodeIndex::new(to);
    for edge in mol.graph.edges(from_idx) {
        let other = if edge.source() == from_idx {
            edge.target()
        } else {
            edge.source()
        };
        if other == to_idx {
            match edge.weight().order {
                BondOrder::Double => smiles.push('='),
                BondOrder::Triple => smiles.push('#'),
                _ => {} // Single and aromatic are implicit
            }
            return;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_canonical_ethanol() {
        let mol = Molecule::from_smiles("CCO").unwrap();
        let smiles = to_canonical_smiles(&mol);
        assert!(!smiles.is_empty());
        // Should parse back to the same molecule
        let mol2 = Molecule::from_smiles(&smiles).unwrap();
        // Same number of heavy atoms
        let heavy1 = mol
            .graph
            .node_indices()
            .filter(|&i| mol.graph[i].element != 1)
            .count();
        let heavy2 = mol2
            .graph
            .node_indices()
            .filter(|&i| mol2.graph[i].element != 1)
            .count();
        assert_eq!(heavy1, heavy2);
    }

    #[test]
    fn test_canonical_deterministic() {
        let mol1 = Molecule::from_smiles("CCO").unwrap();
        let mol2 = Molecule::from_smiles("OCC").unwrap();
        let s1 = to_canonical_smiles(&mol1);
        let s2 = to_canonical_smiles(&mol2);
        // Both should produce the same canonical SMILES
        assert_eq!(
            s1, s2,
            "Canonical SMILES should be deterministic: {} vs {}",
            s1, s2
        );
    }

    #[test]
    fn test_canonical_benzene() {
        let mol = Molecule::from_smiles("c1ccccc1").unwrap();
        let smiles = to_canonical_smiles(&mol);
        assert!(!smiles.is_empty());
        // Parse back should give same number of atoms
        let mol2 = Molecule::from_smiles(&smiles).unwrap();
        let heavy1 = mol
            .graph
            .node_indices()
            .filter(|&i| mol.graph[i].element != 1)
            .count();
        let heavy2 = mol2
            .graph
            .node_indices()
            .filter(|&i| mol2.graph[i].element != 1)
            .count();
        assert_eq!(heavy1, heavy2);
    }

    #[test]
    fn test_canonical_double_bond() {
        let mol = Molecule::from_smiles("C=C").unwrap();
        let smiles = to_canonical_smiles(&mol);
        assert!(
            smiles.contains('='),
            "Should contain double bond: {}",
            smiles
        );
    }

    #[test]
    fn test_canonical_acetic_acid() {
        let mol = Molecule::from_smiles("CC(=O)O").unwrap();
        let smiles = to_canonical_smiles(&mol);
        assert!(!smiles.is_empty());
        assert!(
            smiles.contains('='),
            "Acetic acid should have a double bond: {}",
            smiles
        );
    }
}
