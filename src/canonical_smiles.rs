//! Canonical SMILES generation using the Morgan algorithm.
//!
//! Implements deterministic SMILES output via:
//! 1. Morgan-like atom ranking for canonical numbering
//! 2. DFS traversal of the molecular graph
//! 3. Ring closure handling with minimal digit assignment
//!
//! The output is *canonical* in the sense that it only depends on the
//! molecular graph and not on the order in which the atoms were supplied:
//! the ring-closure digits are handed out in the order in which the ring
//! closures are opened in the output string and every ring-closure bond is
//! written at a deterministic position.

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

    // 2. The starting atom is the one with the lowest canonical rank.  When
    //    several atoms share that rank they are symmetry equivalent, so every
    //    one of them is tried and the lexicographically smallest SMILES is
    //    kept: the result then does not depend on the order in which the
    //    atoms were supplied.
    let heavy: Vec<usize> = (0..n)
        .filter(|&i| mol.graph[NodeIndex::new(i)].element != 1)
        .collect();
    if heavy.is_empty() {
        return String::new();
    }
    let min_key = heavy
        .iter()
        .map(|&i| (ranks[i], mol.graph[NodeIndex::new(i)].element))
        .min()
        .unwrap();
    let starts: Vec<usize> = heavy
        .iter()
        .copied()
        .filter(|&i| (ranks[i], mol.graph[NodeIndex::new(i)].element) == min_key)
        .collect();

    let mut best: Option<String> = None;
    for start in starts {
        let smiles = canonical_from_start(mol, &ranks, start);
        if best.as_ref().is_none_or(|b| smiles < *b) {
            best = Some(smiles);
        }
    }
    best.unwrap_or_default()
}

/// Build the canonical SMILES starting the traversal at `start`.
fn canonical_from_start(mol: &Molecule, ranks: &[u64], start: usize) -> String {
    let n = mol.graph.node_count();

    // Ring-closure bonds: the back-edges of the depth-first traversal.
    let ring_bonds = find_ring_closure_bonds(mol, start, ranks);

    // Build the SMILES string; the ring-closure digits are assigned while
    // writing, in the order in which the closures are opened.
    let mut writer = SmilesWriter {
        mol,
        ranks: ranks.to_vec(),
        ring_bonds,
        ring_numbers: HashMap::new(),
        next_ring_num: 1,
        smiles: String::new(),
    };
    let mut visited = vec![false; n];
    writer.write_node(start, &mut visited);
    writer.smiles
}

/// Compute the canonical atom labelling.
///
/// The Morgan invariants are refined with the Weisfeiler–Leman operator
/// until the partition is stable.  Atoms that remain equivalent (they are
/// related by a molecular symmetry) are separated by *individualisation*:
/// one atom of the first non-singleton cell is given a unique rank, the
/// refinement is repeated and the individualisation that produces the
/// smallest signature is kept.  The resulting labels depend only on the
/// molecular graph, never on the order in which the atoms were supplied.
fn compute_morgan_ranks(mol: &Molecule) -> Vec<u64> {
    let n = mol.graph.node_count();
    let mut ranks = refine_to_stability(mol, initial_ranks(mol));

    for _ in 0..n {
        // The first cell (ordered by rank value) that holds more than one
        // atom; the cell order is a property of the graph, not of the input
        // order.
        let mut cells: std::collections::BTreeMap<u64, Vec<usize>> =
            std::collections::BTreeMap::new();
        for i in 0..n {
            cells.entry(ranks[i]).or_default().push(i);
        }
        let members = match cells.values().find(|v| v.len() > 1) {
            Some(v) => v.clone(),
            None => break,
        };

        let mut best: Option<(Vec<u64>, Vec<u64>)> = None;
        for &atom in &members {
            let mut candidate = ranks.clone();
            candidate[atom] = unused_rank_value(&ranks);
            let candidate = refine_to_stability(mol, candidate);
            let mut signature = candidate.clone();
            signature.sort_unstable();
            let better = match &best {
                None => true,
                Some((sig, _)) => signature < *sig,
            };
            if better {
                best = Some((signature, candidate));
            }
        }
        if let Some((_, candidate)) = best {
            ranks = candidate;
        } else {
            break;
        }
    }
    ranks
}

/// A rank value that does not occur in `ranks`.
fn unused_rank_value(ranks: &[u64]) -> u64 {
    let used: std::collections::HashSet<u64> = ranks.iter().copied().collect();
    let mut candidate = u64::MAX;
    while used.contains(&candidate) {
        candidate -= 1;
    }
    candidate
}

/// Initial invariants: (atomic_number, degree, h_count, formal_charge, is_aromatic).
fn initial_ranks(mol: &Molecule) -> Vec<u64> {
    let n = mol.graph.node_count();
    (0..n)
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
        .collect()
}

/// One Weisfeiler–Leman refinement step: hash the current rank together
/// with the sorted ranks of the neighbours.
fn refine_ranks(mol: &Molecule, ranks: &[u64]) -> Vec<u64> {
    let n = mol.graph.node_count();
    let mut new_ranks = vec![0u64; n];
    for i in 0..n {
        let idx = NodeIndex::new(i);
        let mut neighbor_ranks: Vec<u64> = mol
            .graph
            .neighbors(idx)
            .map(|nb| ranks[nb.index()])
            .collect();
        neighbor_ranks.sort();

        let mut hash: u64 = 0xcbf29ce484222325;
        for &val in std::iter::once(&ranks[i]).chain(neighbor_ranks.iter()) {
            for byte in val.to_le_bytes() {
                hash ^= byte as u64;
                hash = hash.wrapping_mul(0x100000001b3);
            }
        }
        new_ranks[i] = hash;
    }
    new_ranks
}

/// Refine until the partition (the number of distinct ranks) is stable.
fn refine_to_stability(mol: &Molecule, mut ranks: Vec<u64>) -> Vec<u64> {
    let n = mol.graph.node_count();
    for _ in 0..n.min(20) {
        let new_ranks = refine_ranks(mol, &ranks);
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

/// Canonical key of an undirected bond.
fn bond_key(a: usize, b: usize) -> (usize, usize) {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}

/// Find the bonds that close a ring (the back-edges of the depth-first
/// traversal that starts at `start` and follows `sorted_neighbors`).
fn find_ring_closure_bonds(mol: &Molecule, start: usize, ranks: &[u64]) -> HashSet<(usize, usize)> {
    let n = mol.graph.node_count();
    let mut ring_bonds = HashSet::new();
    let mut in_stack = vec![false; n];
    let mut dfs_visited = vec![false; n];
    let mut dfs_parent = vec![usize::MAX; n];

    dfs_visited[start] = true;
    in_stack[start] = true;
    let mut actual_stack = vec![(start, 0usize)]; // (node, neighbor_idx)
    while let Some(&mut (node, ref mut ni)) = actual_stack.last_mut() {
        let neighbors: Vec<usize> = sorted_neighbors(mol, node, ranks)
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
                // Back-edge: this bond closes a ring.
                ring_bonds.insert(bond_key(node, nb));
            }
        } else {
            in_stack[node] = false;
            actual_stack.pop();
        }
    }
    ring_bonds
}

/// State of the SMILES writer.
struct SmilesWriter<'a> {
    mol: &'a Molecule,
    ranks: Vec<u64>,
    ring_bonds: HashSet<(usize, usize)>,
    ring_numbers: HashMap<(usize, usize), u32>,
    next_ring_num: u32,
    smiles: String,
}

impl SmilesWriter<'_> {
    /// Write the symbol of one atom.
    ///
    /// Charged atoms and aromatic heteroatoms with an explicit hydrogen are
    /// written as bracket atoms (`[O-]`, `[N+]`, `[nH]`, …): the charge and
    /// the hydrogen count cannot be inferred from the valence rules and
    /// would otherwise be lost, so the SMILES would describe a different
    /// molecule.
    fn write_atom_symbol(&mut self, idx: NodeIndex) {
        let atom = &self.mol.graph[idx];
        let z = atom.element;
        let aromatic = is_aromatic_organic(self.mol, idx);
        let n_h = self
            .mol
            .graph
            .neighbors(idx)
            .filter(|nb| self.mol.graph[*nb].element == 1)
            .count();
        let charge = atom.formal_charge;
        let needs_bracket = charge != 0 || (aromatic && z != 6 && n_h > 0);

        let symbol = element_symbol(z);
        if !needs_bracket {
            if aromatic {
                self.smiles.push_str(&symbol.to_lowercase());
            } else {
                self.smiles.push_str(symbol);
            }
            return;
        }

        self.smiles.push('[');
        if aromatic {
            self.smiles.push_str(&symbol.to_lowercase());
        } else {
            self.smiles.push_str(symbol);
        }
        if n_h == 1 {
            self.smiles.push('H');
        } else if n_h > 1 {
            self.smiles.push('H');
            self.smiles.push_str(&n_h.to_string());
        }
        if charge > 0 {
            for _ in 0..charge {
                self.smiles.push('+');
            }
        } else if charge < 0 {
            for _ in 0..(-charge) {
                self.smiles.push('-');
            }
        }
        self.smiles.push(']');
    }

    /// Ring-closure bonds incident to `node`, sorted deterministically by the
    /// canonical rank of the partner atom (then by index).
    fn ring_partners(&self, node: usize) -> Vec<usize> {
        let mut partners: Vec<usize> = self
            .ring_bonds
            .iter()
            .filter_map(|&(a, b)| {
                if a == node {
                    Some(b)
                } else if b == node {
                    Some(a)
                } else {
                    None
                }
            })
            .collect();
        partners.sort_by_key(|&p| (self.ranks[p], p));
        partners
    }

    /// Write one atom (and its subtree) to the SMILES string.
    fn write_node(&mut self, node: usize, visited: &mut [bool]) {
        visited[node] = true;
        let idx = NodeIndex::new(node);
        let atom = &self.mol.graph[idx];

        // Skip explicit hydrogens (they're implicit in SMILES)
        if atom.element == 1 {
            return;
        }

        // Write atom symbol
        self.write_atom_symbol(idx);

        // Ring closures: the first end that is written *opens* the closure
        // and takes the next free digit, the second end closes it.
        let partners = self.ring_partners(node);
        let mut digits: Vec<u32> = Vec::with_capacity(partners.len());
        for other in partners {
            let key = bond_key(node, other);
            if visited[other] {
                // Closing end: reuse the digit of the opening end.
                if let Some(&num) = self.ring_numbers.get(&key) {
                    digits.push(num);
                }
            } else {
                // Opening end: hand out the next digit.
                let num = self.next_ring_num;
                self.next_ring_num += 1;
                self.ring_numbers.insert(key, num);
                digits.push(num);
            }
        }
        digits.sort_unstable();
        for num in digits {
            if num < 10 {
                self.smiles.push(char::from(b'0' + num as u8));
            } else {
                self.smiles.push('%');
                self.smiles.push_str(&num.to_string());
            }
        }

        // Get neighbors (excluding visited ones, H atoms, and ring closure edges)
        let neighbors = sorted_neighbors(self.mol, node, &self.ranks);
        let unvisited: Vec<usize> = neighbors
            .iter()
            .filter(|&&nb| {
                if visited[nb] || self.mol.graph[NodeIndex::new(nb)].element == 1 {
                    return false;
                }
                // Don't follow back-edges (ring closure edges) as tree edges
                !self.ring_bonds.contains(&bond_key(node, nb))
            })
            .copied()
            .collect();

        // Write branches and continue DFS
        if unvisited.len() <= 1 {
            // No branching needed
            for &nb in &unvisited {
                // Write bond symbol if not single/aromatic
                self.write_bond_symbol(node, nb);
                self.write_node(nb, visited);
            }
        } else {
            // Multiple branches: main chain = last, branches = parenthesized
            for (i, &nb) in unvisited.iter().enumerate() {
                if i < unvisited.len() - 1 {
                    self.smiles.push('(');
                    self.write_bond_symbol(node, nb);
                    self.write_node(nb, visited);
                    self.smiles.push(')');
                } else {
                    self.write_bond_symbol(node, nb);
                    self.write_node(nb, visited);
                }
            }
        }
    }

    /// Write the bond symbol between two atoms.
    ///
    /// Single bonds are implicit, *except* between two aromatic atoms: an
    /// aromatic–aromatic single bond (e.g. the central bond of biphenyl)
    /// must be written explicitly, otherwise the bond would be re-perceived
    /// as aromatic when the SMILES is parsed again.
    fn write_bond_symbol(&mut self, from: usize, to: usize) {
        let from_idx = NodeIndex::new(from);
        let to_idx = NodeIndex::new(to);
        for edge in self.mol.graph.edges(from_idx) {
            let other = if edge.source() == from_idx {
                edge.target()
            } else {
                edge.source()
            };
            if other == to_idx {
                match edge.weight().order {
                    BondOrder::Double => self.smiles.push('='),
                    BondOrder::Triple => self.smiles.push('#'),
                    BondOrder::Single => {
                        let aromatic_bond = is_aromatic_organic(self.mol, from_idx)
                            && is_aromatic_organic(self.mol, to_idx);
                        if aromatic_bond {
                            self.smiles.push('-');
                        }
                    }
                    _ => {} // Aromatic and unknown are implicit
                }
                return;
            }
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
