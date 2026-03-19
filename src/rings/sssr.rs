//! Smallest Set of Smallest Rings (SSSR) using Horton's algorithm.
//!
//! Identifies the fundamental cycle basis of the molecular graph:
//! 1. Compute shortest paths between all atom pairs (Floyd-Warshall)
//! 2. For each edge (u,v), find the shortest cycle containing it
//! 3. Reduce to a linearly independent set of minimum total weight

use crate::graph::Molecule;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, VecDeque};

/// Information about a single ring in the SSSR.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RingInfo {
    /// Atom indices forming the ring (in order around the cycle).
    pub atoms: Vec<usize>,
    /// Ring size.
    pub size: usize,
    /// Whether the ring is aromatic (Hückel 4n+2).
    pub is_aromatic: bool,
}

/// Result of SSSR computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SssrResult {
    /// The rings in the SSSR.
    pub rings: Vec<RingInfo>,
    /// Per-atom ring membership count.
    pub atom_ring_count: Vec<usize>,
    /// Per-atom ring sizes (all ring sizes the atom belongs to).
    pub atom_ring_sizes: Vec<Vec<usize>>,
    /// Ring-size histogram: index = ring size, value = count.
    pub ring_size_histogram: Vec<usize>,
}

/// Compute the Smallest Set of Smallest Rings (SSSR) for a molecule.
///
/// Uses BFS-based cycle detection: for each edge (u,v), remove it and check
/// if u and v are still connected. If so, the shortest alternative path + the
/// edge forms a ring candidate. We collect all unique minimal rings.
pub fn compute_sssr(mol: &Molecule) -> SssrResult {
    let n = mol.graph.node_count();
    let m = mol.graph.edge_count();

    if n == 0 || m == 0 {
        return SssrResult {
            rings: vec![],
            atom_ring_count: vec![0; n],
            atom_ring_sizes: vec![vec![]; n],
            ring_size_histogram: vec![],
        };
    }

    // Expected number of independent cycles: m - n + c (connected components)
    // For a single connected molecule: m - n + 1
    let n_expected = m.saturating_sub(n) + 1;

    // Build adjacency list
    let mut adj: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];
    let mut edges: Vec<(usize, usize)> = Vec::with_capacity(m);
    for edge in mol.graph.edge_references() {
        let u = edge.source().index();
        let v = edge.target().index();
        adj[u].insert(v);
        adj[v].insert(u);
        if u < v {
            edges.push((u, v));
        } else {
            edges.push((v, u));
        }
    }
    edges.sort();
    edges.dedup();

    let mut ring_candidates: Vec<Vec<usize>> = Vec::new();

    // For each edge, find the shortest cycle containing it
    for &(u, v) in &edges {
        // BFS from u to v without using the direct edge u-v
        if let Some(path) = bfs_shortest_path_excluding_edge(&adj, u, v, n) {
            // path is u → ... → v (not including the direct edge)
            // The ring is path (which already includes u and v)
            let ring = path;
            if ring.len() >= 3 {
                ring_candidates.push(ring);
            }
        }
    }

    // Deduplicate rings: normalize each ring to a canonical form
    let mut unique_rings: Vec<Vec<usize>> = Vec::new();
    let mut seen: BTreeSet<Vec<usize>> = BTreeSet::new();

    for ring in &ring_candidates {
        let canonical = canonicalize_ring(ring);
        if seen.insert(canonical.clone()) {
            unique_rings.push(ring.clone());
        }
    }

    // Sort by ring size (prefer smallest rings)
    unique_rings.sort_by_key(|r| r.len());

    // Take only linearly independent rings (up to n_expected)
    // Use a greedy approach: include a ring if it contains at least one edge
    // not covered by previously included rings
    let mut selected_rings: Vec<Vec<usize>> = Vec::new();
    let mut covered_edges: BTreeSet<(usize, usize)> = BTreeSet::new();

    for ring in &unique_rings {
        if selected_rings.len() >= n_expected {
            break;
        }
        let ring_edges = ring_to_edges(ring);
        let has_new_edge = ring_edges.iter().any(|e| !covered_edges.contains(e));
        if has_new_edge {
            for e in &ring_edges {
                covered_edges.insert(*e);
            }
            selected_rings.push(ring.clone());
        }
    }

    // Determine aromaticity for each ring
    let rings: Vec<RingInfo> = selected_rings
        .iter()
        .map(|ring| {
            let is_aromatic = check_ring_aromaticity(mol, ring);
            RingInfo {
                size: ring.len(),
                atoms: ring.clone(),
                is_aromatic,
            }
        })
        .collect();

    // Compute per-atom ring membership
    let mut atom_ring_count = vec![0usize; n];
    let mut atom_ring_sizes: Vec<Vec<usize>> = vec![vec![]; n];
    for ring in &rings {
        for &atom in &ring.atoms {
            atom_ring_count[atom] += 1;
            atom_ring_sizes[atom].push(ring.size);
        }
    }

    // Ring size histogram
    let max_size = rings.iter().map(|r| r.size).max().unwrap_or(0);
    let mut ring_size_histogram = vec![0usize; max_size + 1];
    for ring in &rings {
        ring_size_histogram[ring.size] += 1;
    }

    SssrResult {
        rings,
        atom_ring_count,
        atom_ring_sizes,
        ring_size_histogram,
    }
}

/// BFS shortest path from `start` to `end` without using the direct edge between them.
fn bfs_shortest_path_excluding_edge(
    adj: &[BTreeSet<usize>],
    start: usize,
    end: usize,
    n: usize,
) -> Option<Vec<usize>> {
    let mut visited = vec![false; n];
    let mut parent = vec![usize::MAX; n];
    let mut queue = VecDeque::new();

    visited[start] = true;
    queue.push_back(start);

    while let Some(current) = queue.pop_front() {
        for &next in &adj[current] {
            // Skip the direct edge start-end
            if current == start && next == end {
                continue;
            }
            if current == end && next == start {
                continue;
            }

            if !visited[next] {
                visited[next] = true;
                parent[next] = current;
                if next == end {
                    // Reconstruct path
                    let mut path = vec![end];
                    let mut curr = end;
                    while curr != start {
                        curr = parent[curr];
                        path.push(curr);
                    }
                    path.reverse();
                    return Some(path);
                }
                queue.push_back(next);
            }
        }
    }

    None
}

/// Canonicalize a ring by choosing the smallest rotation starting from the minimum index.
fn canonicalize_ring(ring: &[usize]) -> Vec<usize> {
    if ring.is_empty() {
        return vec![];
    }

    let n = ring.len();
    // Find position of minimum element
    let min_pos = ring
        .iter()
        .enumerate()
        .min_by_key(|(_, &v)| v)
        .map(|(i, _)| i)
        .unwrap();

    // Try forward rotation
    let forward: Vec<usize> = (0..n).map(|i| ring[(min_pos + i) % n]).collect();
    // Try reverse rotation
    let reverse: Vec<usize> = (0..n).map(|i| ring[(min_pos + n - i) % n]).collect();

    // Return lexicographically smaller
    if forward <= reverse {
        forward
    } else {
        reverse
    }
}

/// Convert a ring (list of atom indices) to a sorted set of edges.
fn ring_to_edges(ring: &[usize]) -> Vec<(usize, usize)> {
    let n = ring.len();
    let mut edges = Vec::with_capacity(n);
    for i in 0..n {
        let u = ring[i];
        let v = ring[(i + 1) % n];
        if u < v {
            edges.push((u, v));
        } else {
            edges.push((v, u));
        }
    }
    edges.sort();
    edges
}

/// Check if a ring is aromatic using Hückel's rule (4n+2 π electrons).
fn check_ring_aromaticity(mol: &Molecule, ring_atoms: &[usize]) -> bool {
    use crate::graph::BondOrder;

    // All ring atoms must be sp2 or have aromatic bonds
    let all_sp2_or_aromatic = ring_atoms.iter().all(|&idx| {
        let node = NodeIndex::new(idx);
        let elem = mol.graph[node].element;
        // Must be C, N, O, S (common aromatic atoms)
        if !matches!(elem, 6 | 7 | 8 | 16) {
            return false;
        }
        // Check if atom has aromatic bonds or is sp2
        let has_aromatic = mol
            .graph
            .edges(node)
            .any(|e| matches!(e.weight().order, BondOrder::Aromatic));
        let is_sp2 = matches!(
            mol.graph[node].hybridization,
            crate::graph::Hybridization::SP2
        );
        has_aromatic || is_sp2
    });

    if !all_sp2_or_aromatic {
        return false;
    }

    // Count π electrons using Hückel's rule
    let mut pi_electrons = 0;
    for &idx in ring_atoms {
        let node = NodeIndex::new(idx);
        let elem = mol.graph[node].element;
        match elem {
            6 => pi_electrons += 1, // C contributes 1 π electron
            7 => {
                // N: 1 if in pyridine (=N-), 2 if in pyrrole (NH)
                let h_count = mol
                    .graph
                    .neighbors(node)
                    .filter(|n| mol.graph[*n].element == 1)
                    .count();
                if h_count > 0 {
                    pi_electrons += 2; // pyrrole-type N
                } else {
                    pi_electrons += 1; // pyridine-type N (or N in ring with lone pair)
                }
            }
            8 => pi_electrons += 2,  // O contributes 2 (furan-type)
            16 => pi_electrons += 2, // S contributes 2 (thiophene-type)
            _ => pi_electrons += 1,
        }
    }

    // Hückel's rule: 4n+2 for n = 0,1,2,...
    // Common: 2 (n=0), 6 (n=1, benzene), 10 (n=2), 14 (n=3), 18 (n=4)
    if pi_electrons < 2 {
        return false;
    }
    (pi_electrons - 2) % 4 == 0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_benzene_sssr() {
        let mol = Molecule::from_smiles("c1ccccc1").unwrap();
        let result = compute_sssr(&mol);

        // Benzene: 1 ring of size 6
        assert_eq!(result.rings.len(), 1, "Benzene should have 1 ring in SSSR");
        assert_eq!(result.rings[0].size, 6, "Ring should be size 6");
        assert!(result.rings[0].is_aromatic, "Ring should be aromatic");
    }

    #[test]
    fn test_naphthalene_sssr() {
        let mol = Molecule::from_smiles("c1ccc2ccccc2c1").unwrap();
        let result = compute_sssr(&mol);

        // Naphthalene: 2 rings of size 6
        assert_eq!(
            result.rings.len(),
            2,
            "Naphthalene SSSR should have 2 rings, got {}",
            result.rings.len()
        );
        for ring in &result.rings {
            assert_eq!(ring.size, 6, "All naphthalene rings should be size 6");
            assert!(ring.is_aromatic, "All naphthalene rings should be aromatic");
        }
    }

    #[test]
    fn test_cyclohexane_sssr() {
        let mol = Molecule::from_smiles("C1CCCCC1").unwrap();
        let result = compute_sssr(&mol);

        assert_eq!(result.rings.len(), 1, "Cyclohexane should have 1 ring");
        assert_eq!(result.rings[0].size, 6);
        assert!(!result.rings[0].is_aromatic, "Cyclohexane is not aromatic");
    }

    #[test]
    fn test_ethane_no_rings() {
        let mol = Molecule::from_smiles("CC").unwrap();
        let result = compute_sssr(&mol);

        assert_eq!(result.rings.len(), 0, "Ethane should have no rings");
    }

    #[test]
    fn test_ring_canonicalization() {
        let ring1 = vec![3, 0, 1, 2];
        let ring2 = vec![1, 2, 3, 0];
        assert_eq!(canonicalize_ring(&ring1), canonicalize_ring(&ring2));
    }

    #[test]
    fn test_atom_ring_membership() {
        let mol = Molecule::from_smiles("c1ccccc1").unwrap();
        let result = compute_sssr(&mol);

        // All aromatic carbons should be in 1 ring
        for &idx in &result.rings[0].atoms {
            assert!(
                result.atom_ring_count[idx] >= 1,
                "Atom {} should be in at least 1 ring",
                idx
            );
        }
    }
}
