//! Substructure matching: match a SmartsPattern against a Molecule graph.

use super::parser::*;
use crate::graph::{BondOrder, ChiralType, Hybridization, Molecule};
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;

/// Find all substructure matches of `pattern` in `mol`.
/// Returns a vec of mappings: each mapping is pattern_atom_idx → molecule_atom_idx.
pub fn substruct_match(mol: &Molecule, pattern: &SmartsPattern) -> Vec<Vec<usize>> {
    if pattern.atoms.is_empty() {
        return vec![];
    }

    // Pre-compute ring info (lazily, only if needed)
    let needs_ring = pattern_needs_ring(pattern);
    let ring_info = if needs_ring {
        Some(compute_ring_info(mol))
    } else {
        None
    };

    substruct_match_inner(mol, pattern, ring_info.as_ref())
}

/// Find all substructure matches using pre-computed ring info.
/// This avoids recomputing SSSR when matching many patterns against the same molecule.
pub fn substruct_match_with_ring_info(
    mol: &Molecule,
    pattern: &SmartsPattern,
    ring_info: &RingInfo,
) -> Vec<Vec<usize>> {
    if pattern.atoms.is_empty() {
        return vec![];
    }
    substruct_match_inner(mol, pattern, Some(ring_info))
}

/// Compute ring info for a molecule (SSSR + ring membership + ring sizes).
/// Call once per molecule, then pass to `substruct_match_with_ring_info`.
pub fn precompute_ring_info(mol: &Molecule) -> RingInfo {
    compute_ring_info(mol)
}

fn substruct_match_inner(
    mol: &Molecule,
    pattern: &SmartsPattern,
    ring_info: Option<&RingInfo>,
) -> Vec<Vec<usize>> {
    let n_pat = pattern.atoms.len();
    let n_mol = mol.graph.node_count();

    // Build adjacency for the pattern
    let mut pat_adj: Vec<Vec<(usize, usize)>> = vec![vec![]; n_pat]; // (neighbor, bond_idx)
    for (bi, bond) in pattern.bonds.iter().enumerate() {
        pat_adj[bond.from].push((bond.to, bi));
        pat_adj[bond.to].push((bond.from, bi));
    }

    let mut results = Vec::new();
    let mut mapping = vec![usize::MAX; n_pat];
    let mut used = vec![false; n_mol];

    // Try each molecule atom as starting point for pattern atom 0
    for start in 0..n_mol {
        mapping[0] = start;
        used[start] = true;
        if atom_matches(
            mol,
            NodeIndex::new(start),
            &pattern.atoms[0].query,
            ring_info,
        ) {
            backtrack(
                mol,
                pattern,
                &pat_adj,
                &mut mapping,
                &mut used,
                1,
                ring_info,
                &mut results,
            );
        }
        used[start] = false;
        mapping[0] = usize::MAX;
    }
    results
}

fn backtrack(
    mol: &Molecule,
    pattern: &SmartsPattern,
    pat_adj: &[Vec<(usize, usize)>],
    mapping: &mut Vec<usize>,
    used: &mut Vec<bool>,
    depth: usize,
    ring_info: Option<&RingInfo>,
    results: &mut Vec<Vec<usize>>,
) {
    if depth == pattern.atoms.len() {
        results.push(mapping.clone());
        return;
    }

    // Find the next unmapped pattern atom that connects to an already-mapped atom.
    // This ensures we extend the match along the pattern's connectivity.
    let (pat_atom, connected_to, bond_idx) = match find_next_atom(pattern, pat_adj, mapping, depth)
    {
        Some(x) => x,
        None => return,
    };

    let mol_anchor = NodeIndex::new(mapping[connected_to]);

    // Try each neighbor of the molecule anchor atom
    for nb in mol.graph.neighbors(mol_anchor) {
        let ni = nb.index();
        if used[ni] {
            continue;
        }

        // Check atom match
        if !atom_matches(mol, nb, &pattern.atoms[pat_atom].query, ring_info) {
            continue;
        }

        // Check bond match
        let edge = mol.graph.find_edge(mol_anchor, nb).unwrap();
        let mol_bond = &mol.graph[edge];
        if !bond_matches(
            mol_bond.order,
            &pattern.bonds[bond_idx].query,
            mol,
            mol_anchor,
            nb,
            ring_info,
        ) {
            continue;
        }

        // Also verify that all other already-mapped neighbors of pat_atom have matching bonds
        let mut all_bonds_ok = true;
        for &(other_pat, bi) in &pat_adj[pat_atom] {
            if other_pat == connected_to {
                continue;
            }
            if mapping[other_pat] == usize::MAX {
                continue;
            }
            // There should be a bond from nb to mapping[other_pat] in the molecule
            let other_mol = NodeIndex::new(mapping[other_pat]);
            if let Some(e) = mol.graph.find_edge(nb, other_mol) {
                if !bond_matches(
                    mol.graph[e].order,
                    &pattern.bonds[bi].query,
                    mol,
                    nb,
                    other_mol,
                    ring_info,
                ) {
                    all_bonds_ok = false;
                    break;
                }
            } else {
                all_bonds_ok = false;
                break;
            }
        }
        if !all_bonds_ok {
            continue;
        }

        mapping[pat_atom] = ni;
        used[ni] = true;
        backtrack(
            mol,
            pattern,
            pat_adj,
            mapping,
            used,
            depth + 1,
            ring_info,
            results,
        );
        used[ni] = false;
        mapping[pat_atom] = usize::MAX;
    }
}

/// Find the next pattern atom to map: must be unmapped but adjacent to a mapped atom.
fn find_next_atom(
    pattern: &SmartsPattern,
    pat_adj: &[Vec<(usize, usize)>],
    mapping: &[usize],
    _depth: usize,
) -> Option<(usize, usize, usize)> {
    // BFS order: find first unmapped atom connected to a mapped atom
    for pat_idx in 0..pattern.atoms.len() {
        if mapping[pat_idx] != usize::MAX {
            continue;
        }
        for &(neighbor, bond_idx) in &pat_adj[pat_idx] {
            if mapping[neighbor] != usize::MAX {
                return Some((pat_idx, neighbor, bond_idx));
            }
        }
    }
    None
}

/// Check if a molecule atom matches an atom query.
fn atom_matches(
    mol: &Molecule,
    atom: NodeIndex,
    query: &AtomQuery,
    ring_info: Option<&RingInfo>,
) -> bool {
    let a = &mol.graph[atom];
    match query {
        AtomQuery::True => true,
        AtomQuery::Element(z) => a.element == *z && !is_aromatic_atom(mol, atom),
        AtomQuery::AromaticElem(z) => a.element == *z && is_aromatic_atom(mol, atom),
        AtomQuery::AnyAromatic => is_aromatic_atom(mol, atom),
        AtomQuery::AnyAliphatic => !is_aromatic_atom(mol, atom),
        AtomQuery::AtomicNum(z) => a.element == *z,
        AtomQuery::NotAtomicNum(z) => a.element != *z,
        AtomQuery::TotalH(n) => count_h(mol, atom) == *n as usize,
        AtomQuery::TotalDegree(n) => mol.graph.neighbors(atom).count() == *n as usize,
        AtomQuery::HeavyDegree(n) => {
            mol.graph
                .neighbors(atom)
                .filter(|&nb| mol.graph[nb].element != 1)
                .count()
                == *n as usize
        }
        AtomQuery::RingBondCount(n) => {
            if let Some(ri) = ring_info {
                count_ring_bonds(mol, atom, ri) == *n as usize
            } else {
                false
            }
        }
        AtomQuery::InRing => {
            if let Some(ri) = ring_info {
                ri.atom_in_ring[atom.index()]
            } else {
                false
            }
        }
        AtomQuery::RingCount(n) => {
            if let Some(ri) = ring_info {
                let count = ri
                    .rings
                    .iter()
                    .filter(|r| r.contains(&atom.index()))
                    .count();
                count == *n as usize
            } else {
                *n == 0
            }
        }
        AtomQuery::RingSize(n) => {
            if let Some(ri) = ring_info {
                ri.atom_ring_sizes[atom.index()].contains(&(*n as usize))
            } else {
                false
            }
        }
        AtomQuery::RingSizeRange(lo, hi) => {
            if let Some(ri) = ring_info {
                ri.atom_ring_sizes[atom.index()]
                    .iter()
                    .any(|&s| s >= *lo as usize && s <= *hi as usize)
            } else {
                false
            }
        }
        AtomQuery::RingSizeMin(lo) => {
            if let Some(ri) = ring_info {
                ri.atom_ring_sizes[atom.index()]
                    .iter()
                    .any(|&s| s >= *lo as usize)
            } else {
                false
            }
        }
        AtomQuery::FormalCharge(c) => a.formal_charge == *c,
        AtomQuery::Hybridization(n) => matches!(
            (n, &a.hybridization),
            (1, Hybridization::SP) | (2, Hybridization::SP2) | (3, Hybridization::SP3)
        ),
        AtomQuery::Chiral(chiral) => matches!(
            (chiral, &a.chiral_tag),
            (ChiralType::TetrahedralCW, ChiralType::TetrahedralCW)
                | (ChiralType::TetrahedralCCW, ChiralType::TetrahedralCCW)
        ),
        AtomQuery::Recursive(inner) => {
            // The atom must match as atom 0 of the inner pattern
            let matches = substruct_match_from(mol, inner, atom, ring_info);
            !matches.is_empty()
        }
        AtomQuery::And(parts) => parts.iter().all(|q| atom_matches(mol, atom, q, ring_info)),
        AtomQuery::Or(parts) => parts.iter().any(|q| atom_matches(mol, atom, q, ring_info)),
        AtomQuery::Not(inner) => !atom_matches(mol, atom, inner, ring_info),
    }
}

/// Check if a bond matches a bond query.
fn bond_matches(
    order: BondOrder,
    query: &BondQuery,
    _mol: &Molecule,
    from: NodeIndex,
    to: NodeIndex,
    ring_info: Option<&RingInfo>,
) -> bool {
    match query {
        BondQuery::Single => order == BondOrder::Single,
        BondQuery::Double => order == BondOrder::Double,
        BondQuery::Triple => order == BondOrder::Triple,
        BondQuery::Aromatic => order == BondOrder::Aromatic,
        BondQuery::Any => true,
        BondQuery::Ring => {
            if let Some(ri) = ring_info {
                is_ring_bond(from, to, ri)
            } else {
                false
            }
        }
        BondQuery::NotRing => {
            if let Some(ri) = ring_info {
                !is_ring_bond(from, to, ri)
            } else {
                true
            }
        }
        BondQuery::Implicit => {
            // Default bond: single or aromatic
            order == BondOrder::Single || order == BondOrder::Aromatic
        }
        BondQuery::And(parts) => parts
            .iter()
            .all(|q| bond_matches(order, q, _mol, from, to, ring_info)),
        BondQuery::Not(inner) => !bond_matches(order, inner, _mol, from, to, ring_info),
    }
}

/// Match a pattern starting from a specific molecule atom (for recursive SMARTS).
fn substruct_match_from(
    mol: &Molecule,
    pattern: &SmartsPattern,
    start_atom: NodeIndex,
    ring_info: Option<&RingInfo>,
) -> Vec<Vec<usize>> {
    if pattern.atoms.is_empty() {
        return vec![];
    }

    let n_pat = pattern.atoms.len();
    let n_mol = mol.graph.node_count();

    let mut pat_adj: Vec<Vec<(usize, usize)>> = vec![vec![]; n_pat];
    for (bi, bond) in pattern.bonds.iter().enumerate() {
        pat_adj[bond.from].push((bond.to, bi));
        pat_adj[bond.to].push((bond.from, bi));
    }

    let mut results = Vec::new();
    let mut mapping = vec![usize::MAX; n_pat];
    let mut used = vec![false; n_mol];

    // Start from the specified atom
    mapping[0] = start_atom.index();
    used[start_atom.index()] = true;
    if atom_matches(mol, start_atom, &pattern.atoms[0].query, ring_info) {
        backtrack(
            mol,
            pattern,
            &pat_adj,
            &mut mapping,
            &mut used,
            1,
            ring_info,
            &mut results,
        );
    }
    used[start_atom.index()] = false;
    results
}

// ── Ring info helpers ──

/// Pre-computed ring information for a molecule, used to speed up SMARTS matching.
///
/// Compute once with [`precompute_ring_info`] and reuse for many pattern matches.
pub struct RingInfo {
    pub atom_in_ring: Vec<bool>,
    pub atom_ring_sizes: Vec<Vec<usize>>,
    pub rings: Vec<Vec<usize>>,
}

fn compute_ring_info(mol: &Molecule) -> RingInfo {
    let n = mol.graph.node_count();
    let rings = crate::distgeom::find_sssr_pub(mol);

    let mut atom_in_ring = vec![false; n];
    let mut atom_ring_sizes: Vec<Vec<usize>> = vec![vec![]; n];

    for ring in &rings {
        for &a in ring {
            atom_in_ring[a] = true;
            let size = ring.len();
            if !atom_ring_sizes[a].contains(&size) {
                atom_ring_sizes[a].push(size);
            }
        }
    }

    // For macrocycle detection (r{9-}), we also need to detect large rings.
    // The SSSR may not include all rings. For atoms not in SSSR but in a ring,
    // check using BFS shortest alternative path.
    // For now, also detect rings up to size 20 for any bond in the molecule.
    for edge in mol.graph.edge_references() {
        let u = edge.source().index();
        let v = edge.target().index();
        // If already in a detected ring, skip
        if atom_in_ring[u] && atom_in_ring[v] {
            // Check if they share a ring (both in same ring)
            let shared = rings.iter().any(|r| r.contains(&u) && r.contains(&v));
            if shared {
                continue;
            }
        }
        // BFS for alternative path (ring detection)
        if let Some(alt_len) = crate::graph::min_path_excluding2(
            mol,
            NodeIndex::new(u),
            NodeIndex::new(v),
            NodeIndex::new(u),
            NodeIndex::new(v),
            19,
        ) {
            let ring_size = alt_len + 1;
            if ring_size >= 3 {
                atom_in_ring[u] = true;
                atom_in_ring[v] = true;
                if !atom_ring_sizes[u].contains(&ring_size) {
                    atom_ring_sizes[u].push(ring_size);
                }
                if !atom_ring_sizes[v].contains(&ring_size) {
                    atom_ring_sizes[v].push(ring_size);
                }
            }
        }
    }

    RingInfo {
        atom_in_ring,
        atom_ring_sizes,
        rings,
    }
}

fn is_ring_bond(a: NodeIndex, b: NodeIndex, ring_info: &RingInfo) -> bool {
    ring_info.rings.iter().any(|ring| {
        let ai = a.index();
        let bi = b.index();
        if !ring.contains(&ai) || !ring.contains(&bi) {
            return false;
        }
        // Check if a and b are adjacent in the ring
        let len = ring.len();
        for i in 0..len {
            let j = (i + 1) % len;
            if (ring[i] == ai && ring[j] == bi) || (ring[i] == bi && ring[j] == ai) {
                return true;
            }
        }
        false
    })
}

fn count_ring_bonds(mol: &Molecule, atom: NodeIndex, ring_info: &RingInfo) -> usize {
    mol.graph
        .neighbors(atom)
        .filter(|&nb| is_ring_bond(atom, nb, ring_info))
        .count()
}

fn is_aromatic_atom(mol: &Molecule, atom: NodeIndex) -> bool {
    mol.graph
        .edges(atom)
        .any(|e| mol.graph[e.id()].order == BondOrder::Aromatic)
}

fn count_h(mol: &Molecule, atom: NodeIndex) -> usize {
    mol.graph
        .neighbors(atom)
        .filter(|&nb| mol.graph[nb].element == 1)
        .count()
}

fn pattern_needs_ring(pattern: &SmartsPattern) -> bool {
    for atom in &pattern.atoms {
        if query_needs_ring(&atom.query) {
            return true;
        }
    }
    for bond in &pattern.bonds {
        if bond_query_needs_ring(&bond.query) {
            return true;
        }
    }
    false
}

fn query_needs_ring(q: &AtomQuery) -> bool {
    match q {
        AtomQuery::InRing
        | AtomQuery::RingSize(_)
        | AtomQuery::RingSizeRange(..)
        | AtomQuery::RingSizeMin(_)
        | AtomQuery::RingBondCount(_)
        | AtomQuery::RingCount(_) => true,
        AtomQuery::And(parts) | AtomQuery::Or(parts) => parts.iter().any(query_needs_ring),
        AtomQuery::Not(inner) => query_needs_ring(inner),
        AtomQuery::Recursive(inner) => pattern_needs_ring(inner),
        _ => false,
    }
}

fn bond_query_needs_ring(q: &BondQuery) -> bool {
    match q {
        BondQuery::Ring | BondQuery::NotRing => true,
        BondQuery::And(parts) => parts.iter().any(bond_query_needs_ring),
        BondQuery::Not(inner) => bond_query_needs_ring(inner),
        _ => false,
    }
}

/// Batch substructure matching: match a single pattern against many molecules.
/// Returns one Vec<Vec<usize>> per molecule (empty if no match).
pub fn substruct_match_batch(
    molecules: &[&Molecule],
    pattern: &SmartsPattern,
) -> Vec<Vec<Vec<usize>>> {
    molecules
        .iter()
        .map(|mol| substruct_match(mol, pattern))
        .collect()
}

/// Batch substructure matching with rayon parallelism.
#[cfg(feature = "parallel")]
pub fn substruct_match_batch_parallel(
    molecules: &[&Molecule],
    pattern: &SmartsPattern,
) -> Vec<Vec<Vec<usize>>> {
    use rayon::prelude::*;
    molecules
        .par_iter()
        .map(|mol| substruct_match(mol, pattern))
        .collect()
}

/// Check if a molecule contains a substructure (boolean, faster than full match).
pub fn has_substruct_match(mol: &Molecule, pattern: &SmartsPattern) -> bool {
    !substruct_match(mol, pattern).is_empty()
}

/// Batch boolean substructure check with rayon parallelism.
#[cfg(feature = "parallel")]
pub fn has_substruct_match_batch_parallel(
    molecules: &[&Molecule],
    pattern: &SmartsPattern,
) -> Vec<bool> {
    use rayon::prelude::*;
    molecules
        .par_iter()
        .map(|mol| has_substruct_match(mol, pattern))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Molecule;

    #[test]
    fn test_tetrahedral_chirality_matches_explicit_query() {
        let mol = Molecule::from_smiles("C[C@H](F)Cl").unwrap();
        let pattern = parse_smarts("[C@H]").unwrap();
        let inverse = parse_smarts("[C@@H]").unwrap();

        assert!(has_substruct_match(&mol, &pattern));
        assert!(!has_substruct_match(&mol, &inverse));
    }
}
