//! Match CSD torsion patterns against a molecule using SMARTS matching.
//! Equivalent to RDKit's getExperimentalTorsions from TorsionPreferences.cpp.

use crate::forcefield::etkdg_3d::M6TorsionContrib;
use crate::graph::Molecule;
use super::matcher::{substruct_match_with_ring_info, precompute_ring_info};
use super::parser::{parse_smarts, SmartsPattern};
use super::torsion_data::{TORSION_PREFS_MACROCYCLES, TORSION_PREFS_V2};
use std::collections::HashSet;

/// A parsed torsion pattern with its Fourier coefficients.
struct TorsionEntry {
    pattern: SmartsPattern,
    /// Which pattern atom indices correspond to :1, :2, :3, :4
    map_idx: [usize; 4],
    signs: [f64; 6],
    v: [f64; 6],
}

fn parse_table(table: &[(&str, [f32; 12])]) -> Vec<TorsionEntry> {
    let mut entries = Vec::with_capacity(table.len());
    for &(smarts_str, ref coeffs) in table {
        let pattern = match parse_smarts(smarts_str) {
            Ok(p) => p,
            Err(_) => continue,
        };
        let mut map_idx = [usize::MAX; 4];
        for (i, atom) in pattern.atoms.iter().enumerate() {
            if let Some(m) = atom.map_idx {
                if m >= 1 && m <= 4 {
                    map_idx[(m - 1) as usize] = i;
                }
            }
        }
        if map_idx.iter().any(|&x| x == usize::MAX) {
            continue;
        }
        let signs = [
            coeffs[0] as f64, coeffs[2] as f64, coeffs[4] as f64,
            coeffs[6] as f64, coeffs[8] as f64, coeffs[10] as f64,
        ];
        let v = [
            coeffs[1] as f64, coeffs[3] as f64, coeffs[5] as f64,
            coeffs[7] as f64, coeffs[9] as f64, coeffs[11] as f64,
        ];
        entries.push(TorsionEntry { pattern, map_idx, signs, v });
    }
    entries
}

/// Pre-parsed pattern library (parsed once, reused across calls).
struct TorsionLibrary {
    v2: Vec<TorsionEntry>,
    macrocycles: Vec<TorsionEntry>,
}

impl TorsionLibrary {
    fn new() -> Self {
        Self {
            v2: parse_table(TORSION_PREFS_V2),
            macrocycles: parse_table(TORSION_PREFS_MACROCYCLES),
        }
    }
}

// Thread-local cache for parsed patterns (avoids re-parsing on every call).
thread_local! {
    static LIBRARY: TorsionLibrary = TorsionLibrary::new();
}

/// Given a molecule, find all matching CSD experimental torsion contributions.
/// Matches ETKDGv3: v2 patterns + macrocycle patterns (no small ring patterns).
pub fn match_experimental_torsions(mol: &Molecule) -> Vec<M6TorsionContrib> {
    LIBRARY.with(|lib| match_with_library(mol, lib))
}

fn match_with_library(mol: &Molecule, lib: &TorsionLibrary) -> Vec<M6TorsionContrib> {
    let ring_bonds = compute_bond_rings(mol);
    let excluded_bonds = compute_excluded_bonds(mol, &ring_bonds);
    // Pre-compute ring info ONCE for the molecule (SSSR + ring sizes)
    let ring_info = precompute_ring_info(mol);

    let num_bonds = mol.graph.edge_count();
    let mut done_bonds = vec![false; num_bonds];
    let mut torsion_contribs = Vec::new();

    // Process v2 + macrocycles (not small rings, matching ETKDGv3)
    for entry in lib.v2.iter().chain(lib.macrocycles.iter()) {
        // Early exit: if all rotatable bonds already matched, skip remaining patterns
        if done_bonds.iter().all(|&d| d) {
            break;
        }
        // Pre-filter: skip patterns with more atoms than the molecule
        if entry.pattern.atoms.len() > mol.graph.node_count() {
            continue;
        }
        let matches = substruct_match_with_ring_info(mol, &entry.pattern, &ring_info);
        for m in &matches {
            // Extract mapped atoms :1-:4
            let aid1 = m[entry.map_idx[0]];
            let aid2 = m[entry.map_idx[1]];
            let aid3 = m[entry.map_idx[2]];
            let aid4 = m[entry.map_idx[3]];

            // Get the central bond between :2 and :3
            let n2 = petgraph::graph::NodeIndex::new(aid2);
            let n3 = petgraph::graph::NodeIndex::new(aid3);
            let edge = match mol.graph.find_edge(n2, n3) {
                Some(e) => e,
                None => continue,
            };
            let bid = edge.index();

            // Exclusion: bridged/fused ring bonds or bonds in 4+ rings
            if excluded_bonds.contains(&bid) || count_bond_rings(bid, &ring_bonds) > 3 {
                done_bonds[bid] = true;
                continue;
            }

            // First-match-wins
            if done_bonds[bid] {
                continue;
            }
            done_bonds[bid] = true;

            torsion_contribs.push(M6TorsionContrib {
                i: aid1,
                j: aid2,
                k: aid3,
                l: aid4,
                signs: entry.signs,
                v: entry.v,
            });
        }
    }

    torsion_contribs
}

/// Compute bond-ring membership: for each ring, which bond indices are in it.
fn compute_bond_rings(mol: &Molecule) -> Vec<Vec<usize>> {
    let sssr = crate::distgeom::find_sssr_pub(mol);
    let mut bond_rings = Vec::new();
    for ring in &sssr {
        let mut ring_bond_indices = Vec::new();
        let len = ring.len();
        for i in 0..len {
            let a = petgraph::graph::NodeIndex::new(ring[i]);
            let b = petgraph::graph::NodeIndex::new(ring[(i + 1) % len]);
            if let Some(e) = mol.graph.find_edge(a, b) {
                ring_bond_indices.push(e.index());
            }
        }
        bond_rings.push(ring_bond_indices);
    }
    bond_rings
}

/// Count how many SSSR rings contain a given bond.
fn count_bond_rings(bid: usize, ring_bonds: &[Vec<usize>]) -> usize {
    ring_bonds.iter().filter(|ring| ring.contains(&bid)).count()
}

/// Compute excluded bonds: bonds in fused/bridged small ring systems.
/// A bond is excluded if it's in a ring that shares >1 bond with another ring,
/// unless the ring is a macrocycle (size >= 9).
fn compute_excluded_bonds(mol: &Molecule, ring_bonds: &[Vec<usize>]) -> HashSet<usize> {
    const MIN_MACROCYCLE_SIZE: usize = 9;
    let sssr = crate::distgeom::find_sssr_pub(mol);
    let mut excluded = HashSet::new();

    let n_rings = ring_bonds.len();
    for i in 0..n_rings {
        for j in (i + 1)..n_rings {
            // Skip if both are macrocycles
            if sssr[i].len() >= MIN_MACROCYCLE_SIZE && sssr[j].len() >= MIN_MACROCYCLE_SIZE {
                continue;
            }
            // Count shared bonds
            let shared = ring_bonds[i].iter().filter(|b| ring_bonds[j].contains(b)).count();
            if shared > 1 {
                // Exclude all bonds of small rings
                if sssr[i].len() < MIN_MACROCYCLE_SIZE {
                    for &b in &ring_bonds[i] {
                        excluded.insert(b);
                    }
                }
                if sssr[j].len() < MIN_MACROCYCLE_SIZE {
                    for &b in &ring_bonds[j] {
                        excluded.insert(b);
                    }
                }
            }
        }
    }
    excluded
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_torsion_matching_first_mol() {
        let mol = Molecule::from_smiles("C#CCC1COCC2(C)CC(NC23CCOC3)C1=O").unwrap();
        let torsions = match_experimental_torsions(&mol);
        assert_eq!(torsions.len(), 1);
        assert_eq!((torsions[0].i, torsions[0].j, torsions[0].k, torsions[0].l), (1, 2, 3, 4));
    }
}
