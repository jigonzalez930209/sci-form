//! Aromaticity perception for Kekulé SMILES input.
//!
//! sci-form keeps the Kekulé bond orders produced by the SMILES parser;
//! RDKit's sanitization instead re-perceives aromaticity (Hückel's 4n+2 rule)
//! and rewrites the ring bonds. Consumers that must match RDKit — the UFF
//! force field above all — call [`perceive_aromaticity`] (or
//! [`perceived_aromaticity`]) before typing, which marks the ring bonds of
//! Hückel-aromatic rings as [`BondOrder::Aromatic`] and promotes the ring
//! atoms to sp², mirroring RDKit's `MolOps::setAromaticity`.
//!
//! The electron bookkeeping follows the RDKit default model for the elements
//! the validated corpus exercises (B, C, N, O, P, S):
//!
//! * carbon contributes 1 π electron when it carries a ring double bond,
//!   0 with an exocyclic double bond (e.g. the C2 of oxazole), 2 as a
//!   carbanion and 0 as a carbocation;
//! * nitrogen/phosphorus contribute 1 with a ring double bond and 2 when
//!   they only donate a lone pair (pyrrole-type), and never when they are
//!   four-coordinate;
//! * oxygen/sulfur contribute 2 (furan/thiophene-type);
//! * boron contributes 0 (empty p orbital).
//!
//! Rings whose bonds are already aromatic (lowercase SMILES input) are left
//! untouched.

use super::sssr::compute_sssr;
use crate::graph::{BondOrder, Hybridization, Molecule};
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;

/// Mark the Hückel-aromatic rings of `mol` in place.
///
/// Ring bonds become [`BondOrder::Aromatic`] and their atoms sp². Molecules
/// parsed from lowercase aromatic SMILES are unchanged.
pub fn perceive_aromaticity(mol: &mut Molecule) {
    let rings = aromatic_rings(mol);
    for ring in rings {
        mark_ring_aromatic(mol, &ring);
    }
}

/// Return an aromatized copy of `mol` when Kekulé aromatic rings are detected.
///
/// `None` means the molecule already carries aromatic bonds (or has no
/// aromatic ring), so callers can keep using the original.
pub fn perceived_aromaticity(mol: &Molecule) -> Option<Molecule> {
    let rings = aromatic_rings(mol);
    if rings.is_empty() {
        return None;
    }
    let mut out = mol.clone();
    for ring in rings {
        mark_ring_aromatic(&mut out, &ring);
    }
    Some(out)
}

/// Atom rings (in SSSR order) that satisfy the Hückel rule.
fn aromatic_rings(mol: &Molecule) -> Vec<Vec<usize>> {
    let sssr = compute_sssr(mol);
    sssr.rings
        .iter()
        .filter(|ring| ring_aromatic(mol, &ring.atoms))
        .map(|ring| ring.atoms.clone())
        .collect()
}

/// True when the ring is aromatic and not already marked as such.
fn ring_aromatic(mol: &Molecule, ring: &[usize]) -> bool {
    // Already-perceived rings (lowercase SMILES) need no rewrite.
    if ring.iter().any(|&a| has_aromatic_bond(mol, a)) {
        return false;
    }

    let mut pi_electrons = 0i32;
    for &atom_idx in ring {
        let node = NodeIndex::new(atom_idx);
        let atom = &mol.graph[node];
        let degree = mol.graph.edges(node).count();

        let mut ring_double = false;
        let mut exocyclic_pi = false;
        for edge in mol.graph.edges(node) {
            let other = if edge.source() == node {
                edge.target().index()
            } else {
                edge.source().index()
            };
            let is_pi = matches!(edge.weight().order, BondOrder::Double | BondOrder::Triple);
            if !is_pi {
                continue;
            }
            if ring.contains(&other) {
                ring_double = true;
            } else {
                exocyclic_pi = true;
            }
        }

        let contribution = match atom.element {
            5 => 0, // boron: empty p orbital
            6 => {
                if atom.formal_charge < 0 {
                    2 // carbanion (cyclopentadienyl-type)
                } else if atom.formal_charge > 0 || exocyclic_pi {
                    // Carbocation, or the π system lies outside the ring
                    // (e.g. the C2 of oxazole, `N=C1OC=CO1`).
                    0
                } else if ring_double {
                    1
                } else {
                    return false; // saturated ring carbon
                }
            }
            7 | 15 => {
                if degree >= 4 {
                    return false; // four-coordinate N/P cannot donate
                }
                if ring_double {
                    1 // pyridine-type
                } else {
                    2 // pyrrole-type lone pair
                }
            }
            8 | 16 => {
                if degree >= 3 {
                    return false;
                }
                2 // furan/thiophene-type lone pair
            }
            _ => return false,
        };
        pi_electrons += contribution;
    }

    pi_electrons >= 2 && (pi_electrons - 2) % 4 == 0
}

fn has_aromatic_bond(mol: &Molecule, atom_idx: usize) -> bool {
    mol.graph
        .edges(NodeIndex::new(atom_idx))
        .any(|e| e.weight().order == BondOrder::Aromatic)
}

/// Mark the ring bonds aromatic and promote the ring atoms to sp².
fn mark_ring_aromatic(mol: &mut Molecule, ring: &[usize]) {
    let edges: Vec<_> = mol
        .graph
        .edge_references()
        .filter(|e| ring.contains(&e.source().index()) && ring.contains(&e.target().index()))
        .map(|e| e.id())
        .collect();
    for edge in edges {
        mol.graph[edge].order = BondOrder::Aromatic;
    }
    for &atom_idx in ring {
        let node = NodeIndex::new(atom_idx);
        if mol.graph[node].hybridization != Hybridization::SP {
            mol.graph[node].hybridization = Hybridization::SP2;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Molecule;

    fn perceived(smiles: &str) -> Molecule {
        let mut mol = Molecule::from_smiles(smiles).unwrap();
        perceive_aromaticity(&mut mol);
        mol
    }

    fn aromatic_bonds(mol: &Molecule) -> usize {
        mol.graph
            .edge_references()
            .filter(|e| e.weight().order == BondOrder::Aromatic)
            .count()
    }

    #[test]
    fn kekule_pyrrole_is_perceived() {
        assert_eq!(aromatic_bonds(&perceived("C1=CC=CN1")), 5);
    }

    #[test]
    fn kekule_furan_is_perceived() {
        assert_eq!(aromatic_bonds(&perceived("C1=COC=C1")), 5);
    }

    #[test]
    fn kekule_benzene_is_perceived() {
        assert_eq!(aromatic_bonds(&perceived("C1=CC=CC=C1")), 6);
    }

    #[test]
    fn saturated_rings_are_not_perceived() {
        assert_eq!(aromatic_bonds(&perceived("C1CCCCC1")), 0);
        assert_eq!(aromatic_bonds(&perceived("C1=CCCC1")), 0);
    }

    #[test]
    fn cyclobutadiene_is_not_perceived() {
        assert_eq!(aromatic_bonds(&perceived("C1=CC=C1")), 0);
    }

    #[test]
    fn lowercase_input_is_untouched() {
        let mol = perceived("c1ccccc1");
        assert_eq!(aromatic_bonds(&mol), 6);
        assert!(perceived_aromaticity(&mol).is_none());
    }
}
