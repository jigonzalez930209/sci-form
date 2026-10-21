//! UFF atom typing compatible with the canonical RDKit implementation
//! (`ForceFieldHelpers/UFF/AtomTyper.cpp`).
//!
//! Labels follow `getAtomLabel`: an element symbol plus a hybridization digit
//! (`1` sp, `2` sp2, `3` sp3, `R` conjugated sp2) and, for sulfur/phosphorus,
//! a formal-valence flag (`S_3+6`, `P_3+5`, ...).

use crate::graph::{BondOrder, Hybridization, Molecule};
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;

/// Elements treated as hypervalent for conjugation/charge-flag purposes.
fn is_hypervalent(z: u8) -> bool {
    z == 15 || z == 16
}

/// Lone-pair donor elements that propagate conjugation (RDKit resonance model).
fn is_lone_pair_donor(z: u8) -> bool {
    z == 7 || z == 8
}

pub(crate) fn order_value(order: &BondOrder) -> f64 {
    match order {
        BondOrder::Single => 1.0,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Aromatic => 1.5,
        BondOrder::Unknown => 1.0,
    }
}

/// True when the given order participates in the pi system.
fn is_pi_order(order: &BondOrder) -> bool {
    matches!(
        order,
        BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
    )
}

/// Determine whether a given atom is part of an aromatic ring by checking if any of its
/// bonds have `BondOrder::Aromatic`.
pub fn is_atom_aromatic(mol: &Molecule, idx: usize) -> bool {
    let node = NodeIndex::new(idx);
    mol.graph
        .edges(node)
        .any(|e| e.weight().order == BondOrder::Aromatic)
}

/// Total valence as used by RDKit's UFF typer (sum of bond orders, aromatic = 1.5).
pub fn atom_total_valence(mol: &Molecule, idx: usize) -> f64 {
    let node = NodeIndex::new(idx);
    mol.graph
        .edges(node)
        .map(|e| order_value(&e.weight().order))
        .sum()
}

/// Conjugation flags for all edges (parallel to `mol.graph.edge_indices()`).
///
/// Approximation of `MolOps::setConjugation` validated against RDKit
/// (F17 fixture, `oracle_rdkit.py labels`):
/// - aromatic bonds are always conjugated;
/// - a single bond is conjugated when one endpoint carries a pi bond to a
///   *different* atom and the other endpoint is a lone-pair donor or also
///   carries a pi bond;
/// - a multiple bond is conjugated when an endpoint has a conjugated single
///   bond or two pi bonds (cumulenes, CO2);
/// - pi bonds at hypervalent P/S do not propagate conjugation.
pub fn conjugated_bond_flags(mol: &Molecule) -> Vec<bool> {
    let n = mol.graph.node_count();
    let mut edges: Vec<(usize, usize, BondOrder)> = Vec::new();
    for edge in mol.graph.edge_references() {
        let i = edge.source().index();
        let j = edge.target().index();
        edges.push((i, j, edge.weight().order));
    }
    let mut adjacency: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n]; // (edge idx, other atom)
    for (eidx, &(i, j, _)) in edges.iter().enumerate() {
        adjacency[i].push((eidx, j));
        adjacency[j].push((eidx, i));
    }
    let elements: Vec<u8> = (0..n)
        .map(|i| mol.graph[NodeIndex::new(i)].element)
        .collect();

    let has_pi = |atom: usize| -> bool {
        !is_hypervalent(elements[atom])
            && adjacency[atom]
                .iter()
                .any(|&(e, _)| is_pi_order(&edges[e].2))
    };
    let partner_capable =
        |atom: usize| -> bool { is_lone_pair_donor(elements[atom]) || has_pi(atom) };

    let mut flags = vec![false; edges.len()];
    // Aromatic bonds.
    for (eidx, (_, _, order)) in edges.iter().enumerate() {
        if matches!(order, BondOrder::Aromatic) {
            flags[eidx] = true;
        }
    }
    // Single bonds.
    for (eidx, &(i, j, ref order)) in edges.iter().enumerate() {
        if !matches!(order, BondOrder::Single) {
            continue;
        }
        let pi_i = adjacency[i].iter().any(|&(e, k)| {
            e != eidx && is_pi_order(&edges[e].2) && !is_hypervalent(elements[i]) && k != j
        });
        let pi_j = adjacency[j].iter().any(|&(e, k)| {
            e != eidx && is_pi_order(&edges[e].2) && !is_hypervalent(elements[j]) && k != i
        });
        flags[eidx] = (pi_i && partner_capable(j)) || (pi_j && partner_capable(i));
    }
    // Multiple bonds.
    for (eidx, &(i, j, ref order)) in edges.iter().enumerate() {
        if !matches!(order, BondOrder::Double | BondOrder::Triple) {
            continue;
        }
        let endpoint_has_conj_single = |atom: usize| -> bool {
            adjacency[atom]
                .iter()
                .any(|&(e, _)| e != eidx && flags[e] && matches!(edges[e].2, BondOrder::Single))
        };
        let atom_pi_count = |atom: usize| -> usize {
            if is_hypervalent(elements[atom]) {
                return 0;
            }
            adjacency[atom]
                .iter()
                .filter(|&&(e, _)| e != eidx && is_pi_order(&edges[e].2))
                .count()
        };
        flags[eidx] = endpoint_has_conj_single(i)
            || endpoint_has_conj_single(j)
            || atom_pi_count(i) >= 1
            || atom_pi_count(j) >= 1;
    }
    flags
}

/// True when any incident bond of `idx` is conjugated.
pub fn atom_has_conjugated_bond(mol: &Molecule, idx: usize) -> bool {
    let node = NodeIndex::new(idx);
    let flags = conjugated_bond_flags(mol);
    mol.graph
        .edge_references()
        .any(|e| flags[e.id().index()] && (e.source() == node || e.target() == node))
}

/// Assign a UFF atom type label following RDKit's `getAtomLabel`.
pub fn assign_uff_type(
    element: u8,
    hyb: &Hybridization,
    is_aromatic: bool,
    has_conjugated_bond: bool,
    total_valence: f64,
) -> &'static str {
    match element {
        1 => "H_",
        5 => match hyb {
            Hybridization::SP2 => "B_2",
            _ => "B_3",
        },
        6 => match hyb {
            Hybridization::SP if !is_aromatic => "C_1",
            Hybridization::SP2 => {
                if is_aromatic || has_conjugated_bond {
                    "C_R"
                } else {
                    "C_2"
                }
            }
            _ => "C_3",
        },
        7 => match hyb {
            Hybridization::SP if !is_aromatic => "N_1",
            Hybridization::SP2 => {
                if is_aromatic || has_conjugated_bond {
                    "N_R"
                } else {
                    "N_2"
                }
            }
            _ => "N_3",
        },
        8 => match hyb {
            Hybridization::SP if !is_aromatic => "O_1",
            Hybridization::SP2 => {
                if is_aromatic || has_conjugated_bond {
                    "O_R"
                } else {
                    "O_2"
                }
            }
            _ => "O_3",
        },
        9 => "F_",
        14 => "Si3",
        15 => match total_valence.round() as i32 {
            3 => "P_3+3",
            5 => "P_3+5",
            _ => "P_3",
        },
        16 => match hyb {
            Hybridization::SP2 => {
                if is_aromatic || has_conjugated_bond {
                    "S_R"
                } else {
                    "S_2"
                }
            }
            _ => match total_valence.round() as i32 {
                6 => "S_3+6",
                4 => "S_3+4",
                2 => "S_3+2",
                _ => "S_3",
            },
        },
        17 => "Cl",
        // Transition metals
        22 => "Ti3+4",
        23 => "V_3+5",
        24 => "Cr6+3",
        25 => "Mn6+2",
        26 => "Fe3+2",
        27 => "Co6+3",
        28 => "Ni4+2",
        29 => "Cu3+1",
        30 => "Zn3+2",
        // Additional main-group elements
        32 => "Ge3",
        33 => "As3",
        34 => "Se3",
        35 => "Br",
        42 => "Mo3+6",
        46 => "Pd4+2",
        47 => "Ag1+1",
        50 => "Sn3",
        51 => "Sb3",
        52 => "Te3",
        53 => "I_",
        78 => "Pt4+2",
        79 => "Au4+3",
        _ => "C_3",
    }
}
