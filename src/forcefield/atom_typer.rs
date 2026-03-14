use crate::graph::Hybridization;

/// Assign a UFF atom type label to an atom given its element, hybridization, and aromaticity.
///
/// Rules mirror RDKit's UFF atom-typing (rdForceFieldHelpers.cpp).
/// Atoms that cannot be typed fall back to "C_3" so the force field always has parameters.
pub fn assign_uff_type(element: u8, hyb: &Hybridization, is_aromatic: bool) -> &'static str {
    match element {
        1 => "H_",
        5 => "B_2",
        6 => {
            if is_aromatic {
                "C_R"
            } else {
                match hyb {
                    Hybridization::SP => "C_1",
                    Hybridization::SP2 => "C_2",
                    _ => "C_3", // SP3 and Unknown
                }
            }
        }
        7 => {
            if is_aromatic {
                "N_R"
            } else {
                match hyb {
                    Hybridization::SP => "N_1",
                    Hybridization::SP2 => "N_2",
                    _ => "N_3",
                }
            }
        }
        8 => {
            if is_aromatic {
                "O_R"
            } else {
                match hyb {
                    Hybridization::SP => "O_1",
                    Hybridization::SP2 => "O_2",
                    _ => "O_3",
                }
            }
        }
        9 => "F_",
        14 => "Si3",
        15 => "P_3",
        16 => {
            if is_aromatic {
                "S_R"
            } else {
                match hyb {
                    Hybridization::SP2 => "S_2",
                    _ => "S_3",
                }
            }
        }
        17 => "Cl",
        35 => "Br",
        53 => "I_",
        _ => "C_3",
    }
}

/// Determine whether a given atom is part of an aromatic ring by checking if any of its
/// bonds have `BondOrder::Aromatic`.
pub fn is_atom_aromatic(mol: &crate::graph::Molecule, idx: usize) -> bool {
    use petgraph::graph::NodeIndex;
    let node = NodeIndex::new(idx);
    for edge in mol.graph.edges(node) {
        if edge.weight().order == crate::graph::BondOrder::Aromatic {
            return true;
        }
    }
    false
}
