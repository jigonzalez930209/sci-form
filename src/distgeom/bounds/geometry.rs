use crate::graph::Molecule;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;

#[derive(PartialEq)]
pub(super) enum AmidePref {
    Cis,
    Trans,
    None,
}

/// Detect amide/ester14 preference: the C=O double bond IS in the 1-4 path.
/// Path: i-n1-n2-j where n2-j is double bond to O/N (or i-n1 is double bond).
/// RDKit's _checkAmideEster14:
///   Pattern: 1—2—3=4 where 3=C, 4=O/N (double bond), 2=N/O (amide/ester)
///     bnd3 is double, a4 is O/N/S, a3 is C, bnd1 is single, a2 is N/O
///   - Secondary amide H on N → TRANS to the O
///   - Others → CIS to the O
pub(super) fn detect_amide_ester14(
    mol: &Molecule,
    i: NodeIndex,
    n1: NodeIndex,
    n2: NodeIndex,
    j: NodeIndex,
) -> AmidePref {
    // Try direction: n2=C(=O), j=O/N (double bond), n1=N/O
    if let Some(pref) = check_amide14_dir(mol, i, n1, n2, j) {
        return pref;
    }
    // Try reverse: n1=C(=O), i=O/N (double bond), n2=N/O
    if let Some(pref) = check_amide14_dir(mol, j, n2, n1, i) {
        return pref;
    }
    AmidePref::None
}

pub(super) fn check_amide14_dir(
    mol: &Molecule,
    outer_atom: NodeIndex, // atom 1 (substituent on amide N)
    amide_n: NodeIndex,    // atom 2 (N or O)
    carbonyl_c: NodeIndex, // atom 3 (C)
    double_o: NodeIndex,   // atom 4 (O/N bonded by double bond to C)
) -> Option<AmidePref> {
    let e_outer = mol.graph[outer_atom].element;
    let e_c = mol.graph[carbonyl_c].element;
    let e_n = mol.graph[amide_n].element;
    let e_o = mol.graph[double_o].element;

    // carbonyl_c must be C(6)
    if e_c != 6 {
        return Option::None;
    }
    // double_o must be O(8) or N(7)
    if e_o != 8 && e_o != 7 {
        return Option::None;
    }
    // amide_n must be N(7) with secondary amide (exactly 1 H neighbor) or O(8)
    if e_n == 7 {
        // Count H neighbors on N
        let h_count = mol
            .graph
            .neighbors(amide_n)
            .filter(|&nb| mol.graph[nb].element == 1)
            .count();
        if h_count != 1 {
            return Option::None;
        }
    } else if e_n != 8 {
        return Option::None;
    }

    // Check that the bond carbonyl_c—double_o is a double bond
    if let Some(edge) = mol.graph.find_edge(carbonyl_c, double_o) {
        if mol.graph[edge].order != crate::graph::BondOrder::Double {
            return Option::None;
        }
    } else {
        return Option::None;
    }

    // Check that bond outer_atom—amide_n is single
    if let Some(edge) = mol.graph.find_edge(outer_atom, amide_n) {
        if mol.graph[edge].order != crate::graph::BondOrder::Single {
            return Option::None;
        }
    } else {
        return Option::None;
    }

    // RDKit CIS/TRANS: H on secondary amide N → TRANS, everything else → CIS
    if e_outer == 1 && e_n == 7 {
        Some(AmidePref::Trans)
    } else {
        Some(AmidePref::Cis)
    }
}

/// Detect amide/ester preference for a 1-4 chain path i-n1-n2-j
/// where n1-n2 is a sp2-sp2 single bond.
/// Matches RDKit's forceTransAmides logic:
/// - If n1 is C bonded to O/N by double bond AND n2 is N or O → amide/ester
/// - For amide N-H, H is TRANS to the carbonyl partner
/// - For other substituents, they are CIS to the carbonyl partner
pub(super) fn detect_amide_ester_preference(
    mol: &Molecule,
    i: NodeIndex,
    n1: NodeIndex,
    n2: NodeIndex,
    j: NodeIndex,
) -> AmidePref {
    // Check pattern: n1 is C(=O) or C(=N), and n2 is N or O
    let e1 = mol.graph[n1].element;
    let e2 = mol.graph[n2].element;

    // Try both directions: n1=C(=X), n2=N/O or n2=C(=X), n1=N/O
    if let Some(pref) = check_amide_dir(mol, i, n1, n2, j, e1, e2) {
        return pref;
    }
    if let Some(pref) = check_amide_dir(mol, j, n2, n1, i, e2, e1) {
        return pref;
    }
    AmidePref::None
}

pub(super) fn check_amide_dir(
    mol: &Molecule,
    atom_near_carbonyl: NodeIndex,
    carbonyl_c: NodeIndex,
    amide_n: NodeIndex,
    atom_on_n: NodeIndex,
    e_c: u8,
    e_n: u8,
) -> Option<AmidePref> {
    // carbonyl_c must be C(6) and amide_n must be N(7) with 1 H or O(8)
    // RDKit's _checkAmideEster15: (a2Num == 8) || (a2Num == 7 && totalNumHs == 1)
    if e_c != 6 {
        return Option::None;
    }
    if e_n == 8 {
        // OK — ester oxygen
    } else if e_n == 7 {
        // Secondary amide: N must have exactly 1 H neighbor
        let h_count = mol
            .graph
            .neighbors(amide_n)
            .filter(|&nb| mol.graph[nb].element == 1)
            .count();
        if h_count != 1 {
            return Option::None;
        }
    } else {
        return Option::None;
    }

    // Check if carbonyl_c has a double bond to O or N (the carbonyl partner)
    let mut has_carbonyl = false;
    for nb in mol.graph.neighbors(carbonyl_c) {
        if nb == amide_n || nb == atom_near_carbonyl {
            continue;
        }
        if let Some(edge) = mol.graph.find_edge(carbonyl_c, nb) {
            let bond = &mol.graph[edge];
            if bond.order == crate::graph::BondOrder::Double {
                let nb_elem = mol.graph[nb].element;
                if nb_elem == 8 || nb_elem == 7 || nb_elem == 16 {
                    has_carbonyl = true;
                    break;
                }
            }
        }
    }

    if !has_carbonyl {
        return Option::None;
    }

    // It's an amide/ester pattern! (_checkAmideEster15: separate C=O branch)
    // RDKit: for this pattern, H on secondary amide N → CIS, others → TRANS
    if mol.graph[atom_on_n].element == 1 && e_n == 7 {
        Some(AmidePref::Cis)
    } else {
        Some(AmidePref::Trans)
    }
}

/// UFF atomic radius by (element, hybridization, conjugated) — from Rappé et al. 1992
/// When conjugated=true and SP2, uses _R (aromatic/conjugated) params instead of _2.
/// This matches RDKit's UFF atom typing: SP2 atoms that are aromatic or have a
/// conjugated bond get the _R type (C_R, N_R, O_R, S_R).
fn uff_radius(element: u8, hyb: &crate::graph::Hybridization, conjugated: bool) -> f64 {
    use crate::graph::Hybridization::*;
    match (element, hyb) {
        (1, _) => 0.354,
        (6, SP) => 0.706,
        (6, SP2) => {
            if conjugated {
                0.729
            } else {
                0.732
            }
        } // C_R vs C_2
        (6, SP3) | (6, _) => 0.757,
        (7, SP) => 0.656,
        (7, SP2) => {
            if conjugated {
                0.699
            } else {
                0.685
            }
        } // N_R vs N_2
        (7, SP3) | (7, _) => 0.700,
        (8, SP) => 0.639, // O_1
        (8, SP2) => {
            if conjugated {
                0.680
            } else {
                0.634
            }
        } // O_R vs O_2
        (8, SP3) | (8, _) => 0.658,
        (9, _) => 0.668,
        (15, _) => 1.075,
        (16, SP2) => {
            if conjugated {
                1.077
            } else {
                0.854
            }
        } // S_R vs S_2
        (16, _) => 1.047,
        (17, _) => 1.044,
        (35, _) => 1.218,
        (53, _) => 1.390,
        _ => crate::graph::get_covalent_radius(element),
    }
}

/// Check if an atom has a conjugated bond (matching RDKit's atomHasConjugatedBond).
/// A bond between two SP2/SP atoms is conjugated if at least one of them has
/// another SP2/SP non-H neighbor — forming a pi conjugation path.
/// An isolated double bond (like C=O in an aldehyde where C's other neighbors
/// are SP3/H) is NOT conjugated.
fn atom_is_conjugated(mol: &Molecule, node: NodeIndex) -> bool {
    use crate::graph::Hybridization::{SP, SP2};
    let atom = &mol.graph[node];
    if atom.hybridization != SP2 && atom.hybridization != SP {
        return false;
    }
    for edge in mol.graph.edges(node) {
        let nb = edge.target();
        let nb_atom = &mol.graph[nb];
        // Both atoms must be SP2 or SP for the bond to be conjugated
        if nb_atom.hybridization != SP2 && nb_atom.hybridization != SP {
            continue;
        }
        // Check if either atom has another SP2/SP non-H neighbor (besides each other)
        let self_has_sp2_nb = mol.graph.neighbors(node).any(|n| {
            n != nb && mol.graph[n].element != 1 && matches!(mol.graph[n].hybridization, SP2 | SP)
        });
        let nb_has_sp2_nb = mol.graph.neighbors(nb).any(|n| {
            n != node && mol.graph[n].element != 1 && matches!(mol.graph[n].hybridization, SP2 | SP)
        });
        if self_has_sp2_nb || nb_has_sp2_nb {
            return true;
        }
    }
    false
}

/// UFF GMP electronegativity (GMP_Xi) by element — from Rappé et al. 1992
fn uff_chi(element: u8) -> f64 {
    match element {
        1 => 4.528,
        5 => 5.110,
        6 => 5.343,
        7 => 6.899,
        8 => 8.741,
        9 => 10.874,
        14 => 4.168,
        15 => 5.463,
        16 => 6.928,
        17 => 8.564,
        35 => 7.790,
        53 => 6.822,
        _ => 5.0, // fallback
    }
}

pub fn get_bond_length(mol: &Molecule, n1: NodeIndex, n2: NodeIndex) -> f64 {
    let a1 = &mol.graph[n1];
    let a2 = &mol.graph[n2];
    let mut e1 = a1.element;
    let mut e2 = a2.element;
    // Detect conjugation for each atom
    let conj1 = atom_is_conjugated(mol, n1);
    let conj2 = atom_is_conjugated(mol, n2);
    let (hyb1, hyb2, c1, c2) = if e1 <= e2 {
        (&a1.hybridization, &a2.hybridization, conj1, conj2)
    } else {
        std::mem::swap(&mut e1, &mut e2);
        (&a2.hybridization, &a1.hybridization, conj2, conj1)
    };

    let mut order = crate::graph::BondOrder::Single;
    if let Some(edge_idx) = mol.graph.find_edge(n1, n2) {
        order = mol.graph[edge_idx].order;
    }

    // Pure UFF formula matching RDKit's set12Bounds:
    // bondLength = ri + rj + rBO - rEN
    // Uses _R radius for conjugated SP2 atoms (C, N, O, S)
    let r1 = uff_radius(e1, hyb1, c1);
    let r2 = uff_radius(e2, hyb2, c2);
    let bo: f64 = match order {
        crate::graph::BondOrder::Single | crate::graph::BondOrder::Unknown => 1.0,
        crate::graph::BondOrder::Double => 2.0,
        crate::graph::BondOrder::Triple => 3.0,
        crate::graph::BondOrder::Aromatic => 1.5,
    };
    let r_bo = if bo > 1.0 {
        -0.1332 * (r1 + r2) * bo.ln()
    } else {
        0.0
    };
    let xi = uff_chi(e1);
    let xj = uff_chi(e2);
    let denom = xi * r1 + xj * r2;
    let r_en = if denom > 1e-6 {
        r1 * r2 * (xi.sqrt() - xj.sqrt()).powi(2) / denom
    } else {
        0.0
    };
    r1 + r2 + r_bo - r_en
}

pub(super) fn compute_14_dist_cis(r12: f64, r23: f64, r34: f64, a123: f64, a234: f64) -> f64 {
    let x0 = r12 * a123.cos();
    let y0 = r12 * a123.sin();
    let x3 = r23 - r34 * a234.cos();
    let y3 = r34 * a234.sin();
    ((x3 - x0).powi(2) + (y3 - y0).powi(2)).sqrt()
}

pub(super) fn compute_14_dist_trans(r12: f64, r23: f64, r34: f64, a123: f64, a234: f64) -> f64 {
    let x0 = r12 * a123.cos();
    let y0 = r12 * a123.sin();
    let x3 = r23 - r34 * a234.cos();
    let y3 = -r34 * a234.sin();
    ((x3 - x0).powi(2) + (y3 - y0).powi(2)).sqrt()
}

#[derive(Clone, Copy, PartialEq)]
pub(super) enum Path14Type {
    Cis,
    Trans,
    Other,
}

pub(super) fn compute13_dist(d1: f64, d2: f64, angle: f64) -> f64 {
    (d1 * d1 + d2 * d2 - 2.0 * d1 * d2 * angle.cos()).sqrt()
}

pub(super) fn compute15_cis_cis(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() - d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / d14).clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 - ang143;
    compute13_dist(d14, d4, ang145)
}

pub(super) fn compute15_cis_trans(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() - d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / d14).clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 + ang143;
    compute13_dist(d14, d4, ang145)
}

pub(super) fn compute15_trans_trans(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() + d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * ang23.cos() + d1 * (ang12 - ang23).cos()) / d14).clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 + ang143;
    compute13_dist(d14, d4, ang145)
}

pub(super) fn compute15_trans_cis(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() + d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * ang23.cos() + d1 * (ang12 - ang23).cos()) / d14).clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 - ang143;
    compute13_dist(d14, d4, ang145)
}
