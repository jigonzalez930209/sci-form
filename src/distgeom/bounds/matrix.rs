use super::geometry::*;
use crate::graph::Molecule;
use nalgebra::DMatrix;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use std::collections::HashSet;

/// Compute the distance-bounds matrix for a molecule using default settings.
///
/// Equivalent to `calculate_bounds_matrix_opts(mol, true)`.
pub fn calculate_bounds_matrix(mol: &Molecule) -> DMatrix<f64> {
    calculate_bounds_matrix_opts(mol, true)
}

/// Compute the distance-bounds matrix for a molecule.
///
/// The matrix `B` encodes lower and upper distance constraints between every pair
/// of atoms. It is populated in four passes:
///
/// 1. **1-2 bounds** — bond lengths from UFF covalent radii.
/// 2. **1-3 bounds** — from bond lengths and ideal angles (law of cosines).
/// 3. **1-4 bounds** — from torsion angle extremes; Set-1 improved bounds if `use_set15 = true`.
/// 4. **VDW bounds** — lower bounds from Van der Waals radii for all other pairs.
///
/// `use_set15` enables improved 1-4 bounds that better match RDKit's ETKDGv2.
pub fn calculate_bounds_matrix_opts(mol: &Molecule, use_set15: bool) -> DMatrix<f64> {
    let n = mol.graph.node_count();
    let mut bounds = DMatrix::from_element(n, n, 0.0f64);
    for i in 0..n {
        for j in 0..n {
            if i < j {
                bounds[(i, j)] = 1000.0;
            }
        }
    }
    let mut top_dist = DMatrix::from_element(n, n, 1000);
    for i in 0..n {
        top_dist[(i, i)] = 0;
    }
    for bond in mol.graph.edge_references() {
        let (i, j) = (bond.source().index(), bond.target().index());
        top_dist[(i, j)] = 1;
        top_dist[(j, i)] = 1;
    }
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                if top_dist[(i, j)] > top_dist[(i, k)] + top_dist[(k, j)] {
                    top_dist[(i, j)] = top_dist[(i, k)] + top_dist[(k, j)];
                }
            }
        }
    }
    fn set_b(b: &mut DMatrix<f64>, i: usize, j: usize, l: f64, u: f64) {
        let (mi, ma) = if i < j { (i, j) } else { (j, i) };
        if b[(ma, mi)] < l {
            b[(ma, mi)] = l;
        }
        if b[(mi, ma)] > u {
            b[(mi, ma)] = u;
        }
    }
    // RDKit-style _checkAndSetBounds: first time → set; subsequent → widen (union)
    fn set_b_union(b: &mut DMatrix<f64>, i: usize, j: usize, l: f64, u: f64) {
        let (mi, ma) = if i < j { (i, j) } else { (j, i) };
        let clb = b[(ma, mi)];
        let cub = b[(mi, ma)];
        // Lower bound: if no prior bounds (≤0.01), set directly; otherwise take smaller
        if (clb <= 0.01) || (l < clb && l > 0.01) {
            b[(ma, mi)] = l;
        }
        // Upper bound: if no prior bounds (≥999), set directly; otherwise take larger
        if (cub >= 999.0) || (u > cub && u < 999.0) {
            b[(mi, ma)] = u;
        }
    }
    // === squishAtoms: flag atoms adjacent to larger heteroatoms in conjugated 5-rings ===
    // RDKit adds extra flexibility for bonds involving these atoms
    let mut squish_atoms = vec![false; n];
    for bond in mol.graph.edge_references() {
        let (si, ti) = (bond.source(), bond.target());
        let sa = &mol.graph[si];
        let ta = &mol.graph[ti];
        // RDKit checks bond->getIsConjugated() && (atomicNum > 10) && inRingOfSize(5)
        let bond_is_conjugated = matches!(
            sa.hybridization,
            crate::graph::Hybridization::SP2 | crate::graph::Hybridization::SP
        ) && matches!(
            ta.hybridization,
            crate::graph::Hybridization::SP2 | crate::graph::Hybridization::SP
        );
        if bond_is_conjugated
            && (sa.element > 10 || ta.element > 10)
            && crate::graph::bond_in_ring_of_size(mol, si, ti, 5)
        {
            squish_atoms[si.index()] = true;
            squish_atoms[ti.index()] = true;
        }
    }
    for bond in mol.graph.edge_references() {
        let d = get_bond_length(mol, bond.source(), bond.target());
        let extra_squish =
            if squish_atoms[bond.source().index()] || squish_atoms[bond.target().index()] {
                0.2
            } else {
                0.0
            };
        set_b(
            &mut bounds,
            bond.source().index(),
            bond.target().index(),
            d - extra_squish - 0.01,
            d + extra_squish + 0.01,
        );
    }
    for i in 0..n {
        let ni = NodeIndex::new(i);
        let nbs: Vec<_> = mol.graph.neighbors(ni).collect();
        for j in 0..nbs.len() {
            for k in (j + 1)..nbs.len() {
                let (n1, n2) = (nbs[j], nbs[k]);

                // RDKit uses DIST13_TOL = 0.04, doubled for each isLargerSP2Atom
                let mut tol: f64 = 0.04;
                // isLargerSP2Atom: atomicNum > 13, SP2, in a ring
                for &atom_idx in &[n1, ni, n2] {
                    let atom = &mol.graph[atom_idx];
                    if atom.element > 13
                        && atom.hybridization == crate::graph::Hybridization::SP2
                        && crate::graph::atom_in_ring(mol, atom_idx)
                    {
                        tol *= 2.0;
                    }
                }

                // Always compute angle-based 1,3 distance, even for 3-ring
                // atoms where n1-n2 are bonded. RDKit computes the distance
                // through the center using the ring angle (60° for 3-ring).
                let a = crate::graph::get_corrected_ideal_angle(mol, ni, n1, n2);
                let (d1, d2) = (get_bond_length(mol, ni, n1), get_bond_length(mol, ni, n2));
                let di = (d1 * d1 + d2 * d2 - 2.0 * d1 * d2 * a.cos()).sqrt();
                let l_val = di - tol;
                let u_val = di + tol;

                // Use union (widen) semantics matching RDKit's _checkAndSetBounds.
                // When a pair has multiple 1,3 paths (through different centers)
                // or is also a 1,2 pair (3-ring), we take the widest range.
                set_b_union(&mut bounds, n1.index(), n2.index(), l_val, u_val);
                top_dist[(n1.index(), n2.index())] = 2;
                top_dist[(n2.index(), n1.index())] = 2;
            }
        }
    }
    // === set14Bounds ===
    // Collect path14 info for set15Bounds
    let mut path14_list: Vec<(usize, usize, usize, usize, Path14Type)> = Vec::new();
    let mut cis_keys: HashSet<(usize, usize, usize, usize)> = HashSet::new();
    let mut trans_keys: HashSet<(usize, usize, usize, usize)> = HashSet::new();
    for i in 0..n {
        for j in (i + 1)..n {
            if top_dist[(i, j)] == 3 {
                let ni = NodeIndex::new(i);
                let nj = NodeIndex::new(j);
                for n1 in mol.graph.neighbors(ni) {
                    for n2 in mol.graph.neighbors(nj) {
                        if let Some(_e) = mol.graph.find_edge(n1, n2) {
                            let (r1, r2, r3) = (
                                get_bond_length(mol, ni, n1),
                                get_bond_length(mol, n1, n2),
                                get_bond_length(mol, n2, nj),
                            );
                            let (a1, a2) = (
                                crate::graph::get_corrected_ideal_angle(mol, n1, ni, n2),
                                crate::graph::get_corrected_ideal_angle(mol, n2, n1, nj),
                            );
                            let (dc, dt) = (
                                compute_14_dist_cis(r1, r2, r3, a1, a2),
                                compute_14_dist_trans(r1, r2, r3, a1, a2),
                            );

                            // RDKit: skip bounds for 1-4 paths in small rings (≤5 members).
                            // Only record the path type for set15Bounds.
                            // Ring size = 3 (path ni→n1→n2→nj) + alt_path_length.
                            // For ring ≤5: alt_path ≤ 2.
                            let in_small_ring =
                                crate::graph::min_path_excluding2(mol, ni, nj, n1, n2, 2).is_some();

                            let hyb_n1 = mol.graph[n1].hybridization;
                            let hyb_n2 = mol.graph[n2].hybridization;
                            let both_sp2 = hyb_n1 == crate::graph::Hybridization::SP2
                                && hyb_n2 == crate::graph::Hybridization::SP2;

                            let stereo = mol.graph[_e].stereo;
                            let central_order = mol.graph[_e].order;

                            let (l, u, path_type) = if stereo == crate::graph::BondStereo::Z {
                                (dc - 0.06, dc + 0.06, Path14Type::Cis)
                            } else if stereo == crate::graph::BondStereo::E {
                                (dt - 0.06, dt + 0.06, Path14Type::Trans)
                            } else if both_sp2 {
                                let ring_back =
                                    crate::graph::min_path_excluding2(mol, ni, nj, n1, n2, 5);
                                if ring_back.is_some() {
                                    (dc - 0.06, dc + 0.06, Path14Type::Cis)
                                } else {
                                    // Check partial ring: n1,n2 share a ring (don't exclude
                                    // ni/nj since ring path may go through them)
                                    let shared_ring =
                                        crate::graph::min_path_excluding2(mol, n1, n2, n1, n1, 6);
                                    if shared_ring.is_some() {
                                        // n1-n2 is a ring bond. Determine if outer atoms are in ring.
                                        // RDKit logic:
                                        //  - Both external (_setShareRingBond14Bounds) → CIS
                                        //  - One external, one in ring (_setTwoInSameRing14Bounds) → TRANS
                                        let ni_in_ring = crate::graph::min_path_excluding2(
                                            mol, ni, n2, n1, n1, 7,
                                        )
                                        .is_some();
                                        let nj_in_ring = crate::graph::min_path_excluding2(
                                            mol, nj, n1, n2, n2, 7,
                                        )
                                        .is_some();
                                        if ni_in_ring || nj_in_ring {
                                            // External substituent is trans to ring atom
                                            (dt - 0.06, dt + 0.06, Path14Type::Trans)
                                        } else {
                                            // Both external: same side → cis
                                            (dc - 0.06, dc + 0.06, Path14Type::Cis)
                                        }
                                    } else if central_order == crate::graph::BondOrder::Double {
                                        // RDKit: for double bonds without stereo, [cis, trans] range
                                        // Only add tolerance if cis ≈ trans
                                        let (mut dl_v, mut du_v) = (dc.min(dt), dc.max(dt));
                                        if (du_v - dl_v).abs() < 0.01 {
                                            dl_v -= 0.06;
                                            du_v += 0.06;
                                        }
                                        (dl_v, du_v, Path14Type::Other)
                                    } else {
                                        // Check _checkAmideEster14 first (C=O in path)
                                        let pref14 = detect_amide_ester14(mol, ni, n1, n2, nj);
                                        if pref14 != AmidePref::None {
                                            match pref14 {
                                                AmidePref::Cis => {
                                                    (dc - 0.06, dc + 0.06, Path14Type::Cis)
                                                }
                                                AmidePref::Trans => {
                                                    (dt - 0.06, dt + 0.06, Path14Type::Trans)
                                                }
                                                AmidePref::None => unreachable!(),
                                            }
                                        } else {
                                            // Then check _checkAmideEster15 (C=O separate branch)
                                            let pref =
                                                detect_amide_ester_preference(mol, ni, n1, n2, nj);
                                            match pref {
                                                AmidePref::Cis => {
                                                    (dc - 0.06, dc + 0.06, Path14Type::Cis)
                                                }
                                                AmidePref::Trans => {
                                                    (dt - 0.06, dt + 0.06, Path14Type::Trans)
                                                }
                                                AmidePref::None => {
                                                    let (mut dl_v, mut du_v) =
                                                        (dc.min(dt), dc.max(dt));
                                                    if (du_v - dl_v).abs() < 0.01 {
                                                        dl_v -= 0.06;
                                                        du_v += 0.06;
                                                    }
                                                    (dl_v, du_v, Path14Type::Other)
                                                }
                                            }
                                        }
                                    }
                                }
                            } else {
                                // Non-sp2: still check amide/ester patterns (RDKit checks for ALL single bonds)
                                let pref14 = detect_amide_ester14(mol, ni, n1, n2, nj);
                                if pref14 != AmidePref::None {
                                    match pref14 {
                                        AmidePref::Cis => (dc - 0.06, dc + 0.06, Path14Type::Cis),
                                        AmidePref::Trans => {
                                            (dt - 0.06, dt + 0.06, Path14Type::Trans)
                                        }
                                        AmidePref::None => unreachable!(),
                                    }
                                } else {
                                    let pref = detect_amide_ester_preference(mol, ni, n1, n2, nj);
                                    match pref {
                                        AmidePref::Cis => (dc - 0.06, dc + 0.06, Path14Type::Cis),
                                        AmidePref::Trans => {
                                            (dt - 0.06, dt + 0.06, Path14Type::Trans)
                                        }
                                        AmidePref::None => {
                                            let (mut dl_v, mut du_v) = (dc.min(dt), dc.max(dt));
                                            if (du_v - dl_v).abs() < 0.01 {
                                                dl_v -= 0.06;
                                                du_v += 0.06;
                                            }
                                            (dl_v, du_v, Path14Type::Other)
                                        }
                                    }
                                }
                            };
                            // RDKit: only set bounds for paths NOT in small rings
                            if !in_small_ring {
                                set_b_union(&mut bounds, i, j, l.max(0.0), u);
                            }
                            // Store path for set15Bounds (both directions)
                            let (n1i, n2i) = (n1.index(), n2.index());
                            path14_list.push((i, n1i, n2i, j, path_type));
                            path14_list.push((j, n2i, n1i, i, path_type));
                            match path_type {
                                Path14Type::Cis => {
                                    cis_keys.insert((i, n1i, n2i, j));
                                    cis_keys.insert((j, n2i, n1i, i));
                                }
                                Path14Type::Trans => {
                                    trans_keys.insert((i, n1i, n2i, j));
                                    trans_keys.insert((j, n2i, n1i, i));
                                }
                                Path14Type::Other => {}
                            }
                        }
                    }
                }
            }
        }
    }
    // === set15Bounds ===
    if use_set15 {
        let mut set15_atoms: HashSet<(usize, usize)> = HashSet::new();
        let dist15_tol: f64 = 0.08;
        for &(a1, a2, a3, a4, path_type) in &path14_list {
            let na4 = NodeIndex::new(a4);
            for a5_node in mol.graph.neighbors(na4) {
                let a5 = a5_node.index();
                if a5 == a1 {
                    continue;
                }
                if top_dist[(a1, a5)] < 4 {
                    continue;
                }
                // Only set if no prior bounds, or already touched by set15
                let (lo, hi) = if a1 < a5 { (a1, a5) } else { (a5, a1) };
                let current_lb = bounds[(hi, lo)];
                if current_lb > 0.01
                    && !set15_atoms.contains(&(a1, a5))
                    && !set15_atoms.contains(&(a5, a1))
                {
                    continue;
                }
                let d1 = get_bond_length(mol, NodeIndex::new(a1), NodeIndex::new(a2));
                let d2 = get_bond_length(mol, NodeIndex::new(a2), NodeIndex::new(a3));
                let d3 = get_bond_length(mol, NodeIndex::new(a3), na4);
                let d4 = get_bond_length(mol, na4, a5_node);
                let ang12 = crate::graph::get_corrected_ideal_angle(
                    mol,
                    NodeIndex::new(a2),
                    NodeIndex::new(a1),
                    NodeIndex::new(a3),
                );
                let ang23 = crate::graph::get_corrected_ideal_angle(
                    mol,
                    NodeIndex::new(a3),
                    NodeIndex::new(a2),
                    na4,
                );
                let ang34 =
                    crate::graph::get_corrected_ideal_angle(mol, na4, NodeIndex::new(a3), a5_node);

                let second_is_cis = cis_keys.contains(&(a2, a3, a4, a5));
                let second_is_trans = trans_keys.contains(&(a2, a3, a4, a5));

                let (dl, mut du): (f64, f64) = match path_type {
                    Path14Type::Cis => {
                        if second_is_cis {
                            let d = compute15_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34);
                            (d - dist15_tol, d + dist15_tol)
                        } else if second_is_trans {
                            let d = compute15_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34);
                            (d - dist15_tol, d + dist15_tol)
                        } else {
                            (
                                compute15_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34) - dist15_tol,
                                compute15_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34)
                                    + dist15_tol,
                            )
                        }
                    }
                    Path14Type::Trans => {
                        if second_is_cis {
                            let d = compute15_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34);
                            (d - dist15_tol, d + dist15_tol)
                        } else if second_is_trans {
                            let d = compute15_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34);
                            (d - dist15_tol, d + dist15_tol)
                        } else {
                            (
                                compute15_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34)
                                    - dist15_tol,
                                compute15_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34)
                                    + dist15_tol,
                            )
                        }
                    }
                    Path14Type::Other => {
                        if second_is_cis {
                            (
                                compute15_cis_cis(d4, d3, d2, d1, ang34, ang23, ang12) - dist15_tol,
                                compute15_cis_trans(d4, d3, d2, d1, ang34, ang23, ang12)
                                    + dist15_tol,
                            )
                        } else if second_is_trans {
                            (
                                compute15_trans_cis(d4, d3, d2, d1, ang34, ang23, ang12)
                                    - dist15_tol,
                                compute15_trans_trans(d4, d3, d2, d1, ang34, ang23, ang12)
                                    + dist15_tol,
                            )
                        } else {
                            let vw1 =
                                crate::graph::get_vdw_radius(mol.graph[NodeIndex::new(a1)].element);
                            let vw5 = crate::graph::get_vdw_radius(mol.graph[a5_node].element);
                            (0.7 * (vw1 + vw5), -1.0)
                        }
                    }
                };
                if du < 0.0 {
                    du = 1000.0;
                }
                set_b_union(&mut bounds, a1, a5, dl, du);
                set15_atoms.insert((a1, a5));
                set15_atoms.insert((a5, a1));
            }
        }
    } // end use_set15
      // VDW lower bounds — only for pairs not already covered by 1-2/1-3/1-4/1-5
    for i in 0..n {
        for j in (i + 1)..n {
            let current_lower = bounds[(j, i)];
            if current_lower > 0.01 {
                continue;
            }
            if top_dist[(i, j)] >= 4 {
                let v = crate::graph::get_vdw_radius(mol.graph[NodeIndex::new(i)].element)
                    + crate::graph::get_vdw_radius(mol.graph[NodeIndex::new(j)].element);
                let l = match top_dist[(i, j)] {
                    4 => v * 0.7,
                    5 => v * 0.85,
                    _ => v,
                };
                set_b(&mut bounds, i, j, l, 1000.0);
            }
        }
    }
    bounds
}
