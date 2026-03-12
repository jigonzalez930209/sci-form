//! ETKDG 3D force field builder — constructs force field terms matching RDKit.

use nalgebra::{DMatrix, Vector3};
use petgraph::visit::EdgeRef;
use super::*;

/// Build the 3D ETKDG force field matching RDKit's construct3DForceField exactly.
///
/// RDKit's order: torsion 1-4 pairs marked first → improper terms → 1-2 → 1-3 → long-range.
/// Torsion 1-4 pairs and all 1-2/1-3 pairs are EXCLUDED from long-range constraints.
/// For 1-3 pairs where the center is an improper-constrained atom, bounds distance is used.
pub fn build_etkdg_3d_ff(
    mol: &crate::graph::Molecule,
    coords: &DMatrix<f64>,
    bounds_matrix: &DMatrix<f64>,
) -> Etkdg3DFF {
    build_etkdg_3d_ff_with_torsions(mol, coords, bounds_matrix, &[])
}

/// Build the 3D ETKDG force field with explicit torsion contributions.
///
/// Matches RDKit's construct3DForceField ordering:
///   1. Mark torsion 1-4 endpoint pairs (from CSD + pattern torsions)
///   2. Mark improper-constrained atoms (SP2 C/N/O with 3 neighbors)
///   3. Add 1-2 bond distance constraints
///   4. Add 1-3 angle distance constraints (using bounds for improper centers)
///   5. Add long-range constraints (excluding all marked pairs)
pub fn build_etkdg_3d_ff_with_torsions(
    mol: &crate::graph::Molecule,
    coords: &DMatrix<f64>,
    bounds_matrix: &DMatrix<f64>,
    torsion_contribs: &[M6TorsionContrib],
) -> Etkdg3DFF {
    let n = mol.graph.node_count();
    let mut dist_12 = Vec::new();
    let mut dist_13 = Vec::new();
    let mut atom_pairs = vec![false; n * n];

    // === Build pattern-based torsion contributions first ===
    let mut torsion_contribs_internal = Vec::new();

    // Basic knowledge: flat ring torsions
    // For 4-6 membered rings where all 4 atoms in a torsion i-j-k-l are SP2,
    // add V2=100.0 with s2=-1 (strong planar preference), matching RDKit.
    // RDKit adds exactly ONE torsion per ring bond (first qualifying atom pair).
    {
        let mut done_bonds = std::collections::HashSet::new();
        for edge in mol.graph.edge_references() {
            let u = edge.source();
            let v = edge.target();
            let key = if u.index() < v.index() { (u.index(), v.index()) } else { (v.index(), u.index()) };
            if done_bonds.contains(&key) { continue; }
            // Check if u-v bond is in a ring of size 4-6.
            // Find the shortest alternative path from u to v NOT using the direct u-v bond.
            // BFS from u's neighbors (excluding v) to reach v, excluding u.
            let alt_path = mol.graph.neighbors(u)
                .filter(|&x| x != v)
                .filter_map(|nu| crate::graph::min_path_excluding(mol, nu, v, u, 5).map(|d| d + 1))
                .min();
            if let Some(alt_len) = alt_path {
                let ring_size = alt_len + 1; // alternative path + direct bond
                if ring_size >= 4 && ring_size <= 6 {
                    // Both u and v must be SP2
                    if mol.graph[u].hybridization != crate::graph::Hybridization::SP2 { continue; }
                    if mol.graph[v].hybridization != crate::graph::Hybridization::SP2 { continue; }
                    // Find the FIRST qualifying SP2 neighbor pair in the ring.
                    // RDKit adds exactly one torsion per ring bond.
                    let nbs_u: Vec<_> = mol.graph.neighbors(u).filter(|&x| x != v).collect();
                    let nbs_v: Vec<_> = mol.graph.neighbors(v).filter(|&x| x != u).collect();
                    let mut added = false;
                    for &nu in &nbs_u {
                        if added { break; }
                        if mol.graph[nu].hybridization != crate::graph::Hybridization::SP2 { continue; }
                        // Check if nu is in the same ring through u-v
                        let ring_path_nu = crate::graph::min_path_excluding(mol, nu, v, u, 5);
                        if ring_path_nu.is_none() { continue; }
                        for &nv in &nbs_v {
                            if nv == nu { continue; }
                            if mol.graph[nv].hybridization != crate::graph::Hybridization::SP2 { continue; }
                            // Check if nv is in the ring through u-v
                            let ring_path_nv = crate::graph::min_path_excluding(mol, u, nv, v, 5);
                            if ring_path_nv.is_none() { continue; }
                            // All 4 atoms are SP2 and in a 4-6 ring → add flat torsion
                            torsion_contribs_internal.push(M6TorsionContrib {
                                i: nu.index(),
                                j: u.index(),
                                k: v.index(),
                                l: nv.index(),
                                signs: [0.0, -1.0, 0.0, 0.0, 0.0, 0.0],
                                v: [0.0, 100.0, 0.0, 0.0, 0.0, 0.0],
                            });
                            added = true;
                            break;
                        }
                    }
                    done_bonds.insert(key);
                }
            }
        }
    }

    // Chain torsion preferences based on RDKit's torsionPreferences_v2.in
    // Simplified atom-type matching for the most impactful generic patterns
    {
        use petgraph::graph::NodeIndex;
        let mut done_chain_bonds = std::collections::HashSet::new();
        for edge in mol.graph.edge_references() {
            let u = edge.source();
            let v = edge.target();
            let bond = &mol.graph[edge.id()];
            
            let key = if u.index() < v.index() { (u.index(), v.index()) } else { (v.index(), u.index()) };
            if done_chain_bonds.contains(&key) { continue; }

            let eu = mol.graph[u].element;
            let ev = mol.graph[v].element;
            let hybu = &mol.graph[u].hybridization;
            let hybv = &mol.graph[v].hybridization;

            // Helper: check if bond is in a ring (path exists from u to v not using this edge)
            let in_ring = crate::graph::min_path_excluding2(mol, u, v, u, v, 7).is_some();

            // Get heavy (non-H) neighbors
            let nbs_u: Vec<_> = mol.graph.neighbors(u).filter(|&x| x != v).collect();
            let nbs_v: Vec<_> = mol.graph.neighbors(v).filter(|&x| x != u).collect();
            if nbs_u.is_empty() || nbs_v.is_empty() { continue; }

            // Helper: check if a node has a double bond to oxygen (carbonyl)
            let has_double_bond_to_o = |node: NodeIndex, exclude: NodeIndex| -> Option<NodeIndex> {
                mol.graph.neighbors(node).find(|&nb| {
                    if nb == exclude { return false; }
                    if mol.graph[nb].element != 8 { return false; }
                    if let Some(e) = mol.graph.find_edge(node, nb) {
                        mol.graph[e].order == crate::graph::BondOrder::Double
                    } else { false }
                })
            };

            // Helper: check if node is aromatic
            let is_aromatic = |node: NodeIndex| -> bool {
                mol.graph.neighbors(node).any(|nb| {
                    if let Some(e) = mol.graph.find_edge(node, nb) {
                        mol.graph[e].order == crate::graph::BondOrder::Aromatic
                    } else { false }
                })
            };

            // ===== DOUBLE BOND PATTERNS =====
            // [*:1][X3,X2:2]=[X3,X2:3][*:4] → V2=100, s2=-1 (planar around double bond)
            if bond.order == crate::graph::BondOrder::Double ||
               bond.order == crate::graph::BondOrder::Aromatic {
                // Already handled by flat ring torsions for ring bonds and OOP terms
                // For non-ring double bonds, enforce planarity
                if !in_ring {
                    for &nu in &nbs_u {
                        for &nv in &nbs_v {
                            if nu == nv { continue; }
                            torsion_contribs_internal.push(M6TorsionContrib {
                                i: nu.index(),
                                j: u.index(),
                                k: v.index(),
                                l: nv.index(),
                                signs: [0.0, -1.0, 0.0, 0.0, 0.0, 0.0],
                                v: [0.0, 100.0, 0.0, 0.0, 0.0, 0.0],
                            });
                        }
                    }
                    done_chain_bonds.insert(key);
                }
                continue;
            }

            // Only single bonds below
            if bond.order != crate::graph::BondOrder::Single { continue; }
            // Skip ring bonds for now (handled separately by flat ring torsions)
            if in_ring { continue; }

            // ===== ESTER/AMIDE: O=C-X patterns =====
            // Check if u or v is a carbonyl carbon bonded to O/N
            let (carbonyl_node, hetero_node) = 
                if eu == 6 && *hybu == crate::graph::Hybridization::SP2 && (ev == 8 || ev == 7) {
                    (u, v)
                } else if ev == 6 && *hybv == crate::graph::Hybridization::SP2 && (eu == 8 || eu == 7) {
                    (v, u)
                } else {
                    (NodeIndex::new(n), NodeIndex::new(n)) // sentinel: no match
                };

            if carbonyl_node.index() < n {
                if let Some(carbonyl_o) = has_double_bond_to_o(carbonyl_node, hetero_node) {
                    // O=C-X-R → trans preference
                    // RDKit: V1=100.0, s1=-1 for ester (O=C-O-C)
                    // RDKit: V1=100.0, s1=-1 for amide (O=C-NH-*) in many patterns
                    let hetero_nbs: Vec<_> = mol.graph.neighbors(hetero_node)
                        .filter(|&x| x != carbonyl_node)
                        .collect();
                    for &hn in &hetero_nbs {
                        torsion_contribs_internal.push(M6TorsionContrib {
                            i: carbonyl_o.index(),
                            j: carbonyl_node.index(),
                            k: hetero_node.index(),
                            l: hn.index(),
                            signs: [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            v: [100.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        });
                    }
                    done_chain_bonds.insert(key);
                    continue;
                }
            }

            // ===== AROMATIC-X SINGLE BONDS =====
            // Ar-O-C, Ar-N-*, Ar-C(=O) patterns
            // Check if one end is aromatic
            let (ar_node, other_node) = if is_aromatic(u) && !is_aromatic(v) {
                (u, v)
            } else if is_aromatic(v) && !is_aromatic(u) {
                (v, u)
            } else if is_aromatic(u) && is_aromatic(v) {
                // Ar-Ar: biaryl
                // [a:1][a:2]-[a:3][a:4] patterns have V2 terms
                // Generic: V2=3-8:
                for &nu in &nbs_u {
                    for &nv in &nbs_v {
                        if nu == nv { continue; }
                        torsion_contribs_internal.push(M6TorsionContrib {
                            i: nu.index(),
                            j: u.index(),
                            k: v.index(),
                            l: nv.index(),
                            signs: [0.0, -1.0, 0.0, 0.0, 0.0, 0.0],
                            v: [0.0, 5.0, 0.0, 0.0, 0.0, 0.0],
                        });
                    }
                }
                done_chain_bonds.insert(key);
                continue;
            } else {
                (NodeIndex::new(n), NodeIndex::new(n))
            };

            if ar_node.index() < n {
                let e_other = mol.graph[other_node].element;
                let hyb_other = &mol.graph[other_node].hybridization;
                
                // Ar-C(=O) pattern: V2=5-10, s2=-1 (planar)
                if e_other == 6 && *hyb_other == crate::graph::Hybridization::SP2 {
                    if has_double_bond_to_o(other_node, ar_node).is_some() {
                        // Ar-C(=O): moderately strong planarity
                        let ar_nbs: Vec<_> = mol.graph.neighbors(ar_node).filter(|&x| x != other_node).collect();
                        let other_nbs: Vec<_> = mol.graph.neighbors(other_node).filter(|&x| x != ar_node).collect();
                        for &na in &ar_nbs {
                            for &no in &other_nbs {
                                if na == no { continue; }
                                torsion_contribs_internal.push(M6TorsionContrib {
                                    i: na.index(),
                                    j: ar_node.index(),
                                    k: other_node.index(),
                                    l: no.index(),
                                    signs: [0.0, -1.0, 0.0, 0.0, 0.0, 0.0],
                                    v: [0.0, 8.0, 0.0, 0.0, 0.0, 0.0],
                                });
                            }
                        }
                        done_chain_bonds.insert(key);
                        continue;
                    }
                }

                // Ar-CX4 (aromatic-sp3): V2=3.6, s2=-1 or V3 depending on substitution
                // Generic fallback: [a:1][c:2]-[CX4H2:3][*:4] → V2=3.6
                if e_other == 6 && *hyb_other == crate::graph::Hybridization::SP3 {
                    let ar_nbs: Vec<_> = mol.graph.neighbors(ar_node).filter(|&x| x != other_node).collect();
                    let other_nbs: Vec<_> = mol.graph.neighbors(other_node).filter(|&x| x != ar_node).collect();
                    for &na in &ar_nbs {
                        for &no in &other_nbs {
                            if na == no { continue; }
                            torsion_contribs_internal.push(M6TorsionContrib {
                                i: na.index(),
                                j: ar_node.index(),
                                k: other_node.index(),
                                l: no.index(),
                                signs: [0.0, -1.0, 0.0, 0.0, 0.0, 0.0],
                                v: [0.0, 3.6, 0.0, 0.0, 0.0, 0.0],
                            });
                        }
                    }
                    done_chain_bonds.insert(key);
                    continue;
                }

                // Ar-O, Ar-N: V2 planarity preference 
                // [a:1][c:2]-[OX2:3][*:4] → V2=7.2, s2=-1
                // [a:1][c:2]-[NX3:3][*:4] → V2=3.0-8.0
                if e_other == 8 || e_other == 7 {
                    let ar_nbs: Vec<_> = mol.graph.neighbors(ar_node).filter(|&x| x != other_node).collect();
                    let other_nbs: Vec<_> = mol.graph.neighbors(other_node).filter(|&x| x != ar_node).collect();
                    let v2_val = if e_other == 8 { 7.2 } else { 5.0 };
                    for &na in &ar_nbs {
                        for &no in &other_nbs {
                            if na == no { continue; }
                            torsion_contribs_internal.push(M6TorsionContrib {
                                i: na.index(),
                                j: ar_node.index(),
                                k: other_node.index(),
                                l: no.index(),
                                signs: [0.0, -1.0, 0.0, 0.0, 0.0, 0.0],
                                v: [0.0, v2_val, 0.0, 0.0, 0.0, 0.0],
                            });
                        }
                    }
                    done_chain_bonds.insert(key);
                    continue;
                }

                done_chain_bonds.insert(key);
                continue;
            }

            // ===== SP3-SP3 and other single bonds =====
            // [!#1:1][CX4:2]-[CX4:3][!#1:4] → V3=7.0, s3=1 (staggered preference)
            if eu == 6 && ev == 6 &&
                *hybu == crate::graph::Hybridization::SP3 && *hybv == crate::graph::Hybridization::SP3 {
                // Get heavy non-H neighbors only
                let heavy_u: Vec<_> = nbs_u.iter().filter(|&&x| mol.graph[x].element != 1).copied().collect();
                let heavy_v: Vec<_> = nbs_v.iter().filter(|&&x| mol.graph[x].element != 1).copied().collect();
                if !heavy_u.is_empty() && !heavy_v.is_empty() {
                    // Use first heavy neighbor pair
                    torsion_contribs_internal.push(M6TorsionContrib {
                        i: heavy_u[0].index(),
                        j: u.index(),
                        k: v.index(),
                        l: heavy_v[0].index(),
                        signs: [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                        v: [0.0, 0.0, 7.0, 0.0, 0.0, 0.0],
                    });
                } else {
                    // At least H neighbors: weaker preference
                    torsion_contribs_internal.push(M6TorsionContrib {
                        i: nbs_u[0].index(),
                        j: u.index(),
                        k: v.index(),
                        l: nbs_v[0].index(),
                        signs: [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                        v: [0.0, 0.0, 4.0, 0.0, 0.0, 0.0],
                    });
                }
                done_chain_bonds.insert(key);
                continue;
            }

            // C(sp3)-O or C(sp3)-N: ether/amine torsion
            // [*:1][CX4:2]-[O:3][CX4:4] → V3=2.5, s3=1
            if (eu == 6 && *hybu == crate::graph::Hybridization::SP3 && (ev == 8 || ev == 7)) ||
               (ev == 6 && *hybv == crate::graph::Hybridization::SP3 && (eu == 8 || eu == 7)) {
                let (c_node, hetero) = if eu == 6 && *hybu == crate::graph::Hybridization::SP3 { (u, v) } else { (v, u) };
                let c_nbs: Vec<_> = mol.graph.neighbors(c_node).filter(|&x| x != hetero).collect();
                let h_nbs: Vec<_> = mol.graph.neighbors(hetero).filter(|&x| x != c_node).collect();
                if !c_nbs.is_empty() && !h_nbs.is_empty() {
                    torsion_contribs_internal.push(M6TorsionContrib {
                        i: c_nbs[0].index(),
                        j: c_node.index(),
                        k: hetero.index(),
                        l: h_nbs[0].index(),
                        signs: [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                        v: [0.0, 0.0, 2.5, 0.0, 0.0, 0.0],
                    });
                    done_chain_bonds.insert(key);
                }
            }
        }
    }

    // === Merge torsion contributions matching RDKit's logic ===
    // RDKit's getExperimentalTorsions:
    //   1. CSD patterns run first, marking covered bonds in doneBonds
    //   2. Basic knowledge runs second:
    //      a. Improper atoms always added (independent of doneBonds)
    //      b. Ring torsions only for bonds NOT covered by CSD
    //
    // Our approach: separate ring torsions from chain torsions in torsion_contribs_internal.
    // When CSD data exists: use CSD + ring torsions (for uncovered bonds only).
    // When no CSD data: use all internal torsions (ring + chain).
    let all_torsions = if !torsion_contribs.is_empty() {
        // Determine which bonds are covered by CSD torsions
        let mut csd_done_bonds = std::collections::HashSet::new();
        for tc in torsion_contribs {
            let (j, k) = if tc.j < tc.k { (tc.j, tc.k) } else { (tc.k, tc.j) };
            csd_done_bonds.insert((j, k));
        }
        // Start with CSD torsions
        let mut merged = torsion_contribs.to_vec();
        // Add ring torsions for bonds NOT covered by CSD
        for tc in &torsion_contribs_internal {
            let (j, k) = if tc.j < tc.k { (tc.j, tc.k) } else { (tc.k, tc.j) };
            // Ring torsions have V2=100.0 — only add ring flatness for uncovered ring bonds
            if tc.v[1] >= 99.0 && tc.signs[1] < 0.0 {
                if !csd_done_bonds.contains(&(j, k)) {
                    merged.push(tc.clone());
                }
            }
        }
        merged
    } else {
        torsion_contribs_internal
    };

    // === Step 1: Mark torsion 1-4 endpoint pairs ===
    for tc in &all_torsions {
        let (i, l) = if tc.i < tc.l { (tc.i, tc.l) } else { (tc.l, tc.i) };
        atom_pairs[i * n + l] = true;
        atom_pairs[l * n + i] = true;
    }

    // === Step 2: Identify improper-constrained atoms and build UFF inversion contribs ===
    // RDKit: SP2 C/N/O with exactly 3 neighbors, 3 permutations each
    let mut is_improper_constrained = vec![false; n];
    let mut inversion_contribs = Vec::new();
    let oob_force_scaling = 10.0f64;
    for i in 0..n {
        let ni = petgraph::graph::NodeIndex::new(i);
        let elem = mol.graph[ni].element;
        if elem != 6 && elem != 7 && elem != 8 { continue; }
        if mol.graph[ni].hybridization != crate::graph::Hybridization::SP2 { continue; }
        let nbs: Vec<_> = mol.graph.neighbors(ni).collect();
        if nbs.len() != 3 { continue; }
        is_improper_constrained[i] = true;

        // Check if central C is bound to SP2 O
        let is_c_bound_to_sp2o = elem == 6 && nbs.iter().any(|&nb| {
            mol.graph[nb].element == 8
                && mol.graph[nb].hybridization == crate::graph::Hybridization::SP2
        });

        // UFF inversion coefficients matching RDKit's calcInversionCoefficientsAndForceConstant
        let (k_base, c0, c1, c2) = (
            if is_c_bound_to_sp2o { 50.0f64 } else { 6.0 },
            1.0f64, -1.0f64, 0.0f64,
        );
        let k = oob_force_scaling * k_base / 3.0;

        // RDKit's 3 permutations: [a0,center,a1,a2], matching addImproperTorsionTerms
        // neighbors stored as nbs[0]=a0, nbs[1]=a1, nbs[2]=a2 (indices 0,2,3 in RDKit's layout)
        let a0 = nbs[0].index();
        let a1 = nbs[1].index();
        let a2 = nbs[2].index();
        // Permutation 0: (a0, center, a1, a2)
        inversion_contribs.push(UFFInversionContrib {
            at1: a0, at2: i, at3: a1, at4: a2,
            force_constant: k, c0, c1, c2,
        });
        // Permutation 1: (a0, center, a2, a1)
        inversion_contribs.push(UFFInversionContrib {
            at1: a0, at2: i, at3: a2, at4: a1,
            force_constant: k, c0, c1, c2,
        });
        // Permutation 2: (a1, center, a2, a0)
        inversion_contribs.push(UFFInversionContrib {
            at1: a1, at2: i, at3: a2, at4: a0,
            force_constant: k, c0, c1, c2,
        });
    }

    // === Step 3: 1-2 constraints (bonds) ===
    for edge in mol.graph.edge_references() {
        let i = edge.source().index();
        let j = edge.target().index();
        atom_pairs[i * n + j] = true;
        atom_pairs[j * n + i] = true;
        let p1 = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
        let p2 = Vector3::new(coords[(j, 0)], coords[(j, 1)], coords[(j, 2)]);
        let d = (p1 - p2).norm();
        dist_12.push(DistConstraint {
            i, j,
            min_len: d - KNOWN_DIST_TOL,
            max_len: d + KNOWN_DIST_TOL,
            k: KNOWN_DIST_K,
        });
    }

    // === Step 4: 1-3 constraints (angle-end-atom distances) ===
    // For SP centers with triple bond or allene geometry: add AngleConstraint (179-180°, k=1)
    // For improper-constrained centers: use bounds matrix distances
    // For others: use current distance ± TOL
    // NOTE: Unlike long-range (step 5), 1-3 constraints are ALWAYS added,
    // even if the pair was already marked by torsion endpoints in step 1.
    // The atomPairs check is only for step 5 (long-range).
    let mut angle_constraints = Vec::new();
    for center in 0..n {
        let nc = petgraph::graph::NodeIndex::new(center);
        let nbs: Vec<_> = mol.graph.neighbors(nc).collect();
        for a in 0..nbs.len() {
            for b in (a + 1)..nbs.len() {
                let i = nbs[a].index();
                let j = nbs[b].index();
                // Always mark and add (matching RDKit's add13Terms)
                atom_pairs[i * n + j] = true;
                atom_pairs[j * n + i] = true;

                // Check if this is an SP (linear) angle:
                // Either bond is triple, or both bonds are double and center has degree 2
                let is_sp_angle = {
                    let bond_a = mol.graph.find_edge(nc, nbs[a])
                        .map(|e| mol.graph[e].order);
                    let bond_b = mol.graph.find_edge(nc, nbs[b])
                        .map(|e| mol.graph[e].order);
                    match (bond_a, bond_b) {
                        (Some(a_ord), Some(b_ord)) => {
                            a_ord == crate::graph::BondOrder::Triple
                                || b_ord == crate::graph::BondOrder::Triple
                                || (a_ord == crate::graph::BondOrder::Double
                                    && b_ord == crate::graph::BondOrder::Double
                                    && nbs.len() == 2)
                        }
                        _ => false,
                    }
                };

                if is_sp_angle {
                    // RDKit: AngleConstraint(179-180°, k=1) for SP/linear centers
                    angle_constraints.push(AngleConstraint {
                        i, j: center, k: j,
                        min_deg: 179.0,
                        max_deg: 180.0,
                        force_k: 1.0,
                    });
                } else if is_improper_constrained[center] {
                    // Use bounds matrix distances (matching RDKit)
                    let (lo, hi) = if i < j {
                        (bounds_matrix[(j, i)], bounds_matrix[(i, j)])
                    } else {
                        (bounds_matrix[(i, j)], bounds_matrix[(j, i)])
                    };
                    dist_13.push(DistConstraint {
                        i, j, min_len: lo, max_len: hi, k: KNOWN_DIST_K,
                    });
                } else {
                    let p1 = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
                    let p2 = Vector3::new(coords[(j, 0)], coords[(j, 1)], coords[(j, 2)]);
                    let d = (p1 - p2).norm();
                    dist_13.push(DistConstraint {
                        i, j,
                        min_len: d - KNOWN_DIST_TOL,
                        max_len: d + KNOWN_DIST_TOL,
                        k: KNOWN_DIST_K,
                    });
                }
            }
        }
    }

    // === Step 5: Long-range constraints (all remaining pairs) ===
    // Matching RDKit's iteration order: i in 1..N, j in 0..i
    // and addContrib(i, j, ...) where i > j
    let long_range_k = 10.0f64;
    let mut dist_long = Vec::new();
    for i in 1..n {
        for j in 0..i {
            if atom_pairs[j * n + i] { continue; }
            let lower = bounds_matrix[(i, j)];
            let upper = bounds_matrix[(j, i)];
            if upper < 1e-6 { continue; }
            dist_long.push(DistConstraint {
                i, j,
                min_len: lower,
                max_len: upper,
                k: long_range_k,
            });
        }
    }

    Etkdg3DFF {
        dist_12,
        dist_13,
        dist_long,
        angle_constraints,
        torsion_contribs: all_torsions,
        inversion_contribs,
        oop_k: 10.0f64,
        torsion_k_omega: 0.0f64,
        bounds_force_scaling: 1.0f64,
        use_m6_torsions: false,
    }
}
