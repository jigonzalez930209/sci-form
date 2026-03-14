//! Validation checks matching RDKit's embedding validation pipeline.
//! Part of the retry-on-failure loop in embedPoints().

use crate::forcefield::bounds_ff::ChiralSet;
use crate::graph::Molecule;
use nalgebra::{DMatrix, Vector3};

const MIN_TETRAHEDRAL_CHIRAL_VOL: f64 = 0.50;
const TETRAHEDRAL_CENTERINVOLUME_TOL: f64 = 0.30;
pub const MAX_MINIMIZED_E_PER_ATOM: f32 = 0.05;

/// Tetrahedral center — any sp3 atom with 4 neighbors (not necessarily chiral)
pub struct TetrahedralCenter {
    pub center: usize,
    pub neighbors: [usize; 4],
    pub in_small_ring: bool,
}

/// Identify tetrahedral centers for volume checks.
/// Matches RDKit's findChiralSets logic for tetrahedralCarbons:
///   - C or N atoms only, with exactly 4 neighbors
///   - Must be in 2+ rings (ring junction atoms)
///   - Must NOT be in any 3-membered ring
pub fn identify_tetrahedral_centers(mol: &Molecule) -> Vec<TetrahedralCenter> {
    let n = mol.graph.node_count();
    // Compute SSSR rings, then derive per-atom ring count and 3-ring membership
    let rings = find_sssr(mol);
    let mut ring_count = vec![0usize; n];
    let mut in_3_ring = vec![false; n];
    let mut small_ring_count = vec![0usize; n]; // rings of size < 5
    for ring in &rings {
        for &atom_idx in ring {
            ring_count[atom_idx] += 1;
            if ring.len() == 3 {
                in_3_ring[atom_idx] = true;
            }
            if ring.len() < 5 {
                small_ring_count[atom_idx] += 1;
            }
        }
    }

    let mut centers = Vec::new();
    for i in 0..n {
        let ni = petgraph::graph::NodeIndex::new(i);
        let atom = &mol.graph[ni];
        // RDKit: only C (6) or N (7) with degree == 4
        let elem = atom.element;
        if elem != 6 && elem != 7 {
            continue;
        }
        let nbs: Vec<_> = mol.graph.neighbors(ni).collect();
        if nbs.len() != 4 {
            continue;
        }

        // RDKit: only add if in 2+ rings AND not in any 3-ring
        if ring_count[i] < 2 || in_3_ring[i] {
            continue;
        }

        centers.push(TetrahedralCenter {
            center: i,
            neighbors: [
                nbs[0].index(),
                nbs[1].index(),
                nbs[2].index(),
                nbs[3].index(),
            ],
            in_small_ring: small_ring_count[i] > 1,
        });
    }
    centers
}

/// Find the Smallest Set of Smallest Rings (SSSR) using Horton's algorithm.
/// Returns a list of rings, each ring being a vector of atom indices.
pub fn find_sssr_pub(mol: &Molecule) -> Vec<Vec<usize>> {
    find_sssr(mol)
}
fn find_sssr(mol: &Molecule) -> Vec<Vec<usize>> {
    use std::collections::VecDeque;
    let n = mol.graph.node_count();
    if n == 0 {
        return vec![];
    }

    // Number of independent cycles = edges - vertices + connected_components
    let num_edges = mol.graph.edge_count();
    // Count connected components via BFS
    let mut visited = vec![false; n];
    let mut num_components = 0;
    for start in 0..n {
        if visited[start] {
            continue;
        }
        num_components += 1;
        let mut queue = VecDeque::new();
        queue.push_back(start);
        visited[start] = true;
        while let Some(curr) = queue.pop_front() {
            for nb in mol.graph.neighbors(petgraph::graph::NodeIndex::new(curr)) {
                if !visited[nb.index()] {
                    visited[nb.index()] = true;
                    queue.push_back(nb.index());
                }
            }
        }
    }
    let cycle_rank = (num_edges + num_components).saturating_sub(n);
    if cycle_rank == 0 {
        return vec![];
    }

    // For each vertex, BFS to compute shortest path tree
    // Then for each non-tree edge found at vertex v, form the candidate ring
    let mut candidates: Vec<Vec<usize>> = Vec::new();

    for root in 0..n {
        let mut dist = vec![usize::MAX; n];
        let mut parent = vec![usize::MAX; n];
        dist[root] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(root);

        while let Some(curr) = queue.pop_front() {
            for nb in mol.graph.neighbors(petgraph::graph::NodeIndex::new(curr)) {
                let nb_idx = nb.index();
                if dist[nb_idx] == usize::MAX {
                    dist[nb_idx] = dist[curr] + 1;
                    parent[nb_idx] = curr;
                    queue.push_back(nb_idx);
                }
            }
        }

        // For each neighbor of root, if they share neighbors that are equidistant or close,
        // form candidate rings. Specifically, look for pairs (u, v) where edge (u,v) exists
        // and dist[u] + dist[v] + 1 gives an odd ring, or dist[u] == dist[v] for even ring.
        for u in 0..n {
            for nb in mol.graph.neighbors(petgraph::graph::NodeIndex::new(u)) {
                let v = nb.index();
                if u >= v {
                    continue;
                } // avoid duplicates
                let ring_len = dist[u] + dist[v] + 1;
                if ring_len > 8 {
                    continue;
                } // skip very large rings
                if dist[u] == usize::MAX || dist[v] == usize::MAX {
                    continue;
                }

                // Build the ring: path from root to u + edge (u,v) + path from v to root
                let path_u = trace_path(&parent, root, u);
                let path_v = trace_path(&parent, root, v);

                // Check that paths don't share intermediate vertices (would make it not a simple cycle)
                let mut ring = path_u.clone();
                // path_v goes root→...→v, we need to reverse it and skip the root
                let mut path_v_rev: Vec<usize> = path_v.into_iter().rev().collect();
                if !path_v_rev.is_empty() && !ring.is_empty() && path_v_rev.last() == ring.first() {
                    path_v_rev.pop(); // remove duplicate root
                }
                ring.extend(path_v_rev);

                // Check it's a valid simple cycle (no repeated vertices)
                let mut seen = std::collections::HashSet::new();
                let is_simple = ring.iter().all(|&x| seen.insert(x));
                if is_simple && ring.len() >= 3 {
                    // Normalize the ring for deduplication
                    let normalized = normalize_ring(&ring);
                    candidates.push(normalized);
                }
            }
        }
    }

    // Deduplicate candidates
    candidates.sort();
    candidates.dedup();

    // Filter: keep only rings that are "relevant" — not the XOR of two smaller rings.
    // For the purposes of ring counting, we keep all unique smallest rings.
    // Sort by size so smallest come first.
    candidates.sort_by_key(|r| r.len());

    // A ring is "relevant" if it cannot be expressed as the symmetric difference of
    // two strictly smaller rings. This matches RDKit's ring perception behavior.
    let edge_sets: Vec<std::collections::HashSet<(usize, usize)>> =
        candidates.iter().map(|r| ring_edges(r).collect()).collect();

    let mut relevant = Vec::new();
    for (i, ring) in candidates.iter().enumerate() {
        let mut is_xor_of_smaller = false;
        // Check all pairs of strictly smaller rings
        for j in 0..i {
            if candidates[j].len() >= ring.len() {
                continue;
            }
            for k in (j + 1)..i {
                if candidates[k].len() >= ring.len() {
                    continue;
                }
                // Symmetric difference of edge sets j and k
                let sym_diff: std::collections::HashSet<(usize, usize)> = edge_sets[j]
                    .symmetric_difference(&edge_sets[k])
                    .copied()
                    .collect();
                if sym_diff == edge_sets[i] {
                    is_xor_of_smaller = true;
                    break;
                }
            }
            if is_xor_of_smaller {
                break;
            }
        }
        if !is_xor_of_smaller {
            relevant.push(ring.clone());
        }
    }

    relevant
}

fn trace_path(parent: &[usize], root: usize, target: usize) -> Vec<usize> {
    let mut path = Vec::new();
    let mut curr = target;
    while curr != root && curr != usize::MAX {
        path.push(curr);
        curr = parent[curr];
    }
    if curr == root {
        path.push(root);
    }
    path.reverse();
    path
}

fn normalize_ring(ring: &[usize]) -> Vec<usize> {
    if ring.is_empty() {
        return vec![];
    }
    // Find minimum element position
    let min_pos = ring.iter().enumerate().min_by_key(|&(_, &v)| v).unwrap().0;
    let n = ring.len();
    // Try both directions (clockwise and counterclockwise)
    let forward: Vec<usize> = (0..n).map(|i| ring[(min_pos + i) % n]).collect();
    let backward: Vec<usize> = (0..n).map(|i| ring[(min_pos + n - i) % n]).collect();
    forward.min(backward)
}

fn ring_edges(ring: &[usize]) -> impl Iterator<Item = (usize, usize)> + '_ {
    let n = ring.len();
    (0..n).map(move |i| {
        let a = ring[i];
        let b = ring[(i + 1) % n];
        (a.min(b), a.max(b))
    })
}

/// Volume test: check that a tetrahedral center has minimum volume.
/// Uses NORMALIZED direction vectors from center to each neighbor.
/// Checks all C(4,3)=4 combinations of 3 vectors.
/// Uses f64 to match RDKit's Point3D (double) precision.
fn volume_test(
    center: usize,
    neighbors: &[usize; 4],
    coords: &DMatrix<f64>,
    relaxed: bool,
) -> bool {
    let dim = coords.ncols().min(3);
    let p0 = Vector3::new(
        coords[(center, 0)],
        coords[(center, 1)],
        if dim >= 3 { coords[(center, 2)] } else { 0.0 },
    );
    let mut vecs = [Vector3::<f64>::zeros(); 4];
    for (k, &nb) in neighbors.iter().enumerate() {
        let pk = Vector3::new(
            coords[(nb, 0)],
            coords[(nb, 1)],
            if dim >= 3 { coords[(nb, 2)] } else { 0.0 },
        );
        let v = p0 - pk; // RDKit: center - neighbor
        let norm = v.norm();
        vecs[k] = if norm > 1e-8 { v / norm } else { v };
    }

    let vol_scale: f64 = if relaxed { 0.25 } else { 1.0 };
    let threshold = vol_scale * MIN_TETRAHEDRAL_CHIRAL_VOL;

    // RDKit checks: (v1×v2)·v3, (v1×v2)·v4, (v1×v3)·v4, (v2×v3)·v4
    let combos: [(usize, usize, usize); 4] = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)];
    for (a, b, c) in combos {
        let cross = vecs[a].cross(&vecs[b]);
        let vol = cross.dot(&vecs[c]).abs();
        if vol < threshold {
            return false;
        }
    }
    true
}

/// Center-in-volume test: center atom must be inside the tetrahedron formed by its 4 neighbors.
fn same_side(
    v1: &Vector3<f64>,
    v2: &Vector3<f64>,
    v3: &Vector3<f64>,
    v4: &Vector3<f64>,
    p0: &Vector3<f64>,
    tol: f64,
) -> bool {
    let normal = (v2 - v1).cross(&(v3 - v1));
    let d1 = normal.dot(&(v4 - v1));
    let d2 = normal.dot(&(p0 - v1));
    if d1.abs() < tol || d2.abs() < tol {
        return false;
    }
    (d1 < 0.0) == (d2 < 0.0)
}

fn center_in_volume(
    center: usize,
    neighbors: &[usize; 4],
    coords: &DMatrix<f64>,
    tol: f64,
) -> bool {
    let dim = coords.ncols().min(3);
    let get_p3d = |idx: usize| -> Vector3<f64> {
        Vector3::new(
            coords[(idx, 0)],
            coords[(idx, 1)],
            if dim >= 3 { coords[(idx, 2)] } else { 0.0 },
        )
    };
    let p0 = get_p3d(center);
    let p = [
        get_p3d(neighbors[0]),
        get_p3d(neighbors[1]),
        get_p3d(neighbors[2]),
        get_p3d(neighbors[3]),
    ];

    same_side(&p[0], &p[1], &p[2], &p[3], &p0, tol)
        && same_side(&p[1], &p[2], &p[3], &p[0], &p0, tol)
        && same_side(&p[2], &p[3], &p[0], &p[1], &p0, tol)
        && same_side(&p[3], &p[0], &p[1], &p[2], &p0, tol)
}

/// Check all tetrahedral centers for minimum volume and center-in-volume.
/// Matches RDKit's checkTetrahedralCenters. Uses f64 coords matching RDKit's Point3D.
pub fn check_tetrahedral_centers(coords: &DMatrix<f64>, centers: &[TetrahedralCenter]) -> bool {
    for tc in centers {
        if !volume_test(tc.center, &tc.neighbors, coords, tc.in_small_ring) {
            return false;
        }
        if !center_in_volume(
            tc.center,
            &tc.neighbors,
            coords,
            TETRAHEDRAL_CENTERINVOLUME_TOL,
        ) {
            return false;
        }
    }
    true
}

/// Check chiral center volumes have correct sign.
/// Matches RDKit's checkChiralCenters — intentionally permissive (allows 20% undershoot if sign matches).
/// Uses f64 coords matching RDKit's Point3D.
pub fn check_chiral_centers(coords: &DMatrix<f64>, chiral_sets: &[ChiralSet]) -> bool {
    for cs in chiral_sets {
        let vol = crate::distgeom::calc_chiral_volume_f64(
            cs.neighbors[0],
            cs.neighbors[1],
            cs.neighbors[2],
            cs.neighbors[3],
            coords,
        );
        let lb = cs.lower_vol as f64;
        let ub = cs.upper_vol as f64;
        if lb > 0.0 && vol < lb && (vol / lb < 0.8 || have_opposite_sign(vol, lb)) {
            return false;
        }
        if ub < 0.0 && vol > ub && (vol / ub < 0.8 || have_opposite_sign(vol, ub)) {
            return false;
        }
    }
    true
}

fn have_opposite_sign(a: f64, b: f64) -> bool {
    (a < 0.0) != (b < 0.0)
}

/// Planarity check: compute OOP (improper torsion) energy for SP2 centers.
/// Reject if energy > n_impropers * tolerance.
/// Matches RDKit's planarity check in minimizeWithExpTorsions.
pub fn check_planarity(mol: &Molecule, coords: &DMatrix<f32>, oop_k: f32, tolerance: f32) -> bool {
    let n = mol.graph.node_count();
    let mut n_impropers = 0usize;
    let mut improper_energy = 0.0f32;

    // SP2 improper (out-of-plane) terms only
    for i in 0..n {
        let ni = petgraph::graph::NodeIndex::new(i);
        if mol.graph[ni].hybridization != crate::graph::Hybridization::SP2 {
            continue;
        }
        let nbs: Vec<_> = mol.graph.neighbors(ni).collect();
        if nbs.len() != 3 {
            continue;
        }
        n_impropers += 1;

        let pc = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
        let p1 = Vector3::new(
            coords[(nbs[0].index(), 0)],
            coords[(nbs[0].index(), 1)],
            coords[(nbs[0].index(), 2)],
        );
        let p2 = Vector3::new(
            coords[(nbs[1].index(), 0)],
            coords[(nbs[1].index(), 1)],
            coords[(nbs[1].index(), 2)],
        );
        let p3 = Vector3::new(
            coords[(nbs[2].index(), 0)],
            coords[(nbs[2].index(), 1)],
            coords[(nbs[2].index(), 2)],
        );
        let v1 = p1 - pc;
        let v2 = p2 - pc;
        let v3 = p3 - pc;
        let vol = v1.dot(&v2.cross(&v3));
        improper_energy += oop_k * vol * vol;
    }

    // SP linearity is enforced by the ETKDG 3D FF distance constraints (k=100),
    // not by the validation check. Including SP angle penalties here caused
    // false rejections for molecules with both SP and SP2 atoms.

    if n_impropers == 0 {
        return true;
    }
    improper_energy <= n_impropers as f32 * tolerance
}

/// Double bond geometry check: reject if substituent-double_bond_atom-other is nearly linear.
/// Matches RDKit's doubleBondGeometryChecks with doubleBondEnds filtering.
/// Uses f64 coords matching RDKit's Point3D.
pub fn check_double_bond_geometry(mol: &Molecule, coords: &DMatrix<f64>) -> bool {
    use petgraph::visit::EdgeRef;
    for edge in mol.graph.edge_references() {
        if mol.graph[edge.id()].order != crate::graph::BondOrder::Double {
            continue;
        }
        let u = edge.source();
        let v = edge.target();

        // Check neighbors of u (substituents around the double bond end)
        let u_deg = mol.graph.neighbors(u).count();
        if u_deg >= 2 {
            for nb in mol.graph.neighbors(u) {
                if nb == v {
                    continue;
                }
                // RDKit filter: skip if bond to neighbor is NOT single and atom has degree 2
                if u_deg == 2 {
                    if let Some(eid) = mol.graph.find_edge(u, nb) {
                        if mol.graph[eid].order != crate::graph::BondOrder::Single {
                            continue;
                        }
                    }
                }
                if !check_linearity(nb.index(), u.index(), v.index(), coords) {
                    return false;
                }
            }
        }
        // Check neighbors of v
        let v_deg = mol.graph.neighbors(v).count();
        if v_deg >= 2 {
            for nb in mol.graph.neighbors(v) {
                if nb == u {
                    continue;
                }
                // RDKit filter: skip if bond to neighbor is NOT single and atom has degree 2
                if v_deg == 2 {
                    if let Some(eid) = mol.graph.find_edge(v, nb) {
                        if mol.graph[eid].order != crate::graph::BondOrder::Single {
                            continue;
                        }
                    }
                }
                if !check_linearity(nb.index(), v.index(), u.index(), coords) {
                    return false;
                }
            }
        }
    }
    true
}

/// Returns false if a0-a1-a2 is nearly linear (angle ≈ 180°).
fn check_linearity(a0: usize, a1: usize, a2: usize, coords: &DMatrix<f64>) -> bool {
    let p0 = Vector3::new(coords[(a0, 0)], coords[(a0, 1)], coords[(a0, 2)]);
    let p1 = Vector3::new(coords[(a1, 0)], coords[(a1, 1)], coords[(a1, 2)]);
    let p2 = Vector3::new(coords[(a2, 0)], coords[(a2, 1)], coords[(a2, 2)]);
    let mut v1 = p1 - p0;
    let n1 = v1.norm();
    if n1 < 1e-8 {
        return true;
    }
    v1 /= n1;
    let mut v2 = p1 - p2;
    let n2 = v2.norm();
    if n2 < 1e-8 {
        return true;
    }
    v2 /= n2;
    // dot ≈ -1 means linear; reject if dot + 1 < 1e-3
    v1.dot(&v2) + 1.0 >= 1e-3
}
