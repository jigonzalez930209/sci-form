use nalgebra::DMatrix;
use petgraph::visit::EdgeRef;
use std::collections::{HashSet, VecDeque};

/// Identifies rotatable bonds and optimizes torsion angles via greedy scanning.
/// This approximates ETKDG's torsion knowledge terms by directly setting
/// preferred torsion angles after DG embedding.

/// A rotatable bond with its torsion atom indices and the atoms to rotate.
struct RotatableBond {
    /// The 4 atom indices defining the dihedral: a-b-c-d
    /// b-c is the rotatable bond
    dihedral: [usize; 4],
    /// Atom indices on the "d side" of the bond (to be rotated)
    mobile_atoms: Vec<usize>,
    /// Preferred torsion angles to try (radians)
    preferred_angles: Vec<f32>,
}

/// Find all rotatable bonds in the molecule and prepare rotation data.
fn find_rotatable_bonds(mol: &crate::graph::Molecule) -> Vec<RotatableBond> {
    let mut bonds = Vec::new();
    let n = mol.graph.node_count();

    for edge in mol.graph.edge_references() {
        let u = edge.source();
        let v = edge.target();

        // Only single bonds
        if edge.weight().order != crate::graph::BondOrder::Single {
            continue;
        }

        // Both must have >= 2 neighbors (not terminal)
        let deg_u = mol.graph.neighbors(u).count();
        let deg_v = mol.graph.neighbors(v).count();
        if deg_u < 2 || deg_v < 2 {
            continue;
        }

        // Skip ring bonds
        if crate::graph::min_path_excluding2(mol, u, v, u, v, 7).is_some() {
            continue;
        }

        // Pick first neighbor of u (not v) as "a" and first neighbor of v (not u) as "d"
        let a = mol.graph.neighbors(u).find(|&x| x != v).unwrap();
        let d = mol.graph.neighbors(v).find(|&x| x != u).unwrap();

        // BFS from v excluding u to find atoms on v's side (mobile atoms)
        let mut mobile = Vec::new();
        let mut visited = HashSet::new();
        visited.insert(u.index());
        visited.insert(v.index());
        let mut queue = VecDeque::new();
        for nb in mol.graph.neighbors(v) {
            if nb != u {
                queue.push_back(nb.index());
                visited.insert(nb.index());
            }
        }
        while let Some(curr) = queue.pop_front() {
            mobile.push(curr);
            let ni = petgraph::graph::NodeIndex::new(curr);
            for nb in mol.graph.neighbors(ni) {
                if !visited.contains(&nb.index()) {
                    visited.insert(nb.index());
                    queue.push_back(nb.index());
                }
            }
        }

        // If mobile side has more atoms than the other side, swap to rotate the smaller side
        let other_count = n - mobile.len() - 2; // exclude u, v
        if mobile.len() > other_count {
            // Rotate the u-side instead
            let mut mobile_u = Vec::new();
            let mut visited_u = HashSet::new();
            visited_u.insert(u.index());
            visited_u.insert(v.index());
            let mut queue_u = VecDeque::new();
            for nb in mol.graph.neighbors(u) {
                if nb != v {
                    queue_u.push_back(nb.index());
                    visited_u.insert(nb.index());
                }
            }
            while let Some(curr) = queue_u.pop_front() {
                mobile_u.push(curr);
                let ni = petgraph::graph::NodeIndex::new(curr);
                for nb in mol.graph.neighbors(ni) {
                    if !visited_u.contains(&nb.index()) {
                        visited_u.insert(nb.index());
                        queue_u.push_back(nb.index());
                    }
                }
            }
            // Use swapped dihedral: d-v-u-a, mobile is u-side
            let preferred = get_preferred_angles(mol, v, u);
            bonds.push(RotatableBond {
                dihedral: [d.index(), v.index(), u.index(), a.index()],
                mobile_atoms: mobile_u,
                preferred_angles: preferred,
            });
        } else {
            let preferred = get_preferred_angles(mol, u, v);
            bonds.push(RotatableBond {
                dihedral: [a.index(), u.index(), v.index(), d.index()],
                mobile_atoms: mobile,
                preferred_angles: preferred,
            });
        }
    }
    bonds
}

/// Get preferred torsion angles based on hybridization of the central bond atoms.
fn get_preferred_angles(
    mol: &crate::graph::Molecule,
    u: petgraph::graph::NodeIndex,
    v: petgraph::graph::NodeIndex,
) -> Vec<f32> {
    use crate::graph::Hybridization::*;
    use std::f32::consts::PI;

    let hyb_u = mol.graph[u].hybridization;
    let hyb_v = mol.graph[v].hybridization;

    match (hyb_u, hyb_v) {
        (SP3, SP3) => {
            // Anti and gauche: 60, 180, 300 degrees
            vec![PI / 3.0, PI, 5.0 * PI / 3.0]
        }
        (SP2, SP2) => {
            // Planar: 0 and 180
            vec![0.0, PI]
        }
        (SP2, SP3) | (SP3, SP2) => {
            // Mixed: 0, 60, 120, 180, 240, 300
            vec![0.0, PI / 3.0, 2.0 * PI / 3.0, PI, 4.0 * PI / 3.0, 5.0 * PI / 3.0]
        }
        _ => {
            // Generic: try all 30-degree increments
            (0..12).map(|i| i as f32 * PI / 6.0).collect()
        }
    }
}

/// Compute dihedral angle (in radians, range [-PI, PI]) from 4 atom positions in a DMatrix.
fn compute_dihedral(coords: &DMatrix<f32>, i: usize, j: usize, k: usize, l: usize) -> f32 {
    let b1 = nalgebra::Vector3::new(
        coords[(j, 0)] - coords[(i, 0)],
        coords[(j, 1)] - coords[(i, 1)],
        coords[(j, 2)] - coords[(i, 2)],
    );
    let b2 = nalgebra::Vector3::new(
        coords[(k, 0)] - coords[(j, 0)],
        coords[(k, 1)] - coords[(j, 1)],
        coords[(k, 2)] - coords[(j, 2)],
    );
    let b3 = nalgebra::Vector3::new(
        coords[(l, 0)] - coords[(k, 0)],
        coords[(l, 1)] - coords[(k, 1)],
        coords[(l, 2)] - coords[(k, 2)],
    );

    let n1 = b1.cross(&b2).normalize();
    let n2 = b2.cross(&b3).normalize();
    let m1 = n1.cross(&b2.normalize());
    let x = n1.dot(&n2);
    let y = m1.dot(&n2);
    y.atan2(x)
}

/// Rotate a set of atoms around the axis defined by points j→k by a given angle.
fn rotate_atoms(coords: &mut DMatrix<f32>, mobile: &[usize], j: usize, k: usize, angle: f32) {
    if angle.abs() < 1e-8 {
        return;
    }
    // Axis of rotation: from j to k
    let axis = nalgebra::Vector3::new(
        coords[(k, 0)] - coords[(j, 0)],
        coords[(k, 1)] - coords[(j, 1)],
        coords[(k, 2)] - coords[(j, 2)],
    );
    let axis_len = axis.norm();
    if axis_len < 1e-8 {
        return;
    }
    let axis = axis / axis_len;

    // Rodrigues' rotation formula
    let cos_a = angle.cos();
    let sin_a = angle.sin();

    // Pivot point is atom j
    let px = coords[(j, 0)];
    let py = coords[(j, 1)];
    let pz = coords[(j, 2)];

    for &idx in mobile {
        let vx = coords[(idx, 0)] - px;
        let vy = coords[(idx, 1)] - py;
        let vz = coords[(idx, 2)] - pz;

        let dot = axis[0] * vx + axis[1] * vy + axis[2] * vz;
        let cx = axis[1] * vz - axis[2] * vy;
        let cy = axis[2] * vx - axis[0] * vz;
        let cz = axis[0] * vy - axis[1] * vx;

        coords[(idx, 0)] = px + vx * cos_a + cx * sin_a + axis[0] * dot * (1.0 - cos_a);
        coords[(idx, 1)] = py + vy * cos_a + cy * sin_a + axis[1] * dot * (1.0 - cos_a);
        coords[(idx, 2)] = pz + vz * cos_a + cz * sin_a + axis[2] * dot * (1.0 - cos_a);
    }
}

/// Snap each rotatable bond's torsion to the nearest preferred angle.
/// No energy scoring — just geometric snapping.
pub fn snap_torsions_to_preferred(
    coords: &mut DMatrix<f32>,
    mol: &crate::graph::Molecule,
) -> usize {
    let rotatable = find_rotatable_bonds(mol);
    let num_rotatable = rotatable.len();

    for rb in &rotatable {
        let [a, b, c, d] = rb.dihedral;
        let current = compute_dihedral(coords, a, b, c, d);

        // Find nearest preferred angle
        let mut best_delta_abs = f32::MAX;
        let mut best_rotation = 0.0f32;
        for &target in &rb.preferred_angles {
            let mut delta = target - current;
            // Normalize to [-PI, PI]
            delta = (delta + std::f32::consts::PI).rem_euclid(2.0 * std::f32::consts::PI)
                - std::f32::consts::PI;
            if delta.abs() < best_delta_abs {
                best_delta_abs = delta.abs();
                best_rotation = delta;
            }
        }

        if best_rotation.abs() > 0.05 {
            rotate_atoms(coords, &rb.mobile_atoms, b, c, best_rotation);
        }
    }
    num_rotatable
}

/// Greedy torsion scan: for each rotatable bond, try preferred angles and pick
/// the one that minimizes the combined energy (bounds + torsion preferences).
/// Repeats for `passes` iterations.
pub fn optimize_torsions_greedy(
    coords: &mut DMatrix<f32>,
    mol: &crate::graph::Molecule,
    bounds: &DMatrix<f64>,
    passes: usize,
) -> usize {
    let rotatable = find_rotatable_bonds(mol);
    let num_rotatable = rotatable.len();
    if rotatable.is_empty() {
        return 0;
    }

    let params = super::energy::FFParams {
        kb: 300.0,
        k_theta: 200.0,
        k_omega: 10.0,
        k_oop: 20.0,
        k_bounds: 100.0,
        k_chiral: 0.0,
        k_vdw: 0.0,
    };

    for _pass in 0..passes {
        for rb in &rotatable {
            let [a, b, c, d] = rb.dihedral;
            let current_angle = compute_dihedral(coords, a, b, c, d);

            let current_energy = super::energy::calculate_total_energy(coords, mol, &params, bounds);
            let mut best_energy = current_energy;
            let mut best_rotation = 0.0f32;

            for &target_angle in &rb.preferred_angles {
                let delta = target_angle - current_angle;
                let delta = ((delta + std::f32::consts::PI) % (2.0 * std::f32::consts::PI))
                    - std::f32::consts::PI;

                rotate_atoms(coords, &rb.mobile_atoms, b, c, delta);

                let e = super::energy::calculate_total_energy(coords, mol, &params, bounds);
                if e < best_energy {
                    best_energy = e;
                    best_rotation = delta;
                }

                rotate_atoms(coords, &rb.mobile_atoms, b, c, -delta);
            }

            if best_rotation.abs() > 1e-6 {
                rotate_atoms(coords, &rb.mobile_atoms, b, c, best_rotation);
            }
        }
    }
    num_rotatable
}

/// Greedy torsion scan using ONLY bounds violation energy as criterion.
/// Better for matching DG-based conformations since it doesn't fight
/// the bounds constraints with torsion preferences.
pub fn optimize_torsions_bounds(
    coords: &mut DMatrix<f32>,
    mol: &crate::graph::Molecule,
    bounds: &DMatrix<f64>,
    passes: usize,
) -> usize {
    let rotatable = find_rotatable_bonds(mol);
    let num_rotatable = rotatable.len();
    if rotatable.is_empty() {
        return 0;
    }

    for _pass in 0..passes {
        for rb in &rotatable {
            let [a, b, c, d] = rb.dihedral;
            let current_angle = compute_dihedral(coords, a, b, c, d);

            let current_energy = super::bounds_ff::bounds_violation_energy(coords, bounds);
            let mut best_energy = current_energy;
            let mut best_rotation = 0.0f32;

            for &target_angle in &rb.preferred_angles {
                let delta = target_angle - current_angle;
                let delta = ((delta + std::f32::consts::PI) % (2.0 * std::f32::consts::PI))
                    - std::f32::consts::PI;

                rotate_atoms(coords, &rb.mobile_atoms, b, c, delta);

                let e = super::bounds_ff::bounds_violation_energy(coords, bounds);
                if e < best_energy {
                    best_energy = e;
                    best_rotation = delta;
                }

                rotate_atoms(coords, &rb.mobile_atoms, b, c, -delta);
            }

            if best_rotation.abs() > 1e-6 {
                rotate_atoms(coords, &rb.mobile_atoms, b, c, best_rotation);
            }
        }
    }
    num_rotatable
}
