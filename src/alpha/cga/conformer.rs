//! CGA-based conformer operations — E1.2
//!
//! Provides dihedral torsion via CGA Motors, subtree rotation,
//! and a CGA refinement pipeline for torsion scans.

use super::motor::{conformal_point, extract_euclidean, Motor};
use super::multivector::Multivector;

/// Lift a 3D point to a CGA conformal null vector.
pub fn embed_point(p: [f64; 3]) -> Multivector {
    conformal_point(p)
}

/// Extract 3D coordinates from a CGA conformal point.
pub fn extract_point(p: &Multivector) -> [f64; 3] {
    extract_euclidean(p)
}

/// Build a Motor that rotates by `angle` radians about the bond axis A→B.
///
/// The rotation axis is the unit vector from atom A to atom B.
/// The translation is chosen so that point A remains fixed
/// (rotate about the line through A along A→B).
pub fn dihedral_motor(a: [f64; 3], b: [f64; 3], angle: f64) -> Motor {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let dz = b[2] - a[2];
    let len = (dx * dx + dy * dy + dz * dz).sqrt();
    if len < 1e-15 {
        return Motor::identity();
    }
    let axis = [dx / len, dy / len, dz / len];

    // To rotate about a line through A (not the origin):
    //   1. Translate so A is at origin: T(-A)
    //   2. Rotate by angle about axis: R
    //   3. Translate back: T(A)
    // Composite: M = T(A) R T(-A)
    let t_neg = Motor::translator([-a[0], -a[1], -a[2]]);
    let r = Motor::rotor(axis, angle);
    let t_pos = Motor::translator(a);
    t_pos.compose(&r.compose(&t_neg))
}

/// Apply a Motor to a subset of atoms in a coordinate array.
///
/// `coords`: flat coordinate array `[x0,y0,z0, x1,y1,z1, ...]`
/// `indices`: atom indices in the subtree to be rotated
/// `motor`: the CGA Motor to apply
///
/// Returns modified coordinate array (input is not mutated).
pub fn apply_motor_to_subtree(coords: &[f64], indices: &[usize], motor: &Motor) -> Vec<f64> {
    let mut result = coords.to_vec();
    for &idx in indices {
        let base = idx * 3;
        let p = [coords[base], coords[base + 1], coords[base + 2]];
        let tp = motor.transform_point(p);
        result[base] = tp[0];
        result[base + 1] = tp[1];
        result[base + 2] = tp[2];
    }
    result
}

/// Find the subtree of atoms to rotate about a rotatable bond.
///
/// `bond_a`, `bond_b`: the two atoms defining the rotatable bond.
/// `adjacency`: adjacency list (for each atom, list of bonded atom indices).
///
/// Returns the set of atom indices on the `bond_b` side (excluding `bond_a`).
pub fn find_rotation_subtree(
    bond_a: usize,
    bond_b: usize,
    adjacency: &[Vec<usize>],
) -> Vec<usize> {
    let n = adjacency.len();
    let mut visited = vec![false; n];
    visited[bond_a] = true; // block traversal through bond_a
    visited[bond_b] = true;

    let mut stack = vec![bond_b];
    let mut subtree = vec![bond_b];

    while let Some(current) = stack.pop() {
        for &neighbor in &adjacency[current] {
            if !visited[neighbor] {
                visited[neighbor] = true;
                subtree.push(neighbor);
                stack.push(neighbor);
            }
        }
    }

    subtree
}

/// Perform a torsion scan using CGA Motors.
///
/// Rotates the subtree about bond A→B through `n_steps` equally-spaced angles
/// from 0 to 2π.  Returns a vector of coordinate snapshots.
///
/// `coords`: flat coordinate array
/// `bond_a`, `bond_b`: atoms defining the rotatable bond
/// `subtree`: atom indices to rotate
/// `n_steps`: number of angular steps
pub fn refine_torsion_cga(
    coords: &[f64],
    bond_a: usize,
    bond_b: usize,
    subtree: &[usize],
    n_steps: usize,
) -> Vec<Vec<f64>> {
    let a = [
        coords[bond_a * 3],
        coords[bond_a * 3 + 1],
        coords[bond_a * 3 + 2],
    ];
    let b = [
        coords[bond_b * 3],
        coords[bond_b * 3 + 1],
        coords[bond_b * 3 + 2],
    ];

    let mut snapshots = Vec::with_capacity(n_steps);
    for step in 0..n_steps {
        let angle = 2.0 * std::f64::consts::PI * step as f64 / n_steps as f64;
        let motor = dihedral_motor(a, b, angle);
        let new_coords = apply_motor_to_subtree(coords, subtree, &motor);
        snapshots.push(new_coords);
    }
    snapshots
}

/// Compute the dihedral angle defined by four points (in radians, range -π..π).
pub fn dihedral_angle(p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], p4: [f64; 3]) -> f64 {
    let b1 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
    let b2 = [p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]];
    let b3 = [p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2]];

    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let m1 = cross(n1, b2_unit(b2));

    let x = dot(n1, n2);
    let y = dot(m1, n2);
    (-y).atan2(-x) + std::f64::consts::PI
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn b2_unit(b: [f64; 3]) -> [f64; 3] {
    let len = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();
    if len < 1e-15 {
        return [0.0, 0.0, 0.0];
    }
    [b[0] / len, b[1] / len, b[2] / len]
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn test_embed_extract_round_trip() {
        let points = [
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 3.0],
            [-5.678, 9.012, -3.456],
            [0.001, -0.001, 100.0],
        ];
        for p in &points {
            let cp = embed_point(*p);
            let ep = extract_point(&cp);
            for i in 0..3 {
                assert!(
                    (ep[i] - p[i]).abs() < 1e-10,
                    "Round-trip error: {:?} → {:?}",
                    p, ep
                );
            }
        }
    }

    #[test]
    fn test_dihedral_motor_preserves_bond() {
        // Bond axis from (0,0,0) to (1,0,0), rotate subtree
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let motor = dihedral_motor(a, b, FRAC_PI_2);

        // Point A should not move
        let ta = motor.transform_point(a);
        for i in 0..3 {
            assert!(
                (ta[i] - a[i]).abs() < 1e-10,
                "Bond atom A moved: {:?} → {:?}",
                a, ta
            );
        }

        // Point B should not move (it's on the rotation axis)
        let tb = motor.transform_point(b);
        for i in 0..3 {
            assert!(
                (tb[i] - b[i]).abs() < 1e-10,
                "Bond atom B moved: {:?} → {:?}",
                b, tb
            );
        }
    }

    #[test]
    fn test_dihedral_motor_rotates_subtree() {
        // Bond along X: A=(0,0,0) B=(1,0,0)
        // Point at (2, 1, 0): after 90° about X, should go to (2, 0, 1)
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let motor = dihedral_motor(a, b, FRAC_PI_2);
        let result = motor.transform_point([2.0, 1.0, 0.0]);
        assert!((result[0] - 2.0).abs() < 1e-10, "x = {}", result[0]);
        assert!((result[1] - 0.0).abs() < 1e-10, "y = {}", result[1]);
        assert!((result[2] - 1.0).abs() < 1e-10, "z = {}", result[2]);
    }

    #[test]
    fn test_apply_motor_to_subtree() {
        // 3 atoms: 0=(0,0,0), 1=(1,0,0), 2=(2,1,0)
        // Bond 0→1, rotate atom 2 by 90° about X
        let coords = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 1.0, 0.0];
        let motor = dihedral_motor([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], FRAC_PI_2);
        let result = apply_motor_to_subtree(&coords, &[2], &motor);

        // Atom 0 unchanged
        assert!((result[0] - 0.0).abs() < 1e-10);
        assert!((result[1] - 0.0).abs() < 1e-10);
        assert!((result[2] - 0.0).abs() < 1e-10);

        // Atom 1 unchanged
        assert!((result[3] - 1.0).abs() < 1e-10);
        assert!((result[4] - 0.0).abs() < 1e-10);
        assert!((result[5] - 0.0).abs() < 1e-10);

        // Atom 2 rotated: (2, 1, 0) → (2, 0, 1)
        assert!((result[6] - 2.0).abs() < 1e-10, "x = {}", result[6]);
        assert!((result[7] - 0.0).abs() < 1e-10, "y = {}", result[7]);
        assert!((result[8] - 1.0).abs() < 1e-10, "z = {}", result[8]);
    }

    #[test]
    fn test_no_gimbal_lock_at_90() {
        // Verify no gimbal lock at ω ≈ ±90°
        let a = [0.0, 0.0, 0.0];
        let b = [0.0, 0.0, 1.0]; // Z axis bond
        let p = [1.0, 0.0, 1.5];

        for &angle in &[FRAC_PI_2, -FRAC_PI_2, PI, -PI] {
            let motor = dihedral_motor(a, b, angle);
            let result = motor.transform_point(p);
            // Distance from axis should be preserved
            let dist_orig = (p[0] * p[0] + p[1] * p[1]).sqrt();
            let dist_new = (result[0] * result[0] + result[1] * result[1]).sqrt();
            assert!(
                (dist_orig - dist_new).abs() < 1e-10,
                "Distance from axis changed at angle {}: {} → {}",
                angle, dist_orig, dist_new
            );
            // Z coordinate preserved
            assert!(
                (result[2] - p[2]).abs() < 1e-10,
                "Z changed at angle {}: {} → {}",
                angle, p[2], result[2]
            );
        }
    }

    #[test]
    fn test_torsion_scan_produces_correct_count() {
        let coords = vec![
            0.0, 0.0, 0.0, // atom 0
            1.0, 0.0, 0.0, // atom 1
            2.0, 1.0, 0.0, // atom 2
        ];
        let snapshots = refine_torsion_cga(&coords, 0, 1, &[2], 12);
        assert_eq!(snapshots.len(), 12);
    }

    #[test]
    fn test_find_rotation_subtree() {
        // Linear chain: 0 - 1 - 2 - 3
        let adj = vec![
            vec![1],       // atom 0
            vec![0, 2],    // atom 1
            vec![1, 3],    // atom 2
            vec![2],       // atom 3
        ];
        // Rotate about bond 1→2, subtree on side of 2
        let subtree = find_rotation_subtree(1, 2, &adj);
        assert!(subtree.contains(&2));
        assert!(subtree.contains(&3));
        assert!(!subtree.contains(&0));
        assert!(!subtree.contains(&1));
    }
}
