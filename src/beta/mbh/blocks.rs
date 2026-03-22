//! MBH Rigid block detection — E9.1
//!
//! Identifies rigid groups in a molecule and decomposes it into
//! rigid blocks (6 DOF) and flexible atoms (3 DOF).

/// A rigid block (group of atoms treated as a rigid body).
#[derive(Debug, Clone)]
pub struct RigidBlock {
    /// Atom indices belonging to this block.
    pub atom_indices: Vec<usize>,
    /// Center of mass of the block.
    pub center_of_mass: [f64; 3],
    /// Whether the block is aromatic.
    pub is_aromatic: bool,
    /// Ring size (if from a ring).
    pub ring_size: usize,
}

/// Decomposition of a molecule into rigid blocks and flexible atoms.
#[derive(Debug, Clone)]
pub struct BlockDecomposition {
    /// Rigid blocks.
    pub blocks: Vec<RigidBlock>,
    /// Indices of flexible atoms (not in any rigid block).
    pub flexible_atoms: Vec<usize>,
    /// Total atoms in the molecule.
    pub n_atoms: usize,
    /// Total reduced DOF: 6 * n_blocks + 3 * n_flexible.
    pub n_dof_reduced: usize,
    /// Full Cartesian DOF: 3 * n_atoms.
    pub n_dof_full: usize,
}

/// Atomic masses for vibrational analysis.
fn atomic_mass(z: u8) -> f64 {
    match z {
        1 => 1.008,
        5 => 10.81,
        6 => 12.011,
        7 => 14.007,
        8 => 15.999,
        9 => 18.998,
        14 => 28.086,
        15 => 30.974,
        16 => 32.065,
        17 => 35.453,
        35 => 79.904,
        53 => 126.904,
        _ => 12.011, // default carbon
    }
}

/// Detect rigid blocks from ring information.
///
/// Uses ring data (atom indices, aromaticity) to identify rigid groups.
/// Atoms not in any ring are classified as flexible.
pub fn detect_rigid_blocks(
    n_atoms: usize,
    elements: &[u8],
    positions: &[[f64; 3]],
    rings: &[(Vec<usize>, bool)], // (atom_indices, is_aromatic)
) -> BlockDecomposition {
    let mut atom_in_block = vec![false; n_atoms];
    let mut blocks = Vec::new();

    for (ring_atoms, is_aromatic) in rings {
        // Filter: only small rings (≤ 8 atoms) are considered rigid
        if ring_atoms.len() > 8 {
            continue;
        }
        let overlap = ring_atoms.iter().any(|&a| atom_in_block[a]);
        if overlap {
            // For fused rings, merge into existing block or skip
            // Simple approach: skip if overlapping
            continue;
        }

        // Compute center of mass
        let mut com = [0.0; 3];
        let mut total_mass = 0.0;
        for &a in ring_atoms {
            let m = atomic_mass(elements[a]);
            for d in 0..3 {
                com[d] += m * positions[a][d];
            }
            total_mass += m;
        }
        if total_mass > 0.0 {
            for d in 0..3 {
                com[d] /= total_mass;
            }
        }

        for &a in ring_atoms {
            atom_in_block[a] = true;
        }

        blocks.push(RigidBlock {
            atom_indices: ring_atoms.clone(),
            center_of_mass: com,
            is_aromatic: *is_aromatic,
            ring_size: ring_atoms.len(),
        });
    }

    let flexible_atoms: Vec<usize> = (0..n_atoms).filter(|&i| !atom_in_block[i]).collect();
    let n_dof_reduced = 6 * blocks.len() + 3 * flexible_atoms.len();
    let n_dof_full = 3 * n_atoms;

    BlockDecomposition {
        blocks,
        flexible_atoms,
        n_atoms,
        n_dof_reduced,
        n_dof_full,
    }
}

/// Build projection matrix L that maps reduced DOF → Cartesian DOF.
///
/// L has shape (3*n_atoms, n_dof_reduced).
/// Columns correspond to:
///   - For each rigid block: 3 translation + 3 rotation vectors
///   - For each flexible atom: 3 Cartesian displacement vectors
pub fn build_projection_matrix(
    decomposition: &BlockDecomposition,
    positions: &[[f64; 3]],
    elements: &[u8],
) -> Vec<Vec<f64>> {
    let n_full = decomposition.n_dof_full;
    let n_reduced = decomposition.n_dof_reduced;
    let mut l_cols: Vec<Vec<f64>> = Vec::with_capacity(n_reduced);

    // For each rigid block: 6 columns (3 translations + 3 rotations)
    for block in &decomposition.blocks {
        let com = block.center_of_mass;

        // 3 translation vectors
        for d in 0..3 {
            let mut col = vec![0.0; n_full];
            let mut norm_sq = 0.0;
            for &a in &block.atom_indices {
                let m = atomic_mass(elements[a]).sqrt();
                col[a * 3 + d] = m;
                norm_sq += m * m;
            }
            // Normalize
            let norm = norm_sq.sqrt();
            if norm > 1e-15 {
                for v in col.iter_mut() {
                    *v /= norm;
                }
            }
            l_cols.push(col);
        }

        // 3 rotation vectors (around COM)
        // Rx: rotation around x-axis → displacement in y,z
        for rot_axis in 0..3 {
            let mut col = vec![0.0; n_full];
            let mut norm_sq = 0.0;

            for &a in &block.atom_indices {
                let m = atomic_mass(elements[a]).sqrt();
                let r = [
                    positions[a][0] - com[0],
                    positions[a][1] - com[1],
                    positions[a][2] - com[2],
                ];

                // Cross product of rotation axis unit vector with r
                let disp = match rot_axis {
                    0 => [0.0, -r[2], r[1]], // x-axis: (0,0,1)×r → (0,-rz,ry) wait, e_x × r = (0*rz - 1*ry, ...) = wrong
                    // e_x × r = (0, r[2], -r[1]) NO...
                    // Actually: e_x = (1,0,0), r = (rx,ry,rz)
                    // e_x × r = |i  j  k | = (0·rz - 0·ry, 0·rx - 1·rz, 1·ry - 0·rx) = (0, -rz, ry)
                    1 => [r[2], 0.0, -r[0]], // e_y × r = (rz, 0, -rx)
                    2 => [-r[1], r[0], 0.0], // e_z × r = (-ry, rx, 0)
                    _ => unreachable!(),
                };

                for d in 0..3 {
                    col[a * 3 + d] = m * disp[d];
                    norm_sq += (m * disp[d]) * (m * disp[d]);
                }
            }

            let norm = norm_sq.sqrt();
            if norm > 1e-15 {
                for v in col.iter_mut() {
                    *v /= norm;
                }
            }
            l_cols.push(col);
        }
    }

    // For each flexible atom: 3 Cartesian displacement vectors
    for &a in &decomposition.flexible_atoms {
        for d in 0..3 {
            let mut col = vec![0.0; n_full];
            col[a * 3 + d] = 1.0;
            l_cols.push(col);
        }
    }

    l_cols
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_blocks_benzene() {
        // 6 ring carbons + 6 hydrogens
        let n = 12;
        let elements: Vec<u8> = vec![6; 6].into_iter().chain(vec![1; 6]).collect();
        let positions: Vec<[f64; 3]> = (0..12)
            .map(|i| {
                if i < 6 {
                    let angle = i as f64 * std::f64::consts::PI / 3.0;
                    [1.4 * angle.cos(), 1.4 * angle.sin(), 0.0]
                } else {
                    let angle = (i - 6) as f64 * std::f64::consts::PI / 3.0;
                    [2.5 * angle.cos(), 2.5 * angle.sin(), 0.0]
                }
            })
            .collect();

        let rings = vec![(vec![0, 1, 2, 3, 4, 5], true)]; // benzene ring
        let decomp = detect_rigid_blocks(n, &elements, &positions, &rings);

        assert_eq!(decomp.blocks.len(), 1);
        assert_eq!(decomp.blocks[0].atom_indices.len(), 6);
        assert!(decomp.blocks[0].is_aromatic);
        assert_eq!(decomp.flexible_atoms.len(), 6); // 6 hydrogens
        assert_eq!(decomp.n_dof_reduced, 6 * 1 + 3 * 6); // 6 + 18 = 24
        assert_eq!(decomp.n_dof_full, 36);
    }

    #[test]
    fn test_detect_blocks_no_rings() {
        let n = 3;
        let elements = vec![6, 1, 1];
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let rings: Vec<(Vec<usize>, bool)> = vec![];
        let decomp = detect_rigid_blocks(n, &elements, &positions, &rings);

        assert_eq!(decomp.blocks.len(), 0);
        assert_eq!(decomp.flexible_atoms.len(), 3);
        assert_eq!(decomp.n_dof_reduced, 9);
    }

    #[test]
    fn test_projection_matrix_dimensions() {
        let n = 9; // 6 ring + 3 flexible
        let elements = vec![6; 9];
        let positions: Vec<[f64; 3]> = (0..9)
            .map(|i| {
                if i < 6 {
                    let a = i as f64 * std::f64::consts::PI / 3.0;
                    [1.4 * a.cos(), 1.4 * a.sin(), 0.0]
                } else {
                    [(i as f64) * 2.0, 0.0, 0.0]
                }
            })
            .collect();

        let rings = vec![(vec![0, 1, 2, 3, 4, 5], true)];
        let decomp = detect_rigid_blocks(n, &elements, &positions, &rings);

        let l = build_projection_matrix(&decomp, &positions, &elements);
        assert_eq!(l.len(), decomp.n_dof_reduced); // number of columns
        for col in &l {
            assert_eq!(col.len(), decomp.n_dof_full); // each column has 3*n_atoms entries
        }
    }

    #[test]
    fn test_projection_columns_normalized() {
        let n = 6;
        let elements = vec![6; 6];
        let positions: Vec<[f64; 3]> = (0..6)
            .map(|i| {
                let a = i as f64 * std::f64::consts::PI / 3.0;
                [1.4 * a.cos(), 1.4 * a.sin(), 0.0]
            })
            .collect();

        let rings = vec![(vec![0, 1, 2, 3, 4, 5], true)];
        let decomp = detect_rigid_blocks(n, &elements, &positions, &rings);
        let l = build_projection_matrix(&decomp, &positions, &elements);

        for (i, col) in l.iter().enumerate() {
            let norm: f64 = col.iter().map(|v| v * v).sum::<f64>().sqrt();
            assert!(
                (norm - 1.0).abs() < 0.1 || norm < 1e-10,
                "Column {} norm should be ~1.0 or 0: {}",
                i,
                norm
            );
        }
    }
}
