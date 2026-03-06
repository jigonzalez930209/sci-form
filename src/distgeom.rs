use crate::graph::Molecule;
use nalgebra::DMatrix;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;

fn compute_14_dist_cis(r12: f32, r23: f32, r34: f32, angle123: f32, angle234: f32) -> f32 {
    let x0 = r12 * angle123.cos();
    let y0 = r12 * angle123.sin();

    let x3 = r23 - r34 * angle234.cos();
    let y3 = r34 * angle234.sin();

    let dx = x3 - x0;
    let dy = y3 - y0;
    (dx * dx + dy * dy).sqrt()
}

fn compute_14_dist_trans(r12: f32, r23: f32, r34: f32, angle123: f32, angle234: f32) -> f32 {
    let x0 = r12 * angle123.cos();
    let y0 = r12 * angle123.sin();

    let x3 = r23 - r34 * angle234.cos();
    let y3 = -r34 * angle234.sin();

    let dx = x3 - x0;
    let dy = y3 - y0;
    (dx * dx + dy * dy).sqrt()
}

fn get_bond_length(mol: &Molecule, n1: NodeIndex, n2: NodeIndex) -> f32 {
    let a1 = &mol.graph[n1];
    let a2 = &mol.graph[n2];

    let mut e1 = a1.element;
    let mut e2 = a2.element;
    if e1 > e2 {
        std::mem::swap(&mut e1, &mut e2);
    }

    let mut order = crate::graph::BondOrder::Single;
    if let Some(edge_idx) = mol.graph.find_edge(n1, n2) {
        order = mol.graph[edge_idx].order.clone();
    }

    match (e1, e2, &order) {
        (1, 6, crate::graph::BondOrder::Single) => return 1.093,
        (1, 7, crate::graph::BondOrder::Single) => return 1.016,
        (1, 8, crate::graph::BondOrder::Single) => return 0.974,
        (6, 6, crate::graph::BondOrder::Single) => return 1.512,
        (6, 6, crate::graph::BondOrder::Double) => return 1.332,
        (6, 6, crate::graph::BondOrder::Triple) => return 1.201,
        (6, 6, crate::graph::BondOrder::Aromatic) => return 1.386,
        (6, 7, crate::graph::BondOrder::Single) => return 1.422,
        (6, 7, crate::graph::BondOrder::Double) => return 1.287,
        (6, 7, crate::graph::BondOrder::Triple) => return 1.161,
        (6, 7, crate::graph::BondOrder::Aromatic) => return 1.349,
        (6, 8, crate::graph::BondOrder::Single) => return 1.417,
        (6, 8, crate::graph::BondOrder::Double) => return 1.224,
        (6, 8, crate::graph::BondOrder::Aromatic) => return 1.361,
        (7, 7, crate::graph::BondOrder::Aromatic) => return 1.344,
        _ => {}
    }

    let mut dist = crate::graph::get_covalent_radius(e1) + crate::graph::get_covalent_radius(e2);

    match order {
        crate::graph::BondOrder::Double => dist -= 0.20,
        crate::graph::BondOrder::Triple => dist -= 0.34,
        crate::graph::BondOrder::Aromatic => dist -= 0.15,
        _ => {}
    }
    dist
}

/// Calculates the initial bounds matrix for a given molecule using RDKit's ETKDG approach.
pub fn calculate_bounds_matrix(mol: &Molecule) -> DMatrix<f32> {
    let n = mol.graph.node_count();
    let mut bounds = DMatrix::from_element(n, n, 0.0);

    let max_dist = 1000.0;
    for i in 0..n {
        for j in 0..n {
            if i != j {
                if i < j {
                    bounds[(i, j)] = max_dist; // Upper triangle
                } else {
                    bounds[(i, j)] = 0.0; // Lower triangle
                }
            }
        }
    }

    // Unweighted topological distances using Floyd-Warshall
    let mut top_dist = DMatrix::from_element(n, n, 1000);
    for i in 0..n {
        top_dist[(i, i)] = 0;
    }
    for bond in mol.graph.edge_references() {
        let i = bond.source().index();
        let j = bond.target().index();
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

    let mut set_bounds = |i: usize, j: usize, lower: f32, upper: f32| {
        let (min_idx, max_idx) = if i < j { (i, j) } else { (j, i) };
        if bounds[(max_idx, min_idx)] < lower {
            bounds[(max_idx, min_idx)] = lower;
        }
        if bounds[(min_idx, max_idx)] > upper {
            bounds[(min_idx, max_idx)] = upper;
        }
    };

    // 1-2 Bounds
    for bond in mol.graph.edge_references() {
        let i = bond.source();
        let j = bond.target();

        let ideal_dist = get_bond_length(mol, i, j);
        set_bounds(i.index(), j.index(), ideal_dist - 0.01, ideal_dist + 0.01);
        top_dist[(i.index(), j.index())] = 1; // Ensure it's marked
        top_dist[(j.index(), i.index())] = 1;
    }

    // 1-3 Bounds
    for i in 0..n {
        let neighbors: Vec<_> = mol.graph.neighbors(NodeIndex::new(i)).collect();

        for j in 0..neighbors.len() {
            for k in (j + 1)..neighbors.len() {
                let n1 = neighbors[j];
                let n2 = neighbors[k];

                let ideal_angle =
                    crate::graph::get_corrected_ideal_angle(mol, NodeIndex::new(i), n1, n2);

                let d1 = get_bond_length(mol, NodeIndex::new(i), n1);
                let d2 = get_bond_length(mol, NodeIndex::new(i), n2);

                let ideal_dist_sq = d1 * d1 + d2 * d2 - 2.0 * d1 * d2 * ideal_angle.cos();
                let ideal_dist = ideal_dist_sq.sqrt();

                let dist = crate::graph::min_path_excluding(mol, n1, n2, NodeIndex::new(i), 3);
                
                let (lower, upper) = if dist == Some(1) {
                    // 3-membered ring case: use the direct 1-2 bond length
                    let d = get_bond_length(mol, n1, n2);
                    (d - 0.01, d + 0.01)
                } else {
                    let d1 = get_bond_length(mol, NodeIndex::new(i), n1);
                    let d2 = get_bond_length(mol, NodeIndex::new(i), n2);
                    let ideal_dist_sq = d1 * d1 + d2 * d2 - 2.0 * d1 * d2 * ideal_angle.cos();
                    let ideal_dist = ideal_dist_sq.sqrt();
                    let (tol_l, tol_u) = match dist {
                        Some(2) => (0.02, 0.02), // 4-ring
                        _ => (0.04, 0.04), // Default
                    };
                    (ideal_dist - tol_l, ideal_dist + tol_u)
                };

                set_bounds(n1.index(), n2.index(), lower, upper);
                top_dist[(n1.index(), n2.index())] = 2; // Mark as 1-3
                top_dist[(n2.index(), n1.index())] = 2;
            }
        }
    }

    // 1-4 Bounds (Dihedrals / Torsions)
    for i in 0..n {
        for j in (i + 1)..n {
            if top_dist[(i, j)] == 3 {
                let mut path_found = false;
                for n1 in mol.graph.neighbors(NodeIndex::new(i)) {
                    for n2 in mol.graph.neighbors(NodeIndex::new(j)) {
                        if mol.graph.contains_edge(n1, n2) {
                            let r12 = get_bond_length(mol, NodeIndex::new(i), n1);
                            let r23 = get_bond_length(mol, n1, n2);
                            let r34 = get_bond_length(mol, n2, NodeIndex::new(j));

                            let angle123 = crate::graph::get_corrected_ideal_angle(mol, n1, NodeIndex::new(i), n2);
                            let angle234 = crate::graph::get_corrected_ideal_angle(mol, n2, n1, NodeIndex::new(j));

                            let d_cis = compute_14_dist_cis(r12, r23, r34, angle123, angle234);
                            let d_trans = compute_14_dist_trans(r12, r23, r34, angle123, angle234);

                            let edge = mol.graph.find_edge(n1, n2).unwrap();
                            let is_cis = mol.graph[edge].stereo == crate::graph::BondStereo::Z;
                            let is_trans = mol.graph[edge].stereo == crate::graph::BondStereo::E;

                            let (mut lower, upper) = if is_cis {
                                (d_cis - 0.06, d_cis + 0.06)
                            } else if is_trans {
                                (d_trans - 0.06, d_trans + 0.06)
                            } else {
                                let min_d = d_cis.min(d_trans);
                                let max_d = d_cis.max(d_trans);
                                if (max_d - min_d) < 0.01 {
                                    (min_d - 0.06, max_d + 0.06)
                                } else {
                                    (min_d, max_d)
                                }
                            };

                            if lower < 0.0 { lower = 0.0; }
                            set_bounds(i, j, lower, upper);
                            path_found = true;
                            break;
                        }
                    }
                    if path_found { break; }
                }
            }
        }
    }

    // VDW Bounds (>= 1-5)
    for i in 0..n {
        for j in (i + 1)..n {
            let dist = top_dist[(i, j)];
            if dist >= 3 { // RDKit: 1-4 and above
                let atom_i = &mol.graph[NodeIndex::new(i)];
                let atom_j = &mol.graph[NodeIndex::new(j)];
                let vdw_i = crate::graph::get_vdw_radius(atom_i.element);
                let vdw_j = crate::graph::get_vdw_radius(atom_j.element);
                let vdw_sum = vdw_i + vdw_j;

                // RDKit exact: VDW_SCALE_15=0.7 for 1-5, 0.85 for 1-6, 1.0 for 1-7+
                let lower = if dist == 3 {
                    vdw_sum * 0.5 // 1-4
                } else if dist == 4 {
                    vdw_sum * 0.7 // 1-5
                } else if dist == 5 {
                    vdw_sum * 0.85 // 1-6
                } else {
                    vdw_sum
                };

                set_bounds(i, j, lower, 1000.0);
            }
        }
    }

    bounds
}

/// Applies Triangle Inequality Smoothing using Floyd-Warshall to update the limits.
/// Converges identically to RDKit's TriangleSmooth implementation.
pub fn triangle_smooth(bounds: &mut DMatrix<f32>) -> bool {
    let n = bounds.nrows();
    for k in 0..n {
        for i in 0..(n - 1) {
            if i == k { continue; }
            
            let u_ik = if i < k { bounds[(i, k)] } else { bounds[(k, i)] };
            let l_ik = if i < k { bounds[(k, i)] } else { bounds[(i, k)] };

            for j in (i + 1)..n {
                if j == k { continue; }
                
                let u_kj = if j < k { bounds[(j, k)] } else { bounds[(k, j)] };
                let l_kj = if j < k { bounds[(k, j)] } else { bounds[(j, k)] };

                // Upper bound: U_ij = min(U_ij, U_ik + U_kj)
                let sum_u = u_ik + u_kj;
                if bounds[(i, j)] > sum_u {
                    bounds[(i, j)] = sum_u;
                }

                // Lower bound: L_ij = max(L_ij, L_ik - U_kj, L_kj - U_ik)
                let diff1 = l_ik - u_kj;
                let diff2 = l_kj - u_ik;
                let mut l_ij = bounds[(j, i)];
                if l_ij < diff1 {
                    l_ij = diff1;
                }
                if l_ij < diff2 {
                    l_ij = diff2;
                }
                bounds[(j, i)] = l_ij;

                if bounds[(j, i)] > bounds[(i, j)] {
                    // Slight violation can happen due to float precision, 
                    // but if it's significant, the bounds are inconsistent.
                    if bounds[(j, i)] - bounds[(i, j)] > 1e-3 {
                        return false; 
                    }
                    bounds[(j, i)] = bounds[(i, j)];
                }
            }
        }
    }
    true
}

use rand::Rng;

pub fn pick_random_distances<R: Rng>(rng: &mut R, bounds: &DMatrix<f32>) -> DMatrix<f32> {
    let n = bounds.nrows();
    let mut dists = DMatrix::from_element(n, n, 0.0);

    for i in 0..n {
        for j in (i + 1)..n {
            let max_val = bounds[(i, j)]; // upper tri: max
            let min_val = bounds[(j, i)]; // lower tri: min

            // Randomly select between min and max uniformly
            let chosen = if max_val > min_val {
                rng.gen_range(min_val..=max_val)
            } else {
                min_val // fallback
            };

            dists[(i, j)] = chosen;
            dists[(j, i)] = chosen; // symmetric matrix
        }
    }

    dists
}

/// Computes Metric Matrix for distance geometry
/// formula: M_ij = (D_i0^2 + D_j0^2 - D_ij^2) / 2
/// We use the centroid as origin. For simplicity we can compute squared distances
/// and center the matrix via the Gram matrix method: M = -1/2 * P * D^2 * P
pub fn compute_metric_matrix(dists: &DMatrix<f32>) -> DMatrix<f32> {
    let n = dists.nrows();
    let mut sq_dists = DMatrix::from_element(n, n, 0.0);

    // Square the distances
    for i in 0..n {
        for j in 0..n {
            sq_dists[(i, j)] = dists[(i, j)].powi(2);
        }
    }

    // Centering Matrix P = I - 1/n * J
    let p =
        DMatrix::from_diagonal_element(n, n, 1.0) - DMatrix::from_element(n, n, 1.0 / (n as f32));

    // M = -0.5 * P * D^2 * P
    -0.5 * (&p * &sq_dists * &p)
}

use nalgebra::SymmetricEigen;

/// Uses nalgebra's Eigen decomposition to project the metric matrix into N dimensions
pub fn generate_nd_coordinates<R: Rng>(rng: &mut R, metric_matrix: &DMatrix<f32>, ndim: usize) -> DMatrix<f32> {
    let n = metric_matrix.nrows();

    // Decompose the symmetric matrix M
    let eigen = SymmetricEigen::new(metric_matrix.clone());

    // Sort eigenvalues and corresponding eigenvectors in descending order
    let mut evals: Vec<(usize, f32)> = eigen.eigenvalues.iter().copied().enumerate().collect();
    evals.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    let mut coords = DMatrix::from_element(n, ndim, 0.0);

    for dim in 0..ndim {
        // We only use the top N eigenvalues
        if dim < evals.len() {
            let (idx, eval) = evals[dim];
            if eval > 0.0 {
                let root_eval = eval.sqrt();
                let evec = eigen.eigenvectors.column(idx);
                for i in 0..n {
                    coords[(i, dim)] = evec[i] * root_eval;
                }
            } else {
                // Initialize with very small noise to break symmetry, but not distort structure
                for i in 0..n {
                    coords[(i, dim)] = (1.0 - 2.0 * rng.gen::<f32>()) * 1e-4;
                }
            }
        } else {
             for i in 0..n {
                coords[(i, dim)] = (1.0 - 2.0 * rng.gen::<f32>()) * 1e-4;
            }
        }
    }

    coords
}

/// Legacy wrapper for 3D
pub fn generate_3d_coordinates<R: Rng>(rng: &mut R, metric_matrix: &DMatrix<f32>) -> DMatrix<f32> {
    generate_nd_coordinates(rng, metric_matrix, 3)
}

/// Calculates the chiral volume for 4 points.
/// This matches RDKit's approach: V = v1 . (v2 x v3)
/// where v_n is the vector from the central atom (idx4 by convention or simply the origin of the chiral center)
/// to the n-th neighbor.
pub fn calc_chiral_volume(
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize, // usually the chiral center in RDKit's chiral volume convention
    coords: &DMatrix<f32>,
) -> f32 {
    let dim = coords.ncols();
    assert!(dim >= 3, "Coordinates must have at least 3 dimensions");

    // v1 = pos1 - pos4
    let v1 = nalgebra::Vector3::new(
        coords[(idx1, 0)] - coords[(idx4, 0)],
        coords[(idx1, 1)] - coords[(idx4, 1)],
        coords[(idx1, 2)] - coords[(idx4, 2)],
    );

    // v2 = pos2 - pos4
    let v2 = nalgebra::Vector3::new(
        coords[(idx2, 0)] - coords[(idx4, 0)],
        coords[(idx2, 1)] - coords[(idx4, 1)],
        coords[(idx2, 2)] - coords[(idx4, 2)],
    );

    // v3 = pos3 - pos4
    let v3 = nalgebra::Vector3::new(
        coords[(idx3, 0)] - coords[(idx4, 0)],
        coords[(idx3, 1)] - coords[(idx4, 1)],
        coords[(idx3, 2)] - coords[(idx4, 2)],
    );

    // Triple scalar product: v1 . (v2 x v3)
    let v2_cross_v3 = v2.cross(&v3);
    v1.dot(&v2_cross_v3)
}

/// Identifies chiral centers and their neighbor sets for volume enforcement.
pub fn identify_chiral_sets(mol: &Molecule) -> Vec<crate::forcefield::bounds_ff::ChiralSet> {
    let mut sets = Vec::new();
    for i in 0..mol.graph.node_count() {
        let ni = petgraph::graph::NodeIndex::new(i);
        let atom = &mol.graph[ni];
        if atom.chiral_tag != crate::graph::ChiralType::Unspecified {
            let mut neighbors: Vec<_> = mol.graph.neighbors(ni).map(|n| n.index()).collect();
            if neighbors.len() == 4 {
                neighbors.sort();
                // RDKit convention: volume sign depends on neighbors permutation.
                // For CW/CCW we just need a non-zero range.
                // RDKit typically uses [2.0, 50.0] or [-50.0, -2.0]
                let (lower, upper) = match atom.chiral_tag {
                    crate::graph::ChiralType::TetrahedralCW => (2.0, 50.0),
                    crate::graph::ChiralType::TetrahedralCCW => (-50.0, -2.0),
                    _ => (0.0, 0.0),
                };
                if lower != 0.0 || upper != 0.0 {
                    sets.push(crate::forcefield::bounds_ff::ChiralSet {
                        center: i,
                        neighbors: [neighbors[0], neighbors[1], neighbors[2], neighbors[3]],
                        lower_vol: lower,
                        upper_vol: upper,
                    });
                }
            }
        }
    }
    sets
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{Atom, Bond, BondOrder, BondStereo, ChiralType, Hybridization};
    use nalgebra::Vector3;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_calculate_bounds_matrix() {
        let mut mol = Molecule::new("Ethane_Dummy");

        let a1 = mol.add_atom(Atom {
            element: 6,
            position: Vector3::zeros(),
            charge: 0.0,
            formal_charge: 0,
            hybridization: Hybridization::SP3,
            chiral_tag: ChiralType::Unspecified,
        });

        let a2 = mol.add_atom(Atom {
            element: 6,
            position: Vector3::zeros(),
            charge: 0.0,
            formal_charge: 0,
            hybridization: Hybridization::SP3,
            chiral_tag: ChiralType::Unspecified,
        });

        mol.add_bond(
            a1,
            a2,
            Bond {
                order: BondOrder::Single,
                stereo: BondStereo::None,
            },
        );

        let bounds = calculate_bounds_matrix(&mol);
        assert_eq!(bounds.nrows(), 2);
        assert_eq!(bounds.ncols(), 2);

        // Check diagonal is zero
        assert_eq!(bounds[(0, 0)], 0.0);
        assert_eq!(bounds[(1, 1)], 0.0);

        // Min bound (lower triangle) for bond (1.52 - 0.01 = 1.51)
        assert!(bounds[(1, 0)] <= 1.51 + 0.001 && bounds[(1, 0)] >= 1.51 - 0.001);

        // Max bound (upper triangle) for bond (1.52 + 0.01 = 1.53)
        assert!(bounds[(0, 1)] >= 1.53 - 0.001 && bounds[(0, 1)] <= 1.53 + 0.001);
    }

    #[test]
    fn test_smooth_bounds_matrix() {
        // Create a 3-atom dummy matrix
        // 0-1 bonded (1.5 max, 1.4 min)
        // 1-2 bonded (1.5 max, 1.4 min)
        // 0-2 no bond (1000 max, 0 min)
        let mut b = DMatrix::from_element(3, 3, 0.0);

        b[(0, 1)] = 1.5;
        b[(1, 0)] = 1.4;
        b[(1, 2)] = 1.5;
        b[(2, 1)] = 1.4;
        b[(0, 2)] = 1000.0;
        b[(2, 0)] = 0.0;

        let smoothed = smooth_bounds_matrix(b);

        // The upper bound for 0-2 should now be at most 1.5 + 1.5 = 3.0 (Triangle inequality)
        assert!(smoothed[(0, 2)] <= 3.0);
    }

    #[test]
    fn test_metric_embedding_pipeline() {
        let mut b = DMatrix::from_element(3, 3, 0.0);
        b[(0, 1)] = 1.6;
        b[(1, 0)] = 1.4;
        b[(1, 2)] = 1.6;
        b[(2, 1)] = 1.4;
        b[(0, 2)] = 3.2;
        b[(2, 0)] = 2.8;

        let mut rng = StdRng::seed_from_u64(42);
        let dists = pick_random_distances(&mut rng, &b);

        // Assert symmetry
        assert_eq!(dists[(0, 1)], dists[(1, 0)]);
        assert!(dists[(0, 1)] >= 1.4 && dists[(0, 1)] <= 1.6);

        let metric = compute_metric_matrix(&dists);
        // Metric matrix should be 3x3
        assert_eq!(metric.nrows(), 3);

        let coords3d = generate_3d_coordinates(&metric);
        // Coordinates should be Nx3
        assert_eq!(coords3d.nrows(), 3);
        assert_eq!(coords3d.ncols(), 3);

        // Ensure not all zeros if distances are non-zero
        let sum_abs: f32 = coords3d.iter().map(|v| v.abs()).sum();
        assert!(sum_abs > 0.0);
    }

    #[test]
    fn test_calc_chiral_volume() {
        // Create a perfect tetrahedron around origin
        let mut coords = DMatrix::from_element(4, 3, 0.0);

        // p1 (1, 1, 1)
        coords[(0, 0)] = 1.0;
        coords[(0, 1)] = 1.0;
        coords[(0, 2)] = 1.0;
        // p2 (-1, -1, 1)
        coords[(1, 0)] = -1.0;
        coords[(1, 1)] = -1.0;
        coords[(1, 2)] = 1.0;
        // p3 (-1, 1, -1)
        coords[(2, 0)] = -1.0;
        coords[(2, 1)] = 1.0;
        coords[(2, 2)] = -1.0;
        // p4 (1, -1, -1)  <-- we'll use this as the center reference for the volume
        coords[(3, 0)] = 1.0;
        coords[(3, 1)] = -1.0;
        coords[(3, 2)] = -1.0;

        let vol = calc_chiral_volume(0, 1, 2, 3, &coords);

        // v1 = (0, 2, 2)
        // v2 = (-2, 0, 2)
        // v3 = (-2, 2, 0)
        // v2 x v3 = (-4, -4, -4)
        // v1 . (v2 x v3) = 0 - 8 - 8 = -16
        assert_eq!(vol, -16.0);

        // Swap 1 and 2 should invert the volume (chirality flip)
        let vol_swapped = calc_chiral_volume(1, 0, 2, 3, &coords);
        assert_eq!(vol_swapped, 16.0);
    }
}
