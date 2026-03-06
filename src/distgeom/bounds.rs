use crate::graph::Molecule;
use nalgebra::DMatrix;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;

pub fn get_bond_length(mol: &Molecule, n1: NodeIndex, n2: NodeIndex) -> f32 {
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

fn compute_14_dist_cis(r12: f32, r23: f32, r34: f32, a123: f32, a234: f32) -> f32 {
    let x0 = r12 * a123.cos();
    let y0 = r12 * a123.sin();
    let x3 = r23 - r34 * a234.cos();
    let y3 = r34 * a234.sin();
    ((x3 - x0).powi(2) + (y3 - y0).powi(2)).sqrt()
}

fn compute_14_dist_trans(r12: f32, r23: f32, r34: f32, a123: f32, a234: f32) -> f32 {
    let x0 = r12 * a123.cos();
    let y0 = r12 * a123.sin();
    let x3 = r23 - r34 * a234.cos();
    let y3 = -r34 * a234.sin();
    ((x3 - x0).powi(2) + (y3 - y0).powi(2)).sqrt()
}

pub fn calculate_bounds_matrix(mol: &Molecule) -> DMatrix<f32> {
    let n = mol.graph.node_count();
    let mut bounds = DMatrix::from_element(n, n, 0.0);
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
    fn set_b(b: &mut DMatrix<f32>, i: usize, j: usize, l: f32, u: f32) {
        let (mi, ma) = if i < j { (i, j) } else { (j, i) };
        if b[(ma, mi)] < l {
            b[(ma, mi)] = l;
        }
        if b[(mi, ma)] > u {
            b[(mi, ma)] = u;
        }
    }
    for bond in mol.graph.edge_references() {
        let d = get_bond_length(mol, bond.source(), bond.target());
        set_b(
            &mut bounds,
            bond.source().index(),
            bond.target().index(),
            d - 0.01,
            d + 0.01,
        );
    }
    for i in 0..n {
        let ni = NodeIndex::new(i);
        let nbs: Vec<_> = mol.graph.neighbors(ni).collect();
        for j in 0..nbs.len() {
            for k in (j + 1)..nbs.len() {
                let (n1, n2) = (nbs[j], nbs[k]);

                // Relax bounds for 1-3 in rings
                let (tol_l, tol_u) = if mol.graph.contains_edge(n1, n2) {
                    (0.1, 0.1)
                } else {
                    (0.02, 0.02)
                };

                let (l_val, u_val);
                if mol.graph.contains_edge(n1, n2) {
                    // If n1 and n2 are bonded (3-ring), use bond-based distance
                    let bond_dist = get_bond_length(mol, n1, n2);
                    l_val = bond_dist - tol_l;
                    u_val = bond_dist + tol_u;
                } else {
                    // Otherwise, use angle-based distance
                    let a = crate::graph::get_corrected_ideal_angle(mol, ni, n1, n2);
                    let (d1, d2) = (get_bond_length(mol, ni, n1), get_bond_length(mol, ni, n2));
                    let di = (d1 * d1 + d2 * d2 - 2.0 * d1 * d2 * a.cos()).sqrt();
                    l_val = di - tol_l;
                    u_val = di + tol_u;
                }
                set_b(&mut bounds, n1.index(), n2.index(), l_val, u_val);
                top_dist[(n1.index(), n2.index())] = 2;
                top_dist[(n2.index(), n1.index())] = 2;
            }
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if top_dist[(i, j)] == 3 {
                for n1 in mol.graph.neighbors(NodeIndex::new(i)) {
                    for n2 in mol.graph.neighbors(NodeIndex::new(j)) {
                        if let Some(e) = mol.graph.find_edge(n1, n2) {
                            let (r1, r2, r3) = (
                                get_bond_length(mol, NodeIndex::new(i), n1),
                                get_bond_length(mol, n1, n2),
                                get_bond_length(mol, n2, NodeIndex::new(j)),
                            );
                            let (a1, a2) = (
                                crate::graph::get_corrected_ideal_angle(
                                    mol,
                                    n1,
                                    NodeIndex::new(i),
                                    n2,
                                ),
                                crate::graph::get_corrected_ideal_angle(
                                    mol,
                                    n2,
                                    n1,
                                    NodeIndex::new(j),
                                ),
                            );
                            let (dc, dt) = (
                                compute_14_dist_cis(r1, r2, r3, a1, a2),
                                compute_14_dist_trans(r1, r2, r3, a1, a2),
                            );
                            let (l, u) = match mol.graph[e].stereo {
                                crate::graph::BondStereo::Z => (dc - 0.06, dc + 0.06),
                                crate::graph::BondStereo::E => (dt - 0.06, dt + 0.06),
                                _ => (dc.min(dt) - 0.06, dc.max(dt) + 0.06),
                            };
                            set_b(&mut bounds, i, j, l.max(0.0), u);
                        }
                    }
                }
            }
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if top_dist[(i, j)] >= 3 {
                let v = crate::graph::get_vdw_radius(mol.graph[NodeIndex::new(i)].element)
                    + crate::graph::get_vdw_radius(mol.graph[NodeIndex::new(j)].element);
                let l = match top_dist[(i, j)] {
                    3 => v * 0.5,
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

pub fn triangle_smooth(b: &mut DMatrix<f32>) -> bool {
    let n = b.nrows();
    for k in 0..n {
        for i in 0..n {
            if i == k {
                continue;
            }
            let (uk, lk) = if i < k {
                (b[(i, k)], b[(k, i)])
            } else {
                (b[(k, i)], b[(i, k)])
            };
            for j in (i + 1)..n {
                if j == k {
                    continue;
                }
                let (ukj, lkj) = if j < k {
                    (b[(j, k)], b[(k, j)])
                } else {
                    (b[(k, j)], b[(j, k)])
                };
                let (su, d1, d2) = (uk + ukj, lk - ukj, lkj - uk);
                if b[(i, j)] > su {
                    b[(i, j)] = su;
                }
                let mut li = b[(j, i)];
                if li < d1 {
                    li = d1;
                }
                if li < d2 {
                    li = d2;
                }
                b[(j, i)] = li;
                if b[(j, i)] > b[(i, j)] + 1e-4 {
                    return false;
                }
                if b[(j, i)] > b[(i, j)] {
                    b[(j, i)] = b[(i, j)];
                }
            }
        }
    }
    true
}

pub fn smooth_bounds_matrix(mut b: DMatrix<f32>) -> DMatrix<f32> {
    triangle_smooth(&mut b);
    b
}
