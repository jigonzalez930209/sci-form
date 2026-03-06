use nalgebra::Vector3;

/// Harmonic bond stretching penalty (Hooke's Law)
/// E = 0.5 * k_b * (r - r_eq)^2
pub fn bond_stretch_energy(p1: &Vector3<f32>, p2: &Vector3<f32>, k_b: f32, r_eq: f32) -> f32 {
    let r = (p1 - p2).norm();
    0.5 * k_b * (r - r_eq).powi(2)
}

/// Harmonic angle bending penalty
/// E = 0.5 * k_theta * (theta - theta_eq)^2
pub fn angle_bend_energy(
    p1: &Vector3<f32>,
    p2: &Vector3<f32>, // central atom
    p3: &Vector3<f32>,
    k_theta: f32,
    theta_eq: f32,
) -> f32 {
    let v1 = (p1 - p2).normalize();
    let v2 = (p3 - p2).normalize();
    let dot = v1.dot(&v2).clamp(-1.0, 1.0);
    let theta = dot.acos();
    0.5 * k_theta * (theta - theta_eq).powi(2)
}

/// Harmonic torsional potential (Dihedral)
pub fn torsional_energy(
    p1: &Vector3<f32>,
    p2: &Vector3<f32>,
    p3: &Vector3<f32>,
    p4: &Vector3<f32>,
    k_omega: f32,
    order: &crate::graph::BondOrder,
) -> f32 {
    let b1 = p2 - p1;
    let b2 = p3 - p2;
    let b3 = p4 - p3;
    if b2.norm() < 1e-4 {
        return 0.0;
    }
    let n1 = b1.cross(&b2).normalize();
    let n2 = b2.cross(&b3).normalize();
    let m = n1.cross(&b2.normalize());
    let x = n1.dot(&n2);
    let y = m.dot(&n2);
    let omega = y.atan2(x);
    match order {
        crate::graph::BondOrder::Double | crate::graph::BondOrder::Aromatic => {
            0.5 * k_omega * 5.0 * (1.0 - (2.0 * omega).cos())
        }
        _ => {
            0.5 * k_omega * (1.0 - (3.0 * omega).cos())
        }
    }
}

/// Harmonic distance bounds penalty (matches RDKit ETKDG Gram rules)
/// If r < lower: E = 0.5 * k_bounds * (lower - r)^2
/// If r > upper: E = 0.5 * k_bounds * (r - upper)^2
pub fn bounds_energy(
    p1: &Vector3<f32>,
    p2: &Vector3<f32>,
    lower: f32,
    upper: f32,
    k_bounds: f32,
) -> f32 {
    let r = (p1 - p2).norm();
    if r < lower {
        0.5 * k_bounds * (lower - r).powi(2)
    } else if r > upper {
        0.5 * k_bounds * (r - upper).powi(2)
    } else {
        0.0
    }
}

/// Out-of-Plane bending energy for sp2 atoms (planarity enforcement)
/// For a central sp2 atom C with 3 neighbors A, B, D:
/// E = k_oop * h^2, where h = distance from D to the plane(A, C, B)
pub fn oop_energy(
    center: &Vector3<f32>,
    n1: &Vector3<f32>,
    n2: &Vector3<f32>,
    n3: &Vector3<f32>,
    k_oop: f32,
) -> f32 {
    let v1 = n1 - center;
    let v2 = n2 - center;
    let normal = v1.cross(&v2);
    let norm_len = normal.norm();
    if norm_len < 1e-6 {
        return 0.0;
    }
    let n_hat = normal / norm_len;
    let v3 = n3 - center;
    let h = v3.dot(&n_hat); // signed distance from n3 to the plane
    0.5 * k_oop * h * h
}

/// Generic container for force field parameters
pub struct FFParams {
    pub kb: f32,
    pub k_theta: f32,
    pub k_omega: f32,
    pub k_oop: f32,
    pub k_bounds: f32,
}

/// Calculates the total potential energy for a molecular system
pub fn calculate_total_energy(
    coords: &nalgebra::DMatrix<f32>,
    mol: &crate::graph::Molecule,
    params: &FFParams,
    bounds_matrix: &nalgebra::DMatrix<f32>,
) -> f32 {
    let mut energy = 0.0;
    use petgraph::visit::EdgeRef;

    // 1. Bond Stretching
    for edge in mol.graph.edge_references() {
        let i = edge.source().index();
        let j = edge.target().index();
        let p1 = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
        let p2 = Vector3::new(coords[(j, 0)], coords[(j, 1)], coords[(j, 2)]);
        let atom_i = &mol.graph[petgraph::graph::NodeIndex::new(i)];
        let atom_j = &mol.graph[petgraph::graph::NodeIndex::new(j)];
        let mut r_eq = crate::graph::get_covalent_radius(atom_i.element)
            + crate::graph::get_covalent_radius(atom_j.element);
        match edge.weight().order {
            crate::graph::BondOrder::Double => r_eq -= 0.20,
            crate::graph::BondOrder::Triple => r_eq -= 0.34,
            crate::graph::BondOrder::Aromatic => r_eq -= 0.15,
            _ => {}
        }
        energy += bond_stretch_energy(&p1, &p2, params.kb, r_eq);
    }

    // 2. Angle Bending
    let n = mol.graph.node_count();
    for i in 0..n {
        let n_idx = petgraph::graph::NodeIndex::new(i);
        let neighbors: Vec<_> = mol.graph.neighbors(n_idx).collect();
        for j in 0..neighbors.len() {
            for k in (j + 1)..neighbors.len() {
                let n1 = neighbors[j];
                let n2 = neighbors[k];
                let ideal_angle = crate::graph::get_corrected_ideal_angle(mol, n_idx, n1, n2);
                let pc = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
                let p1 = Vector3::new(
                    coords[(n1.index(), 0)],
                    coords[(n1.index(), 1)],
                    coords[(n1.index(), 2)],
                );
                let p2 = Vector3::new(
                    coords[(n2.index(), 0)],
                    coords[(n2.index(), 1)],
                    coords[(n2.index(), 2)],
                );
                energy += angle_bend_energy(&p1, &pc, &p2, params.k_theta, ideal_angle);
            }
        }
    }

    // 3. Distance Bounds Enforcement (The "100% Correct Gram Distance" logic)
    for i in 0..n {
        for j in (i + 1)..n {
            // Distance bounds use upper triangle for upper bounds, lower triangle for lower bounds.
            let upper = bounds_matrix[(i, j)];
            let lower = bounds_matrix[(j, i)];
            
            // Only strictly enforce bounds for generic non-bonded pairs (to avoid competing with bond/angle terms directly)
            // But actually RDKit bounds everything. Let's do everything.
            let p1 = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
            let p2 = Vector3::new(coords[(j, 0)], coords[(j, 1)], coords[(j, 2)]);
            energy += bounds_energy(&p1, &p2, lower, upper, params.k_bounds);
        }
    }

    // 4. Torsional energy (1-4)
    for i in 0..n {
        let n1_idx = petgraph::graph::NodeIndex::new(i);
        let n1_neighbors: Vec<_> = mol.graph.neighbors(n1_idx).collect();
        for &n2_idx in &n1_neighbors {
            let n2_neighbors: Vec<_> = mol.graph.neighbors(n2_idx).collect();
            for &n3_idx in &n2_neighbors {
                if n3_idx == n1_idx {
                    continue;
                }
                let n3_neighbors: Vec<_> = mol.graph.neighbors(n3_idx).collect();
                for &n4_idx in &n3_neighbors {
                    if n4_idx == n1_idx || n4_idx == n2_idx {
                        continue;
                    }
                    if n1_idx.index() < n4_idx.index() {
                        let p1 = Vector3::new(
                            coords[(n1_idx.index(), 0)],
                            coords[(n1_idx.index(), 1)],
                            coords[(n1_idx.index(), 2)],
                        );
                        let p2 = Vector3::new(
                            coords[(n2_idx.index(), 0)],
                            coords[(n2_idx.index(), 1)],
                            coords[(n2_idx.index(), 2)],
                        );
                        let p3 = Vector3::new(
                            coords[(n3_idx.index(), 0)],
                            coords[(n3_idx.index(), 1)],
                            coords[(n3_idx.index(), 2)],
                        );
                        let p4 = Vector3::new(
                            coords[(n4_idx.index(), 0)],
                            coords[(n4_idx.index(), 1)],
                            coords[(n4_idx.index(), 2)],
                        );
                        
                        // ETKDG-lite: Infer parameters from hybridization
                        let etkdg_params = super::etkdg_lite::infer_etkdg_parameters(mol, n2_idx.index(), n3_idx.index());
                        energy += super::etkdg_lite::calc_torsion_energy_m6(&p1, &p2, &p3, &p4, &etkdg_params);

                        // Also include standard MMFF torsional energy for bond orders
                        let edge_23 = mol.graph.find_edge(n2_idx, n3_idx).unwrap();
                        let order = &mol.graph[edge_23].order;
                        energy += torsional_energy(&p1, &p2, &p3, &p4, params.k_omega, order);
                    }
                }
            }
        }
    }

    // 5. Out-of-Plane bending (sp2 atoms with exactly 3 neighbors)
    for i in 0..n {
        let ni = petgraph::graph::NodeIndex::new(i);
        let atom = &mol.graph[ni];
        if atom.hybridization != crate::graph::Hybridization::SP2 {
            continue;
        }
        let neighbors: Vec<_> = mol.graph.neighbors(ni).collect();
        if neighbors.len() != 3 {
            continue;
        }
        let pc = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
        for k in 0..3 {
            let a = neighbors[k].index();
            let b = neighbors[(k + 1) % 3].index();
            let c = neighbors[(k + 2) % 3].index();
            let pa = Vector3::new(coords[(a, 0)], coords[(a, 1)], coords[(a, 2)]);
            let pb = Vector3::new(coords[(b, 0)], coords[(b, 1)], coords[(b, 2)]);
            let pn3 = Vector3::new(coords[(c, 0)], coords[(c, 1)], coords[(c, 2)]);
            energy += oop_energy(&pc, &pa, &pb, &pn3, params.k_oop);
        }
    }

    energy
}
