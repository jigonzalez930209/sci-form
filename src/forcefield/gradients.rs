use super::energy::*;
use nalgebra::Vector3;

/// Analytical gradient for harmonic bond stretching
/// dE/dxi = k_b * (r - r_eq) * (xi - xj) / r
pub fn analytical_grad_bond(
    p1: &Vector3<f32>,
    p2: &Vector3<f32>,
    k_b: f32,
    r_eq: f32,
    grad: &mut nalgebra::DMatrix<f32>,
    idx1: usize,
    idx2: usize,
) {
    let diff = p1 - p2;
    let r = diff.norm();
    if r < 1e-4 {
        return;
    }
    let force_mag = k_b * (r - r_eq);
    let g = diff * (force_mag / r);
    grad[(idx1, 0)] += g[0];
    grad[(idx1, 1)] += g[1];
    grad[(idx1, 2)] += g[2];
    grad[(idx2, 0)] -= g[0];
    grad[(idx2, 1)] -= g[1];
    grad[(idx2, 2)] -= g[2];
}

/// Analytical gradient for harmonic angle bending
pub fn analytical_grad_angle(
    p1: &Vector3<f32>,
    p2: &Vector3<f32>,
    p3: &Vector3<f32>,
    k_theta: f32,
    theta_eq: f32,
    grad: &mut nalgebra::DMatrix<f32>,
    idx1: usize,
    idx2: usize,
    idx3: usize,
) {
    let v1 = p1 - p2;
    let v2 = p3 - p2;
    let r1 = v1.norm();
    let r2 = v2.norm();
    if r1 < 1e-4 || r2 < 1e-4 {
        return;
    }
    let v1_n = v1 / r1;
    let v2_n = v2 / r2;
    let cos_theta = v1_n.dot(&v2_n).clamp(-0.9999, 0.9999);
    let theta = cos_theta.acos();
    let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1e-4); // Changed sin_theta calculation and added .max(1e-4)
    if sin_theta < 1e-4 {
        return;
    }
    let pf = k_theta * (theta - theta_eq) / sin_theta;
    let g1 = pf * (cos_theta * v1_n - v2_n) / r1;
    let g3 = pf * (cos_theta * v2_n - v1_n) / r2;
    let g2 = -(g1 + g3);
    grad[(idx1, 0)] += g1[0];
    grad[(idx1, 1)] += g1[1];
    grad[(idx1, 2)] += g1[2];
    grad[(idx2, 0)] += g2[0];
    grad[(idx2, 1)] += g2[1];
    grad[(idx2, 2)] += g2[2];
    grad[(idx3, 0)] += g3[0];
    grad[(idx3, 1)] += g3[1];
    grad[(idx3, 2)] += g3[2];
}

/// Analytical gradient for Distance Bounds penalty
pub fn analytical_grad_bounds(
    p1: &Vector3<f32>,
    p2: &Vector3<f32>,
    lower: f32,
    upper: f32,
    k_bounds: f32,
    grad: &mut nalgebra::DMatrix<f32>,
    idx1: usize,
    idx2: usize,
) {
    let diff = p1 - p2;
    let r = diff.norm();
    if r < 1e-4 {
        return;
    }
    
    let force_mag = if r < lower {
        -k_bounds * (lower - r) // Push away
    } else if r > upper {
        k_bounds * (r - upper) // Pull together
    } else {
        0.0
    };
    
    if force_mag.abs() > 1e-6 {
        let g = diff * (force_mag / r);
        grad[(idx1, 0)] += g[0];
        grad[(idx1, 1)] += g[1];
        grad[(idx1, 2)] += g[2];
        grad[(idx2, 0)] -= g[0];
        grad[(idx2, 1)] -= g[1];
        grad[(idx2, 2)] -= g[2];
    }
}

/// Computes the full analytical gradient of the molecular energy
pub fn compute_analytical_gradient(
    coords: &nalgebra::DMatrix<f32>,
    mol: &crate::graph::Molecule,
    params: &FFParams,
    bounds_matrix: &nalgebra::DMatrix<f32>,
) -> nalgebra::DMatrix<f32> {
    use petgraph::visit::EdgeRef;
    let n = mol.graph.node_count();
    let mut grad = nalgebra::DMatrix::from_element(n, 3, 0.0);

    // 1. Bond Stretching
    for edge in mol.graph.edge_references() {
        let i = edge.source().index();
        let j = edge.target().index();
        let p1 = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
        let p2 = Vector3::new(coords[(j, 0)], coords[(j, 1)], coords[(j, 2)]);
        let ai = &mol.graph[petgraph::graph::NodeIndex::new(i)];
        let aj = &mol.graph[petgraph::graph::NodeIndex::new(j)];
        let mut r_eq = crate::graph::get_covalent_radius(ai.element)
            + crate::graph::get_covalent_radius(aj.element);
        match edge.weight().order {
            crate::graph::BondOrder::Double => r_eq -= 0.20,
            crate::graph::BondOrder::Triple => r_eq -= 0.34,
            crate::graph::BondOrder::Aromatic => r_eq -= 0.15,
            _ => {}
        }
        analytical_grad_bond(&p1, &p2, params.kb, r_eq, &mut grad, i, j);
    }

    // 2. Angle Bending
    for i in 0..n {
        let ni = petgraph::graph::NodeIndex::new(i);
        let neighbors: Vec<_> = mol.graph.neighbors(ni).collect();
        for j in 0..neighbors.len() {
            for k in (j + 1)..neighbors.len() {
                let n1 = neighbors[j];
                let n2 = neighbors[k];
                let ideal = crate::graph::get_corrected_ideal_angle(mol, ni, n1, n2);
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
                analytical_grad_angle(
                    &p1,
                    &pc,
                    &p2,
                    params.k_theta,
                    ideal,
                    &mut grad,
                    n1.index(),
                    i,
                    n2.index(),
                );
            }
        }
    }

    // 3. Distance Bounds Enforcement
    for i in 0..n {
        for j in (i + 1)..n {
            let upper = bounds_matrix[(i, j)];
            let lower = bounds_matrix[(j, i)];
            
            let p1 = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
            let p2 = Vector3::new(coords[(j, 0)], coords[(j, 1)], coords[(j, 2)]);
            
            analytical_grad_bounds(&p1, &p2, lower, upper, params.k_bounds, &mut grad, i, j);
        }
    }

    // 4. Torsional (numerical central-diff fallback)
    let h: f32 = 1e-4;
    for i in 0..n {
        let n1i = petgraph::graph::NodeIndex::new(i);
        for n2i in mol.graph.neighbors(n1i).collect::<Vec<_>>() {
            for n3i in mol.graph.neighbors(n2i).collect::<Vec<_>>() {
                if n3i == n1i {
                    continue;
                }
                for n4i in mol.graph.neighbors(n3i).collect::<Vec<_>>() {
                    if n4i == n1i || n4i == n2i {
                        continue;
                    }
                    if n1i.index() < n4i.index() {
                        let edge = mol.graph.find_edge(n2i, n3i).unwrap();
                        let order = &mol.graph[edge].order;
                        
                        let p1 = Vector3::new(coords[(n1i.index(), 0)], coords[(n1i.index(), 1)], coords[(n1i.index(), 2)]);
                        let p2 = Vector3::new(coords[(n2i.index(), 0)], coords[(n2i.index(), 1)], coords[(n2i.index(), 2)]);
                        let p3 = Vector3::new(coords[(n3i.index(), 0)], coords[(n3i.index(), 1)], coords[(n3i.index(), 2)]);
                        let p4 = Vector3::new(coords[(n4i.index(), 0)], coords[(n4i.index(), 1)], coords[(n4i.index(), 2)]);

                        // 1. Analytical ETKDG-lite M6 gradient
                        let etkdg_params = super::etkdg_lite::infer_etkdg_parameters(mol, n2i.index(), n3i.index());
                        super::etkdg_lite::calc_torsion_grad_m6(
                            &p1, &p2, &p3, &p4, &etkdg_params, &mut grad, 
                            n1i.index(), n2i.index(), n3i.index(), n4i.index()
                        );
                        // 2. Analytical MMFF torsion
                        let mut de_dphi_mmff = 0.0;
                        let phi = {
                            let r12 = p1 - p2;
                            let r23 = p3 - p2;
                            let r34 = p4 - p3;
                            let t1 = r12.cross(&r23);
                            let t2 = r23.cross(&r34);
                            if t1.norm() < 1e-6 || t2.norm() < 1e-6 {
                                0.0
                            } else {
                                (t1.dot(&t2) / (t1.norm() * t2.norm())).clamp(-1.0, 1.0).acos()
                            }
                        };
                        match order {
                            crate::graph::BondOrder::Double | crate::graph::BondOrder::Aromatic => {
                                de_dphi_mmff = 0.5 * params.k_omega * 5.0 * 2.0 * (2.0 * phi).sin();
                            }
                            _ => {
                                de_dphi_mmff = 0.5 * params.k_omega * 3.0 * (3.0 * phi).sin();
                            }
                        }

                        let r12 = p1 - p2;
                        let r23 = p3 - p2;
                        let r34 = p4 - p3;
                        let t1 = r12.cross(&r23);
                        let t2 = r23.cross(&r34);
                        let r23m = r23.norm();
                        if r23m > 1e-6 && t1.norm_squared() > 1e-8 && t2.norm_squared() > 1e-8 {
                            let n1 = t1 / (t1.norm_squared());
                            let n2 = t2 / (t2.norm_squared());
                            let g1 = de_dphi_mmff * r23m * n1;
                            let g4 = -de_dphi_mmff * r23m * n2;
                            let g2 = ((r12.dot(&r23) / (r23m * r23m)) - 1.0) * g1 - (r34.dot(&r23) / (r23m * r23m)) * g4;
                            let g3 = -(g1 + g2 + g4);

                            for d in 0..3 {
                                grad[(n1i.index(), d)] += g1[d];
                                grad[(n2i.index(), d)] += g2[d];
                                grad[(n3i.index(), d)] += g3[d];
                                grad[(n4i.index(), d)] += g4[d];
                            }
                        }
                    }
                }
            }
        }
    }

    // 5. OOP bending (Analytical)
    for i in 0..n {
        let ni = petgraph::graph::NodeIndex::new(i);
        let atom = &mol.graph[ni];
        if atom.hybridization != crate::graph::Hybridization::SP2 { continue; }
        let neighbors: Vec<_> = mol.graph.neighbors(ni).collect();
        if neighbors.len() != 3 { continue; }
        let (n1, n2, n3) = (neighbors[0].index(), neighbors[1].index(), neighbors[2].index());
        
        let pc = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
        let indices = [n1, n2, n3];
        for k in 0..3 {
            let ia = indices[k];
            let ib = indices[(k+1)%3];
            let id = indices[(k+2)%3];
            let pa = Vector3::new(coords[(ia, 0)], coords[(ia, 1)], coords[(ia, 2)]);
            let pb = Vector3::new(coords[(ib, 0)], coords[(ib, 1)], coords[(ib, 2)]);
            let pd = Vector3::new(coords[(id, 0)], coords[(id, 1)], coords[(id, 2)]);
            
            let v1 = pa - pc;
            let v2 = pb - pc;
            let v3 = pd - pc;
            let norm = v1.cross(&v2);
            let nl = norm.norm();
            if nl < 1e-4 { continue; }
            let n_hat = norm / nl;
            let h = v3.dot(&n_hat);
            let de_dh = params.k_oop * h;

            // Simplified analytical dH/dp_i
            let g_d = de_dh * n_hat;
            let g_c = -g_d; // Approximation: move center opposite to d
            
            for dim in 0..3 {
                grad[(id, dim)] += g_d[dim];
                grad[(i, dim)] += g_c[dim];
            }
        }
    }

    // Clamp to prevent explosions
    for i in 0..n {
        for d in 0..3 {
            if !grad[(i, d)].is_finite() {
                grad[(i, d)] = 0.0;
            }
        }
    }
    grad
}
