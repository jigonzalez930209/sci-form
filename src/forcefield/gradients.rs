use super::energy::FFParams;
use nalgebra::{DMatrix, Vector3};
use petgraph::visit::EdgeRef;

/// Analytical gradient for bond stretching
pub fn analytical_grad_bond(
    p1: &Vector3<f32>,
    p2: &Vector3<f32>,
    kb: f32,
    r_eq: f32,
    grad: &mut DMatrix<f32>,
    idx1: usize,
    idx2: usize,
) {
    let diff = p1 - p2;
    let r = diff.norm();
    if r < 1e-4 {
        return;
    }
    let pf = kb * (r - r_eq) / r;
    let g = diff * pf;
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
    pc: &Vector3<f32>,
    p2: &Vector3<f32>,
    k_theta: f32,
    theta_eq: f32,
    g1: &mut Vector3<f32>,
    gc: &mut Vector3<f32>,
    g2: &mut Vector3<f32>,
) {
    let v1 = p1 - pc;
    let v2 = p2 - pc;
    let r1 = v1.norm();
    let r2 = v2.norm();
    if r1 < 1e-4 || r2 < 1e-4 {
        return;
    }
    let u1 = v1 / r1;
    let u2 = v2 / r2;
    let cos_th = u1.dot(&u2).clamp(-0.999999, 0.999999);

    if (theta_eq - std::f32::consts::PI).abs() < 1e-4 {
        // Linear case: E = k * (1 + cos_th)
        // dE/dp1 = k * d(cos_th)/dp1 = k * (u2 - cos_th * u1) / r1
        *g1 = (u2 - cos_th * u1) * (k_theta / r1);
        *g2 = (u1 - cos_th * u2) * (k_theta / r2);
        *gc = -(*g1) - (*g2);
        return;
    }

    let theta = cos_th.acos();
    let sin_th = theta.sin().max(1e-4);
    let pf = k_theta * (theta - theta_eq) / sin_th;

    *g1 = (cos_th * u1 - u2) * (pf / r1);
    *g2 = (cos_th * u2 - u1) * (pf / r2);
    *gc = -(*g1) - (*g2);
}

/// Analytical gradient for distance bounds
pub fn analytical_grad_bounds(
    p1: &Vector3<f32>,
    p2: &Vector3<f32>,
    lower: f32,
    upper: f32,
    k_bounds: f32,
    grad: &mut DMatrix<f32>,
    idx1: usize,
    idx2: usize,
) {
    let diff = p1 - p2;
    let r2 = diff.norm_squared();
    let u2 = upper * upper;
    let l2 = lower * lower;

    if r2 > u2 && u2 > 1e-6 {
        let val = r2 / u2 - 1.0;
        let pre = 4.0 * k_bounds * val / u2;
        grad[(idx1, 0)] += pre * diff[0];
        grad[(idx1, 1)] += pre * diff[1];
        grad[(idx1, 2)] += pre * diff[2];
        grad[(idx2, 0)] -= pre * diff[0];
        grad[(idx2, 1)] -= pre * diff[1];
        grad[(idx2, 2)] -= pre * diff[2];
    } else if r2 < l2 && l2 > 1e-6 {
        // RDKit formula: val = 2L²/(L²+d²) - 1
        let l2d2 = l2 + r2.max(1e-6);
        let pre = 8.0 * k_bounds * l2 * (1.0 - 2.0 * l2 / l2d2) / (l2d2 * l2d2);
        grad[(idx1, 0)] += pre * diff[0];
        grad[(idx1, 1)] += pre * diff[1];
        grad[(idx1, 2)] += pre * diff[2];
        grad[(idx2, 0)] -= pre * diff[0];
        grad[(idx2, 1)] -= pre * diff[1];
        grad[(idx2, 2)] -= pre * diff[2];
    }
}

/// Analytical gradient for torsion energy V*(1 + cos(n*phi - gamma))
/// Uses the standard dihedral gradient formulation
pub fn analytical_grad_torsion(
    p1: &Vector3<f32>,
    p2: &Vector3<f32>,
    p3: &Vector3<f32>,
    p4: &Vector3<f32>,
    v: f32,
    n_fold: f32,
    gamma: f32,
    grad: &mut DMatrix<f32>,
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
) {
    let b1 = p2 - p1;
    let b2 = p3 - p2;
    let b3 = p4 - p3;

    let n1 = b1.cross(&b2);
    let n2 = b2.cross(&b3);

    let n1_sq = n1.norm_squared();
    let n2_sq = n2.norm_squared();
    if n1_sq < 1e-10 || n2_sq < 1e-10 {
        return;
    }

    let b2_len = b2.norm();
    if b2_len < 1e-6 {
        return;
    }

    let m1 = n1.cross(&b2) / b2_len;
    let x = n1.dot(&n2) / (n1_sq.sqrt() * n2_sq.sqrt());
    let y = m1.dot(&n2) / (m1.norm() * n2_sq.sqrt());
    let phi = y.atan2(x);

    // dE/dphi = -V * n * sin(n*phi - gamma)
    let de_dphi = -v * n_fold * (n_fold * phi - gamma).sin();

    // dphi/dr for each atom using standard torsion gradient formulas
    // See Blondel & Karplus, J. Comp. Chem. 17, 1132-1141 (1996)
    // Note: sign convention matches our atan2-based phi definition
    let g1 = b2_len / n1_sq * n1;
    let g4 = -b2_len / n2_sq * n2;
    let b1_dot_b2 = b1.dot(&b2) / (b2_len * b2_len);
    let b3_dot_b2 = b3.dot(&b2) / (b2_len * b2_len);
    let g2 = (-b1_dot_b2 - 1.0) * g1 + b3_dot_b2 * g4;
    let g3 = (-b3_dot_b2 - 1.0) * g4 + b1_dot_b2 * g1;

    for d in 0..3 {
        grad[(idx1, d)] += de_dphi * g1[d];
        grad[(idx2, d)] += de_dphi * g2[d];
        grad[(idx3, d)] += de_dphi * g3[d];
        grad[(idx4, d)] += de_dphi * g4[d];
    }
}

/// Analytical gradient for out-of-plane bending using triple product formulation
/// E = k * V² where V = (p1-pc)·((p2-pc)×(p3-pc))
pub fn analytical_grad_oop(
    pc: &Vector3<f32>,
    p1: &Vector3<f32>,
    p2: &Vector3<f32>,
    p3: &Vector3<f32>,
    k_oop: f32,
    grad: &mut DMatrix<f32>,
    idx_c: usize,
    idx1: usize,
    idx2: usize,
    idx3: usize,
) {
    let v1 = p1 - pc;
    let v2 = p2 - pc;
    let v3 = p3 - pc;
    let vol = v1.dot(&v2.cross(&v3));
    let pre = 2.0 * k_oop * vol;

    let g1 = v2.cross(&v3) * pre;
    let g2 = v3.cross(&v1) * pre;
    let g3 = v1.cross(&v2) * pre;
    let gc = -(g1 + g2 + g3);

    for d in 0..3 {
        grad[(idx_c, d)] += gc[d];
        grad[(idx1, d)] += g1[d];
        grad[(idx2, d)] += g2[d];
        grad[(idx3, d)] += g3[d];
    }
}

/// Analytical gradient for flat-bottom distance constraint
/// Zero gradient when minLen <= d <= maxLen
pub fn analytical_grad_distance_constraint(
    p1: &Vector3<f32>,
    p2: &Vector3<f32>,
    min_len: f32,
    max_len: f32,
    k: f32,
    grad: &mut DMatrix<f32>,
    idx1: usize,
    idx2: usize,
) {
    let diff = p1 - p2;
    let d2 = diff.norm_squared();
    let pre;
    if d2 < min_len * min_len {
        let d = d2.sqrt().max(1e-8);
        // dE/dr = k * (d - minLen), where d < minLen so this is negative
        // dE/dx_i = dE/dr * (x_i - x_j) / d = k * (1 - minLen/d) * (x_i - x_j)
        pre = k * (d - min_len) / d;
    } else if d2 > max_len * max_len {
        let d = d2.sqrt().max(1e-8);
        // dE/dr = k * (d - maxLen), where d > maxLen so this is positive
        pre = k * (d - max_len) / d;
    } else {
        return; // within bounds, zero gradient
    }
    grad[(idx1, 0)] += pre * diff[0];
    grad[(idx1, 1)] += pre * diff[1];
    grad[(idx1, 2)] += pre * diff[2];
    grad[(idx2, 0)] -= pre * diff[0];
    grad[(idx2, 1)] -= pre * diff[1];
    grad[(idx2, 2)] -= pre * diff[2];
}

pub fn compute_analytical_gradient(
    coords: &DMatrix<f32>,
    mol: &crate::graph::Molecule,
    params: &FFParams,
    bounds_matrix: &DMatrix<f64>,
) -> DMatrix<f32> {
    let n = mol.graph.node_count();
    let mut grad = DMatrix::zeros(n, 3);

    // 1. Bonds
    for edge in mol.graph.edge_references() {
        let idx1 = edge.source().index();
        let idx2 = edge.target().index();
        let p1 = Vector3::new(coords[(idx1, 0)], coords[(idx1, 1)], coords[(idx1, 2)]);
        let p2 = Vector3::new(coords[(idx2, 0)], coords[(idx2, 1)], coords[(idx2, 2)]);
        let r_eq = crate::distgeom::get_bond_length(mol, edge.source(), edge.target()) as f32;
        analytical_grad_bond(&p1, &p2, params.kb, r_eq, &mut grad, idx1, idx2);
    }

    // 2. Angles
    for i in 0..n {
        let ni = petgraph::graph::NodeIndex::new(i);
        let nbs: Vec<_> = mol.graph.neighbors(ni).collect();
        for j in 0..nbs.len() {
            for k in (j + 1)..nbs.len() {
                let n1 = nbs[j];
                let n2 = nbs[k];
                let p1 = Vector3::new(
                    coords[(n1.index(), 0)],
                    coords[(n1.index(), 1)],
                    coords[(n1.index(), 2)],
                );
                let pc = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
                let p2 = Vector3::new(
                    coords[(n2.index(), 0)],
                    coords[(n2.index(), 1)],
                    coords[(n2.index(), 2)],
                );
                let theta_eq = crate::graph::get_corrected_ideal_angle(mol, ni, n1, n2) as f32;
                let mut g1 = Vector3::zeros();
                let mut gc = Vector3::zeros();
                let mut g2 = Vector3::zeros();
                analytical_grad_angle(
                    &p1,
                    &pc,
                    &p2,
                    params.k_theta,
                    theta_eq,
                    &mut g1,
                    &mut gc,
                    &mut g2,
                );

                grad[(n1.index(), 0)] += g1.x;
                grad[(n1.index(), 1)] += g1.y;
                grad[(n1.index(), 2)] += g1.z;
                grad[(i, 0)] += gc.x;
                grad[(i, 1)] += gc.y;
                grad[(i, 2)] += gc.z;
                grad[(n2.index(), 0)] += g2.x;
                grad[(n2.index(), 1)] += g2.y;
                grad[(n2.index(), 2)] += g2.z;
            }
        }
    }

    // 3. Distance Bounds
    for i in 0..n {
        for j in (i + 1)..n {
            let p1 = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
            let p2 = Vector3::new(coords[(j, 0)], coords[(j, 1)], coords[(j, 2)]);
            let upper = bounds_matrix[(i, j)] as f32;
            let lower = bounds_matrix[(j, i)] as f32;
            analytical_grad_bounds(&p1, &p2, lower, upper, params.k_bounds, &mut grad, i, j);
        }
    }

    // 3b. Out-of-plane for SP2 with 3 neighbors
    if params.k_oop.abs() > 1e-8 {
        for i in 0..n {
            let ni = petgraph::graph::NodeIndex::new(i);
            if mol.graph[ni].hybridization != crate::graph::Hybridization::SP2 {
                continue;
            }
            let nbs: Vec<_> = mol.graph.neighbors(ni).collect();
            if nbs.len() != 3 {
                continue;
            }
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
            analytical_grad_oop(&pc, &p1, &p2, &p3, params.k_oop,
                &mut grad, i, nbs[0].index(), nbs[1].index(), nbs[2].index());
        }
    }

    // 4. Torsions (UFF-style hybridization-dependent)
    if params.k_omega.abs() > 1e-8 && n >= 4 {
        for edge in mol.graph.edge_references() {
            let u = edge.source();
            let v = edge.target();
            let hyb_u = mol.graph[u].hybridization;
            let hyb_v = mol.graph[v].hybridization;
            if hyb_u == crate::graph::Hybridization::SP
                || hyb_v == crate::graph::Hybridization::SP
            {
                continue;
            }

            let (n_fold, gamma, weight) = crate::forcefield::energy::torsion_params(hyb_u, hyb_v);

            let neighbors_u: Vec<_> = mol.graph.neighbors(u).filter(|&x| x != v).collect();
            let neighbors_v: Vec<_> = mol.graph.neighbors(v).filter(|&x| x != u).collect();

            for &nu in &neighbors_u {
                for &nv in &neighbors_v {
                    let p1 = Vector3::new(
                        coords[(nu.index(), 0)],
                        coords[(nu.index(), 1)],
                        coords[(nu.index(), 2)],
                    );
                    let p2 = Vector3::new(
                        coords[(u.index(), 0)],
                        coords[(u.index(), 1)],
                        coords[(u.index(), 2)],
                    );
                    let p3 = Vector3::new(
                        coords[(v.index(), 0)],
                        coords[(v.index(), 1)],
                        coords[(v.index(), 2)],
                    );
                    let p4 = Vector3::new(
                        coords[(nv.index(), 0)],
                        coords[(nv.index(), 1)],
                        coords[(nv.index(), 2)],
                    );
                    analytical_grad_torsion(
                        &p1, &p2, &p3, &p4,
                        params.k_omega * weight, n_fold, gamma,
                        &mut grad, nu.index(), u.index(), v.index(), nv.index(),
                    );
                }
            }
        }
    }

    // 5. ETKDG-lite M6 torsion gradient (non-ring rotatable bonds)
    if n >= 4 {
        for edge in mol.graph.edge_references() {
            let u = edge.source();
            let v = edge.target();
            if crate::graph::min_path_excluding2(mol, u, v, u, v, 7).is_some() {
                continue;
            }
            let m6 = crate::forcefield::etkdg_lite::infer_etkdg_parameters(mol, u.index(), v.index());
            if m6.v.iter().all(|&x| x.abs() < 1e-6) { continue; }

            let neighbors_u: Vec<_> = mol.graph.neighbors(u).filter(|&x| x != v).collect();
            let neighbors_v: Vec<_> = mol.graph.neighbors(v).filter(|&x| x != u).collect();
            if neighbors_u.is_empty() || neighbors_v.is_empty() { continue; }
            let nu = neighbors_u[0];
            let nv = neighbors_v[0];
            let p1 = Vector3::new(coords[(nu.index(), 0)], coords[(nu.index(), 1)], coords[(nu.index(), 2)]);
            let p2 = Vector3::new(coords[(u.index(), 0)], coords[(u.index(), 1)], coords[(u.index(), 2)]);
            let p3 = Vector3::new(coords[(v.index(), 0)], coords[(v.index(), 1)], coords[(v.index(), 2)]);
            let p4 = Vector3::new(coords[(nv.index(), 0)], coords[(nv.index(), 1)], coords[(nv.index(), 2)]);
            crate::forcefield::etkdg_lite::calc_torsion_grad_m6(
                &p1, &p2, &p3, &p4, &m6,
                &mut grad, nu.index(), u.index(), v.index(), nv.index(),
            );
        }
    }

    // 6. VDW non-bonded gradient (1-4+ pairs)
    if params.k_vdw.abs() > 1e-8 {
        let mut excluded = std::collections::HashSet::new();
        for edge in mol.graph.edge_references() {
            let a = edge.source().index();
            let b = edge.target().index();
            let (lo, hi) = if a < b { (a, b) } else { (b, a) };
            excluded.insert((lo, hi));
        }
        for center in 0..n {
            let nc = petgraph::graph::NodeIndex::new(center);
            let nbs: Vec<_> = mol.graph.neighbors(nc).collect();
            for j in 0..nbs.len() {
                for k in (j + 1)..nbs.len() {
                    let a = nbs[j].index();
                    let b = nbs[k].index();
                    let (lo, hi) = if a < b { (a, b) } else { (b, a) };
                    excluded.insert((lo, hi));
                }
            }
        }
        let mut is_14 = std::collections::HashSet::new();
        for edge in mol.graph.edge_references() {
            let u = edge.source();
            let v = edge.target();
            let neighbors_u: Vec<_> = mol.graph.neighbors(u).filter(|&x| x != v).collect();
            let neighbors_v: Vec<_> = mol.graph.neighbors(v).filter(|&x| x != u).collect();
            for &nu in &neighbors_u {
                for &nv in &neighbors_v {
                    let a = nu.index();
                    let b = nv.index();
                    if a == b { continue; }
                    let (lo, hi) = if a < b { (a, b) } else { (b, a) };
                    if !excluded.contains(&(lo, hi)) {
                        is_14.insert((lo, hi));
                    }
                }
            }
        }

        for i in 0..n {
            let ei = mol.graph[petgraph::graph::NodeIndex::new(i)].element;
            let (xi, di) = crate::forcefield::energy::uff_vdw_params(ei);
            let pi = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
            for j in (i + 1)..n {
                if excluded.contains(&(i, j)) { continue; }
                let ej = mol.graph[petgraph::graph::NodeIndex::new(j)].element;
                let (xj, dj) = crate::forcefield::energy::uff_vdw_params(ej);
                let r_star = (xi + xj) * 0.5;
                let eps_full = (di * dj).sqrt();
                let scale = if is_14.contains(&(i, j)) { 0.5 } else { 1.0 };
                let pj = Vector3::new(coords[(j, 0)], coords[(j, 1)], coords[(j, 2)]);

                let diff = pi - pj;
                let r2 = diff.norm_squared().max(1e-6);
                let r = r2.sqrt();
                if r < 0.5 || r > 8.0 { continue; }

                let u = r_star / r;
                let u6 = u * u * u * u * u * u;
                let u12 = u6 * u6;
                // dE/dr = epsilon * 12 * (u6 - u12) / r
                let de_dr = eps_full * 12.0 * (u6 - u12) / r;
                let pre = params.k_vdw * scale * de_dr / r;
                let gx = pre * diff.x;
                let gy = pre * diff.y;
                let gz = pre * diff.z;
                grad[(i, 0)] += gx;
                grad[(i, 1)] += gy;
                grad[(i, 2)] += gz;
                grad[(j, 0)] -= gx;
                grad[(j, 1)] -= gy;
                grad[(j, 2)] -= gz;
            }
        }
    }

    grad
}
