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
        let r2_eff = r2.max(1e-4);
        let val = l2 / r2_eff - 1.0;
        let pre = -4.0 * k_bounds * val * l2 / (r2_eff * r2_eff);
        grad[(idx1, 0)] += pre * diff[0];
        grad[(idx1, 1)] += pre * diff[1];
        grad[(idx1, 2)] += pre * diff[2];
        grad[(idx2, 0)] -= pre * diff[0];
        grad[(idx2, 1)] -= pre * diff[1];
        grad[(idx2, 2)] -= pre * diff[2];
    }
}

pub fn compute_analytical_gradient(
    coords: &DMatrix<f32>,
    mol: &crate::graph::Molecule,
    params: &FFParams,
    bounds_matrix: &DMatrix<f32>,
) -> DMatrix<f32> {
    let n = mol.graph.node_count();
    let mut grad = DMatrix::zeros(n, 3);

    // 1. Bonds
    for edge in mol.graph.edge_references() {
        let idx1 = edge.source().index();
        let idx2 = edge.target().index();
        let p1 = Vector3::new(coords[(idx1, 0)], coords[(idx1, 1)], coords[(idx1, 2)]);
        let p2 = Vector3::new(coords[(idx2, 0)], coords[(idx2, 1)], coords[(idx2, 2)]);
        let r_eq = crate::distgeom::get_bond_length(mol, edge.source(), edge.target());
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
                let theta_eq = crate::graph::get_corrected_ideal_angle(mol, ni, n1, n2);
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
            let upper = bounds_matrix[(i, j)];
            let lower = bounds_matrix[(j, i)];
            analytical_grad_bounds(&p1, &p2, lower, upper, params.k_bounds, &mut grad, i, j);
        }
    }

    grad
}
