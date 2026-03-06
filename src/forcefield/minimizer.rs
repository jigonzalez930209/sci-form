use super::energy::*;
use super::gradients::compute_analytical_gradient;

/// Builds the topological distance matrix (Floyd-Warshall) for VDW exclusion
pub fn build_topological_distances(mol: &crate::graph::Molecule) -> nalgebra::DMatrix<usize> {
    use petgraph::visit::EdgeRef;
    let n = mol.graph.node_count();
    let mut top_dist = nalgebra::DMatrix::from_element(n, n, 1000);
    for i in 0..n {
        top_dist[(i, i)] = 0;
    }
    for edge in mol.graph.edge_references() {
        let i = edge.source().index();
        let j = edge.target().index();
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
    top_dist
}

/// Steepest Descent minimization using Numerical Gradients
pub fn minimize_energy_steepest_descent(
    coords: &mut nalgebra::DMatrix<f32>,
    mol: &crate::graph::Molecule,
    params: &FFParams,
    bounds_matrix: &nalgebra::DMatrix<f32>,
    max_iters: usize,
    step_size: f32,
    tol: f32,
) {
    let n = coords.nrows();
    let h = 1e-4;
    let mut prev_energy = calculate_total_energy(coords, mol, params, bounds_matrix);

    for _iter in 0..max_iters {
        let mut grad = nalgebra::DMatrix::from_element(n, 3, 0.0);
        for i in 0..n {
            for dim in 0..3 {
                let original = coords[(i, dim)];
                coords[(i, dim)] = original + h;
                let e_plus = calculate_total_energy(coords, mol, params, bounds_matrix);
                coords[(i, dim)] = original - h;
                let e_minus = calculate_total_energy(coords, mol, params, bounds_matrix);
                grad[(i, dim)] = ((e_plus - e_minus) / (2.0 * h)).clamp(-50.0, 50.0);
                coords[(i, dim)] = original;
            }
        }
        for i in 0..n {
            for dim in 0..3 {
                coords[(i, dim)] -= step_size * grad[(i, dim)];
            }
        }
        let current_energy = calculate_total_energy(coords, mol, params, bounds_matrix);
        if (prev_energy - current_energy).abs() < tol {
            break;
        }
        prev_energy = current_energy;
    }
}

/// L-BFGS minimization with Analytical Gradients (custom implementation)
/// Uses the two-loop recursion algorithm for the inverse Hessian approximation.
pub fn minimize_energy_lbfgs(
    coords: &mut nalgebra::DMatrix<f32>,
    mol: &crate::graph::Molecule,
    params: &FFParams,
    bounds_matrix: &nalgebra::DMatrix<f32>,
    max_iters: usize,
) {
    let n = coords.nrows();
    let dim = n * 3;
    let m = 7; // History size

    // Precompute topological terms
    use petgraph::visit::EdgeRef;
    let mut bond_terms = Vec::new();
    for edge in mol.graph.edge_references() {
        let i = edge.source().index();
        let j = edge.target().index();
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
        bond_terms.push((i, j, r_eq, params.kb));
    }

    let mut angle_terms = Vec::new();
    for i in 0..n {
        let ni = petgraph::graph::NodeIndex::new(i);
        let neighbors: Vec<_> = mol.graph.neighbors(ni).collect();
        for j in 0..neighbors.len() {
            for k in (j + 1)..neighbors.len() {
                let n1 = neighbors[j];
                let n2 = neighbors[k];
                let ideal = crate::graph::get_corrected_ideal_angle(mol, ni, n1, n2);
                angle_terms.push((i, n1.index(), n2.index(), ideal, params.k_theta));
            }
        }
    }

    let mut bounds_terms = Vec::new(); // (i, j, lower, upper)
    for i in 0..n {
        for j in (i + 1)..n {
            bounds_terms.push((i, j, bounds_matrix[(j, i)], bounds_matrix[(i, j)]));
        }
    }

    let calc_energy = |c: &[f32]| -> f32 {
        let mut e = 0.0;
        for &(i, j, r_eq, kb) in &bond_terms {
            let dx = c[i * 3] - c[j * 3];
            let dy = c[i * 3 + 1] - c[j * 3 + 1];
            let dz = c[i * 3 + 2] - c[j * 3 + 2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            e += 0.5 * kb * (r - r_eq) * (r - r_eq);
        }
        for &(c_idx, i, j, th_eq, k_th) in &angle_terms {
            let v1x = c[i * 3] - c[c_idx * 3];
            let v1y = c[i * 3 + 1] - c[c_idx * 3 + 1];
            let v1z = c[i * 3 + 2] - c[c_idx * 3 + 2];
            let v2x = c[j * 3] - c[c_idx * 3];
            let v2y = c[j * 3 + 1] - c[c_idx * 3 + 1];
            let v2z = c[j * 3 + 2] - c[c_idx * 3 + 2];
            let r1 = (v1x * v1x + v1y * v1y + v1z * v1z).sqrt();
            let r2 = (v2x * v2x + v2y * v2y + v2z * v2z).sqrt();
            if r1 > 1e-4 && r2 > 1e-4 {
                let dot = (v1x * v2x + v1y * v2y + v1z * v2z) / (r1 * r2);
                let theta = dot.clamp(-1.0, 1.0).acos();
                e += 0.5 * k_th * (theta - th_eq) * (theta - th_eq);
            }
        }
        for &(i, j, lower, upper) in &bounds_terms {
            let dx = c[i * 3] - c[j * 3];
            let dy = c[i * 3 + 1] - c[j * 3 + 1];
            let dz = c[i * 3 + 2] - c[j * 3 + 2];
            let r2 = dx * dx + dy * dy + dz * dz;
            let u2 = upper * upper;
            let l2 = lower * lower;
            if r2 > u2 {
                let val = (r2 / u2) - 1.0;
                e += params.k_bounds * val * val;
            } else if r2 < l2 {
                let val = (2.0 * l2 / (l2 + r2)) - 1.0;
                e += params.k_bounds * val * val;
            }
        }
        e
    };

    let calc_grad = |c: &[f32], grad: &mut [f32]| {
        for v in grad.iter_mut() {
            *v = 0.0;
        }
        for &(i, j, r_eq, kb) in &bond_terms {
            let dx = c[i * 3] - c[j * 3];
            let dy = c[i * 3 + 1] - c[j * 3 + 1];
            let dz = c[i * 3 + 2] - c[j * 3 + 2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r > 1e-4 {
                let pre = kb * (r - r_eq) / r;
                let gx = pre * dx;
                let gy = pre * dy;
                let gz = pre * dz;
                grad[i * 3] += gx;
                grad[i * 3 + 1] += gy;
                grad[i * 3 + 2] += gz;
                grad[j * 3] -= gx;
                grad[j * 3 + 1] -= gy;
                grad[j * 3 + 2] -= gz;
            }
        }
        for &(c_idx, i, j, th_eq, k_th) in &angle_terms {
            let v1x = c[i * 3] - c[c_idx * 3];
            let v1y = c[i * 3 + 1] - c[c_idx * 3 + 1];
            let v1z = c[i * 3 + 2] - c[c_idx * 3 + 2];
            let v2x = c[j * 3] - c[c_idx * 3];
            let v2y = c[j * 3 + 1] - c[c_idx * 3 + 1];
            let v2z = c[j * 3 + 2] - c[c_idx * 3 + 2];
            let r1 = (v1x * v1x + v1y * v1y + v1z * v1z).sqrt();
            let r2 = (v2x * v2x + v2y * v2y + v2z * v2z).sqrt();
            if r1 > 1e-4 && r2 > 1e-4 {
                let u1x = v1x / r1;
                let u1y = v1y / r1;
                let u1z = v1z / r1;
                let u2x = v2x / r2;
                let u2y = v2y / r2;
                let u2z = v2z / r2;
                let dot_val = (u1x * u2x + u1y * u2y + u1z * u2z).clamp(-1.0, 1.0);
                let sin_th = (1.0 - dot_val * dot_val).max(1e-8).sqrt();
                let pre = -k_th * (dot_val.acos() - th_eq) / sin_th;

                let d_th_dv1x = (u2x - dot_val * u1x) / r1;
                let d_th_dv1y = (u2y - dot_val * u1y) / r1;
                let d_th_dv1z = (u2z - dot_val * u1z) / r1;
                let d_th_dv2x = (u1x - dot_val * u2x) / r2;
                let d_th_dv2y = (u1y - dot_val * u2y) / r2;
                let d_th_dv2z = (u1z - dot_val * u2z) / r2;

                grad[i * 3] += pre * d_th_dv1x;
                grad[i * 3 + 1] += pre * d_th_dv1y;
                grad[i * 3 + 2] += pre * d_th_dv1z;
                grad[j * 3] += pre * d_th_dv2x;
                grad[j * 3 + 1] += pre * d_th_dv2y;
                grad[j * 3 + 2] += pre * d_th_dv2z;
                grad[c_idx * 3] -= pre * (d_th_dv1x + d_th_dv2x);
                grad[c_idx * 3 + 1] -= pre * (d_th_dv1y + d_th_dv2y);
                grad[c_idx * 3 + 2] -= pre * (d_th_dv1z + d_th_dv2z);
            }
        }
        for &(i, j, lower, upper) in &bounds_terms {
            let dx = c[i * 3] - c[j * 3];
            let dy = c[i * 3 + 1] - c[j * 3 + 1];
            let dz = c[i * 3 + 2] - c[j * 3 + 2];
            let r2 = dx * dx + dy * dy + dz * dz;
            let u2 = upper * upper;
            let l2 = lower * lower;
            if r2 > u2 {
                let pre = 4.0 * params.k_bounds * ((r2 / u2) - 1.0) / u2;
                grad[i * 3] += pre * dx;
                grad[i * 3 + 1] += pre * dy;
                grad[i * 3 + 2] += pre * dz;
                grad[j * 3] -= pre * dx;
                grad[j * 3 + 1] -= pre * dy;
                grad[j * 3 + 2] -= pre * dz;
            } else if r2 < l2 {
                let l2r2 = l2 + r2;
                let pre = 8.0 * params.k_bounds * l2 * (1.0 - 2.0 * l2 / l2r2) / (l2r2 * l2r2);
                grad[i * 3] += pre * dx;
                grad[i * 3 + 1] += pre * dy;
                grad[i * 3 + 2] += pre * dz;
                grad[j * 3] -= pre * dx;
                grad[j * 3 + 1] -= pre * dy;
                grad[j * 3 + 2] -= pre * dz;
            }
        }
    };

    let dot = |a: &[f32], b: &[f32]| -> f32 {
        let mut sum = 0.0;
        for j in 0..a.len() {
            sum += a[j] * b[j];
        }
        sum
    };

    let mut x = vec![0.0; dim];
    for i in 0..n {
        x[i * 3] = coords[(i, 0)];
        x[i * 3 + 1] = coords[(i, 1)];
        x[i * 3 + 2] = coords[(i, 2)];
    }
    let mut g = vec![0.0; dim];
    calc_grad(&x, &mut g);

    let mut s_hist: Vec<Vec<f32>> = Vec::with_capacity(m);
    let mut y_hist: Vec<Vec<f32>> = Vec::with_capacity(m);
    let mut rho_hist: Vec<f32> = Vec::with_capacity(m);

    let mut q = vec![0.0; dim];
    let mut dir = vec![0.0; dim];
    let mut x_new = vec![0.0; dim];
    let mut g_new = vec![0.0; dim];
    let mut s_k = vec![0.0; dim];
    let mut y_k = vec![0.0; dim];
    let mut alpha = vec![0.0f32; m];

    for _iter in 0..max_iters {
        let k = s_hist.len();
        q.copy_from_slice(&g);

        for i in (0..k).rev() {
            alpha[i] = rho_hist[i] * dot(&s_hist[i], &q);
            for j in 0..dim {
                q[j] -= alpha[i] * y_hist[i][j];
            }
        }

        let gamma = if k > 0 {
            dot(&s_hist[k - 1], &y_hist[k - 1]) / dot(&y_hist[k - 1], &y_hist[k - 1]).max(1e-10)
        } else {
            1.0
        };

        for j in 0..dim {
            dir[j] = q[j] * gamma;
        }
        for i in 0..k {
            let beta = rho_hist[i] * dot(&y_hist[i], &dir);
            for j in 0..dim {
                dir[j] += (alpha[i] - beta) * s_hist[i][j];
            }
        }

        for j in 0..dim {
            dir[j] = -dir[j];
        }

        let mut step = 1.0f32;
        let c1 = 1e-4;
        let dg = dot(&g, &dir);

        let mut dir_norm_sq = 0.0;
        for j in 0..dim {
            dir_norm_sq += dir[j] * dir[j];
        }
        let dir_norm = dir_norm_sq.sqrt();
        if step * dir_norm > 0.5 {
            step = 0.5 / dir_norm;
        }

        let f0 = calc_energy(&x);

        for _ in 0..20 {
            for j in 0..dim {
                x_new[j] = x[j] + step * dir[j];
            }
            let f_new = calc_energy(&x_new);
            if f_new <= f0 + c1 * step * dg {
                break;
            }
            step *= 0.5;
        }

        for j in 0..dim {
            x_new[j] = x[j] + step * dir[j];
        }
        calc_grad(&x_new, &mut g_new);

        for j in 0..dim {
            s_k[j] = x_new[j] - x[j];
            y_k[j] = g_new[j] - g[j];
        }
        let ys = dot(&y_k, &s_k);

        if ys > 1e-10 {
            if s_hist.len() == m {
                s_hist.remove(0);
                y_hist.remove(0);
                rho_hist.remove(0);
            }
            s_hist.push(s_k.clone());
            y_hist.push(y_k.clone());
            rho_hist.push(1.0 / ys);
        }

        let mut g_norm_sq = 0.0;
        for j in 0..dim {
            g_norm_sq += g_new[j] * g_new[j];
        }
        let g_norm = g_norm_sq.sqrt();

        x.copy_from_slice(&x_new);
        g.copy_from_slice(&g_new);

        if g_norm < 1e-4 {
            break;
        }
    }

    // Write back
    for i in 0..n {
        coords[(i, 0)] = x[i * 3];
        coords[(i, 1)] = x[i * 3 + 1];
        coords[(i, 2)] = x[i * 3 + 2];
    }
}
