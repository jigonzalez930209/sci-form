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

    // Helper to flatten / unflatten coords
    let flatten = |c: &nalgebra::DMatrix<f32>| -> Vec<f32> {
        (0..n)
            .flat_map(|i| vec![c[(i, 0)], c[(i, 1)], c[(i, 2)]])
            .collect()
    };
    let unflatten = |v: &[f32], c: &mut nalgebra::DMatrix<f32>| {
        for i in 0..n {
            c[(i, 0)] = v[i * 3];
            c[(i, 1)] = v[i * 3 + 1];
            c[(i, 2)] = v[i * 3 + 2];
        }
    };
    let grad_to_vec = |g: &nalgebra::DMatrix<f32>| -> Vec<f32> {
        (0..n)
            .flat_map(|i| vec![g[(i, 0)], g[(i, 1)], g[(i, 2)]])
            .collect()
    };
    let dot = |a: &[f32], b: &[f32]| -> f32 { a.iter().zip(b).map(|(x, y)| x * y).sum() };

    let mut x = flatten(coords);
    let grad_mat = compute_analytical_gradient(coords, mol, params, bounds_matrix);
    let mut g = grad_to_vec(&grad_mat);

    let mut s_hist: Vec<Vec<f32>> = Vec::with_capacity(m);
    let mut y_hist: Vec<Vec<f32>> = Vec::with_capacity(m);
    let mut rho_hist: Vec<f32> = Vec::with_capacity(m);

    for _iter in 0..max_iters {
        // L-BFGS two-loop recursion
        let k = s_hist.len();
        let mut q = g.clone();
        let mut alpha = vec![0.0f32; k];

        for i in (0..k).rev() {
            alpha[i] = rho_hist[i] * dot(&s_hist[i], &q);
            for j in 0..dim {
                q[j] -= alpha[i] * y_hist[i][j];
            }
        }

        // Initial Hessian approximation: H0 = gamma * I
        let gamma = if k > 0 {
            dot(&s_hist[k - 1], &y_hist[k - 1]) / dot(&y_hist[k - 1], &y_hist[k - 1]).max(1e-10)
        } else {
            1.0
        };
        let mut r: Vec<f32> = q.iter().map(|&v| gamma * v).collect();

        for i in 0..k {
            let beta = rho_hist[i] * dot(&y_hist[i], &r);
            for j in 0..dim {
                r[j] += (alpha[i] - beta) * s_hist[i][j];
            }
        }

        // Direction = -H * g (r already holds H*g)
        let dir: Vec<f32> = r.iter().map(|v| -v).collect();

        // Simple backtracking line search (Armijo condition)
        let mut step = 1.0f32;
        let c1 = 1e-4;
        let dg = dot(&g, &dir);
        
        let dir_norm: f32 = dir.iter().map(|&v| v*v).sum::<f32>().sqrt();
        if step * dir_norm > 0.5 {
            step = 0.5 / dir_norm;
        }

        let f0 = calculate_total_energy(coords, mol, params, bounds_matrix);

        for _ in 0..20 {
            let x_new: Vec<f32> = x.iter().zip(&dir).map(|(xi, di)| xi + step * di).collect();
            unflatten(&x_new, coords);
            let f_new = calculate_total_energy(coords, mol, params, bounds_matrix);
            if f_new <= f0 + c1 * step * dg {
                break;
            }
            step *= 0.5;
        }

        let x_new: Vec<f32> = x.iter().zip(&dir).map(|(xi, di)| xi + step * di).collect();
        unflatten(&x_new, coords);
        let g_mat = compute_analytical_gradient(coords, mol, params, bounds_matrix);
        let g_new = grad_to_vec(&g_mat);

        let s_k: Vec<f32> = x_new.iter().zip(&x).map(|(a, b)| a - b).collect();
        let y_k: Vec<f32> = g_new.iter().zip(&g).map(|(a, b)| a - b).collect();
        let ys = dot(&y_k, &s_k);

        if ys > 1e-10 {
            if s_hist.len() == m {
                s_hist.remove(0);
                y_hist.remove(0);
                rho_hist.remove(0);
            }
            s_hist.push(s_k);
            y_hist.push(y_k);
            rho_hist.push(1.0 / ys);
        }

        // Check convergence
        let g_norm: f32 = g_new.iter().map(|v| v * v).sum::<f32>().sqrt();
        x = x_new;
        g = g_new;

        if g_norm < 1e-4 {
            break;
        }
    }

    unflatten(&x, coords);
}
