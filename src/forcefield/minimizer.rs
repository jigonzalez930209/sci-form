use super::energy::FFParams;
use nalgebra::DMatrix;

/// Optimized L-BFGS minimizer that uses synchronized energy and gradient functions.
pub fn minimize_energy_lbfgs(
    mol: &crate::graph::Molecule,
    initial_coords: &DMatrix<f32>,
    bounds_matrix: &DMatrix<f32>,
    params: &FFParams,
    max_iter: usize,
    tol: f32,
) -> DMatrix<f32> {
    let mut coords = initial_coords.clone();

    // L-BFGS history
    let m = 10;
    let mut s_hist: Vec<DMatrix<f32>> = Vec::with_capacity(m);
    let mut y_hist: Vec<DMatrix<f32>> = Vec::with_capacity(m);
    let mut rho_hist: Vec<f32> = Vec::with_capacity(m);

    let mut g = crate::forcefield::gradients::compute_analytical_gradient(
        &coords,
        mol,
        params,
        bounds_matrix,
    );

    for iter in 0..max_iter {
        let g_norm = g.norm();
        if g_norm < tol {
            println!("L-BFGS converged in {} iterations", iter);
            break;
        }

        // Gradient Scaling for stability (ETKDG requirement)
        let mut g_scaled = g.clone();
        let max_g = g.abs().max();
        if max_g > 10.0 {
            let scale = 10.0 / max_g;
            g_scaled *= scale;
        }

        // Compute search direction p = -H * g_scaled
        let mut p = -g_scaled.clone();
        if !s_hist.is_empty() {
            let k = s_hist.len();
            let mut alphas = vec![0.0; k];
            for i in (0..k).rev() {
                let alpha = rho_hist[i] * s_hist[i].dot(&p);
                p -= alpha * &y_hist[i];
                alphas[i] = alpha;
            }

            // Initial scaling
            let s_last = &s_hist[k - 1];
            let y_last = &y_hist[k - 1];
            let gamma = s_last.dot(y_last) / y_last.dot(y_last).max(1e-10);
            p *= gamma;

            for i in 0..k {
                let beta = rho_hist[i] * y_hist[i].dot(&p);
                p += (alphas[i] - beta) * &s_hist[i];
            }
        }

        // Line search (Simple backtracking with Armijo condition)
        let mut step = 0.1;
        let c1 = 1e-4;
        let e_old =
            crate::forcefield::energy::calculate_total_energy(&coords, mol, params, bounds_matrix);
        let g_dot_p = g.dot(&p);

        let mut found_step = false;
        for _ in 0..15 {
            let next_coords = &coords + step * &p;
            let e_new = crate::forcefield::energy::calculate_total_energy(
                &next_coords,
                mol,
                params,
                bounds_matrix,
            );

            if e_new < e_old + c1 * step * g_dot_p {
                let next_g = crate::forcefield::gradients::compute_analytical_gradient(
                    &next_coords,
                    mol,
                    params,
                    bounds_matrix,
                );

                // Update history
                let s = step * &p;
                let y = &next_g - &g;
                let sy = s.dot(&y);

                if sy > 1e-10 {
                    if s_hist.len() >= m {
                        s_hist.remove(0);
                        y_hist.remove(0);
                        rho_hist.remove(0);
                    }
                    s_hist.push(s);
                    y_hist.push(y);
                    rho_hist.push(1.0 / sy);
                }

                coords = next_coords;
                g = next_g;
                found_step = true;
                break;
            }
            step *= 0.5;
        }

        if !found_step {
            // Restart if line search fails
            s_hist.clear();
            y_hist.clear();
            rho_hist.clear();
            if step < 1e-7 {
                break;
            }
        }
    }

    coords
}

pub fn calculate_rmsd_kabsch(coords: &DMatrix<f32>, reference: &DMatrix<f32>) -> f32 {
    let n = coords.nrows();
    if n == 0 {
        return 0.0;
    }

    let mut c1 = nalgebra::Vector3::zeros();
    let mut c2 = nalgebra::Vector3::zeros();
    for i in 0..n {
        c1 += nalgebra::Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
        c2 += nalgebra::Vector3::new(reference[(i, 0)], reference[(i, 1)], reference[(i, 2)]);
    }
    c1 /= n as f32;
    c2 /= n as f32;

    let mut h = nalgebra::Matrix3::zeros();
    for i in 0..n {
        let p = nalgebra::Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]) - c1;
        let q =
            nalgebra::Vector3::new(reference[(i, 0)], reference[(i, 1)], reference[(i, 2)]) - c2;
        h += p * q.transpose();
    }

    let svd = h.svd(true, true);
    let u = svd.u.unwrap();
    let v_t = svd.v_t.unwrap();
    let v = v_t.transpose();

    let mut d = nalgebra::Matrix3::identity();
    if (v * u.transpose()).determinant() < 0.0 {
        d[(2, 2)] = -1.0;
    }
    let r_mat = v * d * u.transpose();

    let mut sum_sq = 0.0;
    for i in 0..n {
        let p = nalgebra::Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]) - c1;
        let q =
            nalgebra::Vector3::new(reference[(i, 0)], reference[(i, 1)], reference[(i, 2)]) - c2;
        let rotated_p = r_mat * p;
        sum_sq += (rotated_p - q).norm_squared();
    }
    (sum_sq / n as f32).sqrt()
}
