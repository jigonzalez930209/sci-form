//! LBFGS and BFGS optimizers for bounds distance force field.

use super::*;
use nalgebra::DMatrix;

pub fn minimize_bounds_lbfgs(
    coords: &mut DMatrix<f32>,
    bounds: &DMatrix<f64>,
    chiral_sets: &[ChiralSet],
    max_iters: usize,
    force_tol: f32,
) {
    minimize_bounds_lbfgs_ex(
        coords,
        bounds,
        chiral_sets,
        max_iters,
        force_tol,
        1000.0,
        0.1,
    );
}

/// Returns true if converged (gradient norm < force_tol), false if max_iters reached
pub fn minimize_bounds_lbfgs_ex(
    coords: &mut DMatrix<f32>,
    bounds: &DMatrix<f64>,
    chiral_sets: &[ChiralSet],
    max_iters: usize,
    force_tol: f32,
    basin_thresh: f32,
    weight_4d: f32,
) -> bool {
    minimize_bounds_lbfgs_full(
        coords,
        bounds,
        chiral_sets,
        max_iters,
        force_tol,
        basin_thresh,
        weight_4d,
        1.0,
    )
}

/// Full bounds FF minimizer with configurable chiral weight.
/// RDKit uses chiral_weight=1.0 for firstMinimization, chiral_weight=0.2 for minimizeFourthDimension.
pub fn minimize_bounds_lbfgs_full(
    coords: &mut DMatrix<f32>,
    bounds: &DMatrix<f64>,
    chiral_sets: &[ChiralSet],
    max_iters: usize,
    force_tol: f32,
    basin_thresh: f32,
    weight_4d: f32,
    weight_chiral: f32,
) -> bool {
    let n = coords.nrows();
    let dim_coords = coords.ncols();
    let dim_tot = n * dim_coords;
    let m = 7;

    let flatten = |c: &DMatrix<f32>| -> Vec<f32> {
        (0..n)
            .flat_map(|i| (0..dim_coords).map(move |d| c[(i, d)]))
            .collect()
    };
    let unflatten = |v: &[f32], c: &mut DMatrix<f32>| {
        for i in 0..n {
            for d in 0..dim_coords {
                c[(i, d)] = v[i * dim_coords + d];
            }
        }
    };
    let dot = |a: &[f32], b: &[f32]| -> f32 { a.iter().zip(b).map(|(x, y)| x * y).sum() };

    let mut x = flatten(coords);
    let calc_total = |c: &DMatrix<f32>| -> (f32, DMatrix<f32>) {
        let mut e = super::bounds_violation_energy_basin(c, bounds, basin_thresh);
        let mut g = super::bounds_violation_gradient_basin(c, bounds, basin_thresh);
        if !chiral_sets.is_empty() {
            e += weight_chiral * super::chiral_violation_energy(c, chiral_sets);
            let mut cg = DMatrix::from_element(n, dim_coords, 0.0f32);
            super::chiral_violation_gradient(c, chiral_sets, &mut cg);
            g += weight_chiral * cg;
        }
        // 4D penalty: penalize 4th dimension to enable good 3D projection
        if dim_coords == 4 {
            for i in 0..n {
                let x4 = c[(i, 3)];
                e += weight_4d * x4 * x4;
                g[(i, 3)] += weight_4d * x4;
            }
        }
        (e, g)
    };

    let (mut f, g_mat) = calc_total(coords);
    let mut g = vec![0.0; dim_tot];
    for i in 0..n {
        for d in 0..dim_coords {
            g[i * dim_coords + d] = g_mat[(i, d)];
        }
    }

    let mut s_hist: Vec<Vec<f32>> = Vec::with_capacity(m);
    let mut y_hist: Vec<Vec<f32>> = Vec::with_capacity(m);
    let mut rho_hist: Vec<f32> = Vec::with_capacity(m);

    let mut q = vec![0.0; dim_tot];
    let mut dir = vec![0.0; dim_tot];
    let mut x_new = vec![0.0; dim_tot];
    let mut g_new = vec![0.0; dim_tot];
    let mut s_k = vec![0.0; dim_tot];
    let mut y_k = vec![0.0; dim_tot];
    let mut alpha = vec![0.0f32; m];

    let mut converged = false;
    for _iter in 0..max_iters {
        let k = s_hist.len();
        q.copy_from_slice(&g);

        for i in (0..k).rev() {
            alpha[i] = rho_hist[i] * dot(&s_hist[i], &q);
            for j in 0..dim_tot {
                q[j] -= alpha[i] * y_hist[i][j];
            }
        }

        let gamma = if k > 0 {
            dot(&s_hist[k - 1], &y_hist[k - 1]) / dot(&y_hist[k - 1], &y_hist[k - 1]).max(1e-10)
        } else {
            1.0
        };

        for j in 0..dim_tot {
            dir[j] = q[j] * gamma;
        }
        for i in 0..k {
            let beta = rho_hist[i] * dot(&y_hist[i], &dir);
            for j in 0..dim_tot {
                dir[j] += (alpha[i] - beta) * s_hist[i][j];
            }
        }

        for j in 0..dim_tot {
            dir[j] = -dir[j];
        }
        let mut step = 1.0f32;
        let c1 = 1e-4;
        let dg = dot(&g, &dir);
        let mut dir_norm_sq = 0.0;
        for j in 0..dim_tot {
            dir_norm_sq += dir[j] * dir[j];
        }
        let dir_norm = dir_norm_sq.sqrt();
        if step * dir_norm > 0.5 {
            step = 0.5 / dir_norm;
        }

        let f0 = f;
        for _ in 0..20 {
            for j in 0..dim_tot {
                x_new[j] = x[j] + step * dir[j];
            }
            unflatten(&x_new, coords);
            let (f_new, _) = calc_total(coords);
            if f_new <= f0 + c1 * step * dg {
                break;
            }
            step *= 0.5;
        }

        for j in 0..dim_tot {
            x_new[j] = x[j] + step * dir[j];
        }
        unflatten(&x_new, coords);
        let (f_new, g_new_mat) = calc_total(coords);
        for i in 0..n {
            for d in 0..dim_coords {
                g_new[i * dim_coords + d] = g_new_mat[(i, d)];
            }
        }

        for j in 0..dim_tot {
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
        for j in 0..dim_tot {
            g_norm_sq += g_new[j] * g_new[j];
        }
        let g_norm = g_norm_sq.sqrt();

        f = f_new;
        x.copy_from_slice(&x_new);
        g.copy_from_slice(&g_new);

        if g_norm < force_tol {
            converged = true;
            break;
        }
    }
    unflatten(&x, coords);
    converged
}

/// Like minimize_bounds_lbfgs_ex but also applies M6 torsion preferences during DG optimization.
/// Torsion energy/gradient are computed using the first 3 coordinates (even in 4D mode).
pub fn minimize_bounds_with_torsions(
    coords: &mut DMatrix<f32>,
    bounds: &DMatrix<f64>,
    chiral_sets: &[ChiralSet],
    max_iters: usize,
    force_tol: f32,
    basin_thresh: f32,
    weight_4d: f32,
    weight_chiral: f32,
    torsion_contribs: &[crate::forcefield::etkdg_3d::M6TorsionContrib],
) -> bool {
    let n = coords.nrows();
    let dim_coords = coords.ncols();
    let dim_tot = n * dim_coords;
    let m = 7;

    let flatten = |c: &DMatrix<f32>| -> Vec<f32> {
        (0..n)
            .flat_map(|i| (0..dim_coords).map(move |d| c[(i, d)]))
            .collect()
    };
    let unflatten = |v: &[f32], c: &mut DMatrix<f32>| {
        for i in 0..n {
            for d in 0..dim_coords {
                c[(i, d)] = v[i * dim_coords + d];
            }
        }
    };
    let dot = |a: &[f32], b: &[f32]| -> f32 { a.iter().zip(b).map(|(x, y)| x * y).sum() };

    let mut x = flatten(coords);
    let calc_total = |c: &DMatrix<f32>| -> (f32, DMatrix<f32>) {
        let mut e = super::bounds_violation_energy_basin(c, bounds, basin_thresh);
        let mut g = super::bounds_violation_gradient_basin(c, bounds, basin_thresh);
        if !chiral_sets.is_empty() {
            e += weight_chiral * super::chiral_violation_energy(c, chiral_sets);
            let mut cg = DMatrix::from_element(n, dim_coords, 0.0f32);
            super::chiral_violation_gradient(c, chiral_sets, &mut cg);
            g += weight_chiral * cg;
        }
        // 4D penalty
        if dim_coords == 4 {
            for i in 0..n {
                let x4 = c[(i, 3)];
                e += weight_4d * x4 * x4;
                g[(i, 3)] += weight_4d * x4;
            }
        }
        // Torsion preferences (computed in 3D, gradient applied to first 3 dims)
        for tc in torsion_contribs {
            let p1 = nalgebra::Vector3::new(c[(tc.i, 0)], c[(tc.i, 1)], c[(tc.i, 2)]);
            let p2 = nalgebra::Vector3::new(c[(tc.j, 0)], c[(tc.j, 1)], c[(tc.j, 2)]);
            let p3 = nalgebra::Vector3::new(c[(tc.k, 0)], c[(tc.k, 1)], c[(tc.k, 2)]);
            let p4 = nalgebra::Vector3::new(c[(tc.l, 0)], c[(tc.l, 1)], c[(tc.l, 2)]);
            let m6 = crate::forcefield::etkdg_lite::M6Params {
                s: tc.signs.map(|x| x as f32),
                v: tc.v.map(|x| x as f32),
            };
            e += crate::forcefield::etkdg_lite::calc_torsion_energy_m6(&p1, &p2, &p3, &p4, &m6);
            crate::forcefield::etkdg_lite::calc_torsion_grad_m6(
                &p1, &p2, &p3, &p4, &m6, &mut g, tc.i, tc.j, tc.k, tc.l,
            );
        }
        (e, g)
    };

    let (mut f, g_mat) = calc_total(coords);
    let mut g = vec![0.0; dim_tot];
    for i in 0..n {
        for d in 0..dim_coords {
            g[i * dim_coords + d] = g_mat[(i, d)];
        }
    }

    let mut s_hist: Vec<Vec<f32>> = Vec::with_capacity(m);
    let mut y_hist: Vec<Vec<f32>> = Vec::with_capacity(m);
    let mut rho_hist: Vec<f32> = Vec::with_capacity(m);

    let mut q = vec![0.0; dim_tot];
    let mut dir = vec![0.0; dim_tot];
    let mut x_new = vec![0.0; dim_tot];
    let mut g_new = vec![0.0; dim_tot];
    let mut s_k = vec![0.0; dim_tot];
    let mut y_k = vec![0.0; dim_tot];
    let mut alpha = vec![0.0f32; m];

    let mut converged = false;
    for _iter in 0..max_iters {
        let k = s_hist.len();
        q.copy_from_slice(&g);

        for i in (0..k).rev() {
            alpha[i] = rho_hist[i] * dot(&s_hist[i], &q);
            for j in 0..dim_tot {
                q[j] -= alpha[i] * y_hist[i][j];
            }
        }

        let gamma = if k > 0 {
            dot(&s_hist[k - 1], &y_hist[k - 1]) / dot(&y_hist[k - 1], &y_hist[k - 1]).max(1e-10)
        } else {
            1.0
        };

        for j in 0..dim_tot {
            dir[j] = q[j] * gamma;
        }
        for i in 0..k {
            let beta = rho_hist[i] * dot(&y_hist[i], &dir);
            for j in 0..dim_tot {
                dir[j] += (alpha[i] - beta) * s_hist[i][j];
            }
        }

        for j in 0..dim_tot {
            dir[j] = -dir[j];
        }
        let mut step = 1.0f32;
        let c1 = 1e-4;
        let dg = dot(&g, &dir);
        let mut dir_norm_sq = 0.0;
        for j in 0..dim_tot {
            dir_norm_sq += dir[j] * dir[j];
        }
        let dir_norm = dir_norm_sq.sqrt();
        if step * dir_norm > 0.5 {
            step = 0.5 / dir_norm;
        }

        let f0 = f;
        for _ in 0..20 {
            for j in 0..dim_tot {
                x_new[j] = x[j] + step * dir[j];
            }
            unflatten(&x_new, coords);
            let (f_new, _) = calc_total(coords);
            if f_new <= f0 + c1 * step * dg {
                break;
            }
            step *= 0.5;
        }

        unflatten(&x_new, coords);
        let (f_new, g_new_mat) = calc_total(coords);
        for i in 0..n {
            for d in 0..dim_coords {
                g_new[i * dim_coords + d] = g_new_mat[(i, d)];
            }
        }

        for j in 0..dim_tot {
            s_k[j] = x_new[j] - x[j];
            y_k[j] = g_new[j] - g[j];
        }
        let ys = dot(&y_k, &s_k);
        if ys > 1e-10 {
            if s_hist.len() >= m {
                s_hist.remove(0);
                y_hist.remove(0);
                rho_hist.remove(0);
            }
            s_hist.push(s_k.clone());
            y_hist.push(y_k.clone());
            rho_hist.push(1.0 / ys);
        }

        let mut g_norm_sq = 0.0;
        for j in 0..dim_tot {
            g_norm_sq += g_new[j] * g_new[j];
        }
        let g_norm = g_norm_sq.sqrt();

        f = f_new;
        x.copy_from_slice(&x_new);
        g.copy_from_slice(&g_new);

        if g_norm < force_tol {
            converged = true;
            break;
        }
    }
    unflatten(&x, coords);
    converged
}

pub fn minimize_embedding_lbfgs(
    coords: &mut DMatrix<f32>,
    bounds: &DMatrix<f64>,
    chiral_sets: &[ChiralSet],
    torsion_terms: &[EmbedTorsion],
    torsion_weight: f32,
    max_iters: usize,
    force_tol: f32,
    basin_thresh: f32,
    weight_4d: f32,
) {
    let n = coords.nrows();
    let dim_coords = coords.ncols();
    let dim_tot = n * dim_coords;
    let m = 7;

    let flatten = |c: &DMatrix<f32>| -> Vec<f32> {
        (0..n)
            .flat_map(|i| (0..dim_coords).map(move |d| c[(i, d)]))
            .collect()
    };
    let unflatten = |v: &[f32], c: &mut DMatrix<f32>| {
        for i in 0..n {
            for d in 0..dim_coords {
                c[(i, d)] = v[i * dim_coords + d];
            }
        }
    };
    let dot = |a: &[f32], b: &[f32]| -> f32 { a.iter().zip(b).map(|(x, y)| x * y).sum() };

    let mut x = flatten(coords);
    let calc_energy = |c: &DMatrix<f32>| -> f32 {
        let mut e = super::bounds_violation_energy_basin(c, bounds, basin_thresh);
        if !chiral_sets.is_empty() {
            e += super::chiral_violation_energy(c, chiral_sets);
        }
        if !torsion_terms.is_empty() {
            e += torsion_weight * super::energy::torsion_energy_4d(c, torsion_terms);
        }
        if dim_coords == 4 {
            for i in 0..n {
                e += weight_4d * c[(i, 3)] * c[(i, 3)];
            }
        }
        e
    };
    let calc_energy_and_grad = |c: &DMatrix<f32>| -> (f32, DMatrix<f32>) {
        let mut e = super::bounds_violation_energy_basin(c, bounds, basin_thresh);
        let mut g = super::bounds_violation_gradient_basin(c, bounds, basin_thresh);
        if !chiral_sets.is_empty() {
            e += super::chiral_violation_energy(c, chiral_sets);
            super::chiral_violation_gradient(c, chiral_sets, &mut g);
        }
        if !torsion_terms.is_empty() {
            e += torsion_weight * super::energy::torsion_energy_4d(c, torsion_terms);
            let n2 = c.nrows();
            let dim2 = c.ncols();
            let mut tg = DMatrix::zeros(n2, dim2);
            super::energy::torsion_gradient_4d(c, torsion_terms, &mut tg);
            for i in 0..n2 {
                for d in 0..dim2 {
                    g[(i, d)] += torsion_weight * tg[(i, d)];
                }
            }
        }
        if dim_coords == 4 {
            for i in 0..n {
                let x4 = c[(i, 3)];
                e += weight_4d * x4 * x4;
                g[(i, 3)] += weight_4d * x4;
            }
        }
        (e, g)
    };

    let (mut f, g_mat) = calc_energy_and_grad(coords);
    let mut g = vec![0.0; dim_tot];
    for i in 0..n {
        for d in 0..dim_coords {
            g[i * dim_coords + d] = g_mat[(i, d)];
        }
    }

    let mut s_hist: Vec<Vec<f32>> = Vec::with_capacity(m);
    let mut y_hist: Vec<Vec<f32>> = Vec::with_capacity(m);
    let mut rho_hist: Vec<f32> = Vec::with_capacity(m);
    let mut q = vec![0.0; dim_tot];
    let mut dir = vec![0.0; dim_tot];
    let mut x_new = vec![0.0; dim_tot];
    let mut g_new = vec![0.0; dim_tot];
    let mut s_k = vec![0.0; dim_tot];
    let mut y_k = vec![0.0; dim_tot];
    let mut alpha = vec![0.0f32; m];

    for _iter in 0..max_iters {
        let k = s_hist.len();
        q.copy_from_slice(&g);
        for i in (0..k).rev() {
            alpha[i] = rho_hist[i] * dot(&s_hist[i], &q);
            for j in 0..dim_tot {
                q[j] -= alpha[i] * y_hist[i][j];
            }
        }
        let gamma = if k > 0 {
            dot(&s_hist[k - 1], &y_hist[k - 1]) / dot(&y_hist[k - 1], &y_hist[k - 1]).max(1e-10)
        } else {
            1.0
        };
        for j in 0..dim_tot {
            dir[j] = q[j] * gamma;
        }
        for i in 0..k {
            let beta = rho_hist[i] * dot(&y_hist[i], &dir);
            for j in 0..dim_tot {
                dir[j] += (alpha[i] - beta) * s_hist[i][j];
            }
        }
        for j in 0..dim_tot {
            dir[j] = -dir[j];
        }

        let mut step = 1.0f32;
        let c1 = 1e-4;
        let dg = dot(&g, &dir);
        let dir_norm: f32 = dir.iter().map(|x| x * x).sum::<f32>().sqrt();
        if step * dir_norm > 0.5 {
            step = 0.5 / dir_norm;
        }

        let f0 = f;
        for _ in 0..20 {
            for j in 0..dim_tot {
                x_new[j] = x[j] + step * dir[j];
            }
            unflatten(&x_new, coords);
            let f_new = calc_energy(coords);
            if f_new <= f0 + c1 * step * dg {
                break;
            }
            step *= 0.5;
        }

        for j in 0..dim_tot {
            x_new[j] = x[j] + step * dir[j];
        }
        unflatten(&x_new, coords);
        let (f_new, g_new_mat) = calc_energy_and_grad(coords);
        for i in 0..n {
            for d in 0..dim_coords {
                g_new[i * dim_coords + d] = g_new_mat[(i, d)];
            }
        }
        for j in 0..dim_tot {
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

        f = f_new;
        x.copy_from_slice(&x_new);
        g.copy_from_slice(&g_new);

        let g_norm: f32 = g_new.iter().map(|x| x * x).sum::<f32>().sqrt();
        if g_norm < force_tol {
            break;
        }
    }
    unflatten(&x, coords);
}

// ── RDKit-matching full BFGS optimizer ──────────────────────────────────────
