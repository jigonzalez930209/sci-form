use nalgebra::{DMatrix, Vector3};

/// Computes the bounds-violation energy (RDKit's DistViolationContribs).
///
/// For each atom pair (i, j) with lower bound `lb` and upper bound `ub`:
/// - If d² > ub²: val = (d²/ub²) - 1.0; energy += weight * val²
/// - If d² < lb²: val = (2*lb²/(lb²+d²)) - 1.0; energy += weight * val²
/// - Otherwise: no penalty
/// Computes the bounds-violation energy (RDKit's DistViolationContribs).
pub fn bounds_violation_energy(coords: &DMatrix<f32>, bounds: &DMatrix<f32>) -> f32 {
    let n = coords.nrows();
    let dim = coords.ncols();
    let mut energy = 0.0f32;

    for i in 1..n {
        for j in 0..i {
            let ub = bounds[(j, i)]; // upper triangle = upper bound
            let lb = bounds[(i, j)]; // lower triangle = lower bound

            let mut d2 = 0.0f32;
            for d in 0..dim {
                let diff = coords[(i, d)] - coords[(j, d)];
                d2 += diff * diff;
            }

            let ub2 = ub * ub;
            let lb2 = lb * lb;

            let val = if d2 > ub2 {
                (d2 / ub2) - 1.0
            } else if d2 < lb2 {
                (2.0 * lb2 / (lb2 + d2)) - 1.0
            } else {
                0.0
            };

            if val > 0.0 {
                energy += val * val; // weight = 1.0
            }
        }
    }
    energy
}

/// Computes the analytical gradient of the bounds-violation energy.
///
/// Matches RDKit's DistViolationContribs::getGrad exactly.
pub fn bounds_violation_gradient(coords: &DMatrix<f32>, bounds: &DMatrix<f32>) -> DMatrix<f32> {
    let n = coords.nrows();
    let dim = coords.ncols();
    let mut grad = DMatrix::from_element(n, dim, 0.0f32);

    for i in 1..n {
        for j in 0..i {
            let ub = bounds[(j, i)];
            let lb = bounds[(i, j)];

            let mut d2 = 0.0f32;
            let mut diffs = Vec::with_capacity(dim);
            for d_idx in 0..dim {
                let diff = coords[(i, d_idx)] - coords[(j, d_idx)];
                diffs.push(diff);
                d2 += diff * diff;
            }

            let ub2 = ub * ub;
            let lb2 = lb * lb;

            if d2 > ub2 {
                let d = d2.sqrt();
                if d < 1e-8 {
                    continue;
                }
                // RDKit: 2 * (d^2/ub^2 - 1) * (2d/ub^2) = 4d * (d^2/ub^2 - 1) / ub^2
                // dE/dd = 2 * ((d^2/ub^2) - 1) * (2d/ub^2)
                // dE/dx = dE/dd * dd/dx = dE/dd * (x_i - x_j)/d
                // pre_factor = dE/dd / d = 2 * ((d^2/ub^2) - 1) * (2/ub^2)
                let pre_factor = 4.0 * ((d2 / ub2) - 1.0) / ub2;
                for d_idx in 0..dim {
                    let g = pre_factor * diffs[d_idx];
                    grad[(i, d_idx)] += g;
                    grad[(j, d_idx)] -= g;
                }
            } else if d2 < lb2 {
                let d = d2.sqrt();
                if d < 1e-8 {
                    continue;
                }
                // RDKit: 2 * (2*lb^2/(lb^2+d^2) - 1) * (2*lb^2 * (-2d) / (lb^2+d^2)^2)
                // dE/dd = 2 * (2*lb^2/(lb^2+d^2) - 1) * (-4*lb^2*d / (lb^2+d^2)^2)
                // pre_factor = dE/dd / d = 2 * (2*lb^2/(lb^2+d^2) - 1) * (-4*lb^2 / (lb^2+d^2)^2)
                let l2d2 = lb2 + d2;
                let pre_factor = 8.0 * lb2 * (1.0 - 2.0 * lb2 / l2d2) / (l2d2 * l2d2);
                for d_idx in 0..dim {
                    let g = pre_factor * diffs[d_idx];
                    grad[(i, d_idx)] += g;
                    grad[(j, d_idx)] -= g;
                }
            }
        }
    }
    grad
}

/// L-BFGS minimization using bounds-violation energy/gradient only.
/// This matches RDKit's firstMinimization step.
pub fn minimize_bounds_lbfgs(
    coords: &mut DMatrix<f32>,
    bounds: &DMatrix<f32>,
    chiral_sets: &[ChiralSet],
    max_iters: usize,
    force_tol: f32,
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
    let calc_total = |c: &DMatrix<f32>| -> (f32, DMatrix<f32>) {
        let mut e = bounds_violation_energy(c, bounds);
        let mut g = bounds_violation_gradient(c, bounds);
        if !chiral_sets.is_empty() {
            e += 1.0 * chiral_violation_energy(c, chiral_sets);
            chiral_violation_gradient(c, chiral_sets, &mut g);
        }
        // 4D penalty: higher weight for better 3D projection (RDKit uses 0.1, but we need more for stability)
        if dim_coords == 4 {
            let weight_4d = 1.0f32;
            for i in 0..n {
                let x4 = c[(i, 3)];
                e += weight_4d * x4 * x4;
                g[(i, 3)] += 2.0 * weight_4d * x4;
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
            break;
        }
    }
    unflatten(&x, coords);
}

pub struct ChiralSet {
    pub center: usize,
    pub neighbors: [usize; 4],
    pub lower_vol: f32,
    pub upper_vol: f32,
}

pub fn chiral_violation_energy(coords: &DMatrix<f32>, chiral_sets: &[ChiralSet]) -> f32 {
    let mut energy = 0.0f32;
    let dim = coords.ncols();
    for c in chiral_sets {
        let vol = crate::distgeom::calc_chiral_volume(
            c.neighbors[0],
            c.neighbors[1],
            c.neighbors[2],
            c.neighbors[3],
            coords,
        );
        if vol < c.lower_vol {
            let diff = vol - c.lower_vol;
            energy += diff * diff;
        } else if vol > c.upper_vol {
            let diff = vol - c.upper_vol;
            energy += diff * diff;
        }
    }
    energy
}

pub fn chiral_violation_gradient(
    coords: &DMatrix<f32>,
    chiral_sets: &[ChiralSet],
    grad: &mut DMatrix<f32>,
) {
    let dim = coords.ncols();
    for c in chiral_sets {
        let (idx1, idx2, idx3, idx4) = (
            c.neighbors[0],
            c.neighbors[1],
            c.neighbors[2],
            c.neighbors[3],
        );

        // v1 = pos1 - pos4, v2 = pos2 - pos4, v3 = pos3 - pos4
        let v1 = Vector3::new(
            coords[(idx1, 0)] - coords[(idx4, 0)],
            coords[(idx1, 1)] - coords[(idx4, 1)],
            coords[(idx1, 2)] - coords[(idx4, 2)],
        );
        let v2 = Vector3::new(
            coords[(idx2, 0)] - coords[(idx4, 0)],
            coords[(idx2, 1)] - coords[(idx4, 1)],
            coords[(idx2, 2)] - coords[(idx4, 2)],
        );
        let v3 = Vector3::new(
            coords[(idx3, 0)] - coords[(idx4, 0)],
            coords[(idx3, 1)] - coords[(idx4, 1)],
            coords[(idx3, 2)] - coords[(idx4, 2)],
        );

        let v2xv3 = v2.cross(&v3);
        let vol = v1.dot(&v2xv3);

        let pre_factor;
        if vol < c.lower_vol {
            pre_factor = 2.0 * (vol - c.lower_vol);
        } else if vol > c.upper_vol {
            pre_factor = 2.0 * (vol - c.upper_vol);
        } else {
            continue;
        }

        // dV/dpos1 = v2 x v3
        grad[(idx1, 0)] += pre_factor * v2xv3.x;
        grad[(idx1, 1)] += pre_factor * v2xv3.y;
        grad[(idx1, 2)] += pre_factor * v2xv3.z;

        // dV/dpos2 = v3 x v1
        let v3xv1 = v3.cross(&v1);
        grad[(idx2, 0)] += pre_factor * v3xv1.x;
        grad[(idx2, 1)] += pre_factor * v3xv1.y;
        grad[(idx2, 2)] += pre_factor * v3xv1.z;

        // dV/dpos3 = v1 x v2
        let v1xv2 = v1.cross(&v2);
        grad[(idx3, 0)] += pre_factor * v1xv2.x;
        grad[(idx3, 1)] += pre_factor * v1xv2.y;
        grad[(idx3, 2)] += pre_factor * v1xv2.z;

        // dV/dpos4 = -(dV/dpos1 + dV/dpos2 + dV/dpos3)
        grad[(idx4, 0)] -= pre_factor * (v2xv3.x + v3xv1.x + v1xv2.x);
        grad[(idx4, 1)] -= pre_factor * (v2xv3.y + v3xv1.y + v1xv2.y);
        grad[(idx4, 2)] -= pre_factor * (v2xv3.z + v3xv1.z + v1xv2.z);
    }
}
