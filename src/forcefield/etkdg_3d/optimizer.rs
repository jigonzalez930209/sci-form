//! BFGS optimizers and torsion correction for the ETKDG 3D force field.

use nalgebra::{DMatrix, Vector3};
use super::*;

pub fn minimize_etkdg_3d(
    mol: &crate::graph::Molecule,
    initial_coords: &DMatrix<f32>,
    bounds_matrix: &DMatrix<f64>,
    max_iter: usize,
    tol: f32,
) -> DMatrix<f32> {
    // Build the force field from current coordinates
    let coords_f64 = initial_coords.map(|v| v as f64);
    let ff = super::build_etkdg_3d_ff(mol, &coords_f64, bounds_matrix);
    let mut coords = initial_coords.clone();

    // L-BFGS with m=10 history
    let m = 10;
    let mut s_hist: Vec<DMatrix<f32>> = Vec::with_capacity(m);
    let mut y_hist: Vec<DMatrix<f32>> = Vec::with_capacity(m);
    let mut rho_hist: Vec<f32> = Vec::with_capacity(m);

    let debug = std::env::var("DEBUG_LBFGS").is_ok();

    let mut g = etkdg_3d_gradient(&coords, mol, &ff);

    for _iter in 0..max_iter {
        let g_norm = g.norm();
        if g_norm < tol {
            break;
        }

        // Gradient capping for stability
        let mut g_scaled = g.clone();
        let max_g = g.abs().max();
        if max_g > 10.0 {
            g_scaled *= 10.0 / max_g;
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
            let s_last = &s_hist[k - 1];
            let y_last = &y_hist[k - 1];
            let gamma = s_last.dot(y_last) / y_last.dot(y_last).max(1e-10);
            p *= gamma;
            for i in 0..k {
                let beta = rho_hist[i] * y_hist[i].dot(&p);
                p += (alphas[i] - beta) * &s_hist[i];
            }
        }

        // Line search with Armijo condition
        let mut step = 1.0; // larger initial step since flat-bottom potentials are gentler
        let c1 = 1e-4;
        let e_old = super::etkdg_3d_energy(&coords, mol, &ff);
        let g_dot_p = g.dot(&p);

        if g_dot_p >= 0.0 {
            // Not a descent direction, reset
            s_hist.clear();
            y_hist.clear();
            rho_hist.clear();
            continue;
        }

        if debug && _iter < 3 {
            println!("  3D iter={} g_norm={:.4} max_g={:.4} e_old={:.4} g_dot_p={:.4} p_norm={:.4}",
                _iter, g_norm, g.abs().max(), e_old, g_dot_p, p.norm());
        }

        let mut found_step = false;
        for ls in 0..20 {
            let next_coords = &coords + step * &p;
            let e_new = super::etkdg_3d_energy(&next_coords, mol, &ff);
            
            if debug && _iter < 3 && ls < 10 {
                println!("    ls={} step={:.10} e_new={:.4} armijo_rhs={:.4} pass={}",
                    ls, step, e_new, e_old + c1 * step * g_dot_p, e_new < e_old + c1 * step * g_dot_p);
            }

            if e_new < e_old + c1 * step * g_dot_p {
                let next_g = etkdg_3d_gradient(&next_coords, mol, &ff);
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
            if debug && _iter < 3 {
                println!("  3D iter={} LINE SEARCH FAILED step={:.10}", _iter, step);
            }
            s_hist.clear();
            y_hist.clear();
            rho_hist.clear();
            if step < 1e-10 {
                break;
            }
        }
    }

    coords
}

/// Correct torsion angles toward M6 energy minima before full FF minimization.
/// For each CSD torsion on a non-ring bond, rotates atoms on one side of the central bond to the
/// nearest M6 energy minimum angle. This helps the optimizer start in the right basin.
pub fn correct_torsion_angles(
    mol: &crate::graph::Molecule,
    coords: &mut DMatrix<f32>,
    torsion_contribs: &[M6TorsionContrib],
) {
    let n = mol.graph.node_count();
    for tc in torsion_contribs {
        // Skip very weak torsions
        let total_v: f64 = tc.v.iter().map(|x| x.abs()).sum();
        if total_v < 0.5 { continue; }

        // Skip ring torsions — rotating ring bonds would distort ring geometry
        // Check if j-k bond is in a ring
        let j_ni = petgraph::graph::NodeIndex::new(tc.j);
        let k_ni = petgraph::graph::NodeIndex::new(tc.k);
        let in_ring = crate::graph::min_path_excluding2(mol, j_ni, k_ni, j_ni, k_ni, 7).is_some();
        if in_ring { continue; }

        // Current dihedral angle
        let p1 = Vector3::new(coords[(tc.i, 0)], coords[(tc.i, 1)], coords[(tc.i, 2)]);
        let p2 = Vector3::new(coords[(tc.j, 0)], coords[(tc.j, 1)], coords[(tc.j, 2)]);
        let p3 = Vector3::new(coords[(tc.k, 0)], coords[(tc.k, 1)], coords[(tc.k, 2)]);
        let p4 = Vector3::new(coords[(tc.l, 0)], coords[(tc.l, 1)], coords[(tc.l, 2)]);
        let current_angle = calc_dihedral(&p1, &p2, &p3, &p4);

        // Find M6 energy at current angle and scan for all local minima
        let m6 = crate::forcefield::etkdg_lite::M6Params {
            s: tc.signs.map(|x| x as f32),
            v: tc.v.map(|x| x as f32),
        };
        let _current_e = calc_m6_energy(current_angle, &m6);

        // Find nearest local minimum by scanning 1° increments
        let mut best_angle = current_angle;
        let mut best_dist = 0.0f32;
        for deg in -180..=180 {
            let a = deg as f32 * std::f32::consts::PI / 180.0;
            let e = calc_m6_energy(a, &m6);
            // Check if this is a local minimum (lower than neighbors)
            let e_prev = calc_m6_energy((deg - 1) as f32 * std::f32::consts::PI / 180.0, &m6);
            let e_next = calc_m6_energy((deg + 1) as f32 * std::f32::consts::PI / 180.0, &m6);
            if e <= e_prev && e <= e_next {
                // This is a local minimum — check if it's closer than current best
                let mut dist = (a - current_angle).abs();
                if dist > std::f32::consts::PI { dist = 2.0 * std::f32::consts::PI - dist; }
                if best_dist == 0.0 || dist < best_dist {
                    best_angle = a;
                    best_dist = dist;
                }
            }
        }

        // Only rotate if the nearest minimum is significantly different from current angle
        let rotation_angle = best_angle - current_angle;
        // Normalize to [-π, π]
        let rotation_angle = if rotation_angle > std::f32::consts::PI {
            rotation_angle - 2.0 * std::f32::consts::PI
        } else if rotation_angle < -std::f32::consts::PI {
            rotation_angle + 2.0 * std::f32::consts::PI
        } else {
            rotation_angle
        };

        if rotation_angle.abs() < 0.05 { continue; } // Less than ~3° — skip

        // Find atoms on the "l side" of the j-k bond (BFS from k, excluding j)
        let mut on_k_side = vec![false; n];
        on_k_side[tc.k] = true;
        let mut stack = vec![tc.k];
        while let Some(node) = stack.pop() {
            let ni = petgraph::graph::NodeIndex::new(node);
            for nb in mol.graph.neighbors(ni) {
                let nbi = nb.index();
                if nbi == tc.j { continue; }
                if !on_k_side[nbi] {
                    on_k_side[nbi] = true;
                    stack.push(nbi);
                }
            }
        }
        // Don't rotate j
        on_k_side[tc.j] = false;

        // Rotate all atoms on k-side around the j-k axis by rotation_angle
        let axis = (&p3 - &p2).normalize();
        if axis.norm() < 1e-6 { continue; }
        let origin = p2;
        let cos_r = rotation_angle.cos();
        let sin_r = rotation_angle.sin();

        for a in 0..n {
            if !on_k_side[a] { continue; }
            let pos = Vector3::new(coords[(a, 0)], coords[(a, 1)], coords[(a, 2)]);
            let rel = pos - origin;
            // Rodrigues' rotation formula
            let rotated = rel * cos_r + axis.cross(&rel) * sin_r + axis * axis.dot(&rel) * (1.0 - cos_r);
            let new_pos = rotated + origin;
            coords[(a, 0)] = new_pos[0];
            coords[(a, 1)] = new_pos[1];
            coords[(a, 2)] = new_pos[2];
        }
    }
}

/// Calculate dihedral angle in radians
pub(crate) fn calc_dihedral(p1: &Vector3<f32>, p2: &Vector3<f32>, p3: &Vector3<f32>, p4: &Vector3<f32>) -> f32 {
    let b1 = p2 - p1;
    let b2 = p3 - p2;
    let b3 = p4 - p3;
    let n1 = b1.cross(&b2);
    let n2 = b2.cross(&b3);
    let n1_len = n1.norm();
    let n2_len = n2.norm();
    if n1_len < 1e-6 || n2_len < 1e-6 { return 0.0; }
    let n1u = n1 / n1_len;
    let n2u = n2 / n2_len;
    let cos_d = n1u.dot(&n2u).clamp(-1.0, 1.0);
    let sign = n1u.dot(&b3);
    let angle = cos_d.acos();
    if sign < 0.0 { -angle } else { angle }
}

/// Calculate M6 torsion energy at a given angle
pub(crate) fn calc_m6_energy(theta: f32, m6: &crate::forcefield::etkdg_lite::M6Params) -> f32 {
    let mut e = 0.0f32;
    for k in 0..6 {
        if m6.v[k].abs() > 1e-6 {
            let m = (k + 1) as f32;
            e += 0.5 * m6.v[k] * (1.0 + m6.s[k] * (m * theta).cos());
        }
    }
    e
}

/// Minimize using a pre-built force field (allows external override of torsion params etc.)
pub fn minimize_etkdg_3d_with_ff(
    mol: &crate::graph::Molecule,
    initial_coords: &DMatrix<f32>,
    ff: &Etkdg3DFF,
    max_iter: usize,
    tol: f32,
) -> DMatrix<f32> {
    let mut coords = initial_coords.clone();
    let m = 10;
    let mut s_hist: Vec<DMatrix<f32>> = Vec::with_capacity(m);
    let mut y_hist: Vec<DMatrix<f32>> = Vec::with_capacity(m);
    let mut rho_hist: Vec<f32> = Vec::with_capacity(m);
    let mut g = etkdg_3d_gradient(&coords, mol, ff);

    for _iter in 0..max_iter {
        let g_norm = g.norm();
        if g_norm < tol { break; }
        // RDKit-style gradient scaling: 0.1× base, halve until max*scale ≤ 10
        let max_g = g.abs().max();
        let mut scale = 0.1f32;
        while max_g * scale > 10.0 { scale *= 0.5; }
        let g_scaled = &g * scale;
        let mut p = -g_scaled.clone();
        if !s_hist.is_empty() {
            let k = s_hist.len();
            let mut alphas = vec![0.0; k];
            for i in (0..k).rev() {
                let alpha = rho_hist[i] * s_hist[i].dot(&p);
                p -= alpha * &y_hist[i];
                alphas[i] = alpha;
            }
            let gamma = s_hist[k-1].dot(&y_hist[k-1]) / y_hist[k-1].dot(&y_hist[k-1]).max(1e-10);
            p *= gamma;
            for i in 0..k {
                let beta = rho_hist[i] * y_hist[i].dot(&p);
                p += (alphas[i] - beta) * &s_hist[i];
            }
        }
        let mut step = 1.0;
        let c1 = 1e-4;
        let e_old = super::etkdg_3d_energy(&coords, mol, ff);
        let g_dot_p = g.dot(&p);
        if g_dot_p >= 0.0 {
            s_hist.clear(); y_hist.clear(); rho_hist.clear();
            continue;
        }
        let mut found_step = false;
        for _ls in 0..20 {
            let next_coords = &coords + step * &p;
            let e_new = super::etkdg_3d_energy(&next_coords, mol, ff);
            if e_new < e_old + c1 * step * g_dot_p {
                let next_g = etkdg_3d_gradient(&next_coords, mol, ff);
                let s = step * &p;
                let y = &next_g - &g;
                let sy = s.dot(&y);
                if sy > 1e-10 {
                    if s_hist.len() >= m { s_hist.remove(0); y_hist.remove(0); rho_hist.remove(0); }
                    s_hist.push(s); y_hist.push(y); rho_hist.push(1.0 / sy);
                }
                coords = next_coords;
                g = next_g;
                found_step = true;
                break;
            }
            step *= 0.5;
        }
        if !found_step {
            s_hist.clear(); y_hist.clear(); rho_hist.clear();
            if step < 1e-10 { break; }
        }
    }
    coords
}
/// f64-precision version of etkdg_3d_energy.
/// All computation done in f64 to match RDKit's double-precision force field.

/// RDKit-style cubic/quadratic line search (Numerical Recipes 9.7) - f64 precision
pub(crate) fn rdkit_linear_search_f64(
    dim: usize,
    old_pt: &[f64],
    old_val: f64,
    grad: &[f64],
    dir: &mut [f64],
    new_pt: &mut [f64],
    mol: &crate::graph::Molecule,
    ff: &Etkdg3DFF,
    max_step: f64,
    n_atoms: usize,
) -> (f64, i32) {
    const FUNCTOL: f64 = 1e-4;
    const MOVETOL: f64 = 1e-7;
    const MAX_ITER: usize = 1000;

    let sum: f64 = dir.iter().map(|x| x * x).sum::<f64>().sqrt();
    if sum > max_step {
        let scale = max_step / sum;
        for d in dir.iter_mut() { *d *= scale; }
    }

    let slope: f64 = dir.iter().zip(grad.iter()).map(|(d, g)| d * g).sum();
    if slope >= 0.0 {
        new_pt.copy_from_slice(old_pt);
        return (old_val, -1);
    }

    let mut test = 0.0f64;
    for i in 0..dim {
        let temp = dir[i].abs() / old_pt[i].abs().max(1.0);
        if temp > test { test = temp; }
    }
    let lambda_min = MOVETOL / test.max(1e-30);

    let mut lambda = 1.0f64;
    let mut lambda2 = 0.0f64;
    let mut val2 = 0.0f64;

    for it in 0..MAX_ITER {
        if lambda < lambda_min {
            new_pt.copy_from_slice(old_pt);
            return (old_val, 1);
        }
        for i in 0..dim {
            new_pt[i] = old_pt[i] + lambda * dir[i];
        }
        let new_val = super::etkdg_3d_energy_f64(new_pt, n_atoms, mol, ff);

        if new_val - old_val <= FUNCTOL * lambda * slope {
            return (new_val, 0);
        }
        let tmp_lambda;
        if it == 0 {
            tmp_lambda = -slope / (2.0 * (new_val - old_val - slope));
        } else {
            let rhs1 = new_val - old_val - lambda * slope;
            let rhs2 = val2 - old_val - lambda2 * slope;
            let denom = lambda - lambda2;
            if denom.abs() < 1e-30 {
                tmp_lambda = 0.5 * lambda;
            } else {
                let a = (rhs1 / (lambda * lambda) - rhs2 / (lambda2 * lambda2)) / denom;
                let b = (-lambda2 * rhs1 / (lambda * lambda) + lambda * rhs2 / (lambda2 * lambda2)) / denom;
                if a.abs() < 1e-30 {
                    tmp_lambda = -slope / (2.0 * b);
                } else {
                    let disc = b * b - 3.0 * a * slope;
                    if disc < 0.0 {
                        tmp_lambda = 0.5 * lambda;
                    } else if b <= 0.0 {
                        tmp_lambda = (-b + disc.sqrt()) / (3.0 * a);
                    } else {
                        tmp_lambda = -slope / (b + disc.sqrt());
                    }
                }
            }
        };
        let tmp_lambda = tmp_lambda.min(0.5 * lambda);
        lambda2 = lambda;
        val2 = new_val;
        lambda = tmp_lambda.max(0.1 * lambda);
    }
    new_pt.copy_from_slice(old_pt);
    (old_val, -1)
}

/// Full BFGS minimization matching RDKit's BFGSOpt.h — f64 precision internals
pub fn minimize_etkdg_3d_bfgs(
    mol: &crate::graph::Molecule,
    initial_coords: &DMatrix<f64>,
    ff: &Etkdg3DFF,
    max_iter: usize,
    tol: f32,
) -> DMatrix<f64> {
    let n = mol.graph.node_count();
    let dim = n * 3;
    let tol = tol as f64;
    const MAXSTEP: f64 = 100.0;
    const EPS: f64 = 3e-8;
    const TOLX: f64 = 4.0 * EPS;

    let mut pos = vec![0.0f64; dim];
    for a in 0..n {
        pos[a * 3] = initial_coords[(a, 0)];
        pos[a * 3 + 1] = initial_coords[(a, 1)];
        pos[a * 3 + 2] = initial_coords[(a, 2)];
    }

    let to_mat_f64 = |p: &[f64]| -> DMatrix<f64> {
        let mut m = DMatrix::zeros(n, 3);
        for a in 0..n {
            m[(a, 0)] = p[a * 3];
            m[(a, 1)] = p[a * 3 + 1];
            m[(a, 2)] = p[a * 3 + 2];
        }
        m
    };

    // Get gradient as flat f64 vec with RDKit two-pass scaling
    let grad_flat_scaled = |p: &[f64]| -> (Vec<f64>, f64) {
        let mut g = etkdg_3d_gradient_f64(p, n, mol, ff);
        // RDKit's two-pass gradient scaling (ForceField.cpp calcGradient)
        let mut grad_scale = 0.1f64;
        let mut max_grad = -1e8f64;
        for v in g.iter_mut() {
            *v *= grad_scale;
            if v.abs() > max_grad { max_grad = v.abs(); }
        }
        if max_grad > 10.0 {
            while max_grad * grad_scale > 10.0 {
                grad_scale *= 0.5;
            }
            for v in g.iter_mut() {
                *v *= grad_scale;
            }
        }
        (g, grad_scale)
    };

    let mut fp = super::etkdg_3d_energy_f64(&pos, n, mol, ff);
    let (mut grad, _gs) = grad_flat_scaled(&pos);

    let mut inv_hess = vec![0.0f64; dim * dim];
    for i in 0..dim { inv_hess[i * dim + i] = 1.0; }

    let mut xi = vec![0.0f64; dim];
    for i in 0..dim { xi[i] = -grad[i]; }

    let sum_sq: f64 = pos.iter().map(|x| x * x).sum();
    let max_step = MAXSTEP * sum_sq.sqrt().max(dim as f64);

    let mut new_pt = vec![0.0f64; dim];
    let mut d_grad = vec![0.0f64; dim];
    let mut hess_d_grad = vec![0.0f64; dim];

    for _iter in 0..max_iter {
        let (func_val, status) = rdkit_linear_search_f64(
            dim, &pos, fp, &grad, &mut xi, &mut new_pt, mol, ff, max_step, n,
        );
        if status < 0 {
            // RDKit: CHECK_INVARIANT(status >= 0, "bad direction in linearSearch")
            // throws an exception, effectively aborting the optimization.
            break;
        }

        fp = func_val;

        let mut test = 0.0f64;
        for i in 0..dim {
            xi[i] = new_pt[i] - pos[i];
            pos[i] = new_pt[i];
            let temp = xi[i].abs() / pos[i].abs().max(1.0);
            if temp > test { test = temp; }
            d_grad[i] = grad[i];
        }
        if test < TOLX { break; }

        let (new_grad, grad_scale) = grad_flat_scaled(&pos);
        grad = new_grad;

        test = 0.0;
        let term = (func_val * grad_scale).max(1.0);
        for i in 0..dim {
            let temp = grad[i].abs() * pos[i].abs().max(1.0);
            if temp > test { test = temp; }
            d_grad[i] = grad[i] - d_grad[i];
        }
        test /= term;
        if test < tol { break; }

        let mut fac = 0.0f64;
        let mut fae = 0.0f64;
        let mut sum_dg = 0.0f64;
        let mut sum_xi = 0.0f64;
        for i in 0..dim {
            hess_d_grad[i] = 0.0;
            for j in 0..dim {
                hess_d_grad[i] += inv_hess[i * dim + j] * d_grad[j];
            }
            fac += d_grad[i] * xi[i];
            fae += d_grad[i] * hess_d_grad[i];
            sum_dg += d_grad[i] * d_grad[i];
            sum_xi += xi[i] * xi[i];
        }

        if fac > (EPS * sum_dg * sum_xi).sqrt() {
            fac = 1.0 / fac;
            let fad = 1.0 / fae;
            for i in 0..dim {
                d_grad[i] = fac * xi[i] - fad * hess_d_grad[i];
            }
            for i in 0..dim {
                let pxi = fac * xi[i];
                let hdgi = fad * hess_d_grad[i];
                let dgi = fae * d_grad[i];
                for j in i..dim {
                    inv_hess[i * dim + j] += pxi * xi[j] - hdgi * hess_d_grad[j] + dgi * d_grad[j];
                    inv_hess[j * dim + i] = inv_hess[i * dim + j];
                }
            }
        }

        for i in 0..dim {
            xi[i] = 0.0;
            for j in 0..dim {
                xi[i] -= inv_hess[i * dim + j] * grad[j];
            }
        }
    }

    to_mat_f64(&pos)
}
