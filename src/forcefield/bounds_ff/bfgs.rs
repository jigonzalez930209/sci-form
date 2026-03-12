use nalgebra::DMatrix;
use super::*;

const BFGS_FUNCTOL: f64 = 1e-4;
const BFGS_MOVETOL: f64 = 1e-7;
const BFGS_EPS: f64 = 3e-8;
const BFGS_TOLX: f64 = 4.0 * BFGS_EPS;
const BFGS_MAXSTEP: f64 = 100.0;
const BFGS_MAX_ITER_LINEAR_SEARCH: usize = 1000;

pub(crate) fn rdkit_linear_search(
    pos: &[f64], val: f64, grad: &[f64],
    dir: &mut [f64], // may be rescaled in place
    func: &dyn Fn(&[f64]) -> f64,
    max_step: f64,
) -> (Vec<f64>, f64, i32) {
    let dim = pos.len();
    let mut res_code: i32 = -1;

    // get the length of the direction vector
    let mut sum = 0.0f64;
    for i in 0..dim {
        sum += dir[i] * dir[i];
    }
    sum = sum.sqrt();

    // rescale if moving too far
    if sum > max_step {
        for i in 0..dim {
            dir[i] *= max_step / sum;
        }
    }

    // check direction has component along -grad
    let mut slope = 0.0f64;
    for i in 0..dim {
        slope += dir[i] * grad[i];
    }
    if slope >= 0.0 {
        return (pos.to_vec(), val, -1);
    }

    let mut test = 0.0f64;
    for i in 0..dim {
        let temp = dir[i].abs() / pos[i].abs().max(1.0);
        if temp > test {
            test = temp;
        }
    }
    let lambda_min = BFGS_MOVETOL / test;
    let mut lambda = 1.0f64;
    let mut lambda2 = 0.0f64;
    let mut val2 = 0.0f64;
    let mut new_pt = vec![0.0f64; dim];
    let mut new_val = val;

    for _it in 0..BFGS_MAX_ITER_LINEAR_SEARCH {
        if lambda < lambda_min {
            res_code = 1;
            break;
        }
        for i in 0..dim {
            new_pt[i] = pos[i] + lambda * dir[i];
        }
        new_val = func(&new_pt);

        if new_val - val <= BFGS_FUNCTOL * lambda * slope {
            res_code = 0;
            return (new_pt, new_val, res_code);
        }

        // backtrack
        let tmp_lambda = if _it == 0 {
            -slope / (2.0 * (new_val - val - slope))
        } else {
            let rhs1 = new_val - val - lambda * slope;
            let rhs2 = val2 - val - lambda2 * slope;
            let a = (rhs1 / (lambda * lambda) - rhs2 / (lambda2 * lambda2))
                / (lambda - lambda2);
            let b = (-lambda2 * rhs1 / (lambda * lambda)
                + lambda * rhs2 / (lambda2 * lambda2))
                / (lambda - lambda2);
            if a == 0.0 {
                -slope / (2.0 * b)
            } else {
                let disc = b * b - 3.0 * a * slope;
                if disc < 0.0 {
                    0.5 * lambda
                } else if b <= 0.0 {
                    (-b + disc.sqrt()) / (3.0 * a)
                } else {
                    -slope / (b + disc.sqrt())
                }
            }
        };
        let tmp_lambda = if tmp_lambda > 0.5 * lambda { 0.5 * lambda } else { tmp_lambda };
        lambda2 = lambda;
        val2 = new_val;
        lambda = tmp_lambda.max(0.1 * lambda);
    }

    // nothing was done — return original position
    if res_code != 0 {
        return (pos.to_vec(), new_val, res_code);
    }
    (new_pt, new_val, res_code)
}

/// RDKit's gradient scaling: multiply by 0.1, then cap so maxGrad <= 10.0.
/// Returns the total scaling factor (used in convergence test).
pub(crate) fn rdkit_scale_gradient(grad: &mut [f64]) -> f64 {
    let mut grad_scale = 0.1;
    let mut max_grad = -1e8f64;
    for g in grad.iter_mut() {
        *g *= grad_scale;
        if g.abs() > max_grad {
            max_grad = g.abs();
        }
    }
    if max_grad > 10.0 {
        while max_grad * grad_scale > 10.0 {
            grad_scale *= 0.5;
        }
        for g in grad.iter_mut() {
            *g *= grad_scale;
        }
    }
    grad_scale
}

/// Full BFGS minimizer matching RDKit's BFGSOpt::minimize exactly.
/// Uses f64 arithmetic, full N×N inverse Hessian, RDKit's line search,
/// and RDKit's gradient scaling hack.
/// Returns: 0=converged, 1=max iterations reached
pub fn minimize_bfgs_rdkit(
    coords: &mut DMatrix<f64>,
    bounds: &DMatrix<f64>,
    chiral_sets: &[ChiralSet],
    max_iters: usize,
    grad_tol: f64,
    basin_thresh: f32,
    weight_4d: f32,
    weight_chiral: f32,
) -> i32 {
    let n = coords.nrows();
    let dim_coords = coords.ncols();
    let dim = n * dim_coords;

    // Helper: flatten DMatrix<f64> to Vec<f64>
    let flatten = |c: &DMatrix<f64>| -> Vec<f64> {
        (0..n)
            .flat_map(|i| (0..dim_coords).map(move |d| c[(i, d)]))
            .collect()
    };
    let unflatten = |v: &[f64], c: &mut DMatrix<f64>| {
        for i in 0..n {
            for d in 0..dim_coords {
                c[(i, d)] = v[i * dim_coords + d];
            }
        }
    };

    // Energy function computed fully in f64 (matching RDKit's f64 precision)
    let basin_thresh_f64 = basin_thresh as f64;
    let weight_4d_f64 = weight_4d as f64;
    let weight_chiral_f64 = weight_chiral as f64;

    let calc_energy = |pos: &[f64]| -> f64 {
        let mut energy = 0.0f64;
        for i in 1..n {
            for j in 0..i {
                let ub = bounds[(j, i)];
                let lb = bounds[(i, j)];
                if ub - lb > basin_thresh_f64 { continue; }
                let mut d2 = 0.0f64;
                for d in 0..dim_coords {
                    let diff = pos[i * dim_coords + d] - pos[j * dim_coords + d];
                    d2 += diff * diff;
                }
                let ub2 = ub * ub;
                let lb2 = lb * lb;
                let val = if d2 > ub2 {
                    d2 / ub2 - 1.0
                } else if d2 < lb2 {
                    2.0 * lb2 / (lb2 + d2) - 1.0
                } else {
                    0.0
                };
                if val > 0.0 { energy += val * val; }
            }
        }
        if !chiral_sets.is_empty() {
            energy += weight_chiral_f64 * super::chiral_violation_energy_f64(pos, dim_coords, chiral_sets);
        }
        if dim_coords == 4 {
            for i in 0..n {
                let x4 = pos[i * dim_coords + 3];
                energy += weight_4d_f64 * x4 * x4;
            }
        }
        energy
    };

    // Gradient function computed fully in f64
    let calc_gradient_raw = |pos: &[f64]| -> Vec<f64> {
        let mut grad = vec![0.0f64; dim];
        for i in 1..n {
            for j in 0..i {
                let ub = bounds[(j, i)];
                let lb = bounds[(i, j)];
                if ub - lb > basin_thresh_f64 { continue; }
                let mut d2 = 0.0f64;
                let mut diffs = vec![0.0f64; dim_coords];
                for d in 0..dim_coords {
                    let diff = pos[i * dim_coords + d] - pos[j * dim_coords + d];
                    diffs[d] = diff;
                    d2 += diff * diff;
                }
                let ub2 = ub * ub;
                let lb2 = lb * lb;
                if d2 > ub2 {
                    let pre_factor = 4.0 * (d2 / ub2 - 1.0) / ub2;
                    for d in 0..dim_coords {
                        let g = pre_factor * diffs[d];
                        grad[i * dim_coords + d] += g;
                        grad[j * dim_coords + d] -= g;
                    }
                } else if d2 < lb2 {
                    let l2d2 = lb2 + d2;
                    let pre_factor = 8.0 * lb2 * (1.0 - 2.0 * lb2 / l2d2) / (l2d2 * l2d2);
                    for d in 0..dim_coords {
                        let g = pre_factor * diffs[d];
                        grad[i * dim_coords + d] += g;
                        grad[j * dim_coords + d] -= g;
                    }
                }
            }
        }
        if !chiral_sets.is_empty() {
            super::chiral_violation_gradient_f64(pos, dim_coords, chiral_sets, weight_chiral_f64, &mut grad);
        }
        if dim_coords == 4 {
            for i in 0..n {
                let x4 = pos[i * dim_coords + 3];
                grad[i * dim_coords + 3] += weight_4d_f64 * x4;
            }
        }
        grad
    };

    let mut pos = flatten(coords);
    let fp_init = calc_energy(&pos);
    let mut grad = calc_gradient_raw(&pos);
    let _grad_scale = rdkit_scale_gradient(&mut grad);

    // Initialize inverse Hessian to identity, direction to -grad
    let mut inv_hess = vec![0.0f64; dim * dim];
    let mut xi = vec![0.0f64; dim];
    let mut sum_pos = 0.0f64;
    for i in 0..dim {
        inv_hess[i * dim + i] = 1.0;
        xi[i] = -grad[i];
        sum_pos += pos[i] * pos[i];
    }
    let max_step = BFGS_MAXSTEP * sum_pos.sqrt().max(dim as f64);

    let mut fp = fp_init;
    let mut new_pos = vec![0.0f64; dim];
    let mut d_grad = vec![0.0f64; dim];
    let mut hess_d_grad = vec![0.0f64; dim];

    for _iter in 1..=max_iters {
        // line search
        let (found_pos, func_val, status) =
            rdkit_linear_search(&pos, fp, &grad, &mut xi, &calc_energy, max_step);
        debug_assert!(status >= 0, "bad direction in linear search");
        if status < 0 {
            break;
        }

        fp = func_val;
        new_pos.copy_from_slice(&found_pos);

        // compute direction of this line, check position convergence
        let mut test = 0.0f64;
        for i in 0..dim {
            xi[i] = new_pos[i] - pos[i];
            pos[i] = new_pos[i];
            let temp = xi[i].abs() / pos[i].abs().max(1.0);
            if temp > test {
                test = temp;
            }
            d_grad[i] = grad[i];
        }
        if test < BFGS_TOLX {
            unflatten(&pos, coords);
            return 0;
        }

        // update gradient
        grad = calc_gradient_raw(&pos);
        let grad_scale = rdkit_scale_gradient(&mut grad);

        // check gradient convergence
        test = 0.0;
        let term = (func_val * grad_scale).max(1.0);
        for i in 0..dim {
            let temp = grad[i].abs() * pos[i].abs().max(1.0);
            if temp > test {
                test = temp;
            }
            d_grad[i] = grad[i] - d_grad[i];
        }
        test /= term;
        if test < grad_tol {
            unflatten(&pos, coords);
            return 0;
        }

        // compute hessian * dGrad
        let mut fac = 0.0f64;
        let mut fae = 0.0f64;
        let mut sum_d_grad = 0.0f64;
        let mut sum_xi = 0.0f64;
        for i in 0..dim {
            hess_d_grad[i] = 0.0;
            for j in 0..dim {
                hess_d_grad[i] += inv_hess[i * dim + j] * d_grad[j];
            }
            fac += d_grad[i] * xi[i];
            fae += d_grad[i] * hess_d_grad[i];
            sum_d_grad += d_grad[i] * d_grad[i];
            sum_xi += xi[i] * xi[i];
        }

        if fac > (BFGS_EPS * sum_d_grad * sum_xi).sqrt() {
            let fac_inv = 1.0 / fac;
            let fad = 1.0 / fae;
            // compute difference vector for Hessian update
            for i in 0..dim {
                d_grad[i] = fac_inv * xi[i] - fad * hess_d_grad[i];
            }
            // update inverse Hessian (symmetric)
            for i in 0..dim {
                let pxi = fac_inv * xi[i];
                let hdgi = fad * hess_d_grad[i];
                let dgi = fae * d_grad[i];
                for j in i..dim {
                    inv_hess[i * dim + j] +=
                        pxi * xi[j]
                        - hdgi * hess_d_grad[j]
                        + dgi * d_grad[j];
                    inv_hess[j * dim + i] = inv_hess[i * dim + j];
                }
            }
        }

        // generate next direction: xi = -H^{-1} * grad
        for i in 0..dim {
            xi[i] = 0.0;
            for j in 0..dim {
                xi[i] -= inv_hess[i * dim + j] * grad[j];
            }
        }
    }

    unflatten(&pos, coords);
    1 // max iterations reached
}
