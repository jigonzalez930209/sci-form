//! Triangle inequality smoothing for distance bounds matrices.
//!
//! Enforces the triangle inequality on all (i, j, k) triples:
//!   upper(i,j) ≤ upper(i,k) + upper(k,j)
//!   lower(i,j) ≥ lower(i,k) − upper(k,j)
//!
//! Returns `false` if the bounds are infeasible (lower > upper after tightening).

use nalgebra::DMatrix;

pub fn triangle_smooth(b: &mut DMatrix<f64>) -> bool {
    triangle_smooth_tol(b, 0.0)
}

pub fn triangle_smooth_tol(b: &mut DMatrix<f64>, tol: f64) -> bool {
    let n = b.nrows();
    for k in 0..n {
        for i in 0..n {
            if i == k {
                continue;
            }
            let (uk, lk) = if i < k {
                (b[(i, k)], b[(k, i)])
            } else {
                (b[(k, i)], b[(i, k)])
            };
            for j in (i + 1)..n {
                if j == k {
                    continue;
                }
                let (ukj, lkj) = if j < k {
                    (b[(j, k)], b[(k, j)])
                } else {
                    (b[(k, j)], b[(j, k)])
                };
                let (su, d1, d2) = (uk + ukj, lk - ukj, lkj - uk);
                if b[(i, j)] > su {
                    b[(i, j)] = su;
                }
                // Lower bound must satisfy BOTH triangle constraints
                let li = b[(j, i)].max(d1).max(d2);
                b[(j, i)] = li;
                let lower = b[(j, i)];
                let upper = b[(i, j)];
                if lower > upper {
                    if tol > 0.0 && lower > 1e-6 && (lower - upper) / lower < tol {
                        b[(i, j)] = lower; // set upper = lower
                    } else {
                        return false;
                    }
                }
            }
        }
    }
    true
}

pub fn smooth_bounds_matrix(mut b: DMatrix<f64>) -> DMatrix<f64> {
    triangle_smooth(&mut b);
    b
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contradictory_bounds_are_detected() {
        // 3 atoms with an impossible triple: d(0,1)=1, d(1,2)=1 but d(0,2)=5
        // violates the triangle inequality, so smoothing must report failure.
        let mut b = DMatrix::<f64>::zeros(3, 3);
        let set = |b: &mut DMatrix<f64>, i: usize, j: usize, lo: f64, hi: f64| {
            b[(i, j)] = hi;
            b[(j, i)] = lo;
        };
        set(&mut b, 0, 1, 1.0, 1.0);
        set(&mut b, 1, 2, 1.0, 1.0);
        set(&mut b, 0, 2, 5.0, 5.0);
        assert!(!triangle_smooth_tol(&mut b, 0.0));

        // With a tolerance the pair is clamped instead of rejected.
        let mut b2 = DMatrix::<f64>::zeros(3, 3);
        set(&mut b2, 0, 1, 1.0, 1.0);
        set(&mut b2, 1, 2, 1.0, 1.0);
        set(&mut b2, 0, 2, 2.000001, 2.000002);
        assert!(triangle_smooth_tol(&mut b2, 0.05));
    }

    #[test]
    fn test_consistent_bounds_smooth() {
        let mut b = DMatrix::<f64>::zeros(3, 3);
        b[(0, 1)] = 1.5;
        b[(1, 0)] = 1.4;
        b[(1, 2)] = 1.5;
        b[(2, 1)] = 1.4;
        b[(0, 2)] = 2.9;
        b[(2, 0)] = 2.8;
        assert!(triangle_smooth_tol(&mut b, 0.0));
    }
}
