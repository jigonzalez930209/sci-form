use nalgebra::DMatrix;

pub fn triangle_smooth(b: &mut DMatrix<f64>) -> bool {
    triangle_smooth_tol(b, 0.0)
}

pub fn triangle_smooth_tol(b: &mut DMatrix<f64>, tol: f64) -> bool {
    let n = b.nrows();
    let feasible = true;
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
                let mut li = b[(j, i)];
                if li < d1 {
                    li = d1;
                } else if li < d2 {
                    li = d2;
                }
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
    feasible
}

pub fn smooth_bounds_matrix(mut b: DMatrix<f64>) -> DMatrix<f64> {
    triangle_smooth(&mut b);
    b
}
