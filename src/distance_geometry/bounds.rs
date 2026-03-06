use thiserror::Error;

#[derive(Error, Debug)]
pub enum GeometryBoundsError {
    #[error("Infranqueable violación fatal de desigualdad triangular reportada entre los nodos {0} y {1}.")]
    TriangleInequalityViolation(usize, usize),
}

pub fn optimized_triangle_bounds_smoothing(
    num_atoms: usize,
    upper_bounds: &mut [f64],
    lower_bounds: &mut [f64],
    tol_margen: f64,
) -> Result<(), GeometryBoundsError> {
    let n = num_atoms;

    for k in 0..n {
        for i in 0..n {
            if i == k {
                continue;
            }

            let ik_upper = upper_bounds[i * n + k];
            let ik_lower = lower_bounds[i * n + k];

            for j in 0..n {
                if i == j || j == k {
                    continue;
                }

                let idx_ij = i * n + j;
                let idx_kj = k * n + j;

                let kj_upper = upper_bounds[idx_kj];
                let kj_lower = lower_bounds[idx_kj];

                let target_upper = ik_upper + kj_upper;
                if target_upper < upper_bounds[idx_ij] {
                    upper_bounds[idx_ij] = target_upper;
                }

                let target_lower_1 = ik_lower - kj_upper;
                let target_lower_2 = kj_lower - ik_upper;
                let best_lower = target_lower_1.max(target_lower_2);

                if best_lower > lower_bounds[idx_ij] {
                    lower_bounds[idx_ij] = best_lower;
                }

                if lower_bounds[idx_ij] > upper_bounds[idx_ij] + tol_margen {
                    return Err(GeometryBoundsError::TriangleInequalityViolation(i, j));
                }
            }
        }
    }
    Ok(())
}
