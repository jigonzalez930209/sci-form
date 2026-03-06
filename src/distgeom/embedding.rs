use nalgebra::{DMatrix, SymmetricEigen};
use rand::Rng;

pub fn pick_random_distances<R: Rng>(rng: &mut R, bounds: &DMatrix<f32>) -> DMatrix<f32> {
    let n = bounds.nrows();
    let mut dists = DMatrix::from_element(n, n, 0.0);
    for i in 0..n {
        for j in (i + 1)..n {
            let u = bounds[(i, j)];
            let l = bounds[(j, i)];
            let d = if u > l { rng.gen_range(l..=u) } else { l };
            dists[(i, j)] = d;
            dists[(j, i)] = d;
        }
    }
    dists
}

pub fn compute_metric_matrix(dists: &DMatrix<f32>) -> DMatrix<f32> {
    let n = dists.nrows();
    let mut sq = DMatrix::from_element(n, n, 0.0);
    for i in 0..n {
        for j in 0..n {
            sq[(i, j)] = dists[(i, j)] * dists[(i, j)];
        }
    }
    let p =
        DMatrix::from_diagonal_element(n, n, 1.0) - DMatrix::from_element(n, n, 1.0 / (n as f32));
    -0.5 * (&p * &sq * &p)
}

pub fn generate_nd_coordinates<R: Rng>(
    rng: &mut R,
    mm: &DMatrix<f32>,
    ndim: usize,
) -> DMatrix<f32> {
    let n = mm.nrows();
    let eigen = SymmetricEigen::new(mm.clone());
    let mut evals: Vec<(usize, f32)> = eigen.eigenvalues.iter().copied().enumerate().collect();
    evals.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut coords = DMatrix::from_element(n, ndim, 0.0);
    for d in 0..ndim {
        if d < evals.len() {
            let (idx, val) = evals[d];
            if val > 1e-4 {
                let root = val.sqrt();
                let evec = eigen.eigenvectors.column(idx);
                for i in 0..n {
                    coords[(i, d)] = evec[i] * root;
                }
            } else {
                for i in 0..n {
                    let r: f32 = rng.gen_range(-1.0..1.0);
                    coords[(i, d)] = r * 1e-4;
                }
            }
        } else {
            for i in 0..n {
                let r: f32 = rng.gen_range(-1.0..1.0);
                coords[(i, d)] = r * 1e-4;
            }
        }
    }
    coords
}

pub fn generate_3d_coordinates<R: Rng>(rng: &mut R, mm: &DMatrix<f32>) -> DMatrix<f32> {
    generate_nd_coordinates(rng, mm, 3)
}
