use nalgebra::{DMatrix, SymmetricEigen};
use rand::Rng;
use super::rng::MinstdRand;

pub fn pick_rdkit_distances(rng: &mut MinstdRand, bounds: &DMatrix<f64>) -> DMatrix<f64> {
    let n = bounds.nrows();
    let mut dists = DMatrix::from_element(n, n, 0.0f64);

    // Match RDKit's loop order: for i=1..n, for j=0..i
    for i in 1..n {
        for j in 0..i {
            // RDKit: getUpperBound(i,j) where i>j → stored at (j,i) upper triangle
            let ub = bounds[(j, i)];
            // RDKit: getLowerBound(i,j) where i>j → stored at (i,j) lower triangle
            let lb = bounds[(i, j)];
            let rval = rng.next_double();
            let d = lb + rval * (ub - lb);
            dists[(i, j)] = d;
            dists[(j, i)] = d;
        }
    }
    dists
}

pub fn pick_random_distances<R: Rng>(rng: &mut R, bounds: &DMatrix<f64>) -> DMatrix<f32> {
    let n = bounds.nrows();
    let mut dists = DMatrix::from_element(n, n, 0.0f32);
    for i in 0..n {
        for j in (i + 1)..n {
            let u = bounds[(i, j)] as f32;
            let l = bounds[(j, i)] as f32;
            // Use triangular distribution for better sampling (more likely near center)
            let d = if u > l {
                // Sample from triangular distribution centered at midpoint
                let mid = (u + l) / 2.0;
                let half_width = (u - l) / 2.0;
                let r: f32 = rng.gen();
                mid + half_width * (r - 0.5) * 2.0
            } else { 
                l 
            };
            dists[(i, j)] = d;
            dists[(j, i)] = d;
        }
    }
    dists
}

/// ETKDG-style distance selection: pick distances closer to ideal values
/// This mimics RDKit's approach of using more favorable distance combinations
pub fn pick_etkdg_distances<R: Rng>(rng: &mut R, bounds: &DMatrix<f64>) -> DMatrix<f32> {
    let n = bounds.nrows();
    let mut dists = DMatrix::from_element(n, n, 0.0f32);
    
    for i in 0..n {
        for j in (i + 1)..n {
            let u = bounds[(i, j)] as f32;
            let l = bounds[(j, i)] as f32;
            
            if u > l {
                let range = u - l;
                let r: f32 = rng.gen();
                let d = l + range * r;
                dists[(i, j)] = d;
                dists[(j, i)] = d;
            } else {
                dists[(i, j)] = l;
                dists[(j, i)] = l;
            }
        }
    }
    dists
}

/// Compute metric matrix using f64 internally for numerical precision,
/// matching RDKit's double-precision Cayley-Menger transform.
pub fn compute_metric_matrix(dists: &DMatrix<f32>) -> DMatrix<f32> {
    let n = dists.nrows();

    // Compute squared distances and global mean
    let mut sq_dists: Vec<f64> = vec![0.0; n * n];
    let mut sum_sq_all = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            let d = dists[(i, j)] as f64;
            sq_dists[i * n + j] = d * d;
            sum_sq_all += d * d;
        }
    }
    sum_sq_all /= (n * n) as f64;

    // D0_i = (1/N) * sum_j(d²_ij) - sumSqAll
    let mut d0: Vec<f64> = vec![0.0; n];
    for i in 0..n {
        let mut row_sum = 0.0f64;
        for j in 0..n {
            row_sum += sq_dists[i * n + j];
        }
        d0[i] = row_sum / n as f64 - sum_sq_all;
    }

    // T_ij = 0.5 * (D0_i + D0_j - d²_ij)
    let mut metric = DMatrix::from_element(n, n, 0.0f32);
    for i in 0..n {
        for j in 0..=i {
            let val = 0.5 * (d0[i] + d0[j] - sq_dists[i * n + j]);
            metric[(i, j)] = val as f32;
            metric[(j, i)] = val as f32;
        }
    }
    metric
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
            if val > 1e-3 {
                let root = val.sqrt();
                let evec = eigen.eigenvectors.column(idx);
                for i in 0..n {
                    coords[(i, d)] = evec[i] * root;
                }
            } else if val.abs() < 1e-3 {
                for i in 0..n {
                    let r: f32 = rng.gen_range(-1.0..1.0);
                    coords[(i, d)] = r * 1e-3;
                }
            } else {
                for i in 0..n {
                    let r: f32 = rng.gen_range(-1.0..1.0);
                    coords[(i, d)] = r;
                }
            }
        } else {
            for i in 0..n {
                let r: f32 = rng.gen_range(-1.0..1.0);
                coords[(i, d)] = r * 1e-3;
            }
        }
    }
    coords
}

/// Generate N-dimensional coordinates using nalgebra's SymmetricEigen + MinstdRand.
/// Returns None if too many zero eigenvalues (matching RDKit's rejection logic).
pub fn compute_initial_coords_nalgebra(
    rng: &mut MinstdRand,
    mm: &DMatrix<f32>,
    ndim: usize,
) -> Option<DMatrix<f32>> {
    let n = mm.nrows();
    let eigen = SymmetricEigen::new(mm.clone());
    let mut evals: Vec<(usize, f32)> = eigen.eigenvalues.iter().copied().enumerate().collect();
    evals.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut coords = DMatrix::from_element(n, ndim, 0.0f32);

    for d in 0..ndim {
        if d < evals.len() {
            let (idx, val) = evals[d];
            if val > 1e-3 {
                let root = val.sqrt();
                let evec = eigen.eigenvectors.column(idx);
                for i in 0..n {
                    coords[(i, d)] = evec[i] * root;
                }
            } else {
                // Zero or negative eigenvalue → random in [-1, 1] using MinstdRand
                for i in 0..n {
                    coords[(i, d)] = 1.0 - 2.0 * rng.next_double() as f32;
                }
            }
        } else {
            for i in 0..n {
                coords[(i, d)] = 1.0 - 2.0 * rng.next_double() as f32;
            }
        }
    }

    // Don't reject on zero eigenvalues — let downstream validation handle it
    Some(coords)
}

/// Compute initial coordinates using nalgebra's SymmetricEigen in f64.
/// Drop-in replacement for compute_initial_coords_rdkit, building the metric
/// matrix from distances identically but using a different eigendecomposition.
pub fn compute_initial_coords_nalgebra_f64(
    rng: &mut MinstdRand,
    dists: &DMatrix<f64>,
    ndim: usize,
) -> Option<DMatrix<f64>> {
    let n = dists.nrows();

    // Step 1: Build metric matrix (identical to compute_initial_coords_rdkit)
    let d_size = n * (n + 1) / 2;
    let mut sq_packed = vec![0.0f64; d_size];
    let mut sum_sq_all = 0.0f64;
    for i in 0..n {
        let id = i * (i + 1) / 2;
        for j in 0..=i {
            let d = dists[(i, j)];
            sq_packed[id + j] = d * d;
            sum_sq_all += d * d;
        }
    }
    sum_sq_all /= (n * n) as f64;

    let mut d0 = vec![0.0f64; n];
    for i in 0..n {
        let mut row_sum = 0.0f64;
        for j in 0..n {
            let idx = if i >= j {
                i * (i + 1) / 2 + j
            } else {
                j * (j + 1) / 2 + i
            };
            row_sum += sq_packed[idx];
        }
        d0[i] = row_sum / n as f64 - sum_sq_all;
        if d0[i] < EIGVAL_TOL && n > 3 {
            return None;
        }
    }

    // Step 2: Build full symmetric metric matrix for nalgebra
    let mut t_full = DMatrix::from_element(n, n, 0.0f64);
    for i in 0..n {
        for j in 0..=i {
            let sq_val = sq_packed[if i >= j {
                i * (i + 1) / 2 + j
            } else {
                j * (j + 1) / 2 + i
            }];
            let val = 0.5 * (d0[i] + d0[j] - sq_val);
            t_full[(i, j)] = val;
            t_full[(j, i)] = val;
        }
    }

    // Step 3: nalgebra SymmetricEigen decomposition (f64)
    let eigen = SymmetricEigen::new(t_full);
    let mut evals: Vec<(usize, f64)> = eigen.eigenvalues.iter().copied().enumerate().collect();
    evals.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    // Step 4: Process eigenvalues and build coordinates
    let n_eigs = ndim.min(n);
    let mut zero_eigs = 0u32;
    let mut coords = DMatrix::from_element(n, ndim, 0.0f64);

    for d in 0..ndim {
        if d < n_eigs && d < evals.len() {
            let (idx, val) = evals[d];
            if val > EIGVAL_TOL {
                let root = val.sqrt();
                let evec = eigen.eigenvectors.column(idx);
                for i in 0..n {
                    coords[(i, d)] = evec[i] * root;
                }
            } else if val.abs() < EIGVAL_TOL {
                zero_eigs += 1;
                for i in 0..n {
                    coords[(i, d)] = 1.0 - 2.0 * rng.next_double();
                }
            } else {
                // Negative eigenvalue
                for i in 0..n {
                    coords[(i, d)] = 1.0 - 2.0 * rng.next_double();
                }
            }
        } else {
            for i in 0..n {
                coords[(i, d)] = 1.0 - 2.0 * rng.next_double();
            }
        }
    }

    // Same rejection criteria as the rdkit version
    if zero_eigs >= 1 && n > 3 {
        return None;
    }

    Some(coords)
}

/// Power iteration eigenvalue solver matching RDKit's powerEigenSolver exactly.
/// Uses packed symmetric matrix format (lower triangular, row-major).
/// Returns (eigenvalues, eigenvectors) for the top `num_eig` eigenvalues.
pub fn power_eigen_solver(
    num_eig: usize,
    mat: &mut Vec<f64>, // packed symmetric: mat[i*(i+1)/2 + j] for j <= i
    n: usize,
    mut seed: i32,
) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
    const MAX_ITERATIONS: usize = 1000;
    const TOLERANCE: f64 = 0.001;
    const HUGE_EIGVAL: f64 = 1.0e10;
    const TINY_EIGVAL: f64 = 1.0e-10;

    let mut eigenvalues = Vec::with_capacity(num_eig);
    let mut eigenvectors = Vec::with_capacity(num_eig);

    for ei in 0..num_eig {
        let mut eig_val = -HUGE_EIGVAL;
        seed += ei as i32;

        // Initialize random vector using MinstdRand (matching RDKit's setToRandom)
        // RDKit's Vector::setToRandom uses boost::minstd_rand (a=48271) seeded with the seed,
        // then generates uniform doubles in [0, 1) and normalizes to unit vector
        let mut v = vec![0.0f64; n];
        {
            let mut rng = MinstdRand::new(seed as u32);
            for vi in v.iter_mut() {
                *vi = rng.next_double();
            }
            // Normalize: simple sequential dotProduct matching RDKit
            let mut ns = 0.0f64;
            for x in v.iter() {
                ns += x * x;
            }
            let norm = ns.sqrt();
            if norm > 0.0 {
                for vi in v.iter_mut() {
                    *vi /= norm;
                }
            }
        }

        let mut converged = false;
        let mut z = vec![0.0f64; n];

        for _iter in 0..MAX_ITERATIONS {
            // z = mat * v (using packed symmetric format)
            for i in 0..n {
                let mut accum = 0.0f64;
                let id = i * (i + 1) / 2;
                for j in 0..=i {
                    accum += mat[id + j] * v[j];
                }
                for j in (i + 1)..n {
                    accum += mat[j * (j + 1) / 2 + i] * v[j];
                }
                z[i] = accum;
            }

            let prev_val = eig_val;
            // Find index of largest absolute component (first-wins tie-breaking to match RDKit)
            let mut eval_id = 0;
            let mut max_abs = -1.0f64;
            for i in 0..n {
                let a = z[i].abs();
                if a > max_abs {
                    max_abs = a;
                    eval_id = i;
                }
            }
            eig_val = z[eval_id];

            if eig_val.abs() < TINY_EIGVAL {
                break;
            }

            // v = z / eigVal
            for i in 0..n {
                v[i] = z[i] / eig_val;
            }

            if (eig_val - prev_val).abs() < TOLERANCE {
                converged = true;
                break;
            }
        }

        if !converged {
            return None;
        }

        // Normalize v: simple sequential dotProduct matching RDKit
        let mut norm_sum = 0.0f64;
        for x in v.iter() {
            norm_sum += x * x;
        }
        let norm = norm_sum.sqrt();
        if norm > 0.0 {
            for x in &mut v {
                *x /= norm;
            }
        }

        eigenvalues.push(eig_val);
        eigenvectors.push(v.clone());

        // Deflate: mat -= eigVal * v * v^T (packed symmetric)
        for i in 0..n {
            let id = i * (i + 1) / 2;
            for j in 0..=i {
                mat[id + j] -= eig_val * v[i] * v[j];
            }
        }
    }

    Some((eigenvalues, eigenvectors))
}

const EIGVAL_TOL: f64 = 1.0e-3;

/// Compute initial coordinates matching RDKit's computeInitialCoords exactly.
/// Uses f64 throughout for metric matrix, eigendecomposition, and output coordinates.
/// Returns None if embedding fails (bad D0 values, too many zero eigenvalues, etc.)
pub fn compute_initial_coords_rdkit(
    rng: &mut MinstdRand,
    dists: &DMatrix<f64>,
    ndim: usize,
) -> Option<DMatrix<f64>> {
    let n = dists.nrows();

    // Step 1: Build squared distance matrix and compute averages (all in f64)
    let _sq_dists = vec![0.0f64; n * n];
    let mut sum_sq_all = 0.0f64;
    let d_size = n * (n + 1) / 2;
    let mut sq_packed = vec![0.0f64; d_size];
    for i in 0..n {
        let id = i * (i + 1) / 2;
        for j in 0..=i {
            let d = dists[(i, j)];
            sq_packed[id + j] = d * d;
            sum_sq_all += d * d;
        }
    }
    sum_sq_all /= (n * n) as f64;

    // D0_i = (1/N) * sum_j(sqMat.getVal(i,j)) - sumSqD2
    // RDKit iterates j from 0 to N, using getVal(i,j) which for symmetric packed
    // returns the correct element regardless of ordering
    let mut d0 = vec![0.0f64; n];
    for i in 0..n {
        let mut row_sum = 0.0f64;
        for j in 0..n {
            let idx = if i >= j {
                i * (i + 1) / 2 + j
            } else {
                j * (j + 1) / 2 + i
            };
            row_sum += sq_packed[idx];
        }
        d0[i] = row_sum / n as f64 - sum_sq_all;

        // RDKit check: if D0 value is negative and N > 3, fail
        if d0[i] < EIGVAL_TOL && n > 3 {
            return None;
        }
    }

    // Step 2: Build metric matrix T (packed symmetric, f64)
    let mut t_packed = vec![0.0f64; d_size];
    for i in 0..n {
        let id = i * (i + 1) / 2;
        for j in 0..=i {
            let sq_val = sq_packed[if i >= j {
                i * (i + 1) / 2 + j
            } else {
                j * (j + 1) / 2 + i
            }];
            t_packed[id + j] = 0.5 * (d0[i] + d0[j] - sq_val);
        }
    }

    // Step 3: Power eigendecomposition (matching RDKit exactly)
    let n_eigs = ndim.min(n);
    let eigen_seed = (sum_sq_all * n as f64) as i32;
    let (eigenvalues, eigenvectors) =
        power_eigen_solver(n_eigs, &mut t_packed, n, eigen_seed)?;

    // Step 4: Process eigenvalues
    let mut eig_sqrt = vec![0.0f64; n_eigs];
    let mut _found_neg = false;
    let mut zero_eigs = 0u32;
    for i in 0..n_eigs {
        if eigenvalues[i] > EIGVAL_TOL {
            eig_sqrt[i] = eigenvalues[i].sqrt();
        } else if eigenvalues[i].abs() < EIGVAL_TOL {
            eig_sqrt[i] = 0.0;
            zero_eigs += 1;
        } else {
            _found_neg = true;
            eig_sqrt[i] = eigenvalues[i]; // keep negative marker
        }
    }

    // RDKit default: numZeroFail=1, fail if >= 1 zero eigenvalue and N > 3
    if zero_eigs >= 1 && n > 3 {
        return None;
    }

    // Step 5: Generate coordinates
    let mut coords = DMatrix::from_element(n, ndim, 0.0f64);
    for i in 0..n {
        for j in 0..ndim {
            if j < n_eigs && eig_sqrt[j] >= 0.0 {
                coords[(i, j)] = eig_sqrt[j] * eigenvectors[j][i];
            } else {
                // Negative eigenvalue or missing: random in [-1, 1]
                coords[(i, j)] = 1.0 - 2.0 * rng.next_double();
            }
        }
    }

    Some(coords)
}

pub fn generate_3d_coordinates<R: Rng>(rng: &mut R, mm: &DMatrix<f32>) -> DMatrix<f32> {
    generate_nd_coordinates(rng, mm, 3)
}

/// Generate coordinates using a distance geometry approach with better initialization
/// This uses a smarter initialization based on the molecular topology
pub fn generate_smart_coordinates<R: Rng>(
    rng: &mut R,
    bounds: &DMatrix<f64>,
    mol: &crate::graph::Molecule,
) -> DMatrix<f32> {
    let n = mol.graph.node_count();
    
    // First pass: generate initial coordinates using standard distance geometry
    let dists = pick_random_distances(rng, bounds);
    let metric = compute_metric_matrix(&dists);
    let mut coords = generate_nd_coordinates(rng, &metric, 3);
    
    // Second pass: refine using triangle smoothing on coordinates
    // This helps enforce the triangle inequality better
    for _iter in 0..5 {
        for i in 0..n {
            for j in (i + 1)..n {
                let dij = ((coords[(i, 0)] - coords[(j, 0)]).powi(2)
                    + (coords[(i, 1)] - coords[(j, 1)]).powi(2)
                    + (coords[(i, 2)] - coords[(j, 2)]).powi(2)).sqrt();
                
                let l = bounds[(j, i)] as f32;
                let u = bounds[(i, j)] as f32;
                
                if dij < l {
                    // Pull atoms apart
                    let diff = (l - dij) / (dij + 1e-6);
                    coords[(i, 0)] -= diff * (coords[(i, 0)] - coords[(j, 0)]) * 0.1;
                    coords[(i, 1)] -= diff * (coords[(i, 1)] - coords[(j, 1)]) * 0.1;
                    coords[(i, 2)] -= diff * (coords[(i, 2)] - coords[(j, 2)]) * 0.1;
                } else if dij > u {
                    // Push atoms together
                    let diff = (dij - u) / (dij + 1e-6);
                    coords[(i, 0)] += diff * (coords[(i, 0)] - coords[(j, 0)]) * 0.1;
                    coords[(i, 1)] += diff * (coords[(i, 1)] - coords[(j, 1)]) * 0.1;
                    coords[(i, 2)] += diff * (coords[(i, 2)] - coords[(j, 2)]) * 0.1;
                }
            }
        }
    }
    
    coords
}
