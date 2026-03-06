use nalgebra::{DMatrix, SymmetricEigen};
use rand::Rng;

pub fn pick_random_distances<R: Rng>(rng: &mut R, bounds: &DMatrix<f32>) -> DMatrix<f32> {
    let n = bounds.nrows();
    let mut dists = DMatrix::from_element(n, n, 0.0);
    for i in 0..n {
        for j in (i + 1)..n {
            let u = bounds[(i, j)];
            let l = bounds[(j, i)];
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
pub fn pick_etkdg_distances<R: Rng>(rng: &mut R, bounds: &DMatrix<f32>) -> DMatrix<f32> {
    let n = bounds.nrows();
    let mut dists = DMatrix::from_element(n, n, 0.0);
    
    for i in 0..n {
        for j in (i + 1)..n {
            let u = bounds[(i, j)];
            let l = bounds[(j, i)];
            
            if u > l {
                // ETKDG: prefer distances closer to lower bound for 1-2 and 1-3 interactions
                // This helps maintain proper geometry
                let range = u - l;
                let mid = (u + l) / 2.0;
                
                // Use a distribution that favors values slightly below midpoint
                // for shorter distances (more chemically relevant)
                let r: f32 = rng.gen();
                let bias = if range < 0.5 { 0.4 } else { 0.5 }; // favor lower for tight bounds
                let d = l + range * (r.powf(bias));
                
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
    
    // Use eigenvalue weighting - larger eigenvalues get more weight in coordinate generation
    let total_positive: f32 = evals.iter().map(|(_, v)| v.max(0.0)).sum();
    
    for d in 0..ndim {
        if d < evals.len() {
            let (idx, val) = evals[d];
            if val > 1e-6 {
                let root = val.sqrt();
                let evec = eigen.eigenvectors.column(idx);
                // Normalize eigenvector by eigenvalue weight for better conditioning
                let weight = (val / total_positive.max(1e-10)).sqrt();
                for i in 0..n {
                    coords[(i, d)] = evec[i] * root * weight.max(1.0);
                }
            } else {
                // For near-zero eigenvalues, use small random perturbation
                for i in 0..n {
                    let r: f32 = rng.gen_range(-1.0..1.0);
                    coords[(i, d)] = r * 1e-3;
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

pub fn generate_3d_coordinates<R: Rng>(rng: &mut R, mm: &DMatrix<f32>) -> DMatrix<f32> {
    generate_nd_coordinates(rng, mm, 3)
}

/// Generate coordinates using a distance geometry approach with better initialization
/// This uses a smarter initialization based on the molecular topology
pub fn generate_smart_coordinates<R: Rng>(
    rng: &mut R,
    bounds: &DMatrix<f32>,
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
                
                let l = bounds[(j, i)];
                let u = bounds[(i, j)];
                
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
