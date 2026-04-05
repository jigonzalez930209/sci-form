//! Image Dependent Pair Potential (IDPP) initial path generation.
//!
//! Instead of naive linear interpolation in Cartesian space (which causes
//! atom overlaps and unphysical geometries), IDPP generates an initial
//! NEB path where each image has interatomic distances that smoothly
//! interpolate between reactant and product values.
//!
//! Reference: Smidstrup, Pedersen & Jónsson, J. Chem. Phys. 140, 214106 (2014)

/// Generate an IDPP-interpolated initial path between two geometries.
///
/// For each image i at parameter t_i = i/(N-1):
///   1. Compute target pairwise distances: d_ij(t) = (1-t)*d_ij(R) + t*d_ij(P)
///   2. Minimize IDPP objective: Σ_{i<j} w_ij * (|r_i - r_j| - d_ij(t))²
///      where w_ij = 1/d_ij(t)^4 (penalizes short distances more)
///   3. This produces geometries that respect interatomic distances
///      without atom clashes.
///
/// Returns `n_images` coordinate sets (each a flat `Vec<f64>`).
pub fn generate_idpp_path(
    start: &[f64],
    end: &[f64],
    n_images: usize,
    max_steps: usize,
) -> Vec<Vec<f64>> {
    let n_xyz = start.len();
    let n_atoms = n_xyz / 3;
    if n_atoms < 2 || n_images < 2 {
        // Fallback: linear interpolation
        return linear_interpolation(start, end, n_images);
    }

    // Compute pairwise distances for reactant and product
    let d_start = pairwise_distances(start, n_atoms);
    let d_end = pairwise_distances(end, n_atoms);

    // Start from linear interpolation
    let mut images = linear_interpolation(start, end, n_images);

    // Optimize interior images to match IDPP target distances
    let step_size = 0.04;

    for _iter in 0..max_steps {
        let prev = images.clone();
        let mut max_shift = 0.0f64;

        for img_idx in 1..(n_images - 1) {
            let t = img_idx as f64 / (n_images - 1) as f64;

            // Target distances at this parameter value
            let d_target: Vec<f64> = d_start
                .iter()
                .zip(d_end.iter())
                .map(|(&ds, &de)| (1.0 - t) * ds + t * de)
                .collect();

            // Compute IDPP gradient for this image
            let grad = idpp_gradient(&prev[img_idx], n_atoms, &d_target);

            // Steepest descent step
            for k in 0..n_xyz {
                images[img_idx][k] = prev[img_idx][k] - step_size * grad[k];
                let shift = (images[img_idx][k] - prev[img_idx][k]).abs();
                if shift > max_shift {
                    max_shift = shift;
                }
            }
        }

        // Converged if images barely moved
        if max_shift < 1e-4 {
            break;
        }
    }

    images
}

/// Compute pairwise distance matrix (flattened upper triangle, row-major).
fn pairwise_distances(coords: &[f64], n_atoms: usize) -> Vec<f64> {
    let n_pairs = n_atoms * (n_atoms - 1) / 2;
    let mut dists = Vec::with_capacity(n_pairs);
    for i in 0..n_atoms {
        for j in (i + 1)..n_atoms {
            let dx = coords[j * 3] - coords[i * 3];
            let dy = coords[j * 3 + 1] - coords[i * 3 + 1];
            let dz = coords[j * 3 + 2] - coords[i * 3 + 2];
            dists.push((dx * dx + dy * dy + dz * dz).sqrt());
        }
    }
    dists
}

/// Compute gradient of the IDPP objective function.
///
/// IDPP(x) = Σ_{i<j} w_ij * (|r_i - r_j| - d_target_ij)²
/// w_ij = 1 / d_target_ij^4
///
/// ∂IDPP/∂r_i = Σ_{j≠i} 2 * w_ij * (|r_i - r_j| - d_target_ij) * (r_i - r_j) / |r_i - r_j|
fn idpp_gradient(coords: &[f64], n_atoms: usize, d_target: &[f64]) -> Vec<f64> {
    let n_xyz = n_atoms * 3;
    let mut grad = vec![0.0f64; n_xyz];
    let mut pair_idx = 0;

    for i in 0..n_atoms {
        for j in (i + 1)..n_atoms {
            let dx = coords[i * 3] - coords[j * 3];
            let dy = coords[i * 3 + 1] - coords[j * 3 + 1];
            let dz = coords[i * 3 + 2] - coords[j * 3 + 2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();

            let d_tgt = d_target[pair_idx];
            pair_idx += 1;

            if dist < 1e-8 || d_tgt < 1e-8 {
                continue;
            }

            // Weight: penalizes close contacts more heavily
            let w = 1.0 / (d_tgt * d_tgt * d_tgt * d_tgt);
            let factor = 2.0 * w * (dist - d_tgt) / dist;

            grad[i * 3] += factor * dx;
            grad[i * 3 + 1] += factor * dy;
            grad[i * 3 + 2] += factor * dz;
            grad[j * 3] -= factor * dx;
            grad[j * 3 + 1] -= factor * dy;
            grad[j * 3 + 2] -= factor * dz;
        }
    }

    grad
}

/// Simple linear interpolation in Cartesian space (fallback).
fn linear_interpolation(start: &[f64], end: &[f64], n_images: usize) -> Vec<Vec<f64>> {
    let n_xyz = start.len();
    (0..n_images)
        .map(|i| {
            let t = if n_images > 1 {
                i as f64 / (n_images - 1) as f64
            } else {
                0.0
            };
            (0..n_xyz)
                .map(|k| (1.0 - t) * start[k] + t * end[k])
                .collect()
        })
        .collect()
}
