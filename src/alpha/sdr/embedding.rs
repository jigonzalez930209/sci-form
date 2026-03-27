//! SDR Embedding integration — E8.2
//!
//! Coordinate extraction from SDR-relaxed Gram matrices and
//! integration with the conformer pipeline.

use super::projections::{alternating_projections, SdrConfig, SdrConvergence};
use nalgebra::DMatrix;

/// Result from SDR embedding.
#[derive(Debug, Clone)]
pub struct SdrResult {
    /// Extracted 3D coordinates (flat: [x0,y0,z0, x1,y1,z1, ...]).
    pub coords: Vec<f64>,
    /// Number of atoms.
    pub num_atoms: usize,
    /// Convergence information.
    pub convergence: SdrConvergence,
    /// Maximum distance error |d_SDR - d_target| in Å.
    pub max_distance_error: f64,
    /// Number of retries that SDR would have saved.
    pub retries_avoided: usize,
}

/// Build warm-start Gram matrix from distance matrix (double-centering).
///
/// Given squared distances D, compute the Gram matrix:
/// X = -0.5 * J * D * J where J = I - (1/n) * 11^T (centering matrix)
pub fn warm_start_gram(n: usize, distance_pairs: &[(usize, usize, f64)]) -> DMatrix<f64> {
    // Build full distance squared matrix
    let mut d_sq = DMatrix::zeros(n, n);
    for &(i, j, d) in distance_pairs {
        d_sq[(i, j)] = d * d;
        d_sq[(j, i)] = d * d;
    }

    // Double-centering: X = -0.5 * J * D² * J
    let one_n = 1.0 / n as f64;
    let mut gram = DMatrix::zeros(n, n);

    for i in 0..n {
        for j in 0..n {
            let val = d_sq[(i, j)];

            // Row mean
            let row_mean: f64 = (0..n).map(|k| d_sq[(i, k)]).sum::<f64>() * one_n;
            // Column mean
            let col_mean: f64 = (0..n).map(|k| d_sq[(k, j)]).sum::<f64>() * one_n;
            // Grand mean
            let grand_mean: f64 = d_sq.iter().sum::<f64>() * one_n * one_n;

            gram[(i, j)] = -0.5 * (val - row_mean - col_mean + grand_mean);
        }
    }

    gram
}

/// Extract 3D coordinates from a PSD Gram matrix.
///
/// Uses the top-3 eigenvalues/eigenvectors: C = V_3 * Σ_3^{1/2}
pub fn extract_coordinates(gram: &DMatrix<f64>) -> Vec<f64> {
    let n = gram.nrows();
    let eigen = gram.clone().symmetric_eigen();

    // Sort eigenvalues descending
    let mut indexed: Vec<(usize, f64)> = eigen.eigenvalues.iter().cloned().enumerate().collect();
    indexed.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    // Take top 3
    let mut coords = vec![0.0; n * 3];
    for dim in 0..3.min(indexed.len()) {
        let (idx, eval) = indexed[dim];
        let sqrt_eval = if eval > 0.0 { eval.sqrt() } else { 0.0 };

        for i in 0..n {
            coords[i * 3 + dim] = eigen.eigenvectors[(i, idx)] * sqrt_eval;
        }
    }

    coords
}

/// Run full SDR embedding: warm-start → alternating projections → extract coordinates.
pub fn sdr_embed(
    n: usize,
    distance_pairs: &[(usize, usize, f64)],
    config: &SdrConfig,
) -> SdrResult {
    // Warm-start from distance matrix
    let x0 = warm_start_gram(n, distance_pairs);

    // Alternating projections
    let (gram, convergence) = alternating_projections(&x0, distance_pairs, config);

    // Extract 3D coordinates
    let coords = extract_coordinates(&gram);

    // Compute max distance error
    let mut max_err = 0.0;
    let mut retries = 0;

    for &(i, j, d_target) in distance_pairs {
        let dx = coords[i * 3] - coords[j * 3];
        let dy = coords[i * 3 + 1] - coords[j * 3 + 1];
        let dz = coords[i * 3 + 2] - coords[j * 3 + 2];
        let d_actual = (dx * dx + dy * dy + dz * dz).sqrt();
        let err = (d_actual - d_target).abs();
        if err > max_err {
            max_err = err;
        }
    }

    // Count negative eigenvalues in original warm-start
    let eigen_orig = x0.symmetric_eigen();
    for &e in eigen_orig.eigenvalues.iter() {
        if e < -1e-8 {
            retries += 1;
        }
    }

    SdrResult {
        coords,
        num_atoms: n,
        convergence,
        max_distance_error: max_err,
        retries_avoided: retries,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_warm_start_gram_simple() {
        // 3 points forming a right triangle
        let pairs = vec![
            (0, 1, 1.0),           // d=1
            (0, 2, 1.0),           // d=1
            (1, 2, 2.0f64.sqrt()), // d=√2
        ];
        let gram = warm_start_gram(3, &pairs);
        assert_eq!(gram.nrows(), 3);
        // Gram matrix should be approximately PSD for valid distances
    }

    #[test]
    fn test_extract_coordinates_from_known() {
        // Build Gram from known coords
        let known = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let n = 3;
        let mut gram = DMatrix::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                gram[(i, j)] = known[i][0] * known[j][0]
                    + known[i][1] * known[j][1]
                    + known[i][2] * known[j][2];
            }
        }

        let coords = extract_coordinates(&gram);
        assert_eq!(coords.len(), 9);

        // Verify distances are preserved
        let d01 = ((coords[0] - coords[3]).powi(2)
            + (coords[1] - coords[4]).powi(2)
            + (coords[2] - coords[5]).powi(2))
        .sqrt();
        assert!((d01 - 1.0).abs() < 0.1, "d01 should be ~1.0: {}", d01);
    }

    #[test]
    fn test_sdr_embed_triangle() {
        let pairs = vec![(0, 1, 1.5), (0, 2, 1.5), (1, 2, 1.5)];
        let config = SdrConfig::default();
        let result = sdr_embed(3, &pairs, &config);
        assert_eq!(result.num_atoms, 3);
        assert_eq!(result.coords.len(), 9);
    }

    #[test]
    fn test_sdr_embed_tetrahedron() {
        // Regular tetrahedron: all distances equal
        let d = 2.0;
        let pairs = vec![
            (0, 1, d),
            (0, 2, d),
            (0, 3, d),
            (1, 2, d),
            (1, 3, d),
            (2, 3, d),
        ];
        let result = sdr_embed(4, &pairs, &SdrConfig::default());
        assert_eq!(result.num_atoms, 4);

        // Check distances are approximately correct
        for &(i, j, d_target) in &pairs {
            let dx = result.coords[i * 3] - result.coords[j * 3];
            let dy = result.coords[i * 3 + 1] - result.coords[j * 3 + 1];
            let dz = result.coords[i * 3 + 2] - result.coords[j * 3 + 2];
            let d_actual = (dx * dx + dy * dy + dz * dz).sqrt();
            assert!(
                (d_actual - d_target).abs() < 0.5,
                "Distance ({},{}) should be ~{}: got {}",
                i,
                j,
                d_target,
                d_actual
            );
        }
    }
}
