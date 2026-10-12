//! Butina (Taylor-Butina) clustering for conformer ensembles.
//!
//! Groups conformers by RMSD similarity using greedy distance-based clustering.

use crate::alignment::kabsch;
use serde::{Deserialize, Serialize};

/// Result of Butina clustering on a conformer ensemble.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClusterResult {
    /// Number of clusters found.
    pub n_clusters: usize,
    /// Cluster assignment for each conformer (0-indexed cluster ID).
    pub assignments: Vec<usize>,
    /// Indices of cluster centroids (representative conformers).
    pub centroid_indices: Vec<usize>,
    /// Number of members in each cluster.
    pub cluster_sizes: Vec<usize>,
    /// RMSD cutoff used.
    pub rmsd_cutoff: f64,
}

/// Compute the all-pairs RMSD matrix for a set of conformers.
///
/// `conformers`: slice of flat coordinate arrays `[x0,y0,z0, x1,y1,z1, ...]`.
/// All conformers must have the same number of atoms.
pub fn compute_rmsd_matrix(conformers: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = conformers.len();
    let mut matrix = vec![vec![0.0f64; n]; n];

    for i in 0..n {
        for j in (i + 1)..n {
            let rmsd = kabsch::compute_rmsd(&conformers[i], &conformers[j]);
            matrix[i][j] = rmsd;
            matrix[j][i] = rmsd;
        }
    }

    matrix
}

/// Parallel all-pairs RMSD matrix using rayon.
///
/// The O(N²) pairs are computed in parallel since each Kabsch RMSD
/// is independent. Results are collected and assembled into the
/// symmetric matrix.
#[cfg(feature = "parallel")]
pub fn compute_rmsd_matrix_parallel(conformers: &[Vec<f64>]) -> Vec<Vec<f64>> {
    use rayon::prelude::*;

    let n = conformers.len();
    let pairs: Vec<(usize, usize)> = (0..n)
        .flat_map(|i| ((i + 1)..n).map(move |j| (i, j)))
        .collect();

    let results: Vec<(usize, usize, f64)> = pairs
        .into_par_iter()
        .map(|(i, j)| {
            let rmsd = kabsch::compute_rmsd(&conformers[i], &conformers[j]);
            (i, j, rmsd)
        })
        .collect();

    let mut matrix = vec![vec![0.0f64; n]; n];
    for (i, j, rmsd) in results {
        matrix[i][j] = rmsd;
        matrix[j][i] = rmsd;
    }

    matrix
}

/// Perform Butina (Taylor-Butina) clustering on a set of conformers.
///
/// Algorithm:
/// 1. Compute the all-pairs RMSD matrix after Kabsch alignment
/// 2. For each conformer, count neighbors within `rmsd_cutoff`
/// 3. Select the conformer with the most neighbors as a cluster centroid
/// 4. Remove all its neighbors from the pool; repeat until empty
///
/// # Arguments
/// - `conformers`: slice of flat coordinate arrays
/// - `rmsd_cutoff`: RMSD threshold for clustering (Å), typically 1.0
///
/// # Returns
/// `ClusterResult` with cluster assignments, centroids, and sizes.
pub fn butina_cluster(conformers: &[Vec<f64>], rmsd_cutoff: f64) -> ClusterResult {
    let n = conformers.len();

    if n == 0 {
        return ClusterResult {
            n_clusters: 0,
            assignments: vec![],
            centroid_indices: vec![],
            cluster_sizes: vec![],
            rmsd_cutoff,
        };
    }

    if n == 1 {
        return ClusterResult {
            n_clusters: 1,
            assignments: vec![0],
            centroid_indices: vec![0],
            cluster_sizes: vec![1],
            rmsd_cutoff,
        };
    }

    // 1. Compute RMSD matrix (auto-dispatch to parallel when available)
    #[cfg(feature = "parallel")]
    let rmsd_matrix = compute_rmsd_matrix_parallel(conformers);
    #[cfg(not(feature = "parallel"))]
    let rmsd_matrix = compute_rmsd_matrix(conformers);

    // 2. Neighbor lists within the cutoff (self included, like RDKit)
    let neighbor_lists: Vec<Vec<usize>> = (0..n)
        .map(|i| {
            (0..n)
                .filter(|&j| rmsd_matrix[i][j] <= rmsd_cutoff)
                .collect()
        })
        .collect();

    // 3. Candidate seeds ordered by (neighbor count, index) descending.
    let mut sorted_indices: Vec<(usize, usize)> =
        (0..n).map(|i| (neighbor_lists[i].len(), i)).collect();
    sorted_indices.sort_by(|a, b| b.cmp(a));

    let mut seen = vec![false; n];
    let mut assignments = vec![usize::MAX; n];
    let mut centroid_indices = Vec::new();
    let mut cluster_sizes = Vec::new();

    // 4. Reordering variant of Taylor-Butina: seeds with at least one real
    //    neighbor are processed in descending neighbor-count order; after each
    //    cluster the counts of unassigned affected points are refreshed.
    while sorted_indices
        .first()
        .map(|&(count, _)| count > 1)
        .unwrap_or(false)
    {
        let (_, idx) = sorted_indices.remove(0);
        if seen[idx] {
            continue;
        }
        let cluster_id = centroid_indices.len();
        centroid_indices.push(idx);
        seen[idx] = true;
        assignments[idx] = cluster_id;
        let mut size = 1;

        for &nb in &neighbor_lists[idx] {
            if !seen[nb] {
                seen[nb] = true;
                assignments[nb] = cluster_id;
                size += 1;
            }
        }
        cluster_sizes.push(size);

        // Refresh counts for unassigned points adjacent to this cluster.
        let members: Vec<usize> = (0..n).filter(|&i| assignments[i] == cluster_id).collect();
        let mut affected: Vec<usize> = Vec::new();
        for &m in &members {
            for &p in &neighbor_lists[m] {
                if !seen[p] && !affected.contains(&p) {
                    affected.push(p);
                }
            }
        }
        if !affected.is_empty() {
            for entry in sorted_indices.iter_mut() {
                if affected.contains(&entry.1) {
                    entry.0 = neighbor_lists[entry.1]
                        .iter()
                        .filter(|&&p| !seen[p])
                        .count();
                }
            }
            sorted_indices.sort_by(|a, b| b.cmp(a));
        }
    }

    // 5. Remaining points become singleton clusters.
    while let Some((_, idx)) = sorted_indices.first().copied() {
        sorted_indices.remove(0);
        if seen[idx] {
            continue;
        }
        seen[idx] = true;
        let cluster_id = centroid_indices.len();
        centroid_indices.push(idx);
        assignments[idx] = cluster_id;
        cluster_sizes.push(1);
    }

    ClusterResult {
        n_clusters: centroid_indices.len(),
        assignments,
        centroid_indices,
        cluster_sizes,
        rmsd_cutoff,
    }
}

/// Filter a conformer ensemble to keep only cluster centroids.
///
/// Returns the indices of representative conformers.
pub fn filter_diverse_conformers(conformers: &[Vec<f64>], rmsd_cutoff: f64) -> Vec<usize> {
    let result = butina_cluster(conformers, rmsd_cutoff);
    result.centroid_indices
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_conformer() {
        let conformers = vec![vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0]];
        let result = butina_cluster(&conformers, 1.0);
        assert_eq!(result.n_clusters, 1);
        assert_eq!(result.assignments, vec![0]);
    }

    #[test]
    fn test_identical_conformers() {
        let coords = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let conformers = vec![coords.clone(), coords.clone(), coords.clone()];
        let result = butina_cluster(&conformers, 0.5);
        // All identical → 1 cluster
        assert_eq!(result.n_clusters, 1);
        assert_eq!(result.cluster_sizes, vec![3]);
    }

    #[test]
    fn test_distinct_conformers() {
        // Two very different conformers (large RMSD)
        let c1 = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let c2 = vec![0.0, 0.0, 0.0, 100.0, 0.0, 0.0];
        let conformers = vec![c1, c2];
        let result = butina_cluster(&conformers, 0.5);
        assert_eq!(result.n_clusters, 2);
    }

    #[test]
    fn test_empty_input() {
        let conformers: Vec<Vec<f64>> = vec![];
        let result = butina_cluster(&conformers, 1.0);
        assert_eq!(result.n_clusters, 0);
    }

    #[test]
    fn test_rmsd_matrix_symmetry() {
        let c1 = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let c2 = vec![0.0, 0.0, 0.0, 1.5, 0.0, 0.0];
        let c3 = vec![0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        let conformers = vec![c1, c2, c3];
        let matrix = compute_rmsd_matrix(&conformers);

        for i in 0..3 {
            assert!((matrix[i][i]).abs() < 1e-10, "diagonal must be 0");
            for j in 0..3 {
                assert!(
                    (matrix[i][j] - matrix[j][i]).abs() < 1e-10,
                    "matrix must be symmetric"
                );
            }
        }
    }
}
