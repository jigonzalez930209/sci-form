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

    // 1. Compute RMSD matrix
    let rmsd_matrix = compute_rmsd_matrix(conformers);

    // 2. Build neighbor lists
    let mut neighbor_counts: Vec<(usize, usize)> = (0..n)
        .map(|i| {
            let count = (0..n)
                .filter(|&j| j != i && rmsd_matrix[i][j] <= rmsd_cutoff)
                .count();
            (i, count)
        })
        .collect();

    // 3. Greedy clustering
    let mut assignments = vec![usize::MAX; n];
    let mut centroid_indices = Vec::new();
    let mut cluster_sizes = Vec::new();
    let mut assigned = vec![false; n];

    loop {
        // Sort by neighbor count (descending)
        neighbor_counts.sort_by(|a, b| b.1.cmp(&a.1));

        // Find the unassigned conformer with the most neighbors
        let centroid = neighbor_counts
            .iter()
            .find(|(idx, _)| !assigned[*idx])
            .map(|(idx, _)| *idx);

        let centroid = match centroid {
            Some(c) => c,
            None => break,
        };

        let cluster_id = centroid_indices.len();
        centroid_indices.push(centroid);
        assigned[centroid] = true;
        assignments[centroid] = cluster_id;
        let mut size = 1;

        // Assign all unassigned neighbors to this cluster
        for j in 0..n {
            if !assigned[j] && rmsd_matrix[centroid][j] <= rmsd_cutoff {
                assigned[j] = true;
                assignments[j] = cluster_id;
                size += 1;
            }
        }

        cluster_sizes.push(size);

        // Update neighbor counts for remaining unassigned
        for entry in &mut neighbor_counts {
            if assigned[entry.0] {
                entry.1 = 0;
            } else {
                entry.1 = (0..n)
                    .filter(|&j| {
                        !assigned[j] && j != entry.0 && rmsd_matrix[entry.0][j] <= rmsd_cutoff
                    })
                    .count();
            }
        }

        if assigned.iter().all(|&a| a) {
            break;
        }
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
