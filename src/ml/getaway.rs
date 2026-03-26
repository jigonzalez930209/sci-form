//! GETAWAY (GEometry, Topology, and Atom-Weights AssemblY) descriptors.
//!
//! Combine 3D geometry (molecular influence matrix) with topological
//! information and atomic properties. The molecular influence matrix H
//! encodes how each atom influences the overall molecular shape.
//!
//! Reference: Consonni et al., J. Chem. Inf. Model. 42, 682–692 (2002).

use serde::{Deserialize, Serialize};

/// GETAWAY descriptor result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GetawayDescriptors {
    /// Leverage values (diagonal of influence matrix H_ii).
    pub leverages: Vec<f64>,
    /// H autocorrelation of lag k (unweighted): HATk.
    pub hat_autocorrelation: Vec<f64>,
    /// R autocorrelation (geometric distance weighted): Rk.
    pub r_autocorrelation: Vec<f64>,
    /// H total index: HATt = sum of |H_ij| for all bonded pairs.
    pub hat_total: f64,
    /// R total index: Rt.
    pub r_total: f64,
    /// Information content of leverages (Shannon entropy).
    pub h_information: f64,
    /// Maximum leverage.
    pub h_max: f64,
    /// Mean leverage.
    pub h_mean: f64,
    /// Maximum topological distance considered.
    pub max_lag: usize,
}

/// Compute GETAWAY descriptors from 3D coordinates and connectivity.
///
/// # Arguments
/// * `elements` - Atomic numbers
/// * `positions` - 3D coordinates
/// * `bonds` - Bond list as (atom_i, atom_j) pairs
/// * `max_lag` - Maximum topological distance for autocorrelation (default: 8)
pub fn compute_getaway(
    elements: &[u8],
    positions: &[[f64; 3]],
    bonds: &[(usize, usize)],
    max_lag: usize,
) -> GetawayDescriptors {
    let n = elements.len().min(positions.len());
    if n < 2 {
        return GetawayDescriptors {
            leverages: vec![],
            hat_autocorrelation: vec![0.0; max_lag],
            r_autocorrelation: vec![0.0; max_lag],
            hat_total: 0.0,
            r_total: 0.0,
            h_information: 0.0,
            h_max: 0.0,
            h_mean: 0.0,
            max_lag,
        };
    }

    // Compute centered coordinate matrix X (n × 3)
    let mut centroid = [0.0f64; 3];
    for i in 0..n {
        for d in 0..3 {
            centroid[d] += positions[i][d];
        }
    }
    for d in 0..3 {
        centroid[d] /= n as f64;
    }

    // X^T X (3×3 Gram matrix)
    let mut xtx = [[0.0f64; 3]; 3];
    for i in 0..n {
        let dx = [
            positions[i][0] - centroid[0],
            positions[i][1] - centroid[1],
            positions[i][2] - centroid[2],
        ];
        for r in 0..3 {
            for c in 0..3 {
                xtx[r][c] += dx[r] * dx[c];
            }
        }
    }

    // Influence matrix: H = X (X^T X)^{-1} X^T
    // Leverages: h_ii = X_i (X^T X)^{-1} X_i^T
    let xtx_inv = invert_3x3(&xtx);
    let mut leverages = Vec::with_capacity(n);
    for i in 0..n {
        let dx = [
            positions[i][0] - centroid[0],
            positions[i][1] - centroid[1],
            positions[i][2] - centroid[2],
        ];
        // h_ii = dx^T * xtx_inv * dx
        let mut h = 0.0;
        for r in 0..3 {
            for c in 0..3 {
                h += dx[r] * xtx_inv[r][c] * dx[c];
            }
        }
        leverages.push(h);
    }

    // Build topological distance matrix via BFS
    let topo_dist = build_topo_distance(n, bonds);

    // Compute H_ij for off-diagonal elements
    // H_ij = X_i (X^T X)^{-1} X_j^T / sqrt(h_ii * h_jj)
    let mut h_matrix = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        h_matrix[i][i] = leverages[i];
        let dxi = [
            positions[i][0] - centroid[0],
            positions[i][1] - centroid[1],
            positions[i][2] - centroid[2],
        ];
        for j in (i + 1)..n {
            let dxj = [
                positions[j][0] - centroid[0],
                positions[j][1] - centroid[1],
                positions[j][2] - centroid[2],
            ];
            let mut hij = 0.0;
            for r in 0..3 {
                for c in 0..3 {
                    hij += dxi[r] * xtx_inv[r][c] * dxj[c];
                }
            }
            h_matrix[i][j] = hij;
            h_matrix[j][i] = hij;
        }
    }

    // Geometric distance matrix
    let mut geo_dist = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            geo_dist[i][j] = r;
            geo_dist[j][i] = r;
        }
    }

    // Autocorrelation at lag k
    let mut hat_autocorrelation = vec![0.0f64; max_lag];
    let mut r_autocorrelation = vec![0.0f64; max_lag];
    let mut hat_total = 0.0;
    let mut r_total = 0.0;

    for i in 0..n {
        for j in (i + 1)..n {
            let d = topo_dist[i][j];
            if d == 0 || d > max_lag {
                continue;
            }
            let k = d - 1;
            let hi = leverages[i].max(0.0).sqrt();
            let hj = leverages[j].max(0.0).sqrt();
            hat_autocorrelation[k] += hi * hj;

            let rij = geo_dist[i][j];
            if rij > 1e-12 {
                r_autocorrelation[k] += hi * hj / rij;
            }

            hat_total += h_matrix[i][j].abs();
            if rij > 1e-12 {
                r_total += h_matrix[i][j].abs() / rij;
            }
        }
    }

    // Information content of leverages
    let h_sum: f64 = leverages.iter().sum();
    let h_information = if h_sum > 1e-12 {
        let mut entropy = 0.0;
        for &h in &leverages {
            let p = h / h_sum;
            if p > 1e-12 {
                entropy -= p * p.ln();
            }
        }
        entropy
    } else {
        0.0
    };

    let h_max = leverages.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let h_mean = if n > 0 { h_sum / n as f64 } else { 0.0 };

    GetawayDescriptors {
        leverages,
        hat_autocorrelation,
        r_autocorrelation,
        hat_total,
        r_total,
        h_information,
        h_max,
        h_mean,
        max_lag,
    }
}

/// Build topological distance matrix using BFS.
fn build_topo_distance(n: usize, bonds: &[(usize, usize)]) -> Vec<Vec<usize>> {
    let mut adj = vec![vec![]; n];
    for &(a, b) in bonds {
        if a < n && b < n {
            adj[a].push(b);
            adj[b].push(a);
        }
    }

    let mut dist = vec![vec![0usize; n]; n];
    for start in 0..n {
        let mut visited = vec![false; n];
        visited[start] = true;
        let mut queue = std::collections::VecDeque::new();
        queue.push_back((start, 0usize));
        while let Some((node, d)) = queue.pop_front() {
            dist[start][node] = d;
            for &nb in &adj[node] {
                if !visited[nb] {
                    visited[nb] = true;
                    queue.push_back((nb, d + 1));
                }
            }
        }
    }

    dist
}

/// Invert a 3×3 matrix using Cramer's rule.
fn invert_3x3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

    if det.abs() < 1e-12 {
        return [[0.0; 3]; 3];
    }

    let inv_det = 1.0 / det;
    let mut inv = [[0.0f64; 3]; 3];

    inv[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det;
    inv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det;
    inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det;
    inv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det;
    inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det;
    inv[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det;
    inv[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det;
    inv[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det;
    inv[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det;

    inv
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_getaway_ethane() {
        let elements = vec![6, 6, 1, 1, 1, 1, 1, 1];
        let positions = vec![
            [0.0, 0.0, 0.0],   // C
            [1.54, 0.0, 0.0],  // C
            [-0.5, 0.9, 0.0],  // H
            [-0.5, -0.9, 0.0], // H
            [-0.5, 0.0, 0.9],  // H
            [2.04, 0.9, 0.0],  // H
            [2.04, -0.9, 0.0], // H
            [2.04, 0.0, 0.9],  // H
        ];
        let bonds = vec![(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (1, 6), (1, 7)];

        let g = compute_getaway(&elements, &positions, &bonds, 8);
        assert_eq!(g.leverages.len(), 8);
        assert!(g.hat_total > 0.0);
        assert!(g.h_max > 0.0);
    }

    #[test]
    fn test_getaway_empty() {
        let g = compute_getaway(&[], &[], &[], 8);
        assert!(g.leverages.is_empty());
        assert_eq!(g.hat_total, 0.0);
    }
}
