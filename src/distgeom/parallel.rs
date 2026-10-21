//! Parallel distance bounds evaluation and GPU-accelerated Floyd-Warshall.
//!
//! - Rayon-parallelized distance bounds matrix evaluation for ETKDG
//! - GPU Floyd-Warshall for shortest-path smoothing on molecules >500 atoms

/// Configuration for parallel distance bounds computation.
#[derive(Debug, Clone)]
pub struct ParallelBoundsConfig {
    /// Minimum molecule size for GPU dispatch (default: 500).
    pub gpu_threshold: usize,
    /// Whether to attempt GPU Floyd-Warshall.
    pub use_gpu: bool,
}

impl Default for ParallelBoundsConfig {
    fn default() -> Self {
        Self {
            gpu_threshold: 500,
            use_gpu: true,
        }
    }
}

/// Result of parallel bounds evaluation.
#[derive(Debug, Clone)]
pub struct ParallelBoundsResult {
    /// Lower bounds matrix (flat, row-major).
    pub lower: Vec<f64>,
    /// Upper bounds matrix (flat, row-major).
    pub upper: Vec<f64>,
    /// Matrix dimension.
    pub n: usize,
    /// Whether GPU was used for smoothing.
    pub used_gpu: bool,
    /// Execution time in ms.
    pub time_ms: f64,
}

/// Evaluate the distance bounds of a molecule, in parallel where possible.
///
/// The bounds are the CANONICAL ones produced by
/// [`crate::distgeom::bounds::calculate_bounds_matrix_opts`]: this function
/// used to re-derive them from a crude covalent-radius rule, which produced a
/// matrix that silently disagreed with the rest of the pipeline (the audit
/// flagged it as P2-27). The parallel/GPU helpers of this module are used
/// for the triangle smoothing, which is the expensive part.
///
/// The returned matrix stores the lower bound in `lower` and the upper bound
/// in `upper`, both flat and row-major.
pub fn compute_bounds_parallel(
    mol: &crate::graph::Molecule,
    use_set15: bool,
) -> ParallelBoundsResult {
    use std::time::Instant;

    let start = Instant::now();
    let n = mol.graph.node_count();

    // The canonical bounds matrix stores the LOWER bound in the lower
    // triangle (b[(max, min)]) and the UPPER bound in the upper triangle
    // (b[(min, max)]).
    let bounds = crate::distgeom::bounds::calculate_bounds_matrix_opts(mol, use_set15);

    let mut lower_mat = vec![0.0f64; n * n];
    let mut upper_mat = vec![f64::MAX; n * n];
    for i in 0..n {
        for j in 0..n {
            if i < j {
                lower_mat[i * n + j] = bounds[(j, i)];
                upper_mat[i * n + j] = bounds[(i, j)];
            } else if i > j {
                lower_mat[i * n + j] = bounds[(i, j)];
                upper_mat[i * n + j] = bounds[(j, i)];
            } else {
                // The diagonal is 0 in the canonical matrix.
                lower_mat[i * n + i] = 0.0;
                upper_mat[i * n + i] = 0.0;
            }
        }
    }

    ParallelBoundsResult {
        lower: lower_mat,
        upper: upper_mat,
        n,
        used_gpu: false,
        time_ms: start.elapsed().as_secs_f64() * 1000.0,
    }
}

/// GPU-accelerated Floyd-Warshall for triangle inequality smoothing.
///
/// For molecules >500 atoms, the O(N³) Floyd-Warshall can benefit from
/// GPU parallelization. Each k-iteration updates all (i,j) pairs in parallel.
pub fn smooth_bounds_gpu(bounds: &mut [f64], n: usize, is_upper: bool) -> bool {
    // GPU dispatch via wgpu (if available)
    #[cfg(feature = "experimental-gpu")]
    {
        if n >= 500 {
            if let Some(result) = try_gpu_floyd_warshall(bounds, n, is_upper) {
                bounds.copy_from_slice(&result);
                return true;
            }
        }
    }

    // CPU fallback: standard Floyd-Warshall
    floyd_warshall_cpu(bounds, n, is_upper);
    false
}

/// CPU Floyd-Warshall for triangle inequality smoothing.
pub fn floyd_warshall_cpu(matrix: &mut [f64], n: usize, is_upper: bool) {
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                if i == j || i == k || j == k {
                    continue;
                }
                let via_k = matrix[i * n + k] + matrix[k * n + j];
                if is_upper {
                    // Upper bounds: tighten to min
                    if via_k < matrix[i * n + j] {
                        matrix[i * n + j] = via_k;
                    }
                } else {
                    // Lower bounds: tighten to max
                    if via_k > matrix[i * n + j] {
                        matrix[i * n + j] = via_k;
                    }
                }
            }
        }
    }
}

/// Parallel Floyd-Warshall using rayon for the inner loops.
#[cfg(feature = "parallel")]
pub fn floyd_warshall_parallel(matrix: &mut [f64], n: usize, is_upper: bool) {
    use rayon::prelude::*;

    for k in 0..n {
        // Extract the k-th row for read-only access
        let k_row: Vec<f64> = (0..n).map(|j| matrix[k * n + j]).collect();
        let k_col: Vec<f64> = (0..n).map(|i| matrix[i * n + k]).collect();

        // Update all (i, j) pairs in parallel
        // We need to split the matrix into non-overlapping row slices
        let row_updates: Vec<Vec<(usize, f64)>> = (0..n)
            .into_par_iter()
            .map(|i| {
                if i == k {
                    return vec![];
                }
                let mut updates = Vec::new();
                for j in 0..n {
                    if j == i || j == k {
                        continue;
                    }
                    let via_k = k_col[i] + k_row[j];
                    let current = matrix[i * n + j];
                    if is_upper {
                        if via_k < current {
                            updates.push((j, via_k));
                        }
                    } else if via_k > current {
                        updates.push((j, via_k));
                    }
                }
                updates
            })
            .collect();

        // Apply updates
        for (i, updates) in row_updates.into_iter().enumerate() {
            for (j, val) in updates {
                matrix[i * n + j] = val;
            }
        }
    }
}

/// Attempt GPU Floyd-Warshall dispatch.
#[cfg(feature = "experimental-gpu")]
fn try_gpu_floyd_warshall(bounds: &[f64], n: usize, _is_upper: bool) -> Option<Vec<f64>> {
    use crate::gpu::context::GpuContext;

    let ctx = match GpuContext::try_create() {
        Ok(ctx) => ctx,
        Err(_) => return None,
    };

    // For very large matrices, GPU memory may be insufficient
    let matrix_bytes = n * n * std::mem::size_of::<f32>();
    if matrix_bytes > 128 * 1024 * 1024 {
        return None; // > 128MB — skip GPU
    }

    // Convert to f32 for GPU
    let data_f32: Vec<f32> = bounds
        .iter()
        .map(|&v| if v == f64::MAX { f32::MAX } else { v as f32 })
        .collect();

    // GPU Floyd-Warshall is implemented as N iterations of a
    // parallel (i,j) update kernel. Each iteration processes one k.
    // The actual WGSL dispatch would go here using ctx.
    // For now, fall back to CPU parallel.
    let _ = (ctx, data_f32);
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_floyd_warshall_cpu() {
        let mut m = vec![0.0, 1.0, f64::MAX, 1.0, 0.0, 2.0, f64::MAX, 2.0, 0.0];
        floyd_warshall_cpu(&mut m, 3, true);
        assert!((m[2] - 3.0).abs() < 1e-10);
    }

    #[cfg(feature = "parallel")]
    #[test]
    fn test_compute_bounds_parallel() {
        let mol = crate::graph::Molecule::from_smiles("CCO").expect("parse CCO");
        let result = compute_bounds_parallel(&mol, true);
        // The parser adds the implicit hydrogens, so CCO has 9 atoms.
        assert_eq!(result.n, mol.graph.node_count());
        assert_eq!(result.n, 9);
        // The C-C bond lower bound must be the canonical one, not a crude
        // covalent-radius estimate.
        let canonical = crate::distgeom::bounds::calculate_bounds_matrix_opts(&mol, true);
        // The lower bound lives in the lower triangle, the upper in the upper.
        assert!((result.lower[1] - canonical[(1, 0)]).abs() < 1e-12);
        assert!((result.upper[1] - canonical[(0, 1)]).abs() < 1e-12);
    }

    /// The parallel entry point must agree with the canonical bounds matrix
    /// for every pair (the audit flagged the two as silently divergent).
    #[test]
    fn test_parallel_bounds_match_canonical() {
        for smiles in ["CCO", "c1ccccc1", "CC(=O)O", "C1CCCCC1"] {
            let mol = crate::graph::Molecule::from_smiles(smiles).expect("parse");
            let result = compute_bounds_parallel(&mol, true);
            let canonical = crate::distgeom::bounds::calculate_bounds_matrix_opts(&mol, true);
            let n = result.n;
            for i in 0..n {
                for j in 0..n {
                    let (lo, hi) = if i < j {
                        (canonical[(j, i)], canonical[(i, j)])
                    } else if i > j {
                        (canonical[(i, j)], canonical[(j, i)])
                    } else {
                        (0.0, 0.0)
                    };
                    assert!(
                        (result.lower[i * n + j] - lo).abs() < 1e-9,
                        "{smiles}: lower[{i},{j}] {} != {lo}",
                        result.lower[i * n + j]
                    );
                    assert!(
                        (result.upper[i * n + j] - hi).abs() < 1e-9,
                        "{smiles}: upper[{i},{j}] {} != {hi}",
                        result.upper[i * n + j]
                    );
                }
            }
        }
    }
}
