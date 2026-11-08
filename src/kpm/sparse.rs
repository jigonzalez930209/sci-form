//! Sparse CSR matrices for the KPM backend.
//!
//! The dense KPM path materialises an N×N matrix (and, for the exact trace,
//! N vectors). This module adds a minimal CSR representation and sparse
//! Chebyshev recursion so a 1D/2D Hamiltonian with O(N) non-zeros can be
//! expanded without ever allocating the dense matrix.

use super::density::{KpmConfig, KpmDosResult};
use super::{jackson_kernel, ChebyshevExpansion};
use nalgebra::DMatrix;

/// Compressed sparse row matrix.
#[derive(Debug, Clone)]
pub struct SparseMatrix {
    /// Matrix dimension (square).
    pub n: usize,
    /// Row offsets, length `n + 1`.
    pub row_offsets: Vec<usize>,
    /// Column indices for each stored value.
    pub col_indices: Vec<usize>,
    /// Non-zero values (explicit zeros are dropped).
    pub values: Vec<f64>,
}

impl SparseMatrix {
    /// Build a CSR matrix from a dense matrix, dropping entries with
    /// `|value| <= 1e-15`.
    pub fn from_dense(matrix: &DMatrix<f64>) -> Self {
        let n = matrix.nrows();
        let mut row_offsets = Vec::with_capacity(n + 1);
        let mut col_indices = Vec::new();
        let mut values = Vec::new();
        row_offsets.push(0);
        for i in 0..n {
            for j in 0..n {
                let value = matrix[(i, j)];
                if value.abs() > 1e-15 {
                    col_indices.push(j);
                    values.push(value);
                }
            }
            row_offsets.push(col_indices.len());
        }
        Self {
            n,
            row_offsets,
            col_indices,
            values,
        }
    }

    /// Build a CSR matrix from `(row, column, value)` triplets.
    ///
    /// Duplicate entries are summed; the result is not required to be
    /// symmetric, but the KPM caller must pass a Hermitian operator.
    pub fn from_triplets(n: usize, triplets: &[(usize, usize, f64)]) -> Self {
        let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        for &(i, j, value) in triplets {
            if i < n && j < n && value.abs() > 1e-15 {
                rows[i].push((j, value));
            }
        }
        let mut row_offsets = Vec::with_capacity(n + 1);
        let mut col_indices = Vec::new();
        let mut values = Vec::new();
        row_offsets.push(0);
        for row in &mut rows {
            row.sort_by_key(|&(j, _)| j);
            let mut j = 0;
            while j < row.len() {
                let col = row[j].0;
                let mut sum = 0.0;
                while j < row.len() && row[j].0 == col {
                    sum += row[j].1;
                    j += 1;
                }
                if sum.abs() > 1e-15 {
                    col_indices.push(col);
                    values.push(sum);
                }
            }
            row_offsets.push(col_indices.len());
        }
        Self {
            n,
            row_offsets,
            col_indices,
            values,
        }
    }

    /// Number of stored non-zero entries.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// Matrix dimension.
    pub fn nrows(&self) -> usize {
        self.n
    }

    /// Sparse matrix-vector product `y = A x`.
    pub fn matvec(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.n);
        assert_eq!(y.len(), self.n);
        for (i, slot) in y.iter_mut().enumerate() {
            let mut sum = 0.0;
            for k in self.row_offsets[i]..self.row_offsets[i + 1] {
                sum += self.values[k] * x[self.col_indices[k]];
            }
            *slot = sum;
        }
    }

    /// Gershgorin spectral bounds with a 1 % margin.
    pub fn gershgorin_bounds(&self) -> (f64, f64) {
        let mut e_min = f64::INFINITY;
        let mut e_max = f64::NEG_INFINITY;
        for i in 0..self.n {
            let mut diagonal = 0.0;
            let mut radius = 0.0;
            for k in self.row_offsets[i]..self.row_offsets[i + 1] {
                if self.col_indices[k] == i {
                    diagonal += self.values[k];
                } else {
                    radius += self.values[k].abs();
                }
            }
            e_min = e_min.min(diagonal - radius);
            e_max = e_max.max(diagonal + radius);
        }
        let margin = (e_max - e_min).max(1e-12) * 0.01;
        (e_min - margin, e_max + margin)
    }
}

impl ChebyshevExpansion {
    /// Sparse Chebyshev expansion.
    ///
    /// `n_vectors = 0` computes the exact diagonal trace by applying the
    /// recurrence to every unit vector (cost `O(n · nnz · order)`); a positive
    /// value uses the same deterministic ±1/√n stochastic trace as the dense
    /// path.
    pub fn from_sparse(h: &SparseMatrix, order: usize, n_vectors: usize, seed: u64) -> Self {
        let n = h.n;
        let (e_min, e_max) = h.gershgorin_bounds();
        let a = (e_max - e_min) / 2.0;
        let b = (e_max + e_min) / 2.0;

        let n_vectors = if n_vectors == 0 { n } else { n_vectors };
        let sqrt_n = (n as f64).sqrt();
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        use rand::{Rng, SeedableRng};
        let vectors: Vec<Vec<f64>> = if n_vectors == n {
            (0..n)
                .map(|i| {
                    let mut v = vec![0.0; n];
                    v[i] = 1.0;
                    v
                })
                .collect()
        } else {
            (0..n_vectors)
                .map(|_| {
                    (0..n)
                        .map(|_| {
                            if rng.gen_bool(0.5) {
                                1.0 / sqrt_n
                            } else {
                                -1.0 / sqrt_n
                            }
                        })
                        .collect()
                })
                .collect()
        };

        let mut moments = vec![0.0; order];
        for r in &vectors {
            let local = sparse_moments_for_vector(r, h, order, a, b);
            for (sum, value) in moments.iter_mut().zip(local) {
                *sum += value;
            }
        }
        for moment in &mut moments {
            *moment /= n_vectors as f64;
        }

        Self {
            moments,
            a,
            b,
            n,
            order,
        }
    }
}

fn sparse_moments_for_vector(
    r: &[f64],
    h: &SparseMatrix,
    order: usize,
    a: f64,
    b: f64,
) -> Vec<f64> {
    let n = h.n;
    let mut moments = vec![0.0; order];
    let mut t_prev = r.to_vec();
    let mut t_curr = vec![0.0; n];
    sparse_h_tilde(h, r, a, b, &mut t_curr);
    moments[0] = dot(r, &t_prev);
    if order > 1 {
        moments[1] = dot(r, &t_curr);
    }
    let mut t_next = vec![0.0; n];
    let mut ht = vec![0.0; n];
    for k in 2..order {
        h.matvec(&t_curr, &mut ht);
        for i in 0..n {
            let scaled = (ht[i] - b * t_curr[i]) / a;
            t_next[i] = 2.0 * scaled - t_prev[i];
        }
        moments[k] = dot(r, &t_next);
        std::mem::swap(&mut t_prev, &mut t_curr);
        std::mem::swap(&mut t_curr, &mut t_next);
    }
    moments
}

fn sparse_h_tilde(h: &SparseMatrix, x: &[f64], a: f64, b: f64, out: &mut [f64]) {
    let mut hx = vec![0.0; h.n];
    h.matvec(x, &mut hx);
    for i in 0..h.n {
        out[i] = (hx[i] - b * x[i]) / a;
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Sparse KPM density of states.
///
/// Same output contract as [`super::compute_kpm_dos`], but the Chebyshev
/// recurrence runs on a [`SparseMatrix`] and never allocates the dense matrix.
pub fn compute_kpm_dos_sparse(
    h: &SparseMatrix,
    config: &KpmConfig,
    e_min_out: f64,
    e_max_out: f64,
    n_points: usize,
) -> KpmDosResult {
    let expansion = ChebyshevExpansion::from_sparse(h, config.order, config.n_vectors, config.seed);
    let gk = jackson_kernel(config.order);
    let step = (e_max_out - e_min_out) / (n_points - 1).max(1) as f64;
    let energies: Vec<f64> = (0..n_points).map(|i| e_min_out + i as f64 * step).collect();
    let total_dos = energies
        .iter()
        .map(|&e| expansion.dos_at_energy(e, &gk))
        .collect();
    KpmDosResult {
        energies,
        total_dos,
        e_min: expansion.b - expansion.a,
        e_max: expansion.b + expansion.a,
        order: config.order,
    }
}
