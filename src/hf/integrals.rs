//! Electron repulsion integral (ERI) evaluator.
//!
//! Implements the Obara-Saika recurrence for two-electron integrals
//! over Cartesian Gaussian shells. For HF-3c, we need up to (pp|pp).

use super::basis::BasisSet;

/// Compute all electron-repulsion integrals of a basis set.
///
/// Returns the packed triangular array of the `(ij|kl)` integrals in
/// Hartree. Only the symmetry-distinct elements are stored: the array is
/// indexed by [`eri_index`], which exploits the permutation symmetry
/// `(ij|kl) = (ji|kl) = (ij|lk) = (kl|ij)`.
///
/// The integrals themselves are produced by the cross-validated `scf`
/// integral engine (see [`crate::scf::two_electron`]) and repacked into the
/// compact storage used here.
///
/// # Example
///
/// ```
/// use sci_form::hf::basis::build_sto3g_basis;
/// use sci_form::hf::integrals::{compute_eris, get_eri};
///
/// let basis = build_sto3g_basis(&[1, 1], &[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]]);
/// let eris = compute_eris(&basis);
/// // (11|11) in Hartree — positive, since it is a repulsion.
/// let eri = get_eri(&eris, 0, 0, 0, 0, basis.n_basis());
/// assert!(eri > 0.0);
/// ```
pub fn compute_eris(basis: &BasisSet) -> Vec<f64> {
    // Route through the cross-validated `scf` integral engine (F36/F37) and
    // repack the integrals into the compact storage used by `get_eri`.
    let scf_basis = super::basis::to_scf_basis(basis);
    let engine = crate::scf::two_electron::TwoElectronIntegrals::compute(&scf_basis);
    let n = basis.n_basis();
    let mut eris = vec![0.0f64; eri_storage_size(n)];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                for l in 0..n {
                    eris[eri_index(i, j, k, l, n)] = engine.get(i, j, k, l);
                }
            }
        }
    }
    eris
}

fn eri_storage_size(n: usize) -> usize {
    let nn = n * (n + 1) / 2;
    nn * (nn + 1) / 2
}

/// Flat index of the ERI `(ij|kl)` in the packed storage array.
///
/// The four AO indices are first folded into the pair indices
/// `ij = max(i,j)(max(i,j)+1)/2 + min(i,j)` and `kl` likewise, then the
/// pair index `max(ij,kl)(max(ij,kl)+1)/2 + min(ij,kl)` is formed. This
/// stores each of the `n(n+1)/2` unique electron pairs only once.
///
/// `n` (the number of basis functions) is not needed by the index
/// computation and is ignored; it is kept in the signature so that the
/// index can be computed without any extra context by the caller.
///
/// # Panics
///
/// This function never panics: it is pure index arithmetic.
pub fn eri_index(i: usize, j: usize, k: usize, l: usize, _n: usize) -> usize {
    let ij = if i >= j {
        i * (i + 1) / 2 + j
    } else {
        j * (j + 1) / 2 + i
    };
    let kl = if k >= l {
        k * (k + 1) / 2 + l
    } else {
        l * (l + 1) / 2 + k
    };
    if ij >= kl {
        ij * (ij + 1) / 2 + kl
    } else {
        kl * (kl + 1) / 2 + ij
    }
}

/// Read the ERI `(ij|kl)` in Hartree from the packed storage array.
///
/// `eris` must be the array produced by [`compute_eris`] for a basis with
/// `n` basis functions, and `i`, `j`, `k`, `l` must all be smaller than
/// `n`.
///
/// # Panics
///
/// Panics if `eris` is shorter than the packed array for `n` basis
/// functions, or if any of the AO indices is out of range for that basis.
pub fn get_eri(eris: &[f64], i: usize, j: usize, k: usize, l: usize, n: usize) -> f64 {
    eris[eri_index(i, j, k, l, n)]
}

#[cfg(test)]
mod tests {
    use super::super::basis::build_sto3g_basis;
    use super::*;

    #[test]
    fn test_eri_h2() {
        let basis = build_sto3g_basis(&[1, 1], &[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]]);
        let eris = compute_eris(&basis);
        let n = basis.n_basis();
        // (11|11) should be positive (electron-electron repulsion)
        let eri_0000 = get_eri(&eris, 0, 0, 0, 0, n);
        assert!(eri_0000 > 0.0, "(11|11) = {eri_0000}");
    }

    #[test]
    fn test_eri_symmetry() {
        let basis = build_sto3g_basis(&[1, 1], &[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]]);
        let eris = compute_eris(&basis);
        let n = basis.n_basis();
        // Permutational symmetry
        assert_eq!(get_eri(&eris, 0, 1, 0, 1, n), get_eri(&eris, 1, 0, 1, 0, n));
        assert_eq!(get_eri(&eris, 0, 1, 0, 1, n), get_eri(&eris, 0, 1, 1, 0, n));
    }
    #[test]
    fn cholesky_decomposition_reconstructs_a_dense_psd_matrix() {
        let dim = 6usize;
        let mut m = vec![0.0f64; dim * dim];
        for i in 0..dim {
            for j in 0..dim {
                m[i * dim + j] = ((i * 7 + j * 3) % 11) as f64 / 7.0;
            }
        }
        let mut a = vec![vec![0.0f64; dim]; dim];
        for i in 0..dim {
            for j in 0..dim {
                for k in 0..dim {
                    a[i][j] += m[i * dim + k] * m[j * dim + k];
                }
            }
        }
        let mut packed = vec![0.0f64; dim * (dim + 1) / 2];
        for i in 0..dim {
            for j in 0..=i {
                packed[i * (i + 1) / 2 + j] = a[i][j];
            }
        }
        let compressed = CholeskyEris::from_packed(&packed, 3, 0.0);
        let mut worst = 0.0f64;
        for i in 0..3 {
            for j in 0..=i {
                for k in 0..3 {
                    for l in 0..=k {
                        let p = pair_index(i, j);
                        let q = pair_index(k, l);
                        let exact = a[p][q];
                        let approx = compressed.get(i, j, k, l);
                        worst = worst.max((exact - approx).abs());
                    }
                }
            }
        }
        assert!(worst < 1e-10, "synthetic Cholesky error {worst}");
    }

    #[test]
    fn cholesky_compression_reproduces_sto3g_water_eris() {
        use super::super::basis::build_sto3g_basis;
        let positions = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
        let basis = build_sto3g_basis(&[8, 1, 1], &positions);
        let n = basis.n_basis();
        let full = compute_eris(&basis);
        let compressed = compute_eris_cholesky(&basis, 1e-10);
        assert!(compressed.rank() > 0);
        assert!(compressed.rank() <= compressed.n_pairs());
        let mut worst = 0.0f64;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    for l in 0..n {
                        let exact = get_eri(&full, i, j, k, l, n);
                        let approx = compressed.get(i, j, k, l);
                        worst = worst.max((exact - approx).abs());
                    }
                }
            }
        }
        println!(
            "Cholesky water: rank {}/{}, ratio {:.3}, max |ΔERI| {:.3e}",
            compressed.rank(),
            compressed.n_pairs(),
            compressed.compression_ratio(),
            worst
        );
        assert!(worst < 1e-7, "Cholesky ERI error too large: {worst:.3e}");
    }

    #[test]
    fn cholesky_compression_reproduces_benzene_eris_and_compresses() {
        use super::super::basis::build_sto3g_basis;
        // Planar benzene, STO-3G (36 basis functions, 666 pairs).
        let r_cc = 1.397;
        let mut positions = Vec::new();
        for k in 0..6 {
            let angle = std::f64::consts::PI * k as f64 / 3.0;
            positions.push([r_cc * angle.cos(), r_cc * angle.sin(), 0.0]);
        }
        for k in 0..6 {
            let angle = std::f64::consts::PI * k as f64 / 3.0;
            positions.push([
                (r_cc + 1.081) * angle.cos(),
                (r_cc + 1.081) * angle.sin(),
                0.0,
            ]);
        }
        let elements = [6u8, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1];
        let basis = build_sto3g_basis(&elements, &positions);
        let n = basis.n_basis();
        let full = compute_eris(&basis);
        let compressed = compute_eris_cholesky(&basis, 1e-6);
        let mut worst = 0.0f64;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    for l in 0..n {
                        let exact = get_eri(&full, i, j, k, l, n);
                        let approx = compressed.get(i, j, k, l);
                        worst = worst.max((exact - approx).abs());
                    }
                }
            }
        }
        println!(
            "Cholesky benzene: rank {}/{}, ratio {:.3}, max |ΔERI| {:.3e}",
            compressed.rank(),
            compressed.n_pairs(),
            compressed.compression_ratio(),
            worst
        );
        assert!(worst < 1e-5, "Cholesky ERI error too large: {worst:.3e}");
        assert!(
            compressed.compression_ratio() < 0.75,
            "expected meaningful compression, got ratio {}",
            compressed.compression_ratio()
        );
    }
}

// ─── Cholesky-compressed ERI storage ────────────────────────────────────────

/// Cholesky-compressed four-index electron-repulsion integrals.
///
/// The packed ERI array is a symmetric positive-semidefinite matrix in the
/// pair index `(ij)`; a pivoted Cholesky factorization gives
///
/// `(ij|kl) ≈ Σ_r L[(ij),r] L[(kl),r]`
///
/// with rank `r` chosen by a diagonal threshold.  This is the Cholesky
/// alternative to density fitting: no auxiliary basis is required and the
/// approximation error is controlled by `threshold`.
#[derive(Debug, Clone)]
pub struct CholeskyEris {
    n_basis: usize,
    rank: usize,
    factors: Vec<f64>,
    pivots: Vec<usize>,
    threshold: f64,
}

impl CholeskyEris {
    /// Number of basis functions.
    pub fn n_basis(&self) -> usize {
        self.n_basis
    }

    /// Retained Cholesky rank.
    pub fn rank(&self) -> usize {
        self.rank
    }

    /// Number of pair indices (`n_basis (n_basis + 1) / 2`).
    pub fn n_pairs(&self) -> usize {
        let n = self.n_basis;
        n * (n + 1) / 2
    }

    /// Threshold used to stop the factorization.
    pub fn threshold(&self) -> f64 {
        self.threshold
    }

    /// Approximate `(ij|kl)` in Hartree.
    pub fn get(&self, i: usize, j: usize, k: usize, l: usize) -> f64 {
        let p = pair_index(i, j);
        let q = pair_index(k, l);
        let rank = self.rank;
        let mut value = 0.0;
        for r in 0..rank {
            value += self.factors[p * rank + r] * self.factors[q * rank + r];
        }
        value
    }

    /// Bytes occupied by the compressed factors and pivot list.
    pub fn memory_bytes(&self) -> usize {
        self.factors.len() * std::mem::size_of::<f64>()
            + self.pivots.len() * std::mem::size_of::<usize>()
    }

    /// Compression ratio relative to the full packed ERI array.
    pub fn compression_ratio(&self) -> f64 {
        let npair = self.n_pairs();
        let full = npair * (npair + 1) / 2 * std::mem::size_of::<f64>();
        if full == 0 {
            return 1.0;
        }
        self.memory_bytes() as f64 / full as f64
    }

    /// Build the compressed representation from a packed ERI array.
    ///
    /// `threshold` is an absolute cutoff on the remaining pivot diagonal
    /// (Hartree²).  A threshold of `1e-10` typically retains micro-Hartree
    /// accuracy for STO-3G organic molecules.
    pub fn from_packed(eris: &[f64], n_basis: usize, threshold: f64) -> Self {
        let npair = n_basis * (n_basis + 1) / 2;
        assert!(
            eris.len() >= npair * (npair + 1) / 2,
            "packed ERI array is too short for {n_basis} basis functions"
        );
        assert!(threshold >= 0.0, "threshold must be non-negative");

        let mut pivots: Vec<usize> = (0..npair).collect();
        let mut diagonal = vec![0.0f64; npair];
        for (p, value) in diagonal.iter_mut().enumerate() {
            // A[(p,p)] in the packed representation.
            *value = eris[p * (p + 1) / 2 + p];
        }

        let max_rank = npair;
        let mut factors = vec![0.0f64; npair * max_rank];
        let mut rank = 0usize;

        for k in 0..npair {
            // Pivot on the largest remaining diagonal.
            let mut pivot = k;
            for i in (k + 1)..npair {
                if diagonal[i] > diagonal[pivot] {
                    pivot = i;
                }
            }
            if diagonal[pivot] <= threshold || diagonal[pivot] <= 0.0 {
                break;
            }
            if pivot != k {
                pivots.swap(k, pivot);
                diagonal.swap(k, pivot);
                for r in 0..k {
                    factors.swap(k * max_rank + r, pivot * max_rank + r);
                }
            }

            let pivot_value = diagonal[k].sqrt();
            factors[k * max_rank + k] = pivot_value;

            for i in (k + 1)..npair {
                let row_i = pivots[i];
                let row_k = pivots[k];
                let mut sum = eris[packed_pair_index(row_i, row_k)];
                for r in 0..k {
                    sum -= factors[i * max_rank + r] * factors[k * max_rank + r];
                }
                let value = sum / pivot_value;
                factors[i * max_rank + k] = value;
                diagonal[i] -= value * value;
            }
            rank = k + 1;
        }

        // Repack from the max-rank working stride to the compact rank stride.
        let mut compact = vec![0.0f64; npair * rank];
        for row in 0..npair {
            compact[row * rank..(row + 1) * rank]
                .copy_from_slice(&factors[row * max_rank..row * max_rank + rank]);
        }
        // Undo the pivot permutation so `get` can index factor rows by the
        // original pair index.
        let mut reordered = vec![0.0f64; npair * rank];
        for position in 0..npair {
            let pair = pivots[position];
            reordered[pair * rank..(pair + 1) * rank]
                .copy_from_slice(&compact[position * rank..(position + 1) * rank]);
        }
        Self {
            n_basis,
            rank,
            factors: reordered,
            pivots,
            threshold,
        }
    }
}

fn pair_index(i: usize, j: usize) -> usize {
    if i >= j {
        i * (i + 1) / 2 + j
    } else {
        j * (j + 1) / 2 + i
    }
}

fn packed_pair_index(p: usize, q: usize) -> usize {
    if p >= q {
        p * (p + 1) / 2 + q
    } else {
        q * (q + 1) / 2 + p
    }
}

/// Compute the Cholesky-compressed ERI tensor for a basis set.
pub fn compute_eris_cholesky(basis: &BasisSet, threshold: f64) -> CholeskyEris {
    let eris = compute_eris(basis);
    CholeskyEris::from_packed(&eris, basis.n_basis(), threshold)
}
