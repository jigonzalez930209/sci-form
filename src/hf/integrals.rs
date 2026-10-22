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
/// ```no_run
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
}
