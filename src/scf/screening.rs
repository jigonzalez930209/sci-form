//! **ALPHA** — Schwarz pre-screening for electron repulsion integrals.
//!
//! Implements the Schwarz inequality screening to skip negligible ERI shells:
//!   |(μν|λσ)| ≤ Q_μν · Q_λσ
//!
//! where Q_μν = sqrt(|(μν|μν)|).

use super::basis::BasisSet;
use super::obara_saika::ShellPairData;

/// Pre-computed Schwarz screening bounds for all basis function pairs.
pub struct SchwarzScreener {
    /// Q_ij values for each pair (i, j). Stored as upper triangle, row-major.
    q_values: Vec<f64>,
    /// Number of basis functions.
    n_basis: usize,
    /// Screening threshold.
    pub threshold: f64,
}

impl SchwarzScreener {
    /// Build Schwarz bounds from a basis set.
    pub fn build(basis: &BasisSet, threshold: f64) -> Self {
        let n = basis.functions.len();
        let mut q_values = vec![0.0; n * n];

        for i in 0..n {
            for j in i..n {
                // Build shell pair for (i,j)|(i,j)
                let fi = &basis.functions[i];
                let fj = &basis.functions[j];

                // Sum over primitives for contracted basis functions
                let mut q_ij_sq = 0.0;
                for pi in &fi.primitives {
                    for pj in &fj.primitives {
                        let sp = ShellPairData::new(pi.alpha, fi.center, pj.alpha, fj.center);
                        let eri = super::obara_saika::eri_ssss(&sp, &sp);
                        q_ij_sq += pi.coefficient
                            * pj.coefficient
                            * pi.coefficient
                            * pj.coefficient
                            * eri.abs();
                    }
                }

                let q = q_ij_sq.sqrt();
                q_values[i * n + j] = q;
                q_values[j * n + i] = q;
            }
        }

        Self {
            q_values,
            n_basis: n,
            threshold,
        }
    }

    /// Check whether (ij|kl) integral should be computed (passes screening).
    #[inline]
    pub fn is_significant(&self, i: usize, j: usize, k: usize, l: usize) -> bool {
        self.q_values[i * self.n_basis + j] * self.q_values[k * self.n_basis + l] >= self.threshold
    }

    /// Get Q_ij bound value.
    #[inline]
    pub fn q(&self, i: usize, j: usize) -> f64 {
        self.q_values[i * self.n_basis + j]
    }

    /// Number of significant pairs above threshold.
    pub fn n_significant_pairs(&self) -> usize {
        let n = self.n_basis;
        let mut count = 0;
        for i in 0..n {
            for j in i..n {
                if self.q_values[i * n + j] >= self.threshold.sqrt() {
                    count += 1;
                }
            }
        }
        count
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scf::basis::BasisSet;
    use crate::scf::types::MolecularSystem;

    fn h2_basis() -> BasisSet {
        let system =
            MolecularSystem::from_angstrom(&[1, 1], &[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]], 0, 1);
        BasisSet::sto3g(&system.atomic_numbers, &system.positions_bohr)
    }

    #[test]
    fn screener_builds_for_h2() {
        let basis = h2_basis();
        let screener = SchwarzScreener::build(&basis, 1e-10);
        assert_eq!(screener.n_basis, basis.functions.len());
    }

    #[test]
    fn schwarz_values_are_positive() {
        let basis = h2_basis();
        let screener = SchwarzScreener::build(&basis, 1e-10);
        for i in 0..screener.n_basis {
            assert!(
                screener.q(i, i) >= 0.0,
                "diagonal Schwarz bound should be non-negative"
            );
        }
    }

    #[test]
    fn is_significant_for_diagonal() {
        let basis = h2_basis();
        let screener = SchwarzScreener::build(&basis, 1e-10);
        // (ii|ii) should always be significant
        for i in 0..screener.n_basis {
            assert!(screener.is_significant(i, i, i, i));
        }
    }

    #[test]
    fn n_significant_pairs_nonzero() {
        let basis = h2_basis();
        let screener = SchwarzScreener::build(&basis, 1e-10);
        assert!(screener.n_significant_pairs() > 0);
    }
}
