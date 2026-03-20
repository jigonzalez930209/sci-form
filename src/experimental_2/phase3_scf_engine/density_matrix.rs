//! Density matrix construction from MO coefficients.
//!
//! P_μν = 2 Σ_{k ∈ occ} C_μk C_νk
//!
//! The density matrix represents the electron density in the AO basis.
//! The factor of 2 accounts for double occupancy in RHF (closed-shell).

use nalgebra::DMatrix;

/// Build the density matrix from MO coefficients.
///
/// P_μν = 2 Σ_{k=0}^{n_occ-1} C_μk · C_νk
///
/// Only the first `n_occupied` columns of `coefficients` are used.
pub fn build_density_matrix(coefficients: &DMatrix<f64>, n_occupied: usize) -> DMatrix<f64> {
    let n = coefficients.nrows();
    let mut p = DMatrix::zeros(n, n);

    for i in 0..n {
        for j in 0..=i {
            let mut p_ij = 0.0;
            for k in 0..n_occupied {
                p_ij += coefficients[(i, k)] * coefficients[(j, k)];
            }
            p_ij *= 2.0;
            p[(i, j)] = p_ij;
            p[(j, i)] = p_ij;
        }
    }

    p
}

/// Compute the density matrix change (RMS difference).
pub fn density_rms_change(p_new: &DMatrix<f64>, p_old: &DMatrix<f64>) -> f64 {
    let n = p_new.nrows();
    let diff = p_new - p_old;
    let mut sum_sq = 0.0;
    for i in 0..n {
        for j in 0..n {
            sum_sq += diff[(i, j)] * diff[(i, j)];
        }
    }
    (sum_sq / (n * n) as f64).sqrt()
}

/// Compute the number of electrons from the density matrix.
///
/// N_e = Tr(PS) = Σ_μν P_μν S_μν
pub fn electron_count(p: &DMatrix<f64>, s: &DMatrix<f64>) -> f64 {
    let ps = p * s;
    let mut trace = 0.0;
    for i in 0..ps.nrows() {
        trace += ps[(i, i)];
    }
    trace
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_density_matrix_symmetric() {
        let c = DMatrix::from_row_slice(3, 3, &[
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
        ]);
        let p = build_density_matrix(&c, 1);

        for i in 0..3 {
            for j in 0..3 {
                assert!((p[(i, j)] - p[(j, i)]).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn test_density_trace() {
        // For n_occ occupied orbitals: Tr(P) = 2·n_occ (if S = I)
        let c = DMatrix::identity(3, 3);
        let p = build_density_matrix(&c, 2);
        let s = DMatrix::identity(3, 3);
        let n_e = electron_count(&p, &s);
        assert!((n_e - 4.0).abs() < 1e-10); // 2 × 2 = 4 electrons
    }

    #[test]
    fn test_density_change() {
        let p1 = DMatrix::identity(2, 2);
        let mut p2 = DMatrix::identity(2, 2);
        p2[(0, 0)] = 1.1;
        let rms = density_rms_change(&p1, &p2);
        assert!(rms > 0.0);
    }
}
