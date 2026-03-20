//! Mulliken population analysis.
//!
//! The Mulliken population assigns electron density to individual atoms:
//!
//!   q_A = Z_A - Σ_{μ ∈ A} (PS)_{μμ}
//!
//! where P is the density matrix and S is the overlap matrix.

use nalgebra::DMatrix;

/// Mulliken population analysis result.
#[derive(Debug, Clone)]
pub struct MullikenResult {
    /// Gross atomic populations (electrons per atom).
    pub populations: Vec<f64>,
    /// Partial charges (= Z - population).
    pub charges: Vec<f64>,
    /// Per-orbital populations.
    pub orbital_populations: Vec<f64>,
}

/// Compute Mulliken population analysis.
///
/// # Arguments
/// - `density` — Density matrix P (n_basis × n_basis)
/// - `overlap` — Overlap matrix S (n_basis × n_basis)
/// - `basis_to_atom` — Mapping from basis function index to atom index
/// - `atomic_numbers` — Atomic numbers Z for each atom
pub fn mulliken_analysis(
    density: &DMatrix<f64>,
    overlap: &DMatrix<f64>,
    basis_to_atom: &[usize],
    atomic_numbers: &[u8],
) -> MullikenResult {
    let n_basis = density.nrows();
    let n_atoms = atomic_numbers.len();

    // PS matrix product
    let ps = density * overlap;

    // Orbital populations: (PS)_μμ
    let orbital_populations: Vec<f64> = (0..n_basis).map(|i| ps[(i, i)]).collect();

    // Gross atomic populations: sum orbital populations per atom
    let mut populations = vec![0.0; n_atoms];
    for (mu, &atom) in basis_to_atom.iter().enumerate() {
        populations[atom] += orbital_populations[mu];
    }

    // Charges: q_A = Z_A - pop_A
    let charges: Vec<f64> = atomic_numbers
        .iter()
        .enumerate()
        .map(|(a, &z)| z as f64 - populations[a])
        .collect();

    MullikenResult {
        populations,
        charges,
        orbital_populations,
    }
}

/// Löwdin population analysis (orthogonalized basis).
///
/// Uses S^{1/2} P S^{1/2} instead of PS for more symmetric charge partitioning.
pub fn lowdin_analysis(
    density: &DMatrix<f64>,
    overlap: &DMatrix<f64>,
    basis_to_atom: &[usize],
    atomic_numbers: &[u8],
) -> MullikenResult {
    let n_basis = density.nrows();
    let n_atoms = atomic_numbers.len();

    // Compute S^{1/2}
    let eigen = overlap.clone().symmetric_eigen();
    let mut s_half = DMatrix::zeros(n_basis, n_basis);
    for k in 0..n_basis {
        let sqrt_eigenval = eigen.eigenvalues[k].max(0.0).sqrt();
        for i in 0..n_basis {
            for j in 0..n_basis {
                s_half[(i, j)] += sqrt_eigenval
                    * eigen.eigenvectors[(i, k)]
                    * eigen.eigenvectors[(j, k)];
            }
        }
    }

    // Löwdin density: P_L = S^{1/2} P S^{1/2}
    let p_lowdin = &s_half * density * &s_half;

    let orbital_populations: Vec<f64> = (0..n_basis).map(|i| p_lowdin[(i, i)]).collect();

    let mut populations = vec![0.0; n_atoms];
    for (mu, &atom) in basis_to_atom.iter().enumerate() {
        populations[atom] += orbital_populations[mu];
    }

    let charges: Vec<f64> = atomic_numbers
        .iter()
        .enumerate()
        .map(|(a, &z)| z as f64 - populations[a])
        .collect();

    MullikenResult {
        populations,
        charges,
        orbital_populations,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mulliken_charge_neutrality() {
        // For a neutral molecule, sum of charges should be ~0
        let n = 3;
        let p = DMatrix::from_row_slice(n, n, &[
            1.2, 0.1, 0.0,
            0.1, 1.0, 0.05,
            0.0, 0.05, 0.8,
        ]);
        let s = DMatrix::identity(n, n);
        let basis_to_atom = vec![0, 0, 1]; // 2 functions on atom 0, 1 on atom 1
        let atomic_numbers = vec![6, 1]; // C, H

        let result = mulliken_analysis(&p, &s, &basis_to_atom, &atomic_numbers);

        // Populations should sum to Tr(P) when S = I
        let trace_p: f64 = (0..n).map(|i| p[(i, i)]).sum();
        let pop_sum: f64 = result.populations.iter().sum();
        assert!((pop_sum - trace_p).abs() < 1e-12);
    }
}
