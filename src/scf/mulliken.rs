//! Mulliken and Löwdin population analysis.
//!
//! Mulliken: q_A = Z_A - Σ_{μ ∈ A} (PS)_{μμ}
//! Löwdin:   uses S^{1/2} P S^{1/2} for more symmetric partitioning.

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
pub fn mulliken_analysis(
    density: &DMatrix<f64>,
    overlap: &DMatrix<f64>,
    basis_to_atom: &[usize],
    atomic_numbers: &[u8],
) -> MullikenResult {
    let n_basis = density.nrows();
    let n_atoms = atomic_numbers.len();

    let ps = density * overlap;

    let orbital_populations: Vec<f64> = (0..n_basis).map(|i| ps[(i, i)]).collect();

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
                s_half[(i, j)] +=
                    sqrt_eigenval * eigen.eigenvectors[(i, k)] * eigen.eigenvectors[(j, k)];
            }
        }
    }

    // Löwdin density: S^{1/2} P S^{1/2}
    let lowdin_density = &s_half * density * &s_half;

    let orbital_populations: Vec<f64> = (0..n_basis).map(|i| lowdin_density[(i, i)]).collect();

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
    fn test_mulliken_identity_overlap() {
        // With S = I, Mulliken charges = Z - diagonal of P
        let p = DMatrix::from_row_slice(2, 2, &[1.5, 0.2, 0.2, 0.5]);
        let s = DMatrix::identity(2, 2);
        let basis_to_atom = vec![0, 1];
        let elements = [6u8, 1];

        let result = mulliken_analysis(&p, &s, &basis_to_atom, &elements);

        assert!((result.populations[0] - 1.5).abs() < 1e-14);
        assert!((result.populations[1] - 0.5).abs() < 1e-14);
        assert!((result.charges[0] - 4.5).abs() < 1e-14); // C: 6 - 1.5
        assert!((result.charges[1] - 0.5).abs() < 1e-14); // H: 1 - 0.5
    }

    #[test]
    fn test_mulliken_charge_conservation() {
        // Total electrons should be conserved (Tr(PS))
        let p = DMatrix::from_row_slice(2, 2, &[1.5, 0.3, 0.3, 0.5]);
        let s = DMatrix::from_row_slice(2, 2, &[1.0, 0.2, 0.2, 1.0]);
        let basis_to_atom = vec![0, 1];
        let elements = [1u8, 1];

        let result = mulliken_analysis(&p, &s, &basis_to_atom, &elements);

        let total_pop: f64 = result.populations.iter().sum();
        let ps = &p * &s;
        let trace_ps = ps[(0, 0)] + ps[(1, 1)];

        assert!((total_pop - trace_ps).abs() < 1e-12);
    }

    #[test]
    fn test_lowdin_total_population() {
        let p = DMatrix::from_row_slice(2, 2, &[1.0, 0.3, 0.3, 1.0]);
        let s = DMatrix::identity(2, 2);
        let basis_to_atom = vec![0, 1];
        let elements = [1u8, 1];

        let result = lowdin_analysis(&p, &s, &basis_to_atom, &elements);

        // With S = I, Löwdin = Mulliken
        assert!((result.populations[0] - 1.0).abs() < 1e-10);
        assert!((result.populations[1] - 1.0).abs() < 1e-10);
    }
}
