//! Generalized eigenproblem solver for EHT: HC = SCE.
//!
//! Uses Löwdin orthogonalization:
//! 1. Diagonalize S → eigenvalues λ, eigenvectors U
//! 2. Build S^{-1/2} = U diag(1/√λ) U^T
//! 3. Transform H' = S^{-1/2} H S^{-1/2}
//! 4. Diagonalize H' → eigenvalues E, eigenvectors C'
//! 5. Back-transform C = S^{-1/2} C'

use nalgebra::{DMatrix, DVector, SymmetricEigen};
use serde::{Deserialize, Serialize};

/// Result of an EHT calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EhtResult {
    /// Orbital energies (eigenvalues) in eV, sorted ascending.
    pub energies: Vec<f64>,
    /// MO coefficient matrix C (rows = AO index, cols = MO index).
    /// Each column is one molecular orbital.
    pub coefficients: Vec<Vec<f64>>,
    /// Total number of valence electrons.
    pub n_electrons: usize,
    /// Index of the HOMO (0-based).
    pub homo_index: usize,
    /// Index of the LUMO (0-based).
    pub lumo_index: usize,
    /// HOMO energy in eV.
    pub homo_energy: f64,
    /// LUMO energy in eV.
    pub lumo_energy: f64,
    /// HOMO-LUMO gap in eV.
    pub gap: f64,
}

/// Solve the generalized eigenproblem HC = SCE using Löwdin orthogonalization.
///
/// Returns eigenvalues (sorted ascending) and the coefficient matrix C.
pub fn solve_generalized_eigenproblem(
    h: &DMatrix<f64>,
    s: &DMatrix<f64>,
) -> (DVector<f64>, DMatrix<f64>) {
    let n = h.nrows();

    // Step 1: Diagonalize S
    let s_eigen = SymmetricEigen::new(s.clone());
    let s_vals = &s_eigen.eigenvalues;
    let s_vecs = &s_eigen.eigenvectors;

    // Step 2: Build S^{-1/2}
    let mut s_inv_sqrt_diag = DMatrix::zeros(n, n);
    for i in 0..n {
        let val = s_vals[i];
        if val > 1e-10 {
            s_inv_sqrt_diag[(i, i)] = 1.0 / val.sqrt();
        }
    }
    let s_inv_sqrt = s_vecs * &s_inv_sqrt_diag * s_vecs.transpose();

    // Step 3: Transform Hamiltonian: H' = S^{-1/2} H S^{-1/2}
    let h_prime = &s_inv_sqrt * h * &s_inv_sqrt;

    // Step 4: Diagonalize H'
    let h_eigen = SymmetricEigen::new(h_prime);
    let energies = h_eigen.eigenvalues.clone();
    let c_prime = h_eigen.eigenvectors.clone();

    // Step 5: Back-transform C = S^{-1/2} C'
    let c = &s_inv_sqrt * c_prime;

    // Sort by energy (ascending)
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| energies[a].partial_cmp(&energies[b]).unwrap());

    let mut sorted_energies = DVector::zeros(n);
    let mut sorted_c = DMatrix::zeros(n, n);
    for (new_idx, &old_idx) in indices.iter().enumerate() {
        sorted_energies[new_idx] = energies[old_idx];
        for row in 0..n {
            sorted_c[(row, new_idx)] = c[(row, old_idx)];
        }
    }

    (sorted_energies, sorted_c)
}

/// Count valence electrons for a set of atomic numbers.
fn count_valence_electrons(elements: &[u8]) -> usize {
    elements
        .iter()
        .map(|&z| match z {
            1 => 1,  // H
            5 => 3,  // B
            6 => 4,  // C
            7 => 5,  // N
            8 => 6,  // O
            9 => 7,  // F
            14 => 4, // Si
            15 => 5, // P
            16 => 6, // S
            17 => 7, // Cl
            35 => 7, // Br
            53 => 7, // I
            _ => 0,
        })
        .sum()
}

/// Run the full EHT calculation pipeline.
///
/// - `elements`: atomic numbers
/// - `positions`: Cartesian coordinates in Ångström
/// - `k`: Wolfsberg-Helmholtz constant (None = 1.75)
pub fn solve_eht(
    elements: &[u8],
    positions: &[[f64; 3]],
    k: Option<f64>,
) -> Result<EhtResult, String> {
    use super::basis::build_basis;
    use super::hamiltonian::build_hamiltonian;
    use super::overlap::build_overlap_matrix;

    if elements.len() != positions.len() {
        return Err("Element and position arrays must have equal length".to_string());
    }

    let basis = build_basis(elements, positions);
    if basis.is_empty() {
        return Err("No valence orbitals found for given elements".to_string());
    }

    let s = build_overlap_matrix(&basis);
    let h = build_hamiltonian(&basis, &s, k);
    let (energies, c) = solve_generalized_eigenproblem(&h, &s);

    let n_electrons = count_valence_electrons(elements);
    let n_orbitals = basis.len();

    // HOMO is the last occupied orbital (electrons fill in pairs)
    let n_occupied = n_electrons / 2;
    let homo_idx = if n_occupied > 0 { n_occupied - 1 } else { 0 };
    let lumo_idx = if homo_idx + 1 < n_orbitals {
        homo_idx + 1
    } else {
        homo_idx
    };

    let homo_energy = energies[homo_idx];
    let lumo_energy = energies[lumo_idx];

    // Convert nalgebra matrices to Vec<Vec<f64>>
    let coefficients: Vec<Vec<f64>> = (0..n_orbitals)
        .map(|row| (0..n_orbitals).map(|col| c[(row, col)]).collect())
        .collect();

    Ok(EhtResult {
        energies: energies.iter().copied().collect(),
        coefficients,
        n_electrons,
        homo_index: homo_idx,
        lumo_index: lumo_idx,
        homo_energy,
        lumo_energy,
        gap: lumo_energy - homo_energy,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_h2_two_eigenvalues() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let result = solve_eht(&elements, &positions, None).unwrap();
        assert_eq!(result.energies.len(), 2);
        // Bonding orbital lower than anti-bonding
        assert!(result.energies[0] < result.energies[1]);
        // HOMO should be the bonding orbital (index 0)
        assert_eq!(result.homo_index, 0);
        assert_eq!(result.lumo_index, 1);
        assert!(result.gap > 0.0, "H2 HOMO-LUMO gap should be positive");
    }

    #[test]
    fn test_h2_energies_sorted() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let result = solve_eht(&elements, &positions, None).unwrap();
        for i in 1..result.energies.len() {
            assert!(
                result.energies[i] >= result.energies[i - 1],
                "Energies not sorted: E[{}]={} < E[{}]={}",
                i,
                result.energies[i],
                i - 1,
                result.energies[i - 1]
            );
        }
    }

    #[test]
    fn test_h2_coefficients_shape() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let result = solve_eht(&elements, &positions, None).unwrap();
        assert_eq!(result.coefficients.len(), 2);
        assert_eq!(result.coefficients[0].len(), 2);
    }

    #[test]
    fn test_h2o_six_orbitals() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let result = solve_eht(&elements, &positions, None).unwrap();
        // O(2s,2px,2py,2pz) + 2×H(1s) = 6 basis functions
        assert_eq!(result.energies.len(), 6);
        // H2O: 8 valence electrons → 4 occupied orbitals → HOMO index 3
        assert_eq!(result.n_electrons, 8);
        assert_eq!(result.homo_index, 3);
        assert_eq!(result.lumo_index, 4);
    }

    #[test]
    fn test_h2o_gap_positive() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let result = solve_eht(&elements, &positions, None).unwrap();
        assert!(
            result.gap > 0.0,
            "H2O HOMO-LUMO gap = {} should be > 0",
            result.gap
        );
    }

    #[test]
    fn test_lowdin_preserves_orthogonality() {
        // After Löwdin: C^T S C should be identity
        use super::super::basis::build_basis;
        use super::super::hamiltonian::build_hamiltonian;
        use super::super::overlap::build_overlap_matrix;

        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let basis = build_basis(&elements, &positions);
        let s = build_overlap_matrix(&basis);
        let h = build_hamiltonian(&basis, &s, None);
        let (_, c) = solve_generalized_eigenproblem(&h, &s);

        // C^T S C should be approximately identity
        let ct_s_c = c.transpose() * &s * &c;
        let n = ct_s_c.nrows();
        for i in 0..n {
            for j in 0..n {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (ct_s_c[(i, j)] - expected).abs() < 1e-8,
                    "C^T S C[{},{}] = {}, expected {}",
                    i,
                    j,
                    ct_s_c[(i, j)],
                    expected,
                );
            }
        }
    }

    #[test]
    fn test_error_mismatched_arrays() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0]]; // Only 1 position for 2 elements
        assert!(solve_eht(&elements, &positions, None).is_err());
    }

    #[test]
    fn test_valence_electron_count() {
        assert_eq!(count_valence_electrons(&[1, 1]), 2); // H2
        assert_eq!(count_valence_electrons(&[8, 1, 1]), 8); // H2O
        assert_eq!(count_valence_electrons(&[6, 1, 1, 1, 1]), 8); // CH4
        assert_eq!(count_valence_electrons(&[7, 1, 1, 1]), 8); // NH3
    }
}
