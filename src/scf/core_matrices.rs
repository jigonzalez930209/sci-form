//! Core one-electron matrices: S, T, V, and H⁰ = T + V.
//!
//! Groups all one-electron integrals into a single coherent structure
//! for use throughout the SCF pipeline.

use nalgebra::DMatrix;

use super::basis::BasisSet;
use super::kinetic_matrix::build_kinetic_matrix;
use super::nuclear_matrix::build_nuclear_matrix;
use super::overlap_matrix::build_overlap_matrix;

/// All one-electron matrices needed for an SCF calculation.
#[derive(Debug, Clone)]
pub struct CoreMatrices {
    /// Overlap matrix S (n_basis × n_basis).
    pub overlap: DMatrix<f64>,
    /// Kinetic energy matrix T (n_basis × n_basis).
    pub kinetic: DMatrix<f64>,
    /// Nuclear attraction matrix V (n_basis × n_basis).
    pub nuclear: DMatrix<f64>,
    /// Core Hamiltonian H⁰ = T + V (n_basis × n_basis).
    pub core_hamiltonian: DMatrix<f64>,
    /// Number of basis functions.
    pub n_basis: usize,
}

impl CoreMatrices {
    /// Build all one-electron matrices for the molecular system.
    ///
    /// Computes S, T, V, and H⁰ = T + V.
    pub fn build(
        basis: &BasisSet,
        elements: &[u8],
        positions_bohr: &[[f64; 3]],
    ) -> Self {
        let s = build_overlap_matrix(basis);
        let t = build_kinetic_matrix(basis);
        let v = build_nuclear_matrix(basis, elements, positions_bohr);
        let h_core = &t + &v;

        CoreMatrices {
            overlap: s,
            kinetic: t,
            nuclear: v,
            core_hamiltonian: h_core,
            n_basis: basis.n_basis,
        }
    }
}

/// Compute nuclear repulsion energy.
///
/// E_nuc = Σ_{A>B} Z_A · Z_B / |R_A - R_B|
pub fn nuclear_repulsion_energy(elements: &[u8], positions_bohr: &[[f64; 3]]) -> f64 {
    let n_atoms = elements.len();
    let mut e_nuc = 0.0;

    for a in 0..n_atoms {
        for b in (a + 1)..n_atoms {
            let za = elements[a] as f64;
            let zb = elements[b] as f64;
            let dx = positions_bohr[a][0] - positions_bohr[b][0];
            let dy = positions_bohr[a][1] - positions_bohr[b][1];
            let dz = positions_bohr[a][2] - positions_bohr[b][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            e_nuc += za * zb / r;
        }
    }

    e_nuc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_core_hamiltonian_symmetric() {
        let elements = [8u8, 1, 1];
        let positions = [
            [0.0, 0.0, 0.2216],
            [0.0, 1.4313, -0.8864],
            [0.0, -1.4313, -0.8864],
        ];
        let basis = BasisSet::sto3g(&elements, &positions);
        let matrices = CoreMatrices::build(&basis, &elements, &positions);

        let h = &matrices.core_hamiltonian;
        for i in 0..matrices.n_basis {
            for j in 0..matrices.n_basis {
                assert!(
                    (h[(i, j)] - h[(j, i)]).abs() < 1e-12,
                    "H⁰ not symmetric at ({}, {})", i, j
                );
            }
        }
    }

    #[test]
    fn test_h_core_equals_t_plus_v() {
        let basis = BasisSet::sto3g(&[1, 1], &[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]);
        let matrices = CoreMatrices::build(&basis, &[1, 1], &[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]);

        for i in 0..matrices.n_basis {
            for j in 0..matrices.n_basis {
                let expected = matrices.kinetic[(i, j)] + matrices.nuclear[(i, j)];
                assert!(
                    (matrices.core_hamiltonian[(i, j)] - expected).abs() < 1e-14,
                    "H⁰ ≠ T + V at ({}, {})", i, j
                );
            }
        }
    }

    #[test]
    fn test_nuclear_repulsion_h2() {
        let e_nuc = nuclear_repulsion_energy(
            &[1, 1],
            &[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]],
        );
        assert!((e_nuc - 1.0 / 1.4).abs() < 1e-10);
    }

    #[test]
    fn test_overlap_diagonal_positive() {
        let elements = [8u8, 1, 1];
        let positions = [
            [0.0, 0.0, 0.2216],
            [0.0, 1.4313, -0.8864],
            [0.0, -1.4313, -0.8864],
        ];
        let basis = BasisSet::sto3g(&elements, &positions);
        let matrices = CoreMatrices::build(&basis, &elements, &positions);

        for i in 0..matrices.n_basis {
            assert!(matrices.overlap[(i, i)] > 0.0, "S[{i},{i}] should be > 0");
        }
    }

    #[test]
    fn test_core_matrices_water_dimensions() {
        let elements = [8u8, 1, 1];
        let positions = [
            [0.0, 0.0, 0.0],
            [1.43, 1.11, 0.0],
            [-1.43, 1.11, 0.0],
        ];
        let basis = BasisSet::sto3g(&elements, &positions);
        let matrices = CoreMatrices::build(&basis, &elements, &positions);

        assert_eq!(matrices.n_basis, 7);
        assert_eq!(matrices.overlap.nrows(), 7);
        assert_eq!(matrices.kinetic.nrows(), 7);
        assert_eq!(matrices.nuclear.nrows(), 7);
        assert_eq!(matrices.core_hamiltonian.nrows(), 7);
    }
}
