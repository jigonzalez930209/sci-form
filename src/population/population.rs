//! Mulliken and Löwdin population analysis.
//!
//! **Mulliken:** q_A = Z_A − Σ_{μ∈A} (PS)_{μμ}
//! **Löwdin:**  q_A = Z_A − Σ_{μ∈A} (S^{1/2} P S^{1/2})_{μμ}
//!
//! where P is the density matrix P = Σ_i^{occ} n_i c_i c_i^T

use nalgebra::{DMatrix, SymmetricEigen};
use serde::{Deserialize, Serialize};

/// Result of a population analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PopulationResult {
    /// Mulliken partial charges per atom.
    pub mulliken_charges: Vec<f64>,
    /// Löwdin partial charges per atom.
    pub lowdin_charges: Vec<f64>,
    /// Mulliken gross orbital populations per AO.
    pub mulliken_populations: Vec<f64>,
    /// Number of atoms.
    pub num_atoms: usize,
    /// Total charge (should match net formal charge).
    pub total_charge_mulliken: f64,
    /// Total charge from Löwdin.
    pub total_charge_lowdin: f64,
}

/// Build the density matrix P from MO coefficients and occupations.
///
/// P_{μν} = Σ_i n_i C_{μi} C_{νi}
///
/// - `coefficients`: rows = AO, cols = MO (same layout as EhtResult.coefficients)
/// - `n_electrons`: total electrons (fills lowest MOs, 2 per orbital)
fn build_density_matrix(coefficients: &[Vec<f64>], n_electrons: usize) -> DMatrix<f64> {
    let n_ao = coefficients.len();
    let n_occupied = n_electrons / 2;
    let mut p = DMatrix::zeros(n_ao, n_ao);

    for i in 0..n_occupied {
        for mu in 0..n_ao {
            for nu in 0..n_ao {
                p[(mu, nu)] += 2.0 * coefficients[mu][i] * coefficients[nu][i];
            }
        }
    }

    // Handle odd electron (singly occupied HOMO)
    if n_electrons % 2 == 1 {
        let i = n_occupied;
        for mu in 0..n_ao {
            for nu in 0..n_ao {
                p[(mu, nu)] += coefficients[mu][i] * coefficients[nu][i];
            }
        }
    }

    p
}

/// Count valence electrons for an element.
fn valence_electrons(z: u8) -> f64 {
    match z {
        1 => 1.0,
        5 => 3.0,
        6 => 4.0,
        7 => 5.0,
        8 => 6.0,
        9 => 7.0,
        14 => 4.0,
        15 => 5.0,
        16 => 6.0,
        17 => 7.0,
        35 => 7.0,
        53 => 7.0,
        _ => 0.0,
    }
}

/// Map each AO index to its parent atom index.
fn ao_to_atom_map(basis: &[crate::eht::basis::AtomicOrbital]) -> Vec<usize> {
    basis.iter().map(|ao| ao.atom_index).collect()
}

/// Compute Mulliken charges.
///
/// Mulliken population for AO μ: (PS)_{μμ}
/// Atom charge: q_A = Z_A − Σ_{μ∈A} (PS)_{μμ}
pub fn mulliken_charges(
    elements: &[u8],
    basis: &[crate::eht::basis::AtomicOrbital],
    overlap: &DMatrix<f64>,
    coefficients: &[Vec<f64>],
    n_electrons: usize,
) -> Vec<f64> {
    let n_ao = basis.len();
    let n_atoms = elements.len();
    let ao_map = ao_to_atom_map(basis);
    let p = build_density_matrix(coefficients, n_electrons);

    // PS product
    let ps = &p * overlap;

    // Gross AO populations
    let mut atom_pop = vec![0.0; n_atoms];
    for mu in 0..n_ao {
        atom_pop[ao_map[mu]] += ps[(mu, mu)];
    }

    // Charges = nuclear valence - electron population
    (0..n_atoms)
        .map(|a| valence_electrons(elements[a]) - atom_pop[a])
        .collect()
}

/// Compute Löwdin charges.
///
/// Löwdin uses S^{1/2} P S^{1/2} instead of PS.
/// q_A = Z_A − Σ_{μ∈A} (S^{1/2} P S^{1/2})_{μμ}
pub fn lowdin_charges(
    elements: &[u8],
    basis: &[crate::eht::basis::AtomicOrbital],
    overlap: &DMatrix<f64>,
    coefficients: &[Vec<f64>],
    n_electrons: usize,
) -> Vec<f64> {
    let n_ao = basis.len();
    let n_atoms = elements.len();
    let ao_map = ao_to_atom_map(basis);
    let p = build_density_matrix(coefficients, n_electrons);

    // Build S^{1/2}
    let s_eigen = SymmetricEigen::new(overlap.clone());
    let mut s_sqrt_diag = DMatrix::zeros(n_ao, n_ao);
    for i in 0..n_ao {
        let val = s_eigen.eigenvalues[i];
        if val > 1e-10 {
            s_sqrt_diag[(i, i)] = val.sqrt();
        }
    }
    let s_sqrt = &s_eigen.eigenvectors * &s_sqrt_diag * s_eigen.eigenvectors.transpose();

    // S^{1/2} P S^{1/2}
    let sps = &s_sqrt * &p * &s_sqrt;

    let mut atom_pop = vec![0.0; n_atoms];
    for mu in 0..n_ao {
        atom_pop[ao_map[mu]] += sps[(mu, mu)];
    }

    (0..n_atoms)
        .map(|a| valence_electrons(elements[a]) - atom_pop[a])
        .collect()
}

/// Full population analysis: Mulliken + Löwdin in one pass.
pub fn compute_population(
    elements: &[u8],
    positions: &[[f64; 3]],
    coefficients: &[Vec<f64>],
    n_electrons: usize,
) -> PopulationResult {
    let basis = crate::eht::basis::build_basis(elements, positions);
    let overlap = crate::eht::build_overlap_matrix(&basis);
    let n_ao = basis.len();
    let n_atoms = elements.len();
    let ao_map = ao_to_atom_map(&basis);
    let p = build_density_matrix(coefficients, n_electrons);

    // Mulliken: PS
    let ps = &p * &overlap;
    let mut mulliken_pop = vec![0.0; n_atoms];
    let mut mulliken_ao_pop = vec![0.0; n_ao];
    for mu in 0..n_ao {
        mulliken_ao_pop[mu] = ps[(mu, mu)];
        mulliken_pop[ao_map[mu]] += ps[(mu, mu)];
    }
    let mulliken_charges: Vec<f64> = (0..n_atoms)
        .map(|a| valence_electrons(elements[a]) - mulliken_pop[a])
        .collect();
    let total_mulliken: f64 = mulliken_charges.iter().sum();

    // Löwdin: S^{1/2} P S^{1/2}
    let s_eigen = SymmetricEigen::new(overlap.clone());
    let mut s_sqrt_diag = DMatrix::zeros(n_ao, n_ao);
    for i in 0..n_ao {
        let val = s_eigen.eigenvalues[i];
        if val > 1e-10 {
            s_sqrt_diag[(i, i)] = val.sqrt();
        }
    }
    let s_sqrt = &s_eigen.eigenvectors * &s_sqrt_diag * s_eigen.eigenvectors.transpose();
    let sps = &s_sqrt * &p * &s_sqrt;
    let mut lowdin_pop = vec![0.0; n_atoms];
    for mu in 0..n_ao {
        lowdin_pop[ao_map[mu]] += sps[(mu, mu)];
    }
    let lowdin_charges: Vec<f64> = (0..n_atoms)
        .map(|a| valence_electrons(elements[a]) - lowdin_pop[a])
        .collect();
    let total_lowdin: f64 = lowdin_charges.iter().sum();

    PopulationResult {
        mulliken_charges,
        lowdin_charges,
        mulliken_populations: mulliken_ao_pop,
        num_atoms: n_atoms,
        total_charge_mulliken: total_mulliken,
        total_charge_lowdin: total_lowdin,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eht::solve_eht;

    fn h2_molecule() -> (Vec<u8>, Vec<[f64; 3]>) {
        (vec![1, 1], vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]])
    }

    fn water_molecule() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![8, 1, 1],
            vec![[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]],
        )
    }

    #[test]
    fn test_h2_symmetric_charges() {
        let (elems, pos) = h2_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        let pop = compute_population(&elems, &pos, &result.coefficients, result.n_electrons);

        // H₂ is symmetric → both atoms should have ~0 charge
        assert!(
            (pop.mulliken_charges[0] - pop.mulliken_charges[1]).abs() < 1e-6,
            "H₂ Mulliken charges should be symmetric"
        );
        assert!(
            pop.mulliken_charges[0].abs() < 0.01,
            "H₂ charge should be ~0"
        );
    }

    #[test]
    fn test_water_oxygen_negative() {
        let (elems, pos) = water_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        let pop = compute_population(&elems, &pos, &result.coefficients, result.n_electrons);

        // Oxygen should be negative (more electronegative)
        assert!(
            pop.mulliken_charges[0] < 0.0,
            "O in water should have negative Mulliken charge, got {}",
            pop.mulliken_charges[0]
        );
        assert!(
            pop.lowdin_charges[0] < 0.0,
            "O in water should have negative Löwdin charge, got {}",
            pop.lowdin_charges[0]
        );
    }

    #[test]
    fn test_charge_sum_conservation() {
        let (elems, pos) = water_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        let pop = compute_population(&elems, &pos, &result.coefficients, result.n_electrons);

        // Neutral molecule → charges sum ~ 0
        assert!(
            pop.total_charge_mulliken.abs() < 0.01,
            "Mulliken total charge should be ~0, got {}",
            pop.total_charge_mulliken
        );
        assert!(
            pop.total_charge_lowdin.abs() < 0.01,
            "Löwdin total charge should be ~0, got {}",
            pop.total_charge_lowdin
        );
    }

    #[test]
    fn test_hydrogen_symmetry_in_water() {
        let (elems, pos) = water_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        let pop = compute_population(&elems, &pos, &result.coefficients, result.n_electrons);

        // Both H atoms should have identical charges (symmetric geometry)
        assert!(
            (pop.mulliken_charges[1] - pop.mulliken_charges[2]).abs() < 0.01,
            "H charges in water should be symmetric: {} vs {}",
            pop.mulliken_charges[1],
            pop.mulliken_charges[2]
        );
    }

    #[test]
    fn test_lowdin_vs_mulliken_different() {
        let (elems, pos) = water_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        let pop = compute_population(&elems, &pos, &result.coefficients, result.n_electrons);

        // Löwdin and Mulliken typically give different charge magnitudes
        // but same sign for electronegativity ordering
        let m_o = pop.mulliken_charges[0];
        let l_o = pop.lowdin_charges[0];
        assert!(
            m_o.signum() == l_o.signum(),
            "Both methods should agree on sign for O"
        );
    }

    #[test]
    fn test_gross_orbital_populations_sum() {
        let (elems, pos) = h2_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        let pop = compute_population(&elems, &pos, &result.coefficients, result.n_electrons);

        // Sum of gross orbital populations should equal n_electrons
        let total: f64 = pop.mulliken_populations.iter().sum();
        assert!(
            (total - result.n_electrons as f64).abs() < 0.01,
            "AO pop sum {} should equal n_electrons {}",
            total,
            result.n_electrons
        );
    }
}
