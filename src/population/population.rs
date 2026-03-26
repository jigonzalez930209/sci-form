//! Mulliken and Löwdin population analysis.
//!
//! **Mulliken:** q_A = Z_A − Σ_{μ∈A} (PS)_{μμ}
//! **Löwdin:**  q_A = Z_A − Σ_{μ∈A} (S^{1/2} P S^{1/2})_{μμ}
//!
//! where P is the density matrix P = Σ_i^{occ} n_i c_i c_i^T

use nalgebra::{DMatrix, SymmetricEigen};
use serde::{Deserialize, Serialize};

/// Bond-order metrics for one atom pair.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BondOrderEntry {
    /// First atom index.
    pub atom_i: usize,
    /// Second atom index.
    pub atom_j: usize,
    /// Interatomic distance in Å.
    pub distance: f64,
    /// Wiberg-like bond index computed in a Löwdin-orthogonalized AO basis.
    pub wiberg: f64,
    /// Mayer-like bond index computed from the PS matrix.
    pub mayer: f64,
}

/// Full bond-order analysis across all atom pairs.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BondOrderResult {
    /// Bond-order metrics for each unique atom pair.
    pub bonds: Vec<BondOrderEntry>,
    /// Number of atoms in the system.
    pub num_atoms: usize,
    /// Sum of Wiberg-like bond indices touching each atom.
    pub wiberg_valence: Vec<f64>,
    /// Sum of Mayer-like bond indices touching each atom.
    pub mayer_valence: Vec<f64>,
}

/// Result of a population analysis.
///
/// Contains Mulliken and Löwdin charges plus a charge-conservation check.
/// `charge_conservation_error` gives the absolute deviation of the
/// Mulliken total charge from the expected net charge. A value above
/// ~0.01 e suggests a basis-set or occupancy issue.
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
    /// Absolute deviation of Mulliken total charge from integer expectation.
    /// Small values (~0) indicate proper charge conservation.
    pub charge_conservation_error: f64,
}

/// Build the density matrix P from MO coefficients and occupations.
///
/// P_{μν} = Σ_i n_i C_{μi} C_{νi}
///
/// - `coefficients`: rows = AO, cols = MO (same layout as EhtResult.coefficients)
/// - `n_electrons`: total electrons (fills lowest MOs, 2 per orbital)
pub(crate) fn build_density_matrix(coefficients: &[Vec<f64>], n_electrons: usize) -> DMatrix<f64> {
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
/// Effective valence electron count for EHT-based population analysis.
pub(crate) fn valence_electrons(z: u8) -> f64 {
    match z {
        // Period 1
        1 => 1.0, // H
        2 => 2.0, // He
        // Period 2
        3 => 1.0,  // Li
        4 => 2.0,  // Be
        5 => 3.0,  // B
        6 => 4.0,  // C
        7 => 5.0,  // N
        8 => 6.0,  // O
        9 => 7.0,  // F
        10 => 8.0, // Ne
        // Period 3
        11 => 1.0, // Na
        12 => 2.0, // Mg
        13 => 3.0, // Al
        14 => 4.0, // Si
        15 => 5.0, // P
        16 => 6.0, // S
        17 => 7.0, // Cl
        18 => 8.0, // Ar
        // Period 4 main-group
        19 => 1.0, // K
        20 => 2.0, // Ca
        31 => 3.0, // Ga
        32 => 4.0, // Ge
        33 => 5.0, // As
        34 => 6.0, // Se
        35 => 7.0, // Br
        36 => 8.0, // Kr
        // Period 5 main-group
        37 => 1.0, // Rb
        38 => 2.0, // Sr
        49 => 3.0, // In
        50 => 4.0, // Sn
        51 => 5.0, // Sb
        52 => 6.0, // Te
        53 => 7.0, // I
        54 => 8.0, // Xe
        // 3d transition metals (valence = 4s + 3d electrons)
        21 => 3.0,  // Sc
        22 => 4.0,  // Ti
        23 => 5.0,  // V
        24 => 6.0,  // Cr
        25 => 7.0,  // Mn
        26 => 8.0,  // Fe
        27 => 9.0,  // Co
        28 => 10.0, // Ni
        29 => 11.0, // Cu
        30 => 12.0, // Zn
        // 4d transition metals
        39 => 3.0,  // Y
        40 => 4.0,  // Zr
        41 => 5.0,  // Nb
        42 => 6.0,  // Mo
        43 => 7.0,  // Tc
        44 => 8.0,  // Ru
        45 => 9.0,  // Rh
        46 => 10.0, // Pd
        47 => 11.0, // Ag
        48 => 12.0, // Cd
        // 5d transition metals
        72 => 4.0,  // Hf
        73 => 5.0,  // Ta
        74 => 6.0,  // W
        75 => 7.0,  // Re
        76 => 8.0,  // Os
        77 => 9.0,  // Ir
        78 => 10.0, // Pt
        79 => 11.0, // Au
        80 => 12.0, // Hg
        // Period 6 main-group
        81 => 3.0, // Tl
        82 => 4.0, // Pb
        83 => 5.0, // Bi
        _ => 0.0,
    }
}

/// Map each AO index to its parent atom index.
pub(crate) fn ao_to_atom_map(basis: &[crate::eht::basis::AtomicOrbital]) -> Vec<usize> {
    basis.iter().map(|ao| ao.atom_index).collect()
}

fn lowdin_orthogonalized_density(overlap: &DMatrix<f64>, density: &DMatrix<f64>) -> DMatrix<f64> {
    let n_ao = overlap.nrows();
    let s_eigen = SymmetricEigen::new(overlap.clone());
    let mut s_sqrt_diag = DMatrix::zeros(n_ao, n_ao);
    for i in 0..n_ao {
        let val = s_eigen.eigenvalues[i];
        if val > 1e-10 {
            s_sqrt_diag[(i, i)] = val.sqrt();
        }
    }
    let s_sqrt = &s_eigen.eigenvectors * &s_sqrt_diag * s_eigen.eigenvectors.transpose();
    &s_sqrt * density * &s_sqrt
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
    let sps = lowdin_orthogonalized_density(overlap, &p);

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
    let sps = lowdin_orthogonalized_density(&overlap, &p);
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
        charge_conservation_error: (total_mulliken - total_mulliken.round()).abs(),
    }
}

/// Compute Wiberg-like and Mayer-like bond orders from an EHT density matrix.
pub fn compute_bond_orders(
    elements: &[u8],
    positions: &[[f64; 3]],
    coefficients: &[Vec<f64>],
    n_electrons: usize,
) -> BondOrderResult {
    let basis = crate::eht::basis::build_basis(elements, positions);
    let overlap = crate::eht::build_overlap_matrix(&basis);
    let ao_map = ao_to_atom_map(&basis);
    let density = build_density_matrix(coefficients, n_electrons);
    let ps = &density * &overlap;
    let p_orth = lowdin_orthogonalized_density(&overlap, &density);

    let mut atom_aos = vec![Vec::new(); elements.len()];
    for (ao_idx, &atom_idx) in ao_map.iter().enumerate() {
        atom_aos[atom_idx].push(ao_idx);
    }

    let mut bonds = Vec::new();
    let mut wiberg_valence = vec![0.0; elements.len()];
    let mut mayer_valence = vec![0.0; elements.len()];

    for atom_i in 0..elements.len() {
        for atom_j in (atom_i + 1)..elements.len() {
            let mut wiberg = 0.0;
            let mut mayer = 0.0;

            for &mu in &atom_aos[atom_i] {
                for &nu in &atom_aos[atom_j] {
                    let p_orth_mn = p_orth[(mu, nu)];
                    wiberg += p_orth_mn * p_orth_mn;
                    mayer += ps[(mu, nu)] * ps[(nu, mu)];
                }
            }

            let dx = positions[atom_i][0] - positions[atom_j][0];
            let dy = positions[atom_i][1] - positions[atom_j][1];
            let dz = positions[atom_i][2] - positions[atom_j][2];
            let distance = (dx * dx + dy * dy + dz * dz).sqrt();

            wiberg_valence[atom_i] += wiberg;
            wiberg_valence[atom_j] += wiberg;
            mayer_valence[atom_i] += mayer;
            mayer_valence[atom_j] += mayer;

            bonds.push(BondOrderEntry {
                atom_i,
                atom_j,
                distance,
                wiberg,
                mayer,
            });
        }
    }

    BondOrderResult {
        bonds,
        num_atoms: elements.len(),
        wiberg_valence,
        mayer_valence,
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

    #[test]
    fn test_h2_bond_order_is_positive() {
        let (elems, pos) = h2_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        let bond_orders =
            compute_bond_orders(&elems, &pos, &result.coefficients, result.n_electrons);

        assert_eq!(bond_orders.bonds.len(), 1);
        assert!(bond_orders.bonds[0].wiberg > 0.1);
        assert!(bond_orders.bonds[0].mayer > 0.1);
    }

    #[test]
    fn test_water_oh_bonds_exceed_hh_bond_order() {
        let (elems, pos) = water_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        let bond_orders =
            compute_bond_orders(&elems, &pos, &result.coefficients, result.n_electrons);

        let oh_1 = bond_orders
            .bonds
            .iter()
            .find(|bond| bond.atom_i == 0 && bond.atom_j == 1)
            .unwrap();
        let oh_2 = bond_orders
            .bonds
            .iter()
            .find(|bond| bond.atom_i == 0 && bond.atom_j == 2)
            .unwrap();
        let hh = bond_orders
            .bonds
            .iter()
            .find(|bond| bond.atom_i == 1 && bond.atom_j == 2)
            .unwrap();

        assert!(oh_1.wiberg > hh.wiberg);
        assert!(oh_2.wiberg > hh.wiberg);
        assert!(oh_1.mayer > hh.mayer);
        assert!(oh_2.mayer > hh.mayer);
    }
}
