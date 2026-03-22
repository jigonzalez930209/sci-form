//! Natural Population Analysis (NPA) and Natural Bond Orbital (NBO) analysis.
//!
//! NPA provides basis-independent atomic charges by constructing Natural
//! Atomic Orbitals (NAOs) that diagonalize the on-atom density blocks.
//! NBO extends this to identify Lewis structure bonding patterns.
//!
//! Reference: Reed, A. E.; Weinstock, R. B.; Weinhold, F.
//! "Natural population analysis." JCP 83 (1985): 735.

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

/// Natural Population Analysis result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NpaResult {
    /// Natural charges per atom (e).
    pub natural_charges: Vec<f64>,
    /// Natural electron configuration per atom (s, p, d contributions).
    pub natural_config: Vec<NaturalConfig>,
    /// Total Rydberg population (electron count in diffuse/Rydberg NAOs).
    pub rydberg_population: f64,
}

/// Natural electron configuration for one atom.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NaturalConfig {
    /// Atom index.
    pub atom_index: usize,
    /// Element atomic number.
    pub element: u8,
    /// s-orbital population.
    pub s_pop: f64,
    /// p-orbital population.
    pub p_pop: f64,
    /// d-orbital population.
    pub d_pop: f64,
    /// Total natural population.
    pub total: f64,
}

/// Natural Bond Orbital result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NboResult {
    /// NPA charges.
    pub npa: NpaResult,
    /// Lewis-structure bond orbitals.
    pub bond_orbitals: Vec<NboBond>,
    /// Lone pair orbitals.
    pub lone_pairs: Vec<NboLonePair>,
    /// Percentage of electrons in the natural Lewis structure.
    pub lewis_population_pct: f64,
}

/// An NBO bonding orbital.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NboBond {
    /// Atom A index.
    pub atom_a: usize,
    /// Atom B index.
    pub atom_b: usize,
    /// Occupancy (ideally ~2.0).
    pub occupancy: f64,
    /// Polarization coefficient on A.
    pub coeff_a: f64,
    /// Polarization coefficient on B.
    pub coeff_b: f64,
    /// Bond type: "sigma", "pi", "delta".
    pub bond_type: String,
}

/// An NBO lone pair orbital.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NboLonePair {
    /// Atom index.
    pub atom_index: usize,
    /// Occupancy (ideally ~2.0).
    pub occupancy: f64,
    /// Orbital type: "s", "p", "sp", "sp2", "sp3".
    pub orbital_type: String,
}

/// Compute Natural Population Analysis from density and overlap matrices.
///
/// The NPA algorithm:
/// 1. Build density matrix: D = C·n·C^T
/// 2. For each atom A, extract the density sub-block D_AA
/// 3. Diagonalize D_AA to get pre-NAOs with natural occupancies
/// 4. Orthogonalize inter-atom overlaps (occupancy-weighted symmetric)
/// 5. Sum resulting occupancies per atom → natural charges
pub fn compute_npa(
    elements: &[u8],
    overlap: &DMatrix<f64>,
    density: &DMatrix<f64>,
    basis_atom_map: &[usize],
) -> Result<NpaResult, String> {
    let n_atoms = elements.len();
    let n_basis = overlap.nrows();

    if density.nrows() != n_basis || density.ncols() != n_basis {
        return Err("Density matrix dimension mismatch".to_string());
    }
    if basis_atom_map.len() != n_basis {
        return Err("basis_atom_map length must match basis size".to_string());
    }

    // Step 1: Compute S^{-1/2} for Löwdin orthogonalization
    let s_eigen = nalgebra::SymmetricEigen::new(overlap.clone());
    let mut s_inv_sqrt = DMatrix::zeros(n_basis, n_basis);
    for i in 0..n_basis {
        let ev = s_eigen.eigenvalues[i];
        if ev > 1e-10 {
            s_inv_sqrt[(i, i)] = 1.0 / ev.sqrt();
        }
    }
    let s_half_inv = &s_eigen.eigenvectors * &s_inv_sqrt * s_eigen.eigenvectors.transpose();

    // Step 2: Transform density to orthogonal basis
    // D_orth = S^{-1/2} · D · S^{-1/2}
    let d_orth = &s_half_inv * density * &s_half_inv;

    // Step 3: Compute per-atom natural populations
    let mut natural_charges = vec![0.0; n_atoms];
    let mut natural_configs = Vec::with_capacity(n_atoms);
    let mut total_rydberg = 0.0;

    for atom in 0..n_atoms {
        // Find basis functions on this atom
        let atom_basis: Vec<usize> = (0..n_basis)
            .filter(|&mu| basis_atom_map[mu] == atom)
            .collect();

        // Extract atomic density sub-block
        let n_atom_basis = atom_basis.len();
        if n_atom_basis == 0 {
            natural_configs.push(NaturalConfig {
                atom_index: atom,
                element: elements[atom],
                s_pop: 0.0,
                p_pop: 0.0,
                d_pop: 0.0,
                total: 0.0,
            });
            continue;
        }

        let mut d_aa = DMatrix::zeros(n_atom_basis, n_atom_basis);
        for (ii, &mu) in atom_basis.iter().enumerate() {
            for (jj, &nu) in atom_basis.iter().enumerate() {
                d_aa[(ii, jj)] = d_orth[(mu, nu)];
            }
        }

        // Diagonalize atomic block → pre-NAO occupancies
        let eigen_aa = nalgebra::SymmetricEigen::new(d_aa);
        let total_pop: f64 = eigen_aa.eigenvalues.iter().sum();

        // Classify occupancies by angular momentum (simplified)
        let z = elements[atom];
        let (s_pop, p_pop, d_pop) = classify_shell_populations(&eigen_aa.eigenvalues, z);

        // Determine "core" + "valence" vs "Rydberg"
        let expected_valence = expected_valence_electrons(z);
        if total_pop > expected_valence as f64 + 0.5 {
            total_rydberg += total_pop - expected_valence as f64;
        }

        let nuclear_charge = z as f64;
        natural_charges[atom] = nuclear_charge - total_pop;

        natural_configs.push(NaturalConfig {
            atom_index: atom,
            element: z,
            s_pop,
            p_pop,
            d_pop,
            total: total_pop,
        });
    }

    Ok(NpaResult {
        natural_charges,
        natural_config: natural_configs,
        rydberg_population: total_rydberg,
    })
}

/// Compute NBO analysis (NPA + bond/lone-pair identification).
///
/// Starting from NPA results, identifies Lewis-structure bonding patterns
/// by analyzing the inter-atomic density blocks.
pub fn compute_nbo(
    elements: &[u8],
    overlap: &DMatrix<f64>,
    density: &DMatrix<f64>,
    basis_atom_map: &[usize],
    bonds: &[(usize, usize)],
) -> Result<NboResult, String> {
    let npa = compute_npa(elements, overlap, density, basis_atom_map)?;
    let n_basis = overlap.nrows();

    // Identify bond orbitals from inter-atomic density blocks
    let mut bond_orbitals = Vec::new();
    let mut lone_pairs = Vec::new();

    for &(a, b) in bonds {
        let a_basis: Vec<usize> = (0..n_basis).filter(|&mu| basis_atom_map[mu] == a).collect();
        let b_basis: Vec<usize> = (0..n_basis).filter(|&mu| basis_atom_map[mu] == b).collect();

        if a_basis.is_empty() || b_basis.is_empty() {
            continue;
        }

        // Extract AB density block
        let na = a_basis.len();
        let nb = b_basis.len();
        let mut d_ab = DMatrix::zeros(na + nb, na + nb);

        for (ii, &mu) in a_basis.iter().chain(b_basis.iter()).enumerate() {
            for (jj, &nu) in a_basis.iter().chain(b_basis.iter()).enumerate() {
                d_ab[(ii, jj)] = density[(mu, nu)];
            }
        }

        // Diagonalize AB block
        let eigen_ab = nalgebra::SymmetricEigen::new(d_ab);

        // Bond orbitals: eigenvalues near 2.0
        for (k, &occ) in eigen_ab.eigenvalues.iter().enumerate() {
            if occ > 1.5 {
                let col = eigen_ab.eigenvectors.column(k);
                let a_weight: f64 = col.rows(0, na).iter().map(|c| c * c).sum();
                let b_weight: f64 = col.rows(na, nb).iter().map(|c| c * c).sum();
                let total = a_weight + b_weight;

                bond_orbitals.push(NboBond {
                    atom_a: a,
                    atom_b: b,
                    occupancy: occ,
                    coeff_a: (a_weight / (total + 1e-30)).sqrt(),
                    coeff_b: (b_weight / (total + 1e-30)).sqrt(),
                    bond_type: "sigma".to_string(),
                });
            }
        }
    }

    // Identify lone pairs from single-atom density blocks
    let n_atoms = elements.len();
    for atom in 0..n_atoms {
        let atom_basis: Vec<usize> = (0..n_basis)
            .filter(|&mu| basis_atom_map[mu] == atom)
            .collect();

        if atom_basis.is_empty() {
            continue;
        }

        let na = atom_basis.len();
        let mut d_aa = DMatrix::zeros(na, na);
        for (ii, &mu) in atom_basis.iter().enumerate() {
            for (jj, &nu) in atom_basis.iter().enumerate() {
                d_aa[(ii, jj)] = density[(mu, nu)];
            }
        }

        let eigen_aa = nalgebra::SymmetricEigen::new(d_aa);

        // Lone pairs: high occupancy NAOs not involved in bonding
        let bonded_count = bonds
            .iter()
            .filter(|&&(a, b)| a == atom || b == atom)
            .count();

        for (k, &occ) in eigen_aa.eigenvalues.iter().enumerate() {
            if occ > 1.8 && k >= bonded_count {
                let orbital_type = if na == 1 {
                    "s"
                } else if k == 0 {
                    "sp"
                } else {
                    "p"
                };
                lone_pairs.push(NboLonePair {
                    atom_index: atom,
                    occupancy: occ,
                    orbital_type: orbital_type.to_string(),
                });
            }
        }
    }

    // Lewis population percentage
    let total_electrons: f64 = npa.natural_config.iter().map(|c| c.total).sum();
    let lewis_pop: f64 = bond_orbitals.iter().map(|b| b.occupancy).sum::<f64>()
        + lone_pairs.iter().map(|l| l.occupancy).sum::<f64>();
    let lewis_pct = if total_electrons > 0.0 {
        100.0 * lewis_pop / total_electrons
    } else {
        0.0
    };

    Ok(NboResult {
        npa,
        bond_orbitals,
        lone_pairs,
        lewis_population_pct: lewis_pct,
    })
}

/// Classify eigenvalues by angular momentum shell.
fn classify_shell_populations(eigenvalues: &nalgebra::DVector<f64>, z: u8) -> (f64, f64, f64) {
    let mut sorted: Vec<f64> = eigenvalues.iter().copied().collect();
    sorted.sort_by(|a, b| b.partial_cmp(a).unwrap());

    let (n_s, n_p, n_d) = shell_count(z);

    let mut s_pop = 0.0;
    let mut p_pop = 0.0;
    let mut d_pop = 0.0;

    for (i, &occ) in sorted.iter().enumerate() {
        if occ < 0.0 {
            continue;
        }
        if i < n_s {
            s_pop += occ;
        } else if i < n_s + n_p {
            p_pop += occ;
        } else if i < n_s + n_p + n_d {
            d_pop += occ;
        }
    }

    (s_pop, p_pop, d_pop)
}

/// Expected number of valence electrons for an element.
fn expected_valence_electrons(z: u8) -> usize {
    match z {
        1 => 1,
        2 => 2,
        3..=4 => (z - 2) as usize,
        5..=10 => (z - 2) as usize,
        11..=12 => (z - 10) as usize,
        13..=18 => (z - 10) as usize,
        19..=20 => (z - 18) as usize,
        21..=30 => (z - 18) as usize, // TMs
        31..=36 => (z - 28) as usize,
        _ => z.min(8) as usize,
    }
}

/// Expected shell counts (s, p, d) for the minimal basis of element z.
fn shell_count(z: u8) -> (usize, usize, usize) {
    match z {
        1 => (1, 0, 0),
        2 => (1, 0, 0),
        3..=4 => (2, 0, 0),
        5..=10 => (2, 3, 0),
        11..=12 => (3, 0, 0),
        13..=18 => (3, 3, 0),
        19..=20 => (4, 0, 0),
        21..=30 => (4, 3, 5),
        31..=36 => (4, 3, 5),
        _ => (1, 0, 0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_npa_simple() {
        // 2-basis system: one function on each atom
        let elements = vec![1u8, 1];
        let overlap = DMatrix::from_row_slice(2, 2, &[1.0, 0.5, 0.5, 1.0]);
        let density = DMatrix::from_row_slice(2, 2, &[0.5, 0.4, 0.4, 0.5]);
        let basis_atom_map = vec![0, 1];

        let result = compute_npa(&elements, &overlap, &density, &basis_atom_map);
        assert!(result.is_ok());
        let npa = result.unwrap();
        assert_eq!(npa.natural_charges.len(), 2);
        // H2: symmetric → charges should be near zero
        assert!((npa.natural_charges[0] - npa.natural_charges[1]).abs() < 0.01);
    }

    #[test]
    fn test_expected_valence() {
        assert_eq!(expected_valence_electrons(1), 1);
        assert_eq!(expected_valence_electrons(6), 4);
        assert_eq!(expected_valence_electrons(8), 6);
        assert_eq!(expected_valence_electrons(26), 8); // Fe
    }

    #[test]
    fn test_shell_count() {
        assert_eq!(shell_count(1), (1, 0, 0));
        assert_eq!(shell_count(6), (2, 3, 0));
        assert_eq!(shell_count(26), (4, 3, 5));
    }
}
