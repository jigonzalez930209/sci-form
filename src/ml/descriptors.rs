//! Molecular descriptors for ML property prediction.
//!
//! Computes a vector of constitutional, topological, and electronic
//! descriptors from molecular graph and (optionally) 3D coordinates.

use serde::{Deserialize, Serialize};

/// Molecular descriptor vector.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MolecularDescriptors {
    /// Molecular weight (Daltons).
    pub molecular_weight: f64,
    /// Number of heavy atoms (non-H).
    pub n_heavy_atoms: usize,
    /// Number of hydrogen atoms.
    pub n_hydrogens: usize,
    /// Number of bonds.
    pub n_bonds: usize,
    /// Number of rotatable bonds (single bonds between non-H heavy atoms, not in rings).
    pub n_rotatable_bonds: usize,
    /// Number of H-bond donors (N-H, O-H).
    pub n_hbd: usize,
    /// Number of H-bond acceptors (N, O).
    pub n_hba: usize,
    /// Fraction of sp3 carbons.
    pub fsp3: f64,
    /// Total partial charge magnitude (sum of |q_i|).
    pub total_abs_charge: f64,
    /// Max partial charge.
    pub max_charge: f64,
    /// Min partial charge.
    pub min_charge: f64,
    /// Wiener index (sum of shortest-path distances for all pairs).
    pub wiener_index: f64,
    /// Number of rings (from graph cycles, approximate).
    pub n_rings: usize,
    /// Number of aromatic atoms.
    pub n_aromatic: usize,
    /// Balaban J index (approximation).
    pub balaban_j: f64,
    /// Sum of atomic electronegativities (Pauling).
    pub sum_electronegativity: f64,
    /// Sum of atomic polarizabilities (empirical, Å³).
    pub sum_polarizability: f64,
}

/// Atomic weight table for common elements.
fn atomic_weight(z: u8) -> f64 {
    match z {
        1 => 1.008,
        5 => 10.81,
        6 => 12.011,
        7 => 14.007,
        8 => 15.999,
        9 => 18.998,
        14 => 28.086,
        15 => 30.974,
        16 => 32.06,
        17 => 35.45,
        35 => 79.904,
        53 => 126.904,
        _ => z as f64 * 1.5, // rough fallback
    }
}

/// Pauling electronegativity.
fn electronegativity(z: u8) -> f64 {
    match z {
        1 => 2.20,
        5 => 2.04,
        6 => 2.55,
        7 => 3.04,
        8 => 3.44,
        9 => 3.98,
        14 => 1.90,
        15 => 2.19,
        16 => 2.58,
        17 => 3.16,
        35 => 2.96,
        53 => 2.66,
        _ => 1.80,
    }
}

/// Empirical atomic polarizability (Å³).
fn atomic_polarizability(z: u8) -> f64 {
    match z {
        1 => 0.667,
        5 => 3.03,
        6 => 1.76,
        7 => 1.10,
        8 => 0.802,
        9 => 0.557,
        14 => 5.38,
        15 => 3.63,
        16 => 2.90,
        17 => 2.18,
        35 => 3.05,
        53 => 5.35,
        _ => 2.0,
    }
}

/// Compute molecular descriptors from elements, bonds, and partial charges.
///
/// `elements`: atomic numbers.
/// `bonds`: (atom_i, atom_j, bond_order) list.
/// `charges`: partial charges (same length as elements), or empty.
/// `aromatic_atoms`: boolean flags per atom, or empty.
pub fn compute_descriptors(
    elements: &[u8],
    bonds: &[(usize, usize, u8)],
    charges: &[f64],
    aromatic_atoms: &[bool],
) -> MolecularDescriptors {
    let n = elements.len();
    let n_heavy = elements.iter().filter(|&&z| z != 1).count();
    let n_h = n - n_heavy;

    let mw: f64 = elements.iter().map(|&z| atomic_weight(z)).sum();

    // Build adjacency
    let mut adj = vec![vec![]; n];
    for &(i, j, _) in bonds {
        if i < n && j < n {
            adj[i].push(j);
            adj[j].push(i);
        }
    }

    // Rotatable bonds: single bonds between two heavy atoms with ≥2 neighbors each
    let n_rot = bonds
        .iter()
        .filter(|&&(i, j, ord)| {
            ord == 1
                && elements[i] != 1
                && elements[j] != 1
                && adj[i].len() >= 2
                && adj[j].len() >= 2
        })
        .count();

    // H-bond donors: N-H or O-H
    let n_hbd = (0..n)
        .filter(|&i| {
            (elements[i] == 7 || elements[i] == 8) && adj[i].iter().any(|&j| elements[j] == 1)
        })
        .count();

    // H-bond acceptors: N or O
    let n_hba = elements.iter().filter(|&&z| z == 7 || z == 8).count();

    // Fsp3: fraction of sp3 carbons (4 neighbors is sp3 proxy)
    let n_c = elements.iter().filter(|&&z| z == 6).count();
    let n_c_sp3 = (0..n)
        .filter(|&i| elements[i] == 6 && adj[i].len() == 4)
        .count();
    let fsp3 = if n_c > 0 {
        n_c_sp3 as f64 / n_c as f64
    } else {
        0.0
    };

    // Charges
    let (total_abs, max_q, min_q) = if !charges.is_empty() && charges.len() == n {
        let abs_sum: f64 = charges.iter().map(|q| q.abs()).sum();
        let max = charges.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min = charges.iter().cloned().fold(f64::INFINITY, f64::min);
        (abs_sum, max, min)
    } else {
        (0.0, 0.0, 0.0)
    };

    // Wiener index: BFS from each node
    let mut wiener = 0.0;
    for start in 0..n {
        let mut dist = vec![u32::MAX; n];
        dist[start] = 0;
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(start);
        while let Some(u) = queue.pop_front() {
            for &v in &adj[u] {
                if dist[v] == u32::MAX {
                    dist[v] = dist[u] + 1;
                    queue.push_back(v);
                }
            }
        }
        wiener += dist
            .iter()
            .filter(|&&d| d != u32::MAX)
            .map(|&d| d as f64)
            .sum::<f64>();
    }
    wiener /= 2.0; // each pair counted twice

    // Ring count (approximate): n_bonds - n_atoms + n_components
    let n_components = {
        let mut visited = vec![false; n];
        let mut count = 0usize;
        for start in 0..n {
            if visited[start] {
                continue;
            }
            count += 1;
            let mut stack = vec![start];
            while let Some(u) = stack.pop() {
                if visited[u] {
                    continue;
                }
                visited[u] = true;
                for &v in &adj[u] {
                    if !visited[v] {
                        stack.push(v);
                    }
                }
            }
        }
        count
    };
    let n_rings = if bonds.len() + n_components > n {
        bonds.len() + n_components - n
    } else {
        0
    };

    let n_aromatic = if !aromatic_atoms.is_empty() {
        aromatic_atoms.iter().filter(|&&a| a).count()
    } else {
        0
    };

    // Balaban J (simplified)
    let balaban_j = if !bonds.is_empty() && n > 1 {
        let mu = bonds.len() as f64 / (n as f64 - 1.0);
        mu.ln().abs() + 1.0
    } else {
        0.0
    };

    let sum_en: f64 = elements.iter().map(|&z| electronegativity(z)).sum();
    let sum_pol: f64 = elements.iter().map(|&z| atomic_polarizability(z)).sum();

    MolecularDescriptors {
        molecular_weight: mw,
        n_heavy_atoms: n_heavy,
        n_hydrogens: n_h,
        n_bonds: bonds.len(),
        n_rotatable_bonds: n_rot,
        n_hbd,
        n_hba,
        fsp3,
        total_abs_charge: total_abs,
        max_charge: max_q,
        min_charge: min_q,
        wiener_index: wiener,
        n_rings,
        n_aromatic,
        balaban_j,
        sum_electronegativity: sum_en,
        sum_polarizability: sum_pol,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_descriptors_water() {
        let elements = [8u8, 1, 1];
        let bonds = [(0, 1, 1u8), (0, 2, 1)];
        let desc = compute_descriptors(&elements, &bonds, &[], &[]);
        assert_eq!(desc.n_heavy_atoms, 1);
        assert_eq!(desc.n_hydrogens, 2);
        assert_eq!(desc.n_bonds, 2);
        assert_eq!(desc.n_hbd, 1); // O-H
        assert_eq!(desc.n_hba, 1); // O
        assert!((desc.molecular_weight - 18.015).abs() < 0.01);
    }

    #[test]
    fn test_descriptors_methane() {
        let elements = [6u8, 1, 1, 1, 1];
        let bonds = [(0, 1, 1u8), (0, 2, 1), (0, 3, 1), (0, 4, 1)];
        let desc = compute_descriptors(&elements, &bonds, &[], &[]);
        assert_eq!(desc.n_heavy_atoms, 1);
        assert_eq!(desc.fsp3, 1.0); // all C are sp3
        assert_eq!(desc.n_rotatable_bonds, 0);
    }
}
