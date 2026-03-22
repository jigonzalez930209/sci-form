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

/// 3D molecular shape descriptors computed from atomic coordinates.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Descriptors3D {
    /// Radius of gyration (Å).
    pub radius_of_gyration: f64,
    /// Asphericity (0 = sphere, 1 = rod).
    pub asphericity: f64,
    /// Eccentricity (0 = sphere).
    pub eccentricity: f64,
    /// Principal moments of inertia ratios (NPR1 = I1/I3, NPR2 = I2/I3).
    pub npr1: f64,
    pub npr2: f64,
    /// Sphericity (0 = asymmetric, 1 = sphere).
    pub sphericity: f64,
    /// Molecular span: max distance between any two atoms (Å).
    pub span: f64,
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

/// Compute 3D shape descriptors from atomic coordinates and masses.
///
/// `elements`: atomic numbers (used for mass weighting).
/// `positions`: 3D coordinates as flat [x0,y0,z0,x1,y1,z1,...].
pub fn compute_3d_descriptors(elements: &[u8], positions: &[f64]) -> Descriptors3D {
    let n = elements.len();
    if n == 0 || positions.len() < n * 3 {
        return Descriptors3D {
            radius_of_gyration: 0.0,
            asphericity: 0.0,
            eccentricity: 0.0,
            npr1: 0.0,
            npr2: 0.0,
            sphericity: 0.0,
            span: 0.0,
        };
    }

    let masses: Vec<f64> = elements.iter().map(|&z| atomic_weight(z)).collect();
    let total_mass: f64 = masses.iter().sum();

    // Center of mass
    let mut com = [0.0f64; 3];
    for i in 0..n {
        for k in 0..3 {
            com[k] += masses[i] * positions[i * 3 + k];
        }
    }
    for k in 0..3 {
        com[k] /= total_mass;
    }

    // Radius of gyration
    let mut rg2 = 0.0f64;
    for i in 0..n {
        let mut d2 = 0.0;
        for k in 0..3 {
            let d = positions[i * 3 + k] - com[k];
            d2 += d * d;
        }
        rg2 += masses[i] * d2;
    }
    rg2 /= total_mass;
    let rg = rg2.sqrt();

    // Gyration tensor (3x3 symmetric)
    let mut gt = [[0.0f64; 3]; 3];
    for i in 0..n {
        let r = [
            positions[i * 3] - com[0],
            positions[i * 3 + 1] - com[1],
            positions[i * 3 + 2] - com[2],
        ];
        for a in 0..3 {
            for b in 0..3 {
                gt[a][b] += masses[i] * r[a] * r[b];
            }
        }
    }
    for a in 0..3 {
        for b in 0..3 {
            gt[a][b] /= total_mass;
        }
    }

    // Eigenvalues of 3x3 symmetric matrix (analytical Cardano's method)
    let evals = eigenvalues_3x3_symmetric(&gt);
    let mut sorted = evals;
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let (l1, l2, l3) = (sorted[0].max(0.0), sorted[1].max(0.0), sorted[2].max(0.0));
    let sum_l = l1 + l2 + l3;

    let asphericity = if sum_l > 1e-14 {
        ((l1 - l2).powi(2) + (l2 - l3).powi(2) + (l1 - l3).powi(2)) / (2.0 * sum_l * sum_l)
    } else {
        0.0
    };

    let eccentricity = if l3 > 1e-14 {
        (1.0 - l1 / l3).sqrt().clamp(0.0, 1.0)
    } else {
        0.0
    };

    let (npr1, npr2) = if l3 > 1e-14 {
        (l1 / l3, l2 / l3)
    } else {
        (0.0, 0.0)
    };

    let sphericity = if l3 > 1e-14 {
        let prod = (l1 * l2 * l3).powf(1.0 / 3.0);
        (prod / l3).min(1.0)
    } else {
        0.0
    };

    // Molecular span: max pairwise distance
    let mut span = 0.0f64;
    for i in 0..n {
        for j in (i + 1)..n {
            let mut d2 = 0.0;
            for k in 0..3 {
                let d = positions[i * 3 + k] - positions[j * 3 + k];
                d2 += d * d;
            }
            span = span.max(d2.sqrt());
        }
    }

    Descriptors3D {
        radius_of_gyration: rg,
        asphericity,
        eccentricity,
        npr1,
        npr2,
        sphericity,
        span,
    }
}

/// Eigenvalues of a 3x3 symmetric matrix using Cardano's formula.
fn eigenvalues_3x3_symmetric(m: &[[f64; 3]; 3]) -> [f64; 3] {
    let p1 = m[0][1] * m[0][1] + m[0][2] * m[0][2] + m[1][2] * m[1][2];

    if p1.abs() < 1e-30 {
        // Already diagonal
        return [m[0][0], m[1][1], m[2][2]];
    }

    let q = (m[0][0] + m[1][1] + m[2][2]) / 3.0;
    let p2 = (m[0][0] - q).powi(2) + (m[1][1] - q).powi(2) + (m[2][2] - q).powi(2) + 2.0 * p1;
    let p = (p2 / 6.0).sqrt();

    // B = (1/p) * (A - q*I)
    let b = [
        [(m[0][0] - q) / p, m[0][1] / p, m[0][2] / p],
        [m[1][0] / p, (m[1][1] - q) / p, m[1][2] / p],
        [m[2][0] / p, m[2][1] / p, (m[2][2] - q) / p],
    ];

    let det_b = b[0][0] * (b[1][1] * b[2][2] - b[1][2] * b[2][1])
        - b[0][1] * (b[1][0] * b[2][2] - b[1][2] * b[2][0])
        + b[0][2] * (b[1][0] * b[2][1] - b[1][1] * b[2][0]);

    let r = det_b / 2.0;
    let r_clamped = r.clamp(-1.0, 1.0);
    let phi = r_clamped.acos() / 3.0;

    let eig1 = q + 2.0 * p * phi.cos();
    let eig3 = q + 2.0 * p * (phi + 2.0 * std::f64::consts::FRAC_PI_3).cos();
    let eig2 = 3.0 * q - eig1 - eig3;

    [eig1, eig2, eig3]
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
