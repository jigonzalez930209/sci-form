//! EHT band structure calculation for periodic systems.
//!
//! Computes electronic band structure along high-symmetry k-paths
//! in the Brillouin zone for crystalline materials with periodic
//! boundary conditions.

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

/// A k-point in reciprocal space.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KPoint {
    /// Fractional coordinates in reciprocal space.
    pub frac: [f64; 3],
    /// Label (e.g., "Γ", "X", "M", "K").
    pub label: Option<String>,
    /// Linear path distance for plotting.
    pub path_distance: f64,
}

/// Result of a band structure calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BandStructure {
    /// k-points along the path.
    pub kpoints: Vec<KPoint>,
    /// Eigenvalues at each k-point: bands[k_idx][band_idx] in eV.
    pub bands: Vec<Vec<f64>>,
    /// Number of bands (= number of basis functions per cell).
    pub n_bands: usize,
    /// Number of k-points.
    pub n_kpoints: usize,
    /// Fermi energy estimate (eV).
    pub fermi_energy: f64,
    /// Direct band gap (eV), if any.
    pub direct_gap: Option<f64>,
    /// Indirect band gap (eV), if any.
    pub indirect_gap: Option<f64>,
    /// High-symmetry point labels and their k-indices.
    pub high_symmetry_points: Vec<(String, usize)>,
}

/// Configuration for band structure calculation.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct BandStructureConfig {
    /// Number of k-points between high-symmetry points.
    pub n_kpoints_per_segment: usize,
    /// High-symmetry path (pairs of labels).
    pub path: Vec<([f64; 3], String)>,
}

impl Default for BandStructureConfig {
    fn default() -> Self {
        // Default: Γ → X → M → Γ for cubic systems
        Self {
            n_kpoints_per_segment: 50,
            path: vec![
                ([0.0, 0.0, 0.0], "Γ".to_string()),
                ([0.5, 0.0, 0.0], "X".to_string()),
                ([0.5, 0.5, 0.0], "M".to_string()),
                ([0.0, 0.0, 0.0], "Γ".to_string()),
            ],
        }
    }
}

/// Compute electronic band structure using EHT with Bloch's theorem.
///
/// For a periodic system, the EHT Hamiltonian at k-point **k** is:
///   H(k) = Σ_R H(R) exp(i k·R)
///   S(k) = Σ_R S(R) exp(i k·R)
///
/// where R are lattice translation vectors.
pub fn compute_band_structure(
    elements: &[u8],
    positions: &[[f64; 3]],
    lattice: &[[f64; 3]; 3],
    config: &BandStructureConfig,
    n_electrons: usize,
) -> Result<BandStructure, String> {
    if elements.is_empty() {
        return Err("No atoms provided".to_string());
    }

    // Generate k-point path
    let kpoints = generate_kpath(&config.path, config.n_kpoints_per_segment);
    let n_kpts = kpoints.len();

    // Build real-space EHT matrices for the unit cell and nearest neighbors
    let eht_result = crate::eht::solve_eht(elements, positions, None)?;
    let n_basis = eht_result.energies.len();

    // At each k-point, solve the generalized eigenvalue problem
    let mut bands = Vec::with_capacity(n_kpts);
    let mut high_sym = Vec::new();

    for (k_idx, kpt) in kpoints.iter().enumerate() {
        // Build H(k) and S(k) using Bloch phase factors
        let (h_k, s_k) = build_bloch_matrices(elements, positions, lattice, &kpt.frac, n_basis);

        // Solve generalized eigenvalue problem H(k) C = S(k) C ε
        let eigenvalues = solve_generalized_eigen(&h_k, &s_k)?;
        bands.push(eigenvalues);

        if let Some(ref label) = kpt.label {
            high_sym.push((label.clone(), k_idx));
        }
    }

    // Estimate Fermi energy
    let n_occupied = n_electrons / 2;
    let fermi_energy = estimate_fermi_energy(&bands, n_occupied);

    // Compute band gaps
    let (direct_gap, indirect_gap) = compute_band_gaps(&bands, n_occupied);

    Ok(BandStructure {
        kpoints,
        bands,
        n_bands: n_basis,
        n_kpoints: n_kpts,
        fermi_energy,
        direct_gap,
        indirect_gap,
        high_symmetry_points: high_sym,
    })
}

/// Generate a k-point path through the Brillouin zone.
fn generate_kpath(path: &[([f64; 3], String)], n_per_segment: usize) -> Vec<KPoint> {
    let mut kpoints = Vec::new();
    let mut path_dist = 0.0;

    for i in 0..path.len() {
        let (k, label) = &path[i];

        if i == 0 {
            kpoints.push(KPoint {
                frac: *k,
                label: Some(label.clone()),
                path_distance: 0.0,
            });
            continue;
        }

        let (k_prev, _) = &path[i - 1];
        let dk = [k[0] - k_prev[0], k[1] - k_prev[1], k[2] - k_prev[2]];
        let seg_len = (dk[0] * dk[0] + dk[1] * dk[1] + dk[2] * dk[2]).sqrt();

        for j in 1..=n_per_segment {
            let t = j as f64 / n_per_segment as f64;
            let frac = [
                k_prev[0] + t * dk[0],
                k_prev[1] + t * dk[1],
                k_prev[2] + t * dk[2],
            ];
            let is_endpoint = j == n_per_segment;
            path_dist += seg_len / n_per_segment as f64;

            kpoints.push(KPoint {
                frac,
                label: if is_endpoint {
                    Some(label.clone())
                } else {
                    None
                },
                path_distance: path_dist,
            });
        }
    }

    kpoints
}

/// Build Bloch Hamiltonian and overlap matrices at k-point.
fn build_bloch_matrices(
    elements: &[u8],
    positions: &[[f64; 3]],
    lattice: &[[f64; 3]; 3],
    k: &[f64; 3],
    n_basis: usize,
) -> (DMatrix<f64>, DMatrix<f64>) {
    // For the R=0 image, use the standard EHT matrices
    let basis = crate::eht::basis::build_basis(elements, positions);
    let s_0 = crate::eht::overlap::build_overlap_matrix(&basis);
    let h_0 = crate::eht::hamiltonian::build_hamiltonian(&basis, &s_0, None);

    let n = n_basis.min(s_0.nrows());

    // Start with R=0 contribution (phase factor = 1 for k·0 = 0)
    let mut h_k = DMatrix::zeros(n, n);
    let mut s_k = DMatrix::zeros(n, n);

    for i in 0..n {
        for j in 0..n {
            h_k[(i, j)] = h_0[(i, j)];
            s_k[(i, j)] = s_0[(i, j)];
        }
    }

    // Add contributions from nearest-neighbor cells R = ±a, ±b, ±c
    let translations: Vec<[i32; 3]> = vec![
        [1, 0, 0],
        [-1, 0, 0],
        [0, 1, 0],
        [0, -1, 0],
        [0, 0, 1],
        [0, 0, -1],
    ];

    for r in &translations {
        let phase = 2.0
            * std::f64::consts::PI
            * (k[0] * r[0] as f64 + k[1] * r[1] as f64 + k[2] * r[2] as f64);
        let cos_phase = phase.cos();

        // Build translated positions
        let translated: Vec<[f64; 3]> = positions
            .iter()
            .map(|p| {
                [
                    p[0] + r[0] as f64 * lattice[0][0]
                        + r[1] as f64 * lattice[1][0]
                        + r[2] as f64 * lattice[2][0],
                    p[1] + r[0] as f64 * lattice[0][1]
                        + r[1] as f64 * lattice[1][1]
                        + r[2] as f64 * lattice[2][1],
                    p[2] + r[0] as f64 * lattice[0][2]
                        + r[1] as f64 * lattice[1][2]
                        + r[2] as f64 * lattice[2][2],
                ]
            })
            .collect();

        // Build inter-cell overlap using combined basis
        let basis_r = crate::eht::basis::build_basis(elements, &translated);
        // Compute cross-overlap S_{0R} by building overlap of combined [basis, basis_r]
        // and extracting the off-diagonal block
        let mut combined = basis.clone();
        combined.extend_from_slice(&basis_r);
        let s_combined = crate::eht::overlap::build_overlap_matrix(&combined);
        let s_r = s_combined.view((0, n), (n, basis_r.len())).clone_owned();
        let h_r = build_intercell_hamiltonian(&basis, &basis_r, &s_r);

        let nr = n.min(s_r.nrows()).min(s_r.ncols());
        for i in 0..nr {
            for j in 0..nr {
                h_k[(i, j)] += cos_phase * h_r[(i, j)];
                s_k[(i, j)] += cos_phase * s_r[(i, j)];
            }
        }
    }

    (h_k, s_k)
}

/// Build inter-cell Hamiltonian using Wolfsberg-Helmholz approximation.
fn build_intercell_hamiltonian(
    basis_0: &[crate::eht::basis::AtomicOrbital],
    basis_r: &[crate::eht::basis::AtomicOrbital],
    s_0r: &DMatrix<f64>,
) -> DMatrix<f64> {
    let n = basis_0.len().min(s_0r.nrows());
    let m = basis_r.len().min(s_0r.ncols());
    let mut h = DMatrix::zeros(n, m);
    let k_wh = 1.75; // Wolfsberg-Helmholz K parameter

    for i in 0..n {
        let hii = basis_0[i].vsip;
        for j in 0..m {
            let hjj = basis_r[j].vsip;
            h[(i, j)] = 0.5 * k_wh * (hii + hjj) * s_0r[(i, j)];
        }
    }

    h
}

/// Solve generalized eigenvalue problem H C = S C ε via Löwdin orthogonalization.
fn solve_generalized_eigen(h: &DMatrix<f64>, s: &DMatrix<f64>) -> Result<Vec<f64>, String> {
    let n = h.nrows();
    if n == 0 {
        return Ok(vec![]);
    }

    // S^{-1/2} via eigendecomposition
    let s_eigen = nalgebra::SymmetricEigen::new(s.clone());
    let mut s_inv_sqrt = DMatrix::zeros(n, n);

    for (i, &eval) in s_eigen.eigenvalues.iter().enumerate() {
        if eval > 1e-8 {
            let inv_sqrt = 1.0 / eval.sqrt();
            for j in 0..n {
                for k in 0..n {
                    s_inv_sqrt[(j, k)] +=
                        inv_sqrt * s_eigen.eigenvectors[(j, i)] * s_eigen.eigenvectors[(k, i)];
                }
            }
        }
    }

    // H' = S^{-1/2} H S^{-1/2}
    let h_prime = &s_inv_sqrt * h * &s_inv_sqrt;
    let eigen = nalgebra::SymmetricEigen::new(h_prime);

    let mut eigenvalues: Vec<f64> = eigen.eigenvalues.iter().copied().collect();
    eigenvalues.sort_by(|a, b| a.partial_cmp(b).unwrap());

    Ok(eigenvalues)
}

/// Estimate Fermi energy from band eigenvalues.
fn estimate_fermi_energy(bands: &[Vec<f64>], n_occupied: usize) -> f64 {
    if bands.is_empty() || n_occupied == 0 {
        return 0.0;
    }

    // Collect all occupied orbital energies
    let mut all_occupied: Vec<f64> = bands
        .iter()
        .filter_map(|eigenvals| {
            if eigenvals.len() > n_occupied {
                Some((eigenvals[n_occupied - 1] + eigenvals[n_occupied]) / 2.0)
            } else {
                eigenvals.last().copied()
            }
        })
        .collect();

    all_occupied.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if all_occupied.is_empty() {
        return 0.0;
    }
    all_occupied[all_occupied.len() / 2]
}

/// Compute direct and indirect band gaps.
fn compute_band_gaps(bands: &[Vec<f64>], n_occupied: usize) -> (Option<f64>, Option<f64>) {
    if bands.is_empty() || n_occupied == 0 {
        return (None, None);
    }

    let mut min_direct = f64::MAX;
    let mut max_vb = f64::MIN;
    let mut min_cb = f64::MAX;

    for eigenvals in bands {
        if eigenvals.len() <= n_occupied {
            continue;
        }
        let vb_top = eigenvals[n_occupied - 1];
        let cb_bottom = eigenvals[n_occupied];

        let gap = cb_bottom - vb_top;
        if gap < min_direct && gap > 0.0 {
            min_direct = gap;
        }

        if vb_top > max_vb {
            max_vb = vb_top;
        }
        if cb_bottom < min_cb {
            min_cb = cb_bottom;
        }
    }

    let direct = if min_direct < f64::MAX {
        Some(min_direct)
    } else {
        None
    };

    let indirect = if min_cb > max_vb {
        Some(min_cb - max_vb)
    } else {
        None
    };

    (direct, indirect)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_kpath() {
        let config = BandStructureConfig::default();
        let kpoints = generate_kpath(&config.path, 10);
        assert!(!kpoints.is_empty());
        // First point should be Γ
        assert_eq!(kpoints[0].label.as_deref(), Some("Γ"));
    }

    #[test]
    fn test_band_gaps() {
        let bands = vec![vec![-5.0, -3.0, 1.0, 3.0], vec![-4.5, -2.5, 1.5, 3.5]];
        let (direct, indirect) = compute_band_gaps(&bands, 2);
        assert!(direct.is_some());
        assert!(indirect.is_some());
        assert!(indirect.unwrap() > 0.0);
    }
}
