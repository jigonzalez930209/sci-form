//! MBH Hessian construction and frequency computation — E9.2
//!
//! Projects the full Hessian onto the reduced MBH basis and
//! solves the generalized eigenvalue problem for frequencies.

use super::blocks::{build_projection_matrix, detect_rigid_blocks};

/// MBH configuration.
#[derive(Debug, Clone)]
pub struct MbhConfig {
    /// Finite difference step for numerical Hessian (Å).
    pub fd_step: f64,
}

impl Default for MbhConfig {
    fn default() -> Self {
        Self { fd_step: 0.005 }
    }
}

/// Result of MBH vibrational analysis.
#[derive(Debug, Clone)]
pub struct MbhResult {
    /// Frequencies in cm⁻¹ (sorted ascending, may include negatives for imaginary).
    pub frequencies: Vec<f64>,
    /// Number of rigid blocks.
    pub n_blocks: usize,
    /// Number of flexible atoms.
    pub n_flexible: usize,
    /// Reduced DOF count.
    pub n_dof_reduced: usize,
    /// Full Cartesian DOF count.
    pub n_dof_full: usize,
    /// Speedup factor (n_dof_full / n_dof_reduced).
    pub speedup: f64,
}

/// Atomic masses.
fn atomic_mass(z: u8) -> f64 {
    match z {
        1 => 1.008,
        5 => 10.81,
        6 => 12.011,
        7 => 14.007,
        8 => 15.999,
        9 => 18.998,
        14 => 28.086,
        15 => 30.974,
        16 => 32.065,
        17 => 35.453,
        35 => 79.904,
        53 => 126.904,
        _ => 12.011,
    }
}

/// Compute MBH frequencies given an energy function.
///
/// `energy_fn` takes positions (flat coords) and returns energy in kcal/mol.
/// `rings` contains (atom_indices, is_aromatic) tuples from SSSR.
pub fn compute_mbh_frequencies(
    elements: &[u8],
    positions: &[[f64; 3]],
    rings: &[(Vec<usize>, bool)],
    energy_fn: &dyn Fn(&[f64]) -> f64,
    config: &MbhConfig,
) -> MbhResult {
    let n_atoms = elements.len();

    // 1. Detect rigid blocks
    let decomp = detect_rigid_blocks(n_atoms, elements, positions, rings);

    // 2. Build projection matrix L
    let l_cols = build_projection_matrix(&decomp, positions, elements);
    let n_reduced = l_cols.len();

    if n_reduced == 0 {
        return MbhResult {
            frequencies: vec![],
            n_blocks: 0,
            n_flexible: 0,
            n_dof_reduced: 0,
            n_dof_full: 3 * n_atoms,
            speedup: 1.0,
        };
    }

    // 3. Build reduced Hessian by finite differences in reduced coordinates
    let n_full = 3 * n_atoms;
    let h = config.fd_step;

    // Reference coordinates (flat)
    let ref_coords: Vec<f64> = positions.iter().flat_map(|p| p.iter().cloned()).collect();

    // Build H_MBH[i][j] = d²E / dq_i dq_j via central differences
    let mut h_mbh = vec![vec![0.0; n_reduced]; n_reduced];

    for i in 0..n_reduced {
        // Forward displacement in direction i
        let mut coords_p: Vec<f64> = ref_coords.clone();
        let mut coords_m: Vec<f64> = ref_coords.clone();
        for k in 0..n_full {
            coords_p[k] += h * l_cols[i][k];
            coords_m[k] -= h * l_cols[i][k];
        }

        let e_p = energy_fn(&coords_p);
        let e_m = energy_fn(&coords_m);
        let e_0 = energy_fn(&ref_coords);

        // Diagonal: (E+ - 2E0 + E-) / h²
        h_mbh[i][i] = (e_p - 2.0 * e_0 + e_m) / (h * h);

        // Off-diagonal via (E++ - E+- - E-+ + E--) / (4h²)
        for j in (i + 1)..n_reduced {
            let mut cpp = ref_coords.clone();
            let mut cpm = ref_coords.clone();
            let mut cmp = ref_coords.clone();
            let mut cmm = ref_coords.clone();

            for k in 0..n_full {
                cpp[k] += h * l_cols[i][k] + h * l_cols[j][k];
                cpm[k] += h * l_cols[i][k] - h * l_cols[j][k];
                cmp[k] -= h * l_cols[i][k] - h * l_cols[j][k];
                cmm[k] -= h * l_cols[i][k] + h * l_cols[j][k];
            }

            let epp = energy_fn(&cpp);
            let epm = energy_fn(&cpm);
            let emp = energy_fn(&cmp);
            let emm = energy_fn(&cmm);

            h_mbh[i][j] = (epp - epm - emp + emm) / (4.0 * h * h);
            h_mbh[j][i] = h_mbh[i][j];
        }
    }

    // 4. Build mass-weighted Hessian
    // M_MBH = L^T * M * L (but L already has mass-weighting in its columns for blocks)
    // For simplicity in the MBH approach, we mass-weight the Hessian directly

    // Build mass vector
    let mut mass_diag = vec![0.0; n_full];
    for i in 0..n_atoms {
        let m = atomic_mass(elements[i]);
        mass_diag[i * 3] = m;
        mass_diag[i * 3 + 1] = m;
        mass_diag[i * 3 + 2] = m;
    }

    // M_red[i][j] = sum_k L[k][i] * mass[k] * L[k][j]
    let mut m_red = vec![vec![0.0; n_reduced]; n_reduced];
    for i in 0..n_reduced {
        for j in i..n_reduced {
            let mut val = 0.0;
            for k in 0..n_full {
                val += l_cols[i][k] * mass_diag[k] * l_cols[j][k];
            }
            m_red[i][j] = val;
            m_red[j][i] = val;
        }
    }

    // 5. Solve generalized eigenvalue problem H*c = λ*M*c
    // Use M^{-1/2} * H * M^{-1/2} standard form
    // First compute M^{-1/2} via eigendecomposition of M
    let m_mat = nalgebra::DMatrix::from_fn(n_reduced, n_reduced, |i, j| m_red[i][j]);
    let h_mat = nalgebra::DMatrix::from_fn(n_reduced, n_reduced, |i, j| h_mbh[i][j]);

    let m_eigen = m_mat.symmetric_eigen();
    let mut m_inv_sqrt = nalgebra::DMatrix::zeros(n_reduced, n_reduced);
    for i in 0..n_reduced {
        if m_eigen.eigenvalues[i] > 1e-10 {
            let val = 1.0 / m_eigen.eigenvalues[i].sqrt();
            for r in 0..n_reduced {
                for c in 0..n_reduced {
                    m_inv_sqrt[(r, c)] +=
                        val * m_eigen.eigenvectors[(r, i)] * m_eigen.eigenvectors[(c, i)];
                }
            }
        }
    }

    let h_standard = &m_inv_sqrt * &h_mat * &m_inv_sqrt;
    let eigen = h_standard.symmetric_eigen();

    // 6. Convert eigenvalues to cm⁻¹
    // λ in kcal/(mol·amu·Å²) → convert to cm⁻¹
    // Factor: sqrt(λ * 4184 / (6.022e23 * 1e-20)) / (2π * 2.998e10)
    // Simplified: ν(cm⁻¹) = 108.59 * sqrt(|λ|) * sign(λ)
    let conversion = 108.59; // sqrt(kcal/mol / (amu * Å²)) → cm⁻¹

    let mut frequencies: Vec<f64> = eigen
        .eigenvalues
        .iter()
        .map(|&e| {
            if e >= 0.0 {
                conversion * e.sqrt()
            } else {
                -conversion * (-e).sqrt() // imaginary → negative
            }
        })
        .collect();

    frequencies.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let speedup = if decomp.n_dof_reduced > 0 {
        decomp.n_dof_full as f64 / decomp.n_dof_reduced as f64
    } else {
        1.0
    };

    MbhResult {
        frequencies,
        n_blocks: decomp.blocks.len(),
        n_flexible: decomp.flexible_atoms.len(),
        n_dof_reduced: decomp.n_dof_reduced,
        n_dof_full: decomp.n_dof_full,
        speedup,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn harmonic_energy(coords: &[f64]) -> f64 {
        // Simple harmonic potential: sum of 0.5 * k * (r - r0)^2 for bonded pairs
        let k = 100.0; // kcal/(mol·Å²)
        let r0 = 1.4;
        let n = coords.len() / 3;
        let mut e = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = coords[i * 3] - coords[j * 3];
                let dy = coords[i * 3 + 1] - coords[j * 3 + 1];
                let dz = coords[i * 3 + 2] - coords[j * 3 + 2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                if r < 3.0 {
                    // only nearby atoms "bonded"
                    e += 0.5 * k * (r - r0).powi(2);
                }
            }
        }
        e
    }

    #[test]
    fn test_mbh_no_blocks() {
        let elements = vec![6, 1, 1];
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let rings: Vec<(Vec<usize>, bool)> = vec![];
        let config = MbhConfig::default();

        let result =
            compute_mbh_frequencies(&elements, &positions, &rings, &harmonic_energy, &config);
        assert_eq!(result.n_blocks, 0);
        assert_eq!(result.n_flexible, 3);
        assert_eq!(result.n_dof_reduced, 9);
        assert_eq!(result.frequencies.len(), 9);
    }

    #[test]
    fn test_mbh_with_ring() {
        let _n = 9; // 6 ring + 3 external
        let elements = vec![6u8; 9];
        let positions: Vec<[f64; 3]> = (0..9)
            .map(|i| {
                if i < 6 {
                    let a = i as f64 * std::f64::consts::PI / 3.0;
                    [1.4 * a.cos(), 1.4 * a.sin(), 0.0]
                } else {
                    [4.0 + (i - 6) as f64 * 1.5, 0.0, 0.0]
                }
            })
            .collect();

        let rings = vec![(vec![0, 1, 2, 3, 4, 5], true)];
        let result = compute_mbh_frequencies(
            &elements,
            &positions,
            &rings,
            &harmonic_energy,
            &MbhConfig::default(),
        );

        assert_eq!(result.n_blocks, 1);
        assert_eq!(result.n_flexible, 3);
        // 6 DOF for ring + 9 DOF for 3 flexible atoms = 15
        assert_eq!(result.n_dof_reduced, 15);
        // Full would be 27
        assert_eq!(result.n_dof_full, 27);
        assert!(
            result.speedup > 1.0,
            "Should show speedup: {}",
            result.speedup
        );
    }

    #[test]
    fn test_mbh_frequencies_finite() {
        let elements = vec![6u8; 6];
        let positions: Vec<[f64; 3]> = (0..6)
            .map(|i| {
                let a = i as f64 * std::f64::consts::PI / 3.0;
                [1.4 * a.cos(), 1.4 * a.sin(), 0.0]
            })
            .collect();

        let rings = vec![(vec![0, 1, 2, 3, 4, 5], true)];
        let result = compute_mbh_frequencies(
            &elements,
            &positions,
            &rings,
            &harmonic_energy,
            &MbhConfig::default(),
        );

        for &f in &result.frequencies {
            assert!(f.is_finite(), "Frequency should be finite: {}", f);
        }
    }
}
