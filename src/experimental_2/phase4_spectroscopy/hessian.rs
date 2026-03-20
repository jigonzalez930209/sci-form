//! Semi-numerical Hessian matrix computation.
//!
//! Computes the Hessian (second derivatives of energy w.r.t. nuclear coordinates)
//! using central finite differences of analytical or numerical gradients:
//!
//!   H_{ij} = ∂²E / (∂R_i ∂R_j) ≈ [g_i(R_j + δ) - g_i(R_j - δ)] / (2δ)
//!
//! where g_i = ∂E/∂R_i is the gradient component.
//!
//! From the mass-weighted Hessian, eigenvalue decomposition yields:
//! - Normal mode frequencies ω_k = sqrt(λ_k)
//! - Normal mode displacement vectors Q_k
//!
//! The first 5 (linear) or 6 (nonlinear) eigenvalues correspond to
//! translations and rotations and should be near zero.

use nalgebra::DMatrix;

use crate::experimental_2::types::MolecularSystem;
use crate::experimental_2::phase3_scf_engine::gradients::numerical_gradient;
use crate::experimental_2::phase3_scf_engine::scf_loop::ScfConfig;


/// Atomic masses in amu for common elements.
pub fn atomic_mass(z: u8) -> f64 {
    match z {
        1 => 1.00794,
        2 => 4.00260,
        3 => 6.941,
        4 => 9.01218,
        5 => 10.811,
        6 => 12.011,
        7 => 14.007,
        8 => 15.999,
        9 => 18.998,
        10 => 20.180,
        11 => 22.990,
        12 => 24.305,
        13 => 26.982,
        14 => 28.086,
        15 => 30.974,
        16 => 32.065,
        17 => 35.453,
        35 => 79.904,
        53 => 126.904,
        _ => z as f64 * 2.0, // rough estimate
    }
}

/// Result of Hessian calculation.
#[derive(Debug, Clone)]
pub struct HessianResult {
    /// Cartesian Hessian matrix (3N × 3N) in Hartree/Bohr².
    pub hessian: DMatrix<f64>,
    /// Mass-weighted Hessian in Hartree/(amu·Bohr²).
    pub mass_weighted_hessian: DMatrix<f64>,
    /// Normal mode frequencies in cm⁻¹.
    pub frequencies: Vec<f64>,
    /// Normal mode eigenvectors (columns of 3N × 3N matrix).
    pub normal_modes: DMatrix<f64>,
    /// Number of imaginary frequencies.
    pub n_imaginary: usize,
}

/// Build the Cartesian Hessian by semi-numerical differentiation.
///
/// For each coordinate j, displaces by ±δ and computes the gradient,
/// then uses central differences:
///
///   H_{ij} = [g_i(R_j + δ) - g_i(R_j - δ)] / (2δ)
///
/// Cost: 6N gradient evaluations for N atoms.
pub fn compute_hessian(
    system: &MolecularSystem,
    scf_config: &ScfConfig,
    grad_step: f64,
    hess_step: f64,
) -> HessianResult {
    let n_atoms = system.n_atoms();
    let n_coords = n_atoms * 3;

    let mut hessian = DMatrix::zeros(n_coords, n_coords);

    // For each coordinate j, displace and compute gradient
    for atom_j in 0..n_atoms {
        for coord_j in 0..3 {
            let j = atom_j * 3 + coord_j;

            // Forward displacement
            let mut sys_plus = system.clone();
            sys_plus.positions_bohr[atom_j][coord_j] += hess_step;
            let grad_plus = numerical_gradient(&sys_plus, scf_config, grad_step);

            // Backward displacement
            let mut sys_minus = system.clone();
            sys_minus.positions_bohr[atom_j][coord_j] -= hess_step;
            let grad_minus = numerical_gradient(&sys_minus, scf_config, grad_step);

            // Central difference for all i components
            for atom_i in 0..n_atoms {
                for coord_i in 0..3 {
                    let i = atom_i * 3 + coord_i;
                    hessian[(i, j)] = (grad_plus.gradients[atom_i][coord_i]
                        - grad_minus.gradients[atom_i][coord_i])
                        / (2.0 * hess_step);
                }
            }
        }
    }

    // Symmetrize: H = (H + H^T) / 2
    let hessian_sym = (&hessian + hessian.transpose()) * 0.5;

    // Mass-weight the Hessian
    let (mw_hessian, _masses) = mass_weight_hessian(&hessian_sym, &system.atomic_numbers);

    // Diagonalize mass-weighted Hessian
    let eigen = mw_hessian.clone().symmetric_eigen();

    // Convert eigenvalues to frequencies in cm⁻¹
    // ω = sqrt(λ) in atomic units, then convert
    // 1 Hartree/(amu·Bohr²) → frequency conversion
    let hartree_to_cm1: f64 = 219474.63; // 1 Hartree in cm⁻¹ (for E_h / ℏ)
    let amu_to_au: f64 = 1822.888; // 1 amu in atomic mass unit (m_e)

    let mut freq_eigenvalue_pairs: Vec<(f64, usize)> = eigen
        .eigenvalues
        .iter()
        .enumerate()
        .map(|(i, &val)| {
            let freq = if val >= 0.0 {
                val.sqrt() * hartree_to_cm1 / (2.0 * std::f64::consts::PI * amu_to_au.sqrt())
            } else {
                // Imaginary frequency (transition state)
                -((-val).sqrt() * hartree_to_cm1 / (2.0 * std::f64::consts::PI * amu_to_au.sqrt()))
            };
            (freq, i)
        })
        .collect();

    freq_eigenvalue_pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let frequencies: Vec<f64> = freq_eigenvalue_pairs.iter().map(|(f, _)| *f).collect();
    let n_imaginary = frequencies.iter().filter(|&&f| f < -10.0).count();

    // Sort normal modes to match sorted frequencies
    let mut normal_modes = DMatrix::zeros(n_coords, n_coords);
    for (new_idx, &(_, old_idx)) in freq_eigenvalue_pairs.iter().enumerate() {
        for i in 0..n_coords {
            normal_modes[(i, new_idx)] = eigen.eigenvectors[(i, old_idx)];
        }
    }

    HessianResult {
        hessian: hessian_sym,
        mass_weighted_hessian: mw_hessian,
        frequencies,
        normal_modes,
        n_imaginary,
    }
}

/// Mass-weight the Hessian: H'_{ij} = H_{ij} / sqrt(m_i · m_j)
fn mass_weight_hessian(
    hessian: &DMatrix<f64>,
    atomic_numbers: &[u8],
) -> (DMatrix<f64>, Vec<f64>) {
    let n_atoms = atomic_numbers.len();
    let n_coords = n_atoms * 3;

    // Build mass vector (each atom contributes 3 identical entries)
    let masses: Vec<f64> = atomic_numbers
        .iter()
        .flat_map(|&z| {
            let m = atomic_mass(z);
            vec![m, m, m]
        })
        .collect();

    let mut mw = DMatrix::zeros(n_coords, n_coords);
    for i in 0..n_coords {
        for j in 0..n_coords {
            mw[(i, j)] = hessian[(i, j)] / (masses[i] * masses[j]).sqrt();
        }
    }

    (mw, masses)
}

/// Project out translations and rotations from the Hessian.
///
/// For N atoms, there are 3 translational and 3 rotational (or 2 for linear)
/// degrees of freedom that should be projected out before computing
/// vibrational frequencies.
pub fn project_tr(
    hessian: &DMatrix<f64>,
    positions: &[[f64; 3]],
    masses: &[f64],
) -> DMatrix<f64> {
    let n_coords = hessian.nrows();
    let n_atoms = n_coords / 3;

    // Build translation vectors
    let total_mass: f64 = (0..n_atoms).map(|i| masses[i * 3]).sum();
    let mut tr_vectors: Vec<Vec<f64>> = Vec::new();

    // 3 translations
    for axis in 0..3 {
        let mut v = vec![0.0; n_coords];
        for atom in 0..n_atoms {
            let m_sqrt = masses[atom * 3].sqrt();
            v[atom * 3 + axis] = m_sqrt / total_mass.sqrt();
        }
        tr_vectors.push(v);
    }

    // 3 rotations (Eckart conditions)
    let com: [f64; 3] = {
        let mut c = [0.0; 3];
        for atom in 0..n_atoms {
            let m = masses[atom * 3];
            for k in 0..3 {
                c[k] += m * positions[atom][k];
            }
        }
        for k in 0..3 {
            c[k] /= total_mass;
        }
        c
    };

    // Rotation around x: y*dz - z*dy
    let mut rx = vec![0.0; n_coords];
    let mut ry = vec![0.0; n_coords];
    let mut rz = vec![0.0; n_coords];
    for atom in 0..n_atoms {
        let m_sqrt = masses[atom * 3].sqrt();
        let x = positions[atom][0] - com[0];
        let y = positions[atom][1] - com[1];
        let z = positions[atom][2] - com[2];

        // Rx: (0, -z, y)
        rx[atom * 3 + 1] = -z * m_sqrt;
        rx[atom * 3 + 2] = y * m_sqrt;

        // Ry: (z, 0, -x)
        ry[atom * 3] = z * m_sqrt;
        ry[atom * 3 + 2] = -x * m_sqrt;

        // Rz: (-y, x, 0)
        rz[atom * 3] = -y * m_sqrt;
        rz[atom * 3 + 1] = x * m_sqrt;
    }
    tr_vectors.push(rx);
    tr_vectors.push(ry);
    tr_vectors.push(rz);

    // Gram-Schmidt orthonormalization of TR vectors
    let mut ortho: Vec<Vec<f64>> = Vec::new();
    for v in &tr_vectors {
        let mut u = v.clone();
        for prev in &ortho {
            let dot: f64 = u.iter().zip(prev.iter()).map(|(a, b)| a * b).sum();
            for (i, p) in prev.iter().enumerate() {
                u[i] -= dot * p;
            }
        }
        let norm: f64 = u.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm > 1e-8 {
            for x in &mut u {
                *x /= norm;
            }
            ortho.push(u);
        }
    }

    // Build projector P = I - Σ_k |v_k><v_k|
    let mut proj = DMatrix::identity(n_coords, n_coords);
    for v in &ortho {
        for i in 0..n_coords {
            for j in 0..n_coords {
                proj[(i, j)] -= v[i] * v[j];
            }
        }
    }

    // Project: H' = P^T H P
    &proj.transpose() * hessian * &proj
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atomic_masses() {
        assert!((atomic_mass(1) - 1.00794).abs() < 1e-4);
        assert!((atomic_mass(6) - 12.011).abs() < 1e-2);
        assert!((atomic_mass(8) - 15.999).abs() < 1e-2);
    }

    #[test]
    fn test_mass_weighting_symmetry() {
        let n = 6; // 2 atoms
        let h = DMatrix::from_fn(n, n, |i, j| if i == j { 1.0 } else { 0.1 });
        let h_sym = (&h + h.transpose()) * 0.5;

        let z = vec![1, 6]; // H, C
        let (mw, _) = mass_weight_hessian(&h_sym, &z);

        // Mass-weighted Hessian should be symmetric
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (mw[(i, j)] - mw[(j, i)]).abs() < 1e-12,
                    "Mass-weighted Hessian should be symmetric"
                );
            }
        }
    }

    #[test]
    fn test_projection_matrix_is_idempotent() {
        let n = 9; // 3 atoms
        let pos = vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.0, 1.5, 0.0]];
        let masses: Vec<f64> = vec![12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 16.0, 16.0, 16.0];

        let h = DMatrix::identity(n, n);
        let h_proj = project_tr(&h, &pos, &masses);

        // P·P should equal P (idempotent) within the projected subspace
        let h_proj2 = project_tr(&h_proj, &pos, &masses);
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (h_proj[(i, j)] - h_proj2[(i, j)]).abs() < 1e-10,
                    "Projection should be idempotent"
                );
            }
        }
    }
}
