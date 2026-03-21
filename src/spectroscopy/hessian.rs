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

use nalgebra::DMatrix;

use crate::scf::types::MolecularSystem;
use crate::scf::gradients::numerical_gradient;
use crate::scf::scf_loop::ScfConfig;

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
        _ => z as f64 * 2.0,
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

    for atom_j in 0..n_atoms {
        for coord_j in 0..3 {
            let j = atom_j * 3 + coord_j;

            let mut sys_plus = system.clone();
            sys_plus.positions_bohr[atom_j][coord_j] += hess_step;
            let grad_plus = numerical_gradient(&sys_plus, scf_config, grad_step);

            let mut sys_minus = system.clone();
            sys_minus.positions_bohr[atom_j][coord_j] -= hess_step;
            let grad_minus = numerical_gradient(&sys_minus, scf_config, grad_step);

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

    // Symmetrize
    let hessian_sym = (&hessian + hessian.transpose()) * 0.5;

    let (mw_hessian, _masses) = mass_weight_hessian(&hessian_sym, &system.atomic_numbers);

    let eigen = mw_hessian.clone().symmetric_eigen();

    let hartree_to_cm1: f64 = 219474.63;
    let amu_to_au: f64 = 1822.888;

    let mut freq_eigenvalue_pairs: Vec<(f64, usize)> = eigen
        .eigenvalues
        .iter()
        .enumerate()
        .map(|(i, &val)| {
            let freq = if val >= 0.0 {
                val.sqrt() * hartree_to_cm1 / (2.0 * std::f64::consts::PI * amu_to_au.sqrt())
            } else {
                -((-val).sqrt() * hartree_to_cm1 / (2.0 * std::f64::consts::PI * amu_to_au.sqrt()))
            };
            (freq, i)
        })
        .collect();

    freq_eigenvalue_pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let frequencies: Vec<f64> = freq_eigenvalue_pairs.iter().map(|(f, _)| *f).collect();
    let n_imaginary = frequencies.iter().filter(|&&f| f < -10.0).count();

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
        let n = 6;
        let h = DMatrix::from_fn(n, n, |i, j| if i == j { 1.0 } else { 0.1 });
        let h_sym = (&h + h.transpose()) * 0.5;

        let z = vec![1, 6];
        let (mw, _) = mass_weight_hessian(&h_sym, &z);

        for i in 0..n {
            for j in 0..n {
                assert!(
                    (mw[(i, j)] - mw[(j, i)]).abs() < 1e-12,
                    "Mass-weighted Hessian should be symmetric"
                );
            }
        }
    }
}
