//! Analytical and numerical nuclear gradients ∂E/∂R.
//!
//! Implements the Hellmann-Feynman theorem with Pulay force corrections:
//!
//! dE/dR_A = Σ_μν P_μν ∂H_μν/∂R_A
//!         + Σ_μν (Σ_λσ P_μλ P_νσ - P_μν) ∂S_μν/∂R_A × ε_weighted
//!         + ∂V_nn/∂R_A
//!
//! For initial implementation, uses numerical (central finite differences)
//! gradients which are simpler and useful for validation.

use nalgebra::DMatrix;

use crate::experimental_2::types::MolecularSystem;

use super::scf_loop::{run_scf, ScfConfig};

/// Gradient result for all atoms.
#[derive(Debug, Clone)]
pub struct GradientResult {
    /// Nuclear gradients dE/dR in Hartree/Bohr, shape [n_atoms][3].
    pub gradients: Vec<[f64; 3]>,
    /// RMS gradient magnitude.
    pub rms_gradient: f64,
    /// Maximum gradient component.
    pub max_gradient: f64,
    /// Total energy at the evaluated geometry.
    pub energy: f64,
}

impl GradientResult {
    /// Compute gradient norm (RMS).
    pub fn norm(&self) -> f64 {
        self.rms_gradient
    }
}

/// Compute nuclear gradients using central finite differences.
///
/// This is the simplest and most robust approach. For each coordinate,
/// displaces by ±δ and evaluates the SCF energy:
///
///   ∂E/∂R_i ≈ [E(R_i + δ) - E(R_i - δ)] / (2δ)
///
/// Cost: 6N SCF calculations for N atoms. Suitable for small molecules
/// and for validating analytical gradients.
pub fn numerical_gradient(
    system: &MolecularSystem,
    config: &ScfConfig,
    step_size: f64,
) -> GradientResult {
    let n_atoms = system.n_atoms();
    let mut gradients = vec![[0.0; 3]; n_atoms];

    // Get reference energy
    let e_ref = run_scf(system, config).total_energy;

    for atom in 0..n_atoms {
        for coord in 0..3 {
            // Forward displacement
            let mut sys_plus = system.clone();
            sys_plus.positions_bohr[atom][coord] += step_size;
            let e_plus = run_scf(&sys_plus, config).total_energy;

            // Backward displacement
            let mut sys_minus = system.clone();
            sys_minus.positions_bohr[atom][coord] -= step_size;
            let e_minus = run_scf(&sys_minus, config).total_energy;

            // Central difference
            gradients[atom][coord] = (e_plus - e_minus) / (2.0 * step_size);
        }
    }

    // Compute statistics
    let mut sum_sq = 0.0;
    let mut max_g = 0.0f64;
    let n_components = (n_atoms * 3) as f64;

    for g in &gradients {
        for &gc in g {
            sum_sq += gc * gc;
            max_g = max_g.max(gc.abs());
        }
    }

    let rms_gradient = (sum_sq / n_components).sqrt();

    GradientResult {
        gradients,
        rms_gradient,
        max_gradient: max_g,
        energy: e_ref,
    }
}

/// Compute the energy-weighted density matrix W.
///
///   W_μν = 2 Σ_i^occ ε_i C_μi C_νi
///
/// Used in Pulay force computation: the derivative of the
/// orthogonality constraint contributes -Σ_μν W_μν ∂S_μν/∂R_A.
pub fn energy_weighted_density(
    coefficients: &DMatrix<f64>,
    orbital_energies: &[f64],
    n_occupied: usize,
) -> DMatrix<f64> {
    let n_basis = coefficients.nrows();
    let mut w = DMatrix::zeros(n_basis, n_basis);

    for k in 0..n_occupied {
        let eps_k = orbital_energies[k];
        for mu in 0..n_basis {
            for nu in 0..n_basis {
                w[(mu, nu)] += 2.0 * eps_k * coefficients[(mu, k)] * coefficients[(nu, k)];
            }
        }
    }

    w
}

/// Compute nuclear repulsion gradient for atom A.
///
///   ∂V_nn/∂R_A = Z_A Σ_{B≠A} Z_B (R_A - R_B) / |R_A - R_B|³
pub fn nuclear_repulsion_gradient(
    atomic_numbers: &[u8],
    positions: &[[f64; 3]],
    atom_a: usize,
) -> [f64; 3] {
    let z_a = atomic_numbers[atom_a] as f64;
    let r_a = positions[atom_a];
    let mut grad = [0.0; 3];

    for (b, &z_b_u8) in atomic_numbers.iter().enumerate() {
        if b == atom_a {
            continue;
        }
        let z_b = z_b_u8 as f64;
        let r_b = positions[b];

        let dx = r_a[0] - r_b[0];
        let dy = r_a[1] - r_b[1];
        let dz = r_a[2] - r_b[2];
        let r2 = dx * dx + dy * dy + dz * dz;
        let r = r2.sqrt();
        let r3 = r * r2;

        let factor = z_a * z_b / r3;
        grad[0] += factor * dx;
        grad[1] += factor * dy;
        grad[2] += factor * dz;
    }

    grad
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nuclear_repulsion_gradient_h2() {
        let z = vec![1u8, 1];
        let d = 1.4; // ~0.74 Å in bohr
        let pos = vec![[0.0, 0.0, 0.0], [d, 0.0, 0.0]];

        let g0 = nuclear_repulsion_gradient(&z, &pos, 0);
        let g1 = nuclear_repulsion_gradient(&z, &pos, 1);

        // Gradients should be equal and opposite (Newton's 3rd law)
        for i in 0..3 {
            assert!(
                (g0[i] + g1[i]).abs() < 1e-12,
                "Nuclear gradient not antisymmetric"
            );
        }

        // x-component should be negative for atom 0 (repelled toward -x)
        assert!(g0[0] < 0.0, "Atom 0 should be pushed in -x direction");
    }

    #[test]
    fn test_energy_weighted_density_symmetry() {
        let n = 3;
        let mut c = DMatrix::zeros(n, n);
        c[(0, 0)] = 0.7;
        c[(1, 0)] = 0.5;
        c[(2, 0)] = 0.3;
        c[(0, 1)] = -0.5;
        c[(1, 1)] = 0.7;
        c[(2, 1)] = 0.1;

        let energies = vec![-0.5, 0.3, 1.2];
        let w = energy_weighted_density(&c, &energies, 1);

        // W should be symmetric
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (w[(i, j)] - w[(j, i)]).abs() < 1e-12,
                    "Energy-weighted density not symmetric"
                );
            }
        }
    }
}
