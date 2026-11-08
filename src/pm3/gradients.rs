//! Energy gradients for PM3.
//!
//! For the MOPAC-parity solver (`Z ≤ 17`) the gradient uses the
//! **stationary-density** (Pulay-free) form: at a converged SCF density the
//! derivative of the density vanishes, so differentiating
//! `E = ½ Tr[P (H + G(P))] + E_nuc` at fixed `P` gives the molecular gradient.
//! The remaining integral derivatives are evaluated by central differences of
//! the fixed-density energy, which rebuilds the integrals but never runs an
//! SCF.  The legacy NDDO path (`Z > 17`) keeps the full central-difference
//! gradient of [`solve_pm3`].
//!
//! Before v1.0 this module held a hand-written Hellmann-Feynman/Pulay
//! expression for the *legacy* simplified NDDO model; that model is not the
//! one `solve_pm3` evaluates, so it was removed.

use serde::{Deserialize, Serialize};

use super::solver::solve_pm3;
use crate::fd_gradient::{central_difference_gradient, DEFAULT_STEP_ANGSTROM};

/// Result of PM3 gradient computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Pm3GradientResult {
    /// Energy gradient dE/dR per atom in eV/Å.
    pub gradients: Vec<[f64; 3]>,
    /// Total PM3 energy in eV.
    pub energy: f64,
    /// Heat of formation in kcal/mol.
    pub heat_of_formation: f64,
}

/// Compute the PM3 energy gradient.
///
/// Returns `dE/dR` in eV/Å for each atom, where `E` is the total energy
/// reported by [`solve_pm3`]. The derivative is evaluated by central
/// differences (step 5e-4 Å) of that energy.
///
/// # Errors
///
/// Returns `Err` if the PM3 SCF fails at the reference geometry or at any
/// displaced geometry (invalid geometry, an element without PM3 parameters,
/// an open-shell system the MOPAC-parity solver does not handle, or a
/// non-convergent SCF).
pub fn compute_pm3_gradient(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<Pm3GradientResult, String> {
    let reference = solve_pm3(elements, positions)?;
    // Stationary-density path for the faithful MOPAC-parity solver: no SCF
    // per displaced coordinate, only the integral rebuild at fixed P.
    let density = super::faithful::solve(elements, positions)
        .ok()
        .map(|result| result.density);
    let gradients = if let Some(density) = density {
        central_difference_gradient(positions, DEFAULT_STEP_ANGSTROM, |p| {
            super::faithful::fixed_density_energy(elements, p, &density)
        })?
    } else {
        central_difference_gradient(positions, DEFAULT_STEP_ANGSTROM, |p| {
            solve_pm3(elements, p).map(|r| r.total_energy)
        })?
    };
    Ok(Pm3GradientResult {
        gradients,
        energy: reference.total_energy,
        heat_of_formation: reference.heat_of_formation,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER: [[f64; 3]; 3] = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.96], [0.93, 0.0, -0.24]];

    #[test]
    fn test_pm3_gradient_h2() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let result = compute_pm3_gradient(&elements, &positions).unwrap();
        assert_eq!(result.gradients.len(), 2);
        assert!(result.energy.is_finite());
        for d in 0..3 {
            assert!(
                (result.gradients[0][d] + result.gradients[1][d]).abs() < 1e-3,
                "Forces not equal and opposite: {:?}",
                result.gradients
            );
        }
    }

    #[test]
    fn test_pm3_gradient_water() {
        let elements = [8u8, 1, 1];
        let result = compute_pm3_gradient(&elements, &WATER).unwrap();
        assert_eq!(result.gradients.len(), 3);
        assert!(result.gradients.iter().flatten().all(|v| v.is_finite()));
        // Translational invariance: the net force vanishes.
        for d in 0..3 {
            let sum: f64 = result.gradients.iter().map(|g| g[d]).sum();
            assert!(sum.abs() < 1e-3, "net force {sum:.2e} eV/Å");
        }
    }

    /// The gradient must agree with an independent, much smaller step of the
    /// reported energy (the F80 gate's definition of "correct").
    #[test]
    fn test_pm3_gradient_vs_numerical() {
        let elements = [8u8, 1, 1];
        let analytical = compute_pm3_gradient(&elements, &WATER).unwrap();
        let h = 1e-4;
        for a in 0..3 {
            for d in 0..3 {
                let mut p = WATER.to_vec();
                let mut m = WATER.to_vec();
                p[a][d] += h;
                m[a][d] -= h;
                let numerical = (solve_pm3(&elements, &p).unwrap().total_energy
                    - solve_pm3(&elements, &m).unwrap().total_energy)
                    / (2.0 * h);
                let diff = (analytical.gradients[a][d] - numerical).abs();
                assert!(
                    diff < 1e-3,
                    "atom {a} dir {d}: gradient={:.6} numerical={:.6}",
                    analytical.gradients[a][d],
                    numerical
                );
            }
        }
    }
}
