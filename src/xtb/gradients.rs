//! Energy gradients for the public xTB energy (GFN2-xTB).
//!
//! [`crate::compute_xtb`] evaluates GFN2-xTB, so the gradient returned here is
//! the central-difference derivative of the GFN2 total energy: gradient and
//! energy describe one potential-energy surface by construction.
//!
//! Each displaced geometry is warm-started from the converged shell charges
//! and atomic multipoles of the reference point, so the SCC converges in a
//! few iterations instead of from a zero state.  The energy (and therefore the
//! gradient) is unchanged; the cost drops from `6·N` cold SCC solves to
//! `6·N` warm SCC solves.  A fully analytic GFN2 gradient (multipole
//! electrostatics and D4) remains future work.

use serde::{Deserialize, Serialize};

use super::gfn2::{solve_gfn2, solve_gfn2_warm};
use crate::fd_gradient::{central_difference_gradient, DEFAULT_STEP_ANGSTROM};

/// Result of xTB gradient computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XtbGradientResult {
    /// Energy gradient dE/dR per atom in eV/Å.
    pub gradients: Vec<[f64; 3]>,
    /// Total xTB (GFN2) energy in eV.
    pub energy: f64,
}

/// Compute the xTB (GFN2) energy gradient.
///
/// Returns `dE/dR` in eV/Å for each atom, where `E` is the GFN2-xTB total
/// energy (the one [`crate::compute_xtb`] reports), by central differences
/// with a 5e-4 Å step.
///
/// # Errors
///
/// Returns `Err` if the SCC fails at the reference geometry or at any
/// displaced geometry (invalid geometry, an element without parameters, or a
/// non-convergent SCC).
pub fn compute_xtb_gradient(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<XtbGradientResult, String> {
    let reference = solve_gfn2(elements, positions)?;
    let warm = reference.warm_start();
    let gradients = central_difference_gradient(positions, DEFAULT_STEP_ANGSTROM, |p| {
        solve_gfn2_warm(elements, p, &warm)
            .or_else(|_| solve_gfn2(elements, p))
            .map(|r| r.total_energy)
    })?;
    Ok(XtbGradientResult {
        gradients,
        energy: reference.total_energy,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER: [[f64; 3]; 3] = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];

    #[test]
    fn test_xtb_gradient_h2() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let result = compute_xtb_gradient(&elements, &positions).unwrap();
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
    fn test_xtb_gradient_water() {
        let elements = [8u8, 1, 1];
        let result = compute_xtb_gradient(&elements, &WATER).unwrap();
        assert!(result.gradients.iter().flatten().all(|v| v.is_finite()));
        for d in 0..3 {
            let sum: f64 = result.gradients.iter().map(|g| g[d]).sum();
            assert!(sum.abs() < 1e-3, "net force {sum:.2e} eV/Å");
        }
    }

    #[test]
    fn test_xtb_gradient_vs_numerical() {
        let elements = [8u8, 1, 1];
        let analytical = compute_xtb_gradient(&elements, &WATER).unwrap();
        let h = 1e-4;
        for a in 0..3 {
            for d in 0..3 {
                let mut p = WATER.to_vec();
                let mut m = WATER.to_vec();
                p[a][d] += h;
                m[a][d] -= h;
                let numerical = (solve_gfn2(&elements, &p).unwrap().total_energy
                    - solve_gfn2(&elements, &m).unwrap().total_energy)
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
