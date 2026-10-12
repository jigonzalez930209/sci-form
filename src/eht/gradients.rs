//! EHT energy gradients for geometry optimization.
//!
//! Computes the gradient of the EHT geometry surface by central finite
//! differences of
//!
//!   E(R) = Σ_k n_k ε_k(R) + Σ_{i<j} Z_i Z_j / R_ij
//!
//! (occupied band energy plus nuclear repulsion, the latter in eV·Å with the
//! e²/4πε₀ constant). Without the repulsion term the surface has no short-range
//! wall and optimizations collapse onto each other; the repulsion is what makes
//! the minimum physical at the qualitative EHT level.

use serde::{Deserialize, Serialize};

/// Result of EHT gradient calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EhtGradient {
    /// Gradient vectors per atom: [∂E/∂x, ∂E/∂y, ∂E/∂z] in eV/Å.
    pub gradients: Vec<[f64; 3]>,
    /// RMS gradient magnitude (eV/Å).
    pub rms_gradient: f64,
    /// Maximum gradient component (eV/Å).
    pub max_gradient: f64,
    /// Total energy (eV).
    pub energy: f64,
}

/// Result of EHT geometry optimization.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EhtOptResult {
    /// Optimized atomic positions (Å).
    pub positions: Vec<[f64; 3]>,
    /// Final total energy (eV).
    pub energy: f64,
    /// Number of optimization steps taken.
    pub n_steps: usize,
    /// Whether the optimization converged.
    pub converged: bool,
    /// Final RMS gradient (eV/Å).
    pub rms_gradient: f64,
    /// Energy trajectory.
    pub energies: Vec<f64>,
}

/// Configuration for EHT geometry optimization.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EhtOptConfig {
    /// Maximum number of optimization steps.
    pub max_steps: usize,
    /// RMS gradient convergence threshold (eV/Å).
    pub grad_threshold: f64,
    /// Energy convergence threshold (eV).
    pub energy_threshold: f64,
    /// Initial step size for steepest descent (Å).
    pub step_size: f64,
}

impl Default for EhtOptConfig {
    fn default() -> Self {
        Self {
            max_steps: 200,
            grad_threshold: 0.01,
            energy_threshold: 1e-6,
            step_size: 0.05,
        }
    }
}

/// Compute EHT analytical gradients.
///
/// Uses the expression:
///   ∂E/∂R_A = Σ_μν P_μν ∂H_μν/∂R_A − Σ_μν W_μν ∂S_μν/∂R_A
///
/// Gradients of S and H are computed via numerical finite differences
/// on the basis function overlaps and Hamiltonian elements.
/// Coulomb constant e²/(4πε₀) in eV·Å.
const COULOMB_EV_ANGSTROM: f64 = 14.399_645;

/// Nuclear repulsion energy (eV) of the EHT geometry surface.
fn eht_nuclear_repulsion(elements: &[u8], positions: &[[f64; 3]]) -> f64 {
    let n = elements.len();
    let mut energy = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r > 1e-10 {
                energy += COULOMB_EV_ANGSTROM * (elements[i] as f64) * (elements[j] as f64) / r;
            }
        }
    }
    energy
}

/// Occupied-band energy (eV) from an EHT result (SOMO counted once).
fn eht_band_energy(eht: &crate::eht::EhtResult) -> f64 {
    let n_occ = eht.n_electrons.div_ceil(2); // include SOMO for odd electrons
    let is_odd = eht.n_electrons % 2 == 1;
    eht.energies
        .iter()
        .take(n_occ)
        .enumerate()
        .map(|(i, &e)| if is_odd && i == n_occ - 1 { e } else { 2.0 * e })
        .sum()
}

/// EHT geometry-surface energy (eV): occupied band + nuclear repulsion.
fn eht_surface_energy(elements: &[u8], positions: &[[f64; 3]]) -> Result<f64, String> {
    let eht = crate::eht::solve_eht(elements, positions, None)?;
    Ok(eht_band_energy(&eht) + eht_nuclear_repulsion(elements, positions))
}

pub fn compute_eht_gradient(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<EhtGradient, String> {
    let delta = 1e-5; // Finite difference step in Å
    let n_atoms = elements.len();

    // Reference energy (band + nuclear repulsion).
    let e0 = eht_surface_energy(elements, positions)?;

    let mut gradients = vec![[0.0f64; 3]; n_atoms];

    for atom in 0..n_atoms {
        for coord in 0..3 {
            let mut pos_plus = positions.to_vec();
            let mut pos_minus = positions.to_vec();
            pos_plus[atom][coord] += delta;
            pos_minus[atom][coord] -= delta;

            let e_plus = eht_surface_energy(elements, &pos_plus)?;
            let e_minus = eht_surface_energy(elements, &pos_minus)?;

            gradients[atom][coord] = (e_plus - e_minus) / (2.0 * delta);
        }
    }

    let rms = (gradients
        .iter()
        .flat_map(|g| g.iter())
        .map(|x| x * x)
        .sum::<f64>()
        / (3 * n_atoms) as f64)
        .sqrt();
    let max = gradients
        .iter()
        .flat_map(|g| g.iter())
        .map(|x| x.abs())
        .fold(0.0f64, |a, b| a.max(b));

    Ok(EhtGradient {
        gradients,
        rms_gradient: rms,
        max_gradient: max,
        energy: e0,
    })
}

/// Optimize molecular geometry using EHT gradients with steepest descent.
/// Optimize a molecular geometry on the EHT energy surface.
///
/// Uses the shared BFGS optimizer (same engine as the periodic materials
/// optimizer) over `compute_eht_gradient`. The EHT surface is qualitative:
/// bond lengths are not quantitative and polar molecules can dissociate
/// (there is no correlation/repulsion parametrization), so prefer the
/// force-field/ETKDG pipelines for production geometries.
pub fn optimize_geometry_eht(
    elements: &[u8],
    initial_positions: &[[f64; 3]],
    config: Option<EhtOptConfig>,
) -> Result<EhtOptResult, String> {
    let cfg = config.unwrap_or_default();

    // Shared BFGS optimizer with line search (same engine as the periodic
    // materials optimizer); the EHT gradient provides energy and forces.
    let energy_force = |els: &[u8], pos: &[[f64; 3]]| -> Result<(f64, Vec<[f64; 3]>), String> {
        let g = compute_eht_gradient(els, pos)?;
        let forces: Vec<[f64; 3]> = g.gradients.iter().map(|v| [-v[0], -v[1], -v[2]]).collect();
        Ok((g.energy, forces))
    };
    let fw_config = crate::materials::geometry_opt::FrameworkOptConfig {
        method: crate::materials::geometry_opt::OptMethod::Bfgs,
        max_iter: cfg.max_steps,
        force_tol: cfg.grad_threshold,
        energy_tol: cfg.energy_threshold,
        max_step: 0.2,
        fixed_atoms: vec![],
        lattice: None,
    };
    let result = crate::materials::geometry_opt::optimize_framework(
        elements,
        initial_positions,
        &energy_force,
        &fw_config,
    )?;

    let n = elements.len().max(1);
    let rms_gradient = (result
        .forces
        .iter()
        .map(|f| f.iter().map(|x| x * x).sum::<f64>())
        .sum::<f64>()
        / (3.0 * n as f64))
        .sqrt();

    Ok(EhtOptResult {
        positions: result.positions,
        energy: result.energy,
        n_steps: result.n_iterations,
        converged: result.converged,
        rms_gradient,
        energies: result.energy_history,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eht_gradient_h2() {
        let elements = vec![1u8, 1];
        let positions = vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let grad = compute_eht_gradient(&elements, &positions);
        assert!(grad.is_ok());
        let g = grad.unwrap();
        assert_eq!(g.gradients.len(), 2);
        // H2 at equilibrium should have small gradients
        assert!(g.rms_gradient.is_finite());
    }

    #[test]
    fn test_optimize_geometry_eht_uses_shared_bfgs() {
        // The shared BFGS engine must drive the EHT surface downhill from a
        // distorted start. The EHT surface itself is qualitative (see the
        // module docs: it can favour dissociation for polar molecules), so the
        // test checks mechanics, not bond lengths.
        let elements = [8u8, 1, 1];
        let distorted = [[0.0, 0.0, 0.0], [1.25, 0.20, 0.10], [-1.05, -0.30, -0.15]];
        let result = optimize_geometry_eht(&elements, &distorted, None).unwrap();
        assert!(result.n_steps >= 1);
        assert!(result.energy.is_finite());
        assert!(
            result.energies.len() >= 2,
            "energy trajectory should be recorded"
        );
        let first = result.energies[0];
        assert!(
            result.energy <= first,
            "BFGS should not end above the initial energy: {first} -> {}",
            result.energy
        );
        assert!(
            result.rms_gradient.is_finite() && result.rms_gradient >= 0.0,
            "rms gradient must be reported"
        );
    }
}
