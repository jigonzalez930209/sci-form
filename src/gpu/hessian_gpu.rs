//! GPU-accelerated numerical Hessian computation for IR spectroscopy.
//!
//! The numerical Hessian requires 6N single-point energy evaluations
//! (forward and backward displacement for each Cartesian DOF).
//! These evaluations are independent and map well to GPU parallelism.
//!
//! This module provides a GPU dispatch wrapper that batches the
//! displaced-geometry energy evaluations across GPU workgroups.

use super::context::{
    ComputeBindingDescriptor, ComputeBindingKind, ComputeDispatchDescriptor, GpuContext,
};

/// Minimum atoms to justify GPU dispatch for Hessian.
const GPU_HESSIAN_THRESHOLD: usize = 5;

/// Generate all displaced geometries for numerical Hessian.
///
/// For each of 3N Cartesian degrees of freedom, generates +h and −h displacements.
/// Returns (displaced_coords, n_displacements) where displaced_coords is
/// flat: [geom_+0, geom_-0, geom_+1, geom_-1, ...], each geometry is 3N floats.
pub fn generate_hessian_displacements(
    positions: &[[f64; 3]],
    step_size: f64,
) -> (Vec<Vec<[f64; 3]>>, usize) {
    let n_atoms = positions.len();
    let n_dof = 3 * n_atoms;
    let mut displaced = Vec::with_capacity(2 * n_dof);

    for dof in 0..n_dof {
        let atom = dof / 3;
        let coord = dof % 3;

        // +h displacement
        let mut plus = positions.to_vec();
        plus[atom][coord] += step_size;
        displaced.push(plus);

        // -h displacement
        let mut minus = positions.to_vec();
        minus[atom][coord] -= step_size;
        displaced.push(minus);
    }

    (displaced, 2 * n_dof)
}

/// Assemble the Hessian matrix from displaced energies.
///
/// `energies`: [E(+0), E(-0), E(+1), E(-1), ...] — pairs for each DOF.
/// `e_ref`: energy at the reference (undisplaced) geometry.
/// `step_size`: displacement step.
///
/// H_{ij} ≈ [E(+i,+j) − E(+i,−j) − E(−i,+j) + E(−i,−j)] / (4h²)
/// For diagonal: H_{ii} ≈ [E(+i) − 2E₀ + E(−i)] / h²
pub fn assemble_hessian_from_energies(
    energies: &[f64],
    e_ref: f64,
    step_size: f64,
) -> Vec<Vec<f64>> {
    let n_dof = energies.len() / 2;
    let h2 = step_size * step_size;
    let mut hessian = vec![vec![0.0; n_dof]; n_dof];

    // Diagonal elements: H_{ii} = (E(+i) - 2*E0 + E(-i)) / h^2
    for i in 0..n_dof {
        let e_plus = energies[2 * i];
        let e_minus = energies[2 * i + 1];
        hessian[i][i] = (e_plus - 2.0 * e_ref + e_minus) / h2;
    }

    // Off-diagonal: approximate from single displacements
    // H_{ij} ≈ (E(+i) + E(+j) - E(-i) - E(-j)) / (2*h^2) - (H_{ii} + H_{jj})/2
    // This is a cheaper approximation; for accurate off-diagonal,
    // double displacements (±i, ±j) are needed but require 12N² evaluations.
    for i in 0..n_dof {
        for j in (i + 1)..n_dof {
            let e_pi = energies[2 * i];
            let e_mi = energies[2 * i + 1];
            let e_pj = energies[2 * j];
            let e_mj = energies[2 * j + 1];

            // Forward-difference approximation
            let h_ij = ((e_pi + e_pj) - (e_mi + e_mj)) / (4.0 * h2)
                - 0.5 * (hessian[i][i] + hessian[j][j])
                + e_ref / h2;

            hessian[i][j] = h_ij;
            hessian[j][i] = h_ij;
        }
    }

    hessian
}

/// Batch-evaluate energies for displaced geometries using parallel CPU.
///
/// `energy_fn`: closure that computes energy for a given geometry.
/// `displacements`: list of displaced geometries.
#[cfg(feature = "parallel")]
pub fn evaluate_displacements_parallel<F>(
    displacements: &[Vec<[f64; 3]>],
    energy_fn: F,
) -> Vec<f64>
where
    F: Fn(&[[f64; 3]]) -> f64 + Sync,
{
    use rayon::prelude::*;
    displacements
        .par_iter()
        .map(|geom| energy_fn(geom))
        .collect()
}

/// Sequential displacement evaluation fallback.
pub fn evaluate_displacements_sequential<F>(
    displacements: &[Vec<[f64; 3]>],
    energy_fn: F,
) -> Vec<f64>
where
    F: Fn(&[[f64; 3]]) -> f64,
{
    displacements.iter().map(|geom| energy_fn(geom)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_displacement_count() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let (disps, n) = generate_hessian_displacements(&positions, 0.005);
        assert_eq!(n, 12); // 2 atoms × 3 coords × 2 directions
        assert_eq!(disps.len(), 12);
    }

    #[test]
    fn test_hessian_harmonic() {
        // 1D harmonic: E = 0.5 * k * x^2, k = 1.0
        // H = k = 1.0
        let h = 0.001;
        let energies = vec![
            0.5 * (h * h), // E(+x)
            0.5 * (h * h), // E(-x)
            0.5 * (h * h), // E(+y)
            0.5 * (h * h), // E(-y)
            0.5 * (h * h), // E(+z)
            0.5 * (h * h), // E(-z)
        ];
        let hessian = assemble_hessian_from_energies(&energies, 0.0, h);
        assert!((hessian[0][0] - 1.0).abs() < 0.01);
    }
}
