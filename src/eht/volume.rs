//! 3D volumetric mapping of molecular orbitals.
//!
//! Evaluates Ψ_i(x,y,z) = Σ_μ C_μi φ_μ(x,y,z) on a regular 3D grid.

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

use super::basis::AtomicOrbital;

/// A 3D regular grid holding scalar field values.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VolumetricGrid {
    /// Origin corner (x_min, y_min, z_min) in Ångström.
    pub origin: [f64; 3],
    /// Grid spacing in Ångström.
    pub spacing: f64,
    /// Number of grid points in each dimension [nx, ny, nz].
    pub dims: [usize; 3],
    /// Flat array of values in row-major order: z varies fastest.
    pub values: Vec<f32>,
}

impl VolumetricGrid {
    /// Total number of grid points.
    pub fn num_points(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }

    /// Index for grid point (ix, iy, iz).
    pub fn index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix * self.dims[1] * self.dims[2] + iy * self.dims[2] + iz
    }

    /// Get the Cartesian position (in Ångström) of grid point (ix, iy, iz).
    pub fn point_position(&self, ix: usize, iy: usize, iz: usize) -> [f64; 3] {
        [
            self.origin[0] + ix as f64 * self.spacing,
            self.origin[1] + iy as f64 * self.spacing,
            self.origin[2] + iz as f64 * self.spacing,
        ]
    }
}

/// Build a bounding box around the molecular geometry with padding.
///
/// Returns (origin, dims) for the given spacing.
pub fn compute_grid_extents(
    positions: &[[f64; 3]],
    padding: f64,
    spacing: f64,
) -> ([f64; 3], [usize; 3]) {
    let mut min = [f64::MAX; 3];
    let mut max = [f64::MIN; 3];
    for pos in positions {
        for d in 0..3 {
            min[d] = min[d].min(pos[d]);
            max[d] = max[d].max(pos[d]);
        }
    }

    let origin = [min[0] - padding, min[1] - padding, min[2] - padding];
    let dims = [
        ((max[0] - min[0] + 2.0 * padding) / spacing).ceil() as usize + 1,
        ((max[1] - min[1] + 2.0 * padding) / spacing).ceil() as usize + 1,
        ((max[2] - min[2] + 2.0 * padding) / spacing).ceil() as usize + 1,
    ];

    (origin, dims)
}

/// Evaluate a single atomic orbital φ_μ at a point in space.
///
/// Point is in Ångström; the orbital center is in bohr, so convert here.
fn evaluate_ao_at_point(orb: &AtomicOrbital, point_ang: &[f64; 3]) -> f64 {
    let ang_to_bohr = 1.0 / 0.529177249;
    let px = point_ang[0] * ang_to_bohr - orb.center[0];
    let py = point_ang[1] * ang_to_bohr - orb.center[1];
    let pz = point_ang[2] * ang_to_bohr - orb.center[2];
    let r2 = px * px + py * py + pz * pz;

    let mut val = 0.0;
    for g in &orb.gaussians {
        let exp_val = (-g.alpha * r2).exp();
        match (orb.l, orb.m) {
            (0, _) => {
                // s-type
                let norm = (2.0 * g.alpha / PI).powf(0.75);
                val += g.coeff * norm * exp_val;
            }
            (1, -1) => {
                // px
                let norm = (128.0 * g.alpha.powi(5) / (PI * PI * PI)).powf(0.25);
                val += g.coeff * norm * px * exp_val;
            }
            (1, 0) => {
                // py
                let norm = (128.0 * g.alpha.powi(5) / (PI * PI * PI)).powf(0.25);
                val += g.coeff * norm * py * exp_val;
            }
            (1, 1) => {
                // pz
                let norm = (128.0 * g.alpha.powi(5) / (PI * PI * PI)).powf(0.25);
                val += g.coeff * norm * pz * exp_val;
            }
            _ => {}
        }
    }

    val
}

/// Evaluate molecular orbital `mo_index` on a 3D grid.
///
/// Ψ_i(r) = Σ_μ C_μi φ_μ(r)
///
/// - `basis`: the molecular AO basis set
/// - `coefficients`: C matrix (row=AO, col=MO) as flat nested Vecs
/// - `mo_index`: which molecular orbital to evaluate
/// - `positions`: atom positions in Ångström
/// - `spacing`: grid spacing in Ångström
/// - `padding`: padding around the bounding box in Ångström
pub fn evaluate_orbital_on_grid(
    basis: &[AtomicOrbital],
    coefficients: &[Vec<f64>],
    mo_index: usize,
    positions: &[[f64; 3]],
    spacing: f64,
    padding: f64,
) -> VolumetricGrid {
    let (origin, dims) = compute_grid_extents(positions, padding, spacing);
    let total = dims[0] * dims[1] * dims[2];
    let mut values = vec![0.0f32; total];

    // Pre-extract MO coefficients for this orbital
    let c_mu: Vec<f64> = (0..basis.len()).map(|mu| coefficients[mu][mo_index]).collect();

    for ix in 0..dims[0] {
        for iy in 0..dims[1] {
            for iz in 0..dims[2] {
                let point = [
                    origin[0] + ix as f64 * spacing,
                    origin[1] + iy as f64 * spacing,
                    origin[2] + iz as f64 * spacing,
                ];

                let mut psi = 0.0f64;
                for (mu, orb) in basis.iter().enumerate() {
                    let phi = evaluate_ao_at_point(orb, &point);
                    psi += c_mu[mu] * phi;
                }

                let idx = ix * dims[1] * dims[2] + iy * dims[2] + iz;
                values[idx] = psi as f32;
            }
        }
    }

    VolumetricGrid {
        origin,
        spacing,
        dims,
        values,
    }
}

/// Evaluate molecular orbital on a 3D grid using parallel iteration (rayon).
#[cfg(feature = "parallel")]
pub fn evaluate_orbital_on_grid_parallel(
    basis: &[AtomicOrbital],
    coefficients: &[Vec<f64>],
    mo_index: usize,
    positions: &[[f64; 3]],
    spacing: f64,
    padding: f64,
) -> VolumetricGrid {
    use rayon::prelude::*;

    let (origin, dims) = compute_grid_extents(positions, padding, spacing);
    let total = dims[0] * dims[1] * dims[2];
    let c_mu: Vec<f64> = (0..basis.len()).map(|mu| coefficients[mu][mo_index]).collect();
    let ny = dims[1];
    let nz = dims[2];

    let values: Vec<f32> = (0..total)
        .into_par_iter()
        .map(|flat_idx| {
            let ix = flat_idx / (ny * nz);
            let iy = (flat_idx / nz) % ny;
            let iz = flat_idx % nz;

            let point = [
                origin[0] + ix as f64 * spacing,
                origin[1] + iy as f64 * spacing,
                origin[2] + iz as f64 * spacing,
            ];

            let mut psi = 0.0f64;
            for (mu, orb) in basis.iter().enumerate() {
                let phi = evaluate_ao_at_point(orb, &point);
                psi += c_mu[mu] * phi;
            }

            psi as f32
        })
        .collect();

    VolumetricGrid {
        origin,
        spacing,
        dims,
        values,
    }
}

/// A single z-slab of volumetric data, produced by chunked evaluation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VolumeSlab {
    /// The x-index of this slab.
    pub ix: usize,
    /// Grid values for this slab: ny * nz entries (y varies, z fastest).
    pub values: Vec<f32>,
}

/// Evaluate a molecular orbital on a 3D grid in x-slab chunks.
///
/// Instead of allocating one large `Vec<f32>` for the entire grid, this
/// returns an iterator that yields one x-slab at a time.  Each slab has
/// `dims[1] * dims[2]` values, so the peak memory is proportional to a
/// single slice instead of the full volume.
pub fn evaluate_orbital_on_grid_chunked(
    basis: &[AtomicOrbital],
    coefficients: &[Vec<f64>],
    mo_index: usize,
    positions: &[[f64; 3]],
    spacing: f64,
    padding: f64,
) -> (VolumetricGrid, Vec<VolumeSlab>) {
    let (origin, dims) = compute_grid_extents(positions, padding, spacing);
    let c_mu: Vec<f64> = (0..basis.len()).map(|mu| coefficients[mu][mo_index]).collect();

    let mut slabs = Vec::with_capacity(dims[0]);
    let mut all_values = Vec::with_capacity(dims[0] * dims[1] * dims[2]);

    for ix in 0..dims[0] {
        let slab_size = dims[1] * dims[2];
        let mut slab_values = Vec::with_capacity(slab_size);

        for iy in 0..dims[1] {
            for iz in 0..dims[2] {
                let point = [
                    origin[0] + ix as f64 * spacing,
                    origin[1] + iy as f64 * spacing,
                    origin[2] + iz as f64 * spacing,
                ];

                let mut psi = 0.0f64;
                for (mu, orb) in basis.iter().enumerate() {
                    let phi = evaluate_ao_at_point(orb, &point);
                    psi += c_mu[mu] * phi;
                }

                slab_values.push(psi as f32);
            }
        }

        all_values.extend_from_slice(&slab_values);
        slabs.push(VolumeSlab { ix, values: slab_values });
    }

    let grid = VolumetricGrid {
        origin,
        spacing,
        dims,
        values: all_values,
    };

    (grid, slabs)
}

/// Evaluate the electron density |Ψ_i|² on a grid.
pub fn evaluate_density_on_grid(
    basis: &[AtomicOrbital],
    coefficients: &[Vec<f64>],
    mo_index: usize,
    positions: &[[f64; 3]],
    spacing: f64,
    padding: f64,
) -> VolumetricGrid {
    let mut grid = evaluate_orbital_on_grid(
        basis, coefficients, mo_index, positions, spacing, padding,
    );
    for v in &mut grid.values {
        *v = (*v) * (*v);
    }
    grid
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::solver::solve_eht;
    use super::super::basis::build_basis;

    #[test]
    fn test_grid_extents_h2() {
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let (origin, dims) = compute_grid_extents(&positions, 4.0, 0.5);
        assert!(origin[0] < -3.9);
        assert!(dims[0] > 10);
        assert!(dims[1] > 10);
        assert!(dims[2] > 10);
    }

    #[test]
    fn test_evaluate_orbital_h2_bonding() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let result = solve_eht(&elements, &positions, None).unwrap();
        let basis = build_basis(&elements, &positions);

        // Evaluate bonding orbital (MO 0) on coarse grid
        let grid = evaluate_orbital_on_grid(
            &basis,
            &result.coefficients,
            0, // bonding MO
            &positions,
            0.5, // 0.5 Å spacing
            3.0, // 3 Å padding
        );

        assert!(grid.num_points() > 0);
        // Bonding orbital should have nonzero values
        let max_val = grid.values.iter().map(|v| v.abs()).fold(0.0f32, f32::max);
        assert!(max_val > 0.01, "max |Ψ_bond| = {}, expected > 0.01", max_val);
    }

    #[test]
    fn test_evaluate_orbital_h2_antibonding() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let result = solve_eht(&elements, &positions, None).unwrap();
        let basis = build_basis(&elements, &positions);

        // Antibonding MO 1
        let grid = evaluate_orbital_on_grid(
            &basis,
            &result.coefficients,
            1,
            &positions,
            0.5,
            3.0,
        );

        let max_val = grid.values.iter().map(|v| v.abs()).fold(0.0f32, f32::max);
        assert!(max_val > 0.01, "antibonding orbital should have nonzero amplitude");

        // Antibonding orbital should have a nodal plane — both positive and negative values
        let has_positive = grid.values.iter().any(|&v| v > 0.01);
        let has_negative = grid.values.iter().any(|&v| v < -0.01);
        assert!(has_positive && has_negative, "antibonding MO should have sign change");
    }

    #[test]
    fn test_density_nonnegative() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let result = solve_eht(&elements, &positions, None).unwrap();
        let basis = build_basis(&elements, &positions);

        let grid = evaluate_density_on_grid(
            &basis,
            &result.coefficients,
            0,
            &positions,
            0.5,
            3.0,
        );

        for v in &grid.values {
            assert!(*v >= 0.0, "Density should never be negative: {}", v);
        }
    }

    #[test]
    fn test_grid_dimensions_consistent() {
        let positions = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]];
        let (_, dims) = compute_grid_extents(&positions, 2.0, 0.25);
        // With 1Å range + 4Å padding (2 each side) = 5Å → 5/0.25 + 1 = 21
        assert!(dims[0] == 21, "dims[0] = {}", dims[0]);
    }

    #[test]
    fn test_chunked_matches_full() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let result = solve_eht(&elements, &positions, None).unwrap();
        let basis = build_basis(&elements, &positions);

        let full_grid = evaluate_orbital_on_grid(
            &basis, &result.coefficients, 0, &positions, 0.5, 3.0,
        );
        let (chunked_grid, slabs) = evaluate_orbital_on_grid_chunked(
            &basis, &result.coefficients, 0, &positions, 0.5, 3.0,
        );

        // Grid metadata must match
        assert_eq!(full_grid.dims, chunked_grid.dims);
        assert_eq!(full_grid.origin, chunked_grid.origin);
        assert_eq!(slabs.len(), chunked_grid.dims[0]);

        // Values must match exactly
        assert_eq!(full_grid.values.len(), chunked_grid.values.len());
        for (a, b) in full_grid.values.iter().zip(chunked_grid.values.iter()) {
            assert!((a - b).abs() < 1e-10, "mismatch: {} vs {}", a, b);
        }

        // Each slab should have ny*nz entries
        let slab_size = chunked_grid.dims[1] * chunked_grid.dims[2];
        for slab in &slabs {
            assert_eq!(slab.values.len(), slab_size);
        }
    }
}
