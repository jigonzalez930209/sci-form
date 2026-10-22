//! Isosurface mesh generation for molecular orbital visualization.
//!
//! High-level API that combines orbital grid evaluation + marching cubes
//! into a single pipeline with dual-lobe support (positive/negative ψ)
//! and optional normal smoothing.
//!
//! # Contract with callers
//!
//! The CPU fallback always works: the grid evaluation may run on the GPU
//! (when the optional `experimental-gpu` feature is enabled and a device is
//! present) but the mesh extraction itself always runs on the CPU with the
//! marching-cubes implementation. The returned
//! [`IsosurfaceReport`] states
//! which backend produced the grid.
//!
//! # Units
//!
//! Grid spacing, padding and vertex coordinates are in **Bohr**; the
//! isovalue is the orbital wavefunction value (dimensionless, normalised so
//! that ∫|ψ|²dV = 1).

use super::backend_report::IsosurfaceReport;
use super::marching_cubes::{marching_cubes_cpu, smooth_mesh_normals, McOutput};
use super::orbital_grid::{evaluate_orbital_with_report, GridParams};
use crate::scf::basis::BasisSet;
use nalgebra::DMatrix;

/// Configuration for isosurface generation.
#[derive(Debug, Clone)]
pub struct IsosurfaceConfig {
    /// Grid spacing in Bohr (smaller = finer mesh).
    pub spacing: f64,
    /// Padding around the molecule in Bohr.
    pub padding: f64,
    /// Isovalue for the surface. Typical orbital values: 0.02–0.05.
    ///
    /// Dimensionless: the orbital wavefunction is normalised so that
    /// ∫|ψ|²dV = 1.
    pub isovalue: f64,
    /// Whether to generate both positive and negative lobes.
    pub both_lobes: bool,
    /// Whether to smooth vertex normals.
    pub smooth_normals: bool,
}

impl Default for IsosurfaceConfig {
    fn default() -> Self {
        Self {
            spacing: 0.3,
            padding: 5.0,
            isovalue: 0.02,
            both_lobes: true,
            smooth_normals: true,
        }
    }
}

/// Dual-lobe orbital isosurface.
///
/// Holds the meshes of the positive (ψ > 0) and negative (ψ < 0) lobes of one
/// molecular orbital, extracted with the marching-cubes algorithm.
#[derive(Debug, Clone)]
pub struct OrbitalIsosurface {
    /// Mesh of the positive lobe (the region where ψ ≥ isovalue).
    pub positive_lobe: McOutput,
    /// Mesh of the negative lobe (the region where ψ ≤ −isovalue); empty when
    /// `both_lobes` was disabled.
    pub negative_lobe: McOutput,
    /// Index of the orbital in the MO coefficient matrix (0-based).
    pub orbital_index: usize,
    /// Isovalue the surface was extracted at (dimensionless; the orbital
    /// wavefunction is normalised to 1 electron).
    pub isovalue: f64,
    /// Total number of triangles in both lobes combined.
    pub total_triangles: usize,
}

/// Generate an isosurface mesh for a molecular orbital.
///
/// Uses GPU for grid evaluation when available, CPU fallback otherwise.
/// Returns both the mesh and a report describing which backend was used.
///
/// `positions` are the atom coordinates in **Bohr**; `orbital_index` is the
/// 0-based column of `mo_coefficients` to visualise.
pub fn generate_orbital_isosurface(
    basis: &BasisSet,
    mo_coefficients: &DMatrix<f64>,
    orbital_index: usize,
    positions: &[[f64; 3]],
    config: &IsosurfaceConfig,
) -> (OrbitalIsosurface, IsosurfaceReport) {
    let params = GridParams::from_molecule(positions, config.spacing, config.padding);

    let (grid, grid_report) =
        evaluate_orbital_with_report(basis, mo_coefficients, orbital_index, &params);

    let mut positive = marching_cubes_cpu(&grid.values, &grid.params, config.isovalue);

    let mut negative = if config.both_lobes {
        let neg_values: Vec<f64> = grid.values.iter().map(|v| -v).collect();
        marching_cubes_cpu(&neg_values, &grid.params, config.isovalue)
    } else {
        McOutput {
            vertices: Vec::new(),
            normals: Vec::new(),
            indices: Vec::new(),
            n_triangles: 0,
        }
    };

    if config.smooth_normals {
        smooth_mesh_normals(&mut positive);
        if config.both_lobes {
            smooth_mesh_normals(&mut negative);
        }
    }

    let total = positive.n_triangles + negative.n_triangles;

    let iso = OrbitalIsosurface {
        positive_lobe: positive,
        negative_lobe: negative,
        orbital_index,
        isovalue: config.isovalue,
        total_triangles: total,
    };

    let report = IsosurfaceReport {
        grid_backend: grid_report.backend,
        grid_used_gpu: grid_report.used_gpu,
        n_triangles: total,
        isovalue: config.isovalue,
    };

    (iso, report)
}

/// Generate HOMO isosurface.
///
/// The HOMO index is `n_electrons / 2 - 1` for a closed-shell system (0 when
/// there are no electrons). `positions` are in **Bohr**.
pub fn generate_homo_isosurface(
    basis: &BasisSet,
    mo_coefficients: &DMatrix<f64>,
    n_electrons: usize,
    positions: &[[f64; 3]],
    config: &IsosurfaceConfig,
) -> (OrbitalIsosurface, IsosurfaceReport) {
    let n_occ = n_electrons / 2;
    let homo_index = if n_occ > 0 { n_occ - 1 } else { 0 };
    generate_orbital_isosurface(basis, mo_coefficients, homo_index, positions, config)
}

/// Generate LUMO isosurface.
///
/// The LUMO index is `n_electrons / 2` (0-based), i.e. the first virtual
/// orbital of a closed-shell system. `positions` are in **Bohr**.
pub fn generate_lumo_isosurface(
    basis: &BasisSet,
    mo_coefficients: &DMatrix<f64>,
    n_electrons: usize,
    positions: &[[f64; 3]],
    config: &IsosurfaceConfig,
) -> (OrbitalIsosurface, IsosurfaceReport) {
    let lumo_index = n_electrons / 2;
    generate_orbital_isosurface(basis, mo_coefficients, lumo_index, positions, config)
}

/// Estimate mesh complexity for a given grid resolution.
///
/// `spacing` and `padding` are in **Bohr** and `positions` in **Bohr**.
/// Returns `(estimated_triangles, estimated_vertices)`; the triangle count is
/// derived from the number of grid voxels with the empirical rule that ~5 % of
/// the voxels contribute a triangle, and the vertex count is three times the
/// triangle count (unwelded mesh).
pub fn estimate_mesh_complexity(
    positions: &[[f64; 3]],
    spacing: f64,
    padding: f64,
) -> (usize, usize) {
    let params = GridParams::from_molecule(positions, spacing, padding);
    let n_voxels = (params.dimensions[0].saturating_sub(1))
        * (params.dimensions[1].saturating_sub(1))
        * (params.dimensions[2].saturating_sub(1));
    // ~5% of voxels produce triangles for molecular orbitals
    let est_triangles = n_voxels / 20;
    (est_triangles, est_triangles * 3)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scf::basis::BasisSet;

    #[test]
    fn test_isosurface_config_default() {
        let config = IsosurfaceConfig::default();
        assert!(config.spacing > 0.0);
        assert!(config.isovalue > 0.0);
        assert!(config.both_lobes);
        assert!(config.smooth_normals);
    }

    #[test]
    fn test_estimate_mesh_complexity() {
        let positions = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let (tri, vert) = estimate_mesh_complexity(&positions, 0.5, 3.0);
        assert!(tri > 0);
        assert_eq!(vert, tri * 3);
    }

    #[test]
    fn test_generate_orbital_isosurface_h2() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]];
        let basis = BasisSet::sto3g(&elements, &positions);

        let n = basis.n_basis;
        let mut coeffs = DMatrix::zeros(n, n);
        let c = 1.0 / (2.0f64).sqrt();
        coeffs[(0, 0)] = c;
        if n > 1 {
            coeffs[(1, 0)] = c;
        }

        let config = IsosurfaceConfig {
            spacing: 0.4,
            padding: 3.0,
            isovalue: 0.005,
            both_lobes: true,
            smooth_normals: false,
        };

        let (iso, report) = generate_orbital_isosurface(&basis, &coeffs, 0, &positions, &config);
        assert!(!report.grid_backend.is_empty());
        assert!((report.isovalue - 0.005).abs() < 1e-10);
        assert_eq!(iso.orbital_index, 0);
        // Pipeline runs without panic; triangle count depends on grid/isovalue
    }

    #[test]
    fn test_generate_homo_lumo() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]];
        let basis = BasisSet::sto3g(&elements, &positions);

        let n = basis.n_basis;
        let mut coeffs = DMatrix::zeros(n, n);
        coeffs[(0, 0)] = 1.0 / (2.0f64).sqrt();
        if n > 1 {
            coeffs[(1, 0)] = 1.0 / (2.0f64).sqrt();
            coeffs[(0, 1)] = 1.0 / (2.0f64).sqrt();
            coeffs[(1, 1)] = -1.0 / (2.0f64).sqrt();
        }

        let config = IsosurfaceConfig {
            spacing: 0.8,
            padding: 2.0,
            isovalue: 0.01,
            both_lobes: false,
            smooth_normals: false,
        };

        let (homo, _) = generate_homo_isosurface(&basis, &coeffs, 2, &positions, &config);
        assert_eq!(homo.orbital_index, 0); // HOMO for 2 electrons

        let (lumo, _) = generate_lumo_isosurface(&basis, &coeffs, 2, &positions, &config);
        assert_eq!(lumo.orbital_index, 1); // LUMO
    }
}
