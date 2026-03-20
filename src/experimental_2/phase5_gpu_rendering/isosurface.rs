//! Isosurface mesh generation and utilities.
//!
//! Provides high-level API for generating isosurface meshes from
//! molecular orbital or electron density data:
//!
//! 1. Evaluate orbital/density on a 3D grid
//! 2. Run marching cubes at the specified isovalue
//! 3. Smooth normals (optional)
//! 4. Return mesh ready for Three.js / WebGL rendering
//!
//! For orbital visualization, typically two isosurfaces are generated
//! at +isovalue (positive lobe) and -isovalue (negative lobe).

use super::marching_cubes_gpu::{marching_cubes_cpu, McOutput};
use super::orbital_evaluator::{evaluate_orbital, GridParams};
use crate::experimental_2::phase2_quantum_engine::basis_set::BasisSet;
use crate::experimental_2::types::ScfResult;

/// Configuration for isosurface generation.
#[derive(Debug, Clone)]
pub struct IsosurfaceConfig {
    /// Grid spacing in Bohr. Smaller = finer mesh.
    pub spacing: f64,
    /// Padding around the molecule in Bohr.
    pub padding: f64,
    /// Isovalue for the surface. For orbitals, typical: 0.02-0.05.
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

/// Dual-lobe isosurface (positive and negative lobes of an orbital).
#[derive(Debug, Clone)]
pub struct OrbitalIsosurface {
    /// Positive lobe mesh.
    pub positive_lobe: McOutput,
    /// Negative lobe mesh (empty if both_lobes = false).
    pub negative_lobe: McOutput,
    /// Orbital index.
    pub orbital_index: usize,
    /// Isovalue used.
    pub isovalue: f64,
    /// Total number of triangles (both lobes).
    pub total_triangles: usize,
}

/// Generate an isosurface mesh for a molecular orbital.
pub fn generate_orbital_isosurface(
    basis: &BasisSet,
    scf: &ScfResult,
    orbital_index: usize,
    positions: &[[f64; 3]],
    config: &IsosurfaceConfig,
) -> OrbitalIsosurface {
    // Create grid
    let params = GridParams::from_molecule(positions, config.spacing, config.padding);

    // Evaluate orbital on grid
    let grid = evaluate_orbital(basis, &scf.mo_coefficients, orbital_index, &params);

    // Positive lobe
    let mut positive = marching_cubes_cpu(&grid.values, &grid.params, config.isovalue);

    // Negative lobe
    let mut negative = if config.both_lobes {
        // For negative lobe, negate values and use positive isovalue
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

    OrbitalIsosurface {
        positive_lobe: positive,
        negative_lobe: negative,
        orbital_index,
        isovalue: config.isovalue,
        total_triangles: total,
    }
}

/// Generate the HOMO isosurface.
pub fn generate_homo_isosurface(
    basis: &BasisSet,
    scf: &ScfResult,
    positions: &[[f64; 3]],
    config: &IsosurfaceConfig,
) -> OrbitalIsosurface {
    let n_occ = scf.n_electrons / 2;
    let homo_index = if n_occ > 0 { n_occ - 1 } else { 0 };
    generate_orbital_isosurface(basis, scf, homo_index, positions, config)
}

/// Generate the LUMO isosurface.
pub fn generate_lumo_isosurface(
    basis: &BasisSet,
    scf: &ScfResult,
    positions: &[[f64; 3]],
    config: &IsosurfaceConfig,
) -> OrbitalIsosurface {
    let n_occ = scf.n_electrons / 2;
    generate_orbital_isosurface(basis, scf, n_occ, positions, config)
}

/// Simple vertex normal smoothing by averaging face normals at shared positions.
///
/// For each unique vertex position, averages all face normals that share that
/// position. This creates a smoother-looking surface.
fn smooth_mesh_normals(mesh: &mut McOutput) {
    if mesh.vertices.is_empty() {
        return;
    }

    let n_verts = mesh.vertices.len() / 3;

    // Find vertices at approximately the same position and average their normals
    // Using a simple spatial hash for efficiency
    let tolerance = 1e-4f32;

    // For each vertex, accumulate normals from matching vertices
    let mut smoothed_normals = vec![[0.0f32; 3]; n_verts];

    for i in 0..n_verts {
        let pi = [
            mesh.vertices[i * 3],
            mesh.vertices[i * 3 + 1],
            mesh.vertices[i * 3 + 2],
        ];

        let mut nx = 0.0f32;
        let mut ny = 0.0f32;
        let mut nz = 0.0f32;
        let mut count = 0u32;

        for j in 0..n_verts {
            let pj = [
                mesh.vertices[j * 3],
                mesh.vertices[j * 3 + 1],
                mesh.vertices[j * 3 + 2],
            ];

            let dx = pi[0] - pj[0];
            let dy = pi[1] - pj[1];
            let dz = pi[2] - pj[2];
            let dist2 = dx * dx + dy * dy + dz * dz;

            if dist2 < tolerance * tolerance {
                nx += mesh.normals[j * 3];
                ny += mesh.normals[j * 3 + 1];
                nz += mesh.normals[j * 3 + 2];
                count += 1;
            }
        }

        if count > 0 {
            let len = (nx * nx + ny * ny + nz * nz).sqrt();
            if len > 1e-10 {
                smoothed_normals[i] = [nx / len, ny / len, nz / len];
            } else {
                smoothed_normals[i] = [
                    mesh.normals[i * 3],
                    mesh.normals[i * 3 + 1],
                    mesh.normals[i * 3 + 2],
                ];
            }
        }
    }

    // Write back
    for i in 0..n_verts {
        mesh.normals[i * 3] = smoothed_normals[i][0];
        mesh.normals[i * 3 + 1] = smoothed_normals[i][1];
        mesh.normals[i * 3 + 2] = smoothed_normals[i][2];
    }
}

/// Estimate mesh complexity for a given grid resolution.
pub fn estimate_mesh_complexity(
    positions: &[[f64; 3]],
    spacing: f64,
    padding: f64,
) -> (usize, usize) {
    let params = GridParams::from_molecule(positions, spacing, padding);
    let n_voxels = (params.dimensions[0] - 1)
        * (params.dimensions[1] - 1)
        * (params.dimensions[2] - 1);

    // Empirical: ~5% of voxels produce triangles on average for molecular orbitals
    let est_triangles = n_voxels / 20;
    let est_vertices = est_triangles * 3;

    (est_triangles, est_vertices)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_isosurface_config_default() {
        let config = IsosurfaceConfig::default();
        assert!(config.spacing > 0.0);
        assert!(config.isovalue > 0.0);
        assert!(config.both_lobes);
    }

    #[test]
    fn test_estimate_mesh_complexity() {
        let positions = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let (tri, vert) = estimate_mesh_complexity(&positions, 0.5, 3.0);
        assert!(tri > 0, "Should estimate non-zero triangles");
        assert_eq!(vert, tri * 3);
    }
}
