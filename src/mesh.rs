//! Generic orbital mesh generation for all electronic-structure methods.
//!
//! Provides orbital volumetric grids and isosurface meshes that can be used
//! with any implemented method (EHT, PM3, xTB, HF-3c).
//!
//! The orbital shape is evaluated using the EHT basis set (STO-nG Gaussians),
//! which provides a good visual representation of molecular orbitals regardless
//! of the method used to determine the MO coefficients and energies.

use serde::{Deserialize, Serialize};

use crate::eht;

/// Which electronic-structure method to use for orbital mesh generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MeshMethod {
    /// Extended Hückel Theory.
    Eht,
    /// PM3 semi-empirical.
    Pm3,
    /// GFN0-xTB tight-binding.
    Xtb,
    /// HF-3c composite method.
    Hf3c,
}

/// Result of orbital mesh generation including method context.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrbitalMeshResult {
    /// The generated isosurface mesh.
    pub mesh: eht::IsosurfaceMesh,
    /// The volumetric grid (orbital values on a 3D grid).
    pub grid: eht::VolumetricGrid,
    /// Method used.
    pub method: MeshMethod,
    /// MO index that was visualized.
    pub mo_index: usize,
    /// HOMO index for reference.
    pub homo_index: usize,
    /// Orbital energies (eV) from the chosen method.
    pub orbital_energies: Vec<f64>,
    /// HOMO-LUMO gap (eV).
    pub gap: f64,
}

/// Generate an orbital mesh using EHT basis with EHT MO coefficients.
///
/// This is the standard EHT orbital visualization path.
fn compute_eht_mesh(
    elements: &[u8],
    positions: &[[f64; 3]],
    mo_index: usize,
    spacing: f64,
    padding: f64,
    isovalue: f32,
) -> Result<OrbitalMeshResult, String> {
    let result = eht::solve_eht(elements, positions, None)?;
    let basis = eht::basis::build_basis(elements, positions);
    let grid = eht::evaluate_orbital_on_grid(
        &basis,
        &result.coefficients,
        mo_index,
        positions,
        spacing,
        padding,
    );
    let mesh = eht::marching_cubes(&grid, isovalue);

    Ok(OrbitalMeshResult {
        mesh,
        grid,
        method: MeshMethod::Eht,
        mo_index,
        homo_index: result.homo_index,
        orbital_energies: result.energies.clone(),
        gap: result.gap,
    })
}

/// Generate an orbital mesh using PM3 orbital energies for identification,
/// but EHT basis for the spatial orbital shape.
fn compute_pm3_mesh(
    elements: &[u8],
    positions: &[[f64; 3]],
    mo_index: usize,
    spacing: f64,
    padding: f64,
    isovalue: f32,
) -> Result<OrbitalMeshResult, String> {
    let pm3 = crate::pm3::solve_pm3(elements, positions)?;
    let eht_result = eht::solve_eht(elements, positions, None)?;
    let basis = eht::basis::build_basis(elements, positions);
    let grid = eht::evaluate_orbital_on_grid(
        &basis,
        &eht_result.coefficients,
        mo_index,
        positions,
        spacing,
        padding,
    );
    let mesh = eht::marching_cubes(&grid, isovalue);

    let homo_idx = if pm3.n_electrons > 0 {
        pm3.n_electrons / 2 - 1
    } else {
        0
    };

    Ok(OrbitalMeshResult {
        mesh,
        grid,
        method: MeshMethod::Pm3,
        mo_index,
        homo_index: homo_idx,
        orbital_energies: pm3.orbital_energies,
        gap: pm3.gap,
    })
}

/// Generate an orbital mesh using xTB orbital energies for identification,
/// but EHT basis for the spatial orbital shape.
fn compute_xtb_mesh(
    elements: &[u8],
    positions: &[[f64; 3]],
    mo_index: usize,
    spacing: f64,
    padding: f64,
    isovalue: f32,
) -> Result<OrbitalMeshResult, String> {
    let xtb = crate::xtb::solve_xtb(elements, positions)?;
    let eht_result = eht::solve_eht(elements, positions, None)?;
    let basis = eht::basis::build_basis(elements, positions);
    let grid = eht::evaluate_orbital_on_grid(
        &basis,
        &eht_result.coefficients,
        mo_index,
        positions,
        spacing,
        padding,
    );
    let mesh = eht::marching_cubes(&grid, isovalue);

    let homo_idx = if xtb.n_electrons > 0 {
        xtb.n_electrons / 2 - 1
    } else {
        0
    };

    Ok(OrbitalMeshResult {
        mesh,
        grid,
        method: MeshMethod::Xtb,
        mo_index,
        homo_index: homo_idx,
        orbital_energies: xtb.orbital_energies,
        gap: xtb.gap,
    })
}

/// Generate an orbital mesh using HF-3c orbital energies for identification,
/// but EHT basis for the spatial orbital shape.
fn compute_hf3c_mesh(
    elements: &[u8],
    positions: &[[f64; 3]],
    mo_index: usize,
    spacing: f64,
    padding: f64,
    isovalue: f32,
) -> Result<OrbitalMeshResult, String> {
    let config = crate::hf::HfConfig::default();
    let hf = crate::hf::api::solve_hf3c(elements, positions, &config)?;
    let eht_result = eht::solve_eht(elements, positions, None)?;
    let basis = eht::basis::build_basis(elements, positions);
    let grid = eht::evaluate_orbital_on_grid(
        &basis,
        &eht_result.coefficients,
        mo_index,
        positions,
        spacing,
        padding,
    );
    let mesh = eht::marching_cubes(&grid, isovalue);

    // Count electrons for HOMO index
    let n_electrons: usize = elements.iter().map(|&z| z as usize).sum();
    let homo_idx = if n_electrons > 0 {
        n_electrons / 2 - 1
    } else {
        0
    };

    let gap = if hf.orbital_energies.len() > homo_idx + 1 {
        hf.orbital_energies[homo_idx + 1] - hf.orbital_energies[homo_idx]
    } else {
        0.0
    };

    Ok(OrbitalMeshResult {
        mesh,
        grid,
        method: MeshMethod::Hf3c,
        mo_index,
        homo_index: homo_idx,
        orbital_energies: hf.orbital_energies,
        gap,
    })
}

/// Compute an orbital isosurface mesh for the specified method.
///
/// - `elements`: atomic numbers
/// - `positions`: Cartesian coordinates in Å
/// - `method`: which electronic-structure method to use
/// - `mo_index`: which molecular orbital to visualize
/// - `spacing`: grid spacing in Å (default 0.2)
/// - `padding`: padding around molecule in Å (default 3.0)
/// - `isovalue`: isosurface threshold (default 0.02)
pub fn compute_orbital_mesh(
    elements: &[u8],
    positions: &[[f64; 3]],
    method: MeshMethod,
    mo_index: usize,
    spacing: f64,
    padding: f64,
    isovalue: f32,
) -> Result<OrbitalMeshResult, String> {
    match method {
        MeshMethod::Eht => {
            compute_eht_mesh(elements, positions, mo_index, spacing, padding, isovalue)
        }
        MeshMethod::Pm3 => {
            compute_pm3_mesh(elements, positions, mo_index, spacing, padding, isovalue)
        }
        MeshMethod::Xtb => {
            compute_xtb_mesh(elements, positions, mo_index, spacing, padding, isovalue)
        }
        MeshMethod::Hf3c => {
            compute_hf3c_mesh(elements, positions, mo_index, spacing, padding, isovalue)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn water_positions() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![8, 1, 1],
            vec![[0.0, 0.0, 0.0], [0.0, 0.757, 0.587], [0.0, -0.757, 0.587]],
        )
    }

    #[test]
    fn test_eht_mesh_water() {
        let (elements, positions) = water_positions();
        let result =
            compute_orbital_mesh(&elements, &positions, MeshMethod::Eht, 0, 0.4, 3.0, 0.02);
        assert!(result.is_ok());
        let r = result.unwrap();
        assert!(r.grid.num_points() > 0);
        assert_eq!(r.method, MeshMethod::Eht);
    }

    #[test]
    fn test_pm3_mesh_water() {
        let (elements, positions) = water_positions();
        let result =
            compute_orbital_mesh(&elements, &positions, MeshMethod::Pm3, 0, 0.4, 3.0, 0.02);
        assert!(result.is_ok());
        let r = result.unwrap();
        assert!(r.grid.num_points() > 0);
        assert_eq!(r.method, MeshMethod::Pm3);
    }

    #[test]
    fn test_xtb_mesh_water() {
        let (elements, positions) = water_positions();
        let result =
            compute_orbital_mesh(&elements, &positions, MeshMethod::Xtb, 0, 0.4, 3.0, 0.02);
        assert!(result.is_ok());
        let r = result.unwrap();
        assert!(r.grid.num_points() > 0);
        assert_eq!(r.method, MeshMethod::Xtb);
    }

    #[test]
    fn test_hf3c_mesh_water() {
        let (elements, positions) = water_positions();
        let result =
            compute_orbital_mesh(&elements, &positions, MeshMethod::Hf3c, 0, 0.4, 3.0, 0.02);
        assert!(result.is_ok());
        let r = result.unwrap();
        assert!(r.grid.num_points() > 0);
        assert_eq!(r.method, MeshMethod::Hf3c);
    }
}
