//! Universal orbital mesh generation bindings (EHT/PM3/xTB/HF-3c).

use crate::system::coords_to_positions;
use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
pub(crate) struct OrbitalMeshResultPy {
    #[pyo3(get)]
    vertices: Vec<f32>,
    #[pyo3(get)]
    normals: Vec<f32>,
    #[pyo3(get)]
    indices: Vec<u32>,
    #[pyo3(get)]
    num_triangles: usize,
    #[pyo3(get)]
    method: String,
    #[pyo3(get)]
    mo_index: usize,
    #[pyo3(get)]
    homo_index: usize,
    #[pyo3(get)]
    orbital_energies: Vec<f64>,
    #[pyo3(get)]
    gap: f64,
}

/// Generate an orbital mesh for any supported method.
///
/// `elements`: list of atomic numbers.
/// `coords`: flat xyz list in Å.
/// `method`: "eht", "pm3", "xtb", or "hf3c".
/// `mo_index`: molecular orbital index to visualize.
/// `spacing`: grid spacing in Å.
/// `padding`: padding around molecule in Å.
/// `isovalue`: isosurface threshold.
#[pyfunction]
#[pyo3(signature = (elements, coords, method="eht", mo_index=0, spacing=0.2, padding=3.0, isovalue=0.02))]
fn orbital_mesh(
    elements: Vec<u8>,
    coords: Vec<f64>,
    method: &str,
    mo_index: usize,
    spacing: f64,
    padding: f64,
    isovalue: f32,
) -> PyResult<OrbitalMeshResultPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    let positions = coords_to_positions(&coords);
    sci_form_core::compute_orbital_mesh(
        &elements, &positions, method, mo_index, spacing, padding, isovalue,
    )
    .map(|r| OrbitalMeshResultPy {
        vertices: r.mesh.vertices,
        normals: r.mesh.normals,
        indices: r.mesh.indices,
        num_triangles: r.mesh.num_triangles,
        method: format!("{:?}", r.method).to_lowercase(),
        mo_index: r.mo_index,
        homo_index: r.homo_index,
        orbital_energies: r.orbital_energies,
        gap: r.gap,
    })
    .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(orbital_mesh, m)?)?;
    m.add_class::<OrbitalMeshResultPy>()?;
    Ok(())
}
