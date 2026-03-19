//! UV-Vis, graph feature analysis, and topology analysis bindings.

use crate::system::coords_to_positions;
use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
pub(crate) struct UvVisPeakPy {
    #[pyo3(get)]
    energy_ev: f64,
    #[pyo3(get)]
    wavelength_nm: f64,
    #[pyo3(get)]
    intensity: f64,
    #[pyo3(get)]
    from_mo: usize,
    #[pyo3(get)]
    to_mo: usize,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct UvVisSpectrumPy {
    #[pyo3(get)]
    energies_ev: Vec<f64>,
    #[pyo3(get)]
    intensities: Vec<f64>,
    #[pyo3(get)]
    peaks: Vec<UvVisPeakPy>,
    #[pyo3(get)]
    sigma: f64,
    #[pyo3(get)]
    notes: Vec<String>,
}
#[pymethods]
impl UvVisSpectrumPy {
    fn __repr__(&self) -> String {
        format!(
            "UvVisSpectrum(n_points={}, n_peaks={})",
            self.energies_ev.len(),
            self.peaks.len()
        )
    }
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct GraphFeatureAnalysisPy {
    #[pyo3(get)]
    aromatic_atoms: Vec<bool>,
    #[pyo3(get)]
    aromatic_bonds: Vec<(usize, usize)>,
    #[pyo3(get)]
    tagged_tetrahedral_centers: Vec<usize>,
    #[pyo3(get)]
    inferred_tetrahedral_centers: Vec<usize>,
}
#[pymethods]
impl GraphFeatureAnalysisPy {
    fn __repr__(&self) -> String {
        format!(
            "GraphFeatureAnalysis(aromatic_atoms={}, tagged_centers={}, inferred_centers={})",
            self.aromatic_atoms.iter().filter(|v| **v).count(),
            self.tagged_tetrahedral_centers.len(),
            self.inferred_tetrahedral_centers.len()
        )
    }
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct MetalCoordinationCenterPy {
    #[pyo3(get)]
    atom_index: usize,
    #[pyo3(get)]
    element: u8,
    #[pyo3(get)]
    ligand_indices: Vec<usize>,
    #[pyo3(get)]
    coordination_number: usize,
    #[pyo3(get)]
    geometry: String,
    #[pyo3(get)]
    geometry_score: f64,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct TopologyAnalysisResultPy {
    #[pyo3(get)]
    metal_centers: Vec<MetalCoordinationCenterPy>,
    #[pyo3(get)]
    warnings: Vec<String>,
}

#[pyfunction]
#[pyo3(signature = (elements, coords, sigma=0.25, e_min=0.5, e_max=8.0, n_points=600))]
fn uvvis_spectrum(
    elements: Vec<u8>,
    coords: Vec<f64>,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> PyResult<UvVisSpectrumPy> {
    let positions = coords_to_positions(&coords);
    sci_form_core::compute_uv_vis_spectrum(&elements, &positions, sigma, e_min, e_max, n_points)
        .map(|r| UvVisSpectrumPy {
            energies_ev: r.energies_ev,
            intensities: r.intensities,
            peaks: r
                .peaks
                .into_iter()
                .map(|p| UvVisPeakPy {
                    energy_ev: p.energy_ev,
                    wavelength_nm: p.wavelength_nm,
                    intensity: p.intensity,
                    from_mo: p.from_mo,
                    to_mo: p.to_mo,
                })
                .collect(),
            sigma: r.sigma,
            notes: r.notes,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
fn graph_features(smiles: &str) -> PyResult<GraphFeatureAnalysisPy> {
    sci_form_core::analyze_graph_features(smiles)
        .map(|r| GraphFeatureAnalysisPy {
            aromatic_atoms: r.aromaticity.aromatic_atoms,
            aromatic_bonds: r.aromaticity.aromatic_bonds,
            tagged_tetrahedral_centers: r.stereocenters.tagged_tetrahedral_centers,
            inferred_tetrahedral_centers: r.stereocenters.inferred_tetrahedral_centers,
        })
        .map_err(pyo3::exceptions::PyValueError::new_err)
}

#[pyfunction]
fn topology_analysis(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<TopologyAnalysisResultPy> {
    let positions = coords_to_positions(&coords);
    let result = sci_form_core::compute_topology(&elements, &positions);
    Ok(TopologyAnalysisResultPy {
        metal_centers: result
            .metal_centers
            .into_iter()
            .map(|c| MetalCoordinationCenterPy {
                atom_index: c.atom_index,
                element: c.element,
                ligand_indices: c.ligand_indices,
                coordination_number: c.coordination_number,
                geometry: format!("{:?}", c.geometry).to_lowercase(),
                geometry_score: c.geometry_score,
            })
            .collect(),
        warnings: result.warnings,
    })
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(uvvis_spectrum, m)?)?;
    m.add_function(wrap_pyfunction!(graph_features, m)?)?;
    m.add_function(wrap_pyfunction!(topology_analysis, m)?)?;
    m.add_class::<UvVisPeakPy>()?;
    m.add_class::<UvVisSpectrumPy>()?;
    m.add_class::<GraphFeatureAnalysisPy>()?;
    m.add_class::<MetalCoordinationCenterPy>()?;
    m.add_class::<TopologyAnalysisResultPy>()?;
    Ok(())
}
