//! Stereochemistry bindings: CIP priorities, R/S, E/Z assignment.

use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
pub(crate) struct StereocenterPy {
    #[pyo3(get)]
    atom_index: usize,
    #[pyo3(get)]
    element: u8,
    #[pyo3(get)]
    substituent_indices: Vec<usize>,
    #[pyo3(get)]
    priorities: Vec<usize>,
    #[pyo3(get)]
    configuration: Option<String>,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct DoubleBondStereoPy {
    #[pyo3(get)]
    atom1: usize,
    #[pyo3(get)]
    atom2: usize,
    #[pyo3(get)]
    configuration: Option<String>,
    #[pyo3(get)]
    high_priority_sub1: Option<usize>,
    #[pyo3(get)]
    high_priority_sub2: Option<usize>,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct StereoAnalysisPy {
    #[pyo3(get)]
    stereocenters: Vec<StereocenterPy>,
    #[pyo3(get)]
    double_bonds: Vec<DoubleBondStereoPy>,
    #[pyo3(get)]
    n_stereocenters: usize,
    #[pyo3(get)]
    n_double_bonds: usize,
}

/// Analyze stereochemistry: detect R/S stereocenters and E/Z double bonds.
#[pyfunction]
#[pyo3(signature = (smiles, coords=vec![]))]
fn stereo_analysis(smiles: &str, coords: Vec<f64>) -> PyResult<StereoAnalysisPy> {
    sci_form_core::analyze_stereo(smiles, &coords)
        .map(|r| StereoAnalysisPy {
            stereocenters: r
                .stereocenters
                .iter()
                .map(|s| StereocenterPy {
                    atom_index: s.atom_index,
                    element: s.element,
                    substituent_indices: s.substituent_indices.clone(),
                    priorities: s.priorities.clone(),
                    configuration: s.configuration.clone(),
                })
                .collect(),
            double_bonds: r
                .double_bonds
                .iter()
                .map(|d| DoubleBondStereoPy {
                    atom1: d.atom1,
                    atom2: d.atom2,
                    configuration: d.configuration.clone(),
                    high_priority_sub1: d.high_priority_sub1,
                    high_priority_sub2: d.high_priority_sub2,
                })
                .collect(),
            n_stereocenters: r.n_stereocenters,
            n_double_bonds: r.n_double_bonds,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(stereo_analysis, m)?)?;
    m.add_class::<StereocenterPy>()?;
    m.add_class::<DoubleBondStereoPy>()?;
    m.add_class::<StereoAnalysisPy>()?;
    Ok(())
}
