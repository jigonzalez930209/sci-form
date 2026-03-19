//! Method comparison bindings.

use crate::system::{coords_to_positions, plan_to_py, scientific_method_name};
use pyo3::prelude::*;

fn comparison_status_name(status: sci_form_core::MethodComparisonStatus) -> String {
    match status {
        sci_form_core::MethodComparisonStatus::Success => "success",
        sci_form_core::MethodComparisonStatus::Unavailable => "unavailable",
        sci_form_core::MethodComparisonStatus::Error => "error",
    }
    .to_string()
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct MethodComparisonEntryPy {
    #[pyo3(get)]
    method: String,
    #[pyo3(get)]
    status: String,
    #[pyo3(get)]
    available: bool,
    #[pyo3(get)]
    confidence: String,
    #[pyo3(get)]
    confidence_score: f64,
    #[pyo3(get)]
    warnings: Vec<String>,
    #[pyo3(get)]
    limitations: Vec<String>,
    #[pyo3(get)]
    error: Option<String>,
    #[pyo3(get)]
    support_level: Option<String>,
    #[pyo3(get)]
    homo_energy: Option<f64>,
    #[pyo3(get)]
    lumo_energy: Option<f64>,
    #[pyo3(get)]
    gap: Option<f64>,
    #[pyo3(get)]
    energy_kcal_mol: Option<f64>,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct MethodComparisonResultPy {
    #[pyo3(get)]
    plan: crate::system::SystemMethodPlanPy,
    #[pyo3(get)]
    comparisons: Vec<MethodComparisonEntryPy>,
}

#[pyfunction]
#[pyo3(signature = (smiles, elements, coords, allow_experimental_eht=false))]
fn compare_methods(
    smiles: &str,
    elements: Vec<u8>,
    coords: Vec<f64>,
    allow_experimental_eht: bool,
) -> PyResult<MethodComparisonResultPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "coords length must be 3 * len(elements)",
        ));
    }
    let positions = coords_to_positions(&coords);
    let result =
        sci_form_core::compare_methods(smiles, &elements, &positions, allow_experimental_eht);
    let comparisons =
        result
            .comparisons
            .into_iter()
            .map(|entry| {
                let (
                    mut support_level,
                    mut homo_energy,
                    mut lumo_energy,
                    mut gap,
                    mut energy_kcal_mol,
                ) = (None, None, None, None, None);
                if let Some(payload) = entry.payload {
                    match payload {
                        sci_form_core::MethodComparisonPayload::Eht {
                            homo_energy: h,
                            lumo_energy: l,
                            gap: g,
                            support,
                        } => {
                            support_level = Some(format!("{:?}", support.level).to_lowercase());
                            homo_energy = Some(h);
                            lumo_energy = Some(l);
                            gap = Some(g);
                        }
                        sci_form_core::MethodComparisonPayload::Uff { energy_kcal_mol: e } => {
                            energy_kcal_mol = Some(e);
                        }
                    }
                }
                MethodComparisonEntryPy {
                    method: scientific_method_name(entry.method),
                    status: comparison_status_name(entry.status),
                    available: entry.available,
                    confidence: format!("{:?}", entry.confidence).to_lowercase(),
                    confidence_score: entry.confidence_score,
                    warnings: entry.warnings,
                    limitations: entry.limitations,
                    error: entry.error,
                    support_level,
                    homo_energy,
                    lumo_energy,
                    gap,
                    energy_kcal_mol,
                }
            })
            .collect();
    Ok(MethodComparisonResultPy {
        plan: plan_to_py(result.plan),
        comparisons,
    })
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compare_methods, m)?)?;
    m.add_class::<MethodComparisonEntryPy>()?;
    m.add_class::<MethodComparisonResultPy>()?;
    Ok(())
}
