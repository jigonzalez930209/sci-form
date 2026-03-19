//! Fukui function descriptors and reactivity ranking bindings.

use crate::system::coords_to_positions;
use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
pub(crate) struct FukuiDescriptorsPy {
    #[pyo3(get)]
    f_plus: Vec<f64>,
    #[pyo3(get)]
    f_minus: Vec<f64>,
    #[pyo3(get)]
    f_radical: Vec<f64>,
    #[pyo3(get)]
    dual_descriptor: Vec<f64>,
    #[pyo3(get)]
    condensed_atom_indices: Vec<usize>,
    #[pyo3(get)]
    condensed_f_plus: Vec<f64>,
    #[pyo3(get)]
    condensed_f_minus: Vec<f64>,
    #[pyo3(get)]
    condensed_f_radical: Vec<f64>,
    #[pyo3(get)]
    condensed_dual_descriptor: Vec<f64>,
    #[pyo3(get)]
    homo_energy: f64,
    #[pyo3(get)]
    lumo_energy: f64,
    #[pyo3(get)]
    gap: f64,
    #[pyo3(get)]
    num_atoms: usize,
    #[pyo3(get)]
    validity_notes: Vec<String>,
}
#[pymethods]
impl FukuiDescriptorsPy {
    fn __repr__(&self) -> String {
        format!(
            "FukuiDescriptors(n_atoms={}, gap={:.3} eV)",
            self.num_atoms, self.gap
        )
    }
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct ReactivitySiteScorePy {
    #[pyo3(get)]
    atom_index: usize,
    #[pyo3(get)]
    score: f64,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct ReactivityRankingPy {
    #[pyo3(get)]
    nucleophilic_attack_sites: Vec<ReactivitySiteScorePy>,
    #[pyo3(get)]
    electrophilic_attack_sites: Vec<ReactivitySiteScorePy>,
    #[pyo3(get)]
    radical_attack_sites: Vec<ReactivitySiteScorePy>,
    #[pyo3(get)]
    notes: Vec<String>,
}
#[pymethods]
impl ReactivityRankingPy {
    fn __repr__(&self) -> String {
        format!(
            "ReactivityRanking(nuc={}, elec={}, rad={})",
            self.nucleophilic_attack_sites.len(),
            self.electrophilic_attack_sites.len(),
            self.radical_attack_sites.len()
        )
    }
}

#[pyfunction]
fn fukui_descriptors(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<FukuiDescriptorsPy> {
    let positions = coords_to_positions(&coords);
    sci_form_core::compute_fukui_descriptors(&elements, &positions)
        .map(|r| FukuiDescriptorsPy {
            condensed_atom_indices: r.condensed.iter().map(|row| row.atom_index).collect(),
            condensed_f_plus: r.condensed.iter().map(|row| row.f_plus).collect(),
            condensed_f_minus: r.condensed.iter().map(|row| row.f_minus).collect(),
            condensed_f_radical: r.condensed.iter().map(|row| row.f_radical).collect(),
            condensed_dual_descriptor: r.condensed.iter().map(|row| row.dual_descriptor).collect(),
            f_plus: r.f_plus,
            f_minus: r.f_minus,
            f_radical: r.f_radical,
            dual_descriptor: r.dual_descriptor,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap: r.gap,
            num_atoms: r.num_atoms,
            validity_notes: r.validity_notes,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

#[pyfunction]
fn reactivity_ranking(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<ReactivityRankingPy> {
    let positions = coords_to_positions(&coords);
    sci_form_core::compute_reactivity_ranking(&elements, &positions)
        .map(|r| {
            let to_scores = |scores: Vec<sci_form_core::reactivity::ReactivitySiteScore>| {
                scores
                    .into_iter()
                    .map(|s| ReactivitySiteScorePy {
                        atom_index: s.atom_index,
                        score: s.score,
                    })
                    .collect::<Vec<_>>()
            };
            ReactivityRankingPy {
                nucleophilic_attack_sites: to_scores(r.nucleophilic_attack_sites),
                electrophilic_attack_sites: to_scores(r.electrophilic_attack_sites),
                radical_attack_sites: to_scores(r.radical_attack_sites),
                notes: r.notes,
            }
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(fukui_descriptors, m)?)?;
    m.add_function(wrap_pyfunction!(reactivity_ranking, m)?)?;
    m.add_class::<FukuiDescriptorsPy>()?;
    m.add_class::<ReactivitySiteScorePy>()?;
    m.add_class::<ReactivityRankingPy>()?;
    Ok(())
}
