//! Electrochemistry bindings: relative pKa and redox potentials.

use pyo3::prelude::*;

use crate::errors::{engine_error, value_error};

#[pyclass(name = "PkaResult")]
#[derive(Clone)]
pub(crate) struct PkaResultPy {
    #[pyo3(get)]
    pka: f64,
    #[pyo3(get)]
    delta_g_kcal_mol: f64,
    #[pyo3(get)]
    relative: bool,
    #[pyo3(get)]
    acid_smiles: String,
    #[pyo3(get)]
    base_smiles: String,
    #[pyo3(get)]
    temperature_k: f64,
    #[pyo3(get)]
    proton_solvation_kcal_mol: f64,
    #[pyo3(get)]
    acid_energy_kcal_mol: f64,
    #[pyo3(get)]
    base_energy_kcal_mol: f64,
    #[pyo3(get)]
    acid_spread_kcal_mol: f64,
    #[pyo3(get)]
    base_spread_kcal_mol: f64,
    #[pyo3(get)]
    notes: Vec<String>,
}

#[pymethods]
impl PkaResultPy {
    fn __repr__(&self) -> String {
        format!(
            "PkaResult(pKa={:.3}, relative={}, dG={:.3} kcal/mol)",
            self.pka, self.relative, self.delta_g_kcal_mol
        )
    }
}

#[pyclass(name = "RedoxResult")]
#[derive(Clone)]
pub(crate) struct RedoxResultPy {
    #[pyo3(get)]
    potential_vs_reference_v: f64,
    #[pyo3(get)]
    absolute_potential_v: f64,
    #[pyo3(get)]
    delta_g_kcal_mol: f64,
    #[pyo3(get)]
    n_electrons: u32,
    #[pyo3(get)]
    relative: bool,
    #[pyo3(get)]
    oxidized_smiles: String,
    #[pyo3(get)]
    reduced_smiles: String,
    #[pyo3(get)]
    reference_absolute_potential_v: f64,
    #[pyo3(get)]
    oxidized_energy_kcal_mol: f64,
    #[pyo3(get)]
    reduced_energy_kcal_mol: f64,
    #[pyo3(get)]
    notes: Vec<String>,
}

#[pymethods]
impl RedoxResultPy {
    fn __repr__(&self) -> String {
        format!(
            "RedoxResult(E={:.3} V, relative={}, n={})",
            self.potential_vs_reference_v, self.relative, self.n_electrons
        )
    }
}

fn method(value: &str) -> PyResult<sci_form_core::electrochemistry::ElectrochemistryMethod> {
    match value.to_ascii_lowercase().as_str() {
        "uhf" | "hf" => Ok(sci_form_core::electrochemistry::ElectrochemistryMethod::Uhf),
        "xtb" | "gfn2" => Ok(sci_form_core::electrochemistry::ElectrochemistryMethod::Xtb),
        other => Err(value_error(format!(
            "unknown electrochemistry method '{other}'; expected 'uhf' or 'xtb'"
        ))),
    }
}

/// Relative or absolute pKa from an explicit acid/conjugate-base pair.
#[pyfunction]
#[pyo3(signature = (
    acid_smiles,
    base_smiles,
    reference_acid=None,
    reference_base=None,
    reference_pka=None,
    method_name="uhf",
    n_conformers=1,
))]
fn compute_pka(
    acid_smiles: &str,
    base_smiles: &str,
    reference_acid: Option<&str>,
    reference_base: Option<&str>,
    reference_pka: Option<f64>,
    method_name: &str,
    n_conformers: usize,
) -> PyResult<PkaResultPy> {
    let reference = match (reference_acid, reference_base, reference_pka) {
        (Some(acid), Some(base), Some(pka)) => {
            Some(sci_form_core::electrochemistry::PkaReference {
                acid_smiles: acid.to_string(),
                base_smiles: base.to_string(),
                pka,
            })
        }
        (None, None, None) => None,
        _ => {
            return Err(value_error(
                "reference_acid, reference_base and reference_pka must be supplied together",
            ))
        }
    };
    let config = sci_form_core::electrochemistry::ElectrochemistryConfig {
        method: method(method_name)?,
        n_conformers,
        ..Default::default()
    };
    let result = sci_form_core::compute_pka(acid_smiles, base_smiles, reference.as_ref(), &config)
        .map_err(engine_error)?;
    Ok(PkaResultPy {
        pka: result.pka,
        delta_g_kcal_mol: result.delta_g_kcal_mol,
        relative: result.relative,
        acid_smiles: result.acid_smiles,
        base_smiles: result.base_smiles,
        temperature_k: result.temperature_k,
        proton_solvation_kcal_mol: result.proton_solvation_kcal_mol,
        acid_energy_kcal_mol: result.acid.selected_kcal_mol,
        base_energy_kcal_mol: result.base.selected_kcal_mol,
        acid_spread_kcal_mol: result.acid.spread_kcal_mol,
        base_spread_kcal_mol: result.base.spread_kcal_mol,
        notes: result.notes,
    })
}

/// Relative or absolute redox potential from an explicit oxidised/reduced pair.
#[pyfunction]
#[pyo3(signature = (
    oxidized_smiles,
    reduced_smiles,
    n_electrons=1,
    reference_oxidized=None,
    reference_reduced=None,
    reference_potential=None,
    method_name="uhf",
    n_conformers=1,
))]
fn compute_redox_potential(
    oxidized_smiles: &str,
    reduced_smiles: &str,
    n_electrons: u32,
    reference_oxidized: Option<&str>,
    reference_reduced: Option<&str>,
    reference_potential: Option<f64>,
    method_name: &str,
    n_conformers: usize,
) -> PyResult<RedoxResultPy> {
    let reference = match (reference_oxidized, reference_reduced, reference_potential) {
        (Some(oxidized), Some(reduced), Some(potential)) => {
            Some(sci_form_core::electrochemistry::RedoxReference {
                oxidized_smiles: oxidized.to_string(),
                reduced_smiles: reduced.to_string(),
                n_electrons,
                potential_v: potential,
            })
        }
        (None, None, None) => None,
        _ => {
            return Err(value_error(
                "reference_oxidized, reference_reduced and reference_potential must be supplied together",
            ))
        }
    };
    let config = sci_form_core::electrochemistry::ElectrochemistryConfig {
        method: method(method_name)?,
        n_conformers,
        ..Default::default()
    };
    let result = sci_form_core::compute_redox_potential(
        oxidized_smiles,
        reduced_smiles,
        n_electrons,
        reference.as_ref(),
        &config,
    )
    .map_err(engine_error)?;
    Ok(RedoxResultPy {
        potential_vs_reference_v: result.potential_vs_reference_v,
        absolute_potential_v: result.absolute_potential_v,
        delta_g_kcal_mol: result.delta_g_kcal_mol,
        n_electrons: result.n_electrons,
        relative: result.relative,
        oxidized_smiles: result.oxidized_smiles,
        reduced_smiles: result.reduced_smiles,
        reference_absolute_potential_v: result.reference_absolute_potential_v,
        oxidized_energy_kcal_mol: result.oxidized.selected_kcal_mol,
        reduced_energy_kcal_mol: result.reduced.selected_kcal_mol,
        notes: result.notes,
    })
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_pka, m)?)?;
    m.add_function(wrap_pyfunction!(compute_redox_potential, m)?)?;
    m.add_class::<PkaResultPy>()?;
    m.add_class::<RedoxResultPy>()?;
    Ok(())
}
