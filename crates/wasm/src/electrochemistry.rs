//! Electrochemistry WASM bindings: pKa and redox potentials.

use crate::helpers::{json_error, serialize_or_error};
use wasm_bindgen::prelude::*;

fn method(name: &str) -> Result<sci_form::electrochemistry::ElectrochemistryMethod, String> {
    match name.to_ascii_lowercase().as_str() {
        "" | "uhf" | "hf" => Ok(sci_form::electrochemistry::ElectrochemistryMethod::Uhf),
        "xtb" | "gfn2" => Ok(sci_form::electrochemistry::ElectrochemistryMethod::Xtb),
        other => Err(format!(
            "unknown electrochemistry method '{other}'; expected 'uhf' or 'xtb'"
        )),
    }
}

fn config(method_name: &str) -> Result<sci_form::electrochemistry::ElectrochemistryConfig, String> {
    Ok(sci_form::electrochemistry::ElectrochemistryConfig {
        method: method(method_name)?,
        ..Default::default()
    })
}

/// Absolute pKa from an explicit acid/conjugate-base pair.
#[wasm_bindgen]
pub fn compute_pka(acid_smiles: &str, base_smiles: &str, method_name: &str) -> String {
    let config = match config(method_name) {
        Ok(value) => value,
        Err(error) => return json_error(&error),
    };
    match sci_form::compute_pka(acid_smiles, base_smiles, None, &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// Relative pKa against a reference acid with a known experimental value.
#[wasm_bindgen]
pub fn compute_pka_relative(
    acid_smiles: &str,
    base_smiles: &str,
    reference_acid: &str,
    reference_base: &str,
    reference_pka: f64,
    method_name: &str,
) -> String {
    let config = match config(method_name) {
        Ok(value) => value,
        Err(error) => return json_error(&error),
    };
    let reference = sci_form::electrochemistry::PkaReference {
        acid_smiles: reference_acid.to_string(),
        base_smiles: reference_base.to_string(),
        pka: reference_pka,
    };
    match sci_form::compute_pka(acid_smiles, base_smiles, Some(&reference), &config) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// Absolute redox potential from an explicit oxidised/reduced pair.
#[wasm_bindgen]
pub fn compute_redox_potential(
    oxidized_smiles: &str,
    reduced_smiles: &str,
    n_electrons: u32,
    method_name: &str,
) -> String {
    let config = match config(method_name) {
        Ok(value) => value,
        Err(error) => return json_error(&error),
    };
    match sci_form::compute_redox_potential(
        oxidized_smiles,
        reduced_smiles,
        n_electrons,
        None,
        &config,
    ) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}

/// Relative redox potential against a reference couple with a known value.
#[wasm_bindgen]
pub fn compute_redox_potential_relative(
    oxidized_smiles: &str,
    reduced_smiles: &str,
    n_electrons: u32,
    reference_oxidized: &str,
    reference_reduced: &str,
    reference_potential_v: f64,
    method_name: &str,
) -> String {
    let config = match config(method_name) {
        Ok(value) => value,
        Err(error) => return json_error(&error),
    };
    let reference = sci_form::electrochemistry::RedoxReference {
        oxidized_smiles: reference_oxidized.to_string(),
        reduced_smiles: reference_reduced.to_string(),
        n_electrons,
        potential_v: reference_potential_v,
    };
    match sci_form::compute_redox_potential(
        oxidized_smiles,
        reduced_smiles,
        n_electrons,
        Some(&reference),
        &config,
    ) {
        Ok(result) => serialize_or_error(&result),
        Err(error) => json_error(&error),
    }
}
