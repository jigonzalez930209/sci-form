//! Force field and pKa WASM bindings: UFF, MMFF94, empirical pKa.

use crate::helpers::*;
use wasm_bindgen::prelude::*;

/// Compute UFF force field energy.
#[wasm_bindgen]
pub fn compute_uff_energy(smiles: &str, coords: &str) -> String {
    let flat = match parse_flat_coords(coords) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_uff_energy(smiles, &flat) {
        Ok(energy) => format!("{{\"energy\":{:.6},\"unit\":\"kcal/mol\"}}", energy),
        Err(e) => json_error(&e),
    }
}

/// Compute UFF energy with aromaticity-informed heuristic correction metadata.
#[wasm_bindgen]
pub fn compute_uff_energy_with_aromatic_heuristics(smiles: &str, coords: &str) -> String {
    let flat = match parse_flat_coords(coords) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_uff_energy_with_aromatic_heuristics(smiles, &flat) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Compute MMFF94 force field energy.
#[wasm_bindgen]
pub fn compute_mmff94_energy(smiles: &str, coords: &str) -> String {
    let flat = match parse_flat_coords(coords) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_mmff94_energy(smiles, &flat) {
        Ok(energy) => format!("{{\"energy\":{:.6},\"unit\":\"kcal/mol\"}}", energy),
        Err(e) => json_error(&e),
    }
}

/// Estimate acidic/basic pKa sites from graph and charge heuristics.
#[wasm_bindgen]
pub fn compute_empirical_pka(smiles: &str) -> String {
    match sci_form::compute_empirical_pka(smiles) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}
