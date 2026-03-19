//! Reactivity descriptors and molecular analysis WASM bindings.
//! Frontier descriptors, Fukui, reactivity ranking, UV-Vis, graph features, topology.

use crate::helpers::*;
use wasm_bindgen::prelude::*;

/// Compute atom-resolved HOMO/LUMO frontier descriptors.
#[wasm_bindgen]
pub fn compute_frontier_descriptors(elements: &str, coords_flat: &str) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_frontier_descriptors(&elems, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Compute Fukui-function workflows and condensed per-atom descriptors.
#[wasm_bindgen]
pub fn compute_fukui_descriptors(elements: &str, coords_flat: &str) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_fukui_descriptors(&elems, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Build empirical local-reactivity rankings.
#[wasm_bindgen]
pub fn compute_reactivity_ranking(elements: &str, coords_flat: &str) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_reactivity_ranking(&elems, &positions) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Build an exploratory UV-Vis-like spectrum from low-cost EHT transitions.
#[wasm_bindgen]
pub fn compute_uv_vis_spectrum(
    elements: &str,
    coords_flat: &str,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_uv_vis_spectrum(&elems, &positions, sigma, e_min, e_max, n_points) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Analyze aromaticity and graph-level stereocenters from SMILES.
#[wasm_bindgen]
pub fn analyze_graph_features(smiles: &str) -> String {
    match sci_form::analyze_graph_features(smiles) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Detect metal coordination geometry and return topology.
#[wasm_bindgen]
pub fn compute_topology(elements: &str, coords_flat: &str) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    serialize_or_error(&sci_form::compute_topology(&elems, &positions))
}
