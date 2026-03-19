//! WASM bindings for ring detection (SSSR), ECFP fingerprints, Tanimoto
//! similarity, and Butina clustering.

use crate::helpers::{json_error, serialize_or_error};
use wasm_bindgen::prelude::*;

/// Compute the Smallest Set of Smallest Rings (SSSR) for a SMILES string.
///
/// Returns JSON `SssrResult`.
#[wasm_bindgen]
pub fn compute_sssr(smiles: &str) -> String {
    match sci_form::compute_sssr(smiles) {
        Ok(r) => serialize_or_error(&r),
        Err(e) => json_error(&e),
    }
}

/// Compute ECFP fingerprint for a SMILES string.
///
/// radius=2 → ECFP4. Returns JSON `ECFPFingerprint`.
#[wasm_bindgen]
pub fn compute_ecfp(smiles: &str, radius: usize, n_bits: usize) -> String {
    match sci_form::compute_ecfp(smiles, radius, n_bits) {
        Ok(fp) => serialize_or_error(&fp),
        Err(e) => json_error(&e),
    }
}

/// Compute Tanimoto similarity between two ECFP fingerprint JSON objects.
///
/// Returns JSON `{"tanimoto": 0.85}`.
#[wasm_bindgen]
pub fn compute_tanimoto(fp1_json: &str, fp2_json: &str) -> String {
    let fp1: sci_form::rings::ecfp::ECFPFingerprint = match serde_json::from_str(fp1_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad fp1: {e}")),
    };
    let fp2: sci_form::rings::ecfp::ECFPFingerprint = match serde_json::from_str(fp2_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad fp2: {e}")),
    };
    let t = sci_form::compute_tanimoto(&fp1, &fp2);
    format!("{{\"tanimoto\":{t}}}")
}

/// Butina clustering on conformers (flat coordinate arrays in JSON).
///
/// `conformers_json`: JSON array of arrays, e.g. `[[x0,y0,z0,...], [x0,y0,z0,...]]`.
/// Returns JSON `ClusterResult`.
#[wasm_bindgen]
pub fn butina_cluster(conformers_json: &str, rmsd_cutoff: f64) -> String {
    let conformers: Vec<Vec<f64>> = match serde_json::from_str(conformers_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad conformers: {e}")),
    };
    let r = sci_form::butina_cluster(&conformers, rmsd_cutoff);
    serialize_or_error(&r)
}

/// Compute all-pairs RMSD matrix for conformers.
///
/// Returns JSON 2D array.
#[wasm_bindgen]
pub fn compute_rmsd_matrix(conformers_json: &str) -> String {
    let conformers: Vec<Vec<f64>> = match serde_json::from_str(conformers_json) {
        Ok(v) => v,
        Err(e) => return json_error(&format!("bad conformers: {e}")),
    };
    let m = sci_form::compute_rmsd_matrix(&conformers);
    serialize_or_error(&m)
}
