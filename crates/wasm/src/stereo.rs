//! WASM bindings for stereochemistry analysis.

use crate::helpers::serialize_or_error;
use wasm_bindgen::prelude::*;

/// Analyze stereochemistry: detect R/S stereocenters and E/Z double bonds.
///
/// `coords_flat` is a JSON array of flat coordinates (may be empty `"[]"` for topology-only).
/// Returns JSON `StereoAnalysis`.
#[wasm_bindgen]
pub fn analyze_stereo(smiles: &str, coords_flat: &str) -> String {
    let coords: Vec<f64> = match serde_json::from_str(coords_flat) {
        Ok(c) => c,
        Err(e) => return crate::helpers::json_error(&format!("bad coords: {e}")),
    };
    match sci_form::analyze_stereo(smiles, &coords) {
        Ok(r) => serialize_or_error(&r),
        Err(e) => crate::helpers::json_error(&e),
    }
}
