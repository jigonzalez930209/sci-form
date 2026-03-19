//! Transport/streaming WASM bindings: Arrow batch, worker tasks.

use wasm_bindgen::prelude::*;

/// Pack a batch of conformer results into Arrow-compatible columnar format.
#[wasm_bindgen]
pub fn pack_batch_arrow(results_json: &str) -> String {
    let results: Vec<sci_form::ConformerResult> = match serde_json::from_str(results_json) {
        Ok(r) => r,
        Err(e) => return crate::helpers::json_error(&format!("bad JSON: {}", e)),
    };
    let batch = sci_form::transport::pack_conformers(&results);
    crate::helpers::serialize_or_error(&batch)
}

/// Split a batch of SMILES into worker tasks for Web Worker dispatch.
#[wasm_bindgen]
pub fn split_worker_tasks(smiles_json: &str, n_workers: usize, seed: u32) -> String {
    let smiles: Vec<String> = match serde_json::from_str(smiles_json) {
        Ok(s) => s,
        Err(e) => return crate::helpers::json_error(&format!("bad JSON: {}", e)),
    };
    let tasks = sci_form::transport::split_batch(&smiles, n_workers, seed as u64);
    crate::helpers::serialize_or_error(&tasks)
}

/// Estimate optimal number of Web Workers for a batch size.
#[wasm_bindgen]
pub fn estimate_workers(n_items: usize, max_workers: usize) -> usize {
    sci_form::transport::estimate_workers(n_items, max_workers)
}
