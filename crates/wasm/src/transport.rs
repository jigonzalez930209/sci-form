//! Transport/streaming WASM bindings: Arrow batch, worker tasks.

use wasm_bindgen::prelude::*;

/// Pack a batch of conformer results into sci-form's columnar (Arrow-like) format.
///
/// This is not the Apache Arrow IPC format: the payload is a JSON
/// serialisation of sci-form's own columnar layout.
#[wasm_bindgen]
pub fn pack_batch_arrow(results_json: &str) -> String {
    let results: Vec<sci_form::ConformerResult> = match serde_json::from_str(results_json) {
        Ok(r) => r,
        Err(e) => {
            return crate::helpers::json_error(&crate::helpers::bad_json(
                "'results'",
                results_json,
                &e,
            ))
        }
    };
    let batch = sci_form::transport::pack_conformers(&results);
    crate::helpers::serialize_or_error(&batch)
}

/// Pack a batch of conformer results into sci-form's columnar format.
///
/// Alias of [`pack_batch_arrow`], the name of the Rust root function
/// `sci_form::transport::pack_conformers` (blocker C-H1).
#[wasm_bindgen]
pub fn pack_conformers(results_json: &str) -> String {
    pack_batch_arrow(results_json)
}

/// Split a batch of SMILES into worker tasks for Web Worker dispatch.
///
/// `seed` is a 32-bit random seed; use [`split_worker_tasks_seed64`] for the
/// full 64-bit range.
#[wasm_bindgen]
pub fn split_worker_tasks(smiles_json: &str, n_workers: usize, seed: u32) -> String {
    split_worker_tasks_seed64(smiles_json, n_workers, seed as u64)
}

/// Split a batch of SMILES into worker tasks for Web Worker dispatch, 64-bit seed.
///
/// `seed` is a JavaScript `BigInt` (e.g. `42n`).
#[wasm_bindgen]
pub fn split_worker_tasks_seed64(smiles_json: &str, n_workers: usize, seed: u64) -> String {
    let smiles: Vec<String> = match serde_json::from_str(smiles_json) {
        Ok(s) => s,
        Err(e) => {
            return crate::helpers::json_error(&crate::helpers::bad_json(
                "'smiles'",
                smiles_json,
                &e,
            ))
        }
    };
    let tasks = sci_form::transport::split_batch(&smiles, n_workers, seed);
    crate::helpers::serialize_or_error(&tasks)
}

/// Estimate optimal number of Web Workers for a batch size.
#[wasm_bindgen]
pub fn estimate_workers(n_items: usize, max_workers: usize) -> usize {
    sci_form::transport::estimate_workers(n_items, max_workers)
}
