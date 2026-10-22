//! Conformer embedding WASM bindings.
//!
//! # Seed width (blocker C-H7)
//!
//! The Rust core and the Python bindings use a 64-bit seed.  `wasm-bindgen`
//! maps a `u64` parameter to a JavaScript `BigInt`, which is a breaking
//! change for every existing JS caller, so the 32-bit variants are kept and
//! the 64-bit variants are added next to them:
//!
//! | 32-bit (kept, backwards compatible) | 64-bit (full range) |
//! |-------------------------------------|---------------------|
//! | `embed(smiles, seed)`                | `embed_seed64(smiles, seed)` |
//! | `embed_coords(smiles, seed)`         | `embed_coords_seed64(smiles, seed)` |
//! | `embed_coords_typed(smiles, seed)`   | `embed_coords_typed_seed64(smiles, seed)` |
//! | `embed_batch(list, seed)`            | `embed_batch_seed64(list, seed)` |
//!
//! The 32-bit variants are **not** silent about their range: their doc
//! comments state that only the low 32 bits are used.

use wasm_bindgen::prelude::*;

use crate::helpers::{json_error, json_error_list, truncate};

/// Message used when an embedding fails.
///
/// Kept in one place so that the 32-bit and 64-bit variants report exactly
/// the same thing.
fn embed_error(smiles: &str, reason: &str) -> String {
    format!(
        "Embedding failed for SMILES '{}': {}",
        truncate(smiles, 120),
        reason
    )
}

/// Generate a 3D conformer from a SMILES string.
///
/// Returns a JSON string with the full result including atoms, bonds, and
/// coordinates.  When the embedding fails the function returns the canonical
/// WASM error payload `{"error": "..."}` — the same shape every other WASM
/// entry point uses — so a failure is never mistaken for a result object
/// (blocker C-M4).
///
/// `seed` is a 32-bit unsigned integer.  Use [`embed_seed64`] for the full
/// 64-bit seed range of the Rust core.
#[wasm_bindgen]
pub fn embed(smiles: &str, seed: u32) -> String {
    embed_seed64(smiles, seed as u64)
}

/// Generate a 3D conformer from a SMILES string, with a full 64-bit seed.
///
/// `seed` is a JavaScript `BigInt` (e.g. `42n`).  Numbers are also accepted
/// for seeds below 2^53.
#[wasm_bindgen]
pub fn embed_seed64(smiles: &str, seed: u64) -> String {
    let result = sci_form::embed(smiles, seed);
    if let Some(ref e) = result.error {
        return json_error(&embed_error(smiles, e));
    }
    crate::helpers::serialize_or_error(&result)
}

/// Generate 3D coordinates only (compact format).
///
/// Returns JSON: `{"coords": [x0,y0,z0,...], "num_atoms": N}` or `{"error": "..."}`.
///
/// `seed` is a 32-bit unsigned integer; use [`embed_coords_seed64`] for the
/// full 64-bit range.
#[wasm_bindgen]
pub fn embed_coords(smiles: &str, seed: u32) -> String {
    embed_coords_seed64(smiles, seed as u64)
}

/// Generate 3D coordinates only, with a full 64-bit seed.
///
/// `seed` is a JavaScript `BigInt` (e.g. `42n`).
#[wasm_bindgen]
pub fn embed_coords_seed64(smiles: &str, seed: u64) -> String {
    let result = sci_form::embed(smiles, seed);
    if let Some(ref e) = result.error {
        return json_error(&embed_error(smiles, e));
    }
    format!(
        "{{\"coords\":[{}],\"num_atoms\":{}}}",
        result
            .coords
            .iter()
            .map(|v| format!("{:.4}", v))
            .collect::<Vec<_>>()
            .join(","),
        result.num_atoms
    )
}

/// Generate 3D coordinates as a Float64Array (typed-array transfer).
///
/// Returns a Float64Array `[x0,y0,z0, x1,y1,z1,...]` or empty on failure.
///
/// `seed` is a 32-bit unsigned integer; use [`embed_coords_typed_seed64`] for
/// the full 64-bit range.
#[wasm_bindgen]
pub fn embed_coords_typed(smiles: &str, seed: u32) -> js_sys::Float64Array {
    embed_coords_typed_seed64(smiles, seed as u64)
}

/// Generate 3D coordinates as a Float64Array, with a full 64-bit seed.
///
/// `seed` is a JavaScript `BigInt` (e.g. `42n`).
#[wasm_bindgen]
pub fn embed_coords_typed_seed64(smiles: &str, seed: u64) -> js_sys::Float64Array {
    let result = sci_form::embed(smiles, seed);
    if result.error.is_some() || result.coords.is_empty() {
        return js_sys::Float64Array::new_with_length(0);
    }
    let arr = js_sys::Float64Array::new_with_length(result.coords.len() as u32);
    arr.copy_from(&result.coords);
    arr
}

/// Batch-embed multiple molecules from newline-separated SMILES.
///
/// Returns a JSON array of results.  A batch-level failure returns a JSON
/// array holding a single `{"error": "..."}` object so the array shape is
/// preserved.
///
/// `seed` is a 32-bit unsigned integer; use [`embed_batch_seed64`] for the
/// full 64-bit range.
#[wasm_bindgen]
pub fn embed_batch(smiles_list: &str, seed: u32) -> String {
    embed_batch_seed64(smiles_list, seed as u64)
}

/// Batch-embed multiple molecules, with a full 64-bit seed.
///
/// `seed` is a JavaScript `BigInt` (e.g. `42n`).
#[wasm_bindgen]
pub fn embed_batch_seed64(smiles_list: &str, seed: u64) -> String {
    let lines: Vec<&str> = smiles_list
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.trim())
        .collect();
    let config = sci_form::ConformerConfig {
        seed,
        num_threads: 0,
    };
    match sci_form::embed_batch(&lines, &config) {
        Ok(results) => serde_json::to_string(&results)
            .unwrap_or_else(|e| json_error_list(&e.to_string())),
        Err(e) => json_error_list(&e),
    }
}

/// Parse a SMILES string and return molecular info (no 3D).
///
/// Returns JSON:
/// `{"num_atoms": N, "num_bonds": M, "atoms": [{"element": 6, "hybridization": "Sp3", "formal_charge": 0}, ...]}`
///
/// The `atoms` array matches the Python `parse()` payload (blocker C-M10).
#[wasm_bindgen]
pub fn parse_smiles(smiles: &str) -> String {
    match sci_form::parse(smiles) {
        Ok(mol) => {
            let n = mol.graph.node_count();
            let atoms: Vec<serde_json::Value> = (0..n)
                .map(|i| {
                    let idx = sci_form::graph::NodeIndex::new(i);
                    let atom = &mol.graph[idx];
                    serde_json::json!({
                        "element": atom.element,
                        "hybridization": format!("{:?}", atom.hybridization),
                        "formal_charge": atom.formal_charge,
                    })
                })
                .collect();
            serde_json::json!({
                "num_atoms": n,
                "num_bonds": mol.graph.edge_count(),
                "atoms": atoms,
            })
            .to_string()
        }
        Err(e) => json_error(&crate::helpers::smiles_error(&e, smiles)),
    }
}

/// Generate a diverse conformer ensemble and cluster it by Butina RMSD.
///
/// Embeds `n_conformers` times with consecutive seeds and keeps only the
/// cluster centroids (RMSD threshold `rmsd_cutoff` in Å).
///
/// `base_seed` is a full 64-bit seed (a JavaScript `BigInt`).
#[wasm_bindgen]
pub fn embed_diverse(
    smiles: &str,
    n_conformers: usize,
    rmsd_cutoff: f64,
    base_seed: u64,
) -> String {
    let results = sci_form::embed_diverse(smiles, n_conformers, rmsd_cutoff, base_seed);
    serde_json::to_string(&results).unwrap_or_else(|e| json_error(&e.to_string()))
}

/// Return the canonical SMILES representation of a molecule.
#[wasm_bindgen]
pub fn to_canonical_smiles(smiles: &str) -> String {
    match sci_form::to_canonical_smiles(smiles) {
        Ok(s) => serde_json::to_string(&s).unwrap_or_else(|e| json_error(&e.to_string())),
        Err(e) => json_error(&crate::helpers::smiles_error(&e, smiles)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// C-C1: the failure payload of `embed` must be valid JSON.
    #[test]
    fn embed_failure_is_a_valid_json_error_payload() {
        let payload = embed("INVALID_XYZ", 42);
        let parsed: serde_json::Value =
            serde_json::from_str(&payload).expect("the error payload must be valid JSON");
        let message = parsed
            .get("error")
            .and_then(|v| v.as_str())
            .expect("the payload must carry an error");
        assert!(message.contains("INVALID_XYZ"));
    }

    /// C-C1: a successful embed carries the conformer, not an error.
    #[test]
    fn embed_success_carries_the_conformer() {
        let payload = embed("CCO", 42);
        let parsed: serde_json::Value = serde_json::from_str(&payload).expect("valid JSON");
        assert!(
            parsed.get("error").map(|v| v.is_null()).unwrap_or(true),
            "a successful embed must not carry an error"
        );
        assert_eq!(parsed.get("num_atoms").and_then(|v| v.as_u64()), Some(9));
    }

    /// C-M10: `parse_smiles` must expose the same `atoms` array as the
    /// Python `parse()`.
    #[test]
    fn parse_smiles_exposes_the_atoms_array() {
        let payload = parse_smiles("CCO");
        let parsed: serde_json::Value = serde_json::from_str(&payload).expect("valid JSON");
        assert_eq!(parsed.get("num_atoms").and_then(|v| v.as_u64()), Some(9));
        let atoms = parsed
            .get("atoms")
            .and_then(|v| v.as_array())
            .expect("the atoms array");
        assert_eq!(atoms.len(), 9);
        assert_eq!(atoms[0].get("element").and_then(|v| v.as_u64()), Some(6));
        assert!(atoms[0].get("hybridization").is_some());
        assert!(atoms[0].get("formal_charge").is_some());
    }

    /// C-C1: a malformed SMILES must produce a parseable error payload.
    #[test]
    fn parse_smiles_failure_is_a_valid_json_error_payload() {
        let payload = parse_smiles("INVALID_XYZ");
        let parsed: serde_json::Value =
            serde_json::from_str(&payload).expect("the error payload must be valid JSON");
        assert!(parsed.get("error").is_some());
    }
}
