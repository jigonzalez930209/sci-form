use wasm_bindgen::prelude::*;

/// Library version.
#[wasm_bindgen]
pub fn version() -> String {
    sci_form::version()
}

/// Generate a 3D conformer from a SMILES string.
///
/// Returns a JSON string with the full result including atoms, bonds, and coordinates.
#[wasm_bindgen]
pub fn embed(smiles: &str, seed: u32) -> String {
    let result = sci_form::embed(smiles, seed as u64);
    serde_json::to_string(&result).unwrap_or_else(|e| {
        format!("{{\"error\":\"{}\"}}", e)
    })
}

/// Generate 3D coordinates only (compact format).
///
/// Returns JSON: {"coords": [x0,y0,z0,...], "num_atoms": N} or {"error": "..."}
#[wasm_bindgen]
pub fn embed_coords(smiles: &str, seed: u32) -> String {
    let result = sci_form::embed(smiles, seed as u64);
    if let Some(ref e) = result.error {
        return format!("{{\"error\":\"{}\"}}", e.replace('"', "\\'"));
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

/// Batch-embed multiple molecules from newline-separated SMILES.
///
/// Returns a JSON array of results.
#[wasm_bindgen]
pub fn embed_batch(smiles_list: &str, seed: u32) -> String {
    let lines: Vec<&str> = smiles_list
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.trim())
        .collect();
    let config = sci_form::ConformerConfig {
        seed: seed as u64,
        num_threads: 0,
    };
    let results = sci_form::embed_batch(&lines, &config);
    serde_json::to_string(&results).unwrap_or_else(|e| {
        format!("[{{\"error\":\"{}\"}}]", e)
    })
}

/// Parse a SMILES string and return molecular info (no 3D).
#[wasm_bindgen]
pub fn parse_smiles(smiles: &str) -> String {
    match sci_form::parse(smiles) {
        Ok(mol) => {
            let n = mol.graph.node_count();
            let nb = mol.graph.edge_count();
            format!("{{\"num_atoms\":{},\"num_bonds\":{}}}", n, nb)
        }
        Err(e) => format!("{{\"error\":\"{}\"}}", e.replace('"', "\\'")),
    }
}
