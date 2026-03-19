//! ML property prediction WASM bindings.

use wasm_bindgen::prelude::*;

/// Predict molecular properties using ML proxy models.
#[wasm_bindgen]
pub fn compute_ml_properties(smiles: &str) -> String {
    let mol = match sci_form::parse(smiles) {
        Ok(m) => m,
        Err(e) => return crate::helpers::json_error(&e),
    };
    let n = mol.graph.node_count();
    let elements: Vec<u8> = (0..n)
        .map(|i| mol.graph[sci_form::graph::NodeIndex::new(i)].element)
        .collect();
    let bonds: Vec<(usize, usize, u8)> = mol
        .graph
        .edge_indices()
        .map(|e| {
            let (a, b) = mol.graph.edge_endpoints(e).unwrap();
            let order = match mol.graph[e].order {
                sci_form::graph::BondOrder::Single => 1u8,
                sci_form::graph::BondOrder::Double => 2,
                sci_form::graph::BondOrder::Triple => 3,
                sci_form::graph::BondOrder::Aromatic => 2,
                sci_form::graph::BondOrder::Unknown => 1,
            };
            (a.index(), b.index(), order)
        })
        .collect();
    let desc = sci_form::compute_ml_descriptors(&elements, &bonds, &[], &[]);
    let result = sci_form::predict_ml_properties(&desc);
    format!(
        "{{\"logp\":{:.4},\"molar_refractivity\":{:.4},\"log_solubility\":{:.4},\"lipinski_violations\":{},\"lipinski_passes\":{},\"druglikeness\":{:.4}}}",
        result.logp, result.molar_refractivity, result.log_solubility,
        result.lipinski.violations, result.lipinski.passes, result.druglikeness
    )
}

/// Compute molecular descriptors from SMILES.
#[wasm_bindgen]
pub fn compute_molecular_descriptors(smiles: &str) -> String {
    let mol = match sci_form::parse(smiles) {
        Ok(m) => m,
        Err(e) => return crate::helpers::json_error(&e),
    };
    let n = mol.graph.node_count();
    let elements: Vec<u8> = (0..n)
        .map(|i| mol.graph[sci_form::graph::NodeIndex::new(i)].element)
        .collect();
    let bonds: Vec<(usize, usize, u8)> = mol
        .graph
        .edge_indices()
        .map(|e| {
            let (a, b) = mol.graph.edge_endpoints(e).unwrap();
            let order = match mol.graph[e].order {
                sci_form::graph::BondOrder::Single => 1u8,
                sci_form::graph::BondOrder::Double => 2,
                sci_form::graph::BondOrder::Triple => 3,
                sci_form::graph::BondOrder::Aromatic => 2,
                sci_form::graph::BondOrder::Unknown => 1,
            };
            (a.index(), b.index(), order)
        })
        .collect();
    let desc = sci_form::compute_ml_descriptors(&elements, &bonds, &[], &[]);
    crate::helpers::serialize_or_error(&desc)
}
