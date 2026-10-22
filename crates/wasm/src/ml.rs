//! ML property prediction WASM bindings.

use wasm_bindgen::prelude::*;

/// Predict molecular properties using ML proxy models.
///
/// Canonical name matching the Rust root (`sci_form::predict_ml_properties`).
/// Returns JSON:
/// `{"logp":…, "molar_refractivity":…, "log_solubility":…, "lipinski_violations":…,
/// "lipinski_passes":…, "druglikeness":…}`
#[wasm_bindgen]
pub fn predict_ml_properties(smiles: &str) -> String {
    let mol = match sci_form::parse(smiles) {
        Ok(m) => m,
        Err(e) => return crate::helpers::json_error(&crate::helpers::smiles_error(&e, smiles)),
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
                sci_form::graph::BondOrder::Aromatic => 4,
                sci_form::graph::BondOrder::Unknown => 1,
                // `BondOrder` is `#[non_exhaustive]`: a bond order added after
                // this binding was written falls back to a single bond, the
                // same conservative choice as `Unknown`.
                _ => 1,
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

/// Predict molecular properties using ML proxy models.
///
/// **Deprecated alias** of [`predict_ml_properties`]; kept so existing JS
/// callers keep working.  Use `predict_ml_properties`, the name of the Rust
/// root function.
#[wasm_bindgen]
pub fn compute_ml_properties(smiles: &str) -> String {
    predict_ml_properties(smiles)
}

/// Compute molecular descriptors from SMILES.
///
/// Canonical name matching the Rust root (`sci_form::compute_ml_descriptors`).
#[wasm_bindgen]
pub fn compute_ml_descriptors(smiles: &str) -> String {
    let mol = match sci_form::parse(smiles) {
        Ok(m) => m,
        Err(e) => return crate::helpers::json_error(&crate::helpers::smiles_error(&e, smiles)),
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
                sci_form::graph::BondOrder::Aromatic => 4,
                sci_form::graph::BondOrder::Unknown => 1,
                // `BondOrder` is `#[non_exhaustive]`: a bond order added after
                // this binding was written falls back to a single bond, the
                // same conservative choice as `Unknown`.
                _ => 1,
            };
            (a.index(), b.index(), order)
        })
        .collect();
    let desc = sci_form::compute_ml_descriptors(&elements, &bonds, &[], &[]);
    crate::helpers::serialize_or_error(&desc)
}

/// Compute molecular descriptors from SMILES.
///
/// **Deprecated alias** of [`compute_ml_descriptors`]; kept so existing JS
/// callers keep working.  Use `compute_ml_descriptors`, the name of the Rust
/// root function.
#[wasm_bindgen]
pub fn compute_molecular_descriptors(smiles: &str) -> String {
    compute_ml_descriptors(smiles)
}

/// Compute WHIM 3D molecular descriptors.
///
/// `weighting`: `unit`, `mass`, `volume`, `electronegativity` or `polarizability`.
#[wasm_bindgen]
pub fn compute_whim(elements_json: &str, coords_flat_json: &str, weighting: &str) -> String {
    let (elements, positions) =
        match crate::helpers::parse_elements_and_positions(elements_json, coords_flat_json) {
            Ok(v) => v,
            Err(e) => return crate::helpers::json_error(&e),
        };
    let w = match weighting.to_ascii_lowercase().as_str() {
        "unit" | "unweighted" => sci_form::ml::whim::WhimWeighting::Unit,
        "mass" => sci_form::ml::whim::WhimWeighting::Mass,
        "volume" | "vdw" => sci_form::ml::whim::WhimWeighting::Volume,
        "electronegativity" | "en" => sci_form::ml::whim::WhimWeighting::Electronegativity,
        "polarizability" => sci_form::ml::whim::WhimWeighting::Polarizability,
        other => {
            return crate::helpers::json_error(&format!(
                "unknown WHIM weighting '{other}' (expected unit, mass, volume, \
                 electronegativity or polarizability)"
            ))
        }
    };
    crate::helpers::serialize_or_error(&sci_form::ml::whim::compute_whim(&elements, &positions, w))
}

/// Compute radial distribution function (RDF) descriptors.
///
/// `charges_json`: JSON array of partial charges; pass `""` to derive the
/// charges with the Gasteiger solver from geometry-perceived bonds.
#[wasm_bindgen]
pub fn compute_rdf(
    elements_json: &str,
    coords_flat_json: &str,
    charges_json: &str,
    beta: f64,
    r_max: f64,
    n_points: usize,
) -> String {
    let (elements, positions) =
        match crate::helpers::parse_elements_and_positions(elements_json, coords_flat_json) {
            Ok(v) => v,
            Err(e) => return crate::helpers::json_error(&e),
        };
    let charges = if charges_json.is_empty() {
        match crate::helpers::gasteiger_charges_from_geometry(&elements, &positions) {
            Ok(q) => q,
            Err(e) => return crate::helpers::json_error(&e),
        }
    } else {
        match serde_json::from_str::<Vec<f64>>(charges_json) {
            Ok(v) => v,
            Err(e) => {
                return crate::helpers::json_error(&crate::helpers::bad_json(
                    "'charges'",
                    charges_json,
                    &e,
                ))
            }
        }
    };
    crate::helpers::serialize_or_error(&sci_form::ml::rdf_descriptors::compute_rdf(
        &elements, &positions, &charges, beta, r_max, n_points,
    ))
}

/// Compute GETAWAY 3D molecular descriptors.
///
/// `bonds_json`: JSON array of `[i, j]` atom-index pairs; pass `""` to
/// perceive the bonds from the geometry.
#[wasm_bindgen]
pub fn compute_getaway(
    elements_json: &str,
    coords_flat_json: &str,
    bonds_json: &str,
    max_lag: usize,
) -> String {
    let (elements, positions) =
        match crate::helpers::parse_elements_and_positions(elements_json, coords_flat_json) {
            Ok(v) => v,
            Err(e) => return crate::helpers::json_error(&e),
        };
    let bonds = if bonds_json.is_empty() {
        crate::helpers::perceive_bonds(&elements, &positions, 0.3)
    } else {
        match serde_json::from_str::<Vec<(usize, usize)>>(bonds_json) {
            Ok(v) => v,
            Err(e) => {
                return crate::helpers::json_error(&crate::helpers::bad_json(
                    "'bonds'", bonds_json, &e,
                ))
            }
        }
    };
    crate::helpers::serialize_or_error(&sci_form::ml::getaway::compute_getaway(
        &elements, &positions, &bonds, max_lag,
    ))
}

fn rmse(predictions: &[f64], targets: &[f64]) -> f64 {
    if predictions.is_empty() {
        return 0.0;
    }
    let sum: f64 = predictions
        .iter()
        .zip(targets.iter())
        .map(|(p, t)| (p - t) * (p - t))
        .sum();
    (sum / predictions.len() as f64).sqrt()
}

/// Train a Random Forest regressor and return its training-set predictions.
///
/// `features_json`: JSON array of feature vectors; `targets_json`: JSON array
/// of target values.
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn train_random_forest(
    features_json: &str,
    targets_json: &str,
    n_trees: usize,
    max_depth: usize,
    min_samples_split: usize,
    min_samples_leaf: usize,
    sample_fraction: f64,
    seed: u64,
) -> String {
    let features: Vec<Vec<f64>> = match serde_json::from_str(features_json) {
        Ok(v) => v,
        Err(e) => {
            return crate::helpers::json_error(&crate::helpers::bad_json(
                "'features'",
                features_json,
                &e,
            ))
        }
    };
    let targets: Vec<f64> = match serde_json::from_str(targets_json) {
        Ok(v) => v,
        Err(e) => {
            return crate::helpers::json_error(&crate::helpers::bad_json(
                "'targets'",
                targets_json,
                &e,
            ))
        }
    };
    if features.len() != targets.len() {
        return crate::helpers::json_error(&format!(
            "features length {} != targets length {}",
            features.len(),
            targets.len()
        ));
    }
    let config = sci_form::ml::advanced_models::RandomForestConfig {
        n_trees,
        tree_config: sci_form::ml::advanced_models::TreeConfig {
            max_depth,
            min_samples_split,
            min_samples_leaf,
            max_features: 0,
        },
        sample_fraction,
        seed,
    };
    let model = sci_form::ml::advanced_models::train_random_forest(&features, &targets, &config);
    let predictions: Vec<f64> = features.iter().map(|f| model.predict(f)).collect();
    crate::helpers::serialize_or_error(&serde_json::json!({
        "model_type": "random_forest",
        "n_estimators": model.n_trees,
        "n_samples": targets.len(),
        "train_rmse": rmse(&predictions, &targets),
        "predictions": predictions,
    }))
}

/// Train a Gradient Boosting regressor and return its training-set predictions.
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn train_gradient_boosting(
    features_json: &str,
    targets_json: &str,
    n_estimators: usize,
    learning_rate: f64,
    max_depth: usize,
    min_samples_split: usize,
    min_samples_leaf: usize,
    subsample: f64,
    seed: u64,
) -> String {
    let features: Vec<Vec<f64>> = match serde_json::from_str(features_json) {
        Ok(v) => v,
        Err(e) => {
            return crate::helpers::json_error(&crate::helpers::bad_json(
                "'features'",
                features_json,
                &e,
            ))
        }
    };
    let targets: Vec<f64> = match serde_json::from_str(targets_json) {
        Ok(v) => v,
        Err(e) => {
            return crate::helpers::json_error(&crate::helpers::bad_json(
                "'targets'",
                targets_json,
                &e,
            ))
        }
    };
    if features.len() != targets.len() {
        return crate::helpers::json_error(&format!(
            "features length {} != targets length {}",
            features.len(),
            targets.len()
        ));
    }
    let config = sci_form::ml::advanced_models::GradientBoostingConfig {
        n_estimators,
        learning_rate,
        tree_config: sci_form::ml::advanced_models::TreeConfig {
            max_depth,
            min_samples_split,
            min_samples_leaf,
            max_features: 0,
        },
        subsample,
        seed,
    };
    let model =
        sci_form::ml::advanced_models::train_gradient_boosting(&features, &targets, &config);
    let predictions: Vec<f64> = features.iter().map(|f| model.predict(f)).collect();
    crate::helpers::serialize_or_error(&serde_json::json!({
        "model_type": "gradient_boosting",
        "n_estimators": model.n_estimators,
        "n_samples": targets.len(),
        "train_rmse": rmse(&predictions, &targets),
        "predictions": predictions,
    }))
}
