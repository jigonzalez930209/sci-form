//! ML molecular descriptors and property prediction bindings.

use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
pub(crate) struct MolecularDescriptorsPy {
    #[pyo3(get)]
    molecular_weight: f64,
    #[pyo3(get)]
    n_heavy_atoms: usize,
    #[pyo3(get)]
    n_hydrogens: usize,
    #[pyo3(get)]
    n_bonds: usize,
    #[pyo3(get)]
    n_rotatable_bonds: usize,
    #[pyo3(get)]
    n_hbd: usize,
    #[pyo3(get)]
    n_hba: usize,
    #[pyo3(get)]
    fsp3: f64,
    #[pyo3(get)]
    wiener_index: f64,
    #[pyo3(get)]
    n_rings: usize,
    #[pyo3(get)]
    n_aromatic: usize,
    #[pyo3(get)]
    sum_electronegativity: f64,
    #[pyo3(get)]
    sum_polarizability: f64,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct MlPropertyResultPy {
    #[pyo3(get)]
    logp: f64,
    #[pyo3(get)]
    molar_refractivity: f64,
    #[pyo3(get)]
    log_solubility: f64,
    #[pyo3(get)]
    lipinski_violations: u8,
    #[pyo3(get)]
    lipinski_passes: bool,
    #[pyo3(get)]
    druglikeness: f64,
}

/// Elements and (i, j, order) bonds extracted from a SMILES string.
type MolElementsAndBonds = (Vec<u8>, Vec<(usize, usize, u8)>);

#[pyclass]
#[derive(Clone)]
pub(crate) struct WhimDescriptorsPy {
    #[pyo3(get)]
    eigenvalues: [f64; 3],
    #[pyo3(get)]
    lambda_: [f64; 3],
    #[pyo3(get)]
    nu: [f64; 3],
    #[pyo3(get)]
    gamma: [f64; 3],
    #[pyo3(get)]
    eta: [f64; 3],
    #[pyo3(get)]
    total_spread: f64,
    #[pyo3(get)]
    anisotropy: f64,
    #[pyo3(get)]
    directional: f64,
    #[pyo3(get)]
    kurtosis: f64,
}
#[pymethods]
impl WhimDescriptorsPy {
    fn __repr__(&self) -> String {
        format!(
            "WhimDescriptors(total_spread={:.4}, anisotropy={:.4}, kurtosis={:.4})",
            self.total_spread, self.anisotropy, self.kurtosis
        )
    }
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct RdfDescriptorsPy {
    #[pyo3(get)]
    radii: Vec<f64>,
    #[pyo3(get)]
    rdf_unweighted: Vec<f64>,
    #[pyo3(get)]
    rdf_mass: Vec<f64>,
    #[pyo3(get)]
    rdf_electronegativity: Vec<f64>,
    #[pyo3(get)]
    rdf_charge: Vec<f64>,
    #[pyo3(get)]
    beta: f64,
    #[pyo3(get)]
    n_pairs: usize,
}
#[pymethods]
impl RdfDescriptorsPy {
    fn __repr__(&self) -> String {
        format!(
            "RdfDescriptors(n_points={}, n_pairs={}, beta={:.2})",
            self.radii.len(),
            self.n_pairs,
            self.beta
        )
    }
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct GetawayDescriptorsPy {
    #[pyo3(get)]
    leverages: Vec<f64>,
    #[pyo3(get)]
    hat_autocorrelation: Vec<f64>,
    #[pyo3(get)]
    r_autocorrelation: Vec<f64>,
    #[pyo3(get)]
    hat_total: f64,
    #[pyo3(get)]
    r_total: f64,
    #[pyo3(get)]
    h_information: f64,
    #[pyo3(get)]
    h_max: f64,
    #[pyo3(get)]
    h_mean: f64,
    #[pyo3(get)]
    max_lag: usize,
}
#[pymethods]
impl GetawayDescriptorsPy {
    fn __repr__(&self) -> String {
        format!(
            "GetawayDescriptors(n_atoms={}, HATt={:.4}, Rt={:.4})",
            self.leverages.len(),
            self.hat_total,
            self.r_total
        )
    }
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct RegressionModelPy {
    /// Model predictions for the supplied feature rows.
    #[pyo3(get)]
    predictions: Vec<f64>,
    /// Root-mean-square error on the training set.
    #[pyo3(get)]
    train_rmse: f64,
    /// Model family ("random_forest" or "gradient_boosting").
    #[pyo3(get)]
    model_type: String,
    /// Number of trees / boosting stages.
    #[pyo3(get)]
    n_estimators: usize,
    /// Number of training samples.
    #[pyo3(get)]
    n_samples: usize,
}
#[pymethods]
impl RegressionModelPy {
    fn __repr__(&self) -> String {
        format!(
            "RegressionModel(type={}, n_estimators={}, train_rmse={:.6})",
            self.model_type, self.n_estimators, self.train_rmse
        )
    }
}

fn parse_mol_for_ml(smiles: &str) -> PyResult<MolElementsAndBonds> {
    let mol = sci_form_core::parse(smiles).map_err(pyo3::exceptions::PyValueError::new_err)?;
    let n = mol.graph.node_count();
    let elements: Vec<u8> = (0..n)
        .map(|i| mol.graph[sci_form_core::graph::NodeIndex::new(i)].element)
        .collect();
    let bonds: Vec<(usize, usize, u8)> = mol
        .graph
        .edge_indices()
        .map(|e| {
            let (a, b) = mol.graph.edge_endpoints(e).unwrap();
            let order = match mol.graph[e].order {
                sci_form_core::graph::BondOrder::Single => 1u8,
                sci_form_core::graph::BondOrder::Double => 2,
                sci_form_core::graph::BondOrder::Triple => 3,
                sci_form_core::graph::BondOrder::Aromatic => 2,
                sci_form_core::graph::BondOrder::Unknown => 1,
            };
            (a.index(), b.index(), order)
        })
        .collect();
    Ok((elements, bonds))
}

#[pyfunction]
fn ml_descriptors(smiles: &str) -> PyResult<MolecularDescriptorsPy> {
    let (elements, bonds) = parse_mol_for_ml(smiles)?;
    let desc = sci_form_core::compute_ml_descriptors(&elements, &bonds, &[], &[]);
    Ok(MolecularDescriptorsPy {
        molecular_weight: desc.molecular_weight,
        n_heavy_atoms: desc.n_heavy_atoms,
        n_hydrogens: desc.n_hydrogens,
        n_bonds: desc.n_bonds,
        n_rotatable_bonds: desc.n_rotatable_bonds,
        n_hbd: desc.n_hbd,
        n_hba: desc.n_hba,
        fsp3: desc.fsp3,
        wiener_index: desc.wiener_index,
        n_rings: desc.n_rings,
        n_aromatic: desc.n_aromatic,
        sum_electronegativity: desc.sum_electronegativity,
        sum_polarizability: desc.sum_polarizability,
    })
}

#[pyfunction]
fn ml_predict(smiles: &str) -> PyResult<MlPropertyResultPy> {
    let (elements, bonds) = parse_mol_for_ml(smiles)?;
    let desc = sci_form_core::compute_ml_descriptors(&elements, &bonds, &[], &[]);
    let result = sci_form_core::predict_ml_properties(&desc);
    Ok(MlPropertyResultPy {
        logp: result.logp,
        molar_refractivity: result.molar_refractivity,
        log_solubility: result.log_solubility,
        lipinski_violations: result.lipinski.violations,
        lipinski_passes: result.lipinski.passes,
        druglikeness: result.druglikeness,
    })
}

/// Compute WHIM 3D molecular descriptors.
///
/// `weighting` is one of `unit`, `mass`, `volume`, `electronegativity` or
/// `polarizability`.
#[pyfunction]
#[pyo3(signature = (elements, coords, weighting="mass"))]
fn compute_whim(
    elements: Vec<u8>,
    coords: Vec<f64>,
    weighting: &str,
) -> PyResult<WhimDescriptorsPy> {
    let positions = crate::system::coords_to_positions(&coords);
    let w = match weighting.to_ascii_lowercase().as_str() {
        "unit" | "unweighted" => sci_form_core::ml::whim::WhimWeighting::Unit,
        "mass" => sci_form_core::ml::whim::WhimWeighting::Mass,
        "volume" | "vdw" => sci_form_core::ml::whim::WhimWeighting::Volume,
        "electronegativity" | "en" => sci_form_core::ml::whim::WhimWeighting::Electronegativity,
        "polarizability" => sci_form_core::ml::whim::WhimWeighting::Polarizability,
        other => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unknown WHIM weighting '{}' (expected unit, mass, volume, \
                 electronegativity or polarizability)",
                other
            )))
        }
    };
    let d = sci_form_core::ml::whim::compute_whim(&elements, &positions, w);
    Ok(WhimDescriptorsPy {
        eigenvalues: d.eigenvalues,
        lambda_: d.lambda,
        nu: d.nu,
        gamma: d.gamma,
        eta: d.eta,
        total_spread: d.total_spread,
        anisotropy: d.anisotropy,
        directional: d.directional,
        kurtosis: d.kurtosis,
    })
}

/// Compute radial distribution function (RDF) descriptors.
///
/// `charges` is optional; when omitted the partial charges are derived with
/// the Gasteiger solver from geometry-perceived bonds.
#[pyfunction]
#[pyo3(signature = (elements, coords, charges=None, beta=100.0, r_max=12.0, n_points=128))]
fn compute_rdf(
    elements: Vec<u8>,
    coords: Vec<f64>,
    charges: Option<Vec<f64>>,
    beta: f64,
    r_max: f64,
    n_points: usize,
) -> PyResult<RdfDescriptorsPy> {
    let positions = crate::system::coords_to_positions(&coords);
    let q = match charges {
        Some(q) => q,
        None => crate::system::gasteiger_charges_from_geometry(&elements, &positions)
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?,
    };
    let d = sci_form_core::ml::rdf_descriptors::compute_rdf(
        &elements, &positions, &q, beta, r_max, n_points,
    );
    Ok(RdfDescriptorsPy {
        radii: d.radii,
        rdf_unweighted: d.rdf_unweighted,
        rdf_mass: d.rdf_mass,
        rdf_electronegativity: d.rdf_electronegativity,
        rdf_charge: d.rdf_charge,
        beta: d.beta,
        n_pairs: d.n_pairs,
    })
}

/// Compute GETAWAY 3D molecular descriptors.
///
/// `bonds` is an optional list of `[atom_i, atom_j]` index pairs; when
/// omitted the bonds are perceived from the geometry (covalent radii).
#[pyfunction]
#[pyo3(signature = (elements, coords, bonds=None, max_lag=8))]
fn compute_getaway(
    elements: Vec<u8>,
    coords: Vec<f64>,
    bonds: Option<Vec<Vec<usize>>>,
    max_lag: usize,
) -> PyResult<GetawayDescriptorsPy> {
    let positions = crate::system::coords_to_positions(&coords);
    let bond_list = match bonds {
        Some(b) => {
            let mut out = Vec::with_capacity(b.len());
            for (k, pair) in b.iter().enumerate() {
                if pair.len() != 2 {
                    return Err(pyo3::exceptions::PyValueError::new_err(format!(
                        "bonds[{}] must have exactly 2 atom indices, got {}",
                        k,
                        pair.len()
                    )));
                }
                out.push((pair[0], pair[1]));
            }
            out
        }
        None => crate::system::perceive_bonds(&elements, &positions, 0.3),
    };
    let d = sci_form_core::ml::getaway::compute_getaway(&elements, &positions, &bond_list, max_lag);
    Ok(GetawayDescriptorsPy {
        leverages: d.leverages,
        hat_autocorrelation: d.hat_autocorrelation,
        r_autocorrelation: d.r_autocorrelation,
        hat_total: d.hat_total,
        r_total: d.r_total,
        h_information: d.h_information,
        h_max: d.h_max,
        h_mean: d.h_mean,
        max_lag: d.max_lag,
    })
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
#[pyfunction]
#[pyo3(signature = (features, targets, n_trees=100, max_depth=10, min_samples_split=5, min_samples_leaf=2, sample_fraction=1.0, seed=42))]
fn train_random_forest(
    features: Vec<Vec<f64>>,
    targets: Vec<f64>,
    n_trees: usize,
    max_depth: usize,
    min_samples_split: usize,
    min_samples_leaf: usize,
    sample_fraction: f64,
    seed: u64,
) -> PyResult<RegressionModelPy> {
    if features.len() != targets.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "features length {} != targets length {}",
            features.len(),
            targets.len()
        )));
    }
    let config = sci_form_core::ml::advanced_models::RandomForestConfig {
        n_trees,
        tree_config: sci_form_core::ml::advanced_models::TreeConfig {
            max_depth,
            min_samples_split,
            min_samples_leaf,
            max_features: 0,
        },
        sample_fraction,
        seed,
    };
    let model =
        sci_form_core::ml::advanced_models::train_random_forest(&features, &targets, &config);
    let predictions: Vec<f64> = features.iter().map(|f| model.predict(f)).collect();
    Ok(RegressionModelPy {
        train_rmse: rmse(&predictions, &targets),
        predictions,
        model_type: "random_forest".to_string(),
        n_estimators: model.n_trees,
        n_samples: targets.len(),
    })
}

/// Train a Gradient Boosting regressor and return its training-set predictions.
#[pyfunction]
#[pyo3(signature = (features, targets, n_estimators=100, learning_rate=0.1, max_depth=5, min_samples_split=5, min_samples_leaf=2, subsample=0.8, seed=42))]
fn train_gradient_boosting(
    features: Vec<Vec<f64>>,
    targets: Vec<f64>,
    n_estimators: usize,
    learning_rate: f64,
    max_depth: usize,
    min_samples_split: usize,
    min_samples_leaf: usize,
    subsample: f64,
    seed: u64,
) -> PyResult<RegressionModelPy> {
    if features.len() != targets.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "features length {} != targets length {}",
            features.len(),
            targets.len()
        )));
    }
    let config = sci_form_core::ml::advanced_models::GradientBoostingConfig {
        n_estimators,
        learning_rate,
        tree_config: sci_form_core::ml::advanced_models::TreeConfig {
            max_depth,
            min_samples_split,
            min_samples_leaf,
            max_features: 0,
        },
        subsample,
        seed,
    };
    let model =
        sci_form_core::ml::advanced_models::train_gradient_boosting(&features, &targets, &config);
    let predictions: Vec<f64> = features.iter().map(|f| model.predict(f)).collect();
    Ok(RegressionModelPy {
        train_rmse: rmse(&predictions, &targets),
        predictions,
        model_type: "gradient_boosting".to_string(),
        n_estimators: model.n_estimators,
        n_samples: targets.len(),
    })
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(ml_descriptors, m)?)?;
    m.add_function(wrap_pyfunction!(ml_predict, m)?)?;
    m.add_function(wrap_pyfunction!(compute_whim, m)?)?;
    m.add_function(wrap_pyfunction!(compute_rdf, m)?)?;
    m.add_function(wrap_pyfunction!(compute_getaway, m)?)?;
    m.add_function(wrap_pyfunction!(train_random_forest, m)?)?;
    m.add_function(wrap_pyfunction!(train_gradient_boosting, m)?)?;
    m.add_class::<MolecularDescriptorsPy>()?;
    m.add_class::<MlPropertyResultPy>()?;
    m.add_class::<WhimDescriptorsPy>()?;
    m.add_class::<RdfDescriptorsPy>()?;
    m.add_class::<GetawayDescriptorsPy>()?;
    m.add_class::<RegressionModelPy>()?;
    Ok(())
}
