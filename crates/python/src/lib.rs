use pyo3::prelude::*;
use pyo3::types::PyDict;

/// A single conformer result returned to Python.
#[pyclass]
#[derive(Clone)]
struct ConformerResult {
    #[pyo3(get)]
    smiles: String,
    #[pyo3(get)]
    num_atoms: usize,
    #[pyo3(get)]
    coords: Vec<f64>,
    #[pyo3(get)]
    elements: Vec<u8>,
    #[pyo3(get)]
    bonds: Vec<(usize, usize, String)>,
    #[pyo3(get)]
    error: Option<String>,
    #[pyo3(get)]
    time_ms: f64,
}

#[pymethods]
impl ConformerResult {
    /// Get coordinates as list of (x, y, z) tuples.
    fn get_positions(&self) -> Vec<(f64, f64, f64)> {
        self.coords
            .chunks_exact(3)
            .map(|c| (c[0], c[1], c[2]))
            .collect()
    }

    /// Check if generation succeeded.
    fn is_ok(&self) -> bool {
        self.error.is_none()
    }

    fn __repr__(&self) -> String {
        if let Some(ref e) = self.error {
            format!("ConformerResult(smiles='{}', error='{}')", self.smiles, e)
        } else {
            format!(
                "ConformerResult(smiles='{}', atoms={}, time={:.1}ms)",
                self.smiles, self.num_atoms, self.time_ms
            )
        }
    }
}

impl From<sci_form_core::ConformerResult> for ConformerResult {
    fn from(r: sci_form_core::ConformerResult) -> Self {
        ConformerResult {
            smiles: r.smiles,
            num_atoms: r.num_atoms,
            coords: r.coords,
            elements: r.elements,
            bonds: r.bonds,
            error: r.error,
            time_ms: r.time_ms,
        }
    }
}

/// Generate a 3D conformer from a SMILES string.
///
/// Args:
///     smiles: SMILES string of the molecule.
///     seed: RNG seed for reproducibility (default: 42).
///
/// Returns:
///     ConformerResult with 3D coordinates.
#[pyfunction]
#[pyo3(signature = (smiles, seed=42))]
fn embed(smiles: &str, seed: u64) -> ConformerResult {
    sci_form_core::embed(smiles, seed).into()
}

/// Batch-generate 3D conformers for multiple SMILES in parallel.
///
/// Args:
///     smiles_list: List of SMILES strings.
///     seed: RNG seed (default: 42).
///     num_threads: Number of threads, 0 = auto (default: 0).
///
/// Returns:
///     List of ConformerResult objects.
#[pyfunction]
#[pyo3(signature = (smiles_list, seed=42, num_threads=0))]
fn embed_batch(smiles_list: Vec<String>, seed: u64, num_threads: usize) -> Vec<ConformerResult> {
    let refs: Vec<&str> = smiles_list.iter().map(|s| s.as_str()).collect();
    let config = sci_form_core::ConformerConfig { seed, num_threads };
    sci_form_core::embed_batch(&refs, &config)
        .into_iter()
        .map(|r| r.into())
        .collect()
}

/// Parse a SMILES string without generating 3D coordinates.
///
/// Returns a dict with atom/bond info.
#[pyfunction]
fn parse(py: Python<'_>, smiles: &str) -> PyResult<PyObject> {
    match sci_form_core::parse(smiles) {
        Ok(mol) => {
            let dict = PyDict::new_bound(py);
            let n = mol.graph.node_count();
            dict.set_item("num_atoms", n)?;
            dict.set_item("num_bonds", mol.graph.edge_count())?;

            let mut atoms_list: Vec<PyObject> = Vec::new();
            for i in 0..n {
                let idx = sci_form_core::graph::NodeIndex::new(i);
                let atom = &mol.graph[idx];
                let adict = PyDict::new_bound(py);
                adict.set_item("element", atom.element)?;
                adict.set_item("hybridization", format!("{:?}", atom.hybridization))?;
                adict.set_item("formal_charge", atom.formal_charge)?;
                atoms_list.push(adict.into());
            }
            dict.set_item("atoms", atoms_list)?;
            Ok(dict.into())
        }
        Err(e) => Err(pyo3::exceptions::PyValueError::new_err(e)),
    }
}

/// Get the library version string.
#[pyfunction]
fn version() -> String {
    sci_form_core::version()
}

/// sci_form Python module
#[pymodule]
fn sci_form(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(embed, m)?)?;
    m.add_function(wrap_pyfunction!(embed_batch, m)?)?;
    m.add_function(wrap_pyfunction!(parse, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(system_capabilities, m)?)?;
    m.add_function(wrap_pyfunction!(system_method_plan, m)?)?;
    m.add_function(wrap_pyfunction!(compare_methods, m)?)?;
    m.add_function(wrap_pyfunction!(eht_support, m)?)?;
    m.add_function(wrap_pyfunction!(eht_calculate, m)?)?;
    m.add_function(wrap_pyfunction!(eht_or_uff_fallback, m)?)?;
    m.add_function(wrap_pyfunction!(eht_orbital_mesh, m)?)?;
    m.add_function(wrap_pyfunction!(charges, m)?)?;
    m.add_function(wrap_pyfunction!(sasa, m)?)?;
    m.add_function(wrap_pyfunction!(population, m)?)?;
    m.add_function(wrap_pyfunction!(bond_orders, m)?)?;
    m.add_function(wrap_pyfunction!(frontier_descriptors, m)?)?;
    m.add_function(wrap_pyfunction!(fukui_descriptors, m)?)?;
    m.add_function(wrap_pyfunction!(reactivity_ranking, m)?)?;
    m.add_function(wrap_pyfunction!(uvvis_spectrum, m)?)?;
    m.add_function(wrap_pyfunction!(graph_features, m)?)?;
    m.add_function(wrap_pyfunction!(topology_analysis, m)?)?;
    m.add_function(wrap_pyfunction!(dipole, m)?)?;
    m.add_function(wrap_pyfunction!(dos, m)?)?;
    m.add_function(wrap_pyfunction!(rmsd, m)?)?;
    m.add_function(wrap_pyfunction!(uff_energy, m)?)?;
    m.add_function(wrap_pyfunction!(uff_energy_with_aromatic_heuristics, m)?)?;
    m.add_function(wrap_pyfunction!(empirical_pka, m)?)?;
    m.add_function(wrap_pyfunction!(unit_cell, m)?)?;
    m.add_function(wrap_pyfunction!(assemble_framework, m)?)?;
    m.add_function(wrap_pyfunction!(pack_conformers, m)?)?;
    m.add_function(wrap_pyfunction!(split_worker_tasks, m)?)?;
    m.add_function(wrap_pyfunction!(estimate_workers, m)?)?;
    m.add_function(wrap_pyfunction!(mmff94_energy, m)?)?;
    m.add_function(wrap_pyfunction!(pm3_calculate, m)?)?;
    m.add_function(wrap_pyfunction!(xtb_calculate, m)?)?;
    m.add_function(wrap_pyfunction!(ml_descriptors, m)?)?;
    m.add_function(wrap_pyfunction!(ml_predict, m)?)?;
    m.add_function(wrap_pyfunction!(stda_uvvis, m)?)?;
    m.add_function(wrap_pyfunction!(vibrational_analysis, m)?)?;
    m.add_function(wrap_pyfunction!(ir_spectrum, m)?)?;
    m.add_function(wrap_pyfunction!(nmr_shifts, m)?)?;
    m.add_function(wrap_pyfunction!(nmr_couplings, m)?)?;
    m.add_function(wrap_pyfunction!(nmr_spectrum, m)?)?;
    m.add_function(wrap_pyfunction!(hose_codes, m)?)?;
    m.add_class::<ConformerResult>()?;
    m.add_class::<SystemCapabilitiesPy>()?;
    m.add_class::<MethodMetadataPy>()?;
    m.add_class::<PropertyMethodPlanPy>()?;
    m.add_class::<SystemMethodPlanPy>()?;
    m.add_class::<MethodComparisonEntryPy>()?;
    m.add_class::<MethodComparisonResultPy>()?;
    m.add_class::<EhtResultPy>()?;
    m.add_class::<EhtSupportPy>()?;
    m.add_class::<ChargeResultPy>()?;
    m.add_class::<SasaResultPy>()?;
    m.add_class::<PopulationResultPy>()?;
    m.add_class::<BondOrderResultPy>()?;
    m.add_class::<FrontierDescriptorsPy>()?;
    m.add_class::<FukuiDescriptorsPy>()?;
    m.add_class::<ReactivitySiteScorePy>()?;
    m.add_class::<ReactivityRankingPy>()?;
    m.add_class::<UvVisPeakPy>()?;
    m.add_class::<UvVisSpectrumPy>()?;
    m.add_class::<GraphFeatureAnalysisPy>()?;
    m.add_class::<MetalCoordinationCenterPy>()?;
    m.add_class::<TopologyAnalysisResultPy>()?;
    m.add_class::<DipoleResultPy>()?;
    m.add_class::<DosResultPy>()?;
    m.add_class::<AlignmentResultPy>()?;
    m.add_class::<UffHeuristicEnergyPy>()?;
    m.add_class::<EmpiricalPkaSitePy>()?;
    m.add_class::<EmpiricalPkaResultPy>()?;
    m.add_class::<UnitCellPy>()?;
    m.add_class::<CrystalStructurePy>()?;
    m.add_class::<RecordBatchPy>()?;
    m.add_class::<Pm3ResultPy>()?;
    m.add_class::<XtbResultPy>()?;
    m.add_class::<MolecularDescriptorsPy>()?;
    m.add_class::<MlPropertyResultPy>()?;
    m.add_class::<StdaExcitationPy>()?;
    m.add_class::<StdaUvVisSpectrumPy>()?;
    m.add_class::<VibrationalModePy>()?;
    m.add_class::<VibrationalAnalysisPy>()?;
    m.add_class::<IrPeakPy>()?;
    m.add_class::<IrSpectrumPy>()?;
    m.add_class::<ChemicalShiftPy>()?;
    m.add_class::<NmrShiftResultPy>()?;
    m.add_class::<JCouplingPy>()?;
    m.add_class::<NmrPeakPy>()?;
    m.add_class::<NmrSpectrumPy>()?;
    Ok(())
}

#[pyclass]
#[derive(Clone)]
struct MethodCapabilityPy {
    #[pyo3(get)]
    available: bool,
    #[pyo3(get)]
    confidence: String,
    #[pyo3(get)]
    unsupported_elements: Vec<u8>,
    #[pyo3(get)]
    warnings: Vec<String>,
}

impl From<sci_form_core::MethodCapability> for MethodCapabilityPy {
    fn from(value: sci_form_core::MethodCapability) -> Self {
        Self {
            available: value.available,
            confidence: format!("{:?}", value.confidence).to_lowercase(),
            unsupported_elements: value.unsupported_elements,
            warnings: value.warnings,
        }
    }
}

fn scientific_method_name(method: sci_form_core::ScientificMethod) -> String {
    match method {
        sci_form_core::ScientificMethod::Embed => "embed",
        sci_form_core::ScientificMethod::Uff => "uff",
        sci_form_core::ScientificMethod::Eht => "eht",
        sci_form_core::ScientificMethod::Pm3 => "pm3",
        sci_form_core::ScientificMethod::Xtb => "xtb",
        sci_form_core::ScientificMethod::Mmff94 => "mmff94",
    }
    .to_string()
}

fn property_request_name(property: sci_form_core::PropertyRequest) -> String {
    match property {
        sci_form_core::PropertyRequest::Geometry => "geometry",
        sci_form_core::PropertyRequest::ForceFieldEnergy => "force_field_energy",
        sci_form_core::PropertyRequest::Orbitals => "orbitals",
        sci_form_core::PropertyRequest::Population => "population",
        sci_form_core::PropertyRequest::OrbitalGrid => "orbital_grid",
    }
    .to_string()
}

fn comparison_status_name(status: sci_form_core::MethodComparisonStatus) -> String {
    match status {
        sci_form_core::MethodComparisonStatus::Success => "success",
        sci_form_core::MethodComparisonStatus::Unavailable => "unavailable",
        sci_form_core::MethodComparisonStatus::Error => "error",
    }
    .to_string()
}

#[pyclass]
#[derive(Clone)]
struct SystemCapabilitiesPy {
    #[pyo3(get)]
    embed: MethodCapabilityPy,
    #[pyo3(get)]
    uff: MethodCapabilityPy,
    #[pyo3(get)]
    eht: MethodCapabilityPy,
    #[pyo3(get)]
    population: MethodCapabilityPy,
    #[pyo3(get)]
    orbital_grid: MethodCapabilityPy,
}

#[pyclass]
#[derive(Clone)]
struct MethodMetadataPy {
    #[pyo3(get)]
    method: String,
    #[pyo3(get)]
    available: bool,
    #[pyo3(get)]
    confidence: String,
    #[pyo3(get)]
    confidence_score: f64,
    #[pyo3(get)]
    limitations: Vec<String>,
    #[pyo3(get)]
    warnings: Vec<String>,
}

impl From<sci_form_core::MethodMetadata> for MethodMetadataPy {
    fn from(value: sci_form_core::MethodMetadata) -> Self {
        Self {
            method: scientific_method_name(value.method),
            available: value.available,
            confidence: format!("{:?}", value.confidence).to_lowercase(),
            confidence_score: value.confidence_score,
            limitations: value.limitations,
            warnings: value.warnings,
        }
    }
}

#[pyclass]
#[derive(Clone)]
struct PropertyMethodPlanPy {
    #[pyo3(get)]
    property: String,
    #[pyo3(get)]
    recommended: Option<String>,
    #[pyo3(get)]
    fallback: Option<String>,
    #[pyo3(get)]
    rationale: Vec<String>,
    #[pyo3(get)]
    methods: Vec<MethodMetadataPy>,
}

impl From<sci_form_core::PropertyMethodPlan> for PropertyMethodPlanPy {
    fn from(value: sci_form_core::PropertyMethodPlan) -> Self {
        Self {
            property: property_request_name(value.property),
            recommended: value.recommended.map(scientific_method_name),
            fallback: value.fallback.map(scientific_method_name),
            rationale: value.rationale,
            methods: value.methods.into_iter().map(Into::into).collect(),
        }
    }
}

#[pyclass]
#[derive(Clone)]
struct SystemMethodPlanPy {
    #[pyo3(get)]
    capabilities: SystemCapabilitiesPy,
    #[pyo3(get)]
    geometry: PropertyMethodPlanPy,
    #[pyo3(get)]
    force_field_energy: PropertyMethodPlanPy,
    #[pyo3(get)]
    orbitals: PropertyMethodPlanPy,
    #[pyo3(get)]
    population: PropertyMethodPlanPy,
    #[pyo3(get)]
    orbital_grid: PropertyMethodPlanPy,
}

#[pyclass]
#[derive(Clone)]
struct MethodComparisonEntryPy {
    #[pyo3(get)]
    method: String,
    #[pyo3(get)]
    status: String,
    #[pyo3(get)]
    available: bool,
    #[pyo3(get)]
    confidence: String,
    #[pyo3(get)]
    confidence_score: f64,
    #[pyo3(get)]
    warnings: Vec<String>,
    #[pyo3(get)]
    limitations: Vec<String>,
    #[pyo3(get)]
    error: Option<String>,
    #[pyo3(get)]
    support_level: Option<String>,
    #[pyo3(get)]
    homo_energy: Option<f64>,
    #[pyo3(get)]
    lumo_energy: Option<f64>,
    #[pyo3(get)]
    gap: Option<f64>,
    #[pyo3(get)]
    energy_kcal_mol: Option<f64>,
}

#[pyclass]
#[derive(Clone)]
struct MethodComparisonResultPy {
    #[pyo3(get)]
    plan: SystemMethodPlanPy,
    #[pyo3(get)]
    comparisons: Vec<MethodComparisonEntryPy>,
}

/// Query operation capabilities for a given element set.
#[pyfunction]
fn system_capabilities(elements: Vec<u8>) -> SystemCapabilitiesPy {
    let caps = sci_form_core::get_system_capabilities(&elements);
    SystemCapabilitiesPy {
        embed: caps.embed.into(),
        uff: caps.uff.into(),
        eht: caps.eht.into(),
        population: caps.population.into(),
        orbital_grid: caps.orbital_grid.into(),
    }
}

/// Build a structured method plan with recommendations, fallbacks, and confidence scores.
#[pyfunction]
fn system_method_plan(elements: Vec<u8>) -> SystemMethodPlanPy {
    let plan = sci_form_core::get_system_method_plan(&elements);
    SystemMethodPlanPy {
        capabilities: SystemCapabilitiesPy {
            embed: plan.capabilities.embed.into(),
            uff: plan.capabilities.uff.into(),
            eht: plan.capabilities.eht.into(),
            population: plan.capabilities.population.into(),
            orbital_grid: plan.capabilities.orbital_grid.into(),
        },
        geometry: plan.geometry.into(),
        force_field_energy: plan.force_field_energy.into(),
        orbitals: plan.orbitals.into(),
        population: plan.population.into(),
        orbital_grid: plan.orbital_grid.into(),
    }
}

/// Compare available methods on the same system/geometry.
#[pyfunction]
#[pyo3(signature = (smiles, elements, coords, allow_experimental_eht=false))]
fn compare_methods(
    smiles: &str,
    elements: Vec<u8>,
    coords: Vec<f64>,
    allow_experimental_eht: bool,
) -> PyResult<MethodComparisonResultPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "coords length must be 3 * len(elements)",
        ));
    }
    let positions: Vec<[f64; 3]> = coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    let result =
        sci_form_core::compare_methods(smiles, &elements, &positions, allow_experimental_eht);

    let comparisons = result
        .comparisons
        .into_iter()
        .map(|entry| {
            let mut support_level = None;
            let mut homo_energy = None;
            let mut lumo_energy = None;
            let mut gap = None;
            let mut energy_kcal_mol = None;

            if let Some(payload) = entry.payload {
                match payload {
                    sci_form_core::MethodComparisonPayload::Eht {
                        homo_energy: h,
                        lumo_energy: l,
                        gap: g,
                        support,
                    } => {
                        support_level = Some(format!("{:?}", support.level).to_lowercase());
                        homo_energy = Some(h);
                        lumo_energy = Some(l);
                        gap = Some(g);
                    }
                    sci_form_core::MethodComparisonPayload::Uff { energy_kcal_mol: e } => {
                        energy_kcal_mol = Some(e);
                    }
                }
            }

            MethodComparisonEntryPy {
                method: scientific_method_name(entry.method),
                status: comparison_status_name(entry.status),
                available: entry.available,
                confidence: format!("{:?}", entry.confidence).to_lowercase(),
                confidence_score: entry.confidence_score,
                warnings: entry.warnings,
                limitations: entry.limitations,
                error: entry.error,
                support_level,
                homo_energy,
                lumo_energy,
                gap,
                energy_kcal_mol,
            }
        })
        .collect();

    Ok(MethodComparisonResultPy {
        plan: SystemMethodPlanPy {
            capabilities: SystemCapabilitiesPy {
                embed: result.plan.capabilities.embed.into(),
                uff: result.plan.capabilities.uff.into(),
                eht: result.plan.capabilities.eht.into(),
                population: result.plan.capabilities.population.into(),
                orbital_grid: result.plan.capabilities.orbital_grid.into(),
            },
            geometry: result.plan.geometry.into(),
            force_field_energy: result.plan.force_field_energy.into(),
            orbitals: result.plan.orbitals.into(),
            population: result.plan.population.into(),
            orbital_grid: result.plan.orbital_grid.into(),
        },
        comparisons,
    })
}

// ─── EHT (Extended Hückel Theory) API ────────────────────────────────────────

/// EHT calculation result.
#[pyclass]
#[derive(Clone)]
struct EhtResultPy {
    #[pyo3(get)]
    energies: Vec<f64>,
    #[pyo3(get)]
    n_electrons: usize,
    #[pyo3(get)]
    homo_index: usize,
    #[pyo3(get)]
    lumo_index: usize,
    #[pyo3(get)]
    homo_energy: f64,
    #[pyo3(get)]
    lumo_energy: f64,
    #[pyo3(get)]
    gap: f64,
    #[pyo3(get)]
    support_level: String,
    #[pyo3(get)]
    has_transition_metals: bool,
    #[pyo3(get)]
    supported_elements: Vec<u8>,
    #[pyo3(get)]
    provisional_elements: Vec<u8>,
    #[pyo3(get)]
    unsupported_elements: Vec<u8>,
    #[pyo3(get)]
    warnings: Vec<String>,
}

#[pymethods]
impl EhtResultPy {
    fn __repr__(&self) -> String {
        format!(
            "EhtResult(n_mo={}, gap={:.3} eV, HOMO={:.3} eV, LUMO={:.3} eV)",
            self.energies.len(),
            self.gap,
            self.homo_energy,
            self.lumo_energy,
        )
    }
}

/// EHT support and confidence metadata for an element set.
#[pyclass]
#[derive(Clone)]
struct EhtSupportPy {
    #[pyo3(get)]
    level: String,
    #[pyo3(get)]
    has_transition_metals: bool,
    #[pyo3(get)]
    supported_elements: Vec<u8>,
    #[pyo3(get)]
    provisional_elements: Vec<u8>,
    #[pyo3(get)]
    unsupported_elements: Vec<u8>,
    #[pyo3(get)]
    warnings: Vec<String>,
}

impl From<sci_form_core::eht::EhtSupport> for EhtSupportPy {
    fn from(value: sci_form_core::eht::EhtSupport) -> Self {
        Self {
            level: format!("{:?}", value.level).to_lowercase(),
            has_transition_metals: value.has_transition_metals,
            supported_elements: value.supported_elements,
            provisional_elements: value.provisional_elements,
            unsupported_elements: value.unsupported_elements,
            warnings: value.warnings,
        }
    }
}

/// Run an EHT calculation on a molecule.
///
/// Args:
///     elements: List of atomic numbers (e.g. [8, 1, 1] for water).
///     coords: Flat list of xyz coordinates [x0,y0,z0, x1,y1,z1, ...].
///     k: Wolfsberg-Helmholtz constant (default: 1.75).
///
/// Returns:
///     EhtResult with orbital energies, HOMO/LUMO info.
#[pyfunction]
#[pyo3(signature = (elements, coords, k=1.75))]
fn eht_calculate(elements: Vec<u8>, coords: Vec<f64>, k: f64) -> PyResult<EhtResultPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "coords length must be 3 * len(elements)",
        ));
    }
    let positions: Vec<[f64; 3]> = coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    let k_opt = if (k - 1.75).abs() < 1e-10 {
        None
    } else {
        Some(k)
    };

    match sci_form_core::eht::solve_eht(&elements, &positions, k_opt) {
        Ok(r) => Ok(EhtResultPy {
            energies: r.energies,
            n_electrons: r.n_electrons,
            homo_index: r.homo_index,
            lumo_index: r.lumo_index,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap: r.gap,
            support_level: format!("{:?}", r.support.level).to_lowercase(),
            has_transition_metals: r.support.has_transition_metals,
            supported_elements: r.support.supported_elements,
            provisional_elements: r.support.provisional_elements,
            unsupported_elements: r.support.unsupported_elements,
            warnings: r.support.warnings,
        }),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

/// Run EHT or route to UFF fallback based on support/confidence.
#[pyfunction]
#[pyo3(signature = (smiles, elements, coords, allow_experimental_eht=false))]
fn eht_or_uff_fallback(
    py: Python<'_>,
    smiles: &str,
    elements: Vec<u8>,
    coords: Vec<f64>,
    allow_experimental_eht: bool,
) -> PyResult<PyObject> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "coords length must be 3 * len(elements)",
        ));
    }
    let positions: Vec<[f64; 3]> = coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    let result = sci_form_core::compute_eht_or_uff_fallback(
        smiles,
        &elements,
        &positions,
        allow_experimental_eht,
    )
    .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;

    let dict = PyDict::new_bound(py);
    match result {
        sci_form_core::ElectronicWorkflowResult::Eht { result } => {
            dict.set_item("mode", "eht")?;
            dict.set_item("energies", result.energies)?;
            dict.set_item("n_electrons", result.n_electrons)?;
            dict.set_item("homo_index", result.homo_index)?;
            dict.set_item("lumo_index", result.lumo_index)?;
            dict.set_item("homo_energy", result.homo_energy)?;
            dict.set_item("lumo_energy", result.lumo_energy)?;
            dict.set_item("gap", result.gap)?;
            dict.set_item(
                "support_level",
                format!("{:?}", result.support.level).to_lowercase(),
            )?;
            dict.set_item("warnings", result.support.warnings)?;
        }
        sci_form_core::ElectronicWorkflowResult::UffFallback {
            energy_kcal_mol,
            reason,
            support,
        } => {
            dict.set_item("mode", "uff_fallback")?;
            dict.set_item("energy_kcal_mol", energy_kcal_mol)?;
            dict.set_item("reason", reason)?;
            dict.set_item(
                "support_level",
                format!("{:?}", support.level).to_lowercase(),
            )?;
            dict.set_item("warnings", support.warnings)?;
        }
    }
    Ok(dict.into())
}

/// Query EHT support metadata before running a calculation.
#[pyfunction]
fn eht_support(elements: Vec<u8>) -> EhtSupportPy {
    sci_form_core::get_eht_support(&elements).into()
}

/// Generate an orbital isosurface mesh as a dict.
///
/// Args:
///     elements: List of atomic numbers.
///     coords: Flat xyz coordinates.
///     mo_index: Molecular orbital index (0-based).
///     spacing: Grid spacing in Ångström (default: 0.2).
///     isovalue: Isosurface cutoff (default: 0.02).
///
/// Returns:
///     Dict with 'vertices', 'normals', 'indices', 'num_triangles'.
#[pyfunction]
#[pyo3(signature = (elements, coords, mo_index, spacing=0.2, isovalue=0.02))]
fn eht_orbital_mesh(
    py: Python<'_>,
    elements: Vec<u8>,
    coords: Vec<f64>,
    mo_index: usize,
    spacing: f64,
    isovalue: f32,
) -> PyResult<PyObject> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "coords length must be 3 * len(elements)",
        ));
    }
    let positions: Vec<[f64; 3]> = coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();

    let result = sci_form_core::eht::solve_eht(&elements, &positions, None)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e))?;

    let basis = sci_form_core::eht::basis::build_basis(&elements, &positions);
    let grid = sci_form_core::eht::evaluate_orbital_on_grid(
        &basis,
        &result.coefficients,
        mo_index,
        &positions,
        spacing,
        3.0,
    );
    let mesh = sci_form_core::eht::marching_cubes(&grid, isovalue);

    let dict = PyDict::new_bound(py);
    dict.set_item("vertices", mesh.vertices)?;
    dict.set_item("normals", mesh.normals)?;
    dict.set_item("indices", mesh.indices)?;
    dict.set_item("num_triangles", mesh.num_triangles)?;
    dict.set_item("isovalue", mesh.isovalue)?;
    Ok(dict.into())
}

// ─── Charges (Gasteiger-Marsili) API ─────────────────────────────────────────

/// Gasteiger-Marsili charge result.
#[pyclass]
#[derive(Clone)]
struct ChargeResultPy {
    #[pyo3(get)]
    charges: Vec<f64>,
    #[pyo3(get)]
    iterations: usize,
    #[pyo3(get)]
    total_charge: f64,
}

#[pymethods]
impl ChargeResultPy {
    fn __repr__(&self) -> String {
        format!(
            "ChargeResult(n_atoms={}, total_charge={:.4}, iterations={})",
            self.charges.len(),
            self.total_charge,
            self.iterations,
        )
    }
}

/// Compute Gasteiger-Marsili partial charges from a SMILES string.
///
/// Args:
///     smiles: SMILES string of the molecule.
///
/// Returns:
///     ChargeResult with partial charges per atom.
#[pyfunction]
fn charges(smiles: &str) -> PyResult<ChargeResultPy> {
    match sci_form_core::compute_charges(smiles) {
        Ok(r) => Ok(ChargeResultPy {
            charges: r.charges,
            iterations: r.iterations,
            total_charge: r.total_charge,
        }),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

// ─── SASA (Shrake-Rupley) API ────────────────────────────────────────────────

/// SASA calculation result.
#[pyclass]
#[derive(Clone)]
struct SasaResultPy {
    #[pyo3(get)]
    total_sasa: f64,
    #[pyo3(get)]
    atom_sasa: Vec<f64>,
    #[pyo3(get)]
    probe_radius: f64,
    #[pyo3(get)]
    num_points: usize,
}

#[pymethods]
impl SasaResultPy {
    fn __repr__(&self) -> String {
        format!(
            "SasaResult(total={:.2} Ų, n_atoms={}, probe={:.2} Å)",
            self.total_sasa,
            self.atom_sasa.len(),
            self.probe_radius,
        )
    }
}

/// Compute solvent-accessible surface area.
///
/// Args:
///     elements: List of atomic numbers.
///     coords: Flat list of xyz coordinates [x0,y0,z0, x1,y1,z1, ...].
///     probe_radius: Probe radius in Å (default: 1.4).
///
/// Returns:
///     SasaResult with total and per-atom SASA values.
#[pyfunction]
#[pyo3(signature = (elements, coords, probe_radius=1.4))]
fn sasa(elements: Vec<u8>, coords: Vec<f64>, probe_radius: f64) -> PyResult<SasaResultPy> {
    let pr = if (probe_radius - 1.4).abs() < 1e-10 {
        None
    } else {
        Some(probe_radius)
    };
    match sci_form_core::compute_sasa(&elements, &coords, pr) {
        Ok(r) => Ok(SasaResultPy {
            total_sasa: r.total_sasa,
            atom_sasa: r.atom_sasa,
            probe_radius: r.probe_radius,
            num_points: r.num_points,
        }),
        Err(e) => Err(pyo3::exceptions::PyValueError::new_err(e)),
    }
}

// ─── Population Analysis (Mulliken & Löwdin) API ─────────────────────────────

#[pyclass]
#[derive(Clone)]
struct PopulationResultPy {
    #[pyo3(get)]
    mulliken_charges: Vec<f64>,
    #[pyo3(get)]
    lowdin_charges: Vec<f64>,
    #[pyo3(get)]
    num_atoms: usize,
    #[pyo3(get)]
    total_charge_mulliken: f64,
    #[pyo3(get)]
    total_charge_lowdin: f64,
}

#[pymethods]
impl PopulationResultPy {
    fn __repr__(&self) -> String {
        format!(
            "PopulationResult(n_atoms={}, total_mulliken={:.4})",
            self.num_atoms, self.total_charge_mulliken
        )
    }
}

#[pyclass]
#[derive(Clone)]
struct BondOrderResultPy {
    #[pyo3(get)]
    atom_pairs: Vec<(usize, usize)>,
    #[pyo3(get)]
    distances: Vec<f64>,
    #[pyo3(get)]
    wiberg: Vec<f64>,
    #[pyo3(get)]
    mayer: Vec<f64>,
    #[pyo3(get)]
    wiberg_valence: Vec<f64>,
    #[pyo3(get)]
    mayer_valence: Vec<f64>,
    #[pyo3(get)]
    num_atoms: usize,
}

#[pymethods]
impl BondOrderResultPy {
    fn __repr__(&self) -> String {
        format!(
            "BondOrderResult(n_atoms={}, n_pairs={})",
            self.num_atoms,
            self.atom_pairs.len()
        )
    }
}

#[pyclass]
#[derive(Clone)]
struct FrontierDescriptorsPy {
    #[pyo3(get)]
    homo_atom_contributions: Vec<f64>,
    #[pyo3(get)]
    lumo_atom_contributions: Vec<f64>,
    #[pyo3(get)]
    dual_descriptor: Vec<f64>,
    #[pyo3(get)]
    homo_energy: f64,
    #[pyo3(get)]
    lumo_energy: f64,
    #[pyo3(get)]
    gap: f64,
    #[pyo3(get)]
    num_atoms: usize,
}

#[pyclass]
#[derive(Clone)]
struct FukuiDescriptorsPy {
    #[pyo3(get)]
    f_plus: Vec<f64>,
    #[pyo3(get)]
    f_minus: Vec<f64>,
    #[pyo3(get)]
    f_radical: Vec<f64>,
    #[pyo3(get)]
    dual_descriptor: Vec<f64>,
    #[pyo3(get)]
    condensed_atom_indices: Vec<usize>,
    #[pyo3(get)]
    condensed_f_plus: Vec<f64>,
    #[pyo3(get)]
    condensed_f_minus: Vec<f64>,
    #[pyo3(get)]
    condensed_f_radical: Vec<f64>,
    #[pyo3(get)]
    condensed_dual_descriptor: Vec<f64>,
    #[pyo3(get)]
    homo_energy: f64,
    #[pyo3(get)]
    lumo_energy: f64,
    #[pyo3(get)]
    gap: f64,
    #[pyo3(get)]
    num_atoms: usize,
    #[pyo3(get)]
    validity_notes: Vec<String>,
}

#[pymethods]
impl FukuiDescriptorsPy {
    fn __repr__(&self) -> String {
        format!(
            "FukuiDescriptors(n_atoms={}, gap={:.3} eV)",
            self.num_atoms, self.gap
        )
    }
}

#[pyclass]
#[derive(Clone)]
struct ReactivitySiteScorePy {
    #[pyo3(get)]
    atom_index: usize,
    #[pyo3(get)]
    score: f64,
}

#[pyclass]
#[derive(Clone)]
struct ReactivityRankingPy {
    #[pyo3(get)]
    nucleophilic_attack_sites: Vec<ReactivitySiteScorePy>,
    #[pyo3(get)]
    electrophilic_attack_sites: Vec<ReactivitySiteScorePy>,
    #[pyo3(get)]
    radical_attack_sites: Vec<ReactivitySiteScorePy>,
    #[pyo3(get)]
    notes: Vec<String>,
}

#[pymethods]
impl ReactivityRankingPy {
    fn __repr__(&self) -> String {
        format!(
            "ReactivityRanking(nuc={}, elec={}, rad={})",
            self.nucleophilic_attack_sites.len(),
            self.electrophilic_attack_sites.len(),
            self.radical_attack_sites.len()
        )
    }
}

#[pyclass]
#[derive(Clone)]
struct UvVisPeakPy {
    #[pyo3(get)]
    energy_ev: f64,
    #[pyo3(get)]
    wavelength_nm: f64,
    #[pyo3(get)]
    intensity: f64,
    #[pyo3(get)]
    from_mo: usize,
    #[pyo3(get)]
    to_mo: usize,
}

#[pyclass]
#[derive(Clone)]
struct UvVisSpectrumPy {
    #[pyo3(get)]
    energies_ev: Vec<f64>,
    #[pyo3(get)]
    intensities: Vec<f64>,
    #[pyo3(get)]
    peaks: Vec<UvVisPeakPy>,
    #[pyo3(get)]
    sigma: f64,
    #[pyo3(get)]
    notes: Vec<String>,
}

#[pyclass]
#[derive(Clone)]
struct GraphFeatureAnalysisPy {
    #[pyo3(get)]
    aromatic_atoms: Vec<bool>,
    #[pyo3(get)]
    aromatic_bonds: Vec<(usize, usize)>,
    #[pyo3(get)]
    tagged_tetrahedral_centers: Vec<usize>,
    #[pyo3(get)]
    inferred_tetrahedral_centers: Vec<usize>,
}

#[pymethods]
impl GraphFeatureAnalysisPy {
    fn __repr__(&self) -> String {
        format!(
            "GraphFeatureAnalysis(aromatic_atoms={}, tagged_centers={}, inferred_centers={})",
            self.aromatic_atoms.iter().filter(|v| **v).count(),
            self.tagged_tetrahedral_centers.len(),
            self.inferred_tetrahedral_centers.len()
        )
    }
}

#[pymethods]
impl UvVisSpectrumPy {
    fn __repr__(&self) -> String {
        format!(
            "UvVisSpectrum(n_points={}, n_peaks={})",
            self.energies_ev.len(),
            self.peaks.len()
        )
    }
}

#[pyclass]
#[derive(Clone)]
struct MetalCoordinationCenterPy {
    #[pyo3(get)]
    atom_index: usize,
    #[pyo3(get)]
    element: u8,
    #[pyo3(get)]
    ligand_indices: Vec<usize>,
    #[pyo3(get)]
    coordination_number: usize,
    #[pyo3(get)]
    geometry: String,
    #[pyo3(get)]
    geometry_score: f64,
}

#[pyclass]
#[derive(Clone)]
struct TopologyAnalysisResultPy {
    #[pyo3(get)]
    metal_centers: Vec<MetalCoordinationCenterPy>,
    #[pyo3(get)]
    warnings: Vec<String>,
}

#[pyclass]
#[derive(Clone)]
struct UffHeuristicEnergyPy {
    #[pyo3(get)]
    raw_energy_kcal_mol: f64,
    #[pyo3(get)]
    aromatic_stabilization_kcal_mol: f64,
    #[pyo3(get)]
    corrected_energy_kcal_mol: f64,
    #[pyo3(get)]
    aromatic_bond_count: usize,
    #[pyo3(get)]
    notes: Vec<String>,
}

#[pyclass]
#[derive(Clone)]
struct EmpiricalPkaSitePy {
    #[pyo3(get)]
    atom_index: usize,
    #[pyo3(get)]
    site_type: String,
    #[pyo3(get)]
    environment: String,
    #[pyo3(get)]
    estimated_pka: f64,
    #[pyo3(get)]
    confidence: f64,
}

#[pyclass]
#[derive(Clone)]
struct EmpiricalPkaResultPy {
    #[pyo3(get)]
    acidic_sites: Vec<EmpiricalPkaSitePy>,
    #[pyo3(get)]
    basic_sites: Vec<EmpiricalPkaSitePy>,
    #[pyo3(get)]
    notes: Vec<String>,
}

#[pymethods]
impl FrontierDescriptorsPy {
    fn __repr__(&self) -> String {
        format!(
            "FrontierDescriptors(n_atoms={}, gap={:.3} eV)",
            self.num_atoms, self.gap
        )
    }
}

/// Compute Mulliken & Löwdin population analysis.
#[pyfunction]
fn population(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<PopulationResultPy> {
    let positions: Vec<[f64; 3]> = coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    match sci_form_core::compute_population(&elements, &positions) {
        Ok(r) => Ok(PopulationResultPy {
            mulliken_charges: r.mulliken_charges,
            lowdin_charges: r.lowdin_charges,
            num_atoms: r.num_atoms,
            total_charge_mulliken: r.total_charge_mulliken,
            total_charge_lowdin: r.total_charge_lowdin,
        }),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

/// Compute Wiberg-like and Mayer-like bond orders from EHT.
#[pyfunction]
fn bond_orders(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<BondOrderResultPy> {
    let positions: Vec<[f64; 3]> = coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    match sci_form_core::compute_bond_orders(&elements, &positions) {
        Ok(r) => Ok(BondOrderResultPy {
            atom_pairs: r
                .bonds
                .iter()
                .map(|bond| (bond.atom_i, bond.atom_j))
                .collect(),
            distances: r.bonds.iter().map(|bond| bond.distance).collect(),
            wiberg: r.bonds.iter().map(|bond| bond.wiberg).collect(),
            mayer: r.bonds.iter().map(|bond| bond.mayer).collect(),
            wiberg_valence: r.wiberg_valence,
            mayer_valence: r.mayer_valence,
            num_atoms: r.num_atoms,
        }),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

/// Compute atom-resolved HOMO/LUMO frontier descriptors.
#[pyfunction]
fn frontier_descriptors(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<FrontierDescriptorsPy> {
    let positions: Vec<[f64; 3]> = coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    match sci_form_core::compute_frontier_descriptors(&elements, &positions) {
        Ok(r) => Ok(FrontierDescriptorsPy {
            homo_atom_contributions: r.homo_atom_contributions,
            lumo_atom_contributions: r.lumo_atom_contributions,
            dual_descriptor: r.dual_descriptor,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap: r.gap,
            num_atoms: r.num_atoms,
        }),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

/// Compute Fukui-function workflows and condensed per-atom descriptors.
#[pyfunction]
fn fukui_descriptors(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<FukuiDescriptorsPy> {
    let positions: Vec<[f64; 3]> = coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    match sci_form_core::compute_fukui_descriptors(&elements, &positions) {
        Ok(r) => Ok(FukuiDescriptorsPy {
            condensed_atom_indices: r.condensed.iter().map(|row| row.atom_index).collect(),
            condensed_f_plus: r.condensed.iter().map(|row| row.f_plus).collect(),
            condensed_f_minus: r.condensed.iter().map(|row| row.f_minus).collect(),
            condensed_f_radical: r.condensed.iter().map(|row| row.f_radical).collect(),
            condensed_dual_descriptor: r.condensed.iter().map(|row| row.dual_descriptor).collect(),
            f_plus: r.f_plus,
            f_minus: r.f_minus,
            f_radical: r.f_radical,
            dual_descriptor: r.dual_descriptor,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap: r.gap,
            num_atoms: r.num_atoms,
            validity_notes: r.validity_notes,
        }),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

/// Build empirical local-reactivity rankings from Fukui descriptors and Mulliken charges.
#[pyfunction]
fn reactivity_ranking(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<ReactivityRankingPy> {
    let positions: Vec<[f64; 3]> = coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    match sci_form_core::compute_reactivity_ranking(&elements, &positions) {
        Ok(r) => {
            let to_scores = |scores: Vec<sci_form_core::reactivity::ReactivitySiteScore>| {
                scores
                    .into_iter()
                    .map(|s| ReactivitySiteScorePy {
                        atom_index: s.atom_index,
                        score: s.score,
                    })
                    .collect::<Vec<_>>()
            };

            Ok(ReactivityRankingPy {
                nucleophilic_attack_sites: to_scores(r.nucleophilic_attack_sites),
                electrophilic_attack_sites: to_scores(r.electrophilic_attack_sites),
                radical_attack_sites: to_scores(r.radical_attack_sites),
                notes: r.notes,
            })
        }
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

/// Build an exploratory UV-Vis-like spectrum from low-cost EHT transitions.
#[pyfunction]
#[pyo3(signature = (elements, coords, sigma=0.25, e_min=0.5, e_max=8.0, n_points=600))]
fn uvvis_spectrum(
    elements: Vec<u8>,
    coords: Vec<f64>,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> PyResult<UvVisSpectrumPy> {
    let positions: Vec<[f64; 3]> = coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    match sci_form_core::compute_uv_vis_spectrum(
        &elements, &positions, sigma, e_min, e_max, n_points,
    ) {
        Ok(r) => Ok(UvVisSpectrumPy {
            energies_ev: r.energies_ev,
            intensities: r.intensities,
            peaks: r
                .peaks
                .into_iter()
                .map(|p| UvVisPeakPy {
                    energy_ev: p.energy_ev,
                    wavelength_nm: p.wavelength_nm,
                    intensity: p.intensity,
                    from_mo: p.from_mo,
                    to_mo: p.to_mo,
                })
                .collect(),
            sigma: r.sigma,
            notes: r.notes,
        }),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

/// Analyze aromaticity and graph-level stereocenters from a SMILES string.
#[pyfunction]
fn graph_features(smiles: &str) -> PyResult<GraphFeatureAnalysisPy> {
    match sci_form_core::analyze_graph_features(smiles) {
        Ok(r) => Ok(GraphFeatureAnalysisPy {
            aromatic_atoms: r.aromaticity.aromatic_atoms,
            aromatic_bonds: r.aromaticity.aromatic_bonds,
            tagged_tetrahedral_centers: r.stereocenters.tagged_tetrahedral_centers,
            inferred_tetrahedral_centers: r.stereocenters.inferred_tetrahedral_centers,
        }),
        Err(e) => Err(pyo3::exceptions::PyValueError::new_err(e)),
    }
}

/// Detect metal coordination geometry and return structured topology outputs.
#[pyfunction]
fn topology_analysis(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<TopologyAnalysisResultPy> {
    let positions: Vec<[f64; 3]> = coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    let result = sci_form_core::compute_topology(&elements, &positions);
    Ok(TopologyAnalysisResultPy {
        metal_centers: result
            .metal_centers
            .into_iter()
            .map(|center| MetalCoordinationCenterPy {
                atom_index: center.atom_index,
                element: center.element,
                ligand_indices: center.ligand_indices,
                coordination_number: center.coordination_number,
                geometry: format!("{:?}", center.geometry).to_lowercase(),
                geometry_score: center.geometry_score,
            })
            .collect(),
        warnings: result.warnings,
    })
}

// ─── Dipole Moment API ──────────────────────────────────────────────────────

#[pyclass]
#[derive(Clone)]
struct DipoleResultPy {
    #[pyo3(get)]
    vector: [f64; 3],
    #[pyo3(get)]
    magnitude: f64,
    #[pyo3(get)]
    unit: String,
}

#[pymethods]
impl DipoleResultPy {
    fn __repr__(&self) -> String {
        format!("DipoleResult({:.3} {})", self.magnitude, self.unit)
    }
}

/// Compute molecular dipole moment.
#[pyfunction]
fn dipole(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<DipoleResultPy> {
    let positions: Vec<[f64; 3]> = coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    match sci_form_core::compute_dipole(&elements, &positions) {
        Ok(r) => Ok(DipoleResultPy {
            vector: r.vector,
            magnitude: r.magnitude,
            unit: r.unit,
        }),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

// ─── DOS/PDOS API ───────────────────────────────────────────────────────────

#[pyclass]
#[derive(Clone)]
struct DosResultPy {
    #[pyo3(get)]
    energies: Vec<f64>,
    #[pyo3(get)]
    total_dos: Vec<f64>,
    #[pyo3(get)]
    sigma: f64,
}

#[pymethods]
impl DosResultPy {
    fn __repr__(&self) -> String {
        format!(
            "DosResult(n_points={}, sigma={:.2} eV)",
            self.energies.len(),
            self.sigma
        )
    }
}

/// Compute density of states (DOS/PDOS).
#[pyfunction]
#[pyo3(signature = (elements, coords, sigma=0.3, e_min=-30.0, e_max=5.0, n_points=500))]
fn dos(
    elements: Vec<u8>,
    coords: Vec<f64>,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> PyResult<DosResultPy> {
    let positions: Vec<[f64; 3]> = coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    match sci_form_core::compute_dos(&elements, &positions, sigma, e_min, e_max, n_points) {
        Ok(r) => Ok(DosResultPy {
            energies: r.energies,
            total_dos: r.total_dos,
            sigma: r.sigma,
        }),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

// ─── RMSD / Kabsch Alignment API ────────────────────────────────────────────

#[pyclass]
#[derive(Clone)]
struct AlignmentResultPy {
    #[pyo3(get)]
    rmsd: f64,
    #[pyo3(get)]
    aligned_coords: Vec<f64>,
}

#[pymethods]
impl AlignmentResultPy {
    fn __repr__(&self) -> String {
        format!("AlignmentResult(rmsd={:.4} Å)", self.rmsd)
    }
}

/// Compute RMSD between two coordinate sets (Kabsch alignment).
#[pyfunction]
fn rmsd(coords: Vec<f64>, reference: Vec<f64>) -> PyResult<AlignmentResultPy> {
    let result = sci_form_core::alignment::align_coordinates(&coords, &reference);
    Ok(AlignmentResultPy {
        rmsd: result.rmsd,
        aligned_coords: result.aligned_coords,
    })
}

// ─── UFF Force Field API ────────────────────────────────────────────────────

/// Compute UFF force field energy for a molecule.
#[pyfunction]
fn uff_energy(smiles: &str, coords: Vec<f64>) -> PyResult<f64> {
    match sci_form_core::compute_uff_energy(smiles, &coords) {
        Ok(e) => Ok(e),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

/// Compute UFF energy with aromaticity-informed heuristic correction metadata.
#[pyfunction]
fn uff_energy_with_aromatic_heuristics(
    smiles: &str,
    coords: Vec<f64>,
) -> PyResult<UffHeuristicEnergyPy> {
    match sci_form_core::compute_uff_energy_with_aromatic_heuristics(smiles, &coords) {
        Ok(result) => Ok(UffHeuristicEnergyPy {
            raw_energy_kcal_mol: result.raw_energy_kcal_mol,
            aromatic_stabilization_kcal_mol: result.aromatic_stabilization_kcal_mol,
            corrected_energy_kcal_mol: result.corrected_energy_kcal_mol,
            aromatic_bond_count: result.aromatic_bond_count,
            notes: result.notes,
        }),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

/// Estimate acidic/basic pKa sites from graph and charge heuristics.
#[pyfunction]
fn empirical_pka(smiles: &str) -> PyResult<EmpiricalPkaResultPy> {
    match sci_form_core::compute_empirical_pka(smiles) {
        Ok(result) => {
            let map_site = |s: sci_form_core::reactivity::EmpiricalPkaSite| EmpiricalPkaSitePy {
                atom_index: s.atom_index,
                site_type: s.site_type,
                environment: s.environment,
                estimated_pka: s.estimated_pka,
                confidence: s.confidence,
            };

            Ok(EmpiricalPkaResultPy {
                acidic_sites: result.acidic_sites.into_iter().map(map_site).collect(),
                basic_sites: result.basic_sites.into_iter().map(map_site).collect(),
                notes: result.notes,
            })
        }
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

// ─── Materials Assembly API ──────────────────────────────────────────────────

/// Unit cell result with parameters and volume.
#[pyclass]
#[derive(Clone)]
struct UnitCellPy {
    #[pyo3(get)]
    a: f64,
    #[pyo3(get)]
    b: f64,
    #[pyo3(get)]
    c: f64,
    #[pyo3(get)]
    alpha: f64,
    #[pyo3(get)]
    beta: f64,
    #[pyo3(get)]
    gamma: f64,
    #[pyo3(get)]
    volume: f64,
    #[pyo3(get)]
    lattice: Vec<Vec<f64>>,
}

/// Crystal structure result.
#[pyclass]
#[derive(Clone)]
struct CrystalStructurePy {
    #[pyo3(get)]
    num_atoms: usize,
    #[pyo3(get)]
    elements: Vec<u8>,
    #[pyo3(get)]
    frac_coords: Vec<Vec<f64>>,
    #[pyo3(get)]
    cart_coords: Vec<Vec<f64>>,
    #[pyo3(get)]
    labels: Vec<String>,
    #[pyo3(get)]
    lattice: Vec<Vec<f64>>,
}

/// Create a unit cell from crystallographic parameters.
#[pyfunction]
#[pyo3(signature = (a, b, c, alpha=90.0, beta=90.0, gamma=90.0))]
fn unit_cell(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> UnitCellPy {
    let cell = sci_form_core::create_unit_cell(a, b, c, alpha, beta, gamma);
    let p = cell.parameters();
    let vol = cell.volume();
    UnitCellPy {
        a: p.a,
        b: p.b,
        c: p.c,
        alpha: p.alpha,
        beta: p.beta,
        gamma: p.gamma,
        volume: vol,
        lattice: cell.lattice.iter().map(|r| r.to_vec()).collect(),
    }
}

/// Assemble a framework crystal structure.
#[pyfunction]
#[pyo3(signature = (topology="pcu", metal=30, geometry="octahedral", lattice_a=10.0, supercell=1))]
fn assemble_framework(
    topology: &str,
    metal: u8,
    geometry: &str,
    lattice_a: f64,
    supercell: usize,
) -> PyResult<CrystalStructurePy> {
    let geom = match geometry {
        "linear" => sci_form_core::materials::CoordinationGeometry::Linear,
        "trigonal" => sci_form_core::materials::CoordinationGeometry::Trigonal,
        "tetrahedral" => sci_form_core::materials::CoordinationGeometry::Tetrahedral,
        "square_planar" => sci_form_core::materials::CoordinationGeometry::SquarePlanar,
        "octahedral" => sci_form_core::materials::CoordinationGeometry::Octahedral,
        _ => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "Unknown geometry: {}",
                geometry
            )))
        }
    };
    let topo = match topology {
        "pcu" => sci_form_core::materials::Topology::pcu(),
        "dia" => sci_form_core::materials::Topology::dia(),
        "sql" => sci_form_core::materials::Topology::sql(),
        _ => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "Unknown topology: {}",
                topology
            )))
        }
    };
    let node = sci_form_core::materials::Sbu::metal_node(metal, 0.0, geom);
    let linker = sci_form_core::materials::Sbu::linear_linker(&[6, 6], 1.4, "carboxylate");
    let cell = sci_form_core::materials::UnitCell::cubic(lattice_a);
    let mut structure = sci_form_core::assemble_framework(&node, &linker, &topo, &cell);
    if supercell > 1 {
        structure = structure.make_supercell(supercell, supercell, supercell);
    }

    let cart = structure.cartesian_coords();
    Ok(CrystalStructurePy {
        num_atoms: structure.num_atoms(),
        elements: structure.elements(),
        frac_coords: structure
            .atoms
            .iter()
            .map(|a| a.frac_coords.to_vec())
            .collect(),
        cart_coords: cart.iter().map(|c| c.to_vec()).collect(),
        labels: structure.labels,
        lattice: structure.cell.lattice.iter().map(|r| r.to_vec()).collect(),
    })
}

// ─── C9: Transport / Streaming API ──────────────────────────────────────────

/// Arrow-compatible record batch result.
#[pyclass]
#[derive(Clone)]
struct RecordBatchPy {
    #[pyo3(get)]
    num_rows: usize,
    #[pyo3(get)]
    num_columns: usize,
    #[pyo3(get)]
    byte_size: usize,
    #[pyo3(get)]
    column_names: Vec<String>,
    #[pyo3(get)]
    float_data: std::collections::HashMap<String, Vec<f64>>,
    #[pyo3(get)]
    int_data: std::collections::HashMap<String, Vec<i32>>,
    #[pyo3(get)]
    uint8_data: std::collections::HashMap<String, Vec<u8>>,
}

/// Pack conformer results into Arrow-compatible columnar format.
#[pyfunction]
fn pack_conformers(results: Vec<ConformerResult>) -> RecordBatchPy {
    let core_results: Vec<sci_form_core::ConformerResult> = results
        .iter()
        .map(|r| sci_form_core::ConformerResult {
            smiles: r.smiles.clone(),
            num_atoms: r.num_atoms,
            coords: r.coords.clone(),
            elements: r.elements.clone(),
            bonds: vec![],
            error: r.error.clone(),
            time_ms: r.time_ms,
        })
        .collect();
    let batch = sci_form_core::transport::pack_conformers(&core_results);
    let mut float_data = std::collections::HashMap::new();
    let mut int_data = std::collections::HashMap::new();
    let mut uint8_data = std::collections::HashMap::new();
    let column_names: Vec<String> = batch.schema.iter().map(|s| s.name.clone()).collect();
    for c in &batch.float_columns {
        float_data.insert(c.name.clone(), c.values.clone());
    }
    for c in &batch.int_columns {
        int_data.insert(c.name.clone(), c.values.clone());
    }
    for c in &batch.uint8_columns {
        uint8_data.insert(c.name.clone(), c.values.clone());
    }
    RecordBatchPy {
        num_rows: batch.num_rows,
        num_columns: batch.num_columns(),
        byte_size: batch.byte_size(),
        column_names,
        float_data,
        int_data,
        uint8_data,
    }
}

/// Split SMILES batch into worker tasks.
#[pyfunction]
#[pyo3(signature = (smiles, n_workers=4, seed=42))]
fn split_worker_tasks(
    smiles: Vec<String>,
    n_workers: usize,
    seed: u64,
) -> Vec<std::collections::HashMap<String, String>> {
    let tasks = sci_form_core::transport::split_batch(&smiles, n_workers, seed);
    tasks
        .iter()
        .map(|t| {
            let mut m = std::collections::HashMap::new();
            m.insert("id".to_string(), t.id.to_string());
            // Format SMILES as JSON-like array string
            let smiles_str = format!(
                "[{}]",
                t.smiles
                    .iter()
                    .map(|s| format!("\"{}\"", s))
                    .collect::<Vec<_>>()
                    .join(",")
            );
            m.insert("smiles".to_string(), smiles_str);
            m.insert("kind".to_string(), format!("{:?}", t.kind));
            m
        })
        .collect()
}

/// Estimate optimal number of workers.
#[pyfunction]
#[pyo3(signature = (n_items, max_workers=8))]
fn estimate_workers(n_items: usize, max_workers: usize) -> usize {
    sci_form_core::transport::estimate_workers(n_items, max_workers)
}

// ─── MMFF94 Force Field ──────────────────────────────────────────────────────

/// Compute MMFF94 force field energy.
#[pyfunction]
fn mmff94_energy(smiles: &str, coords: Vec<f64>) -> PyResult<f64> {
    sci_form_core::compute_mmff94_energy(smiles, &coords)
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

// ─── PM3 Semi-Empirical Method ───────────────────────────────────────────────

/// PM3 calculation result.
#[pyclass]
#[derive(Clone)]
struct Pm3ResultPy {
    #[pyo3(get)]
    orbital_energies: Vec<f64>,
    #[pyo3(get)]
    electronic_energy: f64,
    #[pyo3(get)]
    nuclear_repulsion: f64,
    #[pyo3(get)]
    total_energy: f64,
    #[pyo3(get)]
    heat_of_formation: f64,
    #[pyo3(get)]
    n_basis: usize,
    #[pyo3(get)]
    n_electrons: usize,
    #[pyo3(get)]
    homo_energy: f64,
    #[pyo3(get)]
    lumo_energy: f64,
    #[pyo3(get)]
    gap: f64,
    #[pyo3(get)]
    mulliken_charges: Vec<f64>,
    #[pyo3(get)]
    scf_iterations: usize,
    #[pyo3(get)]
    converged: bool,
}

/// Run a PM3 semi-empirical calculation.
///
/// `elements`: list of atomic numbers.
/// `coords`: flat xyz list [x0,y0,z0, x1,y1,z1, ...] in Å.
#[pyfunction]
fn pm3_calculate(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<Pm3ResultPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    sci_form_core::compute_pm3(&elements, &positions)
        .map(|r| Pm3ResultPy {
            orbital_energies: r.orbital_energies,
            electronic_energy: r.electronic_energy,
            nuclear_repulsion: r.nuclear_repulsion,
            total_energy: r.total_energy,
            heat_of_formation: r.heat_of_formation,
            n_basis: r.n_basis,
            n_electrons: r.n_electrons,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap: r.gap,
            mulliken_charges: r.mulliken_charges,
            scf_iterations: r.scf_iterations,
            converged: r.converged,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

// ─── xTB Tight-Binding Method ────────────────────────────────────────────────

/// xTB calculation result.
#[pyclass]
#[derive(Clone)]
struct XtbResultPy {
    #[pyo3(get)]
    orbital_energies: Vec<f64>,
    #[pyo3(get)]
    total_energy: f64,
    #[pyo3(get)]
    repulsion_energy: f64,
    #[pyo3(get)]
    electronic_energy: f64,
    #[pyo3(get)]
    n_basis: usize,
    #[pyo3(get)]
    n_electrons: usize,
    #[pyo3(get)]
    homo_energy: f64,
    #[pyo3(get)]
    lumo_energy: f64,
    #[pyo3(get)]
    gap: f64,
    #[pyo3(get)]
    mulliken_charges: Vec<f64>,
    #[pyo3(get)]
    converged: bool,
    #[pyo3(get)]
    scf_iterations: usize,
}

/// Run an xTB tight-binding calculation.
///
/// `elements`: list of atomic numbers.
/// `coords`: flat xyz list [x0,y0,z0, x1,y1,z1, ...] in Å.
#[pyfunction]
fn xtb_calculate(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<XtbResultPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    sci_form_core::compute_xtb(&elements, &positions)
        .map(|r| XtbResultPy {
            orbital_energies: r.orbital_energies,
            total_energy: r.total_energy,
            repulsion_energy: r.repulsive_energy,
            electronic_energy: r.electronic_energy,
            n_basis: r.n_basis,
            n_electrons: r.n_electrons,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap: r.gap,
            mulliken_charges: r.mulliken_charges,
            converged: r.converged,
            scf_iterations: r.scc_iterations,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

// ─── ML Property Proxies ─────────────────────────────────────────────────────

/// Molecular descriptors.
#[pyclass]
#[derive(Clone)]
struct MolecularDescriptorsPy {
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

/// ML-predicted property result.
#[pyclass]
#[derive(Clone)]
struct MlPropertyResultPy {
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

/// Compute molecular descriptors from a SMILES string.
///
/// Returns descriptors usable for ML property prediction.
#[pyfunction]
fn ml_descriptors(smiles: &str) -> PyResult<MolecularDescriptorsPy> {
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

/// Predict ML properties from a SMILES string.
///
/// Returns LogP, molar refractivity, solubility (log S), Lipinski Ro5, and druglikeness.
#[pyfunction]
fn ml_predict(smiles: &str) -> PyResult<MlPropertyResultPy> {
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

// ─── Spectroscopy: UV-Vis sTDA ────────────────────────────────────────────────

/// sTDA UV-Vis excitation.
#[pyclass]
#[derive(Clone)]
struct StdaExcitationPy {
    #[pyo3(get)]
    energy_ev: f64,
    #[pyo3(get)]
    wavelength_nm: f64,
    #[pyo3(get)]
    oscillator_strength: f64,
    #[pyo3(get)]
    from_mo: usize,
    #[pyo3(get)]
    to_mo: usize,
    #[pyo3(get)]
    transition_dipole: f64,
}

/// sTDA UV-Vis spectrum result.
#[pyclass]
#[derive(Clone)]
struct StdaUvVisSpectrumPy {
    #[pyo3(get)]
    energies_ev: Vec<f64>,
    #[pyo3(get)]
    wavelengths_nm: Vec<f64>,
    #[pyo3(get)]
    absorptivity: Vec<f64>,
    #[pyo3(get)]
    excitations: Vec<StdaExcitationPy>,
    #[pyo3(get)]
    sigma: f64,
    #[pyo3(get)]
    broadening: String,
    #[pyo3(get)]
    notes: Vec<String>,
}

/// Compute sTDA UV-Vis spectrum with proper oscillator strengths.
///
/// `elements`: list of atomic numbers.
/// `coords`: flat xyz list in Å.
/// `sigma`: broadening width in eV.
/// `e_min`, `e_max`: energy range in eV.
/// `n_points`: grid resolution.
/// `broadening`: "gaussian" (default) or "lorentzian".
#[pyfunction]
#[pyo3(signature = (elements, coords, sigma=0.3, e_min=1.0, e_max=8.0, n_points=500, broadening="gaussian"))]
fn stda_uvvis(
    elements: Vec<u8>,
    coords: Vec<f64>,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
    broadening: &str,
) -> PyResult<StdaUvVisSpectrumPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    let bt = match broadening {
        "lorentzian" | "Lorentzian" => sci_form_core::reactivity::BroadeningType::Lorentzian,
        _ => sci_form_core::reactivity::BroadeningType::Gaussian,
    };
    sci_form_core::compute_stda_uvvis(&elements, &positions, sigma, e_min, e_max, n_points, bt)
        .map(|r| StdaUvVisSpectrumPy {
            energies_ev: r.energies_ev,
            wavelengths_nm: r.wavelengths_nm,
            absorptivity: r.absorptivity,
            excitations: r
                .excitations
                .iter()
                .map(|e| StdaExcitationPy {
                    energy_ev: e.energy_ev,
                    wavelength_nm: e.wavelength_nm,
                    oscillator_strength: e.oscillator_strength,
                    from_mo: e.from_mo,
                    to_mo: e.to_mo,
                    transition_dipole: e.transition_dipole,
                })
                .collect(),
            sigma: r.sigma,
            broadening: format!("{:?}", r.broadening),
            notes: r.notes,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

// ─── Spectroscopy: IR ─────────────────────────────────────────────────────────

/// Vibrational mode result.
#[pyclass]
#[derive(Clone)]
struct VibrationalModePy {
    #[pyo3(get)]
    frequency_cm1: f64,
    #[pyo3(get)]
    ir_intensity: f64,
    #[pyo3(get)]
    displacement: Vec<f64>,
    #[pyo3(get)]
    is_real: bool,
}

/// Vibrational analysis result.
#[pyclass]
#[derive(Clone)]
struct VibrationalAnalysisPy {
    #[pyo3(get)]
    n_atoms: usize,
    #[pyo3(get)]
    modes: Vec<VibrationalModePy>,
    #[pyo3(get)]
    n_real_modes: usize,
    #[pyo3(get)]
    zpve_ev: f64,
    #[pyo3(get)]
    method: String,
    #[pyo3(get)]
    notes: Vec<String>,
}

/// IR peak.
#[pyclass]
#[derive(Clone)]
struct IrPeakPy {
    #[pyo3(get)]
    frequency_cm1: f64,
    #[pyo3(get)]
    ir_intensity: f64,
    #[pyo3(get)]
    mode_index: usize,
}

/// IR spectrum result.
#[pyclass]
#[derive(Clone)]
struct IrSpectrumPy {
    #[pyo3(get)]
    wavenumbers: Vec<f64>,
    #[pyo3(get)]
    intensities: Vec<f64>,
    #[pyo3(get)]
    peaks: Vec<IrPeakPy>,
    #[pyo3(get)]
    gamma: f64,
    #[pyo3(get)]
    notes: Vec<String>,
}

/// Perform vibrational analysis via numerical Hessian.
///
/// `elements`: list of atomic numbers.
/// `coords`: flat xyz list in Å.
/// `method`: "eht", "pm3", or "xtb".
/// `step_size`: finite-difference step in Å (default 0.005).
#[pyfunction]
#[pyo3(signature = (elements, coords, method="eht", step_size=None))]
fn vibrational_analysis(
    elements: Vec<u8>,
    coords: Vec<f64>,
    method: &str,
    step_size: Option<f64>,
) -> PyResult<VibrationalAnalysisPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    sci_form_core::compute_vibrational_analysis(&elements, &positions, method, step_size)
        .map(|r| VibrationalAnalysisPy {
            n_atoms: r.n_atoms,
            modes: r
                .modes
                .iter()
                .map(|m| VibrationalModePy {
                    frequency_cm1: m.frequency_cm1,
                    ir_intensity: m.ir_intensity,
                    displacement: m.displacement.clone(),
                    is_real: m.is_real,
                })
                .collect(),
            n_real_modes: r.n_real_modes,
            zpve_ev: r.zpve_ev,
            method: r.method.clone(),
            notes: r.notes.clone(),
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

/// Generate Lorentzian-broadened IR spectrum.
///
/// `elements`: list of atomic numbers.
/// `coords`: flat xyz list in Å.
/// `method`: "eht", "pm3", or "xtb".
/// `gamma`: line width in cm⁻¹ (default 15).
/// `wn_min`, `wn_max`: wavenumber range in cm⁻¹.
/// `n_points`: grid resolution.
#[pyfunction]
#[pyo3(signature = (elements, coords, method="eht", gamma=15.0, wn_min=400.0, wn_max=4000.0, n_points=1000, step_size=None))]
fn ir_spectrum(
    elements: Vec<u8>,
    coords: Vec<f64>,
    method: &str,
    gamma: f64,
    wn_min: f64,
    wn_max: f64,
    n_points: usize,
    step_size: Option<f64>,
) -> PyResult<IrSpectrumPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    let analysis =
        sci_form_core::compute_vibrational_analysis(&elements, &positions, method, step_size)
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
    let spectrum = sci_form_core::compute_ir_spectrum(&analysis, gamma, wn_min, wn_max, n_points);
    Ok(IrSpectrumPy {
        wavenumbers: spectrum.wavenumbers,
        intensities: spectrum.intensities,
        peaks: spectrum
            .peaks
            .iter()
            .map(|p| IrPeakPy {
                frequency_cm1: p.frequency_cm1,
                ir_intensity: p.ir_intensity,
                mode_index: p.mode_index,
            })
            .collect(),
        gamma: spectrum.gamma,
        notes: spectrum.notes,
    })
}

// ─── Spectroscopy: NMR ────────────────────────────────────────────────────────

/// NMR chemical shift.
#[pyclass]
#[derive(Clone)]
struct ChemicalShiftPy {
    #[pyo3(get)]
    atom_index: usize,
    #[pyo3(get)]
    element: u8,
    #[pyo3(get)]
    shift_ppm: f64,
    #[pyo3(get)]
    environment: String,
    #[pyo3(get)]
    confidence: f64,
}

/// NMR chemical shift result (¹H + ¹³C).
#[pyclass]
#[derive(Clone)]
struct NmrShiftResultPy {
    #[pyo3(get)]
    h_shifts: Vec<ChemicalShiftPy>,
    #[pyo3(get)]
    c_shifts: Vec<ChemicalShiftPy>,
    #[pyo3(get)]
    notes: Vec<String>,
}

/// J-coupling constant.
#[pyclass]
#[derive(Clone)]
struct JCouplingPy {
    #[pyo3(get)]
    h1_index: usize,
    #[pyo3(get)]
    h2_index: usize,
    #[pyo3(get)]
    j_hz: f64,
    #[pyo3(get)]
    n_bonds: usize,
    #[pyo3(get)]
    coupling_type: String,
}

/// NMR peak.
#[pyclass]
#[derive(Clone)]
struct NmrPeakPy {
    #[pyo3(get)]
    shift_ppm: f64,
    #[pyo3(get)]
    intensity: f64,
    #[pyo3(get)]
    atom_index: usize,
    #[pyo3(get)]
    multiplicity: String,
    #[pyo3(get)]
    environment: String,
}

/// NMR spectrum result.
#[pyclass]
#[derive(Clone)]
struct NmrSpectrumPy {
    #[pyo3(get)]
    ppm_axis: Vec<f64>,
    #[pyo3(get)]
    intensities: Vec<f64>,
    #[pyo3(get)]
    peaks: Vec<NmrPeakPy>,
    #[pyo3(get)]
    nucleus: String,
    #[pyo3(get)]
    gamma: f64,
    #[pyo3(get)]
    notes: Vec<String>,
}

/// Predict ¹H and ¹³C NMR chemical shifts from SMILES.
#[pyfunction]
fn nmr_shifts(smiles: &str) -> PyResult<NmrShiftResultPy> {
    sci_form_core::predict_nmr_shifts(smiles)
        .map(|r| NmrShiftResultPy {
            h_shifts: r
                .h_shifts
                .iter()
                .map(|s| ChemicalShiftPy {
                    atom_index: s.atom_index,
                    element: s.element,
                    shift_ppm: s.shift_ppm,
                    environment: s.environment.clone(),
                    confidence: s.confidence,
                })
                .collect(),
            c_shifts: r
                .c_shifts
                .iter()
                .map(|s| ChemicalShiftPy {
                    atom_index: s.atom_index,
                    element: s.element,
                    shift_ppm: s.shift_ppm,
                    environment: s.environment.clone(),
                    confidence: s.confidence,
                })
                .collect(),
            notes: r.notes,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

/// Predict J-coupling constants.
///
/// `smiles`: SMILES string.
/// `coords`: flat xyz list in Å (or empty for topological estimate).
#[pyfunction]
#[pyo3(signature = (smiles, coords=vec![]))]
fn nmr_couplings(smiles: &str, coords: Vec<f64>) -> PyResult<Vec<JCouplingPy>> {
    let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    sci_form_core::predict_nmr_couplings(smiles, &positions)
        .map(|couplings| {
            couplings
                .iter()
                .map(|c| JCouplingPy {
                    h1_index: c.h1_index,
                    h2_index: c.h2_index,
                    j_hz: c.j_hz,
                    n_bonds: c.n_bonds,
                    coupling_type: c.coupling_type.clone(),
                })
                .collect()
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

/// Generate a complete NMR spectrum from SMILES.
///
/// `smiles`: SMILES string.
/// `nucleus`: "1H" or "13C".
/// `gamma`: Lorentzian line width in ppm.
/// `ppm_min`, `ppm_max`: spectral window.
/// `n_points`: grid resolution.
#[pyfunction]
#[pyo3(signature = (smiles, nucleus="1H", gamma=0.02, ppm_min=0.0, ppm_max=12.0, n_points=1000))]
fn nmr_spectrum(
    smiles: &str,
    nucleus: &str,
    gamma: f64,
    ppm_min: f64,
    ppm_max: f64,
    n_points: usize,
) -> PyResult<NmrSpectrumPy> {
    sci_form_core::compute_nmr_spectrum(smiles, nucleus, gamma, ppm_min, ppm_max, n_points)
        .map(|r| NmrSpectrumPy {
            ppm_axis: r.ppm_axis,
            intensities: r.intensities,
            peaks: r
                .peaks
                .iter()
                .map(|p| NmrPeakPy {
                    shift_ppm: p.shift_ppm,
                    intensity: p.intensity,
                    atom_index: p.atom_index,
                    multiplicity: p.multiplicity.clone(),
                    environment: p.environment.clone(),
                })
                .collect(),
            nucleus: format!("{:?}", r.nucleus),
            gamma: r.gamma,
            notes: r.notes,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

/// Generate HOSE codes for all atoms in a molecule.
///
/// `smiles`: SMILES string.
/// `max_radius`: maximum sphere radius (default 2).
#[pyfunction]
#[pyo3(signature = (smiles, max_radius=2))]
fn hose_codes(smiles: &str, max_radius: usize) -> PyResult<Vec<(usize, u8, String)>> {
    sci_form_core::compute_hose_codes(smiles, max_radius)
        .map(|codes| {
            codes
                .iter()
                .map(|c| (c.atom_index, c.element, c.full_code.clone()))
                .collect()
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}
