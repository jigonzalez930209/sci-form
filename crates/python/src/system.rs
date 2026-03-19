//! System capabilities and method plan bindings.

use pyo3::prelude::*;

pub(crate) fn coords_to_positions(coords: &[f64]) -> Vec<[f64; 3]> {
    coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect()
}

pub(crate) fn scientific_method_name(method: sci_form_core::ScientificMethod) -> String {
    match method {
        sci_form_core::ScientificMethod::Embed => "embed",
        sci_form_core::ScientificMethod::Uff => "uff",
        sci_form_core::ScientificMethod::Eht => "eht",
        sci_form_core::ScientificMethod::Pm3 => "pm3",
        sci_form_core::ScientificMethod::Xtb => "xtb",
        sci_form_core::ScientificMethod::Mmff94 => "mmff94",
        sci_form_core::ScientificMethod::Ani => "ani",
        sci_form_core::ScientificMethod::Hf3c => "hf3c",
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

#[pyclass]
#[derive(Clone)]
pub(crate) struct MethodCapabilityPy {
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
    fn from(v: sci_form_core::MethodCapability) -> Self {
        Self {
            available: v.available,
            confidence: format!("{:?}", v.confidence).to_lowercase(),
            unsupported_elements: v.unsupported_elements,
            warnings: v.warnings,
        }
    }
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct SystemCapabilitiesPy {
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
pub(crate) struct MethodMetadataPy {
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
    fn from(v: sci_form_core::MethodMetadata) -> Self {
        Self {
            method: scientific_method_name(v.method),
            available: v.available,
            confidence: format!("{:?}", v.confidence).to_lowercase(),
            confidence_score: v.confidence_score,
            limitations: v.limitations,
            warnings: v.warnings,
        }
    }
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct PropertyMethodPlanPy {
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
    fn from(v: sci_form_core::PropertyMethodPlan) -> Self {
        Self {
            property: property_request_name(v.property),
            recommended: v.recommended.map(scientific_method_name),
            fallback: v.fallback.map(scientific_method_name),
            rationale: v.rationale,
            methods: v.methods.into_iter().map(Into::into).collect(),
        }
    }
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct SystemMethodPlanPy {
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

pub(crate) fn caps_to_py(caps: sci_form_core::SystemCapabilities) -> SystemCapabilitiesPy {
    SystemCapabilitiesPy {
        embed: caps.embed.into(),
        uff: caps.uff.into(),
        eht: caps.eht.into(),
        population: caps.population.into(),
        orbital_grid: caps.orbital_grid.into(),
    }
}

pub(crate) fn plan_to_py(plan: sci_form_core::SystemMethodPlan) -> SystemMethodPlanPy {
    SystemMethodPlanPy {
        capabilities: caps_to_py(plan.capabilities),
        geometry: plan.geometry.into(),
        force_field_energy: plan.force_field_energy.into(),
        orbitals: plan.orbitals.into(),
        population: plan.population.into(),
        orbital_grid: plan.orbital_grid.into(),
    }
}

#[pyfunction]
fn system_capabilities(elements: Vec<u8>) -> SystemCapabilitiesPy {
    caps_to_py(sci_form_core::get_system_capabilities(&elements))
}

#[pyfunction]
fn system_method_plan(elements: Vec<u8>) -> SystemMethodPlanPy {
    plan_to_py(sci_form_core::get_system_method_plan(&elements))
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(system_capabilities, m)?)?;
    m.add_function(wrap_pyfunction!(system_method_plan, m)?)?;
    m.add_class::<SystemCapabilitiesPy>()?;
    m.add_class::<MethodMetadataPy>()?;
    m.add_class::<PropertyMethodPlanPy>()?;
    m.add_class::<SystemMethodPlanPy>()?;
    Ok(())
}
