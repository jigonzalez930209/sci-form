//! System capabilities and method plan bindings.

use pyo3::prelude::*;

pub(crate) fn coords_to_positions(coords: &[f64]) -> Vec<[f64; 3]> {
    coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect()
}

/// Perceive covalent bonds from 3D coordinates using covalent radii.
///
/// Used by the bindings that need a bond list (GETAWAY descriptors, hapticity
/// detection) or partial charges (RDF descriptors) but only receive a
/// geometry. When the caller supplies an explicit bond list it is used
/// verbatim instead.
pub(crate) fn perceive_bonds(
    elements: &[u8],
    positions: &[[f64; 3]],
    tolerance: f64,
) -> Vec<(usize, usize)> {
    let n = elements.len().min(positions.len());
    let mut bonds = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            let r_sum = sci_form_core::graph::get_covalent_radius(elements[i])
                + sci_form_core::graph::get_covalent_radius(elements[j]);
            if d > 0.4 && d < r_sum + tolerance {
                bonds.push((i, j));
            }
        }
    }
    bonds
}

/// Gasteiger-Marsili charges derived from geometry-perceived bonds.
pub(crate) fn gasteiger_charges_from_geometry(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<Vec<f64>, String> {
    let bonds = perceive_bonds(elements, positions, 0.3);
    let formal_charges = vec![0i8; elements.len()];
    let result = sci_form_core::charges::gasteiger::gasteiger_marsili_charges(
        elements,
        &bonds,
        &formal_charges,
        6,
    )?;
    Ok(result.charges)
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
