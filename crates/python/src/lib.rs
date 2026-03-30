use pyo3::prelude::*;
use pyo3::types::PyList;

mod alpha;
mod analysis;
mod beta;
mod compare;
mod dynamics;
mod eht;
mod electronic;
mod embed;
mod experimental;
mod forcefield;
mod hf_ani_esp;
mod materials;
mod mesh;
mod ml;
mod nmr;
mod population;
mod potentials;
mod properties;
mod reactivity;
mod rings;
mod solvation;
mod spectroscopy;
mod stereo;
mod system;
mod transport;

/// sci_form Python module
#[pymodule]
fn sci_form(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__package__", "sci_form")?;
    m.add("__path__", PyList::empty_bound(m.py()))?;
    if let Ok(spec) = m.getattr("__spec__") {
        if !spec.is_none() {
            spec.setattr("submodule_search_locations", PyList::empty_bound(m.py()))?;
        }
    }

    embed::register(m)?;
    system::register(m)?;
    compare::register(m)?;
    eht::register(m)?;
    properties::register(m)?;
    population::register(m)?;
    reactivity::register(m)?;
    analysis::register(m)?;
    forcefield::register(m)?;
    materials::register(m)?;
    transport::register(m)?;
    electronic::register(m)?;
    ml::register(m)?;
    spectroscopy::register(m)?;
    nmr::register(m)?;
    hf_ani_esp::register(m)?;
    mesh::register(m)?;
    dynamics::register(m)?;
    stereo::register(m)?;
    solvation::register(m)?;
    rings::register(m)?;
    experimental::register(m)?;
    alpha::register(m)?;
    beta::register(m)?;
    potentials::register(m)?;
    Ok(())
}
