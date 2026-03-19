use pyo3::prelude::*;

mod analysis;
mod compare;
mod dynamics;
mod eht;
mod electronic;
mod embed;
mod forcefield;
mod hf_ani_esp;
mod materials;
mod mesh;
mod ml;
mod nmr;
mod population;
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
    Ok(())
}
