//! Solvation bindings: non-polar SASA solvation, Generalized Born.

use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
pub(crate) struct NonPolarSolvationPy {
    #[pyo3(get)]
    energy_kcal_mol: f64,
    #[pyo3(get)]
    atom_contributions: Vec<f64>,
    #[pyo3(get)]
    atom_sasa: Vec<f64>,
    #[pyo3(get)]
    total_sasa: f64,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct GbSolvationPy {
    #[pyo3(get)]
    electrostatic_energy_kcal_mol: f64,
    #[pyo3(get)]
    nonpolar_energy_kcal_mol: f64,
    #[pyo3(get)]
    total_energy_kcal_mol: f64,
    #[pyo3(get)]
    born_radii: Vec<f64>,
    #[pyo3(get)]
    charges: Vec<f64>,
    #[pyo3(get)]
    solvent_dielectric: f64,
    #[pyo3(get)]
    solute_dielectric: f64,
}

/// Non-polar solvation energy from SASA with atomic solvation parameters.
#[pyfunction]
#[pyo3(signature = (elements, coords, probe_radius=1.4))]
fn nonpolar_solvation(
    elements: Vec<u8>,
    coords: Vec<f64>,
    probe_radius: f64,
) -> PyResult<NonPolarSolvationPy> {
    let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    let r = sci_form_core::compute_nonpolar_solvation(&elements, &positions, Some(probe_radius));
    Ok(NonPolarSolvationPy {
        energy_kcal_mol: r.energy_kcal_mol,
        atom_contributions: r.atom_contributions,
        atom_sasa: r.atom_sasa,
        total_sasa: r.total_sasa,
    })
}

/// Generalized Born solvation energy (electrostatic + non-polar).
#[pyfunction]
#[pyo3(signature = (elements, coords, charges, solvent_dielectric=78.5, solute_dielectric=1.0, probe_radius=1.4))]
fn gb_solvation(
    elements: Vec<u8>,
    coords: Vec<f64>,
    charges: Vec<f64>,
    solvent_dielectric: f64,
    solute_dielectric: f64,
    probe_radius: f64,
) -> PyResult<GbSolvationPy> {
    let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    let r = sci_form_core::compute_gb_solvation(
        &elements,
        &positions,
        &charges,
        Some(solvent_dielectric),
        Some(solute_dielectric),
        Some(probe_radius),
    );
    Ok(GbSolvationPy {
        electrostatic_energy_kcal_mol: r.electrostatic_energy_kcal_mol,
        nonpolar_energy_kcal_mol: r.nonpolar_energy_kcal_mol,
        total_energy_kcal_mol: r.total_energy_kcal_mol,
        born_radii: r.born_radii,
        charges: r.charges,
        solvent_dielectric: r.solvent_dielectric,
        solute_dielectric: r.solute_dielectric,
    })
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(nonpolar_solvation, m)?)?;
    m.add_function(wrap_pyfunction!(gb_solvation, m)?)?;
    m.add_class::<NonPolarSolvationPy>()?;
    m.add_class::<GbSolvationPy>()?;
    Ok(())
}
