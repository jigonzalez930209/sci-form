//! PM3 and xTB-family semi-empirical electronic structure bindings.

use crate::system::coords_to_positions;
use pyo3::prelude::*;

#[pyclass(name = "Pm3Result")]
#[derive(Clone)]
pub(crate) struct Pm3ResultPy {
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

#[pyclass(name = "XtbResult")]
#[derive(Clone)]
pub(crate) struct XtbResultPy {
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

#[pyclass(name = "Gfn1Result")]
#[derive(Clone)]
pub(crate) struct Gfn1ResultPy {
    #[pyo3(get)]
    orbital_energies: Vec<f64>,
    #[pyo3(get)]
    total_energy: f64,
    #[pyo3(get)]
    repulsive_energy: f64,
    #[pyo3(get)]
    dispersion_energy: f64,
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
    shell_charges: Vec<Vec<f64>>,
    #[pyo3(get)]
    converged: bool,
    #[pyo3(get)]
    scc_iterations: usize,
}

#[pyclass(name = "Gfn2Result")]
#[derive(Clone)]
pub(crate) struct Gfn2ResultPy {
    #[pyo3(get)]
    orbital_energies: Vec<f64>,
    #[pyo3(get)]
    total_energy: f64,
    #[pyo3(get)]
    repulsive_energy: f64,
    #[pyo3(get)]
    dispersion_energy: f64,
    #[pyo3(get)]
    halogen_bond_energy: f64,
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
    atomic_dipoles: Vec<(f64, f64, f64)>,
    #[pyo3(get)]
    atomic_quadrupoles: Vec<f64>,
    #[pyo3(get)]
    converged: bool,
    #[pyo3(get)]
    scc_iterations: usize,
}

#[pyclass(name = "CisdExcitation")]
#[derive(Clone)]
pub(crate) struct CisdExcitationPy {
    #[pyo3(get)]
    energy: f64,
    #[pyo3(get)]
    energy_ev: f64,
    #[pyo3(get)]
    wavelength_nm: f64,
    #[pyo3(get)]
    oscillator_strength: f64,
    #[pyo3(get)]
    dominant_transition: Vec<[usize; 2]>,
    #[pyo3(get)]
    character: String,
    #[pyo3(get)]
    singles_weight: f64,
}

#[pyclass(name = "CisdResult")]
#[derive(Clone)]
pub(crate) struct CisdResultPy {
    #[pyo3(get)]
    excitations: Vec<CisdExcitationPy>,
    #[pyo3(get)]
    correlation_energy: f64,
    #[pyo3(get)]
    n_csfs: usize,
    #[pyo3(get)]
    approximation: String,
    #[pyo3(get)]
    scf_energy: f64,
    #[pyo3(get)]
    scf_converged: bool,
    #[pyo3(get)]
    n_basis: usize,
    #[pyo3(get)]
    n_occupied: usize,
}
#[pymethods]
impl CisdResultPy {
    fn __repr__(&self) -> String {
        format!(
            "CisdResult(n_exc={}, correlation={:.6} Ha, approx={})",
            self.excitations.len(),
            self.correlation_energy,
            self.approximation
        )
    }
}

fn validate_coords(elements: &[u8], coords: &[f64]) -> PyResult<Vec<[f64; 3]>> {
    if coords.len() != elements.len() * 3 {
        return Err(crate::errors::value_error(format!(
            "coords length {} != elements.len() * 3 = {}",
            coords.len(),
            elements.len() * 3
        )));
    }
    Ok(coords_to_positions(coords))
}

#[pyfunction]
fn compute_pm3(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<Pm3ResultPy> {
    let positions = validate_coords(&elements, &coords)?;
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
        .map_err(crate::errors::engine_error)
}

#[pyfunction]
fn compute_xtb(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<XtbResultPy> {
    let positions = validate_coords(&elements, &coords)?;
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
        .map_err(crate::errors::engine_error)
}

#[pyfunction]
fn compute_gfn1(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<Gfn1ResultPy> {
    let positions = validate_coords(&elements, &coords)?;
    sci_form_core::xtb::solve_gfn1(&elements, &positions)
        .map(|r| Gfn1ResultPy {
            orbital_energies: r.orbital_energies,
            total_energy: r.total_energy,
            repulsive_energy: r.repulsive_energy,
            dispersion_energy: r.dispersion_energy,
            electronic_energy: r.electronic_energy,
            n_basis: r.n_basis,
            n_electrons: r.n_electrons,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap: r.gap,
            mulliken_charges: r.mulliken_charges,
            shell_charges: r.shell_charges,
            converged: r.converged,
            scc_iterations: r.scc_iterations,
        })
        .map_err(crate::errors::engine_error)
}

#[pyfunction]
fn compute_gfn2(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<Gfn2ResultPy> {
    let positions = validate_coords(&elements, &coords)?;
    sci_form_core::xtb::solve_gfn2(&elements, &positions)
        .map(|r| Gfn2ResultPy {
            orbital_energies: r.orbital_energies,
            total_energy: r.total_energy,
            repulsive_energy: r.repulsive_energy,
            dispersion_energy: r.dispersion_energy,
            halogen_bond_energy: r.halogen_bond_energy,
            electronic_energy: r.electronic_energy,
            n_basis: r.n_basis,
            n_electrons: r.n_electrons,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap: r.gap,
            mulliken_charges: r.mulliken_charges,
            atomic_dipoles: r
                .atomic_dipoles
                .into_iter()
                .map(|dipole| (dipole[0], dipole[1], dipole[2]))
                .collect(),
            atomic_quadrupoles: r.atomic_quadrupoles,
            converged: r.converged,
            scc_iterations: r.scc_iterations,
        })
        .map_err(crate::errors::engine_error)
}

/// Compute CISD excited states from an internal STO-3G RHF calculation.
///
/// The binding runs the shared SCF engine (RHF/STO-3G) and the validated
/// two-electron integral engine internally, then feeds the canonical orbitals
/// and ERIs to the CISD driver. `n_states` is the number of excited states
/// to report. Large systems fall back to the perturbative (CIS + MP2)
/// approximation automatically, reported in `approximation`.
#[pyfunction]
#[pyo3(signature = (elements, coords, n_states=5))]
fn compute_cisd(elements: Vec<u8>, coords: Vec<f64>, n_states: usize) -> PyResult<CisdResultPy> {
    let positions = validate_coords(&elements, &coords)?;

    // STO-3G RHF SCF (shared engine).
    let system =
        sci_form_core::scf::types::MolecularSystem::from_angstrom(&elements, &positions, 0, 1);
    let scf_config = sci_form_core::scf::scf_loop::ScfConfig::default();
    let scf = sci_form_core::scf::scf_loop::run_scf(&system, &scf_config);
    if !scf.converged {
        return Err(crate::errors::engine_error(
            "RHF/STO-3G SCF did not converge; CISD requires a converged reference",
        ));
    }

    // AO ERIs in the packed layout the CISD driver expects.
    let basis =
        sci_form_core::scf::basis::BasisSet::sto3g(&system.atomic_numbers, &system.positions_bohr);
    let n_basis = basis.n_basis;
    let engine = sci_form_core::scf::two_electron::TwoElectronIntegrals::compute(&basis);
    let n_pair = n_basis * (n_basis + 1) / 2;
    let mut eris = vec![0.0f64; n_pair * (n_pair + 1) / 2];
    for i in 0..n_basis {
        for j in 0..n_basis {
            for k in 0..n_basis {
                for l in 0..n_basis {
                    eris[sci_form_core::hf::integrals::eri_index(i, j, k, l, n_basis)] =
                        engine.get(i, j, k, l);
                }
            }
        }
    }

    let n_occupied = scf.n_electrons / 2;
    let result = sci_form_core::hf::compute_cisd(
        &scf.orbital_energies,
        &scf.mo_coefficients,
        &eris,
        n_basis,
        n_occupied,
        n_states,
    );
    Ok(CisdResultPy {
        excitations: result
            .excitations
            .iter()
            .map(|e| CisdExcitationPy {
                energy: e.energy,
                energy_ev: e.energy_ev,
                wavelength_nm: e.wavelength_nm,
                oscillator_strength: e.oscillator_strength,
                dominant_transition: e.dominant_transition.iter().map(|t| [t.0, t.1]).collect(),
                character: e.character.clone(),
                singles_weight: e.singles_weight,
            })
            .collect(),
        correlation_energy: result.correlation_energy,
        n_csfs: result.n_csfs,
        approximation: result.approximation,
        scf_energy: scf.total_energy,
        scf_converged: scf.converged,
        n_basis,
        n_occupied,
    })
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_pm3, m)?)?;
    m.add_function(wrap_pyfunction!(compute_xtb, m)?)?;
    m.add_function(wrap_pyfunction!(compute_gfn1, m)?)?;
    m.add_function(wrap_pyfunction!(compute_gfn2, m)?)?;
    m.add_function(wrap_pyfunction!(compute_cisd, m)?)?;
    m.add_class::<CisdExcitationPy>()?;
    m.add_class::<CisdResultPy>()?;
    m.add_class::<Pm3ResultPy>()?;
    m.add_class::<XtbResultPy>()?;
    m.add_class::<Gfn1ResultPy>()?;
    m.add_class::<Gfn2ResultPy>()?;
    Ok(())
}
