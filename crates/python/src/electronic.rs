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

/// Result of an open-shell SCF calculation (UHF or ROHF).
#[pyclass(name = "UhfResult")]
#[derive(Clone)]
pub(crate) struct UhfResultPy {
    /// Alpha orbital energies in Hartree, ascending.
    #[pyo3(get)]
    alpha_orbital_energies: Vec<f64>,
    /// Beta orbital energies in Hartree, ascending.
    #[pyo3(get)]
    beta_orbital_energies: Vec<f64>,
    /// Electronic energy in Hartree.
    #[pyo3(get)]
    electronic_energy: f64,
    /// Nuclear repulsion energy in Hartree.
    #[pyo3(get)]
    nuclear_repulsion: f64,
    /// Total energy in Hartree.
    #[pyo3(get)]
    total_energy: f64,
    /// HOMO energy in Hartree.
    #[pyo3(get)]
    homo_energy: f64,
    /// LUMO energy in Hartree, if the basis has an empty virtual space.
    #[pyo3(get)]
    lumo_energy: Option<f64>,
    /// HOMO-LUMO gap in eV.
    #[pyo3(get)]
    gap_ev: f64,
    /// Mulliken charges per atom, in units of `e`.
    #[pyo3(get)]
    mulliken_charges: Vec<f64>,
    /// Number of SCF iterations.
    #[pyo3(get)]
    scf_iterations: usize,
    /// Whether the SCF converged.
    #[pyo3(get)]
    converged: bool,
    /// Number of basis functions.
    #[pyo3(get)]
    n_basis: usize,
    /// Number of alpha electrons.
    #[pyo3(get)]
    n_alpha: usize,
    /// Number of beta electrons.
    #[pyo3(get)]
    n_beta: usize,
    /// Expectation value <S²>, in units of ħ².
    #[pyo3(get)]
    s2_expectation: f64,
    /// Spin contamination <S²> - S(S+1), in units of ħ².
    #[pyo3(get)]
    spin_contamination: f64,
    /// Alpha MO coefficients, one row per basis function.
    #[pyo3(get)]
    alpha_coefficients: Vec<Vec<f64>>,
    /// Beta MO coefficients, one row per basis function.
    #[pyo3(get)]
    beta_coefficients: Vec<Vec<f64>>,
}

impl UhfResultPy {
    fn from_core(r: sci_form_core::scf::uhf::UhfResult) -> Self {
        Self {
            alpha_orbital_energies: r.alpha_orbital_energies,
            beta_orbital_energies: r.beta_orbital_energies,
            electronic_energy: r.electronic_energy,
            nuclear_repulsion: r.nuclear_repulsion,
            total_energy: r.total_energy,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap_ev: r.gap_ev,
            mulliken_charges: r.mulliken_charges,
            scf_iterations: r.scf_iterations,
            converged: r.converged,
            n_basis: r.n_basis,
            n_alpha: r.n_alpha,
            n_beta: r.n_beta,
            s2_expectation: r.s2_expectation,
            spin_contamination: r.spin_contamination,
            alpha_coefficients: matrix_rows(&r.alpha_coefficients),
            beta_coefficients: matrix_rows(&r.beta_coefficients),
        }
    }
}

/// A dense matrix as a list of rows, so that the bindings do not have to
/// expose `nalgebra` types.
fn matrix_rows(m: &nalgebra::DMatrix<f64>) -> Vec<Vec<f64>> {
    m.row_iter()
        .map(|row| row.iter().copied().collect())
        .collect()
}

/// Result of a PM3 gradient calculation.
#[pyclass(name = "Pm3GradientResult")]
#[derive(Clone)]
pub(crate) struct Pm3GradientResultPy {
    /// Energy gradient dE/dR per atom in eV/Å, one `[x, y, z]` per atom.
    #[pyo3(get)]
    gradients: Vec<[f64; 3]>,
    /// Total PM3 energy in eV.
    #[pyo3(get)]
    energy: f64,
    /// Heat of formation in kcal/mol.
    #[pyo3(get)]
    heat_of_formation: f64,
}

/// Result of an xTB gradient calculation.
#[pyclass(name = "XtbGradientResult")]
#[derive(Clone)]
pub(crate) struct XtbGradientResultPy {
    /// Energy gradient dE/dR per atom in eV/Å, one `[x, y, z]` per atom.
    #[pyo3(get)]
    gradients: Vec<[f64; 3]>,
    /// Total xTB energy in eV.
    #[pyo3(get)]
    energy: f64,
}

/// Run an Unrestricted Hartree-Fock (UHF) calculation for open-shell systems.
///
/// `charge` is the total molecular charge (0 = neutral) and `multiplicity`
/// the spin multiplicity (1 = singlet, 2 = doublet, 3 = triplet).
#[pyfunction]
#[pyo3(signature = (elements, coords, charge=0, multiplicity=1))]
fn compute_uhf(
    elements: Vec<u8>,
    coords: Vec<f64>,
    charge: i32,
    multiplicity: u32,
) -> PyResult<UhfResultPy> {
    let positions = validate_coords(&elements, &coords)?;
    sci_form_core::compute_uhf(&elements, &positions, charge, multiplicity)
        .map(UhfResultPy::from_core)
        .map_err(crate::errors::engine_error)
}

/// Run a Restricted Open-shell Hartree-Fock (ROHF) calculation.
///
/// Same arguments as [`compute_uhf`]. ROHF uses the canonical projection
/// (the same spatial orbitals for both spins), so the spin contamination is
/// zero by construction.
#[pyfunction]
#[pyo3(signature = (elements, coords, charge=0, multiplicity=1))]
fn compute_rohf(
    elements: Vec<u8>,
    coords: Vec<f64>,
    charge: i32,
    multiplicity: u32,
) -> PyResult<UhfResultPy> {
    let positions = validate_coords(&elements, &coords)?;
    sci_form_core::compute_rohf(&elements, &positions, charge, multiplicity)
        .map(UhfResultPy::from_core)
        .map_err(crate::errors::engine_error)
}

/// Compute the PM3 analytical energy gradient (eV/Å per atom).
#[pyfunction]
fn compute_pm3_gradient(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<Pm3GradientResultPy> {
    let positions = validate_coords(&elements, &coords)?;
    sci_form_core::compute_pm3_gradient(&elements, &positions)
        .map(|r| Pm3GradientResultPy {
            gradients: r.gradients,
            energy: r.energy,
            heat_of_formation: r.heat_of_formation,
        })
        .map_err(crate::errors::engine_error)
}

/// Compute the xTB analytical energy gradient (eV/Å per atom).
#[pyfunction]
fn compute_xtb_gradient(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<XtbGradientResultPy> {
    let positions = validate_coords(&elements, &coords)?;
    sci_form_core::compute_xtb_gradient(&elements, &positions)
        .map(|r| XtbGradientResultPy {
            gradients: r.gradients,
            energy: r.energy,
        })
        .map_err(crate::errors::engine_error)
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
    m.add_function(wrap_pyfunction!(compute_uhf, m)?)?;
    m.add_function(wrap_pyfunction!(compute_rohf, m)?)?;
    m.add_function(wrap_pyfunction!(compute_pm3_gradient, m)?)?;
    m.add_function(wrap_pyfunction!(compute_xtb_gradient, m)?)?;
    m.add_class::<CisdExcitationPy>()?;
    m.add_class::<CisdResultPy>()?;
    m.add_class::<Pm3ResultPy>()?;
    m.add_class::<XtbResultPy>()?;
    m.add_class::<Gfn1ResultPy>()?;
    m.add_class::<Gfn2ResultPy>()?;
    m.add_class::<UhfResultPy>()?;
    m.add_class::<Pm3GradientResultPy>()?;
    m.add_class::<XtbGradientResultPy>()?;
    Ok(())
}
