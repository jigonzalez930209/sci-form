//! Property bindings: charges, SASA, dipole, DOS, RMSD.

use crate::errors::{engine_error, invalid_value, value_error};
use crate::system::coords_to_positions;
use pyo3::prelude::*;

#[pyclass(name = "ChargeResult")]
#[derive(Clone)]
pub(crate) struct ChargeResultPy {
    #[pyo3(get)]
    charges: Vec<f64>,
    #[pyo3(get)]
    iterations: usize,
    #[pyo3(get)]
    total_charge: f64,
    #[pyo3(get)]
    converged: bool,
}
#[pymethods]
impl ChargeResultPy {
    fn __repr__(&self) -> String {
        format!(
            "ChargeResult(n_atoms={}, total_charge={:.4}, iterations={}, converged={})",
            self.charges.len(),
            self.total_charge,
            self.iterations,
            self.converged
        )
    }
}

#[pyclass(name = "SasaResult")]
#[derive(Clone)]
pub(crate) struct SasaResultPy {
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
            self.probe_radius
        )
    }
}

#[pyclass(name = "DipoleResult")]
#[derive(Clone)]
pub(crate) struct DipoleResultPy {
    #[pyo3(get)]
    vector: [f64; 3],
    #[pyo3(get)]
    magnitude: f64,
    #[pyo3(get)]
    unit: String,
    #[pyo3(get)]
    nuclear_dipole: Option<[f64; 3]>,
    #[pyo3(get)]
    electronic_dipole: Option<[f64; 3]>,
}
#[pymethods]
impl DipoleResultPy {
    fn __repr__(&self) -> String {
        format!(
            "DipoleResult({:.3} {}, decomposed={})",
            self.magnitude,
            self.unit,
            self.nuclear_dipole.is_some()
        )
    }
}

#[pyclass(name = "DosResult")]
#[derive(Clone)]
pub(crate) struct DosResultPy {
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

#[pyclass(name = "AlignmentResult")]
#[derive(Clone)]
pub(crate) struct AlignmentResultPy {
    /// RMSD after the optimal Kabsch alignment, in Å.
    #[pyo3(get)]
    rmsd: f64,
    /// The aligned coordinates.
    #[pyo3(get)]
    aligned_coords: Vec<f64>,
    /// The optimal rotation matrix (row-major).
    #[pyo3(get)]
    rotation: [[f64; 3]; 3],
    /// The translation applied to the input centroid.
    #[pyo3(get)]
    translation: [f64; 3],
}
#[pymethods]
impl AlignmentResultPy {
    fn __repr__(&self) -> String {
        format!("AlignmentResult(rmsd={:.4} Å)", self.rmsd)
    }
}

#[pyclass(name = "MultiMethodDosResult")]
#[derive(Clone)]
pub(crate) struct MultiMethodDosResultPy {
    #[pyo3(get)]
    energies: Vec<f64>,
    #[pyo3(get)]
    total_dos: Vec<f64>,
    #[pyo3(get)]
    method: String,
    #[pyo3(get)]
    homo_energy: f64,
    #[pyo3(get)]
    lumo_energy: f64,
    #[pyo3(get)]
    gap: f64,
    #[pyo3(get)]
    orbital_energies: Vec<f64>,
    #[pyo3(get)]
    sigma: f64,
}
#[pymethods]
impl MultiMethodDosResultPy {
    fn __repr__(&self) -> String {
        format!(
            "MultiMethodDosResult(method={}, n_points={}, gap={:.3} eV)",
            self.method,
            self.energies.len(),
            self.gap
        )
    }
}

/// Compute Gasteiger-Marsili partial charges from a SMILES string.
///
/// Raises `ValueError` for a malformed SMILES and `RuntimeError` when the
/// charge solver fails.
#[pyfunction]
fn compute_charges(smiles: &str) -> PyResult<ChargeResultPy> {
    sci_form_core::compute_charges(smiles)
        .map(|r| ChargeResultPy {
            charges: r.charges,
            iterations: r.iterations,
            total_charge: r.total_charge,
            converged: r.converged,
        })
        .map_err(|e| {
            if e.contains("SMILES") || e.contains("parse") {
                value_error(crate::errors::smiles_error(&e, smiles))
            } else {
                engine_error(e)
            }
        })
}

/// Compute Gasteiger-Marsili charges with configurable iteration parameters.
///
/// `max_iter` is the number of electronegativity-equalization sweeps,
/// `initial_damping` the first damping factor (halved every iteration) and
/// `convergence_threshold` the max charge-transfer convergence criterion.
#[pyfunction]
#[pyo3(signature = (smiles, max_iter=6, initial_damping=0.5, convergence_threshold=1e-10))]
fn compute_charges_configured(
    smiles: &str,
    max_iter: usize,
    initial_damping: f64,
    convergence_threshold: f64,
) -> PyResult<ChargeResultPy> {
    let config = sci_form_core::charges::gasteiger::GasteigerConfig {
        max_iter,
        initial_damping,
        convergence_threshold,
    };
    sci_form_core::compute_charges_configured(smiles, &config)
        .map(|r| ChargeResultPy {
            charges: r.charges,
            iterations: r.iterations,
            total_charge: r.total_charge,
            converged: r.converged,
        })
        .map_err(|e| {
            if e.contains("SMILES") || e.contains("parse") {
                value_error(crate::errors::smiles_error(&e, smiles))
            } else {
                engine_error(e)
            }
        })
}

/// Compute the solvent-accessible surface area.
#[pyfunction]
#[pyo3(signature = (elements, coords, probe_radius=1.4))]
fn compute_sasa(
    elements: Vec<u8>,
    coords: Vec<f64>,
    probe_radius: f64,
) -> PyResult<SasaResultPy> {
    crate::errors::check_coords_len(&coords, elements.len())?;
    let pr = if (probe_radius - 1.4).abs() < 1e-10 {
        None
    } else {
        Some(probe_radius)
    };
    sci_form_core::compute_sasa(&elements, &coords, pr)
        .map(|r| SasaResultPy {
            total_sasa: r.total_sasa,
            atom_sasa: r.atom_sasa,
            probe_radius: r.probe_radius,
            num_points: r.num_points,
        })
        .map_err(value_error)
}

/// Compute the molecular dipole moment.
#[pyfunction]
fn compute_dipole(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<DipoleResultPy> {
    let positions = coords_to_positions(&coords);
    sci_form_core::compute_dipole(&elements, &positions)
        .map(|r| DipoleResultPy {
            vector: r.vector,
            magnitude: r.magnitude,
            unit: r.unit,
            nuclear_dipole: r.nuclear_dipole,
            electronic_dipole: r.electronic_dipole,
        })
        .map_err(engine_error)
}

/// Compute the density of states from the EHT orbital energies.
#[pyfunction]
#[pyo3(signature = (elements, coords, sigma=0.3, e_min=-30.0, e_max=5.0, n_points=500))]
fn compute_dos(
    elements: Vec<u8>,
    coords: Vec<f64>,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> PyResult<DosResultPy> {
    let positions = coords_to_positions(&coords);
    sci_form_core::compute_dos(&elements, &positions, sigma, e_min, e_max, n_points)
        .map(|r| DosResultPy {
            energies: r.energies,
            total_dos: r.total_dos,
            sigma: r.sigma,
        })
        .map_err(engine_error)
}

/// Compute the RMSD between two coordinate sets after a Kabsch alignment.
///
/// Raises `ValueError` when the two structures cannot be compared (different
/// lengths, or a length that is not a multiple of 3).  A `NaN` RMSD is
/// never returned (blocker C-C4).
#[pyfunction]
fn compute_rmsd(coords: Vec<f64>, reference: Vec<f64>) -> PyResult<AlignmentResultPy> {
    let result = sci_form_core::alignment::try_align_coordinates(&coords, &reference)
        .map_err(value_error)?;
    Ok(AlignmentResultPy {
        rmsd: result.rmsd,
        aligned_coords: result.aligned_coords,
        rotation: result.rotation,
        translation: result.translation,
    })
}

/// Compute DOS from the orbital energies of a chosen electronic-structure method.
///
/// `method` is one of `eht`, `pm3`, `xtb` (GFN2), `gfn1`, `gfn2` or `hf3c`.
#[pyfunction]
#[pyo3(signature = (elements, coords, method="eht", sigma=0.3, e_min=-30.0, e_max=5.0, n_points=500))]
fn compute_dos_multimethod(
    elements: Vec<u8>,
    coords: Vec<f64>,
    method: &str,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> PyResult<MultiMethodDosResultPy> {
    let positions = coords_to_positions(&coords);
    let dos_method = match method.to_ascii_lowercase().as_str() {
        "eht" => sci_form_core::dos::multi_method::DosMethod::Eht,
        "pm3" => sci_form_core::dos::multi_method::DosMethod::Pm3,
        "xtb" | "gfn0" => sci_form_core::dos::multi_method::DosMethod::Xtb,
        "gfn1" => sci_form_core::dos::multi_method::DosMethod::Gfn1,
        "gfn2" => sci_form_core::dos::multi_method::DosMethod::Gfn2,
        "hf3c" | "hf-3c" => sci_form_core::dos::multi_method::DosMethod::Hf3c,
        other => {
            return Err(value_error(invalid_value(
                "'method'",
                other,
                &["eht", "pm3", "xtb", "gfn1", "gfn2", "hf3c"],
            )))
        }
    };
    sci_form_core::dos::multi_method::compute_dos_multimethod(
        &elements, &positions, dos_method, sigma, e_min, e_max, n_points,
    )
    .map(|r| MultiMethodDosResultPy {
        energies: r.energies,
        total_dos: r.total_dos,
        method: format!("{:?}", r.method).to_lowercase(),
        homo_energy: r.homo_energy,
        lumo_energy: r.lumo_energy,
        gap: r.gap,
        orbital_energies: r.orbital_energies,
        sigma: r.sigma,
    })
    .map_err(engine_error)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_charges, m)?)?;
    m.add_function(wrap_pyfunction!(compute_charges_configured, m)?)?;
    m.add_function(wrap_pyfunction!(compute_sasa, m)?)?;
    m.add_function(wrap_pyfunction!(compute_dipole, m)?)?;
    m.add_function(wrap_pyfunction!(compute_dos, m)?)?;
    m.add_function(wrap_pyfunction!(compute_dos_multimethod, m)?)?;
    m.add_function(wrap_pyfunction!(compute_rmsd, m)?)?;
    m.add_class::<ChargeResultPy>()?;
    m.add_class::<SasaResultPy>()?;
    m.add_class::<DipoleResultPy>()?;
    m.add_class::<DosResultPy>()?;
    m.add_class::<MultiMethodDosResultPy>()?;
    m.add_class::<AlignmentResultPy>()?;
    Ok(())
}
