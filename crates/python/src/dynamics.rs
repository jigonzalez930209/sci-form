//! Molecular dynamics, NEB pathways, and conformer search bindings.

use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
pub(crate) struct MdFramePy {
    #[pyo3(get)]
    step: usize,
    #[pyo3(get)]
    time_fs: f64,
    #[pyo3(get)]
    coords: Vec<f64>,
    #[pyo3(get)]
    potential_energy_kcal_mol: f64,
    #[pyo3(get)]
    kinetic_energy_kcal_mol: f64,
    #[pyo3(get)]
    temperature_k: f64,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct MdTrajectoryPy {
    #[pyo3(get)]
    frames: Vec<MdFramePy>,
    #[pyo3(get)]
    dt_fs: f64,
    #[pyo3(get)]
    notes: Vec<String>,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct NebImagePy {
    #[pyo3(get)]
    index: usize,
    #[pyo3(get)]
    coords: Vec<f64>,
    #[pyo3(get)]
    potential_energy_kcal_mol: f64,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct NebPathResultPy {
    #[pyo3(get)]
    images: Vec<NebImagePy>,
    #[pyo3(get)]
    notes: Vec<String>,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct ConformerEnsembleMemberPy {
    #[pyo3(get)]
    seed: u64,
    #[pyo3(get)]
    cluster_id: Option<usize>,
    #[pyo3(get)]
    coords: Vec<f64>,
    #[pyo3(get)]
    energy_kcal_mol: f64,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct ConformerClusterSummaryPy {
    #[pyo3(get)]
    cluster_id: usize,
    #[pyo3(get)]
    representative_seed: u64,
    #[pyo3(get)]
    size: usize,
    #[pyo3(get)]
    member_seeds: Vec<u64>,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct ConformerSearchResultPy {
    #[pyo3(get)]
    generated: usize,
    #[pyo3(get)]
    unique: usize,
    #[pyo3(get)]
    rotatable_bonds: usize,
    #[pyo3(get)]
    conformers: Vec<ConformerEnsembleMemberPy>,
    #[pyo3(get)]
    clusters: Vec<ConformerClusterSummaryPy>,
    #[pyo3(get)]
    notes: Vec<String>,
}

fn map_trajectory(r: sci_form_core::dynamics::MdTrajectory) -> MdTrajectoryPy {
    MdTrajectoryPy {
        frames: r
            .frames
            .into_iter()
            .map(|f| MdFramePy {
                step: f.step,
                time_fs: f.time_fs,
                coords: f.coords,
                potential_energy_kcal_mol: f.potential_energy_kcal_mol,
                kinetic_energy_kcal_mol: f.kinetic_energy_kcal_mol,
                temperature_k: f.temperature_k,
            })
            .collect(),
        dt_fs: r.dt_fs,
        notes: r.notes,
    }
}

/// Run NVE molecular dynamics with Velocity Verlet + UFF.
#[pyfunction]
#[pyo3(signature = (smiles, coords, n_steps=100, dt_fs=1.0, seed=42))]
fn md_trajectory(
    smiles: &str,
    coords: Vec<f64>,
    n_steps: usize,
    dt_fs: f64,
    seed: u64,
) -> PyResult<MdTrajectoryPy> {
    sci_form_core::compute_md_trajectory(smiles, &coords, n_steps, dt_fs, seed)
        .map(map_trajectory)
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

/// Run NVT molecular dynamics with Velocity Verlet + Berendsen thermostat + UFF.
#[pyfunction]
#[pyo3(signature = (smiles, coords, n_steps=100, dt_fs=1.0, seed=42, target_temp_k=300.0, thermostat_tau_fs=100.0))]
fn md_trajectory_nvt(
    smiles: &str,
    coords: Vec<f64>,
    n_steps: usize,
    dt_fs: f64,
    seed: u64,
    target_temp_k: f64,
    thermostat_tau_fs: f64,
) -> PyResult<MdTrajectoryPy> {
    sci_form_core::compute_md_trajectory_nvt(
        smiles,
        &coords,
        n_steps,
        dt_fs,
        seed,
        target_temp_k,
        thermostat_tau_fs,
    )
    .map(map_trajectory)
    .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

/// Compute simplified NEB path between two geometries.
#[pyfunction]
#[pyo3(signature = (smiles, start_coords, end_coords, n_images=8, n_iter=50, spring_k=0.1, step_size=0.01))]
fn simplified_neb_path(
    smiles: &str,
    start_coords: Vec<f64>,
    end_coords: Vec<f64>,
    n_images: usize,
    n_iter: usize,
    spring_k: f64,
    step_size: f64,
) -> PyResult<NebPathResultPy> {
    sci_form_core::compute_simplified_neb_path(
        smiles,
        &start_coords,
        &end_coords,
        n_images,
        n_iter,
        spring_k,
        step_size,
    )
    .map(|r| NebPathResultPy {
        images: r
            .images
            .into_iter()
            .map(|img| NebImagePy {
                index: img.index,
                coords: img.coords,
                potential_energy_kcal_mol: img.potential_energy_kcal_mol,
            })
            .collect(),
        notes: r.notes,
    })
    .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

/// Search conformers with UFF ranking and RMSD clustering.
#[pyfunction]
#[pyo3(signature = (smiles, n_samples=10, seed=42, rmsd_threshold=0.5))]
fn search_conformers(
    smiles: &str,
    n_samples: usize,
    seed: u64,
    rmsd_threshold: f64,
) -> PyResult<ConformerSearchResultPy> {
    sci_form_core::search_conformers_with_uff(smiles, n_samples, seed, rmsd_threshold)
        .map(|r| ConformerSearchResultPy {
            generated: r.generated,
            unique: r.unique,
            rotatable_bonds: r.rotatable_bonds,
            conformers: r
                .conformers
                .into_iter()
                .map(|c| ConformerEnsembleMemberPy {
                    seed: c.seed,
                    cluster_id: c.cluster_id,
                    coords: c.coords,
                    energy_kcal_mol: c.energy_kcal_mol,
                })
                .collect(),
            clusters: r
                .clusters
                .into_iter()
                .map(|c| ConformerClusterSummaryPy {
                    cluster_id: c.cluster_id,
                    representative_seed: c.representative_seed,
                    size: c.size,
                    member_seeds: c.member_seeds,
                })
                .collect(),
            notes: r.notes,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(md_trajectory, m)?)?;
    m.add_function(wrap_pyfunction!(md_trajectory_nvt, m)?)?;
    m.add_function(wrap_pyfunction!(simplified_neb_path, m)?)?;
    m.add_function(wrap_pyfunction!(search_conformers, m)?)?;
    m.add_class::<MdFramePy>()?;
    m.add_class::<MdTrajectoryPy>()?;
    m.add_class::<NebImagePy>()?;
    m.add_class::<NebPathResultPy>()?;
    m.add_class::<ConformerEnsembleMemberPy>()?;
    m.add_class::<ConformerClusterSummaryPy>()?;
    m.add_class::<ConformerSearchResultPy>()?;
    Ok(())
}
