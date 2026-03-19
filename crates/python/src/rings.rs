//! Ring detection (SSSR), ECFP fingerprints, Tanimoto similarity, and
//! Butina clustering bindings.

use pyo3::prelude::*;

// ─── SSSR ──────────────────────────────────────────────────────────────────

#[pyclass]
#[derive(Clone)]
pub(crate) struct RingInfoPy {
    #[pyo3(get)]
    atoms: Vec<usize>,
    #[pyo3(get)]
    size: usize,
    #[pyo3(get)]
    is_aromatic: bool,
}

#[pyclass]
#[derive(Clone)]
pub(crate) struct SssrResultPy {
    #[pyo3(get)]
    rings: Vec<RingInfoPy>,
    #[pyo3(get)]
    atom_ring_count: Vec<usize>,
    #[pyo3(get)]
    atom_ring_sizes: Vec<Vec<usize>>,
    #[pyo3(get)]
    ring_size_histogram: Vec<usize>,
}

/// Compute the Smallest Set of Smallest Rings (SSSR) for a SMILES string.
#[pyfunction]
fn sssr(smiles: &str) -> PyResult<SssrResultPy> {
    sci_form_core::compute_sssr(smiles)
        .map(|r| SssrResultPy {
            rings: r
                .rings
                .iter()
                .map(|ri| RingInfoPy {
                    atoms: ri.atoms.clone(),
                    size: ri.size,
                    is_aromatic: ri.is_aromatic,
                })
                .collect(),
            atom_ring_count: r.atom_ring_count,
            atom_ring_sizes: r.atom_ring_sizes,
            ring_size_histogram: r.ring_size_histogram,
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

// ─── ECFP Fingerprints ─────────────────────────────────────────────────────

#[pyclass]
#[derive(Clone)]
pub(crate) struct ECFPFingerprintPy {
    #[pyo3(get)]
    n_bits: usize,
    #[pyo3(get)]
    on_bits: Vec<usize>,
    #[pyo3(get)]
    radius: usize,
    #[pyo3(get)]
    density: f64,
    /// Keep the inner rust struct for tanimoto computation.
    inner: sci_form_core::rings::ecfp::ECFPFingerprint,
}

/// Compute ECFP fingerprint for a molecule.
///
/// radius=2 → ECFP4, radius=3 → ECFP6. n_bits is the folded vector length.
#[pyfunction]
#[pyo3(signature = (smiles, radius=2, n_bits=2048))]
fn ecfp(smiles: &str, radius: usize, n_bits: usize) -> PyResult<ECFPFingerprintPy> {
    sci_form_core::compute_ecfp(smiles, radius, n_bits)
        .map(|fp| {
            let on_bits: Vec<usize> = fp.on_bits.iter().copied().collect();
            let density = fp.density();
            ECFPFingerprintPy {
                n_bits: fp.n_bits,
                on_bits,
                radius: fp.radius,
                density,
                inner: fp,
            }
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

/// Compute Tanimoto similarity between two ECFP fingerprints.
///
/// Returns a float in [0.0, 1.0].
#[pyfunction]
fn tanimoto(fp1: &ECFPFingerprintPy, fp2: &ECFPFingerprintPy) -> f64 {
    sci_form_core::compute_tanimoto(&fp1.inner, &fp2.inner)
}

// ─── Clustering ─────────────────────────────────────────────────────────────

#[pyclass]
#[derive(Clone)]
pub(crate) struct ClusterResultPy {
    #[pyo3(get)]
    n_clusters: usize,
    #[pyo3(get)]
    assignments: Vec<usize>,
    #[pyo3(get)]
    centroid_indices: Vec<usize>,
    #[pyo3(get)]
    cluster_sizes: Vec<usize>,
    #[pyo3(get)]
    rmsd_cutoff: f64,
}

/// Butina (Taylor-Butina) clustering on conformers by RMSD similarity.
///
/// conformers: list of flat coordinate arrays. rmsd_cutoff in Å.
#[pyfunction]
#[pyo3(signature = (conformers, rmsd_cutoff=1.0))]
fn butina_cluster(conformers: Vec<Vec<f64>>, rmsd_cutoff: f64) -> ClusterResultPy {
    let r = sci_form_core::butina_cluster(&conformers, rmsd_cutoff);
    ClusterResultPy {
        n_clusters: r.n_clusters,
        assignments: r.assignments,
        centroid_indices: r.centroid_indices,
        cluster_sizes: r.cluster_sizes,
        rmsd_cutoff: r.rmsd_cutoff,
    }
}

/// Compute the all-pairs RMSD matrix for a set of conformers.
///
/// Returns a list-of-lists of RMSD values.
#[pyfunction]
fn rmsd_matrix(conformers: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    sci_form_core::compute_rmsd_matrix(&conformers)
}

/// Filter conformer ensemble to keep only diverse cluster centroids.
///
/// Returns indices of representative conformers.
#[pyfunction]
#[pyo3(signature = (conformers, rmsd_cutoff=1.0))]
fn filter_diverse(conformers: Vec<Vec<f64>>, rmsd_cutoff: f64) -> Vec<usize> {
    sci_form_core::clustering::filter_diverse_conformers(&conformers, rmsd_cutoff)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // SSSR
    m.add_function(wrap_pyfunction!(sssr, m)?)?;
    m.add_class::<RingInfoPy>()?;
    m.add_class::<SssrResultPy>()?;
    // ECFP
    m.add_function(wrap_pyfunction!(ecfp, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto, m)?)?;
    m.add_class::<ECFPFingerprintPy>()?;
    // Clustering
    m.add_function(wrap_pyfunction!(butina_cluster, m)?)?;
    m.add_function(wrap_pyfunction!(rmsd_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(filter_diverse, m)?)?;
    m.add_class::<ClusterResultPy>()?;
    Ok(())
}
