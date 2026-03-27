//! Python bindings for all beta-tier experimental modules.
//!
//! Exposed as a sub-module of `sci_form`:
//!
//! ```python
//! from sci_form.beta import (
//!     kpm_dos, kpm_mulliken,
//!     mbh_frequencies,
//!     eht_randnla,
//!     psd_distance, psd_projection,
//!     cpm_charges,
//! )
//! ```

use pyo3::prelude::*;

// ─── B1: KPM ─────────────────────────────────────────────────────────────────

#[cfg(feature = "beta-kpm")]
mod kpm_py {
    use nalgebra::DMatrix;
    use pyo3::prelude::*;
    use sci_form_core::beta::kpm::{compute_kpm_dos, compute_kpm_mulliken, KpmConfig};

    fn build_eht_matrices(
        elements: &[u8],
        positions: &[[f64; 3]],
    ) -> (
        Vec<sci_form_core::eht::basis::AtomicOrbital>,
        DMatrix<f64>,
        DMatrix<f64>,
    ) {
        let basis = sci_form_core::eht::basis::build_basis(elements, positions);
        let overlap = sci_form_core::eht::build_overlap_matrix(&basis);
        let hamiltonian = sci_form_core::eht::build_hamiltonian(&basis, &overlap, None);
        (basis, hamiltonian, overlap)
    }

    #[pyclass]
    #[derive(Clone)]
    pub struct KpmDosResultPy {
        #[pyo3(get)]
        pub energies: Vec<f64>,
        #[pyo3(get)]
        pub total_dos: Vec<f64>,
        #[pyo3(get)]
        pub e_min: f64,
        #[pyo3(get)]
        pub e_max: f64,
        #[pyo3(get)]
        pub order: usize,
    }

    #[pyclass]
    #[derive(Clone)]
    pub struct KpmMullikenResultPy {
        #[pyo3(get)]
        pub mulliken_charges: Vec<f64>,
        #[pyo3(get)]
        pub orbital_populations: Vec<f64>,
    }

    /// Compute EHT KPM Density of States (O(N) scaling).
    ///
    /// Args:
    ///     elements (list[int]): Atomic numbers.
    ///     coords (list[float]): Flat coordinates in Å.
    ///     order (int, default 100): Chebyshev expansion order.
    ///     temperature (float, default 0): Electronic temperature in K.
    ///
    /// Returns:
    ///     KpmDosResultPy
    #[pyfunction]
    #[pyo3(signature = (elements, coords, order = 100, temperature = 0.0))]
    pub fn kpm_dos(
        elements: Vec<u8>,
        coords: Vec<f64>,
        order: usize,
        temperature: f64,
    ) -> PyResult<KpmDosResultPy> {
        let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let (_basis, hamiltonian, _overlap) = build_eht_matrices(&elements, &positions);
        let config = KpmConfig {
            order,
            temperature,
            n_vectors: 0,
            seed: 42,
        };
        let r = compute_kpm_dos(&hamiltonian, &config, -30.0, 5.0, 500);
        Ok(KpmDosResultPy {
            energies: r.energies,
            total_dos: r.total_dos,
            e_min: r.e_min,
            e_max: r.e_max,
            order: r.order,
        })
    }

    /// Compute EHT Mulliken populations using KPM stochastic trace.
    ///
    /// Returns:
    ///     KpmMullikenResultPy
    #[pyfunction]
    pub fn kpm_mulliken(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<KpmMullikenResultPy> {
        let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let (basis, hamiltonian, overlap) = build_eht_matrices(&elements, &positions);
        let eht = sci_form_core::eht::solve_eht(&elements, &positions)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e))?;
        let config = KpmConfig::default();
        let nuclear_charges: Vec<f64> = basis
            .iter()
            .map(|orbital| elements[orbital.atom_index] as f64)
            .collect();
        let r = compute_kpm_mulliken(
            &hamiltonian,
            &overlap,
            eht.n_electrons,
            &nuclear_charges,
            &config,
        );

        let mut mulliken_charges = vec![0.0; elements.len()];
        let mut orbital_populations = vec![0.0; basis.len()];
        for (idx, orbital) in basis.iter().enumerate() {
            let charge = r.charges.get(idx).copied().unwrap_or(0.0);
            mulliken_charges[orbital.atom_index] += charge;
            orbital_populations[idx] = nuclear_charges[idx] - charge;
        }
        Ok(KpmMullikenResultPy {
            mulliken_charges,
            orbital_populations,
        })
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<KpmDosResultPy>()?;
        m.add_class::<KpmMullikenResultPy>()?;
        m.add_function(wrap_pyfunction!(kpm_dos, m)?)?;
        m.add_function(wrap_pyfunction!(kpm_mulliken, m)?)?;
        Ok(())
    }
}

// ─── B2: MBH ─────────────────────────────────────────────────────────────────

#[cfg(feature = "beta-mbh")]
mod mbh_py {
    use pyo3::prelude::*;
    use sci_form_core::beta::mbh::hessian::{compute_mbh_frequencies, MbhConfig};

    #[pyclass]
    #[derive(Clone)]
    pub struct MbhResultPy {
        #[pyo3(get)]
        pub frequencies: Vec<f64>,
        #[pyo3(get)]
        pub n_blocks: usize,
        #[pyo3(get)]
        pub n_flexible: usize,
        #[pyo3(get)]
        pub n_dof_reduced: usize,
        #[pyo3(get)]
        pub n_dof_full: usize,
        #[pyo3(get)]
        pub speedup: f64,
    }

    /// Compute MBH vibrational frequencies using UFF energy.
    ///
    /// Args:
    ///     elements (list[int]): Atomic numbers.
    ///     coords (list[float]): Flat coordinates in Å.
    ///     smiles (str): SMILES string for topology (ring detection and force field).
    ///     fd_step (float, default 0.005): Finite difference step in Å.
    ///
    /// Returns:
    ///     MbhResultPy
    #[pyfunction]
    #[pyo3(signature = (elements, coords, smiles, fd_step = 0.005))]
    pub fn mbh_frequencies(
        elements: Vec<u8>,
        coords: Vec<f64>,
        smiles: &str,
        fd_step: f64,
    ) -> PyResult<MbhResultPy> {
        let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let rings: Vec<(Vec<usize>, bool)> = match sci_form_core::compute_sssr(smiles) {
            Ok(s) => s
                .rings
                .iter()
                .map(|r| (r.atoms.clone(), r.is_aromatic))
                .collect(),
            Err(_) => vec![],
        };
        let config = MbhConfig { fd_step };
        let owned_smiles = smiles.to_string();
        let energy_fn = move |c: &[f64]| -> f64 {
            sci_form_core::compute_uff_energy(&owned_smiles, c).unwrap_or(0.0)
        };
        let r = compute_mbh_frequencies(&elements, &positions, &rings, &energy_fn, &config);
        Ok(MbhResultPy {
            frequencies: r.frequencies,
            n_blocks: r.n_blocks,
            n_flexible: r.n_flexible,
            n_dof_reduced: r.n_dof_reduced,
            n_dof_full: r.n_dof_full,
            speedup: r.speedup,
        })
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<MbhResultPy>()?;
        m.add_function(wrap_pyfunction!(mbh_frequencies, m)?)?;
        Ok(())
    }
}

// ─── B3: RandNLA ─────────────────────────────────────────────────────────────

#[cfg(feature = "beta-randnla")]
mod randnla_py {
    use nalgebra::DMatrix;
    use pyo3::prelude::*;
    use sci_form_core::beta::rand_nla::{solve_eht_randnla, RandNlaConfig};

    fn build_eht_matrices(
        elements: &[u8],
        positions: &[[f64; 3]],
    ) -> (DMatrix<f64>, DMatrix<f64>, usize) {
        let basis = sci_form_core::eht::basis::build_basis(elements, positions);
        let overlap = sci_form_core::eht::build_overlap_matrix(&basis);
        let hamiltonian = sci_form_core::eht::build_hamiltonian(&basis, &overlap, None);
        let n_electrons = sci_form_core::eht::solve_eht(elements, positions)
            .map(|result| result.n_electrons)
            .unwrap_or(elements.len() * 2);
        (hamiltonian, overlap, n_electrons)
    }

    #[pyclass]
    #[derive(Clone)]
    pub struct RandNlaResultPy {
        #[pyo3(get)]
        pub orbital_energies: Vec<f64>,
        #[pyo3(get)]
        pub homo_index: usize,
        #[pyo3(get)]
        pub homo_energy: f64,
        #[pyo3(get)]
        pub lumo_energy: f64,
        #[pyo3(get)]
        pub gap: f64,
        #[pyo3(get)]
        pub k: usize,
        #[pyo3(get)]
        pub residual_error: f64,
        #[pyo3(get)]
        pub used_fallback: bool,
    }

    /// Solve EHT eigenvalue problem using randomized Nyström approximation.
    ///
    /// Scales as O(N k²) instead of O(N³), where k ≈ √N.
    ///
    /// Args:
    ///     elements (list[int]): Atomic numbers.
    ///     coords (list[float]): Flat coordinates in Å.
    ///     sketch_size (int, optional): Nyström rank k. Default: √N.
    ///     max_error (float, default 0.001): Max relative residual.
    ///
    /// Returns:
    ///     RandNlaResultPy
    #[pyfunction]
    #[pyo3(signature = (elements, coords, sketch_size = None, max_error = 0.001))]
    pub fn eht_randnla(
        elements: Vec<u8>,
        coords: Vec<f64>,
        sketch_size: Option<usize>,
        max_error: f64,
    ) -> PyResult<RandNlaResultPy> {
        let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let config = RandNlaConfig {
            sketch_size,
            max_error,
            seed: 42,
            fallback_enabled: true,
        };
        let (hamiltonian, overlap, n_electrons) = build_eht_matrices(&elements, &positions);
        let (eigenvalues, _, info) = solve_eht_randnla(&hamiltonian, &overlap, &config);
        let energies: Vec<f64> = eigenvalues.iter().cloned().collect();
        let n_occ = n_electrons / 2;
        let homo_idx = n_occ.saturating_sub(1);
        let lumo_idx = n_occ.min(energies.len().saturating_sub(1));
        let homo = energies.get(homo_idx).copied().unwrap_or(0.0);
        let lumo = energies.get(lumo_idx).copied().unwrap_or(0.0);
        Ok(RandNlaResultPy {
            orbital_energies: energies,
            homo_index: homo_idx,
            homo_energy: homo,
            lumo_energy: lumo,
            gap: (lumo - homo).max(0.0),
            k: info.k,
            residual_error: info.residual_error,
            used_fallback: info.used_fallback,
        })
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<RandNlaResultPy>()?;
        m.add_function(wrap_pyfunction!(eht_randnla, m)?)?;
        Ok(())
    }
}

// ─── B4: Riemannian ──────────────────────────────────────────────────────────

#[cfg(feature = "beta-riemannian")]
mod riemannian_py {
    use nalgebra::DMatrix;
    use pyo3::prelude::*;
    use sci_form_core::beta::riemannian::{psd_distance, psd_projection};

    /// Compute Riemannian distance between two PSD matrices.
    ///
    /// Args:
    ///     matrix_a (list[float]): Row-major flat matrix.
    ///     matrix_b (list[float]): Row-major flat matrix.
    ///     dim (int): Dimension of both matrices.
    ///
    /// Returns:
    ///     float
    #[pyfunction]
    pub fn beta_psd_distance(matrix_a: Vec<f64>, matrix_b: Vec<f64>, dim: usize) -> PyResult<f64> {
        if matrix_a.len() != dim * dim || matrix_b.len() != dim * dim {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "matrix flat length must equal dim*dim",
            ));
        }
        let a = DMatrix::from_row_slice(dim, dim, &matrix_a);
        let b = DMatrix::from_row_slice(dim, dim, &matrix_b);
        Ok(psd_distance(&a, &b))
    }

    /// Project a symmetric matrix onto the PSD cone.
    ///
    /// Args:
    ///     matrix (list[float]): Row-major flat matrix.
    ///     dim (int): Dimension.
    ///
    /// Returns:
    ///     list[float] — projected flat matrix
    #[pyfunction]
    pub fn beta_psd_projection(matrix: Vec<f64>, dim: usize) -> PyResult<Vec<f64>> {
        if matrix.len() != dim * dim {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "matrix flat length must equal dim*dim",
            ));
        }
        let m = DMatrix::from_row_slice(dim, dim, &matrix);
        let p = psd_projection(&m);
        Ok(p.as_slice().to_vec())
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(beta_psd_distance, m)?)?;
        m.add_function(wrap_pyfunction!(beta_psd_projection, m)?)?;
        Ok(())
    }
}

// ─── B5: CPM ─────────────────────────────────────────────────────────────────

#[cfg(feature = "beta-cpm")]
mod cpm_py {
    use pyo3::prelude::*;
    use sci_form_core::beta::cpm::grand_potential::{compute_cpm_charges, CpmConfig};

    #[pyclass]
    #[derive(Clone)]
    pub struct CpmResultPy {
        #[pyo3(get)]
        pub charges: Vec<f64>,
        #[pyo3(get)]
        pub total_charge: f64,
        #[pyo3(get)]
        pub energy: f64,
    }

    /// Compute constant-potential CPM charges.
    ///
    /// Args:
    ///     elements (list[int]): Atomic numbers.
    ///     coords (list[float]): Flat coordinates in Å.
    ///     potential (float): Electrode potential in eV.
    ///
    /// Returns:
    ///     CpmResultPy
    #[pyfunction]
    #[pyo3(signature = (elements, coords, potential = 0.0))]
    pub fn cpm_charges(
        elements: Vec<u8>,
        coords: Vec<f64>,
        potential: f64,
    ) -> PyResult<CpmResultPy> {
        let config = CpmConfig {
            mu_ev: potential,
            ..Default::default()
        };
        let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let r = compute_cpm_charges(&elements, &positions, &config);
        Ok(CpmResultPy {
            charges: r.charges,
            total_charge: r.total_charge,
            energy: r.grand_potential,
        })
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<CpmResultPy>()?;
        m.add_function(wrap_pyfunction!(cpm_charges, m)?)?;
        Ok(())
    }
}

// ─── Module registration ─────────────────────────────────────────────────────

pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new_bound(parent.py(), "beta")?;
    m.add("__doc__", "Beta-tier experimental modules for sci_form.\n\nThese APIs have passed validation but may still change before stable release.")?;

    #[cfg(feature = "beta-kpm")]
    kpm_py::register(&m)?;
    #[cfg(feature = "beta-mbh")]
    mbh_py::register(&m)?;
    #[cfg(feature = "beta-randnla")]
    randnla_py::register(&m)?;
    #[cfg(feature = "beta-riemannian")]
    riemannian_py::register(&m)?;
    #[cfg(feature = "beta-cpm")]
    cpm_py::register(&m)?;

    parent.add_submodule(&m)?;
    Ok(())
}
