//! Python bindings for all alpha-tier experimental modules.
//!
//! These are exposed as a sub-module of `sci_form`:
//!
//! ```python
//! from sci_form.alpha import (
//!     dft_calculate, reaxff_gradient, reaxff_energy,
//!     eem_charges, compute_mlff, compute_aevs,
//!     boys_function, eri_ssss,
//!     rotate_dihedral_cga, refine_torsion_cga,
//!     gsm_interpolate, gsm_find_ts,
//!     sdr_embed,
//! )
//! ```

use pyo3::prelude::*;

// ─── A1: Kohn-Sham DFT ──────────────────────────────────────────────────────

#[cfg(feature = "alpha-dft")]
mod dft_py {
    use pyo3::prelude::*;
    use sci_form_core::dft::ks_fock::{solve_ks_dft, DftConfig, DftMethod};

    #[pyclass]
    #[derive(Clone)]
    pub struct DftResultPy {
        #[pyo3(get)]
        pub energy: f64,
        #[pyo3(get)]
        pub homo_energy: f64,
        #[pyo3(get)]
        pub lumo_energy: f64,
        #[pyo3(get)]
        pub gap: f64,
        #[pyo3(get)]
        pub converged: bool,
        #[pyo3(get)]
        pub n_basis: usize,
        #[pyo3(get)]
        pub scf_iterations: usize,
        #[pyo3(get)]
        pub mulliken_charges: Vec<f64>,
        #[pyo3(get)]
        pub orbital_energies: Vec<f64>,
        #[pyo3(get)]
        pub nuclear_repulsion: f64,
        #[pyo3(get)]
        pub xc_energy: f64,
    }

    /// Run a Kohn-Sham DFT single-point calculation.
    ///
    /// Args:
    ///     elements (list[int]): Atomic numbers.
    ///     coords (list[float]): Flat 3D coordinates [x0,y0,z0,...] in Å.
    ///     method (str): ``"svwn"`` (LDA) or ``"pbe"`` (GGA, default).
    ///
    /// Returns:
    ///     DftResultPy
    #[pyfunction]
    #[pyo3(signature = (elements, coords, method = "pbe"))]
    pub fn dft_calculate(
        elements: Vec<u8>,
        coords: Vec<f64>,
        method: &str,
    ) -> PyResult<DftResultPy> {
        let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let dft_method = match method.to_lowercase().as_str() {
            "svwn" | "lda" => DftMethod::Svwn,
            _ => DftMethod::Pbe,
        };
        let config = DftConfig {
            method: dft_method,
            ..Default::default()
        };
        let r = solve_ks_dft(&elements, &positions, &config)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e))?;
        Ok(DftResultPy {
            energy: r.total_energy_ev,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap: r.gap,
            converged: r.converged,
            n_basis: r.n_basis,
            scf_iterations: r.scf_iterations,
            mulliken_charges: r.mulliken_charges,
            orbital_energies: r.orbital_energies,
            nuclear_repulsion: r.nuclear_repulsion,
            xc_energy: r.xc_energy,
        })
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<DftResultPy>()?;
        m.add_function(wrap_pyfunction!(dft_calculate, m)?)?;
        Ok(())
    }
}

// ─── A2: ReaxFF ──────────────────────────────────────────────────────────────

#[cfg(feature = "alpha-reaxff")]
mod reaxff_py {
    use pyo3::prelude::*;
    use sci_form_core::forcefield::reaxff::{
        bond_order::compute_bond_orders,
        eem::solve_eem,
        energy::compute_bonded_energy,
        gradients::compute_reaxff_gradient,
        nonbonded::compute_nonbonded_energy,
        params::{ReaxffAtomParams, ReaxffParams},
    };

    fn build_atom_params(elements: &[u8], params: &ReaxffParams) -> Vec<ReaxffAtomParams> {
        let fallback = params
            .atom_params
            .last()
            .cloned()
            .unwrap_or(ReaxffAtomParams {
                element: 1,
                r_sigma: 1.0,
                r_pi: 0.0,
                r_pipi: 0.0,
                p_bo1: 0.0,
                p_bo2: 1.0,
                p_bo3: 0.0,
                p_bo4: 1.0,
                p_bo5: 0.0,
                p_bo6: 1.0,
                valence: 1.0,
            });

        elements
            .iter()
            .map(|&z| {
                params
                    .element_index(z)
                    .map(|idx| params.atom_params[idx].clone())
                    .unwrap_or_else(|| fallback.clone())
            })
            .collect()
    }

    #[pyclass]
    #[derive(Clone)]
    pub struct ReaxffGradientResultPy {
        #[pyo3(get)]
        pub energy_kcal_mol: f64,
        #[pyo3(get)]
        pub gradient: Vec<f64>,
    }

    #[pyclass]
    #[derive(Clone)]
    pub struct ReaxffEnergyResultPy {
        #[pyo3(get)]
        pub bonded: f64,
        #[pyo3(get)]
        pub coulomb: f64,
        #[pyo3(get)]
        pub van_der_waals: f64,
        #[pyo3(get)]
        pub total: f64,
    }

    #[pyclass]
    #[derive(Clone)]
    pub struct EemChargesResultPy {
        #[pyo3(get)]
        pub charges: Vec<f64>,
        #[pyo3(get)]
        pub total_charge: f64,
    }

    /// Compute ReaxFF numerical gradient and total energy.
    ///
    /// Args:
    ///     elements (list[int]): Atomic numbers.
    ///     coords (list[float]): Flat coordinates in Å.
    ///
    /// Returns:
    ///     ReaxffGradientResultPy: energy_kcal_mol, gradient
    #[pyfunction]
    pub fn reaxff_gradient(
        elements: Vec<u8>,
        coords: Vec<f64>,
    ) -> PyResult<ReaxffGradientResultPy> {
        let params = ReaxffParams::default_chon();
        let (energy, gradient) = compute_reaxff_gradient(&coords, &elements, &params)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e))?;
        Ok(ReaxffGradientResultPy {
            energy_kcal_mol: energy,
            gradient,
        })
    }

    /// Compute ReaxFF energy breakdown (no gradient).
    ///
    /// Returns:
    ///     ReaxffEnergyResultPy: bonded, coulomb, van_der_waals, total
    #[pyfunction]
    pub fn reaxff_energy(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<ReaxffEnergyResultPy> {
        let params = ReaxffParams::default_chon();
        let atom_params = build_atom_params(&elements, &params);
        let eem_params: Vec<_> = elements
            .iter()
            .map(|&z| sci_form_core::forcefield::reaxff::eem::default_eem_params(z))
            .collect();
        let charges = solve_eem(&coords, &eem_params, 0.0)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e))?;
        let bo = compute_bond_orders(&coords, &atom_params, params.cutoff);
        let bonded = compute_bonded_energy(&coords, &bo, &params);
        let (van_der_waals, coulomb) =
            compute_nonbonded_energy(&coords, &charges, &elements, params.cutoff);
        Ok(ReaxffEnergyResultPy {
            bonded: bonded.total,
            coulomb,
            van_der_waals,
            total: bonded.total + coulomb + van_der_waals,
        })
    }

    /// Compute EEM equilibrium charges.
    ///
    /// Returns:
    ///     EemChargesResultPy: charges, total_charge
    #[pyfunction]
    pub fn eem_charges(elements: Vec<u8>, coords: Vec<f64>) -> EemChargesResultPy {
        let eem_params: Vec<_> = elements
            .iter()
            .map(|&z| sci_form_core::forcefield::reaxff::eem::default_eem_params(z))
            .collect();
        let charges = solve_eem(&coords, &eem_params, 0.0).unwrap_or_default();
        let total: f64 = charges.iter().sum();
        EemChargesResultPy {
            charges,
            total_charge: total,
        }
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<ReaxffGradientResultPy>()?;
        m.add_class::<ReaxffEnergyResultPy>()?;
        m.add_class::<EemChargesResultPy>()?;
        m.add_function(wrap_pyfunction!(reaxff_gradient, m)?)?;
        m.add_function(wrap_pyfunction!(reaxff_energy, m)?)?;
        m.add_function(wrap_pyfunction!(eem_charges, m)?)?;
        Ok(())
    }
}

// ─── A3: MLFF ────────────────────────────────────────────────────────────────

#[cfg(feature = "alpha-mlff")]
mod mlff_py {
    use pyo3::prelude::*;
    use sci_form_core::mlff::{compute_aevs, compute_mlff, MlffConfig, SymmetryFunctionParams};

    #[pyclass]
    #[derive(Clone)]
    pub struct MlffResultPy {
        #[pyo3(get)]
        pub energy: f64,
        #[pyo3(get)]
        pub atomic_energies: Vec<f64>,
        #[pyo3(get)]
        pub forces: Vec<f64>,
    }

    #[pyclass]
    #[derive(Clone)]
    pub struct AevResultPy {
        #[pyo3(get)]
        pub n_atoms: usize,
        #[pyo3(get)]
        pub aev_length: usize,
        #[pyo3(get)]
        pub aevs: Vec<Vec<f64>>,
    }

    /// Compute MLFF energy and forces using element-specific neural networks.
    ///
    /// `config_dict` — Python dict matching `MlffConfig` (serializable to JSON).
    ///
    /// Returns:
    ///     MlffResultPy: energy, atomic_energies, forces
    #[pyfunction]
    pub fn alpha_compute_mlff(
        elements: Vec<u8>,
        coords: Vec<f64>,
        config_dict: &pyo3::types::PyDict,
    ) -> PyResult<MlffResultPy> {
        let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let config_value = config_dict
            .extract::<serde_json::Value>()
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("bad config: {}", e)))?;
        let config: MlffConfig = serde_json::from_value(config_value)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("bad config: {}", e)))?;
        let r = compute_mlff(&elements, &positions, &config)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e))?;
        Ok(MlffResultPy {
            energy: r.energy,
            atomic_energies: r.atomic_energies,
            forces: r.forces.into_iter().flatten().collect(),
        })
    }

    /// Compute Atomic Environment Vectors for each atom.
    ///
    /// Returns:
    ///     AevResultPy: n_atoms, aev_length, aevs
    #[pyfunction]
    #[pyo3(signature = (elements, coords, radial_cutoff = 6.0, angular_cutoff = 3.5))]
    pub fn alpha_compute_aevs(
        elements: Vec<u8>,
        coords: Vec<f64>,
        radial_cutoff: f64,
        angular_cutoff: f64,
    ) -> AevResultPy {
        let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let params = SymmetryFunctionParams {
            radial_cutoff,
            angular_cutoff,
            ..Default::default()
        };
        let aevs = compute_aevs(&elements, &positions, &params);
        let aev_length = aevs.first().map(|v| v.len()).unwrap_or(0);
        let aev_arrays: Vec<Vec<f64>> = aevs.iter().map(|v| v.to_vec()).collect();
        AevResultPy {
            n_atoms: elements.len(),
            aev_length,
            aevs: aev_arrays,
        }
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<MlffResultPy>()?;
        m.add_class::<AevResultPy>()?;
        m.add_function(wrap_pyfunction!(alpha_compute_mlff, m)?)?;
        m.add_function(wrap_pyfunction!(alpha_compute_aevs, m)?)?;
        Ok(())
    }
}

// ─── A4: Obara-Saika ERIs ───────────────────────────────────────────────────

#[cfg(feature = "alpha-obara-saika")]
mod obara_saika_py {
    use pyo3::prelude::*;
    use sci_form_core::scf::obara_saika::{boys_function, eri_ssss, schwarz_bound, ShellPairData};

    /// Compute the Boys function F_n(x).
    ///
    /// Returns:
    ///     float
    #[pyfunction]
    pub fn alpha_boys_function(n: usize, x: f64) -> f64 {
        boys_function(n, x)
    }

    /// Compute a (ss|ss) two-electron repulsion integral using Obara-Saika.
    ///
    /// Args:
    ///     alpha_a, center_a, alpha_b, center_b: shell A and B
    ///     alpha_c, center_c, alpha_d, center_d: shell C and D
    ///
    /// Returns:
    ///     float — ERI value in Hartree
    #[pyfunction]
    #[allow(clippy::too_many_arguments)]
    pub fn alpha_eri_ssss(
        alpha_a: f64,
        center_a: [f64; 3],
        alpha_b: f64,
        center_b: [f64; 3],
        alpha_c: f64,
        center_c: [f64; 3],
        alpha_d: f64,
        center_d: [f64; 3],
    ) -> f64 {
        let sp_ab = ShellPairData::new(alpha_a, center_a, alpha_b, center_b);
        let sp_cd = ShellPairData::new(alpha_c, center_c, alpha_d, center_d);
        eri_ssss(&sp_ab, &sp_cd)
    }

    /// Compute Schwarz pre-screening bound Q_ab = √〈ab|ab〉.
    ///
    /// Returns:
    ///     float
    #[pyfunction]
    pub fn alpha_schwarz_bound(
        alpha: f64,
        center_a: [f64; 3],
        beta: f64,
        center_b: [f64; 3],
    ) -> f64 {
        schwarz_bound(alpha, center_a, beta, center_b)
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(alpha_boys_function, m)?)?;
        m.add_function(wrap_pyfunction!(alpha_eri_ssss, m)?)?;
        m.add_function(wrap_pyfunction!(alpha_schwarz_bound, m)?)?;
        Ok(())
    }
}

// ─── A5: CGA Motors ─────────────────────────────────────────────────────────

#[cfg(feature = "alpha-cga")]
mod cga_py {
    use pyo3::prelude::*;
    use sci_form_core::alpha::cga;

    fn dihedral_angle(p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], p4: [f64; 3]) -> f64 {
        let b1 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
        let b2 = [p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]];
        let b3 = [p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2]];

        let n1 = [
            b1[1] * b2[2] - b1[2] * b2[1],
            b1[2] * b2[0] - b1[0] * b2[2],
            b1[0] * b2[1] - b1[1] * b2[0],
        ];
        let n2 = [
            b2[1] * b3[2] - b2[2] * b3[1],
            b2[2] * b3[0] - b2[0] * b3[2],
            b2[0] * b3[1] - b2[1] * b3[0],
        ];
        let b2_norm = (b2[0] * b2[0] + b2[1] * b2[1] + b2[2] * b2[2]).sqrt();
        let b2_unit = if b2_norm > 1e-15 {
            [b2[0] / b2_norm, b2[1] / b2_norm, b2[2] / b2_norm]
        } else {
            [0.0, 0.0, 0.0]
        };
        let m1 = [
            n1[1] * b2_unit[2] - n1[2] * b2_unit[1],
            n1[2] * b2_unit[0] - n1[0] * b2_unit[2],
            n1[0] * b2_unit[1] - n1[1] * b2_unit[0],
        ];

        let x = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2];
        let y = m1[0] * n2[0] + m1[1] * n2[1] + m1[2] * n2[2];
        (-y).atan2(-x) + std::f64::consts::PI
    }

    fn rotation_subtree(
        bond_a: usize,
        bond_b: usize,
        bonds: &[(usize, usize)],
        n_atoms: usize,
    ) -> Vec<usize> {
        let mut adjacency = vec![Vec::new(); n_atoms];
        for &(i, j) in bonds {
            if i < n_atoms && j < n_atoms {
                adjacency[i].push(j);
                adjacency[j].push(i);
            }
        }

        let mut visited = vec![false; n_atoms];
        visited[bond_a] = true;
        visited[bond_b] = true;

        let mut stack = vec![bond_b];
        let mut subtree = vec![bond_b];

        while let Some(current) = stack.pop() {
            for &neighbor in &adjacency[current] {
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    subtree.push(neighbor);
                    stack.push(neighbor);
                }
            }
        }

        subtree
    }

    /// Apply a dihedral rotation to a sub-tree of atoms using CGA motors.
    ///
    /// Args:
    ///     coords (list[float]): Flat coordinates [x0,y0,z0,...].
    ///     subtree_indices (list[int]): Atom indices to rotate.
    ///     axis_a, axis_b (list[float]): [x,y,z] of the two axis atoms.
    ///     angle_rad (float): Rotation angle in radians.
    ///
    /// Returns:
    ///     list[float] — new flat coordinates
    #[pyfunction]
    pub fn rotate_dihedral_cga(
        mut coords: Vec<f64>,
        subtree_indices: Vec<usize>,
        axis_a: [f64; 3],
        axis_b: [f64; 3],
        angle_rad: f64,
    ) -> Vec<f64> {
        let motor = cga::dihedral_motor(axis_a, axis_b, angle_rad);
        cga::apply_motor_to_subtree(&coords, &subtree_indices, &motor)
    }

    /// Refine a torsion angle to a target value using CGA.
    ///
    /// Args:
    ///     elements (list[int]): Atomic numbers.
    ///     coords (list[float]): Flat coordinates.
    ///     bonds (list[tuple[int,int]]): Bond pairs.
    ///     torsion_indices (list[int]): 4 atom indices [i,j,k,l].
    ///     target_angle_rad (float): Target dihedral in radians.
    ///
    /// Returns:
    ///     list[float] — new flat coordinates
    #[pyfunction]
    pub fn refine_torsion_cga(
        elements: Vec<u8>,
        coords: Vec<f64>,
        bonds: Vec<(usize, usize)>,
        torsion_indices: Vec<usize>,
        target_angle_rad: f64,
    ) -> PyResult<Vec<f64>> {
        if torsion_indices.len() != 4 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "torsion_indices must have exactly 4 elements",
            ));
        }
        if coords.len() != elements.len() * 3 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "coords length must equal 3 × number of elements",
            ));
        }

        let i = torsion_indices[0];
        let j = torsion_indices[1];
        let k = torsion_indices[2];
        let l = torsion_indices[3];
        let current_angle = dihedral_angle(
            [coords[i * 3], coords[i * 3 + 1], coords[i * 3 + 2]],
            [coords[j * 3], coords[j * 3 + 1], coords[j * 3 + 2]],
            [coords[k * 3], coords[k * 3 + 1], coords[k * 3 + 2]],
            [coords[l * 3], coords[l * 3 + 1], coords[l * 3 + 2]],
        );

        let mut delta = target_angle_rad - current_angle;
        while delta <= -std::f64::consts::PI {
            delta += 2.0 * std::f64::consts::PI;
        }
        while delta > std::f64::consts::PI {
            delta -= 2.0 * std::f64::consts::PI;
        }

        let subtree = rotation_subtree(j, k, &bonds, elements.len());
        let motor = cga::dihedral_motor(
            [coords[j * 3], coords[j * 3 + 1], coords[j * 3 + 2]],
            [coords[k * 3], coords[k * 3 + 1], coords[k * 3 + 2]],
            delta,
        );
        Ok(cga::apply_motor_to_subtree(&coords, &subtree, &motor))
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(rotate_dihedral_cga, m)?)?;
        m.add_function(wrap_pyfunction!(refine_torsion_cga, m)?)?;
        Ok(())
    }
}

// ─── A6: GSM ────────────────────────────────────────────────────────────────

#[cfg(feature = "alpha-gsm")]
mod gsm_py {
    use pyo3::prelude::*;
    use sci_form_core::alpha::gsm::{self, GsmConfig};

    #[pyclass]
    #[derive(Clone)]
    pub struct GsmResultPy {
        #[pyo3(get)]
        pub ts_coords: Vec<f64>,
        #[pyo3(get)]
        pub ts_energy: f64,
        #[pyo3(get)]
        pub path_energies: Vec<f64>,
        #[pyo3(get)]
        pub reverse_barrier: f64,
        #[pyo3(get)]
        pub path_coords: Vec<Vec<f64>>,
        #[pyo3(get)]
        pub ts_node_index: usize,
        #[pyo3(get)]
        pub n_nodes: usize,
        #[pyo3(get)]
        pub energy_evaluations: usize,
    }

    /// Interpolate reaction path node at parameter t ∈ [0, 1].
    ///
    /// Returns:
    ///     list[float] — flat coordinates
    #[pyfunction]
    pub fn gsm_interpolate(reactant: Vec<f64>, product: Vec<f64>, t: f64) -> Vec<f64> {
        gsm::interpolate_node(&reactant, &product, t)
    }

    /// Find transition state using Growing String Method with UFF energy.
    ///
    /// Returns:
    ///     GsmResultPy
    #[pyfunction]
    #[pyo3(signature = (smiles, reactant_coords, product_coords, n_nodes = 9))]
    pub fn gsm_find_ts(
        smiles: &str,
        reactant_coords: Vec<f64>,
        product_coords: Vec<f64>,
        n_nodes: usize,
    ) -> PyResult<GsmResultPy> {
        let config = GsmConfig {
            max_nodes: n_nodes,
            ..Default::default()
        };
        let owned_smiles = smiles.to_string();
        let energy_fn = move |coords: &[f64]| -> f64 {
            sci_form_core::compute_uff_energy(&owned_smiles, coords).unwrap_or(f64::MAX)
        };
        let r = gsm::find_transition_state(&reactant_coords, &product_coords, &energy_fn, &config);
        Ok(GsmResultPy {
            ts_coords: r.ts_coords,
            ts_energy: r.ts_energy,
            path_energies: r.path_energies,
            reverse_barrier: r.reverse_barrier,
            path_coords: r.path_coords,
            ts_node_index: r.ts_node_index,
            n_nodes: r.n_nodes,
            energy_evaluations: r.energy_evaluations,
        })
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<GsmResultPy>()?;
        m.add_function(wrap_pyfunction!(gsm_interpolate, m)?)?;
        m.add_function(wrap_pyfunction!(gsm_find_ts, m)?)?;
        Ok(())
    }
}

// ─── A7: SDR Embedding ──────────────────────────────────────────────────────

#[cfg(feature = "alpha-sdr")]
mod sdr_py {
    use pyo3::prelude::*;
    use sci_form_core::alpha::sdr::{sdr_embed, SdrConfig};

    #[pyclass]
    #[derive(Clone)]
    pub struct SdrConvergencePy {
        #[pyo3(get)]
        pub iterations: usize,
        #[pyo3(get)]
        pub converged: bool,
        #[pyo3(get)]
        pub final_residual: f64,
        #[pyo3(get)]
        pub neg_eigenvalues_removed: usize,
    }

    #[pyclass]
    #[derive(Clone)]
    pub struct SdrResultPy {
        #[pyo3(get)]
        pub coords: Vec<f64>,
        #[pyo3(get)]
        pub num_atoms: usize,
        #[pyo3(get)]
        pub convergence: SdrConvergencePy,
        #[pyo3(get)]
        pub max_distance_error: f64,
        #[pyo3(get)]
        pub retries_avoided: usize,
    }

    /// Embed a molecule from pairwise distances using Semidefinite Relaxation.
    ///
    /// Args:
    ///     distance_pairs (list[tuple[int, int, float]]): (i, j, d_ij) triples.
    ///     n_atoms (int): Total number of atoms.
    ///     max_iterations (int, default 500)
    ///     tolerance (float, default 1e-6)
    ///
    /// Returns:
    ///     SdrResultPy
    #[pyfunction]
    #[pyo3(signature = (distance_pairs, n_atoms, max_iterations = 500, tolerance = 1e-6))]
    pub fn sdr_embed_py(
        distance_pairs: Vec<(usize, usize, f64)>,
        n_atoms: usize,
        max_iterations: usize,
        tolerance: f64,
    ) -> SdrResultPy {
        let config = SdrConfig {
            max_iter: max_iterations,
            tol: tolerance,
            ..Default::default()
        };
        let result = sdr_embed(n_atoms, &distance_pairs, &config);
        SdrResultPy {
            coords: result.coords,
            num_atoms: result.num_atoms,
            convergence: SdrConvergencePy {
                iterations: result.convergence.iterations,
                converged: result.convergence.converged,
                final_residual: result.convergence.final_residual,
                neg_eigenvalues_removed: result.convergence.neg_eigenvalues_removed,
            },
            max_distance_error: result.max_distance_error,
            retries_avoided: result.retries_avoided,
        }
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<SdrConvergencePy>()?;
        m.add_class::<SdrResultPy>()?;
        m.add_function(wrap_pyfunction!(sdr_embed_py, m)?)?;
        Ok(())
    }
}

// ─── Module registration ─────────────────────────────────────────────────────

pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new_bound(parent.py(), "alpha")?;
    m.add("__doc__", "Alpha-tier experimental modules for sci_form.\n\nThese APIs are early-stage and subject to breaking changes.")?;

    #[cfg(feature = "alpha-dft")]
    dft_py::register(&m)?;
    #[cfg(feature = "alpha-reaxff")]
    reaxff_py::register(&m)?;
    #[cfg(feature = "alpha-mlff")]
    mlff_py::register(&m)?;
    #[cfg(feature = "alpha-obara-saika")]
    obara_saika_py::register(&m)?;
    #[cfg(feature = "alpha-cga")]
    cga_py::register(&m)?;
    #[cfg(feature = "alpha-gsm")]
    gsm_py::register(&m)?;
    #[cfg(feature = "alpha-sdr")]
    sdr_py::register(&m)?;

    parent.add_submodule(&m)?;
    Ok(())
}
