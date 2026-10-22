//! Python bindings for experimental modules.
//!
//! Each experimental track is gated behind its own feature flag.
//!
//! Every binding either delegates to a real core implementation in
//! [`sci_form_core`] or has been removed — no binding returns hard-coded
//! placeholder data.

use pyo3::prelude::*;

// ─── E5: EEQ ───────────────────────────────────────────────────────────────

#[cfg(feature = "experimental-eeq")]
mod eeq_py {
    // `#[pyfunction]` expands to `Err(err.into())`, which newer clippy
    // versions flag as a useless `PyErr` -> `PyErr` conversion.
    #![allow(clippy::useless_conversion)]

    use pyo3::prelude::*;

    #[pyclass(name = "EeqChargeResult")]
    #[derive(Clone)]
    pub struct EeqChargeResultPy {
        #[pyo3(get)]
        pub charges: Vec<f64>,
        #[pyo3(get)]
        pub coordination_numbers: Vec<f64>,
        #[pyo3(get)]
        pub total_charge: f64,
    }

    #[pyclass(name = "EeqEnergyResult")]
    #[derive(Clone)]
    pub struct EeqEnergyResultPy {
        #[pyo3(get)]
        pub electrostatic_energy: f64,
        #[pyo3(get)]
        pub charges: Vec<f64>,
        #[pyo3(get)]
        pub coordination_numbers: Vec<f64>,
    }

    /// Compute EEQ geometry-dependent charges.
    #[pyfunction]
    #[pyo3(signature = (elements, coords, total_charge=0.0))]
    pub fn compute_eeq_charges(
        elements: Vec<u8>,
        coords: Vec<f64>,
        total_charge: f64,
    ) -> PyResult<EeqChargeResultPy> {
        let pos: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let config = sci_form_core::charges_eeq::EeqConfig {
            total_charge,
            regularization: 1e-10,
        };
        let r = sci_form_core::charges_eeq::compute_eeq_charges(&elements, &pos, &config);
        Ok(EeqChargeResultPy {
            charges: r.charges,
            coordination_numbers: r.coordination_numbers,
            total_charge: r.total_charge,
        })
    }

    /// Compute EEQ electrostatic energy.
    #[pyfunction]
    #[pyo3(signature = (elements, coords, total_charge=0.0))]
    pub fn compute_eeq_energy(
        elements: Vec<u8>,
        coords: Vec<f64>,
        total_charge: f64,
    ) -> PyResult<EeqEnergyResultPy> {
        let pos: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let config = sci_form_core::charges_eeq::EeqConfig {
            total_charge,
            regularization: 1e-10,
        };
        let r = sci_form_core::charges_eeq::compute_eeq_energy(&elements, &pos, &config);
        Ok(EeqEnergyResultPy {
            electrostatic_energy: r.electrostatic_energy,
            charges: r.charges,
            coordination_numbers: r.coordination_numbers,
        })
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(compute_eeq_charges, m)?)?;
        m.add_function(wrap_pyfunction!(compute_eeq_energy, m)?)?;
        Ok(())
    }
}

// ─── E6: ALPB ──────────────────────────────────────────────────────────────

#[cfg(feature = "experimental-alpb")]
mod alpb_py {
    // `#[pyfunction]` expands to `Err(err.into())`, which newer clippy
    // versions flag as a useless `PyErr` -> `PyErr` conversion.
    #![allow(clippy::useless_conversion)]

    use pyo3::prelude::*;

    #[pyclass(name = "AlpbResult")]
    #[derive(Clone)]
    pub struct AlpbResultPy {
        #[pyo3(get)]
        pub electrostatic_energy: f64,
        #[pyo3(get)]
        pub nonpolar_energy: f64,
        #[pyo3(get)]
        pub total_energy: f64,
        #[pyo3(get)]
        pub born_radii: Vec<f64>,
        #[pyo3(get)]
        pub alpb_factor: f64,
    }

    #[pyclass(name = "AlpbBornRadii")]
    #[derive(Clone)]
    pub struct AlpbBornRadiiPy {
        #[pyo3(get)]
        pub radii: Vec<f64>,
        #[pyo3(get)]
        pub intrinsic: Vec<f64>,
    }

    /// Compute ALPB implicit solvation energy.
    #[pyfunction]
    #[pyo3(signature = (elements, coords, charges, solvent_dielectric=78.5, probe_radius=1.4, surface_tension=0.005))]
    pub fn compute_alpb_solvation(
        elements: Vec<u8>,
        coords: Vec<f64>,
        charges: Vec<f64>,
        solvent_dielectric: f64,
        probe_radius: f64,
        surface_tension: f64,
    ) -> PyResult<AlpbResultPy> {
        let pos: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let config = sci_form_core::solvation_alpb::AlpbConfig {
            solvent_dielectric,
            probe_radius,
            surface_tension,
        };
        let r = sci_form_core::solvation_alpb::compute_alpb_solvation(
            &elements, &pos, &charges, &config,
        );
        Ok(AlpbResultPy {
            electrostatic_energy: r.electrostatic_energy,
            nonpolar_energy: r.nonpolar_energy,
            total_energy: r.total_energy,
            born_radii: r.born_radii,
            alpb_factor: r.alpb_factor,
        })
    }

    /// Compute ALPB-style Born radii.
    #[pyfunction]
    #[pyo3(signature = (elements, coords, probe_radius=1.4))]
    pub fn compute_alpb_born_radii(
        elements: Vec<u8>,
        coords: Vec<f64>,
        probe_radius: f64,
    ) -> PyResult<AlpbBornRadiiPy> {
        let pos: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let r = sci_form_core::solvation_alpb::compute_born_radii(&elements, &pos, probe_radius);
        Ok(AlpbBornRadiiPy {
            radii: r.radii,
            intrinsic: r.intrinsic,
        })
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(compute_alpb_solvation, m)?)?;
        m.add_function(wrap_pyfunction!(compute_alpb_born_radii, m)?)?;
        Ok(())
    }
}

// ─── E7: D4 ────────────────────────────────────────────────────────────────

#[cfg(feature = "experimental-d4")]
mod d4_py {
    // `#[pyfunction]` expands to `Err(err.into())`, which newer clippy
    // versions flag as a useless `PyErr` -> `PyErr` conversion.
    #![allow(clippy::useless_conversion)]

    use pyo3::prelude::*;

    #[pyclass(name = "D4Result")]
    #[derive(Clone)]
    pub struct D4ResultPy {
        #[pyo3(get)]
        pub e2_body: f64,
        #[pyo3(get)]
        pub e3_body: f64,
        #[pyo3(get)]
        pub total_energy: f64,
        #[pyo3(get)]
        pub total_kcal_mol: f64,
        #[pyo3(get)]
        pub coordination_numbers: Vec<f64>,
    }

    /// Compute DFT-D4 dispersion energy.
    #[pyfunction]
    #[pyo3(signature = (elements, coords, s6=1.0, s8=0.95, a1=0.45, a2=4.0, three_body=false))]
    pub fn compute_d4_energy(
        elements: Vec<u8>,
        coords: Vec<f64>,
        s6: f64,
        s8: f64,
        a1: f64,
        a2: f64,
        three_body: bool,
    ) -> PyResult<D4ResultPy> {
        let pos: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let config = sci_form_core::dispersion::D4Config {
            s6,
            s8,
            a1,
            a2,
            three_body,
            s9: 1.0,
        };
        let r = sci_form_core::dispersion::compute_d4_energy(&elements, &pos, &config);
        Ok(D4ResultPy {
            e2_body: r.e2_body,
            e3_body: r.e3_body,
            total_energy: r.total_energy,
            total_kcal_mol: r.total_kcal_mol,
            coordination_numbers: r.coordination_numbers,
        })
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(compute_d4_energy, m)?)?;
        Ok(())
    }
}

// ─── E10: CPM ──────────────────────────────────────────────────────────────

// The CPM charges are exposed whenever the core's CPM module is available.
// `experimental-cpm` is a legacy alias of `beta-cpm`, so the two features are
// equivalent; gating on both makes `sci_form.compute_cpm_charges` (and its
// deprecated alias `sci_form.cpm_charges`) available in every build that
// ships the core CPM code (blocker C-C3).
#[cfg(any(feature = "experimental-cpm", feature = "beta-cpm"))]
mod cpm_py {
    // `#[pyfunction]` expands to `Err(err.into())`, which newer clippy
    // versions flag as a useless `PyErr` -> `PyErr` conversion.
    #![allow(clippy::useless_conversion)]

    use pyo3::prelude::*;

    #[pyclass(name = "CpmResult")]
    #[derive(Clone)]
    pub struct CpmResultPy {
        #[pyo3(get)]
        pub charges: Vec<f64>,
        #[pyo3(get)]
        pub total_charge: f64,
        #[pyo3(get)]
        pub grand_potential: f64,
        #[pyo3(get)]
        pub electrostatic_energy: f64,
        #[pyo3(get)]
        pub mu_ev: f64,
        #[pyo3(get)]
        pub iterations: usize,
        #[pyo3(get)]
        pub converged: bool,
    }

    #[pyclass(name = "CpmSurface")]
    #[derive(Clone)]
    pub struct CpmSurfacePy {
        #[pyo3(get)]
        pub mu_values: Vec<f64>,
        #[pyo3(get)]
        pub total_charge: Vec<f64>,
        #[pyo3(get)]
        pub free_energy: Vec<f64>,
        #[pyo3(get)]
        pub capacitance: Vec<f64>,
        #[pyo3(get)]
        pub all_converged: bool,
    }

    /// Compute CPM charges at a given electrochemical potential.
    #[pyfunction]
    #[pyo3(signature = (elements, coords, mu_ev=-4.44, dielectric=78.5))]
    pub fn compute_cpm_charges(
        elements: Vec<u8>,
        coords: Vec<f64>,
        mu_ev: f64,
        dielectric: f64,
    ) -> PyResult<CpmResultPy> {
        let pos: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let config = sci_form_core::beta::cpm::CpmConfig {
            mu_ev,
            dielectric,
            max_iter: 100,
            charge_tol: 1e-6,
            mu_per_atom: None,
        };
        let r = sci_form_core::beta::cpm::compute_cpm_charges(&elements, &pos, &config);
        Ok(CpmResultPy {
            charges: r.charges,
            total_charge: r.total_charge,
            grand_potential: r.grand_potential,
            electrostatic_energy: r.electrostatic_energy,
            mu_ev: r.mu_ev,
            iterations: r.iterations,
            converged: r.converged,
        })
    }

    /// Scan electrochemical surface over a potential range.
    #[pyfunction]
    #[pyo3(signature = (elements, coords, mu_min=-5.5, mu_max=-3.5, n_points=10, dielectric=78.5))]
    pub fn compute_cpm_surface(
        elements: Vec<u8>,
        coords: Vec<f64>,
        mu_min: f64,
        mu_max: f64,
        n_points: usize,
        dielectric: f64,
    ) -> PyResult<CpmSurfacePy> {
        let pos: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let r = sci_form_core::beta::cpm::compute_cpm_surface(
            &elements, &pos, mu_min, mu_max, n_points, dielectric,
        );
        Ok(CpmSurfacePy {
            mu_values: r.mu_values,
            total_charge: r.total_charge,
            free_energy: r.free_energy,
            capacitance: r.capacitance,
            all_converged: r.all_converged,
        })
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<CpmResultPy>()?;
        m.add_class::<CpmSurfacePy>()?;
        m.add_function(wrap_pyfunction!(compute_cpm_charges, m)?)?;
        m.add_function(wrap_pyfunction!(compute_cpm_surface, m)?)?;
        Ok(())
    }
}

/// Register all enabled experimental bindings.
pub(crate) fn register(_m: &Bound<'_, PyModule>) -> PyResult<()> {
    #[cfg(feature = "experimental-eeq")]
    eeq_py::register(_m)?;

    #[cfg(feature = "experimental-alpb")]
    alpb_py::register(_m)?;

    #[cfg(feature = "experimental-d4")]
    d4_py::register(_m)?;

    #[cfg(any(feature = "experimental-cpm", feature = "beta-cpm"))]
    cpm_py::register(_m)?;

    #[cfg(feature = "alpha-edl")]
    alpha_edl_py::register(_m)?;

    #[cfg(feature = "alpha-periodic-linear")]
    alpha_periodic_py::register(_m)?;

    #[cfg(feature = "alpha-kinetics")]
    alpha_kinetics_py::register(_m)?;

    #[cfg(feature = "alpha-render-bridge")]
    alpha_render_py::register(_m)?;

    #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-edl"))]
    alpha_edl_chart_py::register(_m)?;

    #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-kinetics"))]
    alpha_kinetics_chart_py::register(_m)?;

    #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
    alpha_periodic_chart_py::register(_m)?;

    #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-gsm"))]
    alpha_gsm_chart_py::register(_m)?;

    Ok(())
}

// ─── Alpha: EDL ────────────────────────────────────────────────────────────

#[cfg(feature = "alpha-edl")]
mod alpha_edl_py {
    // `#[pyfunction]` expands to `Err(err.into())`, which newer clippy
    // versions flag as a useless `PyErr` -> `PyErr` conversion.
    #![allow(clippy::useless_conversion)]

    use pyo3::prelude::*;

    #[pyclass(name = "EdlProfileResult")]
    #[derive(Clone)]
    pub struct EdlProfileResultPy {
        #[pyo3(get)]
        pub distance_axis_angstrom: Vec<f64>,
        #[pyo3(get)]
        pub electrostatic_potential_v: Vec<f64>,
        #[pyo3(get)]
        pub field_strength_v_per_m: Vec<f64>,
        #[pyo3(get)]
        pub charge_density_c_per_m3: Vec<f64>,
        #[pyo3(get)]
        pub compact_layer_drop_v: f64,
        #[pyo3(get)]
        pub diffuse_layer_drop_v: f64,
        #[pyo3(get)]
        pub total_interfacial_drop_v: f64,
        #[pyo3(get)]
        pub capacitance_total_f_per_m2: f64,
        #[pyo3(get)]
        pub model_name: String,
        #[pyo3(get)]
        pub converged: bool,
    }

    /// Compute an EDL profile for the given model, surface potential, and electrolyte conditions.
    #[pyfunction]
    #[pyo3(signature = (
        surface_potential_v,
        model = "gouy-chapman",
        ionic_strength_m = 0.1,
        temperature_k = 298.15,
        stern_thickness_angstrom = 3.0,
        compact_dielectric = 6.0,
        bulk_dielectric = 78.5,
        n_points = 128,
        extent_angstrom = 12.0
    ))]
    #[allow(clippy::too_many_arguments)]
    pub fn compute_edl_profile(
        surface_potential_v: f64,
        model: &str,
        ionic_strength_m: f64,
        temperature_k: f64,
        stern_thickness_angstrom: f64,
        compact_dielectric: f64,
        bulk_dielectric: f64,
        n_points: usize,
        extent_angstrom: f64,
    ) -> PyResult<EdlProfileResultPy> {
        use sci_form_core::alpha::edl::*;
        let edl_model = match model {
            "helmholtz" => EdlModel::Helmholtz,
            "gouy-chapman" => EdlModel::GouyChapman,
            "gouy-chapman-stern" | "gcs" => EdlModel::GouyChapmanStern,
            _ => {
                return Err(crate::errors::value_error(format!(
                    "Invalid 'model' '{}'. Supported: 'helmholtz', 'gouy-chapman', 'gouy-chapman-stern', 'gcs'.",
                    model
                )))
            }
        };
        let config = EdlConfig {
            model: edl_model,
            temperature_k,
            ionic_strength_m,
            stern_thickness_angstrom,
            compact_dielectric,
            bulk_dielectric,
            numerics: EdlNumerics {
                n_points,
                extent_angstrom,
            },
            ..Default::default()
        };
        let result = compute_edl_profile(surface_potential_v, &config)
            .map_err(crate::errors::engine_error)?;
        Ok(EdlProfileResultPy {
            distance_axis_angstrom: result.distance_axis_angstrom,
            electrostatic_potential_v: result.electrostatic_potential_v,
            field_strength_v_per_m: result.field_strength_v_per_m,
            charge_density_c_per_m3: result.charge_density_c_per_m3,
            compact_layer_drop_v: result.compact_layer_drop_v,
            diffuse_layer_drop_v: result.diffuse_layer_drop_v,
            total_interfacial_drop_v: result.total_interfacial_drop_v,
            capacitance_total_f_per_m2: result.differential_capacitance.total_f_per_m2,
            model_name: result.model_name,
            converged: result.converged,
        })
    }

    /// Scan EDL capacitance over a range of surface potentials.
    #[pyfunction]
    #[pyo3(signature = (v_min, v_max, n_points, ionic_strength_m = 0.1, model = "gcs"))]
    pub fn compute_edl_capacitance_scan(
        v_min: f64,
        v_max: f64,
        n_points: usize,
        ionic_strength_m: f64,
        model: &str,
    ) -> PyResult<Vec<(f64, f64)>> {
        use sci_form_core::alpha::edl::*;
        let edl_model = match model {
            "helmholtz" => EdlModel::Helmholtz,
            "gouy-chapman" => EdlModel::GouyChapman,
            "gouy-chapman-stern" | "gcs" => EdlModel::GouyChapmanStern,
            _ => {
                return Err(crate::errors::value_error(format!(
                    "Invalid 'model' '{}'. Supported: 'helmholtz', 'gouy-chapman', 'gouy-chapman-stern', 'gcs'.",
                    model
                )))
            }
        };
        let config = EdlConfig {
            model: edl_model,
            ionic_strength_m,
            ..Default::default()
        };
        scan_edl_capacitance(v_min, v_max, n_points, &config).map_err(crate::errors::engine_error)
    }

    /// Compute Helmholtz-only profile (linearcompact compact layer only, no diffuse response).
    #[pyfunction]
    #[pyo3(signature = (surface_potential_v, ionic_strength_m = 0.1, temperature_k = 298.15))]
    pub fn compute_helmholtz_profile(
        surface_potential_v: f64,
        ionic_strength_m: f64,
        temperature_k: f64,
    ) -> PyResult<EdlProfileResultPy> {
        use sci_form_core::alpha::edl::*;
        let config = EdlConfig {
            model: EdlModel::Helmholtz,
            ionic_strength_m,
            temperature_k,
            ..Default::default()
        };
        let result = compute_helmholtz_profile(surface_potential_v, &config)
            .map_err(crate::errors::engine_error)?;
        Ok(EdlProfileResultPy {
            distance_axis_angstrom: result.distance_axis_angstrom,
            electrostatic_potential_v: result.electrostatic_potential_v,
            field_strength_v_per_m: result.field_strength_v_per_m,
            charge_density_c_per_m3: result.charge_density_c_per_m3,
            compact_layer_drop_v: result.compact_layer_drop_v,
            diffuse_layer_drop_v: result.diffuse_layer_drop_v,
            total_interfacial_drop_v: result.total_interfacial_drop_v,
            capacitance_total_f_per_m2: result.differential_capacitance.total_f_per_m2,
            model_name: result.model_name,
            converged: result.converged,
        })
    }

    /// Compute Gouy-Chapman profile (diffuse layer only).
    #[pyfunction]
    #[pyo3(signature = (surface_potential_v, ionic_strength_m = 0.1, temperature_k = 298.15))]
    pub fn compute_gouy_chapman_profile(
        surface_potential_v: f64,
        ionic_strength_m: f64,
        temperature_k: f64,
    ) -> PyResult<EdlProfileResultPy> {
        use sci_form_core::alpha::edl::*;
        let config = EdlConfig {
            model: EdlModel::GouyChapman,
            ionic_strength_m,
            temperature_k,
            ..Default::default()
        };
        let result = compute_gouy_chapman_profile(surface_potential_v, &config)
            .map_err(crate::errors::engine_error)?;
        Ok(EdlProfileResultPy {
            distance_axis_angstrom: result.distance_axis_angstrom,
            electrostatic_potential_v: result.electrostatic_potential_v,
            field_strength_v_per_m: result.field_strength_v_per_m,
            charge_density_c_per_m3: result.charge_density_c_per_m3,
            compact_layer_drop_v: result.compact_layer_drop_v,
            diffuse_layer_drop_v: result.diffuse_layer_drop_v,
            total_interfacial_drop_v: result.total_interfacial_drop_v,
            capacitance_total_f_per_m2: result.differential_capacitance.total_f_per_m2,
            model_name: result.model_name,
            converged: result.converged,
        })
    }

    /// Compute Gouy-Chapman-Stern profile (compact + diffuse layers).
    #[pyfunction]
    #[pyo3(signature = (surface_potential_v, ionic_strength_m = 0.1, temperature_k = 298.15, stern_thickness_a = 3.0))]
    pub fn compute_gcs_profile(
        surface_potential_v: f64,
        ionic_strength_m: f64,
        temperature_k: f64,
        stern_thickness_a: f64,
    ) -> PyResult<EdlProfileResultPy> {
        use sci_form_core::alpha::edl::*;
        let config = EdlConfig {
            model: EdlModel::GouyChapmanStern,
            ionic_strength_m,
            temperature_k,
            stern_thickness_angstrom: stern_thickness_a,
            ..Default::default()
        };
        let result = compute_gcs_profile(surface_potential_v, &config)
            .map_err(crate::errors::engine_error)?;
        Ok(EdlProfileResultPy {
            distance_axis_angstrom: result.distance_axis_angstrom,
            electrostatic_potential_v: result.electrostatic_potential_v,
            field_strength_v_per_m: result.field_strength_v_per_m,
            charge_density_c_per_m3: result.charge_density_c_per_m3,
            compact_layer_drop_v: result.compact_layer_drop_v,
            diffuse_layer_drop_v: result.diffuse_layer_drop_v,
            total_interfacial_drop_v: result.total_interfacial_drop_v,
            capacitance_total_f_per_m2: result.differential_capacitance.total_f_per_m2,
            model_name: result.model_name,
            converged: result.converged,
        })
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(compute_edl_profile, m)?)?;
        m.add_function(wrap_pyfunction!(compute_edl_capacitance_scan, m)?)?;
        m.add_function(wrap_pyfunction!(compute_helmholtz_profile, m)?)?;
        m.add_function(wrap_pyfunction!(compute_gouy_chapman_profile, m)?)?;
        m.add_function(wrap_pyfunction!(compute_gcs_profile, m)?)?;
        Ok(())
    }
}

// ─── Shared: real-space operator blocks ────────────────────────────────────

/// Parsing helpers for the real-space operator blocks consumed by the core
/// k-mesh API (`build_bloch_hamiltonian`, `assemble_periodic_operators`,
/// `compute_periodic_kpm_dos_from_operators`, `solve_periodic_randnla`).
///
/// A block is `{"matrix": [[...], ...], "translation": [i, j, k]}` where the
/// matrix is square and the translation is expressed in lattice-vector units.
#[cfg(feature = "alpha-periodic-linear")]
mod periodic_blocks {
    use nalgebra::DMatrix;

    /// Owned real-space operator blocks: matrix plus lattice translation.
    pub type RealSpaceBlocks = Vec<(DMatrix<f64>, [i32; 3])>;

    /// Borrowed real-space operator blocks in the shape expected by the core API.
    pub type RealSpaceBlockRefs<'a> = Vec<(&'a DMatrix<f64>, [i32; 3])>;

    /// Parse a JSON array of real-space blocks.
    pub fn parse_blocks(json_str: &str) -> Result<RealSpaceBlocks, String> {
        let value: serde_json::Value = serde_json::from_str(json_str)
            .map_err(|e| format!("invalid real-space block JSON: {e}"))?;
        let array = value
            .as_array()
            .ok_or_else(|| "real-space blocks must be a JSON array".to_string())?;
        let mut blocks = Vec::with_capacity(array.len());
        for item in array {
            let matrix = item
                .get("matrix")
                .and_then(|m| m.as_array())
                .ok_or_else(|| "each real-space block needs a 'matrix' array".to_string())?;
            let translation = item
                .get("translation")
                .and_then(|t| t.as_array())
                .ok_or_else(|| "each real-space block needs a 'translation' array".to_string())?;
            if translation.len() != 3 {
                return Err("lattice translations need exactly three components".into());
            }
            let mut trans = [0i32; 3];
            for (slot, value) in trans.iter_mut().zip(translation.iter()) {
                *slot = value
                    .as_i64()
                    .ok_or_else(|| "lattice translations must be integers".to_string())?
                    as i32;
            }
            let n = matrix.len();
            if n == 0 {
                return Err("real-space matrices must be non-empty".into());
            }
            let mut flat = Vec::with_capacity(n * n);
            for row in matrix {
                let row = row
                    .as_array()
                    .ok_or_else(|| "real-space matrix rows must be arrays".to_string())?;
                if row.len() != n {
                    return Err("real-space matrices must be square".into());
                }
                for entry in row {
                    flat.push(
                        entry.as_f64().ok_or_else(|| {
                            "real-space matrix entries must be numbers".to_string()
                        })?,
                    );
                }
            }
            blocks.push((DMatrix::from_row_slice(n, n, &flat), trans));
        }
        if blocks.is_empty() {
            return Err("at least one real-space block is required".into());
        }
        Ok(blocks)
    }

    /// Borrow parsed blocks in the shape expected by the core k-mesh API.
    pub fn as_core_refs(blocks: &RealSpaceBlocks) -> RealSpaceBlockRefs<'_> {
        blocks
            .iter()
            .map(|(matrix, trans)| (matrix, *trans))
            .collect()
    }
}

// ─── Alpha: Periodic Linear ───────────────────────────────────────────────

#[cfg(feature = "alpha-periodic-linear")]
mod alpha_periodic_py {
    // `#[pyfunction]` expands to `Err(err.into())`, which newer clippy
    // versions flag as a useless `PyErr` -> `PyErr` conversion.
    #![allow(clippy::useless_conversion)]

    use pyo3::prelude::*;

    use super::kmesh_support::{matrix_rows, parse_centering, symmetric_eigenvalues};
    use super::periodic_blocks::{as_core_refs, parse_blocks};

    #[pyclass(name = "KMeshResult")]
    #[derive(Clone)]
    pub struct KMeshResultPy {
        #[pyo3(get)]
        pub n_points: usize,
        #[pyo3(get)]
        pub grid: [usize; 3],
        #[pyo3(get)]
        pub fractional_coords: Vec<Vec<f64>>,
        #[pyo3(get)]
        pub weights: Vec<f64>,
    }

    /// Generate a Monkhorst-Pack or Gamma-centered k-mesh.
    #[pyfunction]
    #[pyo3(signature = (grid, centering = "monkhorst-pack"))]
    pub fn compute_kmesh(grid: [usize; 3], centering: &str) -> PyResult<KMeshResultPy> {
        use sci_form_core::alpha::periodic_linear::*;
        let c = match centering {
            "monkhorst-pack" | "mp" => KMeshCentering::MonkhorstPack,
            "gamma" | "gamma-centered" => KMeshCentering::GammaCentered,
            _ => {
                return Err(crate::errors::value_error(format!(
                    "Invalid 'centering' '{}'. Supported: 'monkhorst-pack', 'mp', 'gamma', 'gamma-centered'.",
                    centering
                )))
            }
        };
        let mesh = monkhorst_pack_mesh(&KMeshConfig { grid, centering: c })
            .map_err(crate::errors::engine_error)?;
        Ok(KMeshResultPy {
            n_points: mesh.points.len(),
            grid,
            fractional_coords: mesh.points.iter().map(|p| p.fractional.to_vec()).collect(),
            weights: mesh.points.iter().map(|p| p.weight).collect(),
        })
    }

    /// Compute the k-averaged density of states with the core KPM solver.
    ///
    /// `blocks_json` is a JSON array of real-space operator blocks; the Bloch
    /// Hamiltonians are assembled for every k-point of the requested mesh and
    /// the density of states is evaluated by
    /// `sci_form_core::alpha::periodic_linear::compute_periodic_kpm_dos_from_operators`.
    #[pyfunction]
    #[pyo3(signature = (
        blocks_json,
        grid = [2, 2, 2],
        centering = "monkhorst-pack",
        order = 100,
        e_min = -20.0,
        e_max = 10.0,
        n_points = 512,
        fermi_energy_ev = None
    ))]
    #[allow(clippy::too_many_arguments)]
    pub fn compute_periodic_dos(
        blocks_json: &str,
        grid: [usize; 3],
        centering: &str,
        order: usize,
        e_min: f64,
        e_max: f64,
        n_points: usize,
        fermi_energy_ev: Option<f64>,
    ) -> PyResult<String> {
        use sci_form_core::alpha::periodic_linear::*;
        let blocks = parse_blocks(blocks_json).map_err(crate::errors::value_error)?;
        let centering = parse_centering(centering).map_err(crate::errors::value_error)?;
        let kmesh_config = KMeshConfig { grid, centering };
        let mesh = monkhorst_pack_mesh(&kmesh_config).map_err(crate::errors::engine_error)?;
        let hamiltonian_refs = as_core_refs(&blocks);
        let (hamiltonians, _overlaps) =
            sci_form_core::kmesh::assemble_periodic_operators(&hamiltonian_refs, &[], &mesh);
        let config = PeriodicKpmAdapterConfig {
            kmesh: kmesh_config,
            order,
            e_min_ev: e_min,
            e_max_ev: e_max,
            n_points,
            fermi_energy_ev,
            ..Default::default()
        };
        let result = compute_periodic_kpm_dos_from_operators(&hamiltonians, &config)
            .map_err(crate::errors::engine_error)?;
        Ok(serde_json::json!({
            "energies_ev": result.energies_ev,
            "total_dos": result.total_dos,
            "n_kpoints": result.diagnostics.n_kpoints,
            "order": order,
            "e_min": e_min,
            "e_max": e_max,
            "fermi_energy_ev": fermi_energy_ev,
            "band_edges": result.band_edges,
            "diagnostics": result.diagnostics,
        })
        .to_string())
    }

    /// Solve the periodic generalized eigenproblem with the core RandNLA solver.
    ///
    /// `hamiltonian_blocks_json` and `overlap_blocks_json` are JSON arrays of
    /// real-space operator blocks; the k-point Hamiltonians/overlaps are
    /// assembled by `sci_form_core::kmesh::assemble_periodic_operators` and
    /// solved by `sci_form_core::alpha::periodic_linear::solve_periodic_randnla`.
    #[pyfunction]
    #[pyo3(signature = (
        hamiltonian_blocks_json,
        overlap_blocks_json,
        n_electrons = 2,
        grid = [2, 2, 2],
        centering = "monkhorst-pack",
        sketch_size = None
    ))]
    #[allow(clippy::too_many_arguments)]
    pub fn solve_periodic_randnla(
        hamiltonian_blocks_json: &str,
        overlap_blocks_json: &str,
        n_electrons: usize,
        grid: [usize; 3],
        centering: &str,
        sketch_size: Option<usize>,
    ) -> PyResult<String> {
        use sci_form_core::alpha::periodic_linear::*;
        let hamiltonian_blocks =
            parse_blocks(hamiltonian_blocks_json).map_err(crate::errors::value_error)?;
        let overlap_blocks =
            parse_blocks(overlap_blocks_json).map_err(crate::errors::value_error)?;
        let centering = parse_centering(centering).map_err(crate::errors::value_error)?;
        let kmesh_config = KMeshConfig { grid, centering };
        let mesh = monkhorst_pack_mesh(&kmesh_config).map_err(crate::errors::engine_error)?;
        let hamiltonian_refs = as_core_refs(&hamiltonian_blocks);
        let overlap_refs = as_core_refs(&overlap_blocks);
        let (hamiltonians, overlaps) = sci_form_core::kmesh::assemble_periodic_operators(
            &hamiltonian_refs,
            &overlap_refs,
            &mesh,
        );
        let config = PeriodicRandNlaAdapterConfig {
            kmesh: kmesh_config,
            sketch_size,
            ..Default::default()
        };
        let result = sci_form_core::alpha::periodic_linear::solve_periodic_randnla(
            &hamiltonians,
            &overlaps,
            n_electrons,
            &config,
        )
        .map_err(crate::errors::engine_error)?;
        Ok(serde_json::json!({
            "n_kpoints": result.diagnostics.n_kpoints,
            "sketch_size": sketch_size,
            "homo_energy_ev": result.band_edges.homo_energy_ev,
            "lumo_energy_ev": result.band_edges.lumo_energy_ev,
            "band_gap_ev": result.band_edges.direct_gap_ev,
            "direct_gap_ev": result.band_edges.direct_gap_ev,
            "indirect_gap_ev": result.band_edges.indirect_gap_ev,
            "max_residual": result.diagnostics.max_residual,
            "used_fallback": result.diagnostics.used_fallback,
            "kpoint_results": result.kpoint_results,
        })
        .to_string())
    }

    /// Compute the Bloch phase factor for a lattice translation.
    #[pyfunction]
    pub fn bloch_phase(
        k_x: f64,
        k_y: f64,
        k_z: f64,
        t_x: i32,
        t_y: i32,
        t_z: i32,
    ) -> PyResult<(f64, f64)> {
        let (cos_theta, sin_theta) =
            sci_form_core::kmesh::bloch_phase([k_x, k_y, k_z], [t_x, t_y, t_z]);
        Ok((cos_theta, sin_theta))
    }

    /// Build the Bloch Hamiltonian `H(k)` from real-space operator blocks.
    #[pyfunction]
    #[pyo3(signature = (blocks_json, k_x = 0.0, k_y = 0.0, k_z = 0.0))]
    pub fn build_bloch_hamiltonian(
        blocks_json: &str,
        k_x: f64,
        k_y: f64,
        k_z: f64,
    ) -> PyResult<String> {
        let blocks = parse_blocks(blocks_json).map_err(crate::errors::value_error)?;
        let refs = as_core_refs(&blocks);
        let hamiltonian = sci_form_core::kmesh::build_bloch_hamiltonian(&refs, [k_x, k_y, k_z]);
        let n = hamiltonian.nrows();
        let eigenvalues = symmetric_eigenvalues(&hamiltonian);
        Ok(serde_json::json!({
            "n_basis": n,
            "k": [k_x, k_y, k_z],
            "matrix_shape": [n, n],
            "matrix": matrix_rows(&hamiltonian),
            "eigenvalues": eigenvalues,
        })
        .to_string())
    }

    /// Assemble periodic Hamiltonians and overlaps for every k-point of a mesh.
    #[pyfunction]
    #[pyo3(signature = (
        hamiltonian_blocks_json,
        overlap_blocks_json = "[]",
        grid = [2, 2, 2],
        centering = "monkhorst-pack"
    ))]
    pub fn assemble_periodic_operators(
        hamiltonian_blocks_json: &str,
        overlap_blocks_json: &str,
        grid: [usize; 3],
        centering: &str,
    ) -> PyResult<String> {
        let hamiltonian_blocks =
            parse_blocks(hamiltonian_blocks_json).map_err(crate::errors::value_error)?;
        let overlap_blocks = if overlap_blocks_json.trim().is_empty() {
            Vec::new()
        } else {
            parse_blocks(overlap_blocks_json).map_err(crate::errors::value_error)?
        };
        let centering = parse_centering(centering).map_err(crate::errors::value_error)?;
        let mesh = sci_form_core::kmesh::monkhorst_pack_mesh(&sci_form_core::kmesh::KMeshConfig {
            grid,
            centering,
        })
        .map_err(crate::errors::engine_error)?;
        let hamiltonian_refs = as_core_refs(&hamiltonian_blocks);
        let overlap_refs = as_core_refs(&overlap_blocks);
        let (hamiltonians, overlaps) = sci_form_core::kmesh::assemble_periodic_operators(
            &hamiltonian_refs,
            &overlap_refs,
            &mesh,
        );
        let n_basis = hamiltonians.first().map(|m| m.nrows()).unwrap_or_default();
        Ok(serde_json::json!({
            "n_kpoints": hamiltonians.len(),
            "n_basis": n_basis,
            "assembled": hamiltonians.len() == mesh.points.len(),
            "kpoints": mesh.points.iter().map(|p| serde_json::json!({
                "fractional": p.fractional,
                "index": p.index,
                "weight": p.weight,
            })).collect::<Vec<_>>(),
            "hamiltonians": hamiltonians.iter().map(matrix_rows).collect::<Vec<_>>(),
            "overlaps": overlaps.iter().map(matrix_rows).collect::<Vec<_>>(),
        })
        .to_string())
    }

    /// Integrate a scalar (density, energy) over the Brillouin zone.
    #[pyfunction]
    pub fn bz_integrate_scalar(values: Vec<f64>, weights: Vec<f64>) -> PyResult<f64> {
        if values.len() != weights.len() {
            return Err(crate::errors::value_error(
                "values and weights must have same length",
            ));
        }
        Ok(sci_form_core::kmesh::bz_integrate_scalar(&values, &weights))
    }

    /// Integrate a vector (current, magnetization) over the Brillouin zone.
    #[pyfunction]
    pub fn bz_integrate_vector(vectors: Vec<Vec<f64>>, weights: Vec<f64>) -> PyResult<Vec<f64>> {
        if vectors.len() != weights.len() {
            return Err(crate::errors::value_error(
                "vectors and weights must have same length",
            ));
        }
        Ok(sci_form_core::kmesh::bz_integrate_vector(
            &vectors, &weights,
        ))
    }

    /// Check that the electron count fits into the available band/k-point states.
    #[pyfunction]
    pub fn validate_electron_count(
        n_electrons: usize,
        n_bands: usize,
        n_kpoints: usize,
    ) -> PyResult<bool> {
        Ok(n_electrons <= n_bands * n_kpoints)
    }

    /// Verify that k-mesh weights sum to 1 within a tolerance.
    #[pyfunction]
    pub fn validate_kmesh_weights(weights: Vec<f64>, tolerance: f64) -> PyResult<bool> {
        use sci_form_core::kmesh::{KMesh, KMeshCentering, KPoint};
        let points: Vec<KPoint> = weights
            .iter()
            .map(|&weight| KPoint {
                fractional: [0.0, 0.0, 0.0],
                index: [0, 0, 0],
                weight,
            })
            .collect();
        let mesh = KMesh {
            total_weight: weights.iter().sum(),
            points,
            grid: [1, 1, 1],
            centering: KMeshCentering::MonkhorstPack,
        };
        Ok(sci_form_core::kmesh::validate_kmesh_weights(
            &mesh, tolerance,
        ))
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(compute_kmesh, m)?)?;
        m.add_function(wrap_pyfunction!(compute_periodic_dos, m)?)?;
        m.add_function(wrap_pyfunction!(solve_periodic_randnla, m)?)?;
        m.add_function(wrap_pyfunction!(bloch_phase, m)?)?;
        m.add_function(wrap_pyfunction!(build_bloch_hamiltonian, m)?)?;
        m.add_function(wrap_pyfunction!(assemble_periodic_operators, m)?)?;
        m.add_function(wrap_pyfunction!(bz_integrate_scalar, m)?)?;
        m.add_function(wrap_pyfunction!(bz_integrate_vector, m)?)?;
        m.add_function(wrap_pyfunction!(validate_electron_count, m)?)?;
        m.add_function(wrap_pyfunction!(validate_kmesh_weights, m)?)?;
        Ok(())
    }
}

// ─── Alpha: Kinetics ──────────────────────────────────────────────────────

#[cfg(feature = "alpha-kinetics")]
mod alpha_kinetics_py {
    // `#[pyfunction]` expands to `Err(err.into())`, which newer clippy
    // versions flag as a useless `PyErr` -> `PyErr` conversion.
    #![allow(clippy::useless_conversion)]

    use pyo3::prelude::*;

    use super::kinetics_support::network_config;

    #[pyclass(name = "HtstRateResult")]
    #[derive(Clone)]
    pub struct HtstRateResultPy {
        #[pyo3(get)]
        pub step_id: String,
        #[pyo3(get)]
        pub forward_rate_s_inv: f64,
        #[pyo3(get)]
        pub reverse_rate_s_inv: f64,
        #[pyo3(get)]
        pub equilibrium_constant: f64,
    }

    /// Evaluate an HTST transition rate at a given temperature.
    #[pyfunction]
    #[pyo3(signature = (
        activation_free_energy_ev,
        reaction_free_energy_ev,
        temperature_k = 298.15,
        prefactor_s_inv = None,
        step_id = "step"
    ))]
    pub fn compute_htst_rate(
        activation_free_energy_ev: f64,
        reaction_free_energy_ev: f64,
        temperature_k: f64,
        prefactor_s_inv: Option<f64>,
        step_id: &str,
    ) -> PyResult<HtstRateResultPy> {
        use sci_form_core::alpha::kinetics::*;
        let step = ElementaryStep {
            step_id: step_id.to_string(),
            activation_free_energy_ev,
            reaction_free_energy_ev,
            prefactor_s_inv,
        };
        let state = ThermodynamicState {
            temperature_k,
            pressure_bar: 1.0,
        };
        let result = evaluate_htst_rate(&step, state).map_err(crate::errors::engine_error)?;
        Ok(HtstRateResultPy {
            step_id: result.step_id,
            forward_rate_s_inv: result.forward_rate_s_inv,
            reverse_rate_s_inv: result.reverse_rate_s_inv,
            equilibrium_constant: result.equilibrium_constant,
        })
    }

    /// Evaluate HTST rates over a temperature sweep.
    #[pyfunction]
    #[pyo3(signature = (
        activation_free_energy_ev,
        reaction_free_energy_ev,
        temperatures_k,
        prefactor_s_inv = None,
        step_id = "step"
    ))]
    pub fn compute_htst_sweep(
        activation_free_energy_ev: f64,
        reaction_free_energy_ev: f64,
        temperatures_k: Vec<f64>,
        prefactor_s_inv: Option<f64>,
        step_id: &str,
    ) -> PyResult<Vec<HtstRateResultPy>> {
        use sci_form_core::alpha::kinetics::*;
        let step = ElementaryStep {
            step_id: step_id.to_string(),
            activation_free_energy_ev,
            reaction_free_energy_ev,
            prefactor_s_inv,
        };
        let results = evaluate_htst_temperature_sweep(&step, &temperatures_k, 1.0)
            .map_err(crate::errors::engine_error)?;
        Ok(results
            .into_iter()
            .map(|r| HtstRateResultPy {
                step_id: r.step_id,
                forward_rate_s_inv: r.forward_rate_s_inv,
                reverse_rate_s_inv: r.reverse_rate_s_inv,
                equilibrium_constant: r.equilibrium_constant,
            })
            .collect())
    }

    /// Extract convergence diagnostics from a solved microkinetic trace.
    ///
    /// `network_json` is an optional `MicrokineticNetworkConfig` JSON document;
    /// when omitted the built-in A ⇌ B demo network is solved.
    #[pyfunction]
    #[pyo3(signature = (network_json = None, n_steps = 20))]
    pub fn extract_kinetics_diagnostics(
        network_json: Option<&str>,
        n_steps: usize,
    ) -> PyResult<String> {
        use sci_form_core::alpha::kinetics::*;
        let config = network_config(network_json, n_steps).map_err(crate::errors::value_error)?;
        let trace = solve_microkinetic_network(&config).map_err(crate::errors::engine_error)?;
        let diagnostics = extract_kinetics_diagnostics(&trace);
        Ok(serde_json::json!({
            "max_population": diagnostics.max_population,
            "total_population": diagnostics.total_population,
            "mass_conservation_error": diagnostics.mass_conservation_error,
            "max_rate_of_change": diagnostics.max_rate_of_change,
            "is_steady_state": diagnostics.is_steady_state,
        })
        .to_string())
    }

    /// Solve a microkinetic network and report the resulting trace.
    #[pyfunction]
    #[pyo3(signature = (network_json = None, n_steps = 20))]
    pub fn solve_microkinetic_network(
        network_json: Option<&str>,
        n_steps: usize,
    ) -> PyResult<String> {
        use sci_form_core::alpha::kinetics::*;
        let config = network_config(network_json, n_steps).map_err(crate::errors::value_error)?;
        let trace = solve_microkinetic_network(&config).map_err(crate::errors::engine_error)?;
        let diagnostics = extract_kinetics_diagnostics(&trace);
        let last = trace.frames.last();
        Ok(serde_json::json!({
            "n_steps": trace.frames.len(),
            "n_rates": trace.rates.len(),
            "converged": diagnostics.is_steady_state,
            "final_time_s": last.map(|f| f.time_s),
            "final_populations": last
                .map(|f| f.species.iter().map(|s| s.population).collect::<Vec<_>>())
                .unwrap_or_default(),
            "diagnostics": diagnostics,
        })
        .to_string())
    }

    /// Solve a microkinetic network and propagate it to steady state.
    #[pyfunction]
    #[pyo3(signature = (network_json = None, n_steps = 20, tol = 1e-10))]
    pub fn solve_microkinetic_steady_state(
        network_json: Option<&str>,
        n_steps: usize,
        tol: f64,
    ) -> PyResult<String> {
        use sci_form_core::alpha::kinetics::*;
        let config = network_config(network_json, n_steps).map_err(crate::errors::value_error)?;
        let trace =
            solve_microkinetic_steady_state(&config, tol).map_err(crate::errors::engine_error)?;
        let diagnostics = extract_kinetics_diagnostics(&trace);
        let last = trace.frames.last();
        Ok(serde_json::json!({
            "n_steps": trace.frames.len(),
            "converged": diagnostics.is_steady_state,
            "residual": diagnostics.max_rate_of_change,
            "tol": tol,
            "final_time_s": last.map(|f| f.time_s),
            "final_populations": last
                .map(|f| f.species.iter().map(|s| s.population).collect::<Vec<_>>())
                .unwrap_or_default(),
            "diagnostics": diagnostics,
        })
        .to_string())
    }

    /// Analyze one GSM + MBH + HTST coupled elementary step.
    #[pyfunction]
    #[pyo3(signature = (
        smiles,
        reactant_coords,
        product_coords,
        method = "uff",
        step_id = "elementary-step",
        temperature_k = 298.15,
        pressure_bar = 1.0
    ))]
    #[allow(clippy::too_many_arguments)]
    pub fn analyze_gsm_mbh_htst_step(
        smiles: &str,
        reactant_coords: Vec<f64>,
        product_coords: Vec<f64>,
        method: &str,
        step_id: &str,
        temperature_k: f64,
        pressure_bar: f64,
    ) -> PyResult<String> {
        use sci_form_core::alpha::gsm::GsmEnergyBackend;
        use sci_form_core::alpha::kinetics::*;
        let backend: GsmEnergyBackend = method.parse().map_err(|_e: String| {
            crate::errors::value_error(crate::errors::invalid_value(
                "'method'",
                method,
                &[
                    "uff", "mmff94", "pm3", "xtb", "gfn1", "gfn2", "hf3c", "reaxff",
                ],
            ))
        })?;
        let config = HtstAdapterConfig {
            step_id: step_id.to_string(),
            state: ThermodynamicState {
                temperature_k,
                pressure_bar,
            },
            ..Default::default()
        };
        let result = analyze_gsm_mbh_htst_step_with_backend(
            smiles,
            &reactant_coords,
            &product_coords,
            backend,
            &config,
        )
        .map_err(|e| crate::errors::smiles_or_engine_error(e, smiles))?;
        Ok(serde_json::json!({
            "method": backend.as_str(),
            "step_id": result.rate.step_id,
            "converged": result.rate.forward_rate_s_inv.is_finite(),
            "barrier_kcal_mol": result.activation_energy_kcal_mol,
            "reverse_barrier_kcal_mol": result.reverse_barrier_kcal_mol,
            "ts_energy_kcal_mol": result.ts_energy_kcal_mol,
            "forward_rate_s_inv": result.rate.forward_rate_s_inv,
            "reverse_rate_s_inv": result.rate.reverse_rate_s_inv,
            "equilibrium_constant": result.rate.equilibrium_constant,
            "imaginary_frequency_cm1": result.imaginary_frequency_cm1,
            "n_blocks": result.n_blocks,
            "n_flexible": result.n_flexible,
            "n_nodes": result.n_nodes,
            "energy_evaluations": result.energy_evaluations,
        })
        .to_string())
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(compute_htst_rate, m)?)?;
        m.add_function(wrap_pyfunction!(compute_htst_sweep, m)?)?;
        m.add_function(wrap_pyfunction!(extract_kinetics_diagnostics, m)?)?;
        m.add_function(wrap_pyfunction!(solve_microkinetic_network, m)?)?;
        m.add_function(wrap_pyfunction!(solve_microkinetic_steady_state, m)?)?;
        m.add_function(wrap_pyfunction!(analyze_gsm_mbh_htst_step, m)?)?;
        Ok(())
    }
}

// ─── Shared: microkinetic network configuration ────────────────────────────

/// Shared helper that turns an optional network JSON document into a
/// [`sci_form_core::alpha::kinetics::MicrokineticNetworkConfig`].
#[cfg(feature = "alpha-kinetics")]
mod kinetics_support {
    use sci_form_core::alpha::kinetics::{
        ElementaryStep, MicrokineticNetworkConfig, ThermodynamicState,
    };

    /// Parse a network configuration, falling back to a built-in A ⇌ B demo.
    pub fn network_config(
        network_json: Option<&str>,
        n_steps: usize,
    ) -> Result<MicrokineticNetworkConfig, String> {
        if let Some(text) = network_json {
            let mut config: MicrokineticNetworkConfig = serde_json::from_str(text)
                .map_err(|e| format!("invalid microkinetic network JSON: {e}"))?;
            if config.n_frames < 2 {
                config.n_frames = n_steps.max(2);
            }
            return Ok(config);
        }
        Ok(MicrokineticNetworkConfig {
            species: vec!["A".into(), "B".into()],
            steps: vec![ElementaryStep {
                step_id: "A->B".into(),
                activation_free_energy_ev: 0.5,
                reaction_free_energy_ev: -0.1,
                prefactor_s_inv: None,
            }],
            stoichiometry: vec![vec![(0, -1.0), (1, 1.0)]],
            initial_populations: vec![1.0, 0.0],
            total_time_s: 1e-3,
            n_frames: n_steps.max(2),
            state: ThermodynamicState {
                temperature_k: 298.15,
                pressure_bar: 1.0,
            },
        })
    }
}

// ─── Alpha: Render Bridge ─────────────────────────────────────────────────

#[cfg(feature = "alpha-render-bridge")]
mod alpha_render_py {
    // `#[pyfunction]` expands to `Err(err.into())`, which newer clippy
    // versions flag as a useless `PyErr` -> `PyErr` conversion.
    #![allow(clippy::useless_conversion)]

    use pyo3::prelude::*;

    #[pyclass(name = "ChartSeries")]
    #[derive(Clone)]
    pub struct ChartSeriesPy {
        #[pyo3(get)]
        pub series_id: String,
        #[pyo3(get)]
        pub label: String,
        #[pyo3(get)]
        pub x: Vec<f64>,
        #[pyo3(get)]
        pub y: Vec<f64>,
        #[pyo3(get)]
        pub x_unit: String,
        #[pyo3(get)]
        pub y_unit: String,
    }

    #[pyclass(name = "TrajectoryFrame")]
    #[derive(Clone)]
    pub struct TrajectoryFramePy {
        #[pyo3(get)]
        pub frame_time: f64,
        #[pyo3(get)]
        pub coords: Vec<f64>,
    }

    #[pyclass(name = "ChartPayload")]
    #[derive(Clone)]
    pub struct ChartPayloadPy {
        #[pyo3(get)]
        pub title: String,
        #[pyo3(get)]
        pub series: Vec<ChartSeriesPy>,
        #[pyo3(get)]
        pub trajectory: Vec<TrajectoryFramePy>,
    }

    /// Serialize a chart payload to the canonical wrapper JSON.
    ///
    /// `chart_json` is a `ChartPayload` JSON document; the result keeps the
    /// `{"format": "json", "version": "1.0", "data": {...}}` envelope with the
    /// real chart in `data`.
    #[pyfunction]
    pub fn chart_to_json(chart_json: &str) -> PyResult<String> {
        let chart = parse_chart(chart_json).map_err(crate::errors::value_error)?;
        Ok(serde_json::json!({
            "format": "json",
            "version": "1.0",
            "data": chart,
        })
        .to_string())
    }

    /// Deserialize a chart payload (wrapper or bare chart) to chart JSON.
    #[pyfunction]
    pub fn chart_from_json(json_str: &str) -> PyResult<String> {
        match parse_chart(json_str) {
            Ok(chart) => Ok(serde_json::to_string(&chart)
                .unwrap_or_else(|e| serde_json::json!({"error": e.to_string()}).to_string())),
            Err(e) => Ok(serde_json::json!({"error": e}).to_string()),
        }
    }

    /// Validate the structure of a chart payload (wrapper or bare chart).
    #[pyfunction]
    pub fn validate_chart_schema(json_str: &str) -> PyResult<bool> {
        Ok(parse_chart(json_str).is_ok())
    }

    /// Pack a chart payload into an Arrow-style record batch.
    ///
    /// Returns the packed column layout plus a lossless round-trip check.
    #[pyfunction]
    pub fn pack_chart_payload(chart_json: &str) -> PyResult<String> {
        use sci_form_core::render_bridge::{pack_chart_payload, verify_chart_round_trip};
        let chart = parse_chart(chart_json).map_err(crate::errors::value_error)?;
        let batch = pack_chart_payload(&chart);
        Ok(serde_json::json!({
            "n_columns": batch.float_columns.len(),
            "column_names": batch.float_columns.iter().map(|c| c.name.clone()).collect::<Vec<_>>(),
            "column_lengths": batch.float_columns.iter().map(|c| c.values.len()).collect::<Vec<_>>(),
            "round_trip_ok": verify_chart_round_trip(&chart, 0.0),
        })
        .to_string())
    }

    /// Parse either the wrapper envelope or a bare chart payload.
    pub fn parse_chart(
        json_str: &str,
    ) -> Result<sci_form_core::render_bridge::ChartPayload, String> {
        let value: serde_json::Value =
            serde_json::from_str(json_str).map_err(|e| format!("invalid JSON: {e}"))?;
        let chart_value = value.get("data").cloned().unwrap_or(value);
        let chart: sci_form_core::render_bridge::ChartPayload = serde_json::from_value(chart_value)
            .map_err(|e| format!("invalid chart payload: {e}"))?;
        if chart.series.iter().any(|s| s.x.len() != s.y.len()) {
            return Err("chart series must have matching x and y lengths".into());
        }
        Ok(chart)
    }

    /// Convert a core chart payload into the Python representation.
    #[cfg(any(
        feature = "alpha-edl",
        feature = "alpha-kinetics",
        feature = "alpha-periodic-linear",
        feature = "alpha-gsm"
    ))]
    pub fn to_py_chart(chart: &sci_form_core::render_bridge::ChartPayload) -> ChartPayloadPy {
        ChartPayloadPy {
            title: chart.title.clone(),
            series: chart
                .series
                .iter()
                .map(|s| ChartSeriesPy {
                    series_id: s.series_id.clone(),
                    label: s.label.clone(),
                    x: s.x.clone(),
                    y: s.y.clone(),
                    x_unit: s.x_unit.clone(),
                    y_unit: s.y_unit.clone(),
                })
                .collect(),
            trajectory: Vec::new(),
        }
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(chart_to_json, m)?)?;
        m.add_function(wrap_pyfunction!(chart_from_json, m)?)?;
        m.add_function(wrap_pyfunction!(validate_chart_schema, m)?)?;
        m.add_function(wrap_pyfunction!(pack_chart_payload, m)?)?;
        Ok(())
    }
}

// ─── Alpha: EDL charts ────────────────────────────────────────────────────

#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-edl"))]
mod alpha_edl_chart_py {
    // `#[pyfunction]` expands to `Err(err.into())`, which newer clippy
    // versions flag as a useless `PyErr` -> `PyErr` conversion.
    #![allow(clippy::useless_conversion)]

    use pyo3::prelude::*;

    use super::alpha_render_py::{to_py_chart, ChartPayloadPy};

    /// Build an EDL profile chart payload from a given surface potential and config.
    #[pyfunction]
    #[pyo3(signature = (surface_potential_v, ionic_strength_m = 0.1))]
    pub fn compute_edl_chart(
        surface_potential_v: f64,
        ionic_strength_m: f64,
    ) -> PyResult<ChartPayloadPy> {
        use sci_form_core::alpha::edl::*;
        use sci_form_core::alpha::render_bridge::edl_profile_chart;
        let config = EdlConfig {
            model: EdlModel::GouyChapman,
            ionic_strength_m,
            ..Default::default()
        };
        let profile = compute_gouy_chapman_profile(surface_potential_v, &config)
            .map_err(crate::errors::engine_error)?;
        Ok(to_py_chart(&edl_profile_chart(&profile)))
    }

    /// Build a capacitance scan chart from min/max potential and config.
    ///
    /// The interfacial charge is obtained by integrating the computed charge
    /// density of every EDL profile in the scan, so the "Q(μ)" series holds
    /// real values instead of placeholders.
    #[pyfunction]
    #[pyo3(signature = (v_min, v_max, n_points, ionic_strength_m = 0.1, model = "gcs"))]
    pub fn compute_capacitance_chart(
        v_min: f64,
        v_max: f64,
        n_points: usize,
        ionic_strength_m: f64,
        model: &str,
    ) -> PyResult<ChartPayloadPy> {
        use sci_form_core::alpha::edl::*;
        use sci_form_core::render_bridge::{ChartPayload, ChartSeries};
        let edl_model = match model {
            "helmholtz" => EdlModel::Helmholtz,
            "gouy-chapman" => EdlModel::GouyChapman,
            "gouy-chapman-stern" | "gcs" => EdlModel::GouyChapmanStern,
            _ => {
                return Err(crate::errors::value_error(format!(
                    "Invalid 'model' '{}'. Supported: 'helmholtz', 'gouy-chapman', 'gouy-chapman-stern', 'gcs'.",
                    model
                )))
            }
        };
        let config = EdlConfig {
            model: edl_model,
            ionic_strength_m,
            ..Default::default()
        };
        let scan = scan_edl_capacitance(v_min, v_max, n_points, &config)
            .map_err(crate::errors::engine_error)?;
        let mut potentials = Vec::with_capacity(scan.len());
        let mut capacitances = Vec::with_capacity(scan.len());
        let mut charges = Vec::with_capacity(scan.len());
        for (potential, capacitance) in &scan {
            let profile =
                compute_edl_profile(*potential, &config).map_err(crate::errors::engine_error)?;
            potentials.push(*potential);
            capacitances.push(*capacitance);
            charges.push(integrated_charge_c_per_m2(&profile));
        }
        let chart = ChartPayload {
            title: "EDL Capacitance Scan".to_string(),
            series: vec![
                ChartSeries {
                    series_id: "charge_vs_potential".into(),
                    label: "Q(μ)".into(),
                    x: potentials.clone(),
                    y: charges,
                    x_unit: "V".into(),
                    y_unit: "C/m^2".into(),
                },
                ChartSeries {
                    series_id: "capacitance".into(),
                    label: "C(μ)".into(),
                    x: potentials,
                    y: capacitances,
                    x_unit: "V".into(),
                    y_unit: "F/m^2".into(),
                },
            ],
        };
        Ok(to_py_chart(&chart))
    }

    /// Integrate the charge density of an EDL profile into a surface charge
    /// density in coulombs per square metre.
    fn integrated_charge_c_per_m2(profile: &sci_form_core::alpha::edl::EdlProfileResult) -> f64 {
        const ANGSTROM_TO_M: f64 = 1.0e-10;
        let mut integral = 0.0_f64;
        for i in 1..profile.distance_axis_angstrom.len() {
            let dx_m = (profile.distance_axis_angstrom[i] - profile.distance_axis_angstrom[i - 1])
                * ANGSTROM_TO_M;
            integral += 0.5
                * (profile.charge_density_c_per_m3[i] + profile.charge_density_c_per_m3[i - 1])
                * dx_m;
        }
        integral
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(compute_edl_chart, m)?)?;
        m.add_function(wrap_pyfunction!(compute_capacitance_chart, m)?)?;
        Ok(())
    }
}

// ─── Alpha: Kinetics charts ───────────────────────────────────────────────

#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-kinetics"))]
mod alpha_kinetics_chart_py {
    // `#[pyfunction]` expands to `Err(err.into())`, which newer clippy
    // versions flag as a useless `PyErr` -> `PyErr` conversion.
    #![allow(clippy::useless_conversion)]

    use pyo3::prelude::*;

    use super::alpha_render_py::{to_py_chart, ChartPayloadPy};
    use super::kinetics_support::network_config;

    /// Build an Arrhenius chart from temperature-dependent rates.
    ///
    /// The chart is the core Arrhenius projection `ln(k)` vs `1000/T`; the
    /// forward rates are supplied by the caller and are used verbatim.
    #[pyfunction]
    #[pyo3(signature = (temperatures_k, rates_s_inv))]
    pub fn compute_arrhenius_chart(
        temperatures_k: Vec<f64>,
        rates_s_inv: Vec<f64>,
    ) -> PyResult<ChartPayloadPy> {
        use sci_form_core::render_bridge::{ChartPayload, ChartSeries};
        if temperatures_k.len() != rates_s_inv.len() {
            return Err(crate::errors::value_error(
                "temperatures and rates must have the same length",
            ));
        }
        let x: Vec<f64> = temperatures_k.iter().map(|t| 1000.0 / t).collect();
        let y: Vec<f64> = rates_s_inv.iter().map(|r| r.ln()).collect();
        let chart = ChartPayload {
            title: "Arrhenius Plot".to_string(),
            series: vec![ChartSeries {
                series_id: "arrhenius".into(),
                label: "ln(k) vs 1000/T".into(),
                x,
                y,
                x_unit: "1000/K".into(),
                y_unit: "ln(s⁻¹)".into(),
            }],
        };
        Ok(to_py_chart(&chart))
    }

    /// Build a thermal-properties chart from an elementary step.
    ///
    /// The HTST rates and equilibrium constants are computed by the core
    /// transition-state evaluator over the requested temperature sweep.
    #[pyfunction]
    #[pyo3(signature = (
        activation_free_energy_ev,
        reaction_free_energy_ev,
        temperatures_k,
        prefactor_s_inv = None,
        step_id = "thermal"
    ))]
    pub fn compute_thermal_prop_chart(
        activation_free_energy_ev: f64,
        reaction_free_energy_ev: f64,
        temperatures_k: Vec<f64>,
        prefactor_s_inv: Option<f64>,
        step_id: &str,
    ) -> PyResult<ChartPayloadPy> {
        use sci_form_core::alpha::kinetics::{evaluate_htst_temperature_sweep, ElementaryStep};
        use sci_form_core::alpha::render_bridge::arrhenius_chart as core_arrhenius_chart;
        use sci_form_core::render_bridge::{ChartPayload, ChartSeries};
        let step = ElementaryStep {
            step_id: step_id.to_string(),
            activation_free_energy_ev,
            reaction_free_energy_ev,
            prefactor_s_inv,
        };
        let results = evaluate_htst_temperature_sweep(&step, &temperatures_k, 1.0)
            .map_err(crate::errors::engine_error)?;
        let base = core_arrhenius_chart(&results);
        let x: Vec<f64> = results
            .iter()
            .map(|r| 1000.0 / r.state.temperature_k)
            .collect();
        let y: Vec<f64> = results
            .iter()
            .map(|r| r.equilibrium_constant.ln())
            .collect();
        let chart = ChartPayload {
            title: "Thermal Properties".to_string(),
            series: vec![
                base.series.into_iter().next().unwrap_or_default(),
                ChartSeries {
                    series_id: "equilibrium_constant".into(),
                    label: "ln(K_eq) vs 1000/T".into(),
                    x,
                    y,
                    x_unit: "1000/K".into(),
                    y_unit: "ln".into(),
                },
            ],
        };
        Ok(to_py_chart(&chart))
    }

    /// Build a phase portrait chart from a solved microkinetic network.
    ///
    /// The two axes are the populations of the first two species over the
    /// simulated trajectory.
    #[pyfunction]
    #[pyo3(signature = (network_json = None, n_steps = 50))]
    pub fn compute_phase_portrait_chart(
        network_json: Option<&str>,
        n_steps: usize,
    ) -> PyResult<ChartPayloadPy> {
        use sci_form_core::alpha::kinetics::solve_microkinetic_network;
        use sci_form_core::render_bridge::{ChartPayload, ChartSeries};
        let config = network_config(network_json, n_steps).map_err(crate::errors::value_error)?;
        let trace = solve_microkinetic_network(&config).map_err(crate::errors::engine_error)?;
        if config.species.len() < 2 {
            return Err(crate::errors::value_error(
                "a phase portrait needs a network with at least two species",
            ));
        }
        let mut x = Vec::with_capacity(trace.frames.len());
        let mut y = Vec::with_capacity(trace.frames.len());
        for frame in &trace.frames {
            let first = frame
                .species
                .iter()
                .find(|s| s.species_id == config.species[0])
                .map(|s| s.population)
                .unwrap_or_default();
            let second = frame
                .species
                .iter()
                .find(|s| s.species_id == config.species[1])
                .map(|s| s.population)
                .unwrap_or_default();
            x.push(first);
            y.push(second);
        }
        let chart = ChartPayload {
            title: "Phase Portrait".to_string(),
            series: vec![ChartSeries {
                series_id: "phase_portrait".into(),
                label: format!("[{}] vs [{}]", config.species[0], config.species[1]),
                x,
                y,
                x_unit: "arb.".into(),
                y_unit: "arb.".into(),
            }],
        };
        Ok(to_py_chart(&chart))
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(compute_arrhenius_chart, m)?)?;
        m.add_function(wrap_pyfunction!(compute_thermal_prop_chart, m)?)?;
        m.add_function(wrap_pyfunction!(compute_phase_portrait_chart, m)?)?;
        Ok(())
    }
}

// ─── Alpha: Periodic charts ───────────────────────────────────────────────

#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
mod alpha_periodic_chart_py {
    // `#[pyfunction]` expands to `Err(err.into())`, which newer clippy
    // versions flag as a useless `PyErr` -> `PyErr` conversion.
    #![allow(clippy::useless_conversion)]

    use pyo3::prelude::*;

    use super::alpha_render_py::{to_py_chart, ChartPayloadPy};
    use super::chart_support::{
        fermi_contour, interpolate_path, parse_lattice, positions_from_flat,
    };
    use super::periodic_blocks::{as_core_refs, parse_blocks};

    /// Build a band-structure chart from a periodic EHT calculation.
    ///
    /// `elements` are atomic numbers, `coords` the flat Cartesian coordinates
    /// of the unit cell and `lattice` the nine lattice-vector components.
    /// The band energies come from the core periodic band-structure adapter.
    #[pyfunction]
    #[pyo3(signature = (elements, coords, lattice, n_kpoints_per_segment = 20))]
    pub fn compute_band_structure_chart(
        elements: Vec<u8>,
        coords: Vec<f64>,
        lattice: Vec<f64>,
        n_kpoints_per_segment: usize,
    ) -> PyResult<ChartPayloadPy> {
        use sci_form_core::alpha::periodic_linear::*;
        use sci_form_core::alpha::render_bridge::band_structure_chart as build_chart;
        if !coords.len().is_multiple_of(3) {
            return Err(crate::errors::value_error(
                "coords length must be a multiple of three",
            ));
        }
        let positions = positions_from_flat(&coords);
        let lattice = parse_lattice(&lattice).map_err(crate::errors::value_error)?;
        let n_electrons = sci_form_core::eht::solve_eht(&elements, &positions, None)
            .map_err(crate::errors::engine_error)?
            .n_electrons;
        let config = BandStructureAdapterConfig {
            n_kpoints_per_segment,
            path: BandStructureAdapterConfig::default().path,
        };
        let bands =
            compute_band_structure_adapter(&elements, &positions, &lattice, n_electrons, &config)
                .map_err(crate::errors::engine_error)?;
        Ok(to_py_chart(&build_chart(&bands)))
    }

    /// Build a DOS chart from a periodic KPM result.
    #[pyfunction]
    #[pyo3(signature = (energies, dos_values))]
    pub fn compute_dos_chart(energies: Vec<f64>, dos_values: Vec<f64>) -> PyResult<ChartPayloadPy> {
        use sci_form_core::alpha::periodic_linear::{
            PeriodicBandEdgeSummary, PeriodicKpmDosResult, PeriodicSpectralDiagnostics,
        };
        use sci_form_core::alpha::render_bridge::periodic_dos_chart as build_chart;
        let dos_result = PeriodicKpmDosResult {
            energies_ev: energies,
            total_dos: dos_values,
            kmesh: None,
            band_edges: PeriodicBandEdgeSummary::default(),
            diagnostics: PeriodicSpectralDiagnostics::default(),
        };
        Ok(to_py_chart(&build_chart(&dos_result)))
    }

    /// Build a k-point path chart from the core high-symmetry path.
    ///
    /// The path is the Γ → X → M → Γ default of the core band-structure
    /// adapter, interpolated to `n_points_per_segment` points per segment.
    #[pyfunction]
    #[pyo3(signature = (n_points_per_segment = 20))]
    pub fn compute_kpoint_path_chart(n_points_per_segment: usize) -> PyResult<ChartPayloadPy> {
        use sci_form_core::alpha::periodic_linear::BandStructureAdapterConfig;
        use sci_form_core::render_bridge::{ChartPayload, ChartSeries};
        let path = BandStructureAdapterConfig::default().path;
        let (distance, points) = interpolate_path(&path, n_points_per_segment);
        let chart = ChartPayload {
            title: "K-Point Path".to_string(),
            series: vec![
                ChartSeries {
                    series_id: "k_x".into(),
                    label: "k_x".into(),
                    x: distance.clone(),
                    y: points.iter().map(|k| k[0]).collect(),
                    x_unit: "1/B".into(),
                    y_unit: "1/B".into(),
                },
                ChartSeries {
                    series_id: "k_y".into(),
                    label: "k_y".into(),
                    x: distance.clone(),
                    y: points.iter().map(|k| k[1]).collect(),
                    x_unit: "1/B".into(),
                    y_unit: "1/B".into(),
                },
                ChartSeries {
                    series_id: "k_z".into(),
                    label: "k_z".into(),
                    x: distance,
                    y: points.iter().map(|k| k[2]).collect(),
                    x_unit: "1/B".into(),
                    y_unit: "1/B".into(),
                },
            ],
        };
        Ok(to_py_chart(&chart))
    }

    /// Build a Fermi-surface slice chart from real-space operator blocks.
    ///
    /// The band energies are evaluated on a regular `n_kx` × `n_ky` mesh in the
    /// `kz = 0` plane; the contour where the selected band crosses the Fermi
    /// energy is extracted by bisection.
    #[pyfunction]
    #[pyo3(signature = (blocks_json, fermi_energy_ev = None, n_kx = 24, n_ky = 24))]
    pub fn compute_fermi_surface_chart(
        blocks_json: &str,
        fermi_energy_ev: Option<f64>,
        n_kx: usize,
        n_ky: usize,
    ) -> PyResult<ChartPayloadPy> {
        use sci_form_core::render_bridge::{ChartPayload, ChartSeries};
        let blocks = parse_blocks(blocks_json).map_err(crate::errors::value_error)?;
        let refs = as_core_refs(&blocks);
        let (contour_x, contour_y, fermi) = fermi_contour(&refs, fermi_energy_ev, n_kx, n_ky);
        let chart = ChartPayload {
            title: "Fermi Surface".to_string(),
            series: vec![ChartSeries {
                series_id: "fermi_contour".into(),
                label: format!("E_F = {:.4} eV", fermi),
                x: contour_x,
                y: contour_y,
                x_unit: "1/B".into(),
                y_unit: "1/B".into(),
            }],
        };
        Ok(to_py_chart(&chart))
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(compute_band_structure_chart, m)?)?;
        m.add_function(wrap_pyfunction!(compute_dos_chart, m)?)?;
        m.add_function(wrap_pyfunction!(compute_kpoint_path_chart, m)?)?;
        m.add_function(wrap_pyfunction!(compute_fermi_surface_chart, m)?)?;
        Ok(())
    }
}

// ─── Shared: chart helpers ────────────────────────────────────────────────

/// Helpers shared by the chart bindings.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
mod chart_support {
    use super::kmesh_support::symmetric_eigenvalues;

    /// Convert a flat coordinate array into Cartesian positions.
    pub fn positions_from_flat(coords: &[f64]) -> Vec<[f64; 3]> {
        coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect()
    }

    /// Parse the nine lattice-vector components into a 3×3 lattice.
    pub fn parse_lattice(values: &[f64]) -> Result<[[f64; 3]; 3], String> {
        if values.len() != 9 {
            return Err("lattice needs exactly nine components".into());
        }
        let mut lattice = [[0.0_f64; 3]; 3];
        for (row, chunk) in values.chunks_exact(3).enumerate() {
            lattice[row] = [chunk[0], chunk[1], chunk[2]];
        }
        Ok(lattice)
    }

    /// Interpolate a high-symmetry k-path and return the cumulative path
    /// distance together with the sampled fractional k-points.
    pub fn interpolate_path(
        path: &[([f64; 3], String)],
        n_points_per_segment: usize,
    ) -> (Vec<f64>, Vec<[f64; 3]>) {
        let n = n_points_per_segment.max(1);
        let mut distance = Vec::new();
        let mut points = Vec::new();
        let mut accumulated = 0.0_f64;
        let mut previous: Option<[f64; 3]> = None;
        for index in 0..path.len().saturating_sub(1) {
            let start = path[index].0;
            let end = path[index + 1].0;
            for step in 0..=n {
                if index > 0 && step == 0 {
                    continue;
                }
                let t = step as f64 / n as f64;
                let mut k = [0.0_f64; 3];
                for axis in 0..3 {
                    k[axis] = start[axis] + t * (end[axis] - start[axis]);
                }
                if let Some(prev) = previous {
                    accumulated += ((k[0] - prev[0]).powi(2)
                        + (k[1] - prev[1]).powi(2)
                        + (k[2] - prev[2]).powi(2))
                    .sqrt();
                }
                previous = Some(k);
                distance.push(accumulated);
                points.push(k);
            }
        }
        (distance, points)
    }

    /// Compute the Fermi-surface contour of the selected band in the `kz = 0`
    /// plane and return the contour points plus the Fermi energy used.
    pub fn fermi_contour(
        refs: &super::periodic_blocks::RealSpaceBlockRefs<'_>,
        fermi_energy_ev: Option<f64>,
        n_kx: usize,
        n_ky: usize,
    ) -> (Vec<f64>, Vec<f64>, f64) {
        let n_kx = n_kx.max(2);
        let n_ky = n_ky.max(2);
        let kx_values: Vec<f64> = (0..n_kx).map(|i| -0.5 + i as f64 / n_kx as f64).collect();
        let ky_values: Vec<f64> = (0..n_ky).map(|j| -0.5 + j as f64 / n_ky as f64).collect();

        // Sample the band structure on the 2D mesh.
        let mut samples: Vec<Vec<Vec<f64>>> = Vec::with_capacity(n_ky);
        let mut all_energies: Vec<f64> = Vec::new();
        for &ky in &ky_values {
            let mut row = Vec::with_capacity(n_kx);
            for &kx in &kx_values {
                let hamiltonian =
                    sci_form_core::kmesh::build_bloch_hamiltonian(refs, [kx, ky, 0.0]);
                let eigenvalues = symmetric_eigenvalues(&hamiltonian);
                all_energies.extend(eigenvalues.iter().copied());
                row.push(eigenvalues);
            }
            samples.push(row);
        }

        // Fermi energy: supplied by the caller or estimated as the median of
        // the sampled eigenvalues.
        let fermi = match fermi_energy_ev {
            Some(value) => value,
            None => {
                let mut sorted = all_energies.clone();
                sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                if sorted.is_empty() {
                    0.0
                } else {
                    sorted[sorted.len() / 2]
                }
            }
        };

        // The band that crosses the Fermi energy is the first unoccupied band
        // at the first sampled k-point.
        let band = samples
            .first()
            .and_then(|row| row.first())
            .map(|energies| {
                energies
                    .iter()
                    .filter(|&&e| e < fermi)
                    .count()
                    .min(energies.len().saturating_sub(1))
            })
            .unwrap_or(0);

        let mut contour_x = Vec::new();
        let mut contour_y = Vec::new();
        for (j, &ky) in ky_values.iter().enumerate() {
            for i in 0..n_kx.saturating_sub(1) {
                let e0 = band_energy(&samples[j][i], band);
                let e1 = band_energy(&samples[j][i + 1], band);
                if !e0.is_finite() || !e1.is_finite() {
                    continue;
                }
                let f0 = e0 - fermi;
                let f1 = e1 - fermi;
                if f0 == 0.0 {
                    contour_x.push(kx_values[i]);
                    contour_y.push(ky);
                    continue;
                }
                if f0 * f1 > 0.0 {
                    continue;
                }
                // Bisect on the linearly interpolated band energy.
                let span = kx_values[i + 1] - kx_values[i];
                let energy_at = |x: f64| -> f64 {
                    if span == 0.0 {
                        e0
                    } else {
                        e0 + (x - kx_values[i]) / span * (e1 - e0)
                    }
                };
                let (mut lo, mut hi) = (kx_values[i], kx_values[i + 1]);
                for _ in 0..60 {
                    let mid = 0.5 * (lo + hi);
                    if (energy_at(mid) - fermi) * f0 <= 0.0 {
                        hi = mid;
                    } else {
                        lo = mid;
                    }
                }
                contour_x.push(0.5 * (lo + hi));
                contour_y.push(ky);
            }
        }
        (contour_x, contour_y, fermi)
    }

    /// Energy of a band index at a sampled k-point.
    fn band_energy(energies: &[f64], band: usize) -> f64 {
        energies.get(band).copied().unwrap_or(f64::NAN)
    }
}

// ─── Alpha: GSM charts ────────────────────────────────────────────────────

#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-gsm"))]
mod alpha_gsm_chart_py {
    // `#[pyfunction]` expands to `Err(err.into())`, which newer clippy
    // versions flag as a useless `PyErr` -> `PyErr` conversion.
    #![allow(clippy::useless_conversion)]

    use pyo3::prelude::*;

    use super::alpha_render_py::{to_py_chart, ChartPayloadPy};

    /// Build a trajectory chart from a GSM reaction path.
    ///
    /// The path is computed by the core growing-string method with the
    /// requested energy backend; the chart shows the energy profile and the
    /// accumulated path length of the string.
    #[pyfunction]
    #[pyo3(signature = (smiles, reactant_coords, product_coords, method = "uff"))]
    pub fn compute_trajectory_chart(
        smiles: &str,
        reactant_coords: Vec<f64>,
        product_coords: Vec<f64>,
        method: &str,
    ) -> PyResult<ChartPayloadPy> {
        use sci_form_core::render_bridge::{ChartPayload, ChartSeries};
        let result = gsm_path(smiles, &reactant_coords, &product_coords, method)
            .map_err(|e| crate::errors::smiles_or_engine_error(e, smiles))?;
        let distances = path_distances(&result.path_coords);
        let index: Vec<f64> = (0..result.path_energies.len()).map(|i| i as f64).collect();
        let chart = ChartPayload {
            title: "GSM Trajectory".to_string(),
            series: vec![
                ChartSeries {
                    series_id: "path_energy".into(),
                    label: "E(path)".into(),
                    x: index.clone(),
                    y: result.path_energies.clone(),
                    x_unit: "node".into(),
                    y_unit: "kcal/mol".into(),
                },
                ChartSeries {
                    series_id: "path_distance".into(),
                    label: "path length".into(),
                    x: index,
                    y: distances,
                    x_unit: "node".into(),
                    y_unit: "A".into(),
                },
            ],
        };
        let mut payload = to_py_chart(&chart);
        payload.trajectory = sci_form_core::render_bridge::gsm_path_trajectory(
            &result.path_coords,
            &result.path_energies,
        )
        .into_iter()
        .map(|frame| super::alpha_render_py::TrajectoryFramePy {
            frame_time: frame.frame_time,
            coords: frame.coords,
        })
        .collect();
        Ok(payload)
    }

    /// Build a reaction-coordinate chart from a GSM reaction path.
    ///
    /// The reaction coordinate is the accumulated path length in Å; the chart
    /// shows the potential energy along that coordinate.
    #[pyfunction]
    #[pyo3(signature = (smiles, reactant_coords, product_coords, method = "uff"))]
    pub fn compute_reaction_coordinate_chart(
        smiles: &str,
        reactant_coords: Vec<f64>,
        product_coords: Vec<f64>,
        method: &str,
    ) -> PyResult<ChartPayloadPy> {
        use sci_form_core::render_bridge::{ChartPayload, ChartSeries};
        let result = gsm_path(smiles, &reactant_coords, &product_coords, method)
            .map_err(|e| crate::errors::smiles_or_engine_error(e, smiles))?;
        let distances = path_distances(&result.path_coords);
        let chart = ChartPayload {
            title: "Reaction Coordinate".to_string(),
            series: vec![ChartSeries {
                series_id: "reaction_coordinate".into(),
                label: "E(rc)".into(),
                x: distances,
                y: result.path_energies.clone(),
                x_unit: "A".into(),
                y_unit: "kcal/mol".into(),
            }],
        };
        Ok(to_py_chart(&chart))
    }

    /// Run the core GSM search with the requested energy backend.
    fn gsm_path(
        smiles: &str,
        reactant: &[f64],
        product: &[f64],
        method: &str,
    ) -> Result<sci_form_core::alpha::gsm::GsmResult, String> {
        use sci_form_core::alpha::gsm::{
            find_transition_state_with_backend, GsmConfig, GsmEnergyBackend,
        };
        let backend: GsmEnergyBackend = method.parse().map_err(|_e: String| {
            crate::errors::invalid_value(
                "'method'",
                method,
                &[
                    "uff", "mmff94", "pm3", "xtb", "gfn1", "gfn2", "hf3c", "reaxff",
                ],
            )
        })?;
        find_transition_state_with_backend(
            smiles,
            reactant,
            product,
            backend,
            &GsmConfig::default(),
        )
    }

    /// Accumulated path length of the string nodes.
    fn path_distances(path_coords: &[Vec<f64>]) -> Vec<f64> {
        let mut distances = Vec::with_capacity(path_coords.len());
        let mut accumulated = 0.0_f64;
        let mut previous: Option<&Vec<f64>> = None;
        for coords in path_coords {
            if let Some(prev) = previous {
                accumulated += prev
                    .iter()
                    .zip(coords.iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt();
            }
            previous = Some(coords);
            distances.push(accumulated);
        }
        distances
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(compute_trajectory_chart, m)?)?;
        m.add_function(wrap_pyfunction!(compute_reaction_coordinate_chart, m)?)?;
        Ok(())
    }
}

// ─── Shared: k-mesh helpers ───────────────────────────────────────────────

/// Helpers shared by the periodic bindings.
#[cfg(feature = "alpha-periodic-linear")]
mod kmesh_support {
    use nalgebra::DMatrix;
    use sci_form_core::kmesh::KMeshCentering;

    /// Parse a k-mesh centering label.
    pub fn parse_centering(centering: &str) -> Result<KMeshCentering, String> {
        match centering {
            "monkhorst-pack" | "mp" => Ok(KMeshCentering::MonkhorstPack),
            "gamma" | "gamma-centered" => Ok(KMeshCentering::GammaCentered),
            _ => Err(format!(
                "Invalid 'centering' '{centering}'. Supported: 'monkhorst-pack', 'mp', 'gamma', 'gamma-centered'."
            )),
        }
    }

    /// Row-major representation of a square matrix.
    pub fn matrix_rows(matrix: &DMatrix<f64>) -> Vec<Vec<f64>> {
        let n = matrix.nrows();
        (0..n)
            .map(|i| (0..n).map(|j| matrix[(i, j)]).collect())
            .collect()
    }

    /// Ascending eigenvalues of a symmetric matrix.
    pub fn symmetric_eigenvalues(matrix: &DMatrix<f64>) -> Vec<f64> {
        let eigen = nalgebra::SymmetricEigen::new(matrix.clone());
        let mut values: Vec<f64> = eigen.eigenvalues.iter().copied().collect();
        values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        values
    }
}
