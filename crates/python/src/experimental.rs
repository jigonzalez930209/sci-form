//! Python bindings for experimental modules.
//!
//! Each experimental track is gated behind its own feature flag.

use pyo3::prelude::*;

// ─── E5: EEQ ───────────────────────────────────────────────────────────────

#[cfg(feature = "experimental-eeq")]
mod eeq_py {
    use pyo3::prelude::*;

    #[pyclass]
    #[derive(Clone)]
    pub struct EeqChargeResultPy {
        #[pyo3(get)]
        pub charges: Vec<f64>,
        #[pyo3(get)]
        pub coordination_numbers: Vec<f64>,
        #[pyo3(get)]
        pub total_charge: f64,
    }

    #[pyclass]
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
    pub fn eeq_charges(
        elements: Vec<u8>,
        coords: Vec<f64>,
        total_charge: f64,
    ) -> PyResult<EeqChargeResultPy> {
        let pos: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let config = sci_form_core::experimental::eeq::EeqConfig {
            total_charge,
            regularization: 1e-10,
        };
        let r = sci_form_core::experimental::eeq::compute_eeq_charges(&elements, &pos, &config);
        Ok(EeqChargeResultPy {
            charges: r.charges,
            coordination_numbers: r.coordination_numbers,
            total_charge: r.total_charge,
        })
    }

    /// Compute EEQ electrostatic energy.
    #[pyfunction]
    #[pyo3(signature = (elements, coords, total_charge=0.0))]
    pub fn eeq_energy(
        elements: Vec<u8>,
        coords: Vec<f64>,
        total_charge: f64,
    ) -> PyResult<EeqEnergyResultPy> {
        let pos: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let config = sci_form_core::experimental::eeq::EeqConfig {
            total_charge,
            regularization: 1e-10,
        };
        let r = sci_form_core::experimental::eeq::compute_eeq_energy(&elements, &pos, &config);
        Ok(EeqEnergyResultPy {
            electrostatic_energy: r.electrostatic_energy,
            charges: r.charges,
            coordination_numbers: r.coordination_numbers,
        })
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(eeq_charges, m)?)?;
        m.add_function(wrap_pyfunction!(eeq_energy, m)?)?;
        Ok(())
    }
}

// ─── E6: ALPB ──────────────────────────────────────────────────────────────

#[cfg(feature = "experimental-alpb")]
mod alpb_py {
    use pyo3::prelude::*;

    #[pyclass]
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

    #[pyclass]
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
    pub fn alpb_solvation(
        elements: Vec<u8>,
        coords: Vec<f64>,
        charges: Vec<f64>,
        solvent_dielectric: f64,
        probe_radius: f64,
        surface_tension: f64,
    ) -> PyResult<AlpbResultPy> {
        let pos: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let config = sci_form_core::experimental::alpb::AlpbConfig {
            solvent_dielectric,
            probe_radius,
            surface_tension,
        };
        let r = sci_form_core::experimental::alpb::compute_alpb_solvation(
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
    pub fn alpb_born_radii(
        elements: Vec<u8>,
        coords: Vec<f64>,
        probe_radius: f64,
    ) -> PyResult<AlpbBornRadiiPy> {
        let pos: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let r = sci_form_core::experimental::alpb::compute_born_radii(
            &elements, &pos, probe_radius,
        );
        Ok(AlpbBornRadiiPy {
            radii: r.radii,
            intrinsic: r.intrinsic,
        })
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(alpb_solvation, m)?)?;
        m.add_function(wrap_pyfunction!(alpb_born_radii, m)?)?;
        Ok(())
    }
}

// ─── E7: D4 ────────────────────────────────────────────────────────────────

#[cfg(feature = "experimental-d4")]
mod d4_py {
    use pyo3::prelude::*;

    #[pyclass]
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
    pub fn d4_energy(
        elements: Vec<u8>,
        coords: Vec<f64>,
        s6: f64,
        s8: f64,
        a1: f64,
        a2: f64,
        three_body: bool,
    ) -> PyResult<D4ResultPy> {
        let pos: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let config = sci_form_core::experimental::d4::D4Config {
            s6, s8, a1, a2, three_body, s9: 1.0,
        };
        let r = sci_form_core::experimental::d4::compute_d4_energy(&elements, &pos, &config);
        Ok(D4ResultPy {
            e2_body: r.e2_body,
            e3_body: r.e3_body,
            total_energy: r.total_energy,
            total_kcal_mol: r.total_kcal_mol,
            coordination_numbers: r.coordination_numbers,
        })
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(d4_energy, m)?)?;
        Ok(())
    }
}

// ─── E10: CPM ──────────────────────────────────────────────────────────────

#[cfg(feature = "experimental-cpm")]
mod cpm_py {
    use pyo3::prelude::*;

    #[pyclass]
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

    #[pyclass]
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
    pub fn cpm_charges(
        elements: Vec<u8>,
        coords: Vec<f64>,
        mu_ev: f64,
        dielectric: f64,
    ) -> PyResult<CpmResultPy> {
        let pos: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let config = sci_form_core::experimental::cpm::CpmConfig {
            mu_ev,
            dielectric,
            max_iter: 100,
            charge_tol: 1e-6,
        };
        let r = sci_form_core::experimental::cpm::compute_cpm_charges(&elements, &pos, &config);
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
    pub fn cpm_surface(
        elements: Vec<u8>,
        coords: Vec<f64>,
        mu_min: f64,
        mu_max: f64,
        n_points: usize,
        dielectric: f64,
    ) -> PyResult<CpmSurfacePy> {
        let pos: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let r = sci_form_core::experimental::cpm::compute_cpm_surface(
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
        m.add_function(wrap_pyfunction!(cpm_charges, m)?)?;
        m.add_function(wrap_pyfunction!(cpm_surface, m)?)?;
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

    #[cfg(feature = "experimental-cpm")]
    cpm_py::register(_m)?;

    Ok(())
}
