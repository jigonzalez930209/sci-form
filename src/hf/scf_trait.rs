//! Unified SCF solver trait for HF, PM3, and xTB backends.
//!
//! Provides a common interface for self-consistent field procedures
//! across different quantum chemistry methods.

use nalgebra::DMatrix;

/// Common SCF convergence result.
#[derive(Debug, Clone)]
pub struct ScfOutput {
    /// Converged total energy (method-specific units).
    ///
    /// Hartree for the Hartree-Fock backend, eV for the PM3 and xTB
    /// backends — the same units the underlying solver reports.
    pub energy: f64,
    /// Orbital energies (eigenvalues), in the method-specific unit of
    /// [`ScfOutput::energy`] (Hartree for HF, eV for PM3 and xTB).
    pub orbital_energies: Vec<f64>,
    /// Number of SCF iterations.
    pub iterations: usize,
    /// Whether SCF converged.
    pub converged: bool,
    /// Converged density matrix.
    pub density: Option<DMatrix<f64>>,
    /// Mulliken charges (if available).
    pub mulliken_charges: Option<Vec<f64>>,
}

/// Configuration for SCF convergence.
#[derive(Debug, Clone)]
pub struct ScfConvergenceConfig {
    /// Maximum number of SCF iterations.
    pub max_iter: usize,
    /// Energy convergence threshold in Hartree.
    pub energy_threshold: f64,
    /// Density convergence threshold.
    ///
    /// The RMS change of the density matrix elements between two iterations
    /// must fall below this value.
    pub density_threshold: f64,
    /// DIIS history size (0 to disable).
    pub diis_size: usize,
    /// Level shift for virtual orbitals in Hartree (0.0 to disable).
    pub level_shift: f64,
    /// Enable ADIIS for initial iterations before switching to DIIS.
    pub use_adiis: bool,
    /// Iteration threshold for switching from ADIIS to DIIS.
    pub adiis_switch_iter: usize,
}

impl Default for ScfConvergenceConfig {
    fn default() -> Self {
        ScfConvergenceConfig {
            max_iter: 100,
            energy_threshold: 1e-8,
            density_threshold: 1e-6,
            diis_size: 8,
            level_shift: 0.0,
            use_adiis: false,
            adiis_switch_iter: 5,
        }
    }
}

/// Unified trait for SCF solvers (HF, PM3, xTB).
pub trait ScfSolver {
    /// Run the SCF procedure to convergence.
    fn solve(&self, config: &ScfConvergenceConfig) -> Result<ScfOutput, String>;

    /// Get the method name for display.
    fn method_name(&self) -> &str;

    /// Get the number of basis functions.
    fn n_basis(&self) -> usize;

    /// Get the number of electrons.
    fn n_electrons(&self) -> usize;
}

/// HF SCF solver wrapping the existing Roothaan-Hall implementation.
pub struct HfScfSolver {
    /// Core Hamiltonian matrix `H_core` in the AO basis, in Hartree.
    pub h_core: DMatrix<f64>,
    /// Overlap matrix `S` in the AO basis (dimensionless).
    pub s_mat: DMatrix<f64>,
    /// Packed electron-repulsion integrals in Hartree, in the layout of
    /// [`super::integrals::compute_eris`].
    pub eris: Vec<f64>,
    /// Number of electrons in the system (closed shell: twice the number of
    /// occupied spatial orbitals).
    pub n_elec: usize,
}

impl ScfSolver for HfScfSolver {
    fn solve(&self, config: &ScfConvergenceConfig) -> Result<ScfOutput, String> {
        let hf_config = super::scf::ScfConfig {
            max_iter: config.max_iter,
            energy_threshold: config.energy_threshold,
            density_threshold: config.density_threshold,
            diis_size: config.diis_size,
            level_shift: config.level_shift,
        };

        let result = super::scf::solve_scf(
            &self.h_core,
            &self.s_mat,
            &self.eris,
            None,
            self.n_elec,
            &hf_config,
        );

        Ok(ScfOutput {
            energy: result.energy,
            orbital_energies: result.orbital_energies,
            iterations: result.iterations,
            converged: result.converged,
            density: Some(result.density),
            mulliken_charges: None,
        })
    }

    fn method_name(&self) -> &str {
        "HF"
    }
    fn n_basis(&self) -> usize {
        self.h_core.nrows()
    }
    fn n_electrons(&self) -> usize {
        self.n_elec
    }
}

/// PM3 SCF solver wrapping the existing NDDO implementation.
pub struct Pm3ScfSolver {
    /// Atomic numbers of the molecule.
    pub elements: Vec<u8>,
    /// Cartesian coordinates in Ångström.
    pub positions: Vec<[f64; 3]>,
}

impl ScfSolver for Pm3ScfSolver {
    fn solve(&self, _config: &ScfConvergenceConfig) -> Result<ScfOutput, String> {
        let result = crate::pm3::solver::solve_pm3(&self.elements, &self.positions)?;
        Ok(ScfOutput {
            energy: result.total_energy,
            orbital_energies: result.orbital_energies,
            iterations: result.scf_iterations,
            converged: result.converged,
            density: None,
            mulliken_charges: Some(result.mulliken_charges),
        })
    }

    fn method_name(&self) -> &str {
        "PM3"
    }
    fn n_basis(&self) -> usize {
        0
    } // not directly accessible
    fn n_electrons(&self) -> usize {
        0
    }
}

/// xTB SCC solver wrapping the existing tight-binding implementation.
pub struct XtbScfSolver {
    /// Atomic numbers of the molecule.
    pub elements: Vec<u8>,
    /// Cartesian coordinates in Ångström.
    pub positions: Vec<[f64; 3]>,
}

impl ScfSolver for XtbScfSolver {
    fn solve(&self, _config: &ScfConvergenceConfig) -> Result<ScfOutput, String> {
        let result = crate::xtb::gfn2::solve_gfn2(&self.elements, &self.positions)?;
        Ok(ScfOutput {
            energy: result.total_energy,
            orbital_energies: result.orbital_energies,
            iterations: result.scc_iterations,
            converged: result.converged,
            density: None,
            mulliken_charges: Some(result.mulliken_charges),
        })
    }

    fn method_name(&self) -> &str {
        "xTB"
    }
    fn n_basis(&self) -> usize {
        0
    }
    fn n_electrons(&self) -> usize {
        0
    }
}
