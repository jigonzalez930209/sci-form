//! Types for advanced spectroscopy methods.
//!
//! Self-contained types that avoid dependency on experimental_2::types.

use nalgebra::DMatrix;

/// Converged SCF data needed by spectroscopy methods.
///
/// A lightweight input struct that can be constructed from any SCF solver
/// (HF-3c, PM3, EHT, xTB) without coupling to a specific implementation.
#[derive(Debug, Clone)]
pub struct ScfInput {
    /// Orbital energies sorted ascending (Hartree).
    pub orbital_energies: Vec<f64>,
    /// MO coefficient matrix C (n_basis × n_basis).
    pub mo_coefficients: DMatrix<f64>,
    /// Density matrix P (n_basis × n_basis).
    pub density_matrix: DMatrix<f64>,
    /// Overlap matrix S (n_basis × n_basis).
    pub overlap_matrix: DMatrix<f64>,
    /// Number of basis functions.
    pub n_basis: usize,
    /// Number of electrons.
    pub n_electrons: usize,
}

/// Information about a single electronic transition (UV-Vis).
#[derive(Debug, Clone)]
pub struct TransitionInfo {
    /// Excitation energy (eV).
    pub energy_ev: f64,
    /// Wavelength (nm).
    pub wavelength_nm: f64,
    /// Oscillator strength (dimensionless).
    pub oscillator_strength: f64,
    /// Transition dipole moment [x, y, z] in atomic units.
    pub transition_dipole: [f64; 3],
}

/// Result of sTDA UV-Vis calculation.
#[derive(Debug, Clone)]
pub struct SpectroscopyResult {
    /// Electronic transitions with energies and oscillator strengths.
    pub transitions: Vec<TransitionInfo>,
    /// Method used.
    pub method: String,
}

/// NMR shielding tensor for a single nucleus.
#[derive(Debug, Clone)]
pub struct ShieldingTensor {
    /// Atom index.
    pub atom_index: usize,
    /// Element (atomic number).
    pub element: u8,
    /// 3×3 shielding tensor (ppm).
    pub tensor: [[f64; 3]; 3],
    /// Isotropic shielding σ_iso = Tr(σ)/3 (ppm).
    pub isotropic: f64,
    /// Anisotropic shielding (ppm).
    pub anisotropy: f64,
    /// Chemical shift δ = σ_ref − σ (ppm).
    pub chemical_shift: f64,
}

/// Result of GIAO NMR calculation.
#[derive(Debug, Clone)]
pub struct NmrShieldingResult {
    /// Chemical shifts (ppm) per atom.
    pub chemical_shifts: Vec<f64>,
    /// Element types (atomic numbers).
    pub elements: Vec<u8>,
    /// Number of atoms.
    pub n_atoms: usize,
}
