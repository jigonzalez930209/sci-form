//! Shared types used across all phases of the experimental engine.

use nalgebra::DMatrix;

/// Represents a molecular system with all data needed for quantum calculations.
#[derive(Debug, Clone)]
pub struct MolecularSystem {
    /// Atomic numbers (Z) for each atom.
    pub atomic_numbers: Vec<u8>,
    /// Cartesian coordinates in Bohr, shape [n_atoms][3].
    pub positions_bohr: Vec<[f64; 3]>,
    /// Total molecular charge (0 = neutral).
    pub charge: i32,
    /// Spin multiplicity (1 = singlet, 2 = doublet, etc.).
    pub multiplicity: u32,
}

impl MolecularSystem {
    /// Create from atomic numbers and Angstrom coordinates.
    pub fn from_angstrom(elements: &[u8], positions_ang: &[[f64; 3]], charge: i32, multiplicity: u32) -> Self {
        let positions_bohr = positions_ang
            .iter()
            .map(|p| {
                [
                    p[0] * super::constants::ANGSTROM_TO_BOHR,
                    p[1] * super::constants::ANGSTROM_TO_BOHR,
                    p[2] * super::constants::ANGSTROM_TO_BOHR,
                ]
            })
            .collect();

        Self {
            atomic_numbers: elements.to_vec(),
            positions_bohr,
            charge,
            multiplicity,
        }
    }

    /// Create from flat Angstrom coordinates [x0,y0,z0,x1,y1,z1,...].
    pub fn from_flat_angstrom(elements: &[u8], coords_flat: &[f64], charge: i32, multiplicity: u32) -> Self {
        let positions: Vec<[f64; 3]> = coords_flat
            .chunks(3)
            .map(|c| [c[0], c[1], c[2]])
            .collect();
        Self::from_angstrom(elements, &positions, charge, multiplicity)
    }

    /// Number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.atomic_numbers.len()
    }

    /// Total number of electrons (sum of Z minus charge).
    pub fn n_electrons(&self) -> usize {
        let z_total: i32 = self.atomic_numbers.iter().map(|&z| z as i32).sum();
        (z_total - self.charge) as usize
    }

    /// Number of occupied orbitals (RHF, closed-shell).
    pub fn n_occupied(&self) -> usize {
        self.n_electrons() / 2
    }

    /// Distance between atoms i and j in Bohr.
    pub fn distance_bohr(&self, i: usize, j: usize) -> f64 {
        let a = self.positions_bohr[i];
        let b = self.positions_bohr[j];
        let dx = a[0] - b[0];
        let dy = a[1] - b[1];
        let dz = a[2] - b[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

/// Result of a complete self-consistent field calculation.
#[derive(Debug, Clone)]
pub struct ScfResult {
    /// Orbital energies sorted ascending (Hartree).
    pub orbital_energies: Vec<f64>,
    /// MO coefficient matrix C (n_basis × n_basis), column-major.
    pub mo_coefficients: DMatrix<f64>,
    /// Density matrix P (n_basis × n_basis).
    pub density_matrix: DMatrix<f64>,
    /// Total electronic energy (Hartree).
    pub electronic_energy: f64,
    /// Nuclear repulsion energy (Hartree).
    pub nuclear_repulsion: f64,
    /// Total energy = electronic + nuclear (Hartree).
    pub total_energy: f64,
    /// HOMO energy (Hartree).
    pub homo_energy: f64,
    /// LUMO energy (Hartree), None if all orbitals occupied.
    pub lumo_energy: Option<f64>,
    /// HOMO-LUMO gap (eV).
    pub gap_ev: f64,
    /// Mulliken charges per atom.
    pub mulliken_charges: Vec<f64>,
    /// Number of SCF iterations performed.
    pub scf_iterations: usize,
    /// Whether SCF converged.
    pub converged: bool,
    /// Number of basis functions.
    pub n_basis: usize,
    /// Number of electrons.
    pub n_electrons: usize,
    /// Overlap matrix S.
    pub overlap_matrix: DMatrix<f64>,
    /// Fock matrix F at convergence.
    pub fock_matrix: DMatrix<f64>,
}

/// Result of a geometry optimization.
#[derive(Debug, Clone)]
pub struct OptimizationResult {
    /// Optimized positions in Bohr.
    pub optimized_positions: Vec<[f64; 3]>,
    /// Final total energy (Hartree).
    pub final_energy: f64,
    /// Whether optimization converged.
    pub converged: bool,
    /// Number of optimization steps.
    pub n_steps: usize,
    /// Final gradient norm (Hartree/Bohr).
    pub gradient_norm: f64,
    /// Energy at each optimization step (Hartree).
    pub energy_trajectory: Vec<f64>,
}

/// Result of a spectroscopic calculation.
#[derive(Debug, Clone)]
pub struct SpectroscopyResult {
    /// Electronic transitions with energies and oscillator strengths.
    pub transitions: Vec<TransitionInfo>,
    /// Method used (e.g. "sTDA", "TD-DFT").
    pub method: String,
}

/// Information about a single electronic transition.
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

/// Result of an IR spectroscopy calculation.
#[derive(Debug, Clone)]
pub struct IrResult {
    /// Vibrational frequencies (cm⁻¹).
    pub frequencies: Vec<f64>,
    /// IR intensities (km/mol).
    pub intensities: Vec<f64>,
    /// Number of vibrational modes.
    pub n_modes: usize,
}

/// Result of NMR chemical shift calculation.
#[derive(Debug, Clone)]
pub struct NmrResult {
    /// Chemical shifts (ppm) per atom, relative to reference.
    pub chemical_shifts: Vec<f64>,
    /// Element types (atomic numbers) for each shift.
    pub elements: Vec<u8>,
    /// Number of atoms.
    pub n_atoms: usize,
}
