//! Shared types for SCF calculations.

use nalgebra::DMatrix;

/// Represents a molecular system with all data needed for quantum calculations.
#[derive(Debug, Clone)]
pub struct MolecularSystem {
    /// Atomic numbers (Z) for each atom.
    pub atomic_numbers: Vec<u8>,
    /// Cartesian coordinates in Bohr, shape `[n_atoms][3]`.
    pub positions_bohr: Vec<[f64; 3]>,
    /// Total molecular charge in units of `e` (0 = neutral).
    pub charge: i32,
    /// Spin multiplicity `2S + 1` (1 = singlet, 2 = doublet, etc.).
    pub multiplicity: u32,
}

impl MolecularSystem {
    /// Create from atomic numbers and Angstrom coordinates.
    pub fn from_angstrom(
        elements: &[u8],
        positions_ang: &[[f64; 3]],
        charge: i32,
        multiplicity: u32,
    ) -> Self {
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

    /// Create from atomic numbers and Angstrom coordinates, validating the
    /// result.
    ///
    /// Rejects a length mismatch between `elements` and `positions_ang`, a
    /// non-finite coordinate, an empty system, a `charge` that would make
    /// the electron count negative, and a `multiplicity` that is
    /// incompatible with the resulting number of electrons.
    pub fn try_from_angstrom(
        elements: &[u8],
        positions_ang: &[[f64; 3]],
        charge: i32,
        multiplicity: u32,
    ) -> Result<Self, String> {
        let system = Self::from_angstrom(elements, positions_ang, charge, multiplicity);
        system.validate()?;
        Ok(system)
    }

    /// Create from flat Angstrom coordinates [x0,y0,z0,x1,y1,z1,...],
    /// validating the result.
    pub fn try_from_flat_angstrom(
        elements: &[u8],
        coords_flat: &[f64],
        charge: i32,
        multiplicity: u32,
    ) -> Result<Self, String> {
        crate::validation::validate_flat_geometry(elements, coords_flat)?;
        let positions: Vec<[f64; 3]> = coords_flat
            .chunks_exact(3)
            .map(|c| [c[0], c[1], c[2]])
            .collect();
        Self::try_from_angstrom(elements, &positions, charge, multiplicity)
    }

    /// Validate the geometry, the total charge and the spin multiplicity.
    ///
    /// The charge/multiplicity messages follow the wording already used by
    /// the xTB unrestricted solver (`src/xtb/uks.rs`).
    pub fn validate(&self) -> Result<(), String> {
        self.validate_geometry()?;

        let z_total: i32 = self.atomic_numbers.iter().map(|&z| z as i32).sum();
        let n_el = z_total - self.charge;
        if n_el < 0 {
            return Err(format!(
                "charge={} is incompatible with {} nuclear charges (the electron count would be {}; charge must be between {} and {})",
                self.charge,
                z_total,
                n_el,
                -z_total,
                z_total
            ));
        }
        let n_electrons = n_el as usize;

        if self.multiplicity == 0 {
            return Err("multiplicity must be >= 1 (got 0)".to_string());
        }
        let n_unpaired = (self.multiplicity as usize).saturating_sub(1);
        if n_unpaired > n_electrons || !(n_electrons + n_unpaired).is_multiple_of(2) {
            return Err(format!(
                "multiplicity={} is incompatible with {} electrons (need N±(M-1) even, 0 <= M-1 <= N)",
                self.multiplicity, n_electrons
            ));
        }
        Ok(())
    }

    /// Validate only the geometry (atom count, non-empty system, finite
    /// coordinates).
    ///
    /// The RHF driver uses this because a closed-shell calculation does not
    /// care about the requested spin multiplicity.
    pub fn validate_geometry(&self) -> Result<(), String> {
        crate::validation::validate_geometry(&self.atomic_numbers, &self.positions_bohr)
    }

    /// Number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.atomic_numbers.len()
    }

    /// Total number of electrons (sum of Z minus charge).
    ///
    /// A negative electron count (a `charge` larger than the total nuclear
    /// charge) yields `0` instead of wrapping around to a huge `usize`.
    pub fn n_electrons(&self) -> usize {
        let z_total: i32 = self.atomic_numbers.iter().map(|&z| z as i32).sum();
        usize::try_from(z_total - self.charge).unwrap_or(0)
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
    /// MO coefficient matrix C (n_basis × n_basis), dimensionless.
    pub mo_coefficients: DMatrix<f64>,
    /// Density matrix P (n_basis × n_basis), dimensionless.
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
    /// Mulliken charges per atom, in units of `e`.
    pub mulliken_charges: Vec<f64>,
    /// Number of SCF iterations performed.
    pub scf_iterations: usize,
    /// Whether SCF converged.
    pub converged: bool,
    /// Number of basis functions.
    pub n_basis: usize,
    /// Number of electrons.
    pub n_electrons: usize,
    /// Overlap matrix S (n_basis × n_basis), dimensionless.
    pub overlap_matrix: DMatrix<f64>,
    /// Fock matrix F at convergence, in Hartree.
    pub fock_matrix: DMatrix<f64>,
    /// Error message when the calculation could not be run at all
    /// (invalid geometry). `None` for a completed calculation.
    pub error: Option<String>,
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
    /// Error message when the calculation could not be run at all
    /// (invalid geometry). `None` for a completed calculation.
    pub error: Option<String>,
}
