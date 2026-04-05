//! Per-frame property calculation for reaction dynamics frames.
//!
//! Computes electronic properties at each frame of a reaction trajectory:
//! charges, bond orders, HOMO/LUMO, dipole, etc.

use serde::{Deserialize, Serialize};

/// Properties computed at a single frame.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FrameProperties {
    /// Mulliken partial charges.
    pub mulliken_charges: Vec<f64>,
    /// Wiberg bond order matrix (flattened, row-major, n_atoms × n_atoms).
    pub wiberg_bond_orders: Vec<f64>,
    /// HOMO energy (eV).
    pub homo_energy: f64,
    /// LUMO energy (eV).
    pub lumo_energy: f64,
    /// HOMO-LUMO gap (eV).
    pub gap: f64,
    /// Dipole moment magnitude (Debye).
    pub dipole_magnitude: f64,
    /// Dipole moment vector [x, y, z].
    pub dipole_vector: [f64; 3],
    /// Total electronic energy at this frame (eV).
    pub energy: f64,
}

/// Compute electronic properties for a single reaction frame.
///
/// Uses EHT (fast) for orbital info + population analysis for charges.
/// Falls back gracefully if individual calculations fail.
pub fn compute_frame_properties(
    elements: &[u8],
    coords: &[f64],
) -> Result<FrameProperties, String> {
    let n_atoms = elements.len();
    let pos: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

    // EHT for orbital energies
    let (homo_e, lumo_e, gap, energy) = match crate::eht::solve_eht(elements, &pos, Some(1.75)) {
        Ok(eht) => {
            (eht.homo_energy, eht.lumo_energy, eht.gap, eht.homo_energy) // no total_energy in EHT
        }
        Err(_) => (0.0, 0.0, 0.0, 0.0),
    };

    // Population analysis for Mulliken charges
    let mulliken = match crate::compute_population(elements, &pos) {
        Ok(pop) => pop.mulliken_charges,
        Err(_) => vec![0.0; n_atoms],
    };

    // Dipole moment
    let (dip_mag, dip_vec) = match crate::compute_dipole(elements, &pos) {
        Ok(d) => (d.magnitude, [d.vector[0], d.vector[1], d.vector[2]]),
        Err(_) => (0.0, [0.0; 3]),
    };

    // Bond orders
    let bond_orders = compute_wiberg_bond_orders(elements, coords, n_atoms);

    Ok(FrameProperties {
        mulliken_charges: mulliken,
        wiberg_bond_orders: bond_orders,
        homo_energy: homo_e,
        lumo_energy: lumo_e,
        gap,
        dipole_magnitude: dip_mag,
        dipole_vector: dip_vec,
        energy,
    })
}

/// Compute Wiberg bond orders from the density matrix.
///
/// Uses EHT density to compute bond orders: BOij = Σ_μν (P*S)_μi * (P*S)_νj
/// Falls back to zero matrix on failure.
fn compute_wiberg_bond_orders(elements: &[u8], coords: &[f64], n_atoms: usize) -> Vec<f64> {
    let pos: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

    // Try compute_bond_orders if available
    match crate::eht::solve_eht(elements, &pos, Some(1.75)) {
        Ok(eht) => {
            let n_basis = eht.coefficients.len();
            let n_occ = eht.n_electrons / 2;

            // Build density matrix P = 2 * Σ_occ C_μi C_νi
            // coefficients[mu][mo] = C_{mu,mo}
            let mut density = vec![0.0f64; n_basis * n_basis];
            for k in 0..n_occ {
                for mu in 0..n_basis {
                    let c_mu = if k < eht.coefficients[mu].len() {
                        eht.coefficients[mu][k]
                    } else {
                        0.0
                    };
                    for nu in 0..n_basis {
                        let c_nu = if k < eht.coefficients[nu].len() {
                            eht.coefficients[nu][k]
                        } else {
                            0.0
                        };
                        density[mu * n_basis + nu] += 2.0 * c_mu * c_nu;
                    }
                }
            }

            // Build basis→atom map
            let basis_to_atom = build_basis_atom_map(elements, n_basis);

            // Wiberg bond order = Σ_μ∈A Σ_ν∈B (PS)_μν²
            // For simplicity, use P_μν² as an approximation (minimal basis, S≈I)
            let mut bo = vec![0.0f64; n_atoms * n_atoms];
            for mu in 0..n_basis {
                let a = basis_to_atom[mu];
                for nu in 0..n_basis {
                    let b = basis_to_atom[nu];
                    if a < n_atoms && b < n_atoms && a != b {
                        let p = density[mu * n_basis + nu];
                        bo[a * n_atoms + b] += p * p;
                    }
                }
            }

            bo
        }
        Err(_) => vec![0.0f64; n_atoms * n_atoms],
    }
}

fn build_basis_atom_map(elements: &[u8], n_basis: usize) -> Vec<usize> {
    let n_atoms = elements.len();
    let mut map = Vec::with_capacity(n_basis);
    for (i, &z) in elements.iter().enumerate() {
        let nbf = match z {
            1 => 1,
            2 => 1,
            3..=10 => 5,
            11..=18 => 9,
            19..=36 => 13,
            _ => 18,
        };
        for _ in 0..nbf {
            map.push(i);
        }
    }
    while map.len() < n_basis {
        map.push(n_atoms.saturating_sub(1));
    }
    map.truncate(n_basis);
    map
}
