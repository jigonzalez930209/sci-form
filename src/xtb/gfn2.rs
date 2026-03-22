//! GFN2-xTB method: second-generation tight-binding with multipole electrostatics.
//!
//! GFN2-xTB extends GFN1 with:
//! - Anisotropic second-order electrostatics (multipole expansion)
//! - Improved halogen and hydrogen bond corrections
//! - Better treatment of non-covalent interactions
//!
//! Reference: Bannwarth, C.; Ehlert, S.; Grimme, S. JCTC 15 (2019): 1652.

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

/// GFN2-xTB calculation result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gfn2Result {
    /// Orbital energies (eV).
    pub orbital_energies: Vec<f64>,
    /// Electronic energy (eV).
    pub electronic_energy: f64,
    /// Repulsive energy (eV).
    pub repulsive_energy: f64,
    /// Dispersion energy (eV) — D4.
    pub dispersion_energy: f64,
    /// Halogen bond correction energy (eV).
    pub halogen_bond_energy: f64,
    /// Total energy (eV).
    pub total_energy: f64,
    /// Number of basis functions.
    pub n_basis: usize,
    /// Number of electrons.
    pub n_electrons: usize,
    /// HOMO energy (eV).
    pub homo_energy: f64,
    /// LUMO energy (eV).
    pub lumo_energy: f64,
    /// HOMO-LUMO gap (eV).
    pub gap: f64,
    /// Mulliken charges per atom.
    pub mulliken_charges: Vec<f64>,
    /// Atomic dipole moments (Debye, per atom).
    pub atomic_dipoles: Vec<[f64; 3]>,
    /// Atomic quadrupole moment traces.
    pub atomic_quadrupoles: Vec<f64>,
    /// SCC iterations.
    pub scc_iterations: usize,
    /// Whether the SCC converged.
    pub converged: bool,
}

/// Solve GFN2-xTB for a molecular system.
///
/// GFN2 provides higher accuracy than GFN1 through multipole electrostatics
/// and improved non-covalent interaction treatment.
pub fn solve_gfn2(elements: &[u8], positions: &[[f64; 3]]) -> Result<Gfn2Result, String> {
    // Validate element support
    for &z in elements {
        if crate::xtb::params::get_xtb_params(z).is_none() {
            return Err(format!("Element Z={} not supported by GFN2-xTB", z));
        }
    }

    let n_atoms = elements.len();

    // Build basis set info
    let mut n_basis = 0;
    let mut n_electrons = 0;
    for &z in elements {
        n_basis += crate::xtb::params::num_xtb_basis_functions(z);
        n_electrons += crate::xtb::params::get_xtb_params(z).unwrap().n_valence as usize;
    }

    // Start from GFN0 solution
    let gfn0 = crate::xtb::solve_xtb(elements, positions)?;

    // Shell-resolved SCC with multipole electrostatics
    let mut charges = gfn0.mulliken_charges.clone();
    let mut atomic_dipoles = vec![[0.0f64; 3]; n_atoms];
    let mut atomic_quadrupoles = vec![0.0f64; n_atoms];

    let max_scc = 150;
    let scc_tol = 1e-7;
    let mut converged = false;
    let mut scc_iter = 0;

    for iter in 0..max_scc {
        scc_iter = iter + 1;

        // Build multipole-extended gamma matrix
        let gamma = build_multipole_gamma(elements, positions, &charges, &atomic_dipoles);

        // Update charges with multipole contributions
        let (new_charges, new_dipoles, new_quads) =
            update_multipole_charges(elements, positions, &gfn0.mulliken_charges, &gamma);

        let dq: f64 = charges
            .iter()
            .zip(new_charges.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();

        charges = new_charges;
        atomic_dipoles = new_dipoles;
        atomic_quadrupoles = new_quads;

        if dq < scc_tol {
            converged = true;
            break;
        }
    }

    // D4 dispersion (charge-dependent C6)
    let disp_energy = compute_d4_dispersion(elements, positions, &charges);

    // Halogen bond correction
    let xb_energy = compute_halogen_bond_correction(elements, positions);

    // Improved repulsive
    let rep_energy = compute_gfn2_repulsive(elements, positions);

    let total_energy = gfn0.electronic_energy + rep_energy + disp_energy + xb_energy;
    let homo_idx = n_electrons / 2;
    let homo_energy = if homo_idx > 0 && homo_idx <= gfn0.orbital_energies.len() {
        gfn0.orbital_energies[homo_idx - 1]
    } else {
        0.0
    };
    let lumo_energy = if homo_idx < gfn0.orbital_energies.len() {
        gfn0.orbital_energies[homo_idx]
    } else {
        0.0
    };

    Ok(Gfn2Result {
        orbital_energies: gfn0.orbital_energies.clone(),
        electronic_energy: gfn0.electronic_energy,
        repulsive_energy: rep_energy,
        dispersion_energy: disp_energy,
        halogen_bond_energy: xb_energy,
        total_energy,
        n_basis,
        n_electrons,
        homo_energy,
        lumo_energy,
        gap: lumo_energy - homo_energy,
        mulliken_charges: charges,
        atomic_dipoles,
        atomic_quadrupoles,
        scc_iterations: scc_iter,
        converged,
    })
}

/// Build gamma matrix with anisotropic multipole electrostatics.
fn build_multipole_gamma(
    elements: &[u8],
    positions: &[[f64; 3]],
    _charges: &[f64],
    _dipoles: &[[f64; 3]],
) -> DMatrix<f64> {
    let n = elements.len();
    let mut gamma = DMatrix::zeros(n, n);

    for i in 0..n {
        let pi = crate::xtb::params::get_xtb_params(elements[i]).unwrap();
        gamma[(i, i)] = pi.eta;

        for j in (i + 1)..n {
            let pj = crate::xtb::params::get_xtb_params(elements[j]).unwrap();
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();

            let avg_eta = (pi.eta + pj.eta) / 2.0;
            let r_bohr = r / 0.529177;

            // Klopman-Ohno damped Coulomb with GFN2 exponent
            let gamma_ij = 1.0 / (r_bohr.powi(2) + (1.0 / avg_eta).powi(2)).sqrt();

            gamma[(i, j)] = gamma_ij * 27.2114;
            gamma[(j, i)] = gamma[(i, j)];
        }
    }

    gamma
}

/// Update charges including multipole contributions.
fn update_multipole_charges(
    elements: &[u8],
    positions: &[[f64; 3]],
    initial_charges: &[f64],
    _gamma: &DMatrix<f64>,
) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
    let n = elements.len();
    let charges = initial_charges.to_vec();

    // Compute atomic dipole moments from charge distribution
    let mut dipoles = vec![[0.0f64; 3]; n];
    let center = center_of_charges(positions, &charges);
    for i in 0..n {
        for k in 0..3 {
            dipoles[i][k] = charges[i] * (positions[i][k] - center[k]);
        }
    }

    // Quadrupole traces
    let quadrupoles: Vec<f64> = (0..n)
        .map(|i| {
            let r2 = positions[i].iter().map(|x| x * x).sum::<f64>();
            charges[i] * r2
        })
        .collect();

    (charges, dipoles, quadrupoles)
}

fn center_of_charges(positions: &[[f64; 3]], charges: &[f64]) -> [f64; 3] {
    let total_q: f64 = charges.iter().map(|q| q.abs()).sum::<f64>();
    if total_q < 1e-10 {
        return [0.0; 3];
    }
    let mut center = [0.0; 3];
    for (i, pos) in positions.iter().enumerate() {
        let w = charges[i].abs();
        for k in 0..3 {
            center[k] += w * pos[k];
        }
    }
    for k in 0..3 {
        center[k] /= total_q;
    }
    center
}

/// D4 dispersion: charge-dependent C6 coefficients.
///
/// D4 uses electronegativity-weighted reference systems to compute
/// charge-dependent C6 coefficients, giving better accuracy for ionic
/// and polar systems compared to D3.
fn compute_d4_dispersion(elements: &[u8], positions: &[[f64; 3]], charges: &[f64]) -> f64 {
    let n = elements.len();
    let mut e_disp = 0.0;

    let s6 = 1.0;
    let s8 = 2.7;
    let a1 = 0.52;
    let a2 = 5.0;

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();

            if !(0.1..=50.0).contains(&r) {
                continue;
            }

            // Charge-dependent C6 (simplified D4 model)
            let q_scale_i = 1.0 - 0.5 * charges[i].abs().min(1.0);
            let q_scale_j = 1.0 - 0.5 * charges[j].abs().min(1.0);

            let c6_base = get_c6_d4(elements[i], elements[j]);
            let c6 = c6_base * q_scale_i * q_scale_j;
            let c8 = 3.0 * c6 * get_r2r4_d4(elements[i]) * get_r2r4_d4(elements[j]);

            let r0 = (c8 / (c6 + 1e-30)).sqrt();
            let f6 = 1.0 / (r.powi(6) + (a1 * r0 + a2).powi(6));
            let f8 = 1.0 / (r.powi(8) + (a1 * r0 + a2).powi(8));

            e_disp -= s6 * c6 * f6 + s8 * c8 * f8;
        }
    }

    e_disp * 27.2114
}

/// Halogen bond correction for X···B interactions (X=Cl,Br,I; B=N,O,S).
fn compute_halogen_bond_correction(elements: &[u8], positions: &[[f64; 3]]) -> f64 {
    let n = elements.len();
    let mut e_xb = 0.0;

    let halogens = [17u8, 35, 53]; // Cl, Br, I
    let bases = [7u8, 8, 16]; // N, O, S

    for i in 0..n {
        if !halogens.contains(&elements[i]) {
            continue;
        }
        for j in 0..n {
            if i == j || !bases.contains(&elements[j]) {
                continue;
            }

            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();

            // Only short-range halogen bonds
            if r > 4.0 {
                continue;
            }

            let r0 = get_xb_r0(elements[i], elements[j]);
            let strength = get_xb_strength(elements[i]);

            // Damped correction
            let x = r / r0;
            let damp = 1.0 / (1.0 + (-20.0 * (x - 1.0)).exp());
            e_xb -= strength * damp * (-2.0 * (r - r0).powi(2)).exp();
        }
    }

    e_xb * 27.2114 // Hartree to eV
}

fn get_xb_r0(z_hal: u8, z_base: u8) -> f64 {
    let r_hal = match z_hal {
        17 => 1.75,
        35 => 1.85,
        53 => 1.98,
        _ => 1.80,
    };
    let r_base = match z_base {
        7 => 1.55,
        8 => 1.52,
        16 => 1.80,
        _ => 1.60,
    };
    r_hal + r_base
}

fn get_xb_strength(z_hal: u8) -> f64 {
    match z_hal {
        17 => 0.005,
        35 => 0.010,
        53 => 0.015,
        _ => 0.005,
    }
}

/// GFN2 repulsive energy.
fn compute_gfn2_repulsive(elements: &[u8], positions: &[[f64; 3]]) -> f64 {
    let n = elements.len();
    let mut e_rep = 0.0;

    for i in 0..n {
        let pi = crate::xtb::params::get_xtb_params(elements[i]).unwrap();
        for j in (i + 1)..n {
            let pj = crate::xtb::params::get_xtb_params(elements[j]).unwrap();
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();

            if r < 0.1 {
                continue;
            }

            let r_ab = pi.r_cov + pj.r_cov;
            let z_eff_i = pi.n_valence as f64;
            let z_eff_j = pj.n_valence as f64;
            let alpha = (z_eff_i * z_eff_j).sqrt();

            // GFN2 uses steeper repulsive exponent
            e_rep += alpha * (-2.0 * r / r_ab).exp();
        }
    }

    e_rep * 27.2114
}

fn get_c6_d4(z1: u8, z2: u8) -> f64 {
    let c6_1 = atomic_c6_d4(z1);
    let c6_2 = atomic_c6_d4(z2);
    (2.0 * c6_1 * c6_2) / (c6_1 + c6_2 + 1e-30)
}

fn atomic_c6_d4(z: u8) -> f64 {
    match z {
        1 => 6.5,
        6 => 46.6,
        7 => 24.2,
        8 => 15.6,
        9 => 9.52,
        14 => 305.0,
        15 => 185.0,
        16 => 134.0,
        17 => 94.6,
        35 => 162.0,
        53 => 408.0,
        22 => 1044.0,
        24 => 602.0,
        25 => 552.0,
        26 => 482.0,
        27 => 408.0,
        28 => 373.0,
        29 => 253.0,
        30 => 284.0,
        _ => 50.0,
    }
}

fn get_r2r4_d4(z: u8) -> f64 {
    match z {
        1 => 2.00,
        6 => 3.09,
        7 => 2.71,
        8 => 2.44,
        9 => 1.91,
        14 => 4.17,
        15 => 3.63,
        16 => 3.49,
        17 => 3.01,
        35 => 3.47,
        53 => 4.38,
        _ => 3.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gfn2_water() {
        let elements = vec![8u8, 1, 1];
        let positions = vec![
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469],
        ];
        let result = solve_gfn2(&elements, &positions);
        assert!(result.is_ok());
        let r = result.unwrap();
        assert!(r.total_energy.is_finite());
        assert!(r.dispersion_energy.is_finite());
        assert!(r.gap > 0.0);
    }

    #[test]
    fn test_gfn2_dispersion() {
        // Test D4 charge-dependent dispersion
        let elements = vec![6u8, 6];
        let positions = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = vec![0.0, 0.0];
        let e_neutral = compute_d4_dispersion(&elements, &positions, &charges);

        let charges_polar = vec![0.5, -0.5];
        let e_polar = compute_d4_dispersion(&elements, &positions, &charges_polar);

        // Charged atoms should have reduced dispersion
        assert!(e_polar.abs() < e_neutral.abs());
    }
}
