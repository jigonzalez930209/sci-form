//! Frontier Molecular Orbital (FMO) guided approach direction.
//!
//! Uses the spatial distribution of HOMO/LUMO coefficients from EHT to
//! determine the natural approach direction for a reaction — far superior
//! to X-axis forcing.

/// Result of orbital approach analysis.
pub struct OrbitalApproachInfo {
    /// Unit vector for the optimal approach direction.
    pub direction: [f64; 3],
    /// HOMO energy (eV).
    pub homo_energy: f64,
    /// LUMO energy (eV).
    pub lumo_energy: f64,
    /// Which orbital interaction drives this approach (HOMO-LUMO or HOMO-HOMO).
    pub interaction_type: OrbitalInteraction,
}

/// Which frontier orbital interaction is dominant.
#[derive(Debug, Clone, Copy)]
pub enum OrbitalInteraction {
    /// Nucleophile HOMO → Electrophile LUMO
    HomoLumo,
    /// Electrophile HOMO → Nucleophile LUMO (inverted)
    LumoHomo,
    /// Charge-controlled (HOMO-HOMO overlap)
    ChargeControlled,
}

impl std::fmt::Display for OrbitalInteraction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::HomoLumo => write!(f, "HOMO→LUMO"),
            Self::LumoHomo => write!(f, "LUMO→HOMO"),
            Self::ChargeControlled => write!(f, "charge-controlled"),
        }
    }
}

/// Compute the optimal approach direction using EHT orbital information.
///
/// For a bimolecular reaction, the approach should maximise the overlap
/// between the nucleophile's HOMO and the electrophile's LUMO (or vice versa).
///
/// Returns `None` if EHT calculation fails (falls back to electrostatic approach).
pub fn compute_orbital_approach_direction(
    elements_a: &[u8],
    coords_a: &[f64],
    elements_b: &[u8],
    coords_b: &[f64],
) -> Option<OrbitalApproachInfo> {
    // Run EHT on each fragment to get frontier orbitals
    let pos_a: Vec<[f64; 3]> = coords_a.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    let pos_b: Vec<[f64; 3]> = coords_b.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

    let eht_a = crate::eht::solve_eht(elements_a, &pos_a, Some(1.75)).ok()?;
    let eht_b = crate::eht::solve_eht(elements_b, &pos_b, Some(1.75)).ok()?;

    let _n_a = elements_a.len();
    let _n_b = elements_b.len();

    // Find HOMO/LUMO indices
    let n_occ_a = eht_a.n_electrons / 2;
    let n_occ_b = eht_b.n_electrons / 2;

    if n_occ_a == 0 || n_occ_b == 0 {
        return None;
    }

    let homo_a_e = eht_a.homo_energy;
    let lumo_a_e = eht_a.lumo_energy;
    let homo_b_e = eht_b.homo_energy;
    let lumo_b_e = eht_b.lumo_energy;

    // Determine which fragment is nucleophile (higher HOMO) vs electrophile
    let (nuc_homo_e, elec_lumo_e, interaction) =
        if (homo_a_e - lumo_b_e).abs() < (homo_b_e - lumo_a_e).abs() {
            // A is nucleophile (HOMO_A → LUMO_B)
            (homo_a_e, lumo_b_e, OrbitalInteraction::HomoLumo)
        } else {
            // B is nucleophile (HOMO_B → LUMO_A)
            (homo_b_e, lumo_a_e, OrbitalInteraction::LumoHomo)
        };

    // Compute the orbital-weighted centroid for the HOMO of the nucleophile
    // and the LUMO of the electrophile, using atomic contributions
    let n_basis_a = eht_a.coefficients.len();
    let n_basis_b = eht_b.coefficients.len();

    let orbital_centroid_a = compute_orbital_centroid(
        &eht_a.coefficients,
        eht_a.homo_index,
        n_basis_a,
        elements_a,
        coords_a,
    );
    let orbital_centroid_b = compute_orbital_centroid(
        &eht_b.coefficients,
        eht_b.lumo_index,
        n_basis_b,
        elements_b,
        coords_b,
    );

    // Approach direction: from nucleophile's HOMO centroid toward electrophile's LUMO centroid
    let dx = orbital_centroid_b[0] - orbital_centroid_a[0];
    let dy = orbital_centroid_b[1] - orbital_centroid_a[1];
    let dz = orbital_centroid_b[2] - orbital_centroid_a[2];
    let len = (dx * dx + dy * dy + dz * dz).sqrt();

    if len < 1e-10 {
        // Atoms are co-located — use COM direction
        let com_a = com_flat(coords_a);
        let com_b = com_flat(coords_b);
        let dx2 = com_b[0] - com_a[0];
        let dy2 = com_b[1] - com_a[1];
        let dz2 = com_b[2] - com_a[2];
        let l2 = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2).sqrt().max(1e-12);
        return Some(OrbitalApproachInfo {
            direction: [dx2 / l2, dy2 / l2, dz2 / l2],
            homo_energy: nuc_homo_e,
            lumo_energy: elec_lumo_e,
            interaction_type: interaction,
        });
    }

    Some(OrbitalApproachInfo {
        direction: [dx / len, dy / len, dz / len],
        homo_energy: nuc_homo_e,
        lumo_energy: elec_lumo_e,
        interaction_type: interaction,
    })
}

/// Compute the spatial centroid of an orbital weighted by squared coefficients.
///
/// Coefficients are stored as `coefficients[ao_index][mo_index]`.
/// For each basis function, maps it to its parent atom, and accumulates
/// (c²) × atom_position to build a weighted centroid.
fn compute_orbital_centroid(
    coefficients: &[Vec<f64>],
    orbital_idx: usize,
    n_basis: usize,
    elements: &[u8],
    coords: &[f64],
) -> [f64; 3] {
    let n_atoms = elements.len();
    // Build basis → atom map using a simple heuristic
    let mut basis_to_atom = Vec::with_capacity(n_basis);
    for (i, &z) in elements.iter().enumerate() {
        let n_bf = basis_functions_for_element(z);
        for _ in 0..n_bf {
            basis_to_atom.push(i);
        }
    }
    while basis_to_atom.len() < n_basis {
        basis_to_atom.push(n_atoms.saturating_sub(1));
    }

    let mut centroid = [0.0f64; 3];
    let mut total_wt = 0.0f64;

    for mu in 0..n_basis.min(basis_to_atom.len()) {
        if mu >= coefficients.len() || orbital_idx >= coefficients[mu].len() {
            continue;
        }
        let c = coefficients[mu][orbital_idx];
        let w = c * c;
        let atom = basis_to_atom[mu];
        if atom < n_atoms {
            centroid[0] += w * coords[atom * 3];
            centroid[1] += w * coords[atom * 3 + 1];
            centroid[2] += w * coords[atom * 3 + 2];
            total_wt += w;
        }
    }

    if total_wt > 1e-14 {
        centroid[0] /= total_wt;
        centroid[1] /= total_wt;
        centroid[2] /= total_wt;
    }

    centroid
}

/// Estimate number of minimal basis functions for an element.
fn basis_functions_for_element(z: u8) -> usize {
    match z {
        1 => 1,        // H: 1s
        2 => 1,        // He: 1s
        3..=4 => 5,    // Li-Be: 1s + 2s + 2p
        5..=10 => 5,   // B-Ne: 1s + 2s + 2p
        11..=18 => 9,  // Na-Ar: + 3s + 3p
        19..=36 => 13, // K-Kr: + 3d + 4s
        37..=54 => 18, // Rb-Xe: + 4d + 5s + 4p
        _ => 18,       // heavy elements: estimate
    }
}

fn com_flat(coords: &[f64]) -> [f64; 3] {
    let n = coords.len() / 3;
    if n == 0 {
        return [0.0; 3];
    }
    let mut c = [0.0; 3];
    for i in 0..n {
        c[0] += coords[i * 3];
        c[1] += coords[i * 3 + 1];
        c[2] += coords[i * 3 + 2];
    }
    let nf = n as f64;
    [c[0] / nf, c[1] / nf, c[2] / nf]
}
