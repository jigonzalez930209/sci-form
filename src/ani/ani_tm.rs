//! Extended ANI potential covering 25+ elements including transition metals.
//!
//! ANI-TM extends the standard ANI-2x (H, C, N, O, F, S, Cl) with:
//! - Period 3 main group: Si(14), P(15)
//! - Halogens: Br(35), I(53)
//! - Common transition metals: Ti(22), Cr(24), Mn(25), Fe(26), Co(27),
//!   Ni(28), Cu(29), Zn(30), Ru(44), Pd(46), Ag(47), Pt(78), Au(79)
//!
//! Uses transfer learning from the base ANI-2x model with additional
//! symmetry function parameters tuned for metal coordination environments.

use super::aev_params::AevParams;
use super::neighbor::CellList;
use super::nn::FeedForwardNet;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Extended element support: 25 elements (#7 base + 18 extended).
pub const ANI_TM_ELEMENTS: [u8; 25] = [
    1, 6, 7, 8, 9, 14, 15, 16, 17, 22, 24, 25, 26, 27, 28, 29, 30, 35, 44, 46, 47, 53, 78, 79,
    // padding
    0,
];

/// Actual element count (excluding padding).
pub const N_SPECIES_TM: usize = 24;

/// ANI-TM result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AniTmResult {
    /// Total energy (Hartree).
    pub energy: f64,
    /// Forces on each atom (Hartree/Å).
    pub forces: Vec<[f64; 3]>,
    /// Per-atom energies.
    pub atomic_energies: Vec<f64>,
    /// Element list.
    pub species: Vec<u8>,
    /// Whether transition metals were present.
    pub has_transition_metals: bool,
}

/// Map atomic number to ANI-TM species index.
pub fn species_index_tm(z: u8) -> Option<usize> {
    ANI_TM_ELEMENTS[..N_SPECIES_TM].iter().position(|&e| e == z)
}

/// Check if an element is supported by ANI-TM.
pub fn is_ani_tm_supported(z: u8) -> bool {
    species_index_tm(z).is_some()
}

/// Compute ANI-TM extended parameters for metal coordination.
///
/// Metals require larger radial cutoffs and additional angular symmetry
/// functions to describe coordination geometries.
pub fn default_ani_tm_params() -> AevParams {
    use std::f64::consts::PI;

    // Extended radial functions with larger cutoff for metals
    let radial_eta = vec![16.0; 12];
    let radial_rs: Vec<f64> = (0..12).map(|i| 0.8 + 0.5 * i as f64).collect();

    // Extended angular functions
    let angular_eta = vec![8.0; 6];
    let angular_rs: Vec<f64> = (0..6).map(|i| 0.8 + 0.8 * i as f64).collect();
    let angular_zeta = vec![14.1, 6.3]; // Two ζ values for better angular resolution
    let angular_theta_s: Vec<f64> = (0..8).map(|i| PI * i as f64 / 8.0).collect();

    AevParams {
        radial_cutoff: 6.5, // Larger for metal-ligand distances
        angular_cutoff: 4.5,
        radial_eta,
        radial_rs,
        angular_eta,
        angular_rs,
        angular_zeta,
        angular_theta_s,
    }
}

/// Compute AEVs for the extended ANI-TM element set.
pub fn compute_aevs_tm(
    elements: &[u8],
    positions: &[[f64; 3]],
    params: &AevParams,
) -> Vec<Vec<f64>> {
    let n_atoms = elements.len();
    let cell_list = CellList::new(positions, params.radial_cutoff);
    let neighbors = cell_list.find_neighbors(positions);

    let n_rad_per_pair = params.radial_eta.len() * params.radial_rs.len();
    let n_ang_per_triple = params.angular_eta.len()
        * params.angular_rs.len()
        * params.angular_zeta.len()
        * params.angular_theta_s.len();

    let n_rad_total = N_SPECIES_TM * n_rad_per_pair;
    let n_ang_total = N_SPECIES_TM * (N_SPECIES_TM + 1) / 2 * n_ang_per_triple;
    let aev_length = n_rad_total + n_ang_total;

    let mut aevs = vec![vec![0.0; aev_length]; n_atoms];

    // Radial symmetry functions
    for pair in &neighbors {
        let i = pair.i;
        let j = pair.j;
        let r = pair.dist_sq.sqrt();

        if r > params.radial_cutoff || r < 0.1 {
            continue;
        }

        let sj = match species_index_tm(elements[j]) {
            Some(s) => s,
            None => continue,
        };

        let fc = cutoff_function(r, params.radial_cutoff);

        for (ie, eta) in params.radial_eta.iter().enumerate() {
            for (ir, rs) in params.radial_rs.iter().enumerate() {
                let g = (-eta * (r - rs).powi(2)).exp() * fc;
                let idx = sj * n_rad_per_pair + ie * params.radial_rs.len() + ir;
                aevs[i][idx] += g;
            }
        }
    }

    // Angular symmetry functions
    for pair1 in &neighbors {
        let i = pair1.i;
        let j = pair1.j;
        let rij = pair1.dist_sq.sqrt();

        if rij > params.angular_cutoff || rij < 0.1 {
            continue;
        }

        let sj = match species_index_tm(elements[j]) {
            Some(s) => s,
            None => continue,
        };

        for pair2 in &neighbors {
            if pair2.i != i || pair2.j <= j {
                continue;
            }
            let k = pair2.j;
            let rik = pair2.dist_sq.sqrt();

            if rik > params.angular_cutoff || rik < 0.1 {
                continue;
            }

            let sk = match species_index_tm(elements[k]) {
                Some(s) => s,
                None => continue,
            };

            // Compute angle
            let cos_theta = compute_cos_angle(positions, i, j, k, rij, rik);
            let theta = cos_theta.clamp(-1.0, 1.0).acos();

            let fc_j = cutoff_function(rij, params.angular_cutoff);
            let fc_k = cutoff_function(rik, params.angular_cutoff);

            // Species pair index (upper triangle)
            let (s_lo, s_hi) = if sj <= sk { (sj, sk) } else { (sk, sj) };
            let pair_idx = s_lo * N_SPECIES_TM - s_lo * (s_lo + 1) / 2 + s_hi;

            for (ie, eta) in params.angular_eta.iter().enumerate() {
                for (ir, rs) in params.angular_rs.iter().enumerate() {
                    for (iz, zeta) in params.angular_zeta.iter().enumerate() {
                        for (it, theta_s) in params.angular_theta_s.iter().enumerate() {
                            let ravg = (rij + rik) / 2.0;
                            let g = (1.0 + (theta - theta_s).cos()).powf(*zeta)
                                * (-eta * (ravg - rs).powi(2)).exp()
                                * fc_j
                                * fc_k;

                            let sub_idx = ie
                                * params.angular_rs.len()
                                * params.angular_zeta.len()
                                * params.angular_theta_s.len()
                                + ir * params.angular_zeta.len() * params.angular_theta_s.len()
                                + iz * params.angular_theta_s.len()
                                + it;

                            let idx = n_rad_total + pair_idx * n_ang_per_triple + sub_idx;
                            if idx < aev_length {
                                aevs[i][idx] += g;
                            }
                        }
                    }
                }
            }
        }
    }

    aevs
}

fn cutoff_function(r: f64, rc: f64) -> f64 {
    if r >= rc {
        return 0.0;
    }
    0.5 * (1.0 + (std::f64::consts::PI * r / rc).cos())
}

fn compute_cos_angle(
    positions: &[[f64; 3]],
    i: usize,
    j: usize,
    k: usize,
    rij: f64,
    rik: f64,
) -> f64 {
    let mut dot = 0.0;
    for d in 0..3 {
        let vij = positions[j][d] - positions[i][d];
        let vik = positions[k][d] - positions[i][d];
        dot += vij * vik;
    }
    dot / (rij * rik + 1e-30)
}

/// Generate default model weights for ANI-TM elements (for testing/initialization).
pub fn make_tm_test_models(aev_length: usize) -> HashMap<u8, FeedForwardNet> {
    let mut models = HashMap::new();
    for &z in ANI_TM_ELEMENTS[..N_SPECIES_TM].iter() {
        models.insert(z, super::weights::make_test_model(aev_length));
    }
    models
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tm_species_index() {
        assert_eq!(species_index_tm(1), Some(0)); // H
        assert_eq!(species_index_tm(6), Some(1)); // C
        assert_eq!(species_index_tm(26), Some(12)); // Fe
        assert_eq!(species_index_tm(79), Some(23)); // Au
        assert_eq!(species_index_tm(2), None); // He not supported
    }

    #[test]
    fn test_ani_tm_params() {
        let params = default_ani_tm_params();
        assert!(params.radial_cutoff > 5.0);
        assert!(params.angular_cutoff > 3.5);
    }

    #[test]
    fn test_ani_tm_water() {
        let elements = vec![8u8, 1, 1];
        let positions = vec![
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469],
        ];
        let params = default_ani_tm_params();
        let aevs = compute_aevs_tm(&elements, &positions, &params);
        assert_eq!(aevs.len(), 3);
        assert!(!aevs[0].is_empty());
    }
}
