//! Behler-Parrinello Atomic Environment Vectors (AEVs).
//!
//! Transforms Cartesian coordinates into rotation/translation-invariant
//! descriptors using radial and angular symmetry functions.
//!
//! Radial: $G_R = \sum_{j} e^{-\eta(R_{ij}-R_s)^2} f_c(R_{ij})$
//! Angular: $G_A = 2^{1-\zeta} \sum_{j<k} (1+\cos(\theta-\theta_s))^\zeta
//!          \cdot e^{-\eta((R_{ij}+R_{ik})/2 - R_s)^2} f_c(R_{ij}) f_c(R_{ik})$

use super::aev_params::{species_index, AevParams, N_SPECIES};
use super::cutoff::cosine_cutoff;
use super::neighbor::NeighborPair;

/// Compute AEVs for all atoms in the system.
///
/// Returns a `Vec` of length `n_atoms`, each entry is an AEV vector.
pub fn compute_aevs(
    elements: &[u8],
    positions: &[[f64; 3]],
    neighbors: &[NeighborPair],
    params: &AevParams,
) -> Vec<Vec<f64>> {
    let n = elements.len();
    let aev_len = params.total_aev_length();

    // Build per-atom neighbor lists (both directions)
    let mut atom_neighbors: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for np in neighbors {
        let d = np.dist_sq.sqrt();
        atom_neighbors[np.i].push((np.j, d));
        atom_neighbors[np.j].push((np.i, d));
    }

    let mut aevs = vec![vec![0.0f64; aev_len]; n];

    for i in 0..n {
        let si = match species_index(elements[i]) {
            Some(s) => s,
            None => continue,
        };
        let _ = si; // atom i species (used for symmetry)
        compute_radial_aev(
            i,
            elements,
            positions,
            &atom_neighbors[i],
            params,
            &mut aevs[i],
        );
        compute_angular_aev(
            i,
            elements,
            positions,
            &atom_neighbors[i],
            params,
            &mut aevs[i],
        );
    }

    aevs
}

fn compute_radial_aev(
    _i: usize,
    elements: &[u8],
    _positions: &[[f64; 3]],
    neighbors_i: &[(usize, f64)],
    params: &AevParams,
    aev: &mut [f64],
) {
    let rad_len = params.radial_length();

    for &(j, rij) in neighbors_i {
        if rij >= params.radial_cutoff {
            continue;
        }
        let sj = match species_index(elements[j]) {
            Some(s) => s,
            None => continue,
        };
        let fc = cosine_cutoff(rij, params.radial_cutoff);
        let offset = sj * rad_len;

        let mut k = 0;
        for eta in &params.radial_eta {
            for rs in &params.radial_rs {
                let dr = rij - rs;
                aev[offset + k] += (-eta * dr * dr).exp() * fc;
                k += 1;
            }
        }
    }
}

fn compute_angular_aev(
    i: usize,
    elements: &[u8],
    positions: &[[f64; 3]],
    neighbors_i: &[(usize, f64)],
    params: &AevParams,
    aev: &mut [f64],
) {
    let rad_total = N_SPECIES * params.radial_length();
    let ang_len = params.angular_length();

    // Filter neighbors within angular cutoff
    let ang_neighbors: Vec<(usize, f64)> = neighbors_i
        .iter()
        .filter(|&&(_, d)| d < params.angular_cutoff)
        .copied()
        .collect();

    for a in 0..ang_neighbors.len() {
        let (j, rij) = ang_neighbors[a];
        let sj = match species_index(elements[j]) {
            Some(s) => s,
            None => continue,
        };
        let fc_ij = cosine_cutoff(rij, params.angular_cutoff);

        for b in (a + 1)..ang_neighbors.len() {
            let (k, rik) = ang_neighbors[b];
            let sk = match species_index(elements[k]) {
                Some(s) => s,
                None => continue,
            };
            let fc_ik = cosine_cutoff(rik, params.angular_cutoff);

            let theta = compute_angle(positions, i, j, k);
            let (s_lo, s_hi) = if sj <= sk { (sj, sk) } else { (sk, sj) };
            let pair_idx = s_lo * (2 * N_SPECIES - s_lo - 1) / 2 + (s_hi - s_lo);
            let offset = rad_total + pair_idx * ang_len;

            let r_avg = (rij + rik) / 2.0;
            let mut m = 0;
            for eta in &params.angular_eta {
                for rs in &params.angular_rs {
                    for zeta in &params.angular_zeta {
                        for theta_s in &params.angular_theta_s {
                            let cos_term = 1.0 + (theta - theta_s).cos();
                            let angular = 2.0f64.powf(1.0 - zeta) * cos_term.powf(*zeta);
                            let radial = (-eta * (r_avg - rs).powi(2)).exp();
                            aev[offset + m] += angular * radial * fc_ij * fc_ik;
                            m += 1;
                        }
                    }
                }
            }
        }
    }
}

/// Compute angle ∠jik from positions.
fn compute_angle(positions: &[[f64; 3]], i: usize, j: usize, k: usize) -> f64 {
    let vij = [
        positions[j][0] - positions[i][0],
        positions[j][1] - positions[i][1],
        positions[j][2] - positions[i][2],
    ];
    let vik = [
        positions[k][0] - positions[i][0],
        positions[k][1] - positions[i][1],
        positions[k][2] - positions[i][2],
    ];
    let dot = vij[0] * vik[0] + vij[1] * vik[1] + vij[2] * vik[2];
    let nij = (vij[0] * vij[0] + vij[1] * vij[1] + vij[2] * vij[2]).sqrt();
    let nik = (vik[0] * vik[0] + vik[1] * vik[1] + vik[2] * vik[2]).sqrt();
    let cos_theta = (dot / (nij * nik)).clamp(-1.0, 1.0);
    cos_theta.acos()
}

#[cfg(test)]
mod tests {
    use super::super::aev_params::default_ani2x_params;
    use super::super::neighbor::CellList;
    use super::*;

    #[test]
    fn test_aev_water() {
        let elements = [8u8, 1, 1];
        let positions = [
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469],
        ];
        let params = default_ani2x_params();
        let cl = CellList::new(&positions, params.radial_cutoff);
        let neighbors = cl.find_neighbors(&positions);
        let aevs = compute_aevs(&elements, &positions, &neighbors, &params);

        assert_eq!(aevs.len(), 3);
        // Both H atoms should have symmetric AEVs
        let diff: f64 = aevs[1]
            .iter()
            .zip(aevs[2].iter())
            .map(|(a, b)| (a - b).abs())
            .sum();
        assert!(diff < 1e-10, "H atoms in water should have symmetric AEVs");
    }

    #[test]
    fn test_aev_nonzero() {
        let elements = [6u8, 1, 1, 1, 1];
        let positions = [
            [0.0, 0.0, 0.0],
            [0.63, 0.63, 0.63],
            [-0.63, -0.63, 0.63],
            [-0.63, 0.63, -0.63],
            [0.63, -0.63, -0.63],
        ];
        let params = default_ani2x_params();
        let cl = CellList::new(&positions, params.radial_cutoff);
        let neighbors = cl.find_neighbors(&positions);
        let aevs = compute_aevs(&elements, &positions, &neighbors, &params);

        // Carbon AEV should have nonzero entries
        let sum: f64 = aevs[0].iter().map(|v| v.abs()).sum();
        assert!(sum > 0.0, "Carbon AEV should be nonzero");
    }
}
