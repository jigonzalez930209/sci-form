//! **ALPHA** — Generalized Behler-Parrinello symmetry functions.
//!
//! Extracts and generalizes the AEV computation from `ani/aev.rs` into a
//! standalone descriptor module with configurable radial/angular parameters.

use serde::{Deserialize, Serialize};

/// Symmetry function parameter set.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SymmetryFunctionParams {
    /// Radial cutoff (Å).
    pub radial_cutoff: f64,
    /// Angular cutoff (Å).
    pub angular_cutoff: f64,
    /// Radial shift values (Å).
    pub radial_shifts: Vec<f64>,
    /// Radial eta values (Å⁻²).
    pub radial_etas: Vec<f64>,
    /// Angular shift values (Å).
    pub angular_shifts: Vec<f64>,
    /// Angular eta values (Å⁻²).
    pub angular_etas: Vec<f64>,
    /// Angular zeta values.
    pub angular_zetas: Vec<f64>,
}

impl Default for SymmetryFunctionParams {
    fn default() -> Self {
        // ANI-2x compatible defaults
        let radial_shifts: Vec<f64> = (0..16).map(|i| 0.9 + 0.25 * i as f64).collect();
        let angular_shifts: Vec<f64> = (0..4).map(|i| 0.9 + 0.5 * i as f64).collect();

        Self {
            radial_cutoff: 5.2,
            angular_cutoff: 3.5,
            radial_shifts,
            radial_etas: vec![16.0; 16],
            angular_shifts,
            angular_etas: vec![8.0; 4],
            angular_zetas: vec![1.0, 2.0, 4.0, 8.0],
        }
    }
}

/// Atomic Environment Vector (AEV) for a single atom.
#[derive(Debug, Clone)]
pub struct Aev {
    /// Radial symmetry function values.
    pub radial: Vec<f64>,
    /// Angular symmetry function values.
    pub angular: Vec<f64>,
}

impl Aev {
    /// Total descriptor length.
    pub fn len(&self) -> usize {
        self.radial.len() + self.angular.len()
    }

    /// Whether the AEV is empty.
    pub fn is_empty(&self) -> bool {
        self.radial.is_empty() && self.angular.is_empty()
    }

    /// Concatenate into a single vector.
    pub fn to_vec(&self) -> Vec<f64> {
        let mut v = self.radial.clone();
        v.extend_from_slice(&self.angular);
        v
    }
}

/// Cutoff function f_c(r) = 0.5 * (cos(π r / R_c) + 1) for r < R_c.
#[inline]
fn cutoff_fn(r: f64, rc: f64) -> f64 {
    if r >= rc {
        0.0
    } else {
        0.5 * (1.0 + (std::f64::consts::PI * r / rc).cos())
    }
}

/// Compute AEVs for all atoms.
///
/// With the `parallel` feature enabled, the outer atom loop runs via rayon,
/// making AEV computation scale with available CPU threads.
pub fn compute_aevs(
    elements: &[u8],
    positions: &[[f64; 3]],
    params: &SymmetryFunctionParams,
) -> Vec<Aev> {
    let n = elements.len();
    let n_radial = params.radial_shifts.len() * params.radial_etas.len();
    let n_angular =
        params.angular_shifts.len() * params.angular_etas.len() * params.angular_zetas.len();

    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        (0..n)
            .into_par_iter()
            .map(|i| compute_single_aev(i, n, positions, params, n_radial, n_angular))
            .collect()
    }

    #[cfg(not(feature = "parallel"))]
    {
        (0..n)
            .map(|i| compute_single_aev(i, n, positions, params, n_radial, n_angular))
            .collect()
    }
}

fn compute_single_aev(
    i: usize,
    n: usize,
    positions: &[[f64; 3]],
    params: &SymmetryFunctionParams,
    n_radial: usize,
    n_angular: usize,
) -> Aev {
    let mut radial = vec![0.0; n_radial];
    let mut angular = vec![0.0; n_angular];

    // Radial terms
    for j in 0..n {
        if i == j {
            continue;
        }
        let r_ij = distance(&positions[i], &positions[j]);
        if r_ij >= params.radial_cutoff {
            continue;
        }
        let fc = cutoff_fn(r_ij, params.radial_cutoff);

        for (s_idx, (&rs, &eta)) in params
            .radial_shifts
            .iter()
            .zip(params.radial_etas.iter())
            .enumerate()
        {
            radial[s_idx] += (-eta * (r_ij - rs).powi(2)).exp() * fc;
        }
    }

    // Angular terms
    for j in 0..n {
        if i == j {
            continue;
        }
        let r_ij = distance(&positions[i], &positions[j]);
        if r_ij >= params.angular_cutoff {
            continue;
        }

        for k in (j + 1)..n {
            if i == k {
                continue;
            }
            let r_ik = distance(&positions[i], &positions[k]);
            if r_ik >= params.angular_cutoff {
                continue;
            }

            let theta = angle(&positions[i], &positions[j], &positions[k]);
            let fc_ij = cutoff_fn(r_ij, params.angular_cutoff);
            let fc_ik = cutoff_fn(r_ik, params.angular_cutoff);

            let mut idx = 0;
            for &rs in &params.angular_shifts {
                for &eta in &params.angular_etas {
                    for &zeta in &params.angular_zetas {
                        let r_mid = 0.5 * (r_ij + r_ik);
                        let g = 2.0f64.powf(1.0 - zeta)
                            * (1.0 + theta.cos()).powf(zeta)
                            * (-eta * (r_mid - rs).powi(2)).exp()
                            * fc_ij
                            * fc_ik;
                        angular[idx] += g;
                        idx += 1;
                    }
                }
            }
        }
    }

    Aev { radial, angular }
}

fn distance(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn angle(center: &[f64; 3], a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let va = [a[0] - center[0], a[1] - center[1], a[2] - center[2]];
    let vb = [b[0] - center[0], b[1] - center[1], b[2] - center[2]];
    let dot = va[0] * vb[0] + va[1] * vb[1] + va[2] * vb[2];
    let ma = (va[0] * va[0] + va[1] * va[1] + va[2] * va[2]).sqrt();
    let mb = (vb[0] * vb[0] + vb[1] * vb[1] + vb[2] * vb[2]).sqrt();
    (dot / (ma * mb + 1e-30)).clamp(-1.0, 1.0).acos()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aev_default_params() {
        let p = SymmetryFunctionParams::default();
        assert!(!p.radial_etas.is_empty());
        assert!(!p.radial_shifts.is_empty());
        assert!(p.radial_cutoff > 0.0);
    }

    #[test]
    fn aev_single_atom_is_zero() {
        let params = SymmetryFunctionParams::default();
        let aevs = compute_aevs(&[1], &[[0.0, 0.0, 0.0]], &params);
        assert_eq!(aevs.len(), 1);
        // Lone atom has no neighbors — all AEV components should be zero
        for v in &aevs[0].radial {
            assert!(v.abs() < 1e-15, "lone atom radial AEV should be zero");
        }
    }

    #[test]
    fn aev_h2_has_nonzero_radial() {
        let params = SymmetryFunctionParams::default();
        let aevs = compute_aevs(&[1, 1], &[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]], &params);
        assert_eq!(aevs.len(), 2);
        let radial_sum: f64 = aevs[0].radial.iter().sum();
        assert!(radial_sum > 0.0, "H2 should have nonzero radial AEV");
    }

    #[test]
    fn aev_symmetric_molecule_equal_descriptors() {
        let params = SymmetryFunctionParams::default();
        let aevs = compute_aevs(&[1, 1], &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], &params);
        let v0 = aevs[0].to_vec();
        let v1 = aevs[1].to_vec();
        assert_eq!(v0.len(), v1.len());
        for i in 0..v0.len() {
            assert!(
                (v0[i] - v1[i]).abs() < 1e-12,
                "symmetric atoms should have equal AEVs at index {i}"
            );
        }
    }

    #[test]
    fn aev_length_consistent() {
        let params = SymmetryFunctionParams::default();
        let aevs = compute_aevs(&[1, 6], &[[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], &params);
        assert!(!aevs[0].is_empty());
        assert_eq!(aevs[0].len(), aevs[0].radial.len() + aevs[0].angular.len());
    }
}
