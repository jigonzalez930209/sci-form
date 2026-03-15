//! Shrake-Rupley solvent-accessible surface area (SASA) algorithm.
//!
//! Distributes test points on a sphere around each atom (van der Waals radius
//! plus probe radius), counts how many points are NOT buried inside a neighbor,
//! and scales to area.
//!
//! Reference: A. Shrake & J.A. Rupley, *J. Mol. Biol.* **79**, 351–371, 1973.

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Default probe radius (water) in Ångström.
pub const DEFAULT_PROBE_RADIUS: f64 = 1.4;

/// Default number of test points per atom sphere (92 is a common choice,
/// 960 or 2562 give higher accuracy).
pub const DEFAULT_NUM_POINTS: usize = 960;

/// Van der Waals radii (Å) by atomic number (Bondi radii).
pub fn vdw_radius(z: u8) -> f64 {
    match z {
        1 => 1.20,
        5 => 1.92,
        6 => 1.70,
        7 => 1.55,
        8 => 1.52,
        9 => 1.47,
        14 => 2.10,
        15 => 1.80,
        16 => 1.80,
        17 => 1.75,
        35 => 1.85,
        53 => 1.98,
        34 => 1.90,
        _ => 1.70, // default fallback
    }
}

/// Result of a SASA calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SasaResult {
    /// Total solvent-accessible surface area in Ų.
    pub total_sasa: f64,
    /// Per-atom SASA values in Ų.
    pub atom_sasa: Vec<f64>,
    /// Probe radius used (Å).
    pub probe_radius: f64,
    /// Number of test points per sphere.
    pub num_points: usize,
}

/// Generate approximately uniformly distributed points on a unit sphere using
/// a golden-section spiral (Fibonacci sphere).
fn generate_sphere_points(n: usize) -> Vec<[f64; 3]> {
    let golden_ratio = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let angle_increment = 2.0 * PI / golden_ratio;

    (0..n)
        .map(|i| {
            let y = 1.0 - (2.0 * i as f64) / (n as f64 - 1.0);
            let radius = (1.0 - y * y).max(0.0).sqrt();
            let theta = angle_increment * i as f64;
            [radius * theta.cos(), y, radius * theta.sin()]
        })
        .collect()
}

fn compute_atom_sasa(
    atom_index: usize,
    positions: &[[f64; 3]],
    radii: &[f64],
    unit_points: &[[f64; 3]],
) -> f64 {
    let ri = radii[atom_index];
    let pi = &positions[atom_index];
    let npts = unit_points.len();
    let area_per_point = 4.0 * PI * ri * ri / npts as f64;

    let accessible = unit_points
        .iter()
        .filter(|pt| {
            let test = [pi[0] + ri * pt[0], pi[1] + ri * pt[1], pi[2] + ri * pt[2]];

            !positions.iter().enumerate().any(|(j, pos_j)| {
                if j == atom_index {
                    return false;
                }
                let rj = radii[j];
                let dx = test[0] - pos_j[0];
                let dy = test[1] - pos_j[1];
                let dz = test[2] - pos_j[2];
                let dist_sq = dx * dx + dy * dy + dz * dz;
                dist_sq < rj * rj
            })
        })
        .count();

    accessible as f64 * area_per_point
}

/// Compute solvent-accessible surface area for a set of atoms.
///
/// # Arguments
/// - `elements`: atomic numbers for each atom
/// - `positions`: 3D coordinates in Å, one per atom
/// - `probe_radius`: solvent probe radius in Å (default: 1.4 for water)
/// - `num_points`: number of test points per atom sphere (default: 960)
///
/// # Returns
/// `SasaResult` with total and per-atom SASA.
pub fn compute_sasa(
    elements: &[u8],
    positions: &[[f64; 3]],
    probe_radius: Option<f64>,
    num_points: Option<usize>,
) -> SasaResult {
    let n = elements.len();
    let probe = probe_radius.unwrap_or(DEFAULT_PROBE_RADIUS);
    let npts = num_points.unwrap_or(DEFAULT_NUM_POINTS);

    let unit_points = generate_sphere_points(npts);

    // Precompute radii (vdw + probe)
    let radii: Vec<f64> = elements.iter().map(|&z| vdw_radius(z) + probe).collect();

    let atom_sasa: Vec<f64> = (0..n)
        .map(|i| compute_atom_sasa(i, positions, &radii, &unit_points))
        .collect();

    let total_sasa: f64 = atom_sasa.iter().sum();

    SasaResult {
        total_sasa,
        atom_sasa,
        probe_radius: probe,
        num_points: npts,
    }
}

/// Compute SASA using rayon to evaluate atoms independently.
#[cfg(feature = "parallel")]
pub fn compute_sasa_parallel(
    elements: &[u8],
    positions: &[[f64; 3]],
    probe_radius: Option<f64>,
    num_points: Option<usize>,
) -> SasaResult {
    use rayon::prelude::*;

    let probe = probe_radius.unwrap_or(DEFAULT_PROBE_RADIUS);
    let npts = num_points.unwrap_or(DEFAULT_NUM_POINTS);
    let unit_points = generate_sphere_points(npts);
    let radii: Vec<f64> = elements.iter().map(|&z| vdw_radius(z) + probe).collect();

    let atom_sasa: Vec<f64> = (0..elements.len())
        .into_par_iter()
        .map(|i| compute_atom_sasa(i, positions, &radii, &unit_points))
        .collect();

    let total_sasa: f64 = atom_sasa.iter().sum();

    SasaResult {
        total_sasa,
        atom_sasa,
        probe_radius: probe,
        num_points: npts,
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_atom_sasa() {
        // Single carbon atom: SASA should be approximately 4π(r_vdw + r_probe)²
        let elems = vec![6];
        let pos = vec![[0.0, 0.0, 0.0]];
        let result = compute_sasa(&elems, &pos, Some(1.4), Some(2000));
        let r = vdw_radius(6) + 1.4;
        let expected = 4.0 * PI * r * r;
        let error = (result.total_sasa - expected).abs() / expected;
        assert!(
            error < 0.02,
            "Single atom SASA should be ~{:.1} Ų, got {:.1} (error {:.1}%)",
            expected,
            result.total_sasa,
            error * 100.0
        );
    }

    #[test]
    fn test_two_distant_atoms() {
        // Two atoms far apart: total SASA ≈ sum of individual SASAs
        let elems = vec![6, 6];
        let pos = vec![[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]]; // very far apart
        let result = compute_sasa(&elems, &pos, Some(1.4), Some(2000));
        let r = vdw_radius(6) + 1.4;
        let expected_per_atom = 4.0 * PI * r * r;
        let expected_total = 2.0 * expected_per_atom;
        let error = (result.total_sasa - expected_total).abs() / expected_total;
        assert!(
            error < 0.02,
            "Far-apart atoms: expected {:.1}, got {:.1}",
            expected_total,
            result.total_sasa
        );
    }

    #[test]
    fn test_bonded_atoms_less_sasa() {
        // Two carbon atoms at bonding distance: SASA < 2× individual
        let elems = vec![6, 6];
        let pos_far = vec![[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]];
        let pos_close = vec![[0.0, 0.0, 0.0], [1.54, 0.0, 0.0]]; // C-C bond ~1.54 Å
        let far_result = compute_sasa(&elems, &pos_far, Some(1.4), Some(960));
        let close_result = compute_sasa(&elems, &pos_close, Some(1.4), Some(960));
        assert!(
            close_result.total_sasa < far_result.total_sasa,
            "Bonded atoms should have less SASA ({:.1}) than distant ({:.1})",
            close_result.total_sasa,
            far_result.total_sasa
        );
    }

    #[test]
    fn test_water_sasa() {
        // H₂O: O at origin, 2 H at typical geometry
        let elems = vec![8, 1, 1];
        let pos = vec![[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let result = compute_sasa(&elems, &pos, Some(1.4), Some(960));
        // SASA for water should be reasonable (roughly 50–150 Ų depending on probe)
        assert!(
            result.total_sasa > 30.0 && result.total_sasa < 200.0,
            "Water SASA should be reasonable: {:.1}",
            result.total_sasa
        );
        assert_eq!(result.atom_sasa.len(), 3);
    }

    #[test]
    fn test_zero_probe_radius() {
        // With zero probe: gives van der Waals surface area
        let elems = vec![6];
        let pos = vec![[0.0, 0.0, 0.0]];
        let result = compute_sasa(&elems, &pos, Some(0.0), Some(2000));
        let r = vdw_radius(6);
        let expected = 4.0 * PI * r * r;
        let error = (result.total_sasa - expected).abs() / expected;
        assert!(
            error < 0.02,
            "Zero-probe SASA should be vdW surface area: expected {:.1}, got {:.1}",
            expected,
            result.total_sasa
        );
    }

    #[test]
    fn test_sphere_points_distribution() {
        // Points should be approximately uniformly distributed
        let points = generate_sphere_points(1000);
        assert_eq!(points.len(), 1000);
        // Check all points are on the unit sphere
        for p in &points {
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!(
                (r - 1.0).abs() < 1e-10,
                "Point should be on unit sphere: r={}",
                r
            );
        }
    }

    #[test]
    fn test_sasa_symmetry() {
        // Two identical atoms equidistant from origin: should have same SASA
        let elems = vec![6, 6];
        let pos = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let result = compute_sasa(&elems, &pos, None, None);
        assert!(
            (result.atom_sasa[0] - result.atom_sasa[1]).abs() < 0.5,
            "Symmetric atoms should have similar SASA: {:.1} vs {:.1}",
            result.atom_sasa[0],
            result.atom_sasa[1]
        );
    }
}
