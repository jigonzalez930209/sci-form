//! RDF (Radial Distribution Function) molecular descriptors.
//!
//! These descriptors encode 3D atomic distance distributions as
//! property-weighted functions sampled at discrete radii.
//! Weighted by atomic properties (charges, masses, electronegativities)
//! they capture shape and electronic distribution information.
//!
//! Reference: Hemmer et al., JCICS 39, 1076–1084 (1999).

use serde::{Deserialize, Serialize};

/// RDF descriptor result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RdfDescriptors {
    /// Sampling radii (Å).
    pub radii: Vec<f64>,
    /// Unweighted RDF values at each radius.
    pub rdf_unweighted: Vec<f64>,
    /// Mass-weighted RDF.
    pub rdf_mass: Vec<f64>,
    /// Electronegativity-weighted RDF.
    pub rdf_electronegativity: Vec<f64>,
    /// Charge-weighted RDF (if charges provided).
    pub rdf_charge: Vec<f64>,
    /// Smoothing factor used (β).
    pub beta: f64,
    /// Number of atom pairs contributing.
    pub n_pairs: usize,
}

/// Compute RDF descriptors from 3D coordinates.
///
/// # Arguments
/// * `elements` - Atomic numbers
/// * `positions` - 3D coordinates
/// * `charges` - Optional partial charges (empty = skip charge-weighted RDF)
/// * `beta` - Gaussian smoothing factor (default: 100 Å⁻²)
/// * `r_max` - Maximum radius (default: 12.0 Å)
/// * `n_points` - Number of sampling points (default: 128)
pub fn compute_rdf(
    elements: &[u8],
    positions: &[[f64; 3]],
    charges: &[f64],
    beta: f64,
    r_max: f64,
    n_points: usize,
) -> RdfDescriptors {
    let n = elements.len().min(positions.len());
    let beta = if beta <= 0.0 { 100.0 } else { beta };
    let r_max = if r_max <= 0.0 { 12.0 } else { r_max };

    if n < 2 || n_points == 0 {
        return RdfDescriptors {
            radii: vec![],
            rdf_unweighted: vec![],
            rdf_mass: vec![],
            rdf_electronegativity: vec![],
            rdf_charge: vec![],
            beta,
            n_pairs: 0,
        };
    }

    let dr = r_max / n_points as f64;
    let radii: Vec<f64> = (0..n_points).map(|k| (k as f64 + 0.5) * dr).collect();

    let mut rdf_unweighted = vec![0.0f64; n_points];
    let mut rdf_mass = vec![0.0f64; n_points];
    let mut rdf_en = vec![0.0f64; n_points];
    let mut rdf_charge = vec![0.0f64; n_points];

    let has_charges = charges.len() >= n;
    let mut n_pairs = 0usize;

    for i in 0..n {
        let mi = atomic_mass_rdf(elements[i]);
        let eni = electronegativity_rdf(elements[i]);
        let qi = if has_charges { charges[i] } else { 0.0 };

        for j in (i + 1)..n {
            let mj = atomic_mass_rdf(elements[j]);
            let enj = electronegativity_rdf(elements[j]);
            let qj = if has_charges { charges[j] } else { 0.0 };

            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let rij = (dx * dx + dy * dy + dz * dz).sqrt();

            n_pairs += 1;

            for (k, &rk) in radii.iter().enumerate() {
                let g = (-beta * (rk - rij) * (rk - rij)).exp();
                rdf_unweighted[k] += g;
                rdf_mass[k] += mi * mj * g;
                rdf_en[k] += eni * enj * g;
                if has_charges {
                    rdf_charge[k] += qi * qj * g;
                }
            }
        }
    }

    // Normalize by 4πr² shell volume to get proper RDF
    for k in 0..n_points {
        let r = radii[k];
        let shell = 4.0 * std::f64::consts::PI * r * r;
        if shell > 1e-12 {
            rdf_unweighted[k] /= shell;
            rdf_mass[k] /= shell;
            rdf_en[k] /= shell;
            rdf_charge[k] /= shell;
        }
    }

    RdfDescriptors {
        radii,
        rdf_unweighted,
        rdf_mass,
        rdf_electronegativity: rdf_en,
        rdf_charge,
        beta,
        n_pairs,
    }
}

fn atomic_mass_rdf(z: u8) -> f64 {
    match z {
        1 => 1.008,
        6 => 12.011,
        7 => 14.007,
        8 => 15.999,
        9 => 18.998,
        15 => 30.974,
        16 => 32.065,
        17 => 35.453,
        35 => 79.904,
        53 => 126.904,
        _ => z as f64 * 1.5,
    }
}

fn electronegativity_rdf(z: u8) -> f64 {
    match z {
        1 => 2.20,
        6 => 2.55,
        7 => 3.04,
        8 => 3.44,
        9 => 3.98,
        14 => 1.90,
        15 => 2.19,
        16 => 2.58,
        17 => 3.16,
        35 => 2.96,
        53 => 2.66,
        _ => 2.00,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rdf_water() {
        let elements = vec![8, 1, 1];
        let positions = vec![
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469],
        ];

        let rdf = compute_rdf(&elements, &positions, &[], 100.0, 12.0, 128);
        assert_eq!(rdf.n_pairs, 3);
        assert_eq!(rdf.radii.len(), 128);
        // Should have peaks near O-H (~0.96 Å) and H-H (~1.51 Å)
        let peak_idx = rdf
            .rdf_unweighted
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        assert!(rdf.radii[peak_idx] < 3.0); // Peak should be near bond distances
    }

    #[test]
    fn test_rdf_with_charges() {
        let elements = vec![8, 1, 1];
        let positions = vec![
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469],
        ];
        let charges = vec![-0.8, 0.4, 0.4];

        let rdf = compute_rdf(&elements, &positions, &charges, 100.0, 12.0, 128);
        // Charge-weighted RDF should have negative values (qi*qj < 0 for O-H pairs)
        let has_negative = rdf.rdf_charge.iter().any(|&v| v < 0.0);
        assert!(has_negative);
    }
}
