//! J-coupling constant estimation via Karplus equation and topological rules.
//!
//! Implements ³J(H,H) prediction using the Karplus equation:
//! ³J = A cos²φ + B cosφ + C
//!
//! where φ is the H-C-C-H dihedral angle.

use serde::{Deserialize, Serialize};

/// A predicted spin-spin coupling constant.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JCoupling {
    /// First hydrogen atom index.
    pub h1_index: usize,
    /// Second hydrogen atom index.
    pub h2_index: usize,
    /// Coupling constant in Hz.
    pub j_hz: f64,
    /// Number of bonds separating the two H atoms.
    pub n_bonds: usize,
    /// Coupling type description.
    pub coupling_type: String,
}

// Karplus equation parameters for ³J(H,H)
// Standard Karplus: A=7.76, B=-1.10, C=1.40 (Altona & Sundaralingam, JACS 1972)
const KARPLUS_A: f64 = 7.76;
const KARPLUS_B: f64 = -1.10;
const KARPLUS_C: f64 = 1.40;

/// Compute dihedral angle (in radians) from four 3D positions.
fn dihedral_angle(p1: &[f64; 3], p2: &[f64; 3], p3: &[f64; 3], p4: &[f64; 3]) -> f64 {
    let b1 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
    let b2 = [p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]];
    let b3 = [p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2]];

    // n1 = b1 × b2
    let n1 = cross(&b1, &b2);
    // n2 = b2 × b3
    let n2 = cross(&b2, &b3);

    let m1 = cross(&n1, &normalize(&b2));

    let x = dot(&n1, &n2);
    let y = dot(&m1, &n2);

    (-y).atan2(x)
}

fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn normalize(v: &[f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-10 {
        return [0.0, 0.0, 0.0];
    }
    [v[0] / len, v[1] / len, v[2] / len]
}

/// Karplus equation: ³J = A cos²φ + B cosφ + C
fn karplus_3j(phi_rad: f64) -> f64 {
    let cos_phi = phi_rad.cos();
    KARPLUS_A * cos_phi * cos_phi + KARPLUS_B * cos_phi + KARPLUS_C
}

/// Predict all J-coupling constants for a molecule.
///
/// `mol`: parsed molecular graph
/// `positions`: 3D coordinates (one per atom, or empty for topology-only estimate)
///
/// Currently supports:
/// - ²J (geminal): topological estimate (typically 2–12 Hz)
/// - ³J (vicinal): Karplus equation if 3D coords available, else topological estimate (6–8 Hz)
pub fn predict_j_couplings(
    mol: &crate::graph::Molecule,
    positions: &[[f64; 3]],
) -> Vec<JCoupling> {
    let n = mol.graph.node_count();
    let has_3d = positions.len() == n;
    let mut couplings = Vec::new();

    // Find all hydrogen atoms
    let h_atoms: Vec<usize> = (0..n)
        .filter(|&i| mol.graph[petgraph::graph::NodeIndex::new(i)].element == 1)
        .collect();

    // For each pair of H atoms, determine coupling pathway
    for i in 0..h_atoms.len() {
        for j in (i + 1)..h_atoms.len() {
            let h1 = h_atoms[i];
            let h2 = h_atoms[j];
            let h1_idx = petgraph::graph::NodeIndex::new(h1);
            let h2_idx = petgraph::graph::NodeIndex::new(h2);

            // Find parent heavy atom for each H
            let parent1: Vec<petgraph::graph::NodeIndex> = mol
                .graph
                .neighbors(h1_idx)
                .filter(|n| mol.graph[*n].element != 1)
                .collect();
            let parent2: Vec<petgraph::graph::NodeIndex> = mol
                .graph
                .neighbors(h2_idx)
                .filter(|n| mol.graph[*n].element != 1)
                .collect();

            if parent1.is_empty() || parent2.is_empty() {
                continue;
            }

            let p1 = parent1[0];
            let p2 = parent2[0];

            if p1 == p2 {
                // ²J geminal coupling (both H on same carbon)
                let j_hz = -12.0; // typical geminal coupling (negative)
                couplings.push(JCoupling {
                    h1_index: h1,
                    h2_index: h2,
                    j_hz: (j_hz as f64).abs(), // report magnitude
                    n_bonds: 2,
                    coupling_type: "geminal_2J".to_string(),
                });
            } else if mol.graph.find_edge(p1, p2).is_some() {
                // ³J vicinal coupling (H-C-C-H pathway)
                let j_hz = if has_3d {
                    let phi = dihedral_angle(
                        &positions[h1],
                        &positions[p1.index()],
                        &positions[p2.index()],
                        &positions[h2],
                    );
                    karplus_3j(phi)
                } else {
                    // Fallback: average vicinal value
                    7.0
                };

                couplings.push(JCoupling {
                    h1_index: h1,
                    h2_index: h2,
                    j_hz,
                    n_bonds: 3,
                    coupling_type: "vicinal_3J".to_string(),
                });
            }
            // We skip longer-range couplings (⁴J, ⁵J) as they are typically < 3 Hz
        }
    }

    couplings
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_karplus_equation_values() {
        // At φ = 0° (eclipsed): ³J should be large (~9.06 Hz)
        let j_0 = karplus_3j(0.0);
        assert!(j_0 > 8.0 && j_0 < 10.0, "³J(0°) = {} Hz, expected ~9 Hz", j_0);

        // At φ = 90°: ³J should be small (~1.4 Hz)
        let j_90 = karplus_3j(std::f64::consts::FRAC_PI_2);
        assert!(j_90 > 0.0 && j_90 < 3.0, "³J(90°) = {} Hz, expected ~1.4 Hz", j_90);

        // At φ = 180° (antiperiplanar): ³J should be large (~10.26 Hz with A-S params)
        let j_180 = karplus_3j(std::f64::consts::PI);
        assert!(j_180 > 6.0 && j_180 < 12.0, "³J(180°) = {} Hz, expected ~10 Hz", j_180);
    }

    #[test]
    fn test_dihedral_angle_basic() {
        // Simple planar dihedral
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [0.0, 0.0, 0.0];
        let p3 = [0.0, 1.0, 0.0];
        let p4 = [-1.0, 1.0, 0.0];
        let angle = dihedral_angle(&p1, &p2, &p3, &p4);
        // Should be 0° (all in same plane)
        assert!(angle.abs() < 0.1 || (angle.abs() - std::f64::consts::PI).abs() < 0.1);
    }

    #[test]
    fn test_ethane_j_couplings() {
        let mol = crate::graph::Molecule::from_smiles("CC").unwrap();
        let couplings = predict_j_couplings(&mol, &[]);

        // Ethane should have geminal ²J and vicinal ³J couplings
        assert!(
            !couplings.is_empty(),
            "Ethane should have J-coupling predictions"
        );

        // Check that we find vicinal couplings
        let vicinal: Vec<&JCoupling> = couplings
            .iter()
            .filter(|c| c.n_bonds == 3)
            .collect();
        assert!(
            !vicinal.is_empty(),
            "Ethane should have ³J vicinal couplings"
        );
    }

    #[test]
    fn test_methane_j_couplings() {
        let mol = crate::graph::Molecule::from_smiles("C").unwrap();
        let couplings = predict_j_couplings(&mol, &[]);

        // Methane: all H on same C → geminal ²J only
        for c in &couplings {
            assert_eq!(c.n_bonds, 2, "Methane H-H should be 2-bond (geminal)");
        }
    }
}
