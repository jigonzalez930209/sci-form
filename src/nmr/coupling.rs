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

/// Karplus parameters for a specific coupling pathway.
#[derive(Debug, Clone, Copy)]
pub struct KarplusParams {
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

impl KarplusParams {
    /// Evaluate Karplus equation: ³J = A cos²φ + B cosφ + C
    pub fn evaluate(&self, phi_rad: f64) -> f64 {
        let cos_phi = phi_rad.cos();
        self.a * cos_phi * cos_phi + self.b * cos_phi + self.c
    }
}

/// Get Karplus parameters based on the coupled pathway atoms.
/// Returns pathway-specific parameters for H-X-Y-H coupling.
fn get_karplus_params(x_elem: u8, y_elem: u8) -> KarplusParams {
    match (x_elem, y_elem) {
        // H-C-C-H: Standard Altona & Sundaralingam (1972)
        (6, 6) => KarplusParams {
            a: 7.76,
            b: -1.10,
            c: 1.40,
        },
        // H-C-N-H: Bystrov parameters for peptides
        (6, 7) | (7, 6) => KarplusParams {
            a: 6.40,
            b: -1.40,
            c: 1.90,
        },
        // H-C-O-H: Parameters for glycosidic/hydroxyl couplings
        (6, 8) | (8, 6) => KarplusParams {
            a: 5.80,
            b: -1.20,
            c: 1.50,
        },
        // H-C-S-H: Thiol coupling
        (6, 16) | (16, 6) => KarplusParams {
            a: 6.00,
            b: -1.00,
            c: 1.30,
        },
        // Default: use standard H-C-C-H parameters
        _ => KarplusParams {
            a: KARPLUS_A,
            b: KARPLUS_B,
            c: KARPLUS_C,
        },
    }
}

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

/// Predict all J-coupling constants for a molecule.
///
/// `mol`: parsed molecular graph
/// `positions`: 3D coordinates (one per atom, or empty for topology-only estimate)
///
/// Currently supports:
/// - ²J (geminal): topological estimate (typically 2–12 Hz)
/// - ³J (vicinal): Karplus equation if 3D coords available, else topological estimate (6–8 Hz)
///   Uses pathway-specific parameters (H-C-C-H, H-C-N-H, H-C-O-H, etc.)
pub fn predict_j_couplings(mol: &crate::graph::Molecule, positions: &[[f64; 3]]) -> Vec<JCoupling> {
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
                let j_hz: f64 = -12.0; // typical geminal coupling (negative)
                couplings.push(JCoupling {
                    h1_index: h1,
                    h2_index: h2,
                    j_hz: j_hz.abs(), // report magnitude
                    n_bonds: 2,
                    coupling_type: "geminal_2J".to_string(),
                });
            } else if mol.graph.find_edge(p1, p2).is_some() {
                // ³J vicinal coupling (H-X-Y-H pathway)
                let x_elem = mol.graph[p1].element;
                let y_elem = mol.graph[p2].element;
                let params = get_karplus_params(x_elem, y_elem);

                let j_hz = if has_3d {
                    let phi = dihedral_angle(
                        &positions[h1],
                        &positions[p1.index()],
                        &positions[p2.index()],
                        &positions[h2],
                    );
                    params.evaluate(phi)
                } else {
                    // Fallback: average vicinal value
                    7.0
                };

                couplings.push(JCoupling {
                    h1_index: h1,
                    h2_index: h2,
                    j_hz,
                    n_bonds: 3,
                    coupling_type: format!(
                        "vicinal_3J_H-{}-{}-H",
                        element_symbol(x_elem),
                        element_symbol(y_elem)
                    ),
                });
            } else {
                // ⁴J long-range coupling: check for H-X-Y-Z-H (4-bond path)
                // Includes W-coupling and allylic coupling patterns
                let p1_neighbors: Vec<petgraph::graph::NodeIndex> = mol
                    .graph
                    .neighbors(p1)
                    .filter(|&nb| nb != h1_idx && mol.graph[nb].element != 1)
                    .collect();
                for &mid in &p1_neighbors {
                    if mol.graph.find_edge(mid, p2).is_some() {
                        // Found 4-bond path: H1-p1-mid-p2-H2
                        let mid_elem = mol.graph[mid].element;
                        let is_sp2_mid =
                            mol.graph[mid].hybridization == crate::graph::Hybridization::SP2;
                        let is_sp2_p1 =
                            mol.graph[p1].hybridization == crate::graph::Hybridization::SP2;
                        let is_sp2_p2 =
                            mol.graph[p2].hybridization == crate::graph::Hybridization::SP2;

                        // Allylic ⁴J: sp3-sp2=sp2-sp3 pattern, typically 1-3 Hz
                        // W-coupling: rigid W-shaped path, typically 1-3 Hz
                        let j_hz = if is_sp2_mid && (is_sp2_p1 || is_sp2_p2) {
                            2.0 // allylic coupling
                        } else {
                            1.0 // generic long-range
                        };

                        couplings.push(JCoupling {
                            h1_index: h1,
                            h2_index: h2,
                            j_hz,
                            n_bonds: 4,
                            coupling_type: format!(
                                "long_range_4J_H-{}-{}-{}-H",
                                element_symbol(mol.graph[p1].element),
                                element_symbol(mid_elem),
                                element_symbol(mol.graph[p2].element)
                            ),
                        });
                        break; // only report one 4-bond path per pair
                    }
                }

                // ⁵J long-range coupling: H-X-Y-Z-W-H (5-bond path)
                // Observed in aromatic systems and extended conjugation, typically 0.2–1.0 Hz
                if couplings.last().map_or(true, |c| c.n_bonds != 4 || c.h1_index != h1 || c.h2_index != h2) {
                    'five_bond: for &mid1 in &p1_neighbors {
                        let mid1_neighbors: Vec<petgraph::graph::NodeIndex> = mol
                            .graph
                            .neighbors(mid1)
                            .filter(|&nb| nb != p1 && nb != h1_idx && mol.graph[nb].element != 1)
                            .collect();
                        for &mid2 in &mid1_neighbors {
                            if mol.graph.find_edge(mid2, p2).is_some() {
                                // Found 5-bond path: H1-p1-mid1-mid2-p2-H2
                                let is_aromatic_path = [p1, mid1, mid2, p2].iter().all(|&n| {
                                    mol.graph[n].hybridization == crate::graph::Hybridization::SP2
                                });

                                let j_hz = if is_aromatic_path {
                                    0.7 // aromatic ⁵J
                                } else {
                                    0.3 // generic five-bond
                                };

                                couplings.push(JCoupling {
                                    h1_index: h1,
                                    h2_index: h2,
                                    j_hz,
                                    n_bonds: 5,
                                    coupling_type: format!(
                                        "long_range_5J_H-{}-{}-{}-{}-H",
                                        element_symbol(mol.graph[p1].element),
                                        element_symbol(mol.graph[mid1].element),
                                        element_symbol(mol.graph[mid2].element),
                                        element_symbol(mol.graph[p2].element)
                                    ),
                                });
                                break 'five_bond;
                            }
                        }
                    }
                }
            }
        }
    }

    couplings
}

fn element_symbol(z: u8) -> &'static str {
    match z {
        6 => "C",
        7 => "N",
        8 => "O",
        16 => "S",
        _ => "X",
    }
}

/// Compute Boltzmann-weighted ensemble-averaged J-couplings from multiple conformers.
///
/// Given conformer coordinates and their energies, compute the Boltzmann-weighted
/// average of ³J coupling constants:
///
/// $\langle {}^3J \rangle = \frac{\sum_i J_i \cdot e^{-E_i / k_B T}}{\sum_i e^{-E_i / k_B T}}$
///
/// # Arguments
/// - `mol`: parsed molecular graph
/// - `conformer_positions`: list of 3D coordinate sets (one per conformer)
/// - `energies_kcal`: relative energies in kcal/mol for each conformer
/// - `temperature_k`: temperature in Kelvin (default 298.15)
pub fn ensemble_averaged_j_couplings(
    mol: &crate::graph::Molecule,
    conformer_positions: &[Vec<[f64; 3]>],
    energies_kcal: &[f64],
    temperature_k: f64,
) -> Vec<JCoupling> {
    if conformer_positions.is_empty() {
        return Vec::new();
    }
    if conformer_positions.len() != energies_kcal.len() {
        return predict_j_couplings(mol, &conformer_positions[0]);
    }

    // k_B in kcal/(mol·K)
    const KB_KCAL: f64 = 0.001987204;
    let beta = 1.0 / (KB_KCAL * temperature_k);

    // Find minimum energy for numerical stability
    let e_min = energies_kcal.iter().cloned().fold(f64::INFINITY, f64::min);

    // Compute Boltzmann weights
    let weights: Vec<f64> = energies_kcal
        .iter()
        .map(|&e| (-(e - e_min) * beta).exp())
        .collect();
    let weight_sum: f64 = weights.iter().sum();

    if weight_sum < 1e-30 {
        return predict_j_couplings(mol, &conformer_positions[0]);
    }

    // Compute J-couplings for each conformer
    let all_couplings: Vec<Vec<JCoupling>> = conformer_positions
        .iter()
        .map(|pos| predict_j_couplings(mol, pos))
        .collect();

    // Average over conformers
    if all_couplings.is_empty() {
        return Vec::new();
    }

    let n_couplings = all_couplings[0].len();
    let mut averaged = all_couplings[0].clone();

    for k in 0..n_couplings {
        let mut weighted_j = 0.0;
        for (conf_idx, couplings) in all_couplings.iter().enumerate() {
            if k < couplings.len() {
                weighted_j += couplings[k].j_hz * weights[conf_idx];
            }
        }
        averaged[k].j_hz = weighted_j / weight_sum;
    }

    averaged
}

/// Parallel ensemble-averaged J-couplings via rayon.
///
/// Each conformer's J-coupling prediction is independent, making
/// this embarrassingly parallel over conformers.
#[cfg(feature = "parallel")]
pub fn ensemble_averaged_j_couplings_parallel(
    mol: &crate::graph::Molecule,
    conformer_positions: &[Vec<[f64; 3]>],
    energies_kcal: &[f64],
    temperature_k: f64,
) -> Vec<JCoupling> {
    use rayon::prelude::*;

    if conformer_positions.is_empty() {
        return Vec::new();
    }
    if conformer_positions.len() != energies_kcal.len() {
        return predict_j_couplings(mol, &conformer_positions[0]);
    }

    const KB_KCAL: f64 = 0.001987204;
    let beta = 1.0 / (KB_KCAL * temperature_k);
    let e_min = energies_kcal.iter().cloned().fold(f64::INFINITY, f64::min);

    let weights: Vec<f64> = energies_kcal
        .iter()
        .map(|&e| (-(e - e_min) * beta).exp())
        .collect();
    let weight_sum: f64 = weights.iter().sum();

    if weight_sum < 1e-30 {
        return predict_j_couplings(mol, &conformer_positions[0]);
    }

    // Parallel J-coupling evaluation per conformer
    let all_couplings: Vec<Vec<JCoupling>> = conformer_positions
        .par_iter()
        .map(|pos| predict_j_couplings(mol, pos))
        .collect();

    if all_couplings.is_empty() {
        return Vec::new();
    }

    let n_couplings = all_couplings[0].len();
    let mut averaged = all_couplings[0].clone();

    for k in 0..n_couplings {
        let mut weighted_j = 0.0;
        for (conf_idx, couplings) in all_couplings.iter().enumerate() {
            if k < couplings.len() {
                weighted_j += couplings[k].j_hz * weights[conf_idx];
            }
        }
        averaged[k].j_hz = weighted_j / weight_sum;
    }

    averaged
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Karplus equation: ³J = A cos²φ + B cosφ + C
    fn karplus_3j(phi_rad: f64) -> f64 {
        let cos_phi = phi_rad.cos();
        KARPLUS_A * cos_phi * cos_phi + KARPLUS_B * cos_phi + KARPLUS_C
    }

    #[test]
    fn test_karplus_equation_values() {
        // At φ = 0° (eclipsed): ³J should be large (~9.06 Hz)
        let j_0 = karplus_3j(0.0);
        assert!(
            j_0 > 8.0 && j_0 < 10.0,
            "³J(0°) = {} Hz, expected ~9 Hz",
            j_0
        );

        // At φ = 90°: ³J should be small (~1.4 Hz)
        let j_90 = karplus_3j(std::f64::consts::FRAC_PI_2);
        assert!(
            j_90 > 0.0 && j_90 < 3.0,
            "³J(90°) = {} Hz, expected ~1.4 Hz",
            j_90
        );

        // At φ = 180° (antiperiplanar): ³J should be large (~10.26 Hz with A-S params)
        let j_180 = karplus_3j(std::f64::consts::PI);
        assert!(
            j_180 > 6.0 && j_180 < 12.0,
            "³J(180°) = {} Hz, expected ~10 Hz",
            j_180
        );
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
        let vicinal: Vec<&JCoupling> = couplings.iter().filter(|c| c.n_bonds == 3).collect();
        assert!(
            !vicinal.is_empty(),
            "Ethane should have ³J vicinal couplings"
        );

        // Vicinal coupling type should include pathway info
        for c in &vicinal {
            assert!(
                c.coupling_type.contains("vicinal_3J"),
                "Coupling type should be vicinal_3J, got {}",
                c.coupling_type
            );
        }
    }

    #[test]
    fn test_karplus_pathway_specific() {
        // H-C-C-H vs H-C-N-H should have different parameters
        let cc_params = get_karplus_params(6, 6);
        let cn_params = get_karplus_params(6, 7);

        // At φ = 0, different pathways should give different J values
        let j_cc = cc_params.evaluate(0.0);
        let j_cn = cn_params.evaluate(0.0);
        assert!(
            (j_cc - j_cn).abs() > 0.1,
            "H-C-C-H and H-C-N-H should have different J at φ=0: {} vs {}",
            j_cc,
            j_cn
        );
    }

    #[test]
    fn test_ensemble_averaging() {
        let mol = crate::graph::Molecule::from_smiles("CC").unwrap();
        // Single conformer: result should be same as regular prediction
        let n = mol.graph.node_count();
        let positions = vec![[0.0, 0.0, 0.0]; n];
        let result = ensemble_averaged_j_couplings(&mol, &[positions], &[0.0], 298.15);
        assert!(!result.is_empty());
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
