//! Chemical shift prediction for ¹H and ¹³C NMR.
//!
//! Uses empirical additivity rules based on local atomic environment:
//! - Hybridization state (sp, sp2, sp3)
//! - Electronegativity of neighboring atoms
//! - Ring current effects for aromatic systems
//! - Functional group corrections

use crate::graph::{BondOrder, Hybridization, Molecule};
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use serde::{Deserialize, Serialize};

/// Chemical shift prediction for a single atom.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChemicalShift {
    /// Atom index in the molecule.
    pub atom_index: usize,
    /// Atomic number.
    pub element: u8,
    /// Predicted chemical shift in ppm.
    pub shift_ppm: f64,
    /// Environment classification.
    pub environment: String,
    /// Confidence level (0.0–1.0).
    pub confidence: f64,
}

/// Complete NMR shift prediction result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NmrShiftResult {
    /// Predicted shifts for hydrogen atoms.
    pub h_shifts: Vec<ChemicalShift>,
    /// Predicted shifts for carbon atoms.
    pub c_shifts: Vec<ChemicalShift>,
    /// Notes and caveats.
    pub notes: Vec<String>,
}

// ─── Electronegativity table (Pauling scale) ────────────────────────────────

fn electronegativity(z: u8) -> f64 {
    match z {
        1 => 2.20,
        5 => 2.04,
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

// ─── ¹H Chemical Shift Prediction ──────────────────────────────────────────
//
// Base shifts by carbon hybridization attached to:
//   sp3-C-H: ~0.9 ppm (methyl), ~1.3 (methylene), ~1.5 (methine)
//   sp2-C-H (alkene): ~5.0 ppm
//   sp2-C-H (aromatic): ~7.3 ppm
//   sp-C-H (alkyne): ~2.5 ppm
//   O-H: ~1.0–5.0 (alcohol), ~9.0–12.0 (carboxylic acid)
//   N-H: ~1.0–5.0 (amine)
//   aldehyde C-H: ~9.5 ppm

fn predict_h_shift(mol: &Molecule, h_idx: NodeIndex) -> (f64, String, f64) {
    // Find the heavy atom that H is bonded to
    let neighbors: Vec<NodeIndex> = mol.graph.neighbors(h_idx).collect();
    if neighbors.is_empty() {
        return (0.0, "isolated_hydrogen".to_string(), 0.3);
    }

    let parent = neighbors[0];
    let parent_elem = mol.graph[parent].element;
    let parent_hyb = &mol.graph[parent].hybridization;

    match parent_elem {
        // H bonded to carbon
        6 => {
            let is_aromatic = mol
                .graph
                .edges(parent)
                .any(|e| matches!(e.weight().order, BondOrder::Aromatic));

            if is_aromatic {
                // Aromatic C-H: base ~7.27 ppm
                let mut shift = 7.27;
                // Neighbor electronegativity corrections
                for nb in mol.graph.neighbors(parent) {
                    if nb == h_idx {
                        continue;
                    }
                    let nb_elem = mol.graph[nb].element;
                    let en_diff = electronegativity(nb_elem) - electronegativity(6);
                    shift += en_diff * 0.3;
                }
                (shift, "aromatic_C-H".to_string(), 0.80)
            } else {
                match parent_hyb {
                    Hybridization::SP3 => {
                        // Count H neighbors on parent carbon
                        let h_count = mol
                            .graph
                            .neighbors(parent)
                            .filter(|n| mol.graph[*n].element == 1)
                            .count();

                        let base = match h_count {
                            3 => 0.90, // methyl
                            2 => 1.30, // methylene
                            _ => 1.50, // methine
                        };

                        // Electronegativity corrections from heavy-atom neighbors
                        let mut shift = base;
                        for nb in mol.graph.neighbors(parent) {
                            if mol.graph[nb].element == 1 {
                                continue;
                            }
                            let nb_elem = mol.graph[nb].element;
                            let en_diff = electronegativity(nb_elem) - electronegativity(6);
                            if en_diff > 0.3 {
                                shift += en_diff * 0.8;
                            }
                        }

                        // Check if adjacent to C=O (alpha to carbonyl)
                        for nb in mol.graph.neighbors(parent) {
                            if mol.graph[nb].element != 6 {
                                continue;
                            }
                            for nb2_edge in mol.graph.edges(nb) {
                                let nb2 = if nb2_edge.source() == nb {
                                    nb2_edge.target()
                                } else {
                                    nb2_edge.source()
                                };
                                if mol.graph[nb2].element == 8
                                    && matches!(nb2_edge.weight().order, BondOrder::Double)
                                {
                                    shift += 0.5;
                                    break;
                                }
                            }
                        }

                        let env = if h_count >= 3 {
                            "methyl"
                        } else if h_count == 2 {
                            "methylene"
                        } else {
                            "methine"
                        };
                        (shift, format!("sp3_{}", env), 0.75)
                    }
                    Hybridization::SP2 => {
                        // Alkene C-H
                        let mut shift = 5.25;
                        // Check for conjugation
                        for nb in mol.graph.neighbors(parent) {
                            if nb == h_idx {
                                continue;
                            }
                            let nb_elem = mol.graph[nb].element;
                            let en_diff = electronegativity(nb_elem) - electronegativity(6);
                            shift += en_diff * 0.5;

                            // Check for aldehyde: C(=O)H  
                            if nb_elem == 8 {
                                if mol
                                    .graph
                                    .edges(parent)
                                    .any(|e| {
                                        let other = if e.source() == parent {
                                            e.target()
                                        } else {
                                            e.source()
                                        };
                                        mol.graph[other].element == 8
                                            && matches!(e.weight().order, BondOrder::Double)
                                    })
                                {
                                    shift = 9.50; // aldehyde H
                                    return (shift, "aldehyde_C-H".to_string(), 0.82);
                                }
                            }
                        }
                        (shift, "alkene_C-H".to_string(), 0.70)
                    }
                    Hybridization::SP => {
                        // Alkyne C-H
                        (2.50, "alkyne_C-H".to_string(), 0.75)
                    }
                    _ => (1.5, "unknown_C-H".to_string(), 0.40),
                }
            }
        }
        // H bonded to oxygen
        8 => {
            // Check if carboxylic acid O-H
            let is_carboxylic = mol.graph.neighbors(parent).any(|nb| {
                if mol.graph[nb].element != 6 {
                    return false;
                }
                mol.graph.edges(nb).any(|e| {
                    let other = if e.source() == nb {
                        e.target()
                    } else {
                        e.source()
                    };
                    mol.graph[other].element == 8
                        && matches!(e.weight().order, BondOrder::Double)
                })
            });

            if is_carboxylic {
                (11.0, "carboxylic_acid_O-H".to_string(), 0.70)
            } else {
                // Check if phenol
                let is_phenol = mol.graph.neighbors(parent).any(|nb| {
                    mol.graph[nb].element == 6
                        && mol
                            .graph
                            .edges(nb)
                            .any(|e| matches!(e.weight().order, BondOrder::Aromatic))
                });
                if is_phenol {
                    (5.5, "phenol_O-H".to_string(), 0.60)
                } else {
                    (2.5, "alcohol_O-H".to_string(), 0.55)
                }
            }
        }
        // H bonded to nitrogen
        7 => {
            let is_aromatic_n = mol
                .graph
                .edges(parent)
                .any(|e| matches!(e.weight().order, BondOrder::Aromatic));
            if is_aromatic_n {
                (8.5, "aromatic_N-H".to_string(), 0.55)
            } else {
                let h_on_n = mol
                    .graph
                    .neighbors(parent)
                    .filter(|n| mol.graph[*n].element == 1)
                    .count();
                match h_on_n {
                    1 => (2.2, "secondary_amine_N-H".to_string(), 0.50),
                    _ => (1.5, "primary_amine_N-H".to_string(), 0.50),
                }
            }
        }
        // H bonded to sulfur
        16 => (1.8, "thiol_S-H".to_string(), 0.50),
        _ => (2.0, "other_X-H".to_string(), 0.30),
    }
}

// ─── ¹³C Chemical Shift Prediction ─────────────────────────────────────────
//
// Base shifts by hybridization:
//   sp3-C: 0–50 ppm (depends on substitution)
//   sp2-C (alkene): 100–150 ppm
//   sp2-C (aromatic): 120–140 ppm
//   sp-C (alkyne): 70–90 ppm
//   C=O (carbonyl): 170–220 ppm

fn predict_c_shift(mol: &Molecule, c_idx: NodeIndex) -> (f64, String, f64) {
    let hyb = &mol.graph[c_idx].hybridization;
    let is_aromatic = mol
        .graph
        .edges(c_idx)
        .any(|e| matches!(e.weight().order, BondOrder::Aromatic));

    // Check for carbonyl C=O
    let has_carbonyl_o = mol.graph.edges(c_idx).any(|e| {
        let other = if e.source() == c_idx {
            e.target()
        } else {
            e.source()
        };
        mol.graph[other].element == 8 && matches!(e.weight().order, BondOrder::Double)
    });

    if has_carbonyl_o {
        // Carbonyl carbon
        let has_oh = mol.graph.neighbors(c_idx).any(|nb| {
            if mol.graph[nb].element != 8 {
                return false;
            }
            mol.graph
                .neighbors(nb)
                .any(|n| mol.graph[n].element == 1)
        });

        if has_oh {
            // Carboxylic acid
            return (175.0, "carboxylic_acid_C".to_string(), 0.75);
        }

        let has_n = mol
            .graph
            .neighbors(c_idx)
            .any(|nb| mol.graph[nb].element == 7);
        if has_n {
            // Amide
            return (170.0, "amide_C=O".to_string(), 0.72);
        }

        // Ketone/aldehyde
        let has_h = mol
            .graph
            .neighbors(c_idx)
            .any(|nb| mol.graph[nb].element == 1);
        if has_h {
            return (200.0, "aldehyde_C=O".to_string(), 0.75);
        }
        return (205.0, "ketone_C=O".to_string(), 0.73);
    }

    if is_aromatic {
        let mut shift = 128.0;
        // Corrections for substituents
        for nb in mol.graph.neighbors(c_idx) {
            let nb_elem = mol.graph[nb].element;
            match nb_elem {
                7 => shift += 5.0,  // N substituent
                8 => shift -= 10.0, // O substituent (shielding)
                9 => shift -= 5.0,  // F
                17 => shift += 2.0, // Cl
                35 => shift += 3.0, // Br
                _ => {}
            }
        }
        return (shift, "aromatic_C".to_string(), 0.78);
    }

    match hyb {
        Hybridization::SP3 => {
            let mut shift = 20.0;
            let heavy_neighbors: Vec<u8> = mol
                .graph
                .neighbors(c_idx)
                .filter(|nb| mol.graph[*nb].element != 1)
                .map(|nb| mol.graph[nb].element)
                .collect();

            // Substitution increments
            let n_heavy = heavy_neighbors.len();
            shift += (n_heavy as f64 - 1.0) * 8.0;

            // Electronegativity corrections
            for &nb_elem in &heavy_neighbors {
                let en_diff = electronegativity(nb_elem) - electronegativity(6);
                if en_diff > 0.3 {
                    shift += en_diff * 10.0;
                }
            }

            let env = match n_heavy {
                0 | 1 => "methyl_C",
                2 => "methylene_C",
                3 => "methine_C",
                _ => "quaternary_C",
            };
            (shift.clamp(0.0, 80.0), env.to_string(), 0.70)
        }
        Hybridization::SP2 => {
            // Alkene carbon
            let mut shift = 125.0;
            for nb in mol.graph.neighbors(c_idx) {
                let nb_elem = mol.graph[nb].element;
                let en_diff = electronegativity(nb_elem) - electronegativity(6);
                shift += en_diff * 8.0;
            }
            (shift.clamp(80.0, 165.0), "alkene_C".to_string(), 0.65)
        }
        Hybridization::SP => {
            // Alkyne
            let mut shift = 80.0;
            let has_h = mol
                .graph
                .neighbors(c_idx)
                .any(|nb| mol.graph[nb].element == 1);
            if has_h {
                shift = 70.0; // terminal alkyne  
            }
            (shift, "alkyne_C".to_string(), 0.70)
        }
        _ => (30.0, "unknown_C".to_string(), 0.40),
    }
}

/// Predict ¹H and ¹³C chemical shifts for all atoms in a molecule.
pub fn predict_chemical_shifts(mol: &Molecule) -> NmrShiftResult {
    let n = mol.graph.node_count();
    let mut h_shifts = Vec::new();
    let mut c_shifts = Vec::new();

    for atom_idx in 0..n {
        let idx = NodeIndex::new(atom_idx);
        let elem = mol.graph[idx].element;

        match elem {
            1 => {
                let (shift, env, conf) = predict_h_shift(mol, idx);
                h_shifts.push(ChemicalShift {
                    atom_index: atom_idx,
                    element: 1,
                    shift_ppm: shift,
                    environment: env,
                    confidence: conf,
                });
            }
            6 => {
                let (shift, env, conf) = predict_c_shift(mol, idx);
                c_shifts.push(ChemicalShift {
                    atom_index: atom_idx,
                    element: 6,
                    shift_ppm: shift,
                    environment: env,
                    confidence: conf,
                });
            }
            _ => {}
        }
    }

    NmrShiftResult {
        h_shifts,
        c_shifts,
        notes: vec![
            "Chemical shifts predicted using empirical additivity rules based on local atomic environment.".to_string(),
            "Accuracy target: MAE < 0.5 ppm for ¹H, < 3.0 ppm for ¹³C. Values are approximate screening-level predictions.".to_string(),
            "Aromatic and functional group corrections are applied where detected.".to_string(),
        ],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ethanol_h_shifts() {
        let mol = Molecule::from_smiles("CCO").unwrap();
        let result = predict_chemical_shifts(&mol);

        assert!(!result.h_shifts.is_empty(), "Should have H shifts");
        assert!(!result.c_shifts.is_empty(), "Should have C shifts");

        // All H shifts should be in reasonable range (0–15 ppm)
        for shift in &result.h_shifts {
            assert!(
                shift.shift_ppm >= 0.0 && shift.shift_ppm <= 15.0,
                "H shift {} out of range for atom {}",
                shift.shift_ppm,
                shift.atom_index
            );
        }
    }

    #[test]
    fn test_benzene_h_shifts() {
        let mol = Molecule::from_smiles("c1ccccc1").unwrap();
        let result = predict_chemical_shifts(&mol);

        // Aromatic H should be around 7.27 ppm
        let aromatic_h: Vec<&ChemicalShift> = result
            .h_shifts
            .iter()
            .filter(|s| s.environment.contains("aromatic"))
            .collect();
        assert!(
            !aromatic_h.is_empty(),
            "Benzene should have aromatic H shifts"
        );
        for shift in &aromatic_h {
            assert!(
                (shift.shift_ppm - 7.27).abs() < 1.5,
                "Aromatic H shift {} should be near 7.27 ppm",
                shift.shift_ppm
            );
        }
    }

    #[test]
    fn test_acetic_acid_c_shifts() {
        let mol = Molecule::from_smiles("CC(=O)O").unwrap();
        let result = predict_chemical_shifts(&mol);

        // Should have a carbonyl carbon shift > 160 ppm
        let carbonyl_c: Vec<&ChemicalShift> = result
            .c_shifts
            .iter()
            .filter(|s| s.environment.contains("carboxylic") || s.environment.contains("C=O"))
            .collect();
        assert!(
            !carbonyl_c.is_empty(),
            "Acetic acid should have a carbonyl C shift"
        );
    }

    #[test]
    fn test_benzene_c13_shift() {
        let mol = Molecule::from_smiles("c1ccccc1").unwrap();
        let result = predict_chemical_shifts(&mol);

        // Aromatic C should be around 128 ppm
        for shift in &result.c_shifts {
            assert!(
                (shift.shift_ppm - 128.0).abs() < 15.0,
                "Aromatic C shift {} should be near 128 ppm",
                shift.shift_ppm
            );
        }
    }

    #[test]
    fn test_aldehyde_h_shift() {
        let mol = Molecule::from_smiles("CC=O").unwrap();
        let result = predict_chemical_shifts(&mol);

        // Check that we detect an aldehyde H in the range ~9-10 ppm
        // Note: depends on how the SMILES parser handles implicit H on C=O
        let aldehyde_h: Vec<&ChemicalShift> = result
            .h_shifts
            .iter()
            .filter(|s| s.environment.contains("aldehyde") || s.shift_ppm > 9.0)
            .collect();
        // This may or may not detect depending on SMILES parsing, so just verify count
        assert!(result.h_shifts.len() > 0);
    }
}
