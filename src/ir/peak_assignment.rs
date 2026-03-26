//! Automatic IR peak assignment using group frequency tables.
//!
//! Maps computed vibrational frequencies to chemical functional groups
//! by matching against standard IR group frequency ranges and analyzing
//! the normal mode character (which atoms move).

use serde::{Deserialize, Serialize};

/// A group frequency range for functional group identification.
#[derive(Debug, Clone)]
pub struct GroupFrequency {
    /// Functional group name.
    pub group: &'static str,
    /// Lower frequency bound (cm⁻¹).
    pub freq_min: f64,
    /// Upper frequency bound (cm⁻¹).
    pub freq_max: f64,
    /// Expected intensity: "strong", "medium", "weak", "variable".
    pub intensity: &'static str,
    /// Description of the vibration type.
    pub vibration_type: &'static str,
}

/// Assigned IR peak with functional group identification.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PeakAssignment {
    /// Frequency (cm⁻¹).
    pub frequency: f64,
    /// Computed intensity.
    pub intensity: f64,
    /// Assigned functional group(s).
    pub assignments: Vec<FunctionalGroupMatch>,
    /// Confidence score (0-1).
    pub confidence: f64,
}

/// A functional group match for a peak.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FunctionalGroupMatch {
    /// Group name.
    pub group: String,
    /// Vibration type (stretch, bend, etc.).
    pub vibration_type: String,
    /// How well the frequency matches (0-1).
    pub match_quality: f64,
}

/// IR peak assignment result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssignmentResult {
    /// Assigned peaks.
    pub assignments: Vec<PeakAssignment>,
    /// Identified functional groups in the molecule.
    pub functional_groups: Vec<String>,
}

/// Standard IR group frequency table.
const GROUP_FREQUENCIES: &[GroupFrequency] = &[
    // O-H stretches
    GroupFrequency {
        group: "O-H (free)",
        freq_min: 3590.0,
        freq_max: 3650.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "O-H (H-bonded)",
        freq_min: 3200.0,
        freq_max: 3570.0,
        intensity: "strong, broad",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "O-H (carboxylic)",
        freq_min: 2500.0,
        freq_max: 3300.0,
        intensity: "strong, broad",
        vibration_type: "stretch",
    },
    // N-H stretches
    GroupFrequency {
        group: "N-H (primary amine)",
        freq_min: 3300.0,
        freq_max: 3500.0,
        intensity: "medium",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "N-H (secondary amine)",
        freq_min: 3310.0,
        freq_max: 3350.0,
        intensity: "medium",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "N-H (amide)",
        freq_min: 3180.0,
        freq_max: 3360.0,
        intensity: "medium",
        vibration_type: "stretch",
    },
    // C-H stretches
    GroupFrequency {
        group: "C-H (sp3)",
        freq_min: 2850.0,
        freq_max: 2960.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C-H (sp2)",
        freq_min: 3010.0,
        freq_max: 3100.0,
        intensity: "medium",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C-H (sp, alkyne)",
        freq_min: 3300.0,
        freq_max: 3320.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C-H (aldehyde)",
        freq_min: 2700.0,
        freq_max: 2830.0,
        intensity: "weak",
        vibration_type: "stretch",
    },
    // Triple bonds
    GroupFrequency {
        group: "C≡C",
        freq_min: 2100.0,
        freq_max: 2260.0,
        intensity: "weak",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C≡N",
        freq_min: 2210.0,
        freq_max: 2260.0,
        intensity: "medium",
        vibration_type: "stretch",
    },
    // Double bonds
    GroupFrequency {
        group: "C=O (ketone)",
        freq_min: 1705.0,
        freq_max: 1725.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C=O (aldehyde)",
        freq_min: 1720.0,
        freq_max: 1740.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C=O (carboxylic acid)",
        freq_min: 1700.0,
        freq_max: 1725.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C=O (ester)",
        freq_min: 1735.0,
        freq_max: 1750.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C=O (amide)",
        freq_min: 1630.0,
        freq_max: 1680.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C=C (alkene)",
        freq_min: 1620.0,
        freq_max: 1680.0,
        intensity: "variable",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C=C (aromatic)",
        freq_min: 1450.0,
        freq_max: 1600.0,
        intensity: "medium",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C=N",
        freq_min: 1600.0,
        freq_max: 1660.0,
        intensity: "medium",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "N=O (nitro)",
        freq_min: 1515.0,
        freq_max: 1560.0,
        intensity: "strong",
        vibration_type: "asymm. stretch",
    },
    GroupFrequency {
        group: "N=O (nitro)",
        freq_min: 1340.0,
        freq_max: 1385.0,
        intensity: "strong",
        vibration_type: "symm. stretch",
    },
    GroupFrequency {
        group: "C=O (acid chloride)",
        freq_min: 1780.0,
        freq_max: 1815.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C=O (anhydride)",
        freq_min: 1760.0,
        freq_max: 1810.0,
        intensity: "strong",
        vibration_type: "asymm. stretch",
    },
    GroupFrequency {
        group: "C=O (anhydride)",
        freq_min: 1740.0,
        freq_max: 1780.0,
        intensity: "strong",
        vibration_type: "symm. stretch",
    },
    GroupFrequency {
        group: "N=C=O (isocyanate)",
        freq_min: 2250.0,
        freq_max: 2275.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C=C=O (ketene)",
        freq_min: 2100.0,
        freq_max: 2140.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    // Bending
    GroupFrequency {
        group: "C-H (methyl)",
        freq_min: 1370.0,
        freq_max: 1380.0,
        intensity: "medium",
        vibration_type: "bend",
    },
    GroupFrequency {
        group: "C-H (methylene)",
        freq_min: 1440.0,
        freq_max: 1465.0,
        intensity: "medium",
        vibration_type: "bend",
    },
    GroupFrequency {
        group: "O-H",
        freq_min: 1200.0,
        freq_max: 1440.0,
        intensity: "strong",
        vibration_type: "bend",
    },
    GroupFrequency {
        group: "N-H",
        freq_min: 1550.0,
        freq_max: 1640.0,
        intensity: "medium",
        vibration_type: "bend",
    },
    GroupFrequency {
        group: "H-O-H bend",
        freq_min: 1590.0,
        freq_max: 1665.0,
        intensity: "strong",
        vibration_type: "bend",
    },
    // Single bond stretches
    GroupFrequency {
        group: "C-O (alcohol)",
        freq_min: 1000.0,
        freq_max: 1100.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C-O (ether)",
        freq_min: 1070.0,
        freq_max: 1150.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C-N",
        freq_min: 1020.0,
        freq_max: 1250.0,
        intensity: "medium",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C-F",
        freq_min: 1000.0,
        freq_max: 1350.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C-Cl",
        freq_min: 600.0,
        freq_max: 800.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C-Br",
        freq_min: 500.0,
        freq_max: 680.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "C-I",
        freq_min: 485.0,
        freq_max: 600.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    // Fingerprint region
    GroupFrequency {
        group: "S=O (sulfoxide)",
        freq_min: 1030.0,
        freq_max: 1070.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "S=O (sulfone)",
        freq_min: 1120.0,
        freq_max: 1160.0,
        intensity: "strong",
        vibration_type: "asymm. stretch",
    },
    GroupFrequency {
        group: "P=O",
        freq_min: 1100.0,
        freq_max: 1300.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "P-O-C",
        freq_min: 900.0,
        freq_max: 1050.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    // Metal-ligand and coordination bands
    GroupFrequency {
        group: "M-H",
        freq_min: 1600.0,
        freq_max: 2100.0,
        intensity: "medium",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "M-CO (bridging carbonyl)",
        freq_min: 1750.0,
        freq_max: 1850.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "M-CO (terminal carbonyl)",
        freq_min: 1850.0,
        freq_max: 2120.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "M=O (terminal oxo)",
        freq_min: 850.0,
        freq_max: 1030.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "M-C",
        freq_min: 400.0,
        freq_max: 650.0,
        intensity: "medium",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "M-N",
        freq_min: 400.0,
        freq_max: 650.0,
        intensity: "medium",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "M-O",
        freq_min: 350.0,
        freq_max: 700.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "M-P",
        freq_min: 350.0,
        freq_max: 500.0,
        intensity: "medium",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "M-S",
        freq_min: 250.0,
        freq_max: 500.0,
        intensity: "medium",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "M-Cl",
        freq_min: 250.0,
        freq_max: 400.0,
        intensity: "strong",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "M-Br",
        freq_min: 180.0,
        freq_max: 320.0,
        intensity: "medium",
        vibration_type: "stretch",
    },
    GroupFrequency {
        group: "M-I",
        freq_min: 150.0,
        freq_max: 260.0,
        intensity: "medium",
        vibration_type: "stretch",
    },
];

/// Assign peaks from vibrational analysis to functional groups.
pub fn assign_peaks(
    frequencies: &[f64],
    intensities: &[f64],
    elements: &[u8],
    mode_displacements: Option<&[Vec<[f64; 3]>]>,
) -> AssignmentResult {
    let mut assignments = Vec::new();
    let mut found_groups = std::collections::BTreeSet::new();

    for (i, &freq) in frequencies.iter().enumerate() {
        if freq < 150.0 {
            continue; // Skip lattice-like modes below the useful far-IR window
        }

        let intensity = if i < intensities.len() {
            intensities[i]
        } else {
            0.0
        };

        let mut matches = Vec::new();

        for gf in GROUP_FREQUENCIES {
            if freq >= gf.freq_min && freq <= gf.freq_max {
                // Check if the molecule contains elements consistent with this group
                let mode_matches = mode_displacements
                    .and_then(|all_modes| all_modes.get(i))
                    .map(|mode| group_matches_mode(gf.group, elements, mode))
                    .unwrap_or(true);

                if group_is_possible(gf.group, elements) && mode_matches {
                    let center = (gf.freq_min + gf.freq_max) / 2.0;
                    let range = gf.freq_max - gf.freq_min;
                    let match_quality = 1.0 - (freq - center).abs() / (range / 2.0 + 1.0);

                    matches.push(FunctionalGroupMatch {
                        group: gf.group.to_string(),
                        vibration_type: gf.vibration_type.to_string(),
                        match_quality: match_quality.clamp(0.0, 1.0),
                    });

                    found_groups.insert(gf.group.to_string());
                }
            }
        }

        let confidence = if matches.is_empty() {
            0.0
        } else {
            matches.iter().map(|m| m.match_quality).fold(0.0, f64::max)
        };

        assignments.push(PeakAssignment {
            frequency: freq,
            intensity,
            assignments: matches,
            confidence,
        });
    }

    AssignmentResult {
        assignments,
        functional_groups: found_groups.into_iter().collect(),
    }
}

fn dominant_mode_elements(
    elements: &[u8],
    mode_displacement: &[[f64; 3]],
) -> std::collections::BTreeSet<u8> {
    let amplitudes: Vec<f64> = mode_displacement
        .iter()
        .map(|vector| {
            (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]).sqrt()
        })
        .collect();
    let max_amplitude = amplitudes.iter().copied().fold(0.0_f64, f64::max);
    if max_amplitude <= 1e-8 {
        return std::collections::BTreeSet::new();
    }

    amplitudes
        .iter()
        .enumerate()
        .filter(|(_, amplitude)| **amplitude >= max_amplitude * 0.35)
        .filter_map(|(index, _)| elements.get(index).copied())
        .collect()
}

fn group_matches_mode(group: &str, elements: &[u8], mode_displacement: &[[f64; 3]]) -> bool {
    let dominant = dominant_mode_elements(elements, mode_displacement);
    if dominant.is_empty() {
        return true;
    }

    let has = |z: u8| dominant.contains(&z);
    let has_metal = dominant.iter().copied().any(is_metal_atomic_number);

    if group.starts_with("M-") || group.starts_with("M=") {
        if !has_metal {
            return false;
        }
        if group.contains("-H") {
            return has(1);
        }
        if group.contains("-CO") {
            return has(6) && has(8);
        }
        if group.contains("=O") || group.contains("-O") {
            return has(8);
        }
        if group.contains("-N") {
            return has(7);
        }
        if group.contains("-S") {
            return has(16);
        }
        if group.contains("-P") {
            return has(15);
        }
        if group.contains("-Cl") {
            return has(17);
        }
        if group.contains("-Br") {
            return has(35);
        }
        if group.contains("-I") {
            return has(53);
        }
        if group == "M-C" {
            return has(6);
        }
        return true;
    }

    if group.contains("O-H") || group == "H-O-H bend" {
        return has(8) && has(1);
    }
    if group.contains("N-H") {
        return has(7) && has(1);
    }
    if group.contains("C-H") {
        return has(6) && has(1);
    }
    if group.contains("C=O") {
        return has(6) && has(8);
    }
    if group.contains("C-O") {
        return has(6) && has(8);
    }
    if group.contains("C-N") || group.contains("C=N") {
        return has(6) && has(7);
    }
    if group.contains("C≡N") {
        return has(6) && has(7);
    }
    if group.contains("C≡C") || group.contains("C=C") {
        return has(6);
    }
    if group.contains("S=O") {
        return has(16) && has(8);
    }
    if group.contains("P=O") {
        return has(15) && has(8);
    }
    if group.contains("P-O-C") {
        return has(15) && has(8) && has(6);
    }
    if group.contains("C-F") {
        return has(6) && has(9);
    }
    if group.contains("C-Cl") {
        return has(6) && has(17);
    }
    if group.contains("C-Br") {
        return has(6) && has(35);
    }
    if group.contains("C-I") {
        return has(6) && has(53);
    }

    true
}

/// Check if a functional group is possible given the elements present.
fn group_is_possible(group: &str, elements: &[u8]) -> bool {
    let has = |z: u8| elements.contains(&z);
    let has_metal = elements.iter().copied().any(is_metal_atomic_number);

    if group.starts_with("M-") || group.starts_with("M=") {
        if !has_metal {
            return false;
        }
        if group.contains("-H") {
            return has(1);
        }
        if group.contains("-CO") {
            return has(6) && has(8);
        }
        if group.contains("=O") || group.contains("-O") {
            return has(8);
        }
        if group.contains("-N") {
            return has(7);
        }
        if group.contains("-S") {
            return has(16);
        }
        if group.contains("-P") {
            return has(15);
        }
        if group.contains("-Cl") {
            return has(17);
        }
        if group.contains("-Br") {
            return has(35);
        }
        if group.contains("-I") {
            return has(53);
        }
        if group == "M-C" {
            return has(6);
        }

        return true;
    }

    if group.contains("O-H") || group.contains("C-O") || group.contains("C=O") {
        return has(8); // Oxygen
    }
    if group.contains("N-H")
        || group.contains("C-N")
        || group.contains("C=N")
        || group.contains("C≡N")
        || group.contains("N=O")
    {
        return has(7); // Nitrogen
    }
    if group.contains("C-F") {
        return has(9); // Fluorine
    }
    if group.contains("C-Cl") {
        return has(17); // Chlorine
    }
    if group.contains("C-Br") {
        return has(35); // Bromine
    }
    if group.contains("C-I") {
        return has(53); // Iodine
    }
    if group.contains("S=O") {
        return has(16); // Sulfur
    }
    if group.contains("P=O") || group.contains("P-O") {
        return has(15); // Phosphorus
    }
    if group.contains("C-H") || group.contains("C=C") || group.contains("C≡C") {
        return has(6); // Carbon (almost always true)
    }

    true // Default: possible
}

fn is_metal_atomic_number(z: u8) -> bool {
    matches!(z, 3 | 4 | 11 | 12 | 13 | 19 | 20 | 21..=31 | 37..=50 | 55..=84 | 87..=103)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_peak_assignment_carbonyl() {
        let freqs = vec![1720.0, 2950.0, 1450.0];
        let intensities = vec![1.0, 0.5, 0.3];
        let elements = vec![6u8, 6, 8, 1, 1, 1, 1]; // acetaldehyde-like

        let result = assign_peaks(&freqs, &intensities, &elements, None);
        assert!(!result.assignments.is_empty());

        // 1720 cm⁻¹ should match C=O
        let carbonyl = &result.assignments[0];
        assert!(carbonyl.assignments.iter().any(|a| a.group.contains("C=O")));
    }

    #[test]
    fn test_group_possible() {
        assert!(group_is_possible("O-H (free)", &[6, 8, 1]));
        assert!(!group_is_possible("C-Cl", &[6, 8, 1]));
        assert!(group_is_possible("C-H (sp3)", &[6, 1]));
    }

    #[test]
    fn test_peak_assignment_metal_ligand_bands() {
        let freqs = vec![295.0, 520.0, 1985.0];
        let intensities = vec![0.4, 0.7, 1.0];
        let elements = vec![26u8, 17, 8, 6];

        let result = assign_peaks(&freqs, &intensities, &elements, None);

        assert!(result.assignments[0]
            .assignments
            .iter()
            .any(|a| a.group == "M-Cl"));
        assert!(result.assignments[1]
            .assignments
            .iter()
            .any(|a| a.group == "M-O"));
        assert!(result.assignments[2]
            .assignments
            .iter()
            .any(|a| a.group == "M-CO (terminal carbonyl)"));
    }
}
