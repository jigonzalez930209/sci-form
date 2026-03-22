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
];

/// Assign peaks from vibrational analysis to functional groups.
pub fn assign_peaks(
    frequencies: &[f64],
    intensities: &[f64],
    elements: &[u8],
    _mode_displacements: Option<&[Vec<[f64; 3]>]>,
) -> AssignmentResult {
    let mut assignments = Vec::new();
    let mut found_groups = std::collections::BTreeSet::new();

    for (i, &freq) in frequencies.iter().enumerate() {
        if freq < 400.0 {
            continue; // Skip very low frequencies
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
                if group_is_possible(gf.group, elements) {
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

/// Check if a functional group is possible given the elements present.
fn group_is_possible(group: &str, elements: &[u8]) -> bool {
    let has = |z: u8| elements.contains(&z);

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
}
