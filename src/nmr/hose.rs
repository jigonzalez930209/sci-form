//! HOSE (Hierarchically Ordered Spherical description of Environment) code generation.
//!
//! Generates unique environment descriptors for each atom by traversing
//! the molecular graph in concentric spheres (radius 0 to max_radius).
//! Includes a compile-time HOSE→chemical-shift lookup database for
//! ¹H and ¹³C prediction with fallback from higher to lower radius.

use crate::graph::{BondOrder, Molecule};
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

/// A HOSE code descriptor for a specific atom.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HoseCode {
    /// Atom index in the molecule.
    pub atom_index: usize,
    /// Atomic number of the center atom.
    pub element: u8,
    /// HOSE string at each sphere radius (0..=max_radius).
    pub spheres: Vec<String>,
    /// Full concatenated HOSE code.
    pub full_code: String,
}

fn bond_symbol(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Single => "",
        BondOrder::Double => "=",
        BondOrder::Triple => "#",
        BondOrder::Aromatic => "*",
        BondOrder::Unknown => "",
    }
}

fn element_symbol(z: u8) -> &'static str {
    match z {
        1 => "H",
        2 => "He",
        3 => "Li",
        4 => "Be",
        5 => "B",
        6 => "C",
        7 => "N",
        8 => "O",
        9 => "F",
        11 => "Na",
        12 => "Mg",
        13 => "Al",
        14 => "Si",
        15 => "P",
        16 => "S",
        17 => "Cl",
        19 => "K",
        20 => "Ca",
        21 => "Sc",
        22 => "Ti",
        23 => "V",
        24 => "Cr",
        25 => "Mn",
        26 => "Fe",
        27 => "Co",
        28 => "Ni",
        29 => "Cu",
        30 => "Zn",
        31 => "Ga",
        32 => "Ge",
        33 => "As",
        34 => "Se",
        35 => "Br",
        37 => "Rb",
        38 => "Sr",
        40 => "Zr",
        41 => "Nb",
        42 => "Mo",
        44 => "Ru",
        45 => "Rh",
        46 => "Pd",
        47 => "Ag",
        48 => "Cd",
        49 => "In",
        50 => "Sn",
        51 => "Sb",
        52 => "Te",
        53 => "I",
        54 => "Xe",
        55 => "Cs",
        56 => "Ba",
        74 => "W",
        78 => "Pt",
        79 => "Au",
        80 => "Hg",
        81 => "Tl",
        82 => "Pb",
        83 => "Bi",
        _ => "X",
    }
}

/// Generate HOSE codes for all atoms in a molecule.
///
/// `mol`: parsed molecular graph
/// `max_radius`: maximum sphere radius (typically 3–5)
pub fn generate_hose_codes(mol: &Molecule, max_radius: usize) -> Vec<HoseCode> {
    let n = mol.graph.node_count();
    let mut codes = Vec::with_capacity(n);

    for atom_idx in 0..n {
        let center = NodeIndex::new(atom_idx);
        let center_element = mol.graph[center].element;
        let center_sym = element_symbol(center_element);

        let mut spheres = Vec::with_capacity(max_radius + 1);
        spheres.push(center_sym.to_string());

        // BFS by sphere depth
        let mut visited = vec![false; n];
        visited[atom_idx] = true;
        let mut current_frontier: Vec<(NodeIndex, BondOrder)> = Vec::new();

        // Initial frontier: direct neighbors
        for edge in mol.graph.edges(center) {
            let neighbor = if edge.source() == center {
                edge.target()
            } else {
                edge.source()
            };
            current_frontier.push((neighbor, edge.weight().order));
        }

        for _radius in 1..=max_radius {
            if current_frontier.is_empty() {
                spheres.push(String::new());
                continue;
            }

            // Sort frontier entries for deterministic output
            let mut sphere_parts: BTreeSet<String> = BTreeSet::new();
            let mut next_frontier: Vec<(NodeIndex, BondOrder)> = Vec::new();

            for (node, bond_order) in &current_frontier {
                let node_idx = node.index();
                if visited[node_idx] {
                    continue;
                }
                visited[node_idx] = true;

                let elem = mol.graph[*node].element;
                let sym = element_symbol(elem);
                let bond_sym = bond_symbol(*bond_order);
                sphere_parts.insert(format!("{}{}", bond_sym, sym));

                // Collect next frontier
                for edge in mol.graph.edges(*node) {
                    let next = if edge.source() == *node {
                        edge.target()
                    } else {
                        edge.source()
                    };
                    if !visited[next.index()] {
                        next_frontier.push((next, edge.weight().order));
                    }
                }
            }

            spheres.push(sphere_parts.into_iter().collect::<Vec<_>>().join(","));
            current_frontier = next_frontier;
        }

        let full_code = format!("{}/{}", spheres[0], spheres[1..].join("/"));

        codes.push(HoseCode {
            atom_index: atom_idx,
            element: center_element,
            spheres,
            full_code,
        });
    }

    codes
}

// ─── HOSE → Chemical Shift Database ─────────────────────────────────────────
//
// Empirical database mapping HOSE code prefixes to chemical shifts (ppm).
// Based on NMRShiftDB2 reference data patterns for common organic environments.
// Uses fallback: try radius 4 → 3 → 2 → 1 match.

/// Result of a HOSE-based chemical shift lookup.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HoseShiftLookup {
    /// Atom index.
    pub atom_index: usize,
    /// Element (atomic number).
    pub element: u8,
    /// Predicted shift in ppm.
    pub shift_ppm: f64,
    /// HOSE code used for matching.
    pub matched_hose: String,
    /// Radius at which the match was found.
    pub match_radius: usize,
    /// Confidence (higher radius = higher confidence).
    pub confidence: f64,
}

/// ¹H HOSE shift database: (HOSE prefix pattern, shift_ppm).
/// Keys are simplified HOSE sphere-1 environment strings.
fn h1_hose_database() -> Vec<(&'static str, f64)> {
    vec![
        // sp3 C-H environments
        ("C/H,H,H,C", 0.90),  // methyl next to alkyl
        ("C/H,H,C,C", 1.30),  // methylene
        ("C/H,C,C,C", 1.50),  // methine
        ("C/H,H,H,O", 3.40),  // methyl next to O (methoxy)
        ("C/H,H,O", 3.50),    // CH₂ next to O
        ("C/H,H,H,N", 2.30),  // methyl next to N
        ("C/H,H,N", 2.60),    // CH₂ next to N
        ("C/H,H,H,=O", 2.10), // methyl next to carbonyl
        ("C/H,H,=O", 2.50),   // methylene next to carbonyl
        ("C/H,H,H,*C", 2.30), // methyl on aromatic
        // sp2 C-H (aromatic)
        ("C/*C,*C,H", 7.27),  // benzene-type ArH
        ("C/*C,*C,*C", 7.50), // fused aromatic ArH
        ("C/*C,*N,H", 7.80),  // pyridine-adjacent ArH
        // sp2 C-H (alkene)
        ("C/=C,H,H", 5.25),  // vinyl =CH₂
        ("C/=C,H,C", 5.40),  // internal alkene
        ("C/=C,=C,H", 6.30), // conjugated diene
        // sp C-H (alkyne)
        ("C/#C,H", 2.50), // terminal alkyne
        // Aldehyde
        ("C/=O,H,C", 9.50), // aldehyde CHO
        ("C/=O,H,H", 9.60), // formaldehyde
        // O-H
        ("O/H,C", 2.50),  // alcohol OH
        ("O/H,*C", 5.50), // phenol OH
        // N-H
        ("N/H,H,C", 1.50), // primary amine
        ("N/H,C,C", 2.20), // secondary amine
    ]
}

/// ¹³C HOSE shift database.
fn c13_hose_database() -> Vec<(&'static str, f64)> {
    vec![
        // sp3 carbon
        ("C/H,H,H,C", 15.0),  // methyl
        ("C/H,H,C,C", 25.0),  // methylene
        ("C/H,C,C,C", 35.0),  // methine
        ("C/C,C,C,C", 40.0),  // quaternary
        ("C/H,H,H,O", 55.0),  // methoxy
        ("C/H,H,O,C", 65.0),  // C-O-C
        ("C/H,H,H,N", 32.0),  // N-methyl
        ("C/H,H,N", 45.0),    // C-N
        ("C/H,H,H,=O", 30.0), // methyl next to carbonyl
        // sp2 carbon (aromatic)
        ("C/*C,*C,H", 128.0),  // unsubstituted ArC
        ("C/*C,*C,C", 137.0),  // alkyl-substituted ArC
        ("C/*C,*C,O", 155.0),  // O-substituted ArC
        ("C/*C,*C,N", 148.0),  // N-substituted ArC
        ("C/*C,*C,F", 163.0),  // F-substituted ArC
        ("C/*C,*C,Cl", 134.0), // Cl-substituted ArC
        // sp2 carbon (alkene)
        ("C/=C,H,H", 115.0), // vinyl =CH₂
        ("C/=C,H,C", 130.0), // internal alkene
        ("C/=C,C,C", 140.0), // trisubstituted alkene
        // Carbonyl
        ("C/=O,O,C", 175.0), // ester / carboxylic acid
        ("C/=O,N,C", 170.0), // amide
        ("C/=O,C,C", 205.0), // ketone
        ("C/=O,H,C", 200.0), // aldehyde
        // sp carbon
        ("C/#C,H", 70.0), // terminal alkyne
        ("C/#C,C", 85.0), // internal alkyne
    ]
}

/// Predict chemical shift from HOSE code using database lookup with fallback.
///
/// Searches the HOSE code against the database, trying the most specific
/// match first (full code), then falling back to shorter prefixes.
pub fn predict_shift_from_hose(hose_code: &HoseCode, nucleus: u8) -> Option<HoseShiftLookup> {
    let database = match nucleus {
        1 => h1_hose_database(),
        6 => c13_hose_database(),
        _ => return None,
    };

    // Try matching from highest radius down to 1
    for radius in (1..=hose_code.spheres.len().saturating_sub(1)).rev() {
        // Build the HOSE prefix up to this radius
        let prefix = format!(
            "{}/{}",
            hose_code.spheres[0],
            hose_code.spheres[1..=radius].join("/")
        );

        for &(pattern, shift) in &database {
            if prefix.contains(pattern)
                || pattern.contains(&prefix)
                || fuzzy_hose_match(&prefix, pattern)
            {
                return Some(HoseShiftLookup {
                    atom_index: hose_code.atom_index,
                    element: nucleus,
                    shift_ppm: shift,
                    matched_hose: pattern.to_string(),
                    match_radius: radius,
                    confidence: 0.5 + 0.1 * radius as f64,
                });
            }
        }
    }

    None
}

// ─── External HOSE database (NMRShiftDB2 / SDBS exports) ────────────────────

/// One entry of an external HOSE→shift database.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HoseDatabaseEntry {
    /// Full HOSE code string, e.g. `C/H,H,H,C`.
    pub hose: String,
    /// Element (atomic number) the shift refers to.
    pub element: u8,
    /// Chemical shift in ppm.
    pub shift_ppm: f64,
    /// Optional number of experimental records behind this entry.
    #[serde(default)]
    pub count: usize,
}

/// External HOSE→chemical-shift database.
///
/// The crate does not bundle NMRShiftDB2 or SDBS data (their licenses require
/// separate review).  This type reads a database exported by the user and
/// performs exact longest-prefix lookup; the built-in tables remain the
/// default when no external database is supplied.
///
/// Supported formats:
///
/// * JSON: either an array of [`HoseDatabaseEntry`] or an object with an
///   `entries` array.
/// * CSV/TSV: `hose,element,shift_ppm[,count]`, with `#` comment lines and an
///   optional header.  `element` can be an atomic number or a symbol.
#[derive(Debug, Clone, Default)]
pub struct HoseDatabase {
    entries: Vec<HoseDatabaseEntry>,
}

impl HoseDatabase {
    /// Create an empty database.
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of stored records.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// True when no record has been loaded.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Add one record.
    pub fn push(&mut self, entry: HoseDatabaseEntry) {
        self.entries.push(entry);
    }

    /// Load a database from a JSON or CSV/TSV file.
    ///
    /// The format is selected from the file extension; unknown extensions are
    /// tried as JSON first and then as CSV.
    pub fn from_path(path: impl AsRef<std::path::Path>) -> Result<Self, String> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .map_err(|error| format!("cannot read {}: {error}", path.display()))?;
        match path
            .extension()
            .and_then(|extension| extension.to_str())
            .map(str::to_ascii_lowercase)
            .as_deref()
        {
            Some("json") => Self::from_json(&text),
            Some("csv") | Some("tsv") => Self::from_csv(&text),
            _ => Self::from_json(&text).or_else(|_| Self::from_csv(&text)),
        }
    }

    /// Parse the JSON export format.
    pub fn from_json(text: &str) -> Result<Self, String> {
        #[derive(Deserialize)]
        #[serde(untagged)]
        enum Payload {
            Array(Vec<HoseDatabaseEntry>),
            Wrapped { entries: Vec<HoseDatabaseEntry> },
        }
        let payload: Payload =
            serde_json::from_str(text).map_err(|error| format!("invalid HOSE JSON: {error}"))?;
        let entries = match payload {
            Payload::Array(entries) => entries,
            Payload::Wrapped { entries } => entries,
        };
        Ok(Self { entries })
    }

    /// Parse a CSV/TSV export (`hose,element,shift_ppm[,count]`).
    pub fn from_csv(text: &str) -> Result<Self, String> {
        let mut entries = Vec::new();
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let tab_separated = line.contains('\t');
            let fields: Vec<&str> = if tab_separated {
                line.split('\t').map(str::trim).collect()
            } else {
                // HOSE codes themselves contain commas, so parse the numeric
                // tail from the right and keep the remaining fields as the code.
                line.split(',').map(str::trim).collect()
            };
            if fields.len() < 3 || fields[0].eq_ignore_ascii_case("hose") {
                continue;
            }

            let (hose, element_field, shift_field, count_field) = if tab_separated {
                (
                    fields[0].to_string(),
                    fields[1],
                    fields[2],
                    fields.get(3).copied(),
                )
            } else {
                let last = fields.len() - 1;
                let last_is_count = fields.len() >= 4
                    && fields[last].parse::<usize>().is_ok()
                    && fields[last - 1].parse::<f64>().is_ok();
                if last_is_count {
                    (
                        fields[..last - 2].join(","),
                        fields[last - 2],
                        fields[last - 1],
                        Some(fields[last]),
                    )
                } else {
                    (
                        fields[..last - 1].join(","),
                        fields[last - 1],
                        fields[last],
                        None,
                    )
                }
            };

            let element = parse_element(element_field)?;
            let shift_ppm: f64 = shift_field
                .parse()
                .map_err(|error| format!("invalid shift '{shift_field}': {error}"))?;
            let count = count_field
                .and_then(|value| value.parse().ok())
                .unwrap_or(1);
            entries.push(HoseDatabaseEntry {
                hose,
                element,
                shift_ppm,
                count,
            });
        }
        Ok(Self { entries })
    }

    /// Look up a HOSE code, trying the full code and then successively shorter
    /// sphere prefixes.  Returns the most specific match.
    pub fn lookup(&self, hose_code: &HoseCode, nucleus: u8) -> Option<HoseShiftLookup> {
        for radius in (1..=hose_code.spheres.len().saturating_sub(1)).rev() {
            let prefix = if radius == 0 {
                hose_code.spheres[0].clone()
            } else {
                format!(
                    "{}/{}",
                    hose_code.spheres[0],
                    hose_code.spheres[1..=radius].join("/")
                )
            };
            if let Some(entry) = self
                .entries
                .iter()
                .find(|entry| entry.element == nucleus && entry.hose == prefix)
            {
                return Some(HoseShiftLookup {
                    atom_index: hose_code.atom_index,
                    element: nucleus,
                    shift_ppm: entry.shift_ppm,
                    matched_hose: entry.hose.clone(),
                    match_radius: radius,
                    confidence: 0.5 + 0.1 * radius as f64,
                });
            }
        }
        None
    }
}

fn parse_element(value: &str) -> Result<u8, String> {
    if let Ok(z) = value.parse::<u8>() {
        return Ok(z);
    }
    match value {
        "H" => Ok(1),
        "C" => Ok(6),
        "N" => Ok(7),
        "O" => Ok(8),
        "F" => Ok(9),
        "Si" => Ok(14),
        "P" => Ok(15),
        "S" => Ok(16),
        "Cl" => Ok(17),
        "Se" => Ok(34),
        "Br" => Ok(35),
        "I" => Ok(53),
        _ => Err(format!("unknown element symbol '{value}'")),
    }
}

/// Predict a shift from an external HOSE database.
pub fn predict_shift_from_database(
    hose_code: &HoseCode,
    nucleus: u8,
    database: &HoseDatabase,
) -> Option<HoseShiftLookup> {
    database.lookup(hose_code, nucleus)
}

/// Fuzzy matching: check if the sphere-1 fragment of the HOSE code
/// matches the sphere-1 fragment of the database pattern.
fn fuzzy_hose_match(hose: &str, pattern: &str) -> bool {
    let hose_parts: Vec<&str> = hose.split('/').collect();
    let pat_parts: Vec<&str> = pattern.split('/').collect();

    if hose_parts.len() < 2 || pat_parts.len() < 2 {
        return false;
    }

    // Match center element
    if hose_parts[0] != pat_parts[0] {
        return false;
    }

    // Fuzzy match sphere-1: check if sorted neighbor sets overlap
    let hose_neighbors: BTreeSet<&str> = hose_parts[1].split(',').collect();
    let pat_neighbors: BTreeSet<&str> = pat_parts[1].split(',').collect();

    let intersection = hose_neighbors.intersection(&pat_neighbors).count();
    let union = hose_neighbors.union(&pat_neighbors).count();

    if union == 0 {
        return false;
    }

    // Require at least 50% overlap
    (intersection as f64 / union as f64) > 0.5
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hose_codes_ethanol() {
        let mol = Molecule::from_smiles("CCO").unwrap();
        let codes = generate_hose_codes(&mol, 3);

        // Should have codes for all atoms (including H)
        assert_eq!(codes.len(), mol.graph.node_count());

        // All codes should have non-empty first sphere
        for code in &codes {
            assert!(!code.spheres[0].is_empty());
            assert!(!code.full_code.is_empty());
        }
    }

    #[test]
    fn test_hose_codes_benzene() {
        let mol = Molecule::from_smiles("c1ccccc1").unwrap();
        let codes = generate_hose_codes(&mol, 3);

        assert_eq!(codes.len(), mol.graph.node_count());

        // All carbon atoms should have similar HOSE codes due to symmetry
        let carbon_codes: Vec<&HoseCode> = codes.iter().filter(|c| c.element == 6).collect();
        assert!(!carbon_codes.is_empty());
    }

    #[test]
    fn test_hose_codes_deterministic() {
        let mol = Molecule::from_smiles("CC(=O)O").unwrap();
        let codes1 = generate_hose_codes(&mol, 3);
        let codes2 = generate_hose_codes(&mol, 3);

        for (c1, c2) in codes1.iter().zip(codes2.iter()) {
            assert_eq!(c1.full_code, c2.full_code);
        }
    }

    #[test]
    fn external_hose_database_json_and_prefix_fallback() {
        let mol = Molecule::from_smiles("CCO").unwrap();
        let codes = generate_hose_codes(&mol, 3);
        let carbon = codes.iter().find(|code| code.element == 6).unwrap();
        let full = carbon.full_code.clone();
        let json = format!(r#"[{{"hose":"{full}","element":6,"shift_ppm":42.0,"count":3}}]"#);
        let db = HoseDatabase::from_json(&json).unwrap();
        assert_eq!(db.len(), 1);
        let hit = db.lookup(carbon, 6).expect("exact HOSE match");
        assert_eq!(hit.shift_ppm, 42.0);
        assert_eq!(hit.match_radius, 3);

        // A database entry at sphere 1 must still match through the fallback.
        let sphere1 = format!("{}/{}", carbon.spheres[0], carbon.spheres[1]);
        let db = HoseDatabase::from_json(&format!(
            r#"[{{"hose":"{sphere1}","element":6,"shift_ppm":17.5}}]"#
        ))
        .unwrap();
        let hit = db.lookup(carbon, 6).expect("sphere-1 fallback");
        assert_eq!(hit.shift_ppm, 17.5);
        assert_eq!(hit.match_radius, 1);
    }

    #[test]
    fn external_hose_database_csv() {
        let csv = "# NMRShiftDB2 export\nhose,element,shift_ppm,count\nC/H,H,H,C,1,0.9,10\nC/*C,*C,H,C,7.27,25\n";
        let db = HoseDatabase::from_csv(csv).unwrap();
        assert_eq!(db.len(), 2);
        assert_eq!(db.entries[0].element, 1);
        assert_eq!(db.entries[0].shift_ppm, 0.9);
        assert_eq!(db.entries[1].element, 6);
        assert_eq!(db.entries[1].count, 25);
    }
}
