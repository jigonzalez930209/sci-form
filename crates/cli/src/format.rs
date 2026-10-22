//! Output formatting and shared parse helpers.
//!
//! # Error message style
//!
//! Every user-facing error follows one template (blocker C-H9):
//!
//! ```text
//! <What failed> '<actual value>'. Supported: 'a', 'b', 'c'.
//! ```
//!
//! All literals are quoted, the offending input is always included and the
//! valid values are always listed when there is a fixed set.  Rust internal
//! names never leak into a CLI message.
//!
//! # Optional-argument convention
//!
//! An optional argument is **absent**, never an empty string, an empty JSON
//! array or a magic number (blocker C-M5).  `--flag` omitted means "not
//! given"; the engine default applies.

/// Format the canonical "invalid value" error message.
///
/// `<what> '<value>'. Supported: 'a', 'b', 'c'.`
pub fn invalid_value(what: &str, value: &str, supported: &[&str]) -> String {
    if supported.is_empty() {
        return format!("Invalid {what} '{value}'.");
    }
    let list = supported
        .iter()
        .map(|s| format!("'{s}'"))
        .collect::<Vec<_>>()
        .join(", ");
    format!("Invalid {what} '{value}'. Supported: {list}.")
}

/// Truncate a long string for inclusion in an error message.
pub fn truncate(value: &str, max_chars: usize) -> String {
    let trimmed = value.trim();
    if trimmed.chars().count() <= max_chars {
        return trimmed.to_string();
    }
    let head: String = trimmed.chars().take(max_chars).collect();
    format!("{head}…")
}

/// Format the canonical "malformed JSON argument" error message.
pub fn bad_json(what: &str, value: &str, err: &serde_json::Error) -> String {
    let fragment = truncate(value, 80);
    format!("Invalid {what} JSON '{fragment}'. ({err})")
}

/// Print an error message to stderr and exit with a non-zero status.
pub fn fail(message: impl AsRef<str>) -> ! {
    eprintln!("{}", message.as_ref());
    std::process::exit(1);
}

pub fn format_xyz(result: &sci_form::ConformerResult) -> String {
    if result.error.is_some() || result.coords.is_empty() {
        return format!(
            "0\nERROR: {}\n",
            result.error.as_deref().unwrap_or("unknown")
        );
    }
    let n = result.num_atoms;
    let mut out = format!("{}\n{}\n", n, result.smiles);
    for i in 0..n {
        let elem = element_symbol(result.elements[i]);
        let (x, y, z) = (
            result.coords[i * 3],
            result.coords[i * 3 + 1],
            result.coords[i * 3 + 2],
        );
        out.push_str(&format!("{:2} {:12.6} {:12.6} {:12.6}\n", elem, x, y, z));
    }
    out
}

pub fn format_sdf(result: &sci_form::ConformerResult) -> String {
    if result.error.is_some() || result.coords.is_empty() {
        return format!(
            "\n  sci-form\n\nERROR: {}\n$$$$\n",
            result.error.as_deref().unwrap_or("unknown")
        );
    }
    let n = result.num_atoms;
    let nb = result.bonds.len();
    let mut out = format!("{}\n  sci-form\n\n", result.smiles);
    out.push_str(&format!(
        "{:3}{:3}  0  0  0  0  0  0  0  0999 V2000\n",
        n, nb
    ));
    for i in 0..n {
        let (x, y, z) = (
            result.coords[i * 3],
            result.coords[i * 3 + 1],
            result.coords[i * 3 + 2],
        );
        let elem = element_symbol(result.elements[i]);
        out.push_str(&format!(
            "{:10.4}{:10.4}{:10.4} {:3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
            x, y, z, elem
        ));
    }
    for (a, b, order) in &result.bonds {
        let bond_type = match order.as_str() {
            "SINGLE" => 1,
            "DOUBLE" => 2,
            "TRIPLE" => 3,
            "AROMATIC" => 4,
            _ => 1,
        };
        out.push_str(&format!("{:3}{:3}{:3}  0\n", a + 1, b + 1, bond_type));
    }
    out.push_str("M  END\n$$$$\n");
    out
}

/// IUPAC element symbol for an atomic number (1–118).
///
/// Returns `"X"` for atomic numbers outside the periodic table; the previous
/// 13-element table produced invalid SDF/XYZ output for every metal.
pub fn element_symbol(atomic_num: u8) -> &'static str {
    const SYMBOLS: [&str; 119] = [
        "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P",
        "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
        "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re",
        "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db",
        "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
    ];
    SYMBOLS.get(atomic_num as usize).copied().unwrap_or("X")
}

/// Parse JSON elements + coords, exit on error.
pub fn parse_elems_coords(elements: &str, coords: &str) -> (Vec<u8>, Vec<f64>, Vec<[f64; 3]>) {
    let elems: Vec<u8> = serde_json::from_str(elements).unwrap_or_else(|e| {
        fail(bad_json("'elements'", elements, &e));
    });
    let flat: Vec<f64> = serde_json::from_str(coords).unwrap_or_else(|e| {
        fail(bad_json("'coords'", coords, &e));
    });
    let positions: Vec<[f64; 3]> = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    (elems, flat, positions)
}

/// Perceive covalent bonds from 3D coordinates using covalent radii.
///
/// Used by the commands that need a bond list (GETAWAY, hapticity) or
/// partial charges (RDF) but only receive a geometry.
pub fn perceive_bonds(
    elements: &[u8],
    positions: &[[f64; 3]],
    tolerance: f64,
) -> Vec<(usize, usize)> {
    let n = elements.len().min(positions.len());
    let mut bonds = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            let r_sum = sci_form::graph::get_covalent_radius(elements[i])
                + sci_form::graph::get_covalent_radius(elements[j]);
            if d > 0.4 && d < r_sum + tolerance {
                bonds.push((i, j));
            }
        }
    }
    bonds
}

/// Gasteiger-Marsili charges derived from geometry-perceived bonds.
pub fn gasteiger_charges_from_geometry(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<Vec<f64>, String> {
    let bonds = perceive_bonds(elements, positions, 0.3);
    let formal_charges = vec![0i8; elements.len()];
    let result = sci_form::charges::gasteiger::gasteiger_marsili_charges(
        elements,
        &bonds,
        &formal_charges,
        6,
    )?;
    Ok(result.charges)
}

/// Parse an optional JSON array of `[i, j]` bond pairs.
///
/// `None` means "not given": the bonds are perceived from the geometry.
pub fn parse_bond_list(bonds: Option<&str>) -> Option<Vec<(usize, usize)>> {
    let raw = bonds?;
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    Some(serde_json::from_str(trimmed).unwrap_or_else(|e| {
        fail(bad_json("'bonds'", trimmed, &e));
    }))
}

/// Parse an optional JSON array of numbers.
///
/// `None` means "not given".
pub fn parse_f64_list(values: Option<&str>, what: &str) -> Option<Vec<f64>> {
    let raw = values?;
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    Some(serde_json::from_str(trimmed).unwrap_or_else(|e| {
        fail(bad_json(what, trimmed, &e));
    }))
}

/// Parse an optional JSON array of integers.
///
/// `None` means "not given".
pub fn parse_usize_list(values: Option<&str>, what: &str) -> Option<Vec<usize>> {
    let raw = values?;
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    Some(serde_json::from_str(trimmed).unwrap_or_else(|e| {
        fail(bad_json(what, trimmed, &e));
    }))
}

/// Parse an optional JSON array of strings.
///
/// `None` means "not given": every available backend is used.
#[allow(dead_code)]
pub fn parse_string_list(values: Option<&str>, what: &str) -> Option<Vec<String>> {
    let raw = values?;
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    Some(serde_json::from_str(trimmed).unwrap_or_else(|e| {
        fail(bad_json(what, trimmed, &e));
    }))
}
