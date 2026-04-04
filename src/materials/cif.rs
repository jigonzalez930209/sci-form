//! CIF (Crystallographic Information File) import and export.
//!
//! Supports reading and writing CIF 1.1 format for crystal structures:
//! - Cell parameters (`_cell_length_a`, `_cell_angle_alpha`, etc.)
//! - Space group (`_symmetry_space_group_name_H-M`, `_symmetry_Int_Tables_number`)
//! - Atom sites (`_atom_site_label`, `_atom_site_fract_x/y/z`, `_atom_site_type_symbol`)

use super::cell::{CellParameters, UnitCell};

/// A parsed CIF crystal structure.
#[derive(Debug, Clone)]
pub struct CifStructure {
    /// Data block name (from `data_` line).
    pub data_block: String,
    /// Unit cell.
    pub cell: UnitCell,
    /// Cell parameters (a, b, c, α, β, γ).
    pub cell_params: CellParameters,
    /// Space group Hermann-Mauguin symbol (if present).
    pub space_group_hm: Option<String>,
    /// Space group ITC number (if present).
    pub space_group_number: Option<u16>,
    /// Atom sites with labels, elements, and fractional coordinates.
    pub atom_sites: Vec<CifAtomSite>,
}

/// A single atom site from a CIF file.
#[derive(Debug, Clone)]
pub struct CifAtomSite {
    /// Atom label (e.g., "Na1", "Cl2").
    pub label: String,
    /// Element symbol (e.g., "Na", "Cl").
    pub type_symbol: String,
    /// Atomic number.
    pub atomic_number: u8,
    /// Fractional coordinates.
    pub frac_x: f64,
    pub frac_y: f64,
    pub frac_z: f64,
    /// Occupancy (default 1.0).
    pub occupancy: f64,
}

/// Parse a CIF-format string into a crystal structure.
///
/// Handles standard CIF 1.1 tags for cell parameters, space group, and atom sites.
/// Numbers with uncertainty notation (e.g., "5.640(1)") are stripped of the parenthetical.
pub fn parse_cif(input: &str) -> Result<CifStructure, String> {
    let mut data_block = String::new();
    let mut a = None;
    let mut b = None;
    let mut c = None;
    let mut alpha = None;
    let mut beta = None;
    let mut gamma = None;
    let mut sg_hm: Option<String> = None;
    let mut sg_number: Option<u16> = None;
    let mut atom_sites = Vec::new();

    let lines: Vec<&str> = input.lines().collect();
    let mut i = 0;

    while i < lines.len() {
        let line = lines[i].trim();

        // Skip comments and empty lines
        if line.is_empty() || line.starts_with('#') {
            i += 1;
            continue;
        }

        // Data block
        if let Some(name) = line.strip_prefix("data_") {
            data_block = name.trim().to_string();
            i += 1;
            continue;
        }

        // Cell parameters — single-value tags
        if line.starts_with("_cell_length_a") {
            a = parse_cif_float(extract_value(line, "_cell_length_a"));
            i += 1;
            continue;
        }
        if line.starts_with("_cell_length_b") {
            b = parse_cif_float(extract_value(line, "_cell_length_b"));
            i += 1;
            continue;
        }
        if line.starts_with("_cell_length_c") {
            c = parse_cif_float(extract_value(line, "_cell_length_c"));
            i += 1;
            continue;
        }
        if line.starts_with("_cell_angle_alpha") {
            alpha = parse_cif_float(extract_value(line, "_cell_angle_alpha"));
            i += 1;
            continue;
        }
        if line.starts_with("_cell_angle_beta") {
            beta = parse_cif_float(extract_value(line, "_cell_angle_beta"));
            i += 1;
            continue;
        }
        if line.starts_with("_cell_angle_gamma") {
            gamma = parse_cif_float(extract_value(line, "_cell_angle_gamma"));
            i += 1;
            continue;
        }

        // Space group
        if line.starts_with("_symmetry_space_group_name_H-M")
            || line.starts_with("_space_group_name_H-M_alt")
        {
            let val = extract_value_quoted(line);
            if !val.is_empty() {
                sg_hm = Some(val);
            }
            i += 1;
            continue;
        }
        if line.starts_with("_symmetry_Int_Tables_number")
            || line.starts_with("_space_group_IT_number")
        {
            let val = extract_value_simple(line);
            sg_number = val.parse::<u16>().ok();
            i += 1;
            continue;
        }

        // Atom site loop
        if line == "loop_" {
            // Check if next lines define atom_site columns
            let mut col_names = Vec::new();
            let mut j = i + 1;
            while j < lines.len() {
                let col_line = lines[j].trim();
                if col_line.starts_with("_atom_site_") {
                    col_names.push(col_line.to_string());
                    j += 1;
                } else {
                    break;
                }
            }

            if !col_names.is_empty() {
                // Map column indices
                let label_col = col_names.iter().position(|c| c == "_atom_site_label");
                let type_col = col_names.iter().position(|c| c == "_atom_site_type_symbol");
                let fx_col = col_names.iter().position(|c| c == "_atom_site_fract_x");
                let fy_col = col_names.iter().position(|c| c == "_atom_site_fract_y");
                let fz_col = col_names.iter().position(|c| c == "_atom_site_fract_z");
                let occ_col = col_names.iter().position(|c| c == "_atom_site_occupancy");
                let cx_col = col_names.iter().position(|c| c == "_atom_site_Cartn_x");
                let cy_col = col_names.iter().position(|c| c == "_atom_site_Cartn_y");
                let cz_col = col_names.iter().position(|c| c == "_atom_site_Cartn_z");

                let has_frac = fx_col.is_some() && fy_col.is_some() && fz_col.is_some();
                let has_cart = cx_col.is_some() && cy_col.is_some() && cz_col.is_some();

                // Read data rows
                while j < lines.len() {
                    let data_line = lines[j].trim();
                    if data_line.is_empty()
                        || data_line.starts_with('#')
                        || data_line.starts_with('_')
                        || data_line == "loop_"
                    {
                        break;
                    }

                    let fields: Vec<&str> = split_cif_fields(data_line);
                    if fields.len() < col_names.len() {
                        j += 1;
                        continue;
                    }

                    let label = label_col.map(|c| fields[c].to_string()).unwrap_or_default();
                    let type_sym = type_col
                        .map(|c| fields[c].to_string())
                        .or_else(|| {
                            // Extract element from label (e.g., "Na1" → "Na")
                            Some(extract_element_from_label(&label))
                        })
                        .unwrap_or_default();

                    let (fx, fy, fz) = if has_frac {
                        (
                            parse_cif_float(fields[fx_col.unwrap()]).unwrap_or(0.0),
                            parse_cif_float(fields[fy_col.unwrap()]).unwrap_or(0.0),
                            parse_cif_float(fields[fz_col.unwrap()]).unwrap_or(0.0),
                        )
                    } else if has_cart {
                        // Will convert below after cell is built
                        (
                            parse_cif_float(fields[cx_col.unwrap()]).unwrap_or(0.0),
                            parse_cif_float(fields[cy_col.unwrap()]).unwrap_or(0.0),
                            parse_cif_float(fields[cz_col.unwrap()]).unwrap_or(0.0),
                        )
                    } else {
                        (0.0, 0.0, 0.0)
                    };

                    let occ = occ_col
                        .and_then(|c| parse_cif_float(fields[c]))
                        .unwrap_or(1.0);

                    let z = symbol_to_atomic_number(&type_sym);

                    atom_sites.push((
                        CifAtomSite {
                            label,
                            type_symbol: type_sym,
                            atomic_number: z,
                            frac_x: fx,
                            frac_y: fy,
                            frac_z: fz,
                            occupancy: occ,
                        },
                        has_cart && !has_frac,
                    ));

                    j += 1;
                }

                i = j;
                continue;
            }

            i += 1;
            continue;
        }

        i += 1;
    }

    // Validate cell parameters
    let a = a.ok_or("Missing _cell_length_a")?;
    let b = b.ok_or("Missing _cell_length_b")?;
    let c = c.ok_or("Missing _cell_length_c")?;
    let alpha = alpha.unwrap_or(90.0);
    let beta = beta.unwrap_or(90.0);
    let gamma = gamma.unwrap_or(90.0);

    let cell_params = CellParameters {
        a,
        b,
        c,
        alpha,
        beta,
        gamma,
    };
    let cell = UnitCell::from_parameters(&cell_params);

    // Convert Cartesian atom sites to fractional if needed
    let atom_sites: Vec<CifAtomSite> = atom_sites
        .into_iter()
        .map(|(mut site, is_cart)| {
            if is_cart {
                let frac = cell.cart_to_frac([site.frac_x, site.frac_y, site.frac_z]);
                site.frac_x = frac[0];
                site.frac_y = frac[1];
                site.frac_z = frac[2];
            }
            site
        })
        .collect();

    Ok(CifStructure {
        data_block,
        cell,
        cell_params,
        space_group_hm: sg_hm,
        space_group_number: sg_number,
        atom_sites,
    })
}

/// Write a CIF-format string from a crystal structure.
pub fn write_cif(structure: &CifStructure) -> String {
    let mut out = String::new();

    out.push_str(&format!("data_{}\n", structure.data_block));
    out.push_str("#\n");

    // Cell parameters
    let p = &structure.cell_params;
    out.push_str(&format!("_cell_length_a                    {:.6}\n", p.a));
    out.push_str(&format!("_cell_length_b                    {:.6}\n", p.b));
    out.push_str(&format!("_cell_length_c                    {:.6}\n", p.c));
    out.push_str(&format!(
        "_cell_angle_alpha                 {:.4}\n",
        p.alpha
    ));
    out.push_str(&format!(
        "_cell_angle_beta                  {:.4}\n",
        p.beta
    ));
    out.push_str(&format!(
        "_cell_angle_gamma                 {:.4}\n",
        p.gamma
    ));
    out.push_str("#\n");

    // Space group
    if let Some(ref hm) = structure.space_group_hm {
        out.push_str(&format!("_symmetry_space_group_name_H-M    '{}'\n", hm));
    }
    if let Some(num) = structure.space_group_number {
        out.push_str(&format!("_symmetry_Int_Tables_number       {}\n", num));
    }
    out.push_str("#\n");

    // Atom sites
    if !structure.atom_sites.is_empty() {
        out.push_str("loop_\n");
        out.push_str("_atom_site_label\n");
        out.push_str("_atom_site_type_symbol\n");
        out.push_str("_atom_site_fract_x\n");
        out.push_str("_atom_site_fract_y\n");
        out.push_str("_atom_site_fract_z\n");
        out.push_str("_atom_site_occupancy\n");
        for site in &structure.atom_sites {
            out.push_str(&format!(
                "{:<8} {:<4} {:.6} {:.6} {:.6} {:.4}\n",
                site.label, site.type_symbol, site.frac_x, site.frac_y, site.frac_z, site.occupancy
            ));
        }
    }

    out.push_str("#END\n");
    out
}

// ─── Helpers ──────────────────────────────────────────────────────────────────

/// Extract value part from a CIF tag line: `_tag  value` → `value`.
fn extract_value<'a>(line: &'a str, tag: &str) -> &'a str {
    line[tag.len()..].trim()
}

/// Extract simple (non-quoted) value from a tag line.
fn extract_value_simple(line: &str) -> &str {
    // Find first whitespace after the tag
    if let Some(pos) = line.find(char::is_whitespace) {
        line[pos..].trim()
    } else {
        ""
    }
}

/// Extract possibly-quoted value (strip single quotes).
fn extract_value_quoted(line: &str) -> String {
    let val = extract_value_simple(line);
    val.trim_matches('\'').trim_matches('"').trim().to_string()
}

/// Parse a CIF float, stripping uncertainty notation like "5.640(1)" → 5.640.
fn parse_cif_float(s: &str) -> Option<f64> {
    let s = s.trim();
    if s == "?" || s == "." {
        return None;
    }
    // Strip parenthetical uncertainty: "5.640(1)" → "5.640"
    let clean = if let Some(paren) = s.find('(') {
        &s[..paren]
    } else {
        s
    };
    clean.parse::<f64>().ok()
}

/// Split a CIF data line into fields, respecting single-quoted strings.
fn split_cif_fields(line: &str) -> Vec<&str> {
    let mut fields = Vec::new();
    let bytes = line.as_bytes();
    let mut i = 0;
    let len = bytes.len();

    while i < len {
        // Skip whitespace
        while i < len && bytes[i] == b' ' {
            i += 1;
        }
        if i >= len {
            break;
        }

        if bytes[i] == b'\'' {
            // Quoted field
            let start = i + 1;
            i = start;
            while i < len && bytes[i] != b'\'' {
                i += 1;
            }
            fields.push(&line[start..i]);
            if i < len {
                i += 1; // skip closing quote
            }
        } else {
            // Unquoted field
            let start = i;
            while i < len && bytes[i] != b' ' {
                i += 1;
            }
            fields.push(&line[start..i]);
        }
    }

    fields
}

/// Extract element symbol from an atom label (e.g., "Na1" → "Na", "Fe2+" → "Fe").
fn extract_element_from_label(label: &str) -> String {
    let mut sym = String::new();
    let mut chars = label.chars();
    if let Some(c) = chars.next() {
        if c.is_ascii_uppercase() {
            sym.push(c);
            if let Some(c2) = chars.next() {
                if c2.is_ascii_lowercase() {
                    sym.push(c2);
                }
            }
        }
    }
    if sym.is_empty() {
        "X".to_string()
    } else {
        sym
    }
}

/// Convert element symbol to atomic number.
fn symbol_to_atomic_number(sym: &str) -> u8 {
    match sym.trim() {
        "H" => 1,
        "He" => 2,
        "Li" => 3,
        "Be" => 4,
        "B" => 5,
        "C" => 6,
        "N" => 7,
        "O" => 8,
        "F" => 9,
        "Ne" => 10,
        "Na" => 11,
        "Mg" => 12,
        "Al" => 13,
        "Si" => 14,
        "P" => 15,
        "S" => 16,
        "Cl" => 17,
        "Ar" => 18,
        "K" => 19,
        "Ca" => 20,
        "Sc" => 21,
        "Ti" => 22,
        "V" => 23,
        "Cr" => 24,
        "Mn" => 25,
        "Fe" => 26,
        "Co" => 27,
        "Ni" => 28,
        "Cu" => 29,
        "Zn" => 30,
        "Ga" => 31,
        "Ge" => 32,
        "As" => 33,
        "Se" => 34,
        "Br" => 35,
        "Kr" => 36,
        "Rb" => 37,
        "Sr" => 38,
        "Y" => 39,
        "Zr" => 40,
        "Nb" => 41,
        "Mo" => 42,
        "Tc" => 43,
        "Ru" => 44,
        "Rh" => 45,
        "Pd" => 46,
        "Ag" => 47,
        "Cd" => 48,
        "In" => 49,
        "Sn" => 50,
        "Sb" => 51,
        "Te" => 52,
        "I" => 53,
        "Xe" => 54,
        "Cs" => 55,
        "Ba" => 56,
        "La" => 57,
        "Ce" => 58,
        "Pr" => 59,
        "Nd" => 60,
        "Pm" => 61,
        "Sm" => 62,
        "Eu" => 63,
        "Gd" => 64,
        "Tb" => 65,
        "Dy" => 66,
        "Ho" => 67,
        "Er" => 68,
        "Tm" => 69,
        "Yb" => 70,
        "Lu" => 71,
        "Hf" => 72,
        "Ta" => 73,
        "W" => 74,
        "Re" => 75,
        "Os" => 76,
        "Ir" => 77,
        "Pt" => 78,
        "Au" => 79,
        "Hg" => 80,
        "Tl" => 81,
        "Pb" => 82,
        "Bi" => 83,
        _ => 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_nacl_cif() {
        let cif = r#"
data_NaCl
_cell_length_a                    5.640
_cell_length_b                    5.640
_cell_length_c                    5.640
_cell_angle_alpha                 90.000
_cell_angle_beta                  90.000
_cell_angle_gamma                 90.000
_symmetry_space_group_name_H-M    'F m -3 m'
_symmetry_Int_Tables_number       225
#
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Na1      Na   0.000000 0.000000 0.000000 1.0000
Cl1      Cl   0.500000 0.500000 0.500000 1.0000
"#;
        let structure = parse_cif(cif).unwrap();
        assert_eq!(structure.data_block, "NaCl");
        assert!((structure.cell_params.a - 5.640).abs() < 1e-6);
        assert_eq!(structure.space_group_hm, Some("F m -3 m".to_string()));
        assert_eq!(structure.space_group_number, Some(225));
        assert_eq!(structure.atom_sites.len(), 2);
        assert_eq!(structure.atom_sites[0].type_symbol, "Na");
        assert_eq!(structure.atom_sites[0].atomic_number, 11);
        assert_eq!(structure.atom_sites[1].type_symbol, "Cl");
        assert_eq!(structure.atom_sites[1].atomic_number, 17);
        assert!((structure.atom_sites[1].frac_x - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_roundtrip_cif() {
        let original = CifStructure {
            data_block: "test".to_string(),
            cell: UnitCell::from_parameters(&CellParameters {
                a: 10.0,
                b: 10.0,
                c: 10.0,
                alpha: 90.0,
                beta: 90.0,
                gamma: 90.0,
            }),
            cell_params: CellParameters {
                a: 10.0,
                b: 10.0,
                c: 10.0,
                alpha: 90.0,
                beta: 90.0,
                gamma: 90.0,
            },
            space_group_hm: Some("P m -3 m".to_string()),
            space_group_number: Some(221),
            atom_sites: vec![CifAtomSite {
                label: "Fe1".to_string(),
                type_symbol: "Fe".to_string(),
                atomic_number: 26,
                frac_x: 0.0,
                frac_y: 0.0,
                frac_z: 0.0,
                occupancy: 1.0,
            }],
        };

        let cif_text = write_cif(&original);
        let parsed = parse_cif(&cif_text).unwrap();
        assert_eq!(parsed.data_block, "test");
        assert!((parsed.cell_params.a - 10.0).abs() < 1e-4);
        assert_eq!(parsed.atom_sites.len(), 1);
        assert_eq!(parsed.atom_sites[0].type_symbol, "Fe");
        assert_eq!(parsed.atom_sites[0].atomic_number, 26);
    }

    #[test]
    fn test_parse_cif_float_uncertainty() {
        assert_eq!(parse_cif_float("5.640(1)"), Some(5.640));
        assert_eq!(parse_cif_float("90.000"), Some(90.0));
        assert_eq!(parse_cif_float("?"), None);
    }

    #[test]
    fn test_extract_element_from_label() {
        assert_eq!(extract_element_from_label("Na1"), "Na");
        assert_eq!(extract_element_from_label("Fe2+"), "Fe");
        assert_eq!(extract_element_from_label("O3"), "O");
        assert_eq!(extract_element_from_label("C12"), "C");
    }
}
