//! ReaxFF force-field file parser.
//!
//! Parses the standard `ffield.reax.*` layout: general parameters, atom
//! entries, bonds, off-diagonal terms, angles, torsions and hydrogen bonds.
//! The resulting [`Ffield`] keeps the raw published values (typed by section)
//! so future energy terms can consume them without re-reading the file.
//!
//! # Scope
//!
//! This is a syntactic parser: it does not decide which sections the current
//! energy model uses. The existing C/H/O energy path consumes the transcribed
//! pair table in [`super::params`]; the remaining sections are preserved for
//! the analytic-gradient work.

use std::path::Path;

/// One atom entry: element symbol plus the 31 numeric fields of the standard
/// four-line atom record.
#[derive(Debug, Clone, PartialEq)]
pub struct FfieldAtom {
    /// Element symbol as written in the file.
    pub symbol: String,
    /// Numeric fields after the symbol (31 values in the standard layout).
    pub values: Vec<f64>,
}

/// One non-atom record: leading integer atom types followed by the numeric
/// parameters of that section.
#[derive(Debug, Clone, PartialEq)]
pub struct FfieldRecord {
    /// Atom types at the start of the record.
    pub atom_types: Vec<usize>,
    /// Numeric parameters after the atom types.
    pub values: Vec<f64>,
}

/// Parsed ReaxFF force-field file.
#[derive(Debug, Clone, PartialEq)]
pub struct Ffield {
    /// General parameters (39 values in the published C/H/O file).
    pub general: Vec<f64>,
    /// Atom entries.
    pub atoms: Vec<FfieldAtom>,
    /// Bond entries (at1, at2 + 16 parameters).
    pub bonds: Vec<FfieldRecord>,
    /// Off-diagonal entries (at1, at2 + 6 parameters).
    pub off_diagonal: Vec<FfieldRecord>,
    /// Angle entries (at1, at2, at3 + 7 parameters).
    pub angles: Vec<FfieldRecord>,
    /// Torsion entries (at1, at2, at3, at4 + 7 parameters).
    pub torsions: Vec<FfieldRecord>,
    /// Hydrogen-bond entries (at1, at2, at3 + 4 parameters).
    pub hydrogen_bonds: Vec<FfieldRecord>,
}

impl Ffield {
    /// Parse a ReaxFF force-field file from disk.
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self, String> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .map_err(|err| format!("cannot read {}: {err}", path.display()))?;
        Self::parse(&text)
    }

    /// Parse a ReaxFF force-field file from a string.
    ///
    /// Comments start at `!` and run to the end of the line. The header
    /// (`DATE:`, `CONTRIBUTOR:`, ...) is skipped; the first standalone integer
    /// is the general-parameter count.
    pub fn parse(text: &str) -> Result<Self, String> {
        let tokens = tokenize(text);
        let mut pos = 0usize;
        while pos < tokens.len() {
            if tokens[pos].parse::<usize>().is_ok() {
                break;
            }
            pos += 1;
        }
        if pos >= tokens.len() {
            return Err("no section counts found in ffield".to_string());
        }

        let n_general = take_usize(&tokens, &mut pos, "general parameter count")?;
        let general = take_floats(&tokens, &mut pos, n_general, "general parameters")?;

        let n_atoms = take_usize(&tokens, &mut pos, "atom count")?;
        let mut atoms = Vec::with_capacity(n_atoms);
        for index in 0..n_atoms {
            let symbol = tokens
                .get(pos)
                .ok_or_else(|| format!("truncated atom {index}"))?
                .clone();
            if symbol.parse::<f64>().is_ok() {
                return Err(format!("atom {index} is missing its element symbol"));
            }
            pos += 1;
            let values = take_floats(&tokens, &mut pos, 32, &format!("atom {index}"))?;
            atoms.push(FfieldAtom { symbol, values });
        }

        let bonds = take_records(&tokens, &mut pos, 2, 16, "bond")?;
        let off_diagonal = take_records(&tokens, &mut pos, 2, 6, "off-diagonal")?;
        let angles = take_records(&tokens, &mut pos, 3, 7, "angle")?;
        let torsions = take_records(&tokens, &mut pos, 4, 7, "torsion")?;
        let hydrogen_bonds = take_records(&tokens, &mut pos, 3, 4, "hydrogen-bond")?;

        Ok(Self {
            general,
            atoms,
            bonds,
            off_diagonal,
            angles,
            torsions,
            hydrogen_bonds,
        })
    }
}

fn tokenize(text: &str) -> Vec<String> {
    text.lines()
        .map(|line| line.split('!').next().unwrap_or(""))
        .flat_map(str::split_whitespace)
        // The standard atom section carries descriptive column headers
        // without a leading `!`; drop those tokens, not the element symbols.
        .filter(|token| !token.contains(';') && !token.contains('(') && !token.contains(')'))
        .map(str::to_string)
        .collect()
}

fn take_usize(tokens: &[String], pos: &mut usize, what: &str) -> Result<usize, String> {
    let token = tokens
        .get(*pos)
        .ok_or_else(|| format!("truncated {what}"))?;
    *pos += 1;
    token
        .parse::<usize>()
        .map_err(|_| format!("invalid {what}: '{token}'"))
}

fn take_floats(
    tokens: &[String],
    pos: &mut usize,
    count: usize,
    what: &str,
) -> Result<Vec<f64>, String> {
    if *pos + count > tokens.len() {
        return Err(format!("truncated {what}"));
    }
    let mut values = Vec::with_capacity(count);
    for token in &tokens[*pos..*pos + count] {
        values.push(
            token
                .parse::<f64>()
                .map_err(|_| format!("invalid number in {what}: '{token}'"))?,
        );
    }
    *pos += count;
    Ok(values)
}

fn take_records(
    tokens: &[String],
    pos: &mut usize,
    type_count: usize,
    value_count: usize,
    what: &str,
) -> Result<Vec<FfieldRecord>, String> {
    let count = take_usize(tokens, pos, &format!("{what} count"))?;
    let mut records = Vec::with_capacity(count);
    for index in 0..count {
        let mut atom_types = Vec::with_capacity(type_count);
        for _ in 0..type_count {
            atom_types.push(take_usize(
                tokens,
                pos,
                &format!("{what} {index} atom type"),
            )?);
        }
        let values = take_floats(tokens, pos, value_count, &format!("{what} {index}"))?;
        records.push(FfieldRecord { atom_types, values });
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[allow(clippy::vec_init_then_push)]
    fn sample() -> String {
        let mut lines = Vec::new();
        lines.push("header: DATE: 2026-11-08".to_string());
        lines.push("2 ! Nr of general parameters".to_string());
        lines.push("1.5 !p1".to_string());
        lines.push("-2.5 !p2".to_string());
        lines.push("1 ! Nr of atoms".to_string());
        lines.push("C 1.0 2.0".to_string());
        let mut atom = vec!["3.0".to_string()];
        for i in 4..=32 {
            atom.push(format!("{i}.0"));
        }
        lines.push(atom.join(" "));
        lines.push("1 ! Nr of bonds".to_string());
        lines.push("1 2 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0".to_string());
        lines.push("18.0 19.0 20.0 21.0 22.0 23.0 24.0 25.0".to_string());
        lines.push("1 ! Nr of off-diagonal terms".to_string());
        lines.push("1 2 1.0 2.0 3.0 4.0 5.0 6.0".to_string());
        lines.push("1 ! Nr of angles".to_string());
        lines.push("1 2 1 1.0 2.0 3.0 4.0 5.0 6.0 7.0".to_string());
        lines.push("1 ! Nr of torsions".to_string());
        lines.push("1 2 1 2 1.0 2.0 3.0 4.0 5.0 6.0 7.0".to_string());
        lines.push("1 ! Nr of hydrogen bonds".to_string());
        lines.push("1 2 1 1.0 2.0 3.0 4.0".to_string());
        lines.join("\n")
    }

    #[test]
    fn parses_every_section() {
        let ffield = Ffield::parse(&sample()).expect("parse");
        assert_eq!(ffield.general, vec![1.5, -2.5]);
        assert_eq!(ffield.atoms.len(), 1);
        assert_eq!(ffield.atoms[0].symbol, "C");
        assert_eq!(ffield.atoms[0].values.len(), 32);
        assert_eq!(ffield.bonds.len(), 1);
        assert_eq!(ffield.bonds[0].atom_types, vec![1, 2]);
        assert_eq!(ffield.bonds[0].values.len(), 16);
        assert_eq!(ffield.off_diagonal.len(), 1);
        assert_eq!(ffield.angles.len(), 1);
        assert_eq!(ffield.torsions.len(), 1);
        assert_eq!(ffield.hydrogen_bonds.len(), 1);
        assert_eq!(ffield.hydrogen_bonds[0].values, vec![1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn rejects_truncated_files() {
        let error = Ffield::parse("1\n2.0\n1\nC 1.0\n").unwrap_err();
        assert!(error.contains("truncated"), "{error}");
    }

    #[test]
    fn parses_reference_ffield_when_available() {
        let path = std::env::var("SCIFORM_FFIELD")
            .ok()
            .map(std::path::PathBuf::from)
            .or_else(|| {
                let candidate = std::path::PathBuf::from(
                    ".cache/refdata/lammps-src/potentials/ffield.reax.cho",
                );
                candidate.exists().then_some(candidate)
            });
        let Some(path) = path else {
            eprintln!("reference ffield not available; skipping");
            return;
        };
        let ffield = Ffield::from_path(&path).expect("parse reference ffield");
        assert_eq!(ffield.general.len(), 39);
        assert_eq!(ffield.atoms.len(), 3);
        assert_eq!(ffield.bonds.len(), 6);
        assert_eq!(ffield.off_diagonal.len(), 3);
        assert_eq!(ffield.angles.len(), 18);
        assert_eq!(ffield.torsions.len(), 26);
        assert_eq!(ffield.hydrogen_bonds.len(), 1);
        assert_eq!(ffield.atoms[0].symbol, "C");
        assert_eq!(ffield.atoms[1].symbol, "H");
        assert_eq!(ffield.atoms[2].symbol, "O");
    }
}
