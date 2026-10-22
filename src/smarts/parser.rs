//! SMARTS parser: converts a SMARTS string into a pattern graph.
//!
//! The parser covers the subset of SMARTS used by the CSD torsion-preference
//! patterns: atom primitives (element, aromaticity, hydrogen count, degree,
//! ring membership and size, formal charge, hybridisation, chirality and
//! recursive `$(...)` queries), bond primitives (single, double, triple,
//! aromatic, any, ring/not-ring) and the logical combinators AND, OR and
//! NOT.
//!
//! The result is a small pattern graph ([`SmartsPattern`]) of atom queries
//! and bond queries, which [`crate::smarts::matcher`] then matches against a
//! molecule.

use crate::graph::ChiralType;

/// A parsed SMARTS pattern as a small graph of atom/bond queries.
///
/// The atoms and bonds are stored in parallel vectors; a bond refers to its
/// endpoints by index into [`SmartsPattern::atoms`].
#[derive(Debug, Clone)]
pub struct SmartsPattern {
    /// The atom queries of the pattern, in the order they were written.
    pub atoms: Vec<SmartsAtom>,
    /// The bond queries of the pattern.
    pub bonds: Vec<SmartsBond>,
}

/// One atom of a SMARTS pattern.
#[derive(Debug, Clone)]
pub struct SmartsAtom {
    /// The primitive (or combination of primitives) the atom must satisfy.
    pub query: AtomQuery,
    /// The SMARTS atom-map number (`:N`) if the pattern declares one.
    pub map_idx: Option<u8>,
}

/// One bond of a SMARTS pattern.
#[derive(Debug, Clone)]
pub struct SmartsBond {
    /// Index of the first atom in [`SmartsPattern::atoms`].
    pub from: usize,
    /// Index of the second atom in [`SmartsPattern::atoms`].
    pub to: usize,
    /// The primitive (or combination of primitives) the bond must satisfy.
    pub query: BondQuery,
}

/// An atom-level SMARTS primitive.
///
/// The combinator variants [`AtomQuery::And`], [`AtomQuery::Or`] and
/// [`AtomQuery::Not`] mirror the SMARTS logical operators; the primitive
/// variants correspond one-to-one to the SMARTS atom primitives and carry
/// the argument of that primitive (an atomic number, a count, a ring size,
/// …).
#[derive(Debug, Clone)]
// The variant set may grow in a future release; downstream crates must
// match with a wildcard arm.
#[non_exhaustive]
pub enum AtomQuery {
    /// Matches any atom; the SMARTS `*` primitive.
    True,
    /// Matches an aliphatic atom of the given element, by atomic number
    /// (SMARTS `C`, `N`, `Cl`, …).
    Element(u8),
    /// Matches an aromatic atom of the given element, by atomic number
    /// (SMARTS `c`, `n`, `o`, `s`, `p`, `se`, `as`).
    AromaticElem(u8),
    /// Matches any aromatic atom (SMARTS `a`).
    AnyAromatic,
    /// Matches any aliphatic atom (SMARTS `A`).
    AnyAliphatic,
    /// Matches an atom of the given atomic number regardless of aromaticity
    /// (SMARTS `#N`).
    AtomicNum(u8),
    /// Matches an atom whose atomic number differs from the given one
    /// (SMARTS `!#N`).
    ///
    /// The parser does not emit this variant — a leading `!` becomes
    /// [`AtomQuery::Not`] — but the matcher understands it.
    NotAtomicNum(u8),
    /// Matches an atom with exactly the given total hydrogen count
    /// (SMARTS `HN`).
    TotalH(u8),
    /// Matches an atom with exactly the given total connection count,
    /// hydrogens included (SMARTS `XN`).
    TotalDegree(u8),
    /// Matches an atom with exactly the given number of connections to
    /// non-hydrogen atoms (SMARTS `DN`).
    HeavyDegree(u8),
    /// Matches an atom that has exactly the given number of ring bonds
    /// (SMARTS `xN`).
    RingBondCount(u8),
    /// Matches an atom that is part of any ring (SMARTS `R`).
    InRing,
    /// Matches an atom that is part of a ring of exactly the given size
    /// (SMARTS `rN`).
    RingSize(u8),
    /// Matches an atom that is part of a ring whose size lies between the
    /// two bounds, inclusive (SMARTS `r{N-M}`).
    RingSizeRange(u8, u8),
    /// Matches an atom that is part of a ring of at least the given size
    /// (SMARTS `r{N-}`).
    RingSizeMin(u8),
    /// Matches an atom with exactly the given formal charge, in units of the
    /// elementary charge (SMARTS `+N` / `-N`).
    FormalCharge(i8),
    /// Matches an atom with the given hybridisation code: 1 = sp, 2 = sp²,
    /// 3 = sp³ (SMARTS `^N`).
    Hybridization(u8),
    /// Matches an atom that belongs to exactly the given number of
    /// smallest-set-of-smallest-rings (SMARTS `RN`).
    RingCount(u8),
    /// Matches an atom with the given tetrahedral configuration
    /// (SMARTS `@` / `@@`).
    Chiral(ChiralType),
    /// Matches an atom that is itself the anchor of an embedded pattern
    /// (SMARTS `$(...)`).
    Recursive(Box<SmartsPattern>),
    /// Matches an atom that satisfies every one of the sub-queries
    /// (juxtaposition or explicit `&`).
    And(Vec<AtomQuery>),
    /// Matches an atom that satisfies at least one of the sub-queries
    /// (SMARTS `,`).
    Or(Vec<AtomQuery>),
    /// Matches an atom that does *not* satisfy the sub-query (SMARTS `!`).
    Not(Box<AtomQuery>),
}

/// A bond-level SMARTS primitive.
///
/// As with [`AtomQuery`], the combinator variants mirror the SMARTS logical
/// operators applied to bond primitives.
#[derive(Debug, Clone)]
// The variant set may grow in a future release; downstream crates must
// match with a wildcard arm.
#[non_exhaustive]
pub enum BondQuery {
    /// Matches a single bond (SMARTS `-`).
    Single,
    /// Matches a double bond (SMARTS `=`).
    Double,
    /// Matches a triple bond (SMARTS `#`).
    Triple,
    /// Matches an aromatic bond (SMARTS `:`).
    Aromatic,
    /// Matches any bond (SMARTS `~`).
    Any,
    /// Matches a bond that is part of a ring (SMARTS `@`).
    Ring,
    /// Matches a bond that is *not* part of a ring (SMARTS `!@`).
    NotRing,
    /// The default bond of an unspecified SMARTS bond: matches a single or
    /// an aromatic bond.
    Implicit,
    /// Matches a bond that satisfies every one of the sub-queries
    /// (several bond primitives juxtaposed, e.g. `!@;-`).
    And(Vec<BondQuery>),
    /// Matches a bond that does *not* satisfy the sub-query (SMARTS `!`).
    Not(Box<BondQuery>),
}

/// Parse a SMARTS string into a SmartsPattern.
pub fn parse_smarts(smarts: &str) -> Result<SmartsPattern, String> {
    let mut parser = SmartsParser::new(smarts);
    parser.parse_chain(None)?;
    Ok(SmartsPattern {
        atoms: parser.atoms,
        bonds: parser.bonds,
    })
}

/// Element symbols indexed by atomic number (H = 1 .. Lr = 103).
const ELEMENT_SYMBOLS: [&str; 103] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
    "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr",
];

/// Atomic number for a case-sensitive element symbol (`Mg` → 12).
fn element_from_symbol(symbol: &str) -> Option<u8> {
    ELEMENT_SYMBOLS
        .iter()
        .position(|&s| s == symbol)
        .map(|idx| (idx + 1) as u8)
}

struct SmartsParser<'a> {
    input: &'a [u8],
    pos: usize,
    atoms: Vec<SmartsAtom>,
    bonds: Vec<SmartsBond>,
    ring_opens: [Option<usize>; 10], // ring closure digits 0-9
}

impl<'a> SmartsParser<'a> {
    fn new(s: &'a str) -> Self {
        Self {
            input: s.as_bytes(),
            pos: 0,
            atoms: Vec::new(),
            bonds: Vec::new(),
            ring_opens: [None; 10],
        }
    }

    fn peek(&self) -> Option<u8> {
        self.input.get(self.pos).copied()
    }

    fn advance(&mut self) -> Option<u8> {
        let c = self.input.get(self.pos).copied();
        if c.is_some() {
            self.pos += 1;
        }
        c
    }

    fn expect(&mut self, ch: u8) -> Result<(), String> {
        if self.advance() == Some(ch) {
            Ok(())
        } else {
            Err(format!("expected '{}' at pos {}", ch as char, self.pos - 1))
        }
    }

    /// Parse a chain of atoms/bonds, optionally connected to `prev_atom`.
    fn parse_chain(&mut self, prev_atom: Option<usize>) -> Result<(), String> {
        let mut prev = prev_atom;
        while self.pos < self.input.len() {
            let c = match self.peek() {
                Some(c) => c,
                None => break,
            };

            match c {
                b')' => break, // end of branch
                b'(' => {
                    // Branch
                    self.advance();
                    self.parse_chain(prev)?;
                    self.expect(b')')?;
                }
                b'[' | b'*' | b'c' | b'n' | b'o' | b's' | b'p' | b'C' | b'N' | b'O' | b'S'
                | b'P' | b'F' | b'B' | b'I' | b'a' | b'A' | b'H' => {
                    // Parse optional bond before atom
                    let bond_q = self.parse_bond_if_present();
                    let atom_idx = self.parse_atom()?;
                    if let Some(p) = prev {
                        self.bonds.push(SmartsBond {
                            from: p,
                            to: atom_idx,
                            query: bond_q.unwrap_or(BondQuery::Implicit),
                        });
                    }
                    prev = Some(atom_idx);
                }
                b'-' | b'=' | b'#' | b'~' | b'/' | b'\\' | b':' | b'!' | b'@' => {
                    // Bond followed by atom
                    let bond_q = self.parse_bond()?;
                    let atom_idx = self.parse_atom()?;
                    if let Some(p) = prev {
                        self.bonds.push(SmartsBond {
                            from: p,
                            to: atom_idx,
                            query: bond_q,
                        });
                    }
                    prev = Some(atom_idx);
                }
                b'0'..=b'9' => {
                    // Ring closure
                    let digit = (self.advance().unwrap() - b'0') as usize;
                    if let Some(open_atom) = self.ring_opens[digit] {
                        self.bonds.push(SmartsBond {
                            from: open_atom,
                            to: prev.unwrap_or(0),
                            query: BondQuery::Implicit,
                        });
                        self.ring_opens[digit] = None;
                    } else {
                        self.ring_opens[digit] = prev;
                    }
                }
                _ => break,
            }
        }
        Ok(())
    }

    /// Try to parse a bond query if one is present (without consuming atom chars).
    fn parse_bond_if_present(&mut self) -> Option<BondQuery> {
        match self.peek() {
            Some(b'-') | Some(b'=') | Some(b'#') | Some(b'~') | Some(b'!') | Some(b'@')
            | Some(b':') => self.parse_bond().ok(),
            _ => None,
        }
    }

    /// Parse a bond query.
    fn parse_bond(&mut self) -> Result<BondQuery, String> {
        let mut parts = Vec::new();
        loop {
            match self.peek() {
                Some(b'-') => {
                    self.advance();
                    parts.push(BondQuery::Single);
                }
                Some(b'=') => {
                    self.advance();
                    parts.push(BondQuery::Double);
                }
                Some(b'#') => {
                    self.advance();
                    parts.push(BondQuery::Triple);
                }
                Some(b'~') => {
                    self.advance();
                    parts.push(BondQuery::Any);
                }
                Some(b':') => {
                    self.advance();
                    parts.push(BondQuery::Aromatic);
                }
                Some(b'@') => {
                    self.advance();
                    parts.push(BondQuery::Ring);
                }
                Some(b'!') => {
                    self.advance();
                    if self.peek() == Some(b'@') {
                        self.advance();
                        parts.push(BondQuery::NotRing);
                    } else {
                        let inner = self.parse_bond()?;
                        parts.push(BondQuery::Not(Box::new(inner)));
                    }
                }
                Some(b';') => {
                    self.advance();
                } // AND separator, continue
                Some(b',') => {
                    // OR — not typically used in torsion SMARTS bonds
                    self.advance();
                }
                _ => break,
            }
        }
        match parts.len() {
            0 => Ok(BondQuery::Implicit),
            1 => Ok(parts.pop().unwrap()),
            _ => Ok(BondQuery::And(parts)),
        }
    }

    /// Parse an atom (either bracket or organic subset).
    fn parse_atom(&mut self) -> Result<usize, String> {
        let atom = match self.peek() {
            Some(b'[') => self.parse_bracket_atom()?,
            Some(b'*') => {
                self.advance();
                SmartsAtom {
                    query: AtomQuery::True,
                    map_idx: None,
                }
            }
            _ => self.parse_organic_atom()?,
        };
        let idx = self.atoms.len();
        self.atoms.push(atom);
        Ok(idx)
    }

    /// Parse an organic subset atom (single uppercase letter, possibly followed by lowercase).
    fn parse_organic_atom(&mut self) -> Result<SmartsAtom, String> {
        let c = self.advance().ok_or("unexpected end")?;
        let query = match c {
            b'C' if self.peek() == Some(b'l') => {
                self.advance();
                AtomQuery::Element(17)
            }
            b'B' if self.peek() == Some(b'r') => {
                self.advance();
                AtomQuery::Element(35)
            }
            b'C' => AtomQuery::Element(6),
            b'N' => AtomQuery::Element(7),
            b'O' => AtomQuery::Element(8),
            b'S' => AtomQuery::Element(16),
            b'P' => AtomQuery::Element(15),
            b'F' => AtomQuery::Element(9),
            b'B' => AtomQuery::Element(5),
            b'I' => AtomQuery::Element(53),
            b'H' => AtomQuery::Element(1),
            b'c' => AtomQuery::AromaticElem(6),
            b'n' => AtomQuery::AromaticElem(7),
            b'o' => AtomQuery::AromaticElem(8),
            b's' => AtomQuery::AromaticElem(16),
            b'p' => AtomQuery::AromaticElem(15),
            b'a' => AtomQuery::AnyAromatic,
            b'A' => AtomQuery::AnyAliphatic,
            _ => {
                return Err(format!(
                    "unexpected atom char '{}' at pos {}",
                    c as char,
                    self.pos - 1
                ))
            }
        };
        Ok(SmartsAtom {
            query,
            map_idx: None,
        })
    }

    /// Parse a bracket atom [...]
    fn parse_bracket_atom(&mut self) -> Result<SmartsAtom, String> {
        self.expect(b'[')?;
        let query = self.parse_atom_spec()?;
        // Check for map class :N before closing bracket
        let map_idx = if self.peek() == Some(b':') {
            self.advance();
            Some(self.parse_number()? as u8)
        } else {
            None
        };
        self.expect(b']')?;
        Ok(SmartsAtom { query, map_idx })
    }

    /// Top-level atom spec: semicolon-separated low-priority AND groups.
    /// Precedence (lowest to highest): ; (low AND) < , (OR) < implicit/& (high AND) < ! (NOT)
    fn parse_atom_spec(&mut self) -> Result<AtomQuery, String> {
        let mut parts = vec![self.parse_atom_query_or()?];
        while self.peek() == Some(b';') {
            self.advance();
            parts.push(self.parse_atom_query_or()?);
        }
        if parts.len() == 1 {
            Ok(parts.pop().unwrap())
        } else {
            Ok(AtomQuery::And(parts))
        }
    }

    /// Parse an atom query with OR (comma-separated).
    fn parse_atom_query_or(&mut self) -> Result<AtomQuery, String> {
        let mut parts = vec![self.parse_atom_query_and()?];
        while self.peek() == Some(b',') {
            self.advance();
            parts.push(self.parse_atom_query_and()?);
        }
        if parts.len() == 1 {
            Ok(parts.pop().unwrap())
        } else {
            Ok(AtomQuery::Or(parts))
        }
    }

    /// Parse an atom query with high-priority AND (implicit juxtaposition or &).
    fn parse_atom_query_and(&mut self) -> Result<AtomQuery, String> {
        let mut parts = Vec::new();
        loop {
            match self.peek() {
                Some(b']') | Some(b',') | Some(b':') | Some(b';') | None => break,
                Some(b'&') => {
                    self.advance();
                } // explicit high-priority AND
                _ => parts.push(self.parse_atom_primitive()?),
            }
        }
        match parts.len() {
            0 => Ok(AtomQuery::True),
            1 => Ok(parts.pop().unwrap()),
            _ => Ok(AtomQuery::And(parts)),
        }
    }

    /// Parse a single atom primitive.
    fn parse_atom_primitive(&mut self) -> Result<AtomQuery, String> {
        let c = self.peek().ok_or("unexpected end in atom spec")?;

        // Two-letter element symbols (`Mg`, `Br`, `As`, …): when the pair
        // resolves to an element, consume it and return the element query;
        // otherwise fall through to the single-letter query handling below.
        if c.is_ascii_uppercase() {
            if let Some(&next) = self.input.get(self.pos + 1) {
                if next.is_ascii_lowercase() {
                    let symbol = [c as char, next as char].iter().collect::<String>();
                    if let Some(z) = element_from_symbol(&symbol) {
                        self.advance();
                        self.advance();
                        return Ok(AtomQuery::Element(z));
                    }
                }
            }
        }
        // Aromatic two-letter elements: `se` (34) and `as` (33).
        if c == b's' && self.input.get(self.pos + 1) == Some(&b'e') {
            self.advance();
            self.advance();
            return Ok(AtomQuery::AromaticElem(34));
        }
        if c == b'a' && self.input.get(self.pos + 1) == Some(&b's') {
            self.advance();
            self.advance();
            return Ok(AtomQuery::AromaticElem(33));
        }

        match c {
            b'!' => {
                self.advance();
                let inner = self.parse_atom_primitive()?;
                Ok(AtomQuery::Not(Box::new(inner)))
            }
            b'#' => {
                self.advance();
                let n = self.parse_number()? as u8;
                Ok(AtomQuery::AtomicNum(n))
            }
            b'@' => {
                self.advance();
                let chiral = if self.peek() == Some(b'@') {
                    self.advance();
                    ChiralType::TetrahedralCW
                } else {
                    ChiralType::TetrahedralCCW
                };
                Ok(AtomQuery::Chiral(chiral))
            }
            b'$' => {
                // Recursive SMARTS: $(smarts)
                self.advance();
                self.expect(b'(')?;
                let start = self.pos;
                // Find matching closing paren, handling nested parens
                let mut depth = 1;
                while depth > 0 && self.pos < self.input.len() {
                    match self.input[self.pos] {
                        b'(' => depth += 1,
                        b')' => depth -= 1,
                        _ => {}
                    }
                    if depth > 0 {
                        self.pos += 1;
                    }
                }
                let inner_str = std::str::from_utf8(&self.input[start..self.pos])
                    .map_err(|_| "invalid utf8 in recursive SMARTS")?;
                self.expect(b')')?;
                let inner = parse_smarts(inner_str)?;
                Ok(AtomQuery::Recursive(Box::new(inner)))
            }
            b'X' => {
                self.advance();
                let n = self.parse_number()? as u8;
                Ok(AtomQuery::TotalDegree(n))
            }
            b'x' => {
                self.advance();
                let n = self.parse_number()? as u8;
                Ok(AtomQuery::RingBondCount(n))
            }
            b'H' => {
                self.advance();
                // H followed by digit = hydrogen count; otherwise H0 is "no H"
                if self.peek().is_some_and(|c| c.is_ascii_digit()) {
                    let n = self.parse_number()? as u8;
                    Ok(AtomQuery::TotalH(n))
                } else {
                    // H without number means H >= 1 (at least one hydrogen)
                    Ok(AtomQuery::TotalH(1))
                }
            }
            b'D' => {
                self.advance();
                if self.peek().is_some_and(|c| c.is_ascii_digit()) {
                    let n = self.parse_number()? as u8;
                    Ok(AtomQuery::HeavyDegree(n))
                } else {
                    Ok(AtomQuery::HeavyDegree(1))
                }
            }
            b'R' => {
                self.advance();
                if self.peek().is_some_and(|c| c.is_ascii_digit()) {
                    let n = self.parse_number()? as u8;
                    Ok(AtomQuery::RingCount(n))
                } else {
                    Ok(AtomQuery::InRing)
                }
            }
            b'r' => {
                self.advance();
                if self.peek() == Some(b'{') {
                    // r{N-M} or r{N-}
                    self.advance(); // '{'
                    if self.peek() == Some(b'-') {
                        // r{-M} → ring size ≤ M
                        self.advance();
                        let m = self.parse_number()? as u8;
                        self.expect(b'}')?;
                        Ok(AtomQuery::RingSizeRange(3, m))
                    } else {
                        let n = self.parse_number()? as u8;
                        if self.peek() == Some(b'-') {
                            self.advance();
                            if self.peek() == Some(b'}') {
                                self.advance();
                                Ok(AtomQuery::RingSizeMin(n))
                            } else {
                                let m = self.parse_number()? as u8;
                                self.expect(b'}')?;
                                Ok(AtomQuery::RingSizeRange(n, m))
                            }
                        } else {
                            self.expect(b'}')?;
                            Ok(AtomQuery::RingSize(n))
                        }
                    }
                } else if self.peek().is_some_and(|c| c.is_ascii_digit()) {
                    let n = self.parse_number()? as u8;
                    Ok(AtomQuery::RingSize(n))
                } else {
                    Ok(AtomQuery::InRing)
                }
            }
            b'+' => {
                self.advance();
                if self.peek().is_some_and(|c| c.is_ascii_digit()) {
                    let n = self.parse_number()? as i8;
                    Ok(AtomQuery::FormalCharge(n))
                } else {
                    Ok(AtomQuery::FormalCharge(1))
                }
            }
            b'-' => {
                // Careful: '-' can also be a bond. Inside bracket atom, it's charge.
                self.advance();
                if self.peek().is_some_and(|c| c.is_ascii_digit()) {
                    let n = self.parse_number()? as i8;
                    Ok(AtomQuery::FormalCharge(-n))
                } else {
                    Ok(AtomQuery::FormalCharge(-1))
                }
            }
            b'^' => {
                self.advance();
                let n = self.parse_number()? as u8;
                Ok(AtomQuery::Hybridization(n))
            }
            b'*' => {
                self.advance();
                Ok(AtomQuery::True)
            }
            b'a' => {
                self.advance();
                Ok(AtomQuery::AnyAromatic)
            }
            b'A' => {
                self.advance();
                // Check if followed by more letters (element symbol like Al, As, etc.)
                // For CSD patterns, 'A' alone means aliphatic
                Ok(AtomQuery::AnyAliphatic)
            }
            b'C' => {
                self.advance();
                if self.peek() == Some(b'l') {
                    self.advance();
                    Ok(AtomQuery::Element(17))
                } else {
                    Ok(AtomQuery::Element(6))
                }
            }
            b'N' => {
                self.advance();
                Ok(AtomQuery::Element(7))
            }
            b'O' => {
                self.advance();
                Ok(AtomQuery::Element(8))
            }
            b'S' => {
                self.advance();
                Ok(AtomQuery::Element(16))
            }
            b'P' => {
                self.advance();
                Ok(AtomQuery::Element(15))
            }
            b'F' => {
                self.advance();
                Ok(AtomQuery::Element(9))
            }
            b'B' => {
                self.advance();
                if self.peek() == Some(b'r') {
                    self.advance();
                    Ok(AtomQuery::Element(35))
                } else {
                    Ok(AtomQuery::Element(5))
                }
            }
            b'I' => {
                self.advance();
                Ok(AtomQuery::Element(53))
            }
            b'c' => {
                self.advance();
                Ok(AtomQuery::AromaticElem(6))
            }
            b'n' => {
                self.advance();
                Ok(AtomQuery::AromaticElem(7))
            }
            b'o' => {
                self.advance();
                Ok(AtomQuery::AromaticElem(8))
            }
            b's' => {
                self.advance();
                Ok(AtomQuery::AromaticElem(16))
            }
            b'p' => {
                self.advance();
                Ok(AtomQuery::AromaticElem(15))
            }
            _ => Err(format!(
                "unexpected '{}' at pos {} in atom spec",
                c as char, self.pos
            )),
        }
    }

    fn parse_number(&mut self) -> Result<i32, String> {
        let start = self.pos;
        while self.pos < self.input.len() && self.input[self.pos].is_ascii_digit() {
            self.pos += 1;
        }
        if self.pos == start {
            return Err(format!("expected number at pos {}", self.pos));
        }
        let s = std::str::from_utf8(&self.input[start..self.pos]).map_err(|_| "invalid utf8")?;
        s.parse::<i32>().map_err(|e| e.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_pattern() {
        let p = parse_smarts("[O:1]=[C:2]!@;-[O:3]~[CH0:4]").unwrap();
        assert_eq!(p.atoms.len(), 4);
        assert_eq!(p.bonds.len(), 3);
        assert_eq!(p.atoms[0].map_idx, Some(1));
        assert_eq!(p.atoms[3].map_idx, Some(4));
    }

    #[test]
    fn test_two_letter_element_symbols() {
        let p = parse_smarts("[Mg:1][Br:2][As:3]").unwrap();
        assert_eq!(p.atoms.len(), 3);
        assert!(matches!(p.atoms[0].query, AtomQuery::Element(12)));
        assert!(matches!(p.atoms[1].query, AtomQuery::Element(35)));
        assert!(matches!(p.atoms[2].query, AtomQuery::Element(33)));

        // Aromatic two-letter elements.
        let p = parse_smarts("[se:1][as:2]").unwrap();
        assert!(matches!(p.atoms[0].query, AtomQuery::AromaticElem(34)));
        assert!(matches!(p.atoms[1].query, AtomQuery::AromaticElem(33)));

        // Unknown symbols still fail.
        assert!(parse_smarts("[Xx:1]").is_err());
    }

    #[test]
    fn test_recursive_smarts() {
        let p = parse_smarts("[$([CX3]=O):1][NX3H1:2]!@;-[c:3][cH1:4]").unwrap();
        assert_eq!(p.atoms.len(), 4);
        assert!(
            matches!(p.atoms[0].query, AtomQuery::Recursive(_)),
            "expected a recursive query, got {:?}",
            p.atoms[0].query
        );
        if let AtomQuery::Recursive(ref inner) = p.atoms[0].query {
            assert_eq!(inner.atoms.len(), 2); // CX3 and O
        }
    }

    #[test]
    fn test_branch() {
        let p = parse_smarts("[a:1][c:2]([a])!@;-[O:3][C:4]").unwrap();
        assert_eq!(p.atoms.len(), 5); // a, c, a_branch, O, C
        assert_eq!(p.bonds.len(), 4);
    }

    #[test]
    fn test_chiral_atom_query() {
        let p = parse_smarts("[C@@H:1]").unwrap();
        assert_eq!(p.atoms.len(), 1);
        assert_eq!(p.atoms[0].map_idx, Some(1));
        assert!(
            matches!(p.atoms[0].query, AtomQuery::And(_)),
            "expected an AND query for a chiral atom, got {:?}",
            p.atoms[0].query
        );
        if let AtomQuery::And(ref parts) = p.atoms[0].query {
            assert!(parts
                .iter()
                .any(|q| matches!(q, AtomQuery::Chiral(ChiralType::TetrahedralCW))));
            assert!(parts.iter().any(|q| matches!(q, AtomQuery::TotalH(1))));
        }
    }

    #[test]
    fn test_ring_size_range() {
        let p = parse_smarts("[c;r{9-}:2]").unwrap();
        assert_eq!(p.atoms.len(), 1);
        if let AtomQuery::And(ref parts) = p.atoms[0].query {
            assert!(parts.iter().any(|q| matches!(q, AtomQuery::RingSizeMin(9))));
        }
    }

    #[test]
    fn test_parse_all_csd_patterns() {
        let data = include_str!("../../tests/fixtures/smarts_patterns.txt");
        let mut ok = 0;
        let mut fail = 0;
        let mut failures = Vec::new();
        for line in data.lines() {
            let smarts = line.split('\t').next().unwrap().trim();
            if smarts.is_empty() {
                continue;
            }
            match parse_smarts(smarts) {
                Ok(p) => {
                    ok += 1;
                    let mapped: Vec<_> = p.atoms.iter().filter(|a| a.map_idx.is_some()).collect();
                    if mapped.len() != 4 {
                        failures.push(format!("WARN mapped={}: {}", mapped.len(), smarts));
                    }
                }
                Err(e) => {
                    fail += 1;
                    failures.push(format!("FAIL: {} → {}", smarts, e));
                }
            }
        }
        for f in &failures {
            eprintln!("{}", f);
        }
        eprintln!("\nParsed: {} ok, {} failed out of {}", ok, fail, ok + fail);
        assert_eq!(fail, 0, "{} patterns failed to parse", fail);
    }
}
