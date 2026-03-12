//! SMARTS pattern parsing and substructure matching.
//!
//! Supports the subset of SMARTS needed for CSD torsion preference matching:
//! atom primitives (element, aromaticity, H count, degree, ring membership, charge, recursive),
//! bond primitives (single, double, any, ring/not-ring), and logical operators (AND, OR, NOT).

mod parser;
mod matcher;
pub mod torsion_data;
pub mod torsion_matcher;

pub use parser::{parse_smarts, SmartsPattern, SmartsAtom, SmartsBond, AtomQuery, BondQuery};
pub use matcher::substruct_match;
pub use torsion_matcher::match_experimental_torsions;
