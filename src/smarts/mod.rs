//! SMARTS pattern parsing and substructure matching.
//!
//! Supports the subset of SMARTS needed for CSD torsion preference matching:
//! atom primitives (element, aromaticity, H count, degree, ring membership, charge, recursive),
//! bond primitives (single, double, any, ring/not-ring), and logical operators (AND, OR, NOT).

mod matcher;
mod parser;
pub mod torsion_data;
pub mod torsion_matcher;

pub use matcher::substruct_match;
pub use matcher::{
    has_substruct_match, precompute_ring_info, substruct_match_batch,
    substruct_match_with_ring_info,
};
#[cfg(feature = "parallel")]
pub use matcher::{has_substruct_match_batch_parallel, substruct_match_batch_parallel};
pub use parser::{parse_smarts, AtomQuery, BondQuery, SmartsAtom, SmartsBond, SmartsPattern};
pub use torsion_matcher::match_experimental_torsions;
