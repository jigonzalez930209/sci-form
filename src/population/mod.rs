//! Population analysis: Mulliken, Löwdin, and Natural (NPA/NBO) partial charges.
//!
//! Given MO coefficients C, overlap matrix S, and electron occupations,
//! computes per-atom partial charges via multiple partitioning schemes.

pub mod npa;
#[allow(clippy::module_inception)]
pub mod population;
pub use npa::{
    compute_nbo, compute_npa, NaturalConfig, NboBond, NboLonePair, NboResult, NpaResult,
};
pub use population::{
    compute_bond_orders, compute_population, compute_population_with_orbitals, lowdin_charges,
    mulliken_charges, BondOrderEntry, BondOrderResult, PopulationResult,
};
