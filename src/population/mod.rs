//! Population analysis: Mulliken and Löwdin partial charges from EHT results.
//!
//! Given MO coefficients C, overlap matrix S, and electron occupations,
//! computes per-atom partial charges via two standard partitioning schemes.

#[allow(clippy::module_inception)]
pub mod population;
pub use population::{
    mulliken_charges, lowdin_charges, compute_population, PopulationResult,
};
