//! **ALPHA** — ReaxFF reactive force field.
//!
//! Implements continuous bond-order based reactive force fields that allow
//! bond breaking and formation during molecular dynamics.
//!
//! # Modules
//!
//! - [`bond_order`] — Continuous bond-order calculation from distances
//! - [`taper`] — Hermite tapering polynomial for smooth cutoff
//! - [`energy`] — Bonded energy terms as functions of continuous BO
//! - [`nonbonded`] — van der Waals + shielded Coulomb
//! - [`eem`] — Electronegativity Equalization Method for charges
//! - [`params`] — ReaxFF parameter parsing
//! - [`gradients`] — Analytical gradients for all terms
//!
//! # Feature gate
//!
//! Requires `alpha-reaxff`:
//! ```toml
//! sci-form = { version = "0.11", features = ["alpha-reaxff"] }
//! ```

pub mod bond_order;
pub mod eem;
pub mod energy;
pub mod gradients;
pub mod nonbonded;
pub mod params;
pub mod taper;

pub use bond_order::compute_bond_orders;
pub use energy::ReaxffEnergy;
pub use gradients::compute_reaxff_gradient;
pub use params::ReaxffParams;
