//! HF-3c: Minimal Hartree-Fock with composite corrections.
//!
//! Implements restricted Hartree-Fock with a minimal basis set (MINIX/STO-3G)
//! plus three empirical corrections:
//! - D3: Grimme dispersion for van der Waals interactions
//! - gCP: Geometric counterpoise for basis set superposition error
//! - SRB: Short-range basis correction for basis incompleteness
//!
//! Also includes CIS (Configuration Interaction Singles) for UV-Vis excitations.
//!
//! Reference: Sure, R.; Grimme, S. "Corrected small basis set Hartree-Fock method
//! for large systems." J. Comput. Chem. 34 (2013): 1672–1685.

pub mod api;
pub mod basis;
pub mod cis;
pub mod cisd;
pub mod d3;
pub mod fock;
pub mod gcp;
pub mod integrals;
pub mod nuclear;
pub mod overlap_kin;
pub mod scf;
pub mod scf_trait;
pub mod srb;

pub use api::{solve_hf3c, Hf3cResult, HfConfig};
pub use cisd::{compute_cisd, CisdExcitation, CisdResult};
pub use scf_trait::{
    HfScfSolver, Pm3ScfSolver, ScfConvergenceConfig, ScfOutput, ScfSolver, XtbScfSolver,
};
