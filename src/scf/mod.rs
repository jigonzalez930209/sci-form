//! Core SCF (Self-Consistent Field) engine.
//!
//! This module provides shared, production-quality implementations of the
//! fundamental SCF building blocks used across the library:
//!
//! - **DIIS** — Pulay's Direct Inversion in the Iterative Subspace convergence accelerator
//! - **Orthogonalization** — Löwdin S^{-1/2} symmetric orthogonalization
//! - **Density matrix** — RHF density matrix construction from MO coefficients
//! - **Energy** — Electronic, nuclear repulsion, and total energy evaluation
//! - **Mulliken** — Mulliken and Löwdin population analysis
//! - **Fock matrix** — Modular Fock matrix construction (HF and DFTB)
//! - **Two-electron integrals** — (μν|λσ) ERI with sequential and parallel paths
//! - **Core matrices** — Grouped one-electron matrices (S, T, V, H⁰)
//! - **Basis set** — STO-3G contracted Gaussian basis functions
//! - **Gaussian integrals** — Primitive integral evaluation (Obara-Saika)
//! - **Overlap matrix** — S_μν overlap integrals
//! - **Kinetic matrix** — T_μν kinetic energy integrals
//! - **Nuclear matrix** — V_μν nuclear attraction integrals

pub mod basis;
pub mod constants;
pub mod core_matrices;
pub mod density_matrix;
pub mod diis;
pub mod energy;
pub mod extended_basis;
pub mod fock_matrix;
pub mod gaussian_integrals;
pub mod gradients;
pub mod kinetic_matrix;
pub mod mulliken;
pub mod nuclear_matrix;
pub mod orthogonalization;
pub mod overlap_matrix;
pub mod scf_loop;
pub mod two_electron;
pub mod types;
pub mod validation;
