//! Phase 2: Parallel Quantum Engine O(N²)
//!
//! Implements the core quantum chemistry primitives: basis sets, overlap
//! integrals, kinetic/nuclear integrals, and core Hamiltonian assembly.
//! Each routine has a CPU reference implementation and is structured for
//! future GPU dispatch via WebGPU compute shaders.
//!
//! # Components
//!
//! - [`basis_set`] — STO-3G contracted Gaussian basis with full parametrization
//! - [`gaussian_integrals`] — Primitive Gaussian overlap via Obara-Saika recursion
//! - [`overlap_matrix`] — Full S_μν matrix construction
//! - [`kinetic_matrix`] — Kinetic energy integrals T_μν
//! - [`nuclear_matrix`] — Nuclear attraction integrals V_μν (Boys function)
//! - [`core_hamiltonian`] — H⁰ = T + V assembly
//! - [`two_electron`] — (μν|λσ) electron repulsion integrals
//! - [`validation`] — CPU vs GPU result comparison utilities

pub mod basis_set;
pub mod core_hamiltonian;
pub mod gaussian_integrals;
pub mod kinetic_matrix;
pub mod nuclear_matrix;
pub mod overlap_matrix;
pub mod two_electron;
pub mod validation;
