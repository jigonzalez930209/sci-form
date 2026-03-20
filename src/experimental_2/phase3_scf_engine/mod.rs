//! Phase 3: Self-Consistent Field (SCF) Engine and Geometry Optimization
//!
//! Implements the full Roothaan-Hall SCF procedure with DIIS acceleration,
//! Mulliken population analysis, analytical gradients, and L-BFGS geometry
//! optimization.
//!
//! # SCF Algorithm
//!
//! 1. Build one-electron matrices (S, T, V) → H⁰
//! 2. Compute orthogonalization matrix X = S^{-1/2}
//! 3. Initial guess: diagonalize X†H⁰X
//! 4. Build density matrix P from occupied MOs
//! 5. Build Fock matrix F = H⁰ + G(P)
//! 6. DIIS extrapolation for convergence acceleration
//! 7. Transform and diagonalize F → new MO coefficients
//! 8. Check convergence (energy + density)
//! 9. Repeat 4–8 until converged
//!
//! # Components
//!
//! - [`orthogonalization`] — Löwdin and canonical S^{-1/2}
//! - [`density_matrix`] — P_μν from MO coefficients
//! - [`fock_matrix`] — F = H⁰ + G(P) construction
//! - [`scf_loop`] — Main SCF iteration driver
//! - [`diis`] — Direct Inversion in Iterative Subspace
//! - [`energy`] — Total energy evaluation
//! - [`mulliken`] — Mulliken population analysis
//! - [`gradients`] — Analytical nuclear gradients ∂E/∂R
//! - [`geometry_optimizer`] — L-BFGS/steepest descent optimizer

pub mod density_matrix;
pub mod diis;
pub mod energy;
pub mod fock_matrix;
pub mod geometry_optimizer;
pub mod gradients;
pub mod mulliken;
pub mod orthogonalization;
pub mod scf_loop;
