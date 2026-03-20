//! # Experimental GPU-Accelerated Quantum Chemistry Engine
//!
//! This module implements a next-generation SCC-DFTB/xTB quantum engine
//! designed for GPU acceleration via WebGPU (wgpu). All code is fully
//! isolated from the production library — no existing functionality
//! is modified or affected.
//!
//! ## Architecture Overview
//!
//! The implementation follows five sequential phases:
//!
//! 1. **GPU Infrastructure** — Aligned memory types, wgpu context, buffer management
//! 2. **Quantum Engine** — STO basis sets, analytic Gaussian integrals, Hamiltonian assembly
//! 3. **SCF Engine** — Self-consistent field loop, DIIS acceleration, geometry optimization
//! 4. **Spectroscopy** — sTDA UV-Vis, semi-numerical Hessian IR, GIAO NMR
//! 5. **GPU Rendering** — Orbital grid evaluation, GPU marching cubes isosurfaces
//!
//! ## Design Principles
//!
//! - All data structures use `#[repr(C)]` alignment for zero-copy GPU transfer
//! - Matrix operations use `nalgebra::DMatrix<f64>` for CPU reference
//! - GPU kernels are written in WGSL (WebGPU Shading Language)
//! - Every GPU computation has a CPU reference for validation
//! - Feature-gated: GPU code behind `experimental-gpu` feature flag

pub mod constants;
pub mod types;

pub mod phase1_gpu_infrastructure;
pub mod phase2_quantum_engine;
pub mod phase3_scf_engine;
pub mod phase4_spectroscopy;
pub mod phase5_gpu_rendering;
