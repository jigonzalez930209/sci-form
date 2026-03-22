//! GPU-accelerated compute infrastructure for quantum chemistry and visualization.
//!
//! This module provides a production-quality GPU pipeline for:
//!
//! - **Context** — wgpu device initialization with explicit CPU fallback
//! - **Backend report** — Every execution reports which backend was used and why
//! - **Orbital grid** — Evaluate ψ_i(r) on 3D grids (GPU or CPU)
//! - **Marching cubes** — Isosurface extraction from scalar fields
//! - **Isosurface** — High-level orbital mesh generation (dual-lobe, normal smoothing)
//! - **Memory budget** — WebGPU/WASM memory limit enforcement and pre-flight checks
//! - **Shader registry** — Centralized WGSL shader catalogue with tier classification
//! - **Pipeline coordinator** — Multi-kernel dispatch coordination with memory-aware scheduling
//! - **Two-electron GPU** — O(N⁴) ERI computation on GPU (Tier 1)
//! - **Fock build GPU** — Fock matrix F = H + G(P) construction on GPU (Tier 1)
//! - **One-electron GPU** — Overlap/kinetic/nuclear matrix build on GPU (Tier 2)
//! - **Density grid GPU** — Electron density ρ(r) on 3D grids (Tier 1)
//! - **ESP grid GPU** — Point-charge ESP evaluation on 3D grids (Tier 1)
//! - **D4 GPU** — Pairwise two-body D4 dispersion accumulation (Tier 3)
//! - **EEQ GPU** — Damped Coulomb matrix build for electronegativity equilibration (Tier 3)
//! - **CPM GPU** — Coulomb matrix build for constant-potential charges (Tier 3)
//!
//! # Feature Gates
//!
//! - `experimental-gpu-rendering` — Enables `#[repr(C)]` aligned types (bytemuck)
//! - `experimental-gpu` — Full wgpu runtime (Vulkan/Metal/DX12), includes rendering
//!
//! When `experimental-gpu` is disabled, all functions fall back to CPU
//! and the backend report makes the fallback explicit.

pub mod aligned_types;
pub mod alpb_born_gpu;
pub mod ani_aev_gpu;
pub mod backend_report;
pub mod buffer_manager;
pub mod context;
#[cfg(feature = "beta-cpm")]
pub mod cpm_gpu;
pub mod d4_dispersion_gpu;
pub mod density_grid_gpu;
pub mod eeq_gpu;
pub mod esp_grid_gpu;
pub mod fock_build_gpu;
pub mod gamma_matrix_gpu;
pub mod hct_born_gpu;
pub mod isosurface;
pub mod marching_cubes;
pub mod memory_budget;
pub mod mmff94_gpu;
pub mod one_electron_gpu;
pub mod orbital_grid;
pub mod pipeline_coordinator;
pub mod shader_registry;
pub mod two_electron_gpu;
