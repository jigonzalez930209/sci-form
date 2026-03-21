//! GPU-accelerated compute infrastructure for orbital grid and visualization.
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
//!
//! # Feature Gates
//!
//! - `experimental-gpu-rendering` — Enables `#[repr(C)]` aligned types (bytemuck)
//! - `experimental-gpu` — Full wgpu runtime (Vulkan/Metal/DX12), includes rendering
//!
//! When `experimental-gpu` is disabled, all functions fall back to CPU
//! and the backend report makes the fallback explicit.

pub mod context;
pub mod backend_report;
pub mod orbital_grid;
pub mod marching_cubes;
pub mod isosurface;
pub mod memory_budget;
pub mod shader_registry;
pub mod pipeline_coordinator;
pub mod aligned_types;
pub mod buffer_manager;
