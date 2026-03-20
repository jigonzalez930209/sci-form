//! Phase 5: GPU Rendering Module
//!
//! Provides GPU-accelerated molecular orbital and density visualization:
//!
//! - **Orbital evaluator**: Evaluate ψ_i(x,y,z) on dense 3D grids
//! - **Marching cubes**: GPU-accelerated isosurface extraction (WGSL-ready)
//! - **Isosurface**: Mesh generation for visualization (Three.js / WebGL)

pub mod isosurface;
pub mod marching_cubes_gpu;
pub mod orbital_evaluator;
