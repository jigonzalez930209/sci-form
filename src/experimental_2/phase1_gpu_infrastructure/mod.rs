//! Phase 1: GPU Infrastructure — Memory Alignment and WebGPU Bridge
//!
//! This module provides the foundational data structures and GPU context
//! needed to transfer molecular data to GPU memory with zero-copy semantics.
//!
//! # Components
//!
//! - [`aligned_types`] — `#[repr(C)]` structures padded to GPU alignment rules
//! - [`gpu_context`] — wgpu device/adapter initialization (headless + WASM)
//! - [`buffer_manager`] — Storage buffer allocation, upload, and readback

pub mod aligned_types;
pub mod buffer_manager;
pub mod gpu_context;
