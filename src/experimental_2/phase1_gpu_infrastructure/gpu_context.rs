//! WebGPU context initialization for headless compute.
//!
//! Provides a GPU compute context that works in three environments:
//! - Desktop (Vulkan/Metal/DX12) via native wgpu
//! - Browser (WebGPU API) via wasm-bindgen
//! - CPU fallback when no GPU is available
//!
//! # Architecture
//!
//! The GPU context follows a three-tier approach:
//! 1. Request adapter (physical GPU)
//! 2. Request device + queue (logical connection)
//! 3. Create compute pipelines from WGSL shaders
//!
//! Since `wgpu` is not yet a dependency, this module provides the
//! interface specification and a CPU-only fallback. Enable the
//! `experimental-gpu` feature to activate real GPU dispatch.

/// Capabilities of the current compute backend.
#[derive(Debug, Clone)]
pub struct ComputeCapabilities {
    /// Backend name (e.g., "Vulkan", "Metal", "WebGPU", "CPU-fallback").
    pub backend: String,
    /// Maximum workgroup size in X dimension.
    pub max_workgroup_size_x: u32,
    /// Maximum workgroup size in Y dimension.
    pub max_workgroup_size_y: u32,
    /// Maximum total invocations per workgroup.
    pub max_workgroup_invocations: u32,
    /// Maximum storage buffer binding size in bytes.
    pub max_storage_buffer_size: u64,
    /// Whether f64 (double precision) is supported on GPU.
    pub supports_f64: bool,
    /// Whether the GPU is actually available.
    pub gpu_available: bool,
}

impl Default for ComputeCapabilities {
    fn default() -> Self {
        Self {
            backend: "CPU-fallback".to_string(),
            max_workgroup_size_x: 256,
            max_workgroup_size_y: 256,
            max_workgroup_invocations: 256,
            max_storage_buffer_size: u64::MAX,
            supports_f64: true,
            gpu_available: false,
        }
    }
}

/// Handle to a GPU compute context.
///
/// In the CPU-fallback mode, this simply tracks capabilities.
/// With `experimental-gpu`, this wraps wgpu::Device + Queue.
#[derive(Debug)]
pub struct GpuContext {
    /// Detected compute capabilities.
    pub capabilities: ComputeCapabilities,
}

impl GpuContext {
    /// Create a CPU-fallback context (no GPU).
    ///
    /// This is always available and enables running all experimental
    /// algorithms in pure CPU mode for validation and testing.
    pub fn cpu_fallback() -> Self {
        Self {
            capabilities: ComputeCapabilities::default(),
        }
    }

    /// Check whether the GPU backend is available.
    pub fn is_gpu_available(&self) -> bool {
        self.capabilities.gpu_available
    }

    /// Calculate optimal workgroup dispatch for a matrix operation.
    ///
    /// Given a matrix of `n × n`, returns (dispatch_x, dispatch_y) workgroup counts.
    pub fn compute_dispatch_size(&self, n: usize, workgroup_size: u32) -> (u32, u32) {
        let wg = workgroup_size as usize;
        let dispatch = (n + wg - 1) / wg;
        (dispatch as u32, dispatch as u32)
    }

    /// Estimate GPU memory needed for an SCF calculation.
    ///
    /// Returns bytes needed for: S, H, F, P matrices (4 × n² × 8 bytes for f64)
    /// plus atom data and basis function descriptors.
    pub fn estimate_memory_bytes(&self, n_basis: usize, n_atoms: usize) -> u64 {
        let matrix_bytes = 4 * n_basis * n_basis * 8; // S, H, F, P
        let atom_bytes = n_atoms * 32; // GpuAtom = 32 bytes
        let basis_bytes = n_basis * 48; // GpuBasisFunction = 48 bytes
        (matrix_bytes + atom_bytes + basis_bytes) as u64
    }
}

/// WGSL shader source for overlap integral computation.
///
/// This shader computes S_μν for all basis function pairs in parallel.
/// Each workgroup thread handles one (μ, ν) pair.
pub const OVERLAP_SHADER_WGSL: &str = r#"
// Overlap integral compute shader
// Dispatch: ceil(n_basis/16) × ceil(n_basis/16) workgroups

struct Params {
    n_basis: u32,
    n_atoms: u32,
    n_occupied: u32,
    _pad: u32,
};

struct BasisFn {
    center_atom: vec4<f32>,    // xyz = center, w = atom_index
    quantum: vec4<f32>,        // x=n, y=l, z=m, w=zeta
    energy_contraction: vec4<f32>, // x=vsip, y=offset, z=n_prims, w=0
};

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> basis: array<BasisFn>;
@group(0) @binding(2) var<storage, read_write> overlap: array<f32>;

// STO s-s overlap approximation
fn sto_overlap(zeta_a: f32, zeta_b: f32, r: f32) -> f32 {
    let p = 0.5 * (zeta_a + zeta_b) * r;
    if (p > 20.0) { return 0.0; }
    return exp(-p) * (1.0 + p + p * p / 3.0);
}

@compute @workgroup_size(16, 16)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let n = params.n_basis;

    if (i >= n || j >= n) { return; }

    if (i == j) {
        overlap[i * n + j] = 1.0;
        return;
    }

    let bi = basis[i];
    let bj = basis[j];

    let dx = bi.center_atom.x - bj.center_atom.x;
    let dy = bi.center_atom.y - bj.center_atom.y;
    let dz = bi.center_atom.z - bj.center_atom.z;
    let r = sqrt(dx*dx + dy*dy + dz*dz);

    let zeta_a = bi.quantum.w;
    let zeta_b = bj.quantum.w;
    let l_a = u32(bi.quantum.y);
    let l_b = u32(bj.quantum.y);

    var s = sto_overlap(zeta_a, zeta_b, r);

    // Angular momentum scaling
    if (l_a == l_b && l_a > 0u) {
        s *= 0.5;
    } else if (l_a != l_b) {
        s *= 0.6;
    }

    overlap[i * n + j] = s;
}
"#;

/// WGSL shader for Wolfsberg-Helmholtz Hamiltonian construction.
pub const HAMILTONIAN_SHADER_WGSL: &str = r#"
// Hamiltonian H_ij = 0.5 * K * S_ij * (H_ii + H_jj)

struct Params {
    n_basis: u32,
    k_wh: f32,
    _pad1: u32,
    _pad2: u32,
};

struct BasisFn {
    center_atom: vec4<f32>,
    quantum: vec4<f32>,
    energy_contraction: vec4<f32>,
};

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> basis: array<BasisFn>;
@group(0) @binding(2) var<storage, read> overlap: array<f32>;
@group(0) @binding(3) var<storage, read_write> hamiltonian: array<f32>;

@compute @workgroup_size(16, 16)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let n = params.n_basis;

    if (i >= n || j >= n) { return; }

    let vsip_i = basis[i].energy_contraction.x;
    let vsip_j = basis[j].energy_contraction.x;

    if (i == j) {
        hamiltonian[i * n + j] = vsip_i;
        return;
    }

    let s_ij = overlap[i * n + j];
    hamiltonian[i * n + j] = 0.5 * params.k_wh * s_ij * (vsip_i + vsip_j);
}
"#;

/// WGSL shader for density matrix construction P = 2 * C_occ * C_occ^T.
pub const DENSITY_SHADER_WGSL: &str = r#"
// Density matrix P_ij = 2 * sum_k C_ik * C_jk (k over occupied orbitals)

struct Params {
    n_basis: u32,
    n_occupied: u32,
    _pad1: u32,
    _pad2: u32,
};

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> coefficients: array<f32>;
@group(0) @binding(2) var<storage, read_write> density: array<f32>;

@compute @workgroup_size(16, 16)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let n = params.n_basis;
    let n_occ = params.n_occupied;

    if (i >= n || j >= n) { return; }

    var p: f32 = 0.0;
    for (var k: u32 = 0u; k < n_occ; k = k + 1u) {
        p += coefficients[i * n + k] * coefficients[j * n + k];
    }

    density[i * n + j] = 2.0 * p;
}
"#;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cpu_fallback_creation() {
        let ctx = GpuContext::cpu_fallback();
        assert!(!ctx.is_gpu_available());
        assert_eq!(ctx.capabilities.backend, "CPU-fallback");
    }

    #[test]
    fn test_dispatch_size() {
        let ctx = GpuContext::cpu_fallback();
        let (dx, dy) = ctx.compute_dispatch_size(100, 16);
        assert_eq!(dx, 7); // ceil(100/16) = 7
        assert_eq!(dy, 7);
    }

    #[test]
    fn test_memory_estimation() {
        let ctx = GpuContext::cpu_fallback();
        let bytes = ctx.estimate_memory_bytes(100, 20);
        // 4 matrices × 100² × 8 + 20 × 32 + 100 × 48
        assert_eq!(bytes, 4 * 100 * 100 * 8 + 20 * 32 + 100 * 48);
    }
}
