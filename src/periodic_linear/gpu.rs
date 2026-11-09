//! GPU-accelerated periodic operator assembly and KPM workloads.
//!
//! Offloads k-point operator assembly and KPM matrix-vector products to GPU
//! compute shaders when the workload is large enough.
//!
//! NOTE: this is currently a **CPU fallback only** — the WGSL kernel it
//! advertises does not exist yet, so the "GPU" path executes exactly the same
//! CPU code and the parity test compares CPU against CPU. The audit
//! (reports/auditoria-total.md, P2-32) flagged this; it is kept as a stable
//! call site for the kernel to be added.
//!
//! # Stability tier
//!
//! **Alpha** — gated behind `alpha-periodic-linear` (plus `experimental-gpu`
//! for the GPU context); the validation gate is
//! `tests/alpha_periodic_validation.rs` (F88).
//!
//! # Units
//!
//! Operators and eigenvalues in **eV**, k-points in **Å⁻¹**.

use nalgebra::DMatrix;

use super::{assemble_periodic_operators, KMesh};
use crate::gpu::context::{
    bytes_to_f32_vec, ceil_div_u32, f32_slice_to_bytes, pack_uniform_values,
    ComputeBindingDescriptor, ComputeBindingKind, ComputeDispatchDescriptor, GpuContext,
    UniformValue,
};

/// Minimum k-mesh size to justify GPU dispatch.
const GPU_PERIODIC_THRESHOLD: usize = 64;

/// GPU-accelerated periodic Bloch Hamiltonian assembly across k-points.
///
/// Each k-point's Bloch sum `H(k) = Σ_R H_R cos(2π k·R)` is independent, so
/// the kernel dispatches one invocation per matrix element per k-point. The
/// CPU path is used below [`GPU_PERIODIC_THRESHOLD`] k-points, when the device
/// is unavailable, or when the blocks are malformed (non-square or ragged), in
/// which case the CPU implementation defines the result.
pub fn assemble_periodic_operators_gpu(
    ctx: &GpuContext,
    h_real_space: &[(&DMatrix<f64>, [i32; 3])],
    s_real_space: &[(&DMatrix<f64>, [i32; 3])],
    kmesh: &KMesh,
) -> (Vec<DMatrix<f64>>, Vec<DMatrix<f64>>) {
    let cpu = || assemble_periodic_operators(h_real_space, s_real_space, kmesh);
    if !ctx.capabilities.gpu_available || kmesh.points.len() < GPU_PERIODIC_THRESHOLD {
        return cpu();
    }
    let Some(n) = matrix_dimension(h_real_space, s_real_space) else {
        return cpu();
    };

    let h = dispatch_bloch_sum(ctx, h_real_space, kmesh, n);
    let s = dispatch_bloch_sum(ctx, s_real_space, kmesh, n);
    match (h, s) {
        (Some(h), Some(s)) => (h, s),
        _ => cpu(),
    }
}

/// Common matrix dimension, or `None` when the blocks are ragged/non-square.
fn matrix_dimension(
    h_real_space: &[(&DMatrix<f64>, [i32; 3])],
    s_real_space: &[(&DMatrix<f64>, [i32; 3])],
) -> Option<usize> {
    let mut n = None;
    for (matrix, _) in h_real_space.iter().chain(s_real_space) {
        if matrix.nrows() != matrix.ncols() {
            return None;
        }
        match n {
            Some(dim) if dim != matrix.nrows() => return None,
            None => n = Some(matrix.nrows()),
            _ => {}
        }
    }
    n.filter(|&n| n > 0)
}

fn dispatch_bloch_sum(
    ctx: &GpuContext,
    blocks: &[(&DMatrix<f64>, [i32; 3])],
    kmesh: &KMesh,
    n: usize,
) -> Option<Vec<DMatrix<f64>>> {
    if blocks.is_empty() {
        return None;
    }
    let n_blocks = blocks.len();
    let n_k = kmesh.points.len();
    let d2 = n * n;
    let total = n_k.checked_mul(d2)?;

    let mut block_bytes = Vec::with_capacity(n_blocks * d2 * 4);
    let mut translations = Vec::with_capacity(n_blocks * 3);
    for (matrix, translation) in blocks {
        // Row-major, matching the WGSL element index and from_row_slice.
        for i in 0..n {
            for j in 0..n {
                block_bytes.extend_from_slice(&(matrix[(i, j)] as f32).to_ne_bytes());
            }
        }
        translations.extend_from_slice(translation);
    }
    let mut translation_bytes = Vec::with_capacity(translations.len() * 4);
    for value in translations {
        translation_bytes.extend_from_slice(&value.to_ne_bytes());
    }
    let kpoints: Vec<f32> = kmesh
        .points
        .iter()
        .flat_map(|point| point.fractional.iter().map(|&x| x as f32))
        .collect();

    let params = pack_uniform_values(&[
        UniformValue::U32(n_blocks as u32),
        UniformValue::U32(n as u32),
        UniformValue::U32(n_k as u32),
        UniformValue::U32(0),
    ]);
    let descriptor = ComputeDispatchDescriptor {
        label: "periodic Bloch sum".to_string(),
        shader_source: PERIODIC_BLOCH_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [ceil_div_u32(total, 64), 1, 1],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "blocks".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: block_bytes,
            },
            ComputeBindingDescriptor {
                label: "translations".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: translation_bytes,
            },
            ComputeBindingDescriptor {
                label: "kpoints".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&kpoints),
            },
            ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: params,
            },
            ComputeBindingDescriptor {
                label: "output".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: f32_slice_to_bytes(&vec![0.0f32; total]),
            },
        ],
    };

    let mut outputs = ctx.run_compute(&descriptor).ok()?.outputs;
    let values = bytes_to_f32_vec(&outputs.pop()?);
    if values.len() != total {
        return None;
    }
    let mut matrices = Vec::with_capacity(n_k);
    for k in 0..n_k {
        let start = k * d2;
        let matrix = DMatrix::from_row_slice(
            n,
            n,
            &values[start..start + d2]
                .iter()
                .map(|&v| v as f64)
                .collect::<Vec<_>>(),
        );
        matrices.push(matrix);
    }
    Some(matrices)
}

/// WGSL source of the Bloch-sum kernel.
///
/// One invocation per `(k-point, matrix element)`. The real part
/// `cos(2π k·R)` is the same convention as
/// [`super::build_bloch_hamiltonian`].
pub const PERIODIC_BLOCH_SHADER: &str = r#"
struct Params {
    n_blocks: u32,
    n: u32,
    n_k: u32,
    _pad: u32,
};

@group(0) @binding(0) var<storage, read> blocks: array<f32>;
@group(0) @binding(1) var<storage, read> translations: array<i32>;
@group(0) @binding(2) var<storage, read> kpoints: array<f32>;
@group(0) @binding(3) var<uniform> params: Params;
@group(0) @binding(4) var<storage, read_write> output: array<f32>;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let n = params.n;
    let d2 = n * n;
    let idx = gid.x;
    if (idx >= params.n_k * d2) { return; }

    let k = idx / d2;
    let element = idx % d2;
    let kx = kpoints[k * 3u];
    let ky = kpoints[k * 3u + 1u];
    let kz = kpoints[k * 3u + 2u];

    var sum: f32 = 0.0;
    for (var b: u32 = 0u; b < params.n_blocks; b = b + 1u) {
        let tx = f32(translations[b * 3u]);
        let ty = f32(translations[b * 3u + 1u]);
        let tz = f32(translations[b * 3u + 2u]);
        let phase = 6.283185307179586 * (kx * tx + ky * ty + kz * tz);
        sum += blocks[b * d2 + element] * cos(phase);
    }
    output[idx] = sum;
}
"#;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::periodic_linear::{KMeshCentering, KMeshConfig};

    #[test]
    fn gpu_assembly_falls_back_to_cpu() {
        let ctx = GpuContext::cpu_fallback();
        let h = DMatrix::<f64>::identity(4, 4);
        let h_pairs: Vec<(&DMatrix<f64>, [i32; 3])> = vec![(&h, [0, 0, 0])];
        let s_pairs: Vec<(&DMatrix<f64>, [i32; 3])> = vec![(&h, [0, 0, 0])];
        let mesh = crate::periodic_linear::monkhorst_pack_mesh(&KMeshConfig {
            grid: [3, 3, 3],
            centering: KMeshCentering::MonkhorstPack,
        })
        .unwrap();

        let (hk_gpu, sk_gpu) = assemble_periodic_operators_gpu(&ctx, &h_pairs, &s_pairs, &mesh);
        let (hk_cpu, sk_cpu) = assemble_periodic_operators(&h_pairs, &s_pairs, &mesh);
        assert_eq!(hk_gpu.len(), hk_cpu.len());
        assert_eq!(sk_gpu.len(), sk_cpu.len());
        for (g, c) in hk_gpu.iter().zip(&hk_cpu) {
            assert!((g - c).norm() < 1e-12, "GPU/CPU H(k) parity violated");
        }
    }
}
