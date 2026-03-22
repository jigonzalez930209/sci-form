//! GPU-accelerated two-electron repulsion integrals (ERIs).
//!
//! Offloads the O(N⁴) ERI computation to the GPU using a WGSL compute shader.
//! Each workgroup thread computes one unique quartet (i,j,k,l) of contracted ERIs,
//! exploiting 8-fold permutation symmetry to reduce work.
//!
//! For STO-3G basis sets (max 3 primitives/function), this provides significant
//! speedup over CPU for molecules with > ~10 basis functions.

use super::context::{
    ComputeBindingDescriptor, ComputeBindingKind, ComputeDispatchDescriptor, GpuContext,
};
use crate::scf::basis::BasisSet;
use crate::scf::two_electron::TwoElectronIntegrals;

/// Maximum primitives per basis function (STO-3G = 3).
const MAX_PRIMITIVES: usize = 3;

/// Minimum basis functions to justify GPU dispatch overhead.
const GPU_DISPATCH_THRESHOLD: usize = 4;

/// Pack basis set data for the ERI GPU shader.
///
/// Basis buffer layout per function (32 bytes = 8 × f32):
///   center: [f32; 3], lx: u32, ly: u32, lz: u32, n_primitives: u32, _pad: u32
///
/// Primitives buffer layout per function (24 bytes = 6 × f32):
///   3 × (alpha: f32, norm_coeff: f32), padded with zeros if fewer primitives
fn pack_basis_eri(basis: &BasisSet) -> (Vec<u8>, Vec<u8>) {
    let mut basis_bytes = Vec::with_capacity(basis.n_basis * 32);
    let mut prim_bytes = Vec::with_capacity(basis.n_basis * 24);

    for bf in &basis.functions {
        // center xyz (12 bytes)
        basis_bytes.extend_from_slice(&(bf.center[0] as f32).to_ne_bytes());
        basis_bytes.extend_from_slice(&(bf.center[1] as f32).to_ne_bytes());
        basis_bytes.extend_from_slice(&(bf.center[2] as f32).to_ne_bytes());
        // lx, ly, lz (12 bytes)
        basis_bytes.extend_from_slice(&bf.angular[0].to_ne_bytes());
        basis_bytes.extend_from_slice(&bf.angular[1].to_ne_bytes());
        basis_bytes.extend_from_slice(&bf.angular[2].to_ne_bytes());
        // n_primitives + padding (8 bytes)
        basis_bytes.extend_from_slice(&(bf.primitives.len() as u32).to_ne_bytes());
        basis_bytes.extend_from_slice(&0u32.to_ne_bytes());

        // Primitives: 3 × (alpha, norm*coeff)
        for i in 0..MAX_PRIMITIVES {
            if i < bf.primitives.len() {
                let norm = crate::scf::basis::BasisFunction::normalization(
                    bf.primitives[i].alpha,
                    bf.angular[0],
                    bf.angular[1],
                    bf.angular[2],
                );
                prim_bytes.extend_from_slice(&(bf.primitives[i].alpha as f32).to_ne_bytes());
                prim_bytes.extend_from_slice(
                    &((bf.primitives[i].coefficient * norm) as f32).to_ne_bytes(),
                );
            } else {
                prim_bytes.extend_from_slice(&0.0f32.to_ne_bytes());
                prim_bytes.extend_from_slice(&0.0f32.to_ne_bytes());
            }
        }
    }

    (basis_bytes, prim_bytes)
}

/// Enumerate unique quartets with 8-fold symmetry: i≥j, k≥l, ij≥kl.
/// Returns (n_quartets, quartet_indices) where quartet_indices is a flat
/// array of [i, j, k, l] u32 tuples.
fn enumerate_unique_quartets(n: usize) -> (usize, Vec<u8>) {
    let mut quartets = Vec::new();
    let mut count = 0usize;

    for i in 0..n {
        for j in 0..=i {
            let ij = i * n + j;
            for k in 0..n {
                for l in 0..=k {
                    let kl = k * n + l;
                    if ij < kl {
                        continue;
                    }
                    quartets.extend_from_slice(&(i as u32).to_ne_bytes());
                    quartets.extend_from_slice(&(j as u32).to_ne_bytes());
                    quartets.extend_from_slice(&(k as u32).to_ne_bytes());
                    quartets.extend_from_slice(&(l as u32).to_ne_bytes());
                    count += 1;
                }
            }
        }
    }

    (count, quartets)
}

/// Compute two-electron integrals on the GPU.
///
/// Returns `TwoElectronIntegrals` with the full symmetrized tensor,
/// or an error if GPU dispatch fails. Falls back to CPU if the
/// basis is too small to benefit from GPU.
pub fn compute_eris_gpu(
    ctx: &GpuContext,
    basis: &BasisSet,
) -> Result<TwoElectronIntegrals, String> {
    let n = basis.n_basis;

    if n < GPU_DISPATCH_THRESHOLD {
        return Err("Basis too small for GPU dispatch".to_string());
    }

    // Check memory: output is n⁴ × 4 bytes (f32)
    let output_size = (n * n * n * n * 4) as u64;
    if output_size > ctx.capabilities.max_storage_buffer_size {
        return Err(format!(
            "ERI tensor ({} bytes) exceeds GPU storage limit ({} bytes)",
            output_size, ctx.capabilities.max_storage_buffer_size
        ));
    }

    let (basis_bytes, prim_bytes) = pack_basis_eri(basis);
    let (n_quartets, quartet_bytes) = enumerate_unique_quartets(n);

    // Params uniform: [n_basis: u32, n_quartets: u32, pad: u32, pad: u32] = 16 bytes
    let mut params = Vec::with_capacity(16);
    params.extend_from_slice(&(n as u32).to_ne_bytes());
    params.extend_from_slice(&(n_quartets as u32).to_ne_bytes());
    params.extend_from_slice(&0u32.to_ne_bytes());
    params.extend_from_slice(&0u32.to_ne_bytes());

    // Output buffer (f32, n⁴ elements)
    let output_seed = vec![0u8; n * n * n * n * 4];

    let workgroup_count = [(n_quartets as u32).div_ceil(64), 1, 1];

    let descriptor = ComputeDispatchDescriptor {
        label: "two-electron ERI".to_string(),
        shader_source: TWO_ELECTRON_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count,
        bindings: vec![
            ComputeBindingDescriptor {
                label: "basis".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: basis_bytes,
            },
            ComputeBindingDescriptor {
                label: "primitives".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: prim_bytes,
            },
            ComputeBindingDescriptor {
                label: "quartets".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: quartet_bytes,
            },
            ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: params,
            },
            ComputeBindingDescriptor {
                label: "output".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: output_seed,
            },
        ],
    };

    let mut result = ctx.run_compute(&descriptor)?;
    let bytes = result.outputs.pop().ok_or("No output from ERI kernel")?;

    // Convert f32 output to f64
    if bytes.len() != n * n * n * n * 4 {
        return Err(format!(
            "ERI output size mismatch: expected {}, got {}",
            n * n * n * n * 4,
            bytes.len()
        ));
    }

    let data: Vec<f64> = bytes
        .chunks_exact(4)
        .map(|c| f32::from_ne_bytes([c[0], c[1], c[2], c[3]]) as f64)
        .collect();

    Ok(TwoElectronIntegrals::from_raw(data, n))
}

/// WGSL compute shader for two-electron repulsion integrals.
///
/// Each thread processes one unique quartet (i,j,k,l) and writes
/// all 8 symmetry-related elements to the output tensor.
///
/// Computes contracted ERIs: (μν|λσ) = Σ_{pqrs} Nₚcₚ·Nqcq·Nrcr·Nscs · [pq|rs]
/// where [pq|rs] = 2π^{5/2}/(p·q·√(p+q)) · Kab·Kcd · F₀(αpq·|PQ|²)
pub const TWO_ELECTRON_SHADER: &str = r#"
struct BasisFunc {
    cx: f32, cy: f32, cz: f32,
    lx: u32, ly: u32, lz: u32,
    n_prims: u32, _pad: u32,
};

struct Params {
    n_basis: u32,
    n_quartets: u32,
    _pad0: u32,
    _pad1: u32,
};

@group(0) @binding(0) var<storage, read> basis: array<BasisFunc>;
@group(0) @binding(1) var<storage, read> primitives: array<vec2<f32>>;  // (alpha, norm*coeff)
@group(0) @binding(2) var<storage, read> quartets: array<vec4<u32>>;
@group(0) @binding(3) var<uniform> params: Params;
@group(0) @binding(4) var<storage, read_write> output: array<f32>;

// Boys function F₀(x) via series expansion
fn boys_f0(x: f32) -> f32 {
    if (x < 1e-7) {
        return 1.0;
    }
    if (x > 30.0) {
        return 0.8862269 / sqrt(x);  // √π / (2√x)
    }
    var sum: f32 = 1.0;
    var term: f32 = 1.0;
    for (var k: u32 = 1u; k < 50u; k = k + 1u) {
        term *= 2.0 * x / f32(2u * k + 1u);
        sum += term;
        if (abs(term) < 1e-6 * abs(sum)) {
            break;
        }
    }
    return exp(-x) * sum;
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let idx = gid.x;
    if (idx >= params.n_quartets) {
        return;
    }

    let q = quartets[idx];
    let i = q.x;
    let j = q.y;
    let k = q.z;
    let l = q.w;

    let bf_i = basis[i];
    let bf_j = basis[j];
    let bf_k = basis[k];
    let bf_l = basis[l];

    var eri: f32 = 0.0;

    // Contract over primitives (max 3 per function for STO-3G)
    for (var pi: u32 = 0u; pi < bf_i.n_prims; pi = pi + 1u) {
        let prim_a = primitives[i * 3u + pi];
        let alpha = prim_a.x;
        let ca = prim_a.y;

        for (var pj: u32 = 0u; pj < bf_j.n_prims; pj = pj + 1u) {
            let prim_b = primitives[j * 3u + pj];
            let beta = prim_b.x;
            let cb = prim_b.y;

            let p = alpha + beta;
            let mu_ab = alpha * beta / p;
            let ab2 = (bf_i.cx - bf_j.cx) * (bf_i.cx - bf_j.cx)
                     + (bf_i.cy - bf_j.cy) * (bf_i.cy - bf_j.cy)
                     + (bf_i.cz - bf_j.cz) * (bf_i.cz - bf_j.cz);
            let k_ab = exp(-mu_ab * ab2);

            let px = (alpha * bf_i.cx + beta * bf_j.cx) / p;
            let py = (alpha * bf_i.cy + beta * bf_j.cy) / p;
            let pz = (alpha * bf_i.cz + beta * bf_j.cz) / p;

            for (var pk: u32 = 0u; pk < bf_k.n_prims; pk = pk + 1u) {
                let prim_c = primitives[k * 3u + pk];
                let gamma = prim_c.x;
                let cc = prim_c.y;

                for (var pl: u32 = 0u; pl < bf_l.n_prims; pl = pl + 1u) {
                    let prim_d = primitives[l * 3u + pl];
                    let delta = prim_d.x;
                    let cd = prim_d.y;

                    let qq = gamma + delta;
                    let mu_cd = gamma * delta / qq;
                    let cd2 = (bf_k.cx - bf_l.cx) * (bf_k.cx - bf_l.cx)
                             + (bf_k.cy - bf_l.cy) * (bf_k.cy - bf_l.cy)
                             + (bf_k.cz - bf_l.cz) * (bf_k.cz - bf_l.cz);
                    let k_cd = exp(-mu_cd * cd2);

                    let qx = (gamma * bf_k.cx + delta * bf_l.cx) / qq;
                    let qy = (gamma * bf_k.cy + delta * bf_l.cy) / qq;
                    let qz = (gamma * bf_k.cz + delta * bf_l.cz) / qq;

                    let pq2 = (px - qx) * (px - qx)
                            + (py - qy) * (py - qy)
                            + (pz - qz) * (pz - qz);
                    let alpha_pq = p * qq / (p + qq);

                    // prefactor = 2π^{5/2} / (p · q · √(p+q))
                    let prefactor = 2.0 * pow(3.14159265, 2.5) / (p * qq * sqrt(p + qq));

                    eri += ca * cb * cc * cd * prefactor * k_ab * k_cd * boys_f0(alpha_pq * pq2);
                }
            }
        }
    }

    let n = params.n_basis;
    let n2 = n * n;

    // Store with 8-fold symmetry
    output[i * n * n2 + j * n2 + k * n + l] = eri;
    output[j * n * n2 + i * n2 + k * n + l] = eri;
    output[i * n * n2 + j * n2 + l * n + k] = eri;
    output[j * n * n2 + i * n2 + l * n + k] = eri;
    output[k * n * n2 + l * n2 + i * n + j] = eri;
    output[l * n * n2 + k * n2 + i * n + j] = eri;
    output[k * n * n2 + l * n2 + j * n + i] = eri;
    output[l * n * n2 + k * n2 + j * n + i] = eri;
}
"#;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_enumerate_quartets_h2() {
        // H₂ STO-3G: 2 basis functions → unique quartets with i≥j, k≥l, ij≥kl
        let (count, bytes) = enumerate_unique_quartets(2);
        assert!(count > 0);
        assert_eq!(bytes.len(), count * 16); // 4 × u32 per quartet
    }

    #[test]
    fn test_enumerate_quartets_single() {
        // 1 basis function → exactly 1 quartet: (0,0,0,0)
        let (count, _) = enumerate_unique_quartets(1);
        assert_eq!(count, 1);
    }

    #[test]
    fn test_pack_basis_eri() {
        let basis =
            crate::scf::basis::BasisSet::sto3g(&[1, 1], &[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]);
        let (basis_bytes, prim_bytes) = pack_basis_eri(&basis);
        // Basis: 32 bytes per function
        assert_eq!(basis_bytes.len(), basis.n_basis * 32);
        // Primitives: 24 bytes per function (3 × 8)
        assert_eq!(prim_bytes.len(), basis.n_basis * 24);
    }

    #[test]
    fn test_gpu_threshold() {
        let ctx = GpuContext::cpu_fallback();
        let basis = crate::scf::basis::BasisSet::sto3g(&[1], &[[0.0, 0.0, 0.0]]);
        let result = compute_eris_gpu(&ctx, &basis);
        assert!(result.is_err()); // Too small for GPU
    }
}
