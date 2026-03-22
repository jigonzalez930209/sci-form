//! GPU-accelerated one-electron matrix construction.
//!
//! Computes overlap (S), kinetic (T), and nuclear attraction (V) matrices
//! on the GPU in a single dispatch. Each thread handles one (μ,ν) pair
//! and outputs three matrix values.
//!
//! For STO-3G basis sets, the O(N²) contraction with max 3×3 primitive
//! pairs per element provides ~9x work per matrix element.

use super::context::{
    ComputeBindingDescriptor, ComputeBindingKind, ComputeDispatchDescriptor, GpuContext,
};
use crate::scf::basis::BasisSet;
use crate::scf::kinetic_matrix::build_kinetic_matrix;
use crate::scf::nuclear_matrix::build_nuclear_matrix;
use crate::scf::overlap_matrix::build_overlap_matrix;

/// Maximum primitives per basis function (STO-3G = 3).
const MAX_PRIMITIVES: usize = 3;

/// Minimum basis functions to justify GPU dispatch.
const GPU_DISPATCH_THRESHOLD: usize = 6;

/// Result of GPU one-electron matrix computation.
pub struct OneElectronResult {
    /// Overlap matrix S, row-major flat (N×N).
    pub overlap: Vec<f64>,
    /// Kinetic energy matrix T, row-major flat (N×N).
    pub kinetic: Vec<f64>,
    /// Nuclear attraction matrix V, row-major flat (N×N).
    pub nuclear: Vec<f64>,
    /// Number of basis functions.
    pub n_basis: usize,
    /// Whether the matrices were actually produced by a GPU dispatch.
    pub used_gpu: bool,
    /// Backend that produced the result.
    pub backend: String,
    /// Human-readable execution note.
    pub note: String,
}

fn flatten_matrix_row_major(matrix: &nalgebra::DMatrix<f64>) -> Vec<f64> {
    (0..matrix.nrows())
        .flat_map(|i| (0..matrix.ncols()).map(move |j| matrix[(i, j)]))
        .collect()
}

fn compute_one_electron_cpu_exact(
    basis: &BasisSet,
    elements: &[u8],
    positions_bohr: &[[f64; 3]],
    note: impl Into<String>,
) -> OneElectronResult {
    let overlap = build_overlap_matrix(basis);
    let kinetic = build_kinetic_matrix(basis);
    let nuclear = build_nuclear_matrix(basis, elements, positions_bohr);

    OneElectronResult {
        overlap: flatten_matrix_row_major(&overlap),
        kinetic: flatten_matrix_row_major(&kinetic),
        nuclear: flatten_matrix_row_major(&nuclear),
        n_basis: basis.n_basis,
        used_gpu: false,
        backend: "CPU-exact".to_string(),
        note: note.into(),
    }
}

fn basis_supports_exact_gpu_kernel(basis: &BasisSet) -> bool {
    basis.functions.iter().all(|bf| bf.l_total == 0)
}

/// Pack basis set and nuclear data for the one-electron GPU shader.
///
/// Basis buffer: per function (32 bytes):
///   center [f32;3], lx: u32, ly: u32, lz: u32, n_prims: u32, _pad: u32
///
/// Primitives buffer: per function (24 bytes):
///   3 × (alpha: f32, norm*coeff: f32)
///
/// Atoms buffer: per atom (16 bytes):
///   position [f32;3], atomic_number: f32
fn pack_one_electron_data(
    basis: &BasisSet,
    elements: &[u8],
    positions_bohr: &[[f64; 3]],
) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    let mut basis_bytes = Vec::with_capacity(basis.n_basis * 32);
    let mut prim_bytes = Vec::with_capacity(basis.n_basis * MAX_PRIMITIVES * 8);

    for bf in &basis.functions {
        basis_bytes.extend_from_slice(&(bf.center[0] as f32).to_ne_bytes());
        basis_bytes.extend_from_slice(&(bf.center[1] as f32).to_ne_bytes());
        basis_bytes.extend_from_slice(&(bf.center[2] as f32).to_ne_bytes());
        basis_bytes.extend_from_slice(&bf.angular[0].to_ne_bytes());
        basis_bytes.extend_from_slice(&bf.angular[1].to_ne_bytes());
        basis_bytes.extend_from_slice(&bf.angular[2].to_ne_bytes());
        basis_bytes.extend_from_slice(&(bf.primitives.len() as u32).to_ne_bytes());
        basis_bytes.extend_from_slice(&0u32.to_ne_bytes());

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

    // Pack atom positions and charges
    let mut atom_bytes = Vec::with_capacity(elements.len() * 16);
    for (i, &z) in elements.iter().enumerate() {
        atom_bytes.extend_from_slice(&(positions_bohr[i][0] as f32).to_ne_bytes());
        atom_bytes.extend_from_slice(&(positions_bohr[i][1] as f32).to_ne_bytes());
        atom_bytes.extend_from_slice(&(positions_bohr[i][2] as f32).to_ne_bytes());
        atom_bytes.extend_from_slice(&(z as f32).to_ne_bytes());
    }

    (basis_bytes, prim_bytes, atom_bytes)
}

/// Compute one-electron matrices (S, T, V) on the GPU.
///
/// Returns overlap, kinetic, and nuclear matrices packed as flat f64 arrays.
pub fn compute_one_electron_gpu(
    ctx: &GpuContext,
    basis: &BasisSet,
    elements: &[u8],
    positions_bohr: &[[f64; 3]],
) -> Result<OneElectronResult, String> {
    let n = basis.n_basis;
    let n_atoms = elements.len();

    if n < GPU_DISPATCH_THRESHOLD {
        return Err("Basis too small for GPU dispatch".to_string());
    }

    if !basis_supports_exact_gpu_kernel(basis) {
        return Ok(compute_one_electron_cpu_exact(
            basis,
            elements,
            positions_bohr,
            "Fell back to exact CPU one-electron builders because the current WGSL kernel only supports pure s-type basis functions.",
        ));
    }

    let (basis_bytes, prim_bytes, atom_bytes) =
        pack_one_electron_data(basis, elements, positions_bohr);

    // Params: n_basis, n_atoms, pad, pad (16 bytes)
    let mut params = Vec::with_capacity(16);
    params.extend_from_slice(&(n as u32).to_ne_bytes());
    params.extend_from_slice(&(n_atoms as u32).to_ne_bytes());
    params.extend_from_slice(&0u32.to_ne_bytes());
    params.extend_from_slice(&0u32.to_ne_bytes());

    // Output: 3 matrices (S, T, V) packed sequentially, each N×N f32
    let output_size = 3 * n * n;
    let output_seed = vec![0.0f32; output_size];
    let output_bytes: Vec<u8> = output_seed.iter().flat_map(|v| v.to_ne_bytes()).collect();

    let wg_size = 16u32;
    let wg_x = (n as u32).div_ceil(wg_size);
    let wg_y = wg_x;

    let descriptor = ComputeDispatchDescriptor {
        label: "one-electron matrices".to_string(),
        shader_source: ONE_ELECTRON_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [wg_x, wg_y, 1],
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
                label: "atoms".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: atom_bytes,
            },
            ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: params,
            },
            ComputeBindingDescriptor {
                label: "output".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: output_bytes,
            },
        ],
    };

    let mut result = ctx.run_compute(&descriptor)?;
    let bytes = result
        .outputs
        .pop()
        .ok_or("No output from one-electron kernel")?;

    if bytes.len() != output_size * 4 {
        return Err(format!(
            "Output size mismatch: expected {}, got {}",
            output_size * 4,
            bytes.len()
        ));
    }

    let all_f64: Vec<f64> = bytes
        .chunks_exact(4)
        .map(|c| f32::from_ne_bytes([c[0], c[1], c[2], c[3]]) as f64)
        .collect();

    let n2 = n * n;
    Ok(OneElectronResult {
        overlap: all_f64[..n2].to_vec(),
        kinetic: all_f64[n2..2 * n2].to_vec(),
        nuclear: all_f64[2 * n2..].to_vec(),
        n_basis: n,
        used_gpu: true,
        backend: ctx.capabilities.backend.clone(),
        note: "Executed WGSL one-electron kernel for a pure s-type basis.".to_string(),
    })
}

/// WGSL compute shader for one-electron matrices (S, T, V).
///
/// Each thread computes S(μ,ν), T(μ,ν), and V(μ,ν) for one basis pair.
/// Exploits symmetric property to only compute upper triangle + reflect.
///
/// For s-type functions:
///   S = (π/p)^{3/2} · exp(-μ·|AB|²)
///   T = β(2l+3)S - 2β² S(a,b+2) [simplified for s-type]
///   V = -2π/p · exp(-μ·|AB|²) · Σ_C Z_C · F₀(p·|PC|²)
pub const ONE_ELECTRON_SHADER: &str = r#"
struct BasisFunc {
    cx: f32, cy: f32, cz: f32,
    lx: u32, ly: u32, lz: u32,
    n_prims: u32, _pad: u32,
};

struct Atom {
    x: f32, y: f32, z: f32,
    charge: f32,
};

struct Params {
    n_basis: u32,
    n_atoms: u32,
    _pad0: u32,
    _pad1: u32,
};

@group(0) @binding(0) var<storage, read> basis: array<BasisFunc>;
@group(0) @binding(1) var<storage, read> primitives: array<vec2<f32>>;
@group(0) @binding(2) var<storage, read> atoms: array<Atom>;
@group(0) @binding(3) var<uniform> params: Params;
@group(0) @binding(4) var<storage, read_write> output: array<f32>;
// Output layout: [S_00..S_{nn}] [T_00..T_{nn}] [V_00..V_{nn}]

fn boys_f0(x: f32) -> f32 {
    if (x < 1e-7) {
        return 1.0;
    }
    if (x > 30.0) {
        return 0.8862269 / sqrt(x);
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

@compute @workgroup_size(16, 16, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let mu = gid.x;
    let nu = gid.y;
    let n = params.n_basis;
    let n2 = n * n;

    if (mu >= n || nu >= n) {
        return;
    }

    // Only compute upper triangle (mu <= nu), then mirror
    let compute_mu = min(mu, nu);
    let compute_nu = max(mu, nu);

    let bf_a = basis[compute_mu];
    let bf_b = basis[compute_nu];

    var s_total: f32 = 0.0;
    var t_total: f32 = 0.0;
    var v_total: f32 = 0.0;

    // Contract over primitives
    for (var pa: u32 = 0u; pa < bf_a.n_prims; pa = pa + 1u) {
        let prim_a = primitives[compute_mu * 3u + pa];
        let alpha = prim_a.x;
        let ca = prim_a.y;

        for (var pb: u32 = 0u; pb < bf_b.n_prims; pb = pb + 1u) {
            let prim_b = primitives[compute_nu * 3u + pb];
            let beta = prim_b.x;
            let cb = prim_b.y;

            let p = alpha + beta;
            let mu_ab = alpha * beta / p;
            let ab2 = (bf_a.cx - bf_b.cx) * (bf_a.cx - bf_b.cx)
                     + (bf_a.cy - bf_b.cy) * (bf_a.cy - bf_b.cy)
                     + (bf_a.cz - bf_b.cz) * (bf_a.cz - bf_b.cz);
            let k_ab = exp(-mu_ab * ab2);

            // Overlap: S = (π/p)^{3/2} · K_ab
            let pi = 3.14159265;
            let s_prim = pow(pi / p, 1.5) * k_ab;
            s_total += ca * cb * s_prim;

            // Kinetic: T = μ(3 - 2μ·|AB|²) · S / p  [s-type simplification]
            let t_prim = mu_ab * (3.0 - 2.0 * mu_ab * ab2) * s_prim;
            t_total += ca * cb * t_prim;

            // Nuclear attraction: V = -Σ_C Z_C · 2π/p · K_ab · F₀(p·|PC|²)
            let px = (alpha * bf_a.cx + beta * bf_b.cx) / p;
            let py = (alpha * bf_a.cy + beta * bf_b.cy) / p;
            let pz = (alpha * bf_a.cz + beta * bf_b.cz) / p;

            for (var c: u32 = 0u; c < params.n_atoms; c = c + 1u) {
                let atom = atoms[c];
                let pc2 = (px - atom.x) * (px - atom.x)
                        + (py - atom.y) * (py - atom.y)
                        + (pz - atom.z) * (pz - atom.z);
                let v_prim = -atom.charge * 2.0 * pi / p * k_ab * boys_f0(p * pc2);
                v_total += ca * cb * v_prim;
            }
        }
    }

    // Write to output: S at offset 0, T at offset n², V at offset 2n²
    output[mu * n + nu] = s_total;
    output[n2 + mu * n + nu] = t_total;
    output[2u * n2 + mu * n + nu] = v_total;
}
"#;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pack_one_electron() {
        let basis =
            crate::scf::basis::BasisSet::sto3g(&[1, 1], &[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]);
        let (basis_bytes, prim_bytes, atom_bytes) =
            pack_one_electron_data(&basis, &[1, 1], &[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]);
        assert_eq!(basis_bytes.len(), basis.n_basis * 32);
        assert_eq!(prim_bytes.len(), basis.n_basis * MAX_PRIMITIVES * 8);
        assert_eq!(atom_bytes.len(), 2 * 16);
    }

    #[test]
    fn test_gpu_threshold() {
        let ctx = GpuContext::cpu_fallback();
        let basis = crate::scf::basis::BasisSet::sto3g(&[1], &[[0.0, 0.0, 0.0]]);
        let result = compute_one_electron_gpu(&ctx, &basis, &[1], &[[0.0, 0.0, 0.0]]);
        assert!(result.is_err());
    }

    #[test]
    fn test_mixed_angular_momentum_falls_back_to_exact_cpu() {
        let elements = [8u8, 1, 1];
        let positions = [
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469],
        ];
        let basis = crate::scf::basis::BasisSet::sto3g(&elements, &positions);
        let ctx = GpuContext::cpu_fallback();

        let result = compute_one_electron_gpu(&ctx, &basis, &elements, &positions)
            .expect("mixed-angular basis should fall back to CPU");

        let overlap = build_overlap_matrix(&basis);
        let overlap_flat = flatten_matrix_row_major(&overlap);
        assert!(!result.used_gpu);
        assert_eq!(result.backend, "CPU-exact");
        assert!(result.note.contains("pure s-type"));
        assert_eq!(result.overlap.len(), overlap_flat.len());
        for (lhs, rhs) in result.overlap.iter().zip(overlap_flat.iter()) {
            assert!((lhs - rhs).abs() < 1e-12);
        }
    }
}
