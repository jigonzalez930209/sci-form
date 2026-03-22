//! GPU-accelerated MMFF94 force field evaluation.
//!
//! Computes MMFF94 non-bonded interactions (van der Waals + electrostatic)
//! on the GPU. The O(N²) pairwise terms are the bottleneck for large
//! molecules and parallelize naturally on GPU hardware.
//!
//! The bonded terms (bond stretching, angle bending, torsion) remain on CPU
//! as they are O(N) and memory-access-heavy.

#[cfg(feature = "experimental-gpu")]
use crate::gpu::context::{
    ComputeBindingDescriptor, ComputeBindingKind, ComputeDispatchDescriptor, GpuContext,
};

use serde::{Deserialize, Serialize};

/// Result of GPU MMFF94 evaluation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Mmff94GpuResult {
    /// Total energy (kcal/mol).
    pub total_energy: f64,
    /// Van der Waals energy (kcal/mol).
    pub vdw_energy: f64,
    /// Electrostatic energy (kcal/mol).
    pub electrostatic_energy: f64,
    /// Number of atom pairs evaluated.
    pub n_pairs: usize,
    /// Whether GPU was actually used.
    pub used_gpu: bool,
    /// Backend description.
    pub backend: String,
}

/// Minimum atoms to justify GPU dispatch for MMFF94.
#[allow(dead_code)]
const GPU_DISPATCH_THRESHOLD: usize = 50;

/// Compute MMFF94 non-bonded energy on GPU.
///
/// Packs atom coordinates, charges, VdW parameters, and 1-4 exclusion
/// data into GPU buffers. A single WGSL kernel evaluates all N(N-1)/2
/// pair interactions in parallel.
#[cfg(feature = "experimental-gpu")]
pub fn compute_mmff94_nonbonded_gpu(
    ctx: &GpuContext,
    coords: &[f64],
    charges: &[f64],
    vdw_radii: &[f64],
    vdw_epsilon: &[f64],
    exclusions_14: &[(usize, usize)],
) -> Result<Mmff94GpuResult, String> {
    let n_atoms = charges.len();

    if n_atoms < GPU_DISPATCH_THRESHOLD {
        return Ok(compute_mmff94_nonbonded_cpu(
            coords,
            charges,
            vdw_radii,
            vdw_epsilon,
            exclusions_14,
            "Atom count below GPU threshold",
        ));
    }

    // Pack atom data: [x, y, z, charge, r_vdw, eps_vdw, pad, pad] per atom (32 bytes)
    let mut atom_bytes = Vec::with_capacity(n_atoms * 32);
    for i in 0..n_atoms {
        atom_bytes.extend_from_slice(&(coords[i * 3] as f32).to_ne_bytes());
        atom_bytes.extend_from_slice(&(coords[i * 3 + 1] as f32).to_ne_bytes());
        atom_bytes.extend_from_slice(&(coords[i * 3 + 2] as f32).to_ne_bytes());
        atom_bytes.extend_from_slice(&(charges[i] as f32).to_ne_bytes());
        atom_bytes.extend_from_slice(&(vdw_radii[i] as f32).to_ne_bytes());
        atom_bytes.extend_from_slice(&(vdw_epsilon[i] as f32).to_ne_bytes());
        atom_bytes.extend_from_slice(&0.0f32.to_ne_bytes());
        atom_bytes.extend_from_slice(&0.0f32.to_ne_bytes());
    }

    // Pack exclusions as a bitmap for fast lookup
    let excl_size = (n_atoms * n_atoms + 31) / 32;
    let mut excl_bits = vec![0u32; excl_size];
    for &(i, j) in exclusions_14 {
        let bit_idx = i * n_atoms + j;
        excl_bits[bit_idx / 32] |= 1 << (bit_idx % 32);
        let bit_idx2 = j * n_atoms + i;
        excl_bits[bit_idx2 / 32] |= 1 << (bit_idx2 % 32);
    }
    let excl_bytes: Vec<u8> = excl_bits.iter().flat_map(|b| b.to_ne_bytes()).collect();

    // Params: n_atoms, pad, pad, pad
    let mut params = Vec::with_capacity(16);
    params.extend_from_slice(&(n_atoms as u32).to_ne_bytes());
    params.extend_from_slice(&0u32.to_ne_bytes());
    params.extend_from_slice(&0u32.to_ne_bytes());
    params.extend_from_slice(&0u32.to_ne_bytes());

    // Output: per-pair energies (vdw + elec), reduced on GPU
    let n_output = 2; // [total_vdw, total_elec]
    let output_bytes = vec![0u8; n_output * 4];

    let wg_size = 64u32;
    let n_pairs = n_atoms * (n_atoms - 1) / 2;
    let wg_count = ((n_pairs as u32) + wg_size - 1) / wg_size;

    let descriptor = ComputeDispatchDescriptor {
        label: "MMFF94 non-bonded".to_string(),
        shader_source: MMFF94_NB_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [wg_count, 1, 1],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "atoms".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: atom_bytes,
            },
            ComputeBindingDescriptor {
                label: "exclusions".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: excl_bytes,
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
    let bytes = result.outputs.pop().ok_or("No output from MMFF94 kernel")?;

    if bytes.len() < 8 {
        return Err("Insufficient output from MMFF94 kernel".to_string());
    }

    let vdw = f32::from_ne_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]) as f64;
    let elec = f32::from_ne_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]) as f64;

    Ok(Mmff94GpuResult {
        total_energy: vdw + elec,
        vdw_energy: vdw,
        electrostatic_energy: elec,
        n_pairs,
        used_gpu: true,
        backend: ctx.capabilities.backend.clone(),
    })
}

/// CPU fallback for MMFF94 non-bonded energy.
pub fn compute_mmff94_nonbonded_cpu(
    coords: &[f64],
    charges: &[f64],
    vdw_radii: &[f64],
    vdw_epsilon: &[f64],
    exclusions_14: &[(usize, usize)],
    note: &str,
) -> Mmff94GpuResult {
    let n_atoms = charges.len();
    let mut vdw_energy = 0.0;
    let mut elec_energy = 0.0;
    let mut n_pairs = 0;

    // Build exclusion set
    let mut excl_set = std::collections::HashSet::new();
    for &(i, j) in exclusions_14 {
        excl_set.insert((i.min(j), i.max(j)));
    }

    for i in 0..n_atoms {
        let xi = coords[i * 3];
        let yi = coords[i * 3 + 1];
        let zi = coords[i * 3 + 2];

        for j in (i + 1)..n_atoms {
            if excl_set.contains(&(i, j)) {
                continue;
            }

            let dx = xi - coords[j * 3];
            let dy = yi - coords[j * 3 + 1];
            let dz = zi - coords[j * 3 + 2];
            let r2 = dx * dx + dy * dy + dz * dz;
            let r = r2.sqrt();

            if !(0.1..=15.0).contains(&r) {
                continue;
            }

            // MMFF94 Buffered 14-7 van der Waals
            let r_star = vdw_radii[i] + vdw_radii[j];
            let eps = (vdw_epsilon[i] * vdw_epsilon[j]).sqrt();
            let rho = r / r_star;
            let rho7 = rho.powi(7);
            let e_vdw = eps * (1.07 / (rho + 0.07)).powi(7) * ((1.12 / (rho7 + 0.12)) - 2.0);
            vdw_energy += e_vdw;

            // Coulomb with dielectric screening
            let e_elec = 332.0716 * charges[i] * charges[j] / (r + 0.05);
            elec_energy += e_elec;

            n_pairs += 1;
        }
    }

    Mmff94GpuResult {
        total_energy: vdw_energy + elec_energy,
        vdw_energy,
        electrostatic_energy: elec_energy,
        n_pairs,
        used_gpu: false,
        backend: format!("CPU ({})", note),
    }
}

/// WGSL compute shader for MMFF94 non-bonded interactions.
///
/// Each workgroup processes a batch of atom pairs, computing buffered-14-7
/// VdW and Coulomb electrostatics. Exclusions are checked via bitmap.
#[cfg(feature = "experimental-gpu")]
pub const MMFF94_NB_SHADER: &str = r#"
struct Atom {
    x: f32, y: f32, z: f32,
    charge: f32,
    r_vdw: f32, eps_vdw: f32,
    _pad0: f32, _pad1: f32,
};

struct Params {
    n_atoms: u32,
    _pad0: u32, _pad1: u32, _pad2: u32,
};

@group(0) @binding(0) var<storage, read> atoms: array<Atom>;
@group(0) @binding(1) var<storage, read> exclusions: array<u32>;
@group(0) @binding(2) var<uniform> params: Params;
@group(0) @binding(3) var<storage, read_write> output: array<atomic<u32>>;

fn pair_to_ij(pair_idx: u32, n: u32) -> vec2<u32> {
    // Convert linear pair index to (i, j) with i < j
    var i: u32 = 0u;
    var remaining: u32 = pair_idx;
    loop {
        let row_size = n - 1u - i;
        if (remaining < row_size) {
            break;
        }
        remaining -= row_size;
        i += 1u;
    }
    let j = i + 1u + remaining;
    return vec2<u32>(i, j);
}

@compute @workgroup_size(64, 1, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let n = params.n_atoms;
    let n_pairs = n * (n - 1u) / 2u;
    let pair_idx = gid.x;

    if (pair_idx >= n_pairs) {
        return;
    }

    let ij = pair_to_ij(pair_idx, n);
    let i = ij.x;
    let j = ij.y;

    // Check exclusion bitmap
    let bit_idx = i * n + j;
    let word = exclusions[bit_idx / 32u];
    if ((word >> (bit_idx % 32u)) & 1u) != 0u {
        return;
    }

    let a = atoms[i];
    let b = atoms[j];

    let dx = a.x - b.x;
    let dy = a.y - b.y;
    let dz = a.z - b.z;
    let r = sqrt(dx * dx + dy * dy + dz * dz);

    if (r < 0.1 || r > 15.0) {
        return;
    }

    // Buffered 14-7 VdW
    let r_star = a.r_vdw + b.r_vdw;
    let eps = sqrt(a.eps_vdw * b.eps_vdw);
    let rho = r / r_star;
    let rho7 = pow(rho, 7.0);
    let e_vdw = eps * pow(1.07 / (rho + 0.07), 7.0) * ((1.12 / (rho7 + 0.12)) - 2.0);

    // Coulomb with distance-dependent dielectric
    let e_elec = 332.0716 * a.charge * b.charge / (r + 0.05);

    // Atomic add to output (using integer representation)
    let vdw_bits = bitcast<u32>(e_vdw);
    let elec_bits = bitcast<u32>(e_elec);
    atomicAdd(&output[0], vdw_bits);
    atomicAdd(&output[1], elec_bits);
}
"#;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mmff94_cpu_fallback() {
        // Simple 2-atom system
        let coords = vec![0.0, 0.0, 0.0, 3.0, 0.0, 0.0];
        let charges = vec![0.3, -0.3];
        let vdw_radii = vec![1.5, 1.7];
        let vdw_epsilon = vec![0.05, 0.06];

        let result =
            compute_mmff94_nonbonded_cpu(&coords, &charges, &vdw_radii, &vdw_epsilon, &[], "test");
        assert_eq!(result.n_pairs, 1);
        assert!(result.total_energy.is_finite());
        assert!(result.electrostatic_energy < 0.0); // opposite charges attract
        assert!(!result.used_gpu);
    }
}
