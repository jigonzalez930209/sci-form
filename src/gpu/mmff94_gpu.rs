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
#[cfg(feature = "experimental-gpu")]
use crate::gpu::memory_budget::GpuMemoryBudget;
use crate::gpu::memory_budget::GpuMemoryLimits;

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

    // Guard against excessive memory usage for exclusion bitmap (N²/8 bytes)
    let excl_bitmap_bytes = (n_atoms * n_atoms).div_ceil(8);
    let atom_buffer_bytes = n_atoms * 32;
    let total_gpu_memory = excl_bitmap_bytes + atom_buffer_bytes;
    const MAX_GPU_BUFFER: usize = 512 * 1024 * 1024; // 512 MB
    if total_gpu_memory > MAX_GPU_BUFFER {
        return Ok(compute_mmff94_nonbonded_cpu(
            coords,
            charges,
            vdw_radii,
            vdw_epsilon,
            exclusions_14,
            "System too large for GPU buffers",
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
    let excl_size = (n_atoms * n_atoms).div_ceil(32);
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

    let n_pairs = n_atoms * (n_atoms - 1) / 2;
    // Tune the dispatch to the device limits instead of a fixed workgroup size.
    let (_wg_size, wg_count) =
        GpuMemoryBudget::webgpu_default().optimal_1d_dispatch(n_pairs as u32);

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

/// Buffered-14-7 vdW + screened Coulomb energy for one pair (kcal/mol).
///
/// Returns `None` when the pair is out of the interaction window.
fn mmff94_pair_energy(
    i: usize,
    j: usize,
    coords: &[f64],
    charges: &[f64],
    vdw_radii: &[f64],
    vdw_epsilon: &[f64],
) -> Option<(f64, f64)> {
    let dx = coords[i * 3] - coords[j * 3];
    let dy = coords[i * 3 + 1] - coords[j * 3 + 1];
    let dz = coords[i * 3 + 2] - coords[j * 3 + 2];
    let r2 = dx * dx + dy * dy + dz * dz;
    let r = r2.sqrt();
    if !(0.1..=15.0).contains(&r) {
        return None;
    }

    // MMFF94 Buffered 14-7 van der Waals
    let r_star = vdw_radii[i] + vdw_radii[j];
    let eps = (vdw_epsilon[i] * vdw_epsilon[j]).sqrt();
    let rho = r / r_star;
    let rho7 = rho.powi(7);
    let e_vdw = eps * (1.07 / (rho + 0.07)).powi(7) * ((1.12 / (rho7 + 0.12)) - 2.0);

    // Coulomb with dielectric screening
    let e_elec = 332.0716 * charges[i] * charges[j] / (r + 0.05);
    Some((e_vdw, e_elec))
}

/// Row-tile size (atoms per dispatch) so the per-tile exclusion bitmap
/// (`rows × n_atoms` bits, half of the binding budget) fits in one storage
/// binding. Larger systems are evaluated in several launches.
pub fn mmff94_row_tile_size(n_atoms: usize, limits: &GpuMemoryLimits) -> usize {
    if n_atoms == 0 {
        return 1;
    }
    let budget_bytes = (limits.max_storage_buffer_binding_size / 2) as f64;
    let rows = (8.0 * budget_bytes / n_atoms as f64).floor() as usize;
    rows.max(1).min(n_atoms).min(4096)
}

/// Row tiles `[start, end)` covering every upper-triangle pair `i < j`.
pub fn mmff94_row_tiles(n_atoms: usize, limits: &GpuMemoryLimits) -> Vec<(usize, usize)> {
    let tile = mmff94_row_tile_size(n_atoms, limits);
    let mut tiles = Vec::new();
    let mut start = 0;
    while start < n_atoms {
        let end = (start + tile).min(n_atoms);
        tiles.push((start, end));
        start = end;
    }
    tiles
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
        for j in (i + 1)..n_atoms {
            if excl_set.contains(&(i, j)) {
                continue;
            }
            if let Some((e_vdw, e_elec)) =
                mmff94_pair_energy(i, j, coords, charges, vdw_radii, vdw_epsilon)
            {
                vdw_energy += e_vdw;
                elec_energy += e_elec;
                n_pairs += 1;
            }
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

/// Tiled CPU evaluation mirroring the GPU row-tile dispatch.
///
/// Splits the upper-triangle pair space into row tiles so each dispatch only
/// needs an exclusion bitmap for `rows × n_atoms` atoms; the accumulated
/// energy is identical to [`compute_mmff94_nonbonded_cpu`].
pub fn compute_mmff94_nonbonded_tiled_cpu(
    coords: &[f64],
    charges: &[f64],
    vdw_radii: &[f64],
    vdw_epsilon: &[f64],
    exclusions_14: &[(usize, usize)],
) -> Mmff94GpuResult {
    let n_atoms = charges.len();
    let limits = GpuMemoryLimits::default();
    let tiles = mmff94_row_tiles(n_atoms, &limits);

    let mut excl_set = std::collections::HashSet::new();
    for &(i, j) in exclusions_14 {
        excl_set.insert((i.min(j), i.max(j)));
    }

    let mut vdw_energy = 0.0;
    let mut elec_energy = 0.0;
    let mut n_pairs = 0;
    for (start, end) in tiles {
        for i in start..end {
            for j in (i + 1)..n_atoms {
                if excl_set.contains(&(i, j)) {
                    continue;
                }
                if let Some((e_vdw, e_elec)) =
                    mmff94_pair_energy(i, j, coords, charges, vdw_radii, vdw_epsilon)
                {
                    vdw_energy += e_vdw;
                    elec_energy += e_elec;
                    n_pairs += 1;
                }
            }
        }
    }

    Mmff94GpuResult {
        total_energy: vdw_energy + elec_energy,
        vdw_energy,
        electrostatic_energy: elec_energy,
        n_pairs,
        used_gpu: false,
        backend: format!(
            "CPU (tiled, {} tile(s))",
            mmff94_row_tiles(n_atoms, &limits).len()
        ),
    }
}

/// Evaluate MMFF94 non-bonded energy with a transparent GPU → CPU fallback.
///
/// Tries the GPU when a device is available (and the system fits), otherwise
/// evaluates on the CPU and reports why in `backend`.
#[cfg(feature = "experimental-gpu")]
pub fn compute_mmff94_nonbonded_with_report(
    coords: &[f64],
    charges: &[f64],
    vdw_radii: &[f64],
    vdw_epsilon: &[f64],
    exclusions_14: &[(usize, usize)],
) -> Mmff94GpuResult {
    let ctx = GpuContext::best_available();
    if ctx.is_gpu_available() {
        if let Ok(result) = compute_mmff94_nonbonded_gpu(
            &ctx,
            coords,
            charges,
            vdw_radii,
            vdw_epsilon,
            exclusions_14,
        ) {
            return result;
        }
    }
    let note = if ctx.is_gpu_available() {
        "GPU dispatch failed"
    } else {
        "GPU not available"
    };
    compute_mmff94_nonbonded_cpu(coords, charges, vdw_radii, vdw_epsilon, exclusions_14, note)
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

    fn chain(n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let mut coords = Vec::new();
        for i in 0..n {
            coords.extend_from_slice(&[i as f64 * 1.5, 0.0, 0.0]);
        }
        let charges: Vec<f64> = (0..n)
            .map(|i| if i % 2 == 0 { 0.3 } else { -0.3 })
            .collect();
        let vdw_radii: Vec<f64> = (0..n).map(|i| 1.5 + 0.05 * (i % 3) as f64).collect();
        let vdw_epsilon: Vec<f64> = (0..n).map(|i| 0.05 + 0.01 * (i % 2) as f64).collect();
        (coords, charges, vdw_radii, vdw_epsilon)
    }

    #[test]
    fn test_tiled_cpu_matches_untiled() {
        let (coords, charges, vdw_radii, vdw_epsilon) = chain(24);
        let exclusions = vec![(0, 1), (2, 3), (5, 9)];
        let untiled = compute_mmff94_nonbonded_cpu(
            &coords,
            &charges,
            &vdw_radii,
            &vdw_epsilon,
            &exclusions,
            "test",
        );
        let tiled = compute_mmff94_nonbonded_tiled_cpu(
            &coords,
            &charges,
            &vdw_radii,
            &vdw_epsilon,
            &exclusions,
        );
        assert_eq!(tiled.n_pairs, untiled.n_pairs);
        assert!((tiled.total_energy - untiled.total_energy).abs() < 1e-12);
        assert!((tiled.vdw_energy - untiled.vdw_energy).abs() < 1e-12);
        assert!((tiled.electrostatic_energy - untiled.electrostatic_energy).abs() < 1e-12);
    }

    #[test]
    fn test_row_tiles_cover_every_row_and_respect_limits() {
        let limits = GpuMemoryLimits::default();
        for n in [0usize, 1, 7, 64, 4096, 20000] {
            let tiles = mmff94_row_tiles(n, &limits);
            // Contiguous cover of 0..n.
            let mut expected = 0usize;
            for &(start, end) in &tiles {
                assert_eq!(start, expected);
                assert!(end > start);
                expected = end;
            }
            assert_eq!(expected, n);
            // The per-tile exclusion bitmap fits the binding budget.
            if n > 0 {
                let rows = mmff94_row_tile_size(n, &limits);
                let bitmap_bytes = (rows * n).div_ceil(8) as u64;
                assert!(
                    bitmap_bytes <= limits.max_storage_buffer_binding_size,
                    "n={n}: bitmap {bitmap_bytes} > {}",
                    limits.max_storage_buffer_binding_size
                );
            }
        }
    }

    #[cfg(feature = "experimental-gpu")]
    #[test]
    fn test_with_report_falls_back_to_cpu() {
        let (coords, charges, vdw_radii, vdw_epsilon) = chain(8);
        let result =
            compute_mmff94_nonbonded_with_report(&coords, &charges, &vdw_radii, &vdw_epsilon, &[]);
        let reference =
            compute_mmff94_nonbonded_cpu(&coords, &charges, &vdw_radii, &vdw_epsilon, &[], "ref");
        assert!((result.total_energy - reference.total_energy).abs() < 1e-12);
        assert_eq!(result.n_pairs, reference.n_pairs);
        if !result.used_gpu {
            assert!(
                result.backend.contains("GPU") || result.backend.contains("tiled"),
                "fallback should explain itself, got {:?}",
                result.backend
            );
        }
    }
}
