//! GPU-accelerated two-body D4 dispersion accumulation.
//!
//! Evaluates the pairwise (two-body) part of the D4 dispersion model
//! (Grimme et al., *J. Chem. Phys.* 2017, 147, 061701) on the GPU.
//!
//! # Contract with callers
//!
//! The CPU implementation in [`crate::dispersion::compute_d4_energy`] always
//! works; this module is the optional accelerator. Only the **two-body** term
//! is computed on the GPU — the optional three-body ATM term is always taken
//! from the CPU implementation.
//!
//! # Units
//!
//! `positions` are in **Å** (the shader converts them to **Bohr**), `c6` and
//! `c8` are in **Hartree·Bohr⁶** and **Hartree·Bohr⁸**, the returned energies
//! in **Hartree** and `total_kcal_mol` in **kcal/mol**.

use super::context::{
    bytes_to_f32_vec, ceil_div_u32, f32_slice_to_bytes, pack_uniform_values,
    pack_vec3_positions_f32, ComputeBindingDescriptor, ComputeBindingKind,
    ComputeDispatchDescriptor, GpuContext, UniformValue,
};
use crate::dispersion::{c8_from_c6, d4_coordination_number, dynamic_c6, D4Config, D4Result};

const ANG_TO_BOHR: f64 = 1.0 / 0.529177;

/// Computes the two-body D4 dispersion energy with a GPU pair kernel.
///
/// `elements` holds the atomic numbers and `positions` the Cartesian
/// coordinates in **Å**. The returned [`D4Result`] holds the two-body energy
/// (Hartree), the three-body ATM energy (Hartree, always from the CPU path),
/// the total (Hartree and kcal/mol) and the per-atom coordination numbers
/// (dimensionless).
///
/// # Errors
///
/// - `"elements/position length mismatch"` when the two slices differ in
///   length.
/// - `"System too small for GPU dispatch"` when fewer than two atoms are
///   supplied.
/// - `"No output from D4 dispersion kernel"` when the dispatch returned no
///   read-back payload.
/// - `"Output size mismatch: expected {n}, got {m}"` when the kernel produced
///   the wrong number of pair energies.
/// - Any error from [`GpuContext::run_compute`].
pub fn compute_d4_energy_gpu(
    ctx: &GpuContext,
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &D4Config,
) -> Result<D4Result, String> {
    let n = elements.len();
    if positions.len() != n {
        return Err("elements/position length mismatch".to_string());
    }
    if n < 2 {
        return Err("System too small for GPU dispatch".to_string());
    }

    let cn = d4_coordination_number(elements, positions);
    let mut pair_params = vec![0.0f32; n * n * 4];
    for i in 0..n {
        for j in (i + 1)..n {
            let c6 = dynamic_c6(elements[i], elements[j], cn[i], cn[j]);
            let c8 = c8_from_c6(c6, elements[i], elements[j]);
            let r0 = if c6 > 1e-10 { (c8 / c6).sqrt() } else { 5.0 };
            let r_cut = config.a1 * r0 + config.a2;
            let base = (i * n + j) * 4;
            pair_params[base] = c6 as f32;
            pair_params[base + 1] = c8 as f32;
            pair_params[base + 2] = r_cut as f32;
            pair_params[base + 3] = r_cut as f32;
        }
    }

    let params_bytes = pack_uniform_values(&[
        UniformValue::U32(n as u32),
        UniformValue::U32(0),
        UniformValue::U32(0),
        UniformValue::U32(0),
        UniformValue::F32(config.s6 as f32),
        UniformValue::F32(config.s8 as f32),
        UniformValue::F32(ANG_TO_BOHR as f32),
        UniformValue::F32(0.0),
    ]);

    let descriptor = ComputeDispatchDescriptor {
        label: "d4 dispersion".to_string(),
        shader_source: D4_DISPERSION_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [ceil_div_u32(n, 16), ceil_div_u32(n, 16), 1],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "positions".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: pack_vec3_positions_f32(positions),
            },
            ComputeBindingDescriptor {
                label: "pair_params".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&pair_params),
            },
            ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: params_bytes,
            },
            ComputeBindingDescriptor {
                label: "output".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: f32_slice_to_bytes(&vec![0.0f32; n * n]),
            },
        ],
    };

    let mut outputs = ctx.run_compute(&descriptor)?.outputs;
    let bytes = outputs.pop().ok_or("No output from D4 dispersion kernel")?;
    let pair_energies = bytes_to_f32_vec(&bytes);
    if pair_energies.len() != n * n {
        return Err(format!(
            "Output size mismatch: expected {}, got {}",
            n * n,
            pair_energies.len()
        ));
    }

    let mut e2 = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            e2 += pair_energies[i * n + j] as f64;
        }
    }

    let e3 = if config.three_body {
        crate::dispersion::compute_d4_energy(elements, positions, config).e3_body
    } else {
        0.0
    };
    let total = e2 + e3;

    Ok(D4Result {
        e2_body: e2,
        e3_body: e3,
        total_energy: total,
        total_kcal_mol: total * 627.509,
        coordination_numbers: cn,
    })
}

/// WGSL source of the two-body D4 dispersion kernel.
///
/// One invocation per atom pair (`@workgroup_size(16, 16, 1)`); the kernel
/// writes the damped pair energy
/// `−s6·C6/r⁶·f₆ − s8·C8/r⁸·f₈` (in **Hartree**, `r` in **Bohr**) for `j > i`
/// and zero for the lower triangle. The uniform block carries `n_atoms`, the
/// scaling factors `s6` and `s8` (dimensionless), the Å→Bohr conversion
/// factor and one padding word.
pub const D4_DISPERSION_SHADER: &str = r#"
struct AtomPos {
    x: f32, y: f32, z: f32, _pad: f32,
};

struct Params {
    n_atoms: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
    s6: f32,
    s8: f32,
    ang_to_bohr: f32,
    _pad3: f32,
};

@group(0) @binding(0) var<storage, read> positions: array<AtomPos>;
@group(0) @binding(1) var<storage, read> pair_params: array<f32>;
@group(0) @binding(2) var<uniform> params: Params;
@group(0) @binding(3) var<storage, read_write> output: array<f32>;

@compute @workgroup_size(16, 16, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let n = params.n_atoms;
    if (i >= n || j >= n) { return; }
    if (j <= i) {
        output[i * n + j] = 0.0;
        return;
    }

    let base = (i * n + j) * 4u;
    let c6 = pair_params[base];
    if (c6 <= 1e-10) {
        output[i * n + j] = 0.0;
        return;
    }
    let c8 = pair_params[base + 1u];
    let r_cut6 = pair_params[base + 2u];
    let r_cut8 = pair_params[base + 3u];

    let pi = positions[i];
    let pj = positions[j];
    let dx = (pi.x - pj.x) * params.ang_to_bohr;
    let dy = (pi.y - pj.y) * params.ang_to_bohr;
    let dz = (pi.z - pj.z) * params.ang_to_bohr;
    let r = sqrt(dx * dx + dy * dy + dz * dz);
    if (r < 1e-10) {
        output[i * n + j] = 0.0;
        return;
    }

    let r2 = r * r;
    let r6 = r2 * r2 * r2;
    let damp6 = r6 / (r6 + pow(r_cut6, 6.0));
    let r8 = r6 * r2;
    let damp8 = r8 / (r8 + pow(r_cut8, 8.0));

    output[i * n + j] = -params.s6 * c6 / r6 * damp6 - params.s8 * c8 / r8 * damp8;
}
"#;
