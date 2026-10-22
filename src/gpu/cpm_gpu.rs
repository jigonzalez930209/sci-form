//! GPU-assisted CPM charges: GPU Coulomb matrix + CPU SCF update.
//!
//! Builds the dielectrically screened Coulomb coupling matrix **J** on the
//! GPU and then equilibrates the constant-potential charges on the CPU.
//!
//! # Contract with callers
//!
//! The CPU implementation in [`crate::cpm::compute_cpm_charges`] always
//! works; this module is the optional accelerator. The GPU path only replaces
//! the O(N²) construction of **J**.
//!
//! # Units
//!
//! `positions` are in **Å**; the coupling entries are in **eV·Å** and the
//! charges, grand potential and electrostatic energy in **e** and **eV**.
//!
//! # Numerical differences w.r.t. the CPU implementation
//!
//! Two deliberate differences exist:
//!
//! 1. The charge equilibration uses a **damped fixed-point iteration**
//!    (`q ← (q + q_new)/2`) instead of the exact LU solve of the CPU path.
//!    For electrode-sized conducting systems the fixed-point iteration can
//!    fail to converge, so `converged` must be checked by the caller.
//! 2. `config.mu_per_atom` (multi-electrode setups) is **ignored**; the
//!    global `config.mu_ev` is applied to every atom.

use super::context::{
    bytes_to_f32_vec, ceil_div_u32, pack_uniform_values, pack_vec3_positions_f32,
    ComputeBindingDescriptor, ComputeBindingKind, ComputeDispatchDescriptor, GpuContext,
    UniformValue,
};
use crate::cpm::grand_potential::{chi_element, eta_element};
use crate::cpm::{CpmConfig, CpmResult};

const EV_PER_ANGSTROM: f64 = 14.3996;

/// Computes constant-potential equilibrium charges using a GPU Coulomb matrix.
///
/// `elements` holds the atomic numbers and `positions` the Cartesian
/// coordinates in **Å**. The returned [`CpmResult`] holds the partial
/// charges (**e**), the total charge (**e**), the grand potential Ω and the
/// electrostatic energy (both **eV**), the applied potential (**eV**), the
/// number of iterations and the convergence flag.
///
/// # Errors
///
/// - `"elements/position length mismatch"` when the two slices differ in
///   length.
/// - `"System too small for GPU dispatch"` when fewer than two atoms are
///   supplied.
/// - `"No output from CPM Coulomb kernel"` when the dispatch returned no
///   read-back payload.
/// - `"Output size mismatch: expected {n}, got {m}"` when the kernel produced
///   the wrong number of matrix elements.
/// - Any error from [`GpuContext::run_compute`].
pub fn compute_cpm_charges_gpu(
    ctx: &GpuContext,
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &CpmConfig,
) -> Result<CpmResult, String> {
    let n = elements.len();
    if positions.len() != n {
        return Err("elements/position length mismatch".to_string());
    }
    if n < 2 {
        return Err("System too small for GPU dispatch".to_string());
    }

    let params_bytes = pack_uniform_values(&[
        UniformValue::U32(n as u32),
        UniformValue::U32(0),
        UniformValue::U32(0),
        UniformValue::U32(0),
        UniformValue::F32(config.dielectric as f32),
        UniformValue::F32(EV_PER_ANGSTROM as f32),
        UniformValue::F32(0.0),
        UniformValue::F32(0.0),
    ]);

    let descriptor = ComputeDispatchDescriptor {
        label: "cpm coulomb".to_string(),
        shader_source: CPM_COULOMB_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [ceil_div_u32(n, 16), ceil_div_u32(n, 16), 1],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "positions".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: pack_vec3_positions_f32(positions),
            },
            ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: params_bytes,
            },
            ComputeBindingDescriptor {
                label: "output".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: vec![0u8; n * n * 4],
            },
        ],
    };

    let mut outputs = ctx.run_compute(&descriptor)?.outputs;
    let bytes = outputs.pop().ok_or("No output from CPM Coulomb kernel")?;
    let j_matrix = bytes_to_f32_vec(&bytes);
    if j_matrix.len() != n * n {
        return Err(format!(
            "Output size mismatch: expected {}, got {}",
            n * n,
            j_matrix.len()
        ));
    }

    let mut charges = vec![0.0; n];
    let mut converged = false;
    let mut iterations = 0;

    for iteration in 0..config.max_iter {
        iterations = iteration + 1;
        let mut max_change = 0.0f64;

        for i in 0..n {
            let mut coupling = 0.0;
            for j in 0..n {
                if i != j {
                    coupling += j_matrix[i * n + j] as f64 * charges[j];
                }
            }

            let new_q =
                (config.mu_ev - chi_element(elements[i]) - coupling) / eta_element(elements[i]);
            let change = (new_q - charges[i]).abs();
            max_change = max_change.max(change);
            charges[i] = 0.5 * charges[i] + 0.5 * new_q;
        }

        if max_change < config.charge_tol {
            converged = true;
            break;
        }
    }

    let total_charge: f64 = charges.iter().sum();
    let mut electrostatic_energy = 0.0;
    for i in 0..n {
        electrostatic_energy += chi_element(elements[i]) * charges[i];
        electrostatic_energy += 0.5 * eta_element(elements[i]) * charges[i] * charges[i];
    }
    for i in 0..n {
        for j in (i + 1)..n {
            electrostatic_energy += charges[i] * j_matrix[i * n + j] as f64 * charges[j];
        }
    }

    Ok(CpmResult {
        grand_potential: electrostatic_energy - config.mu_ev * total_charge,
        electrostatic_energy,
        total_charge,
        mu_ev: config.mu_ev,
        charges,
        iterations,
        converged,
    })
}

/// WGSL source of the CPM Coulomb-matrix kernel.
///
/// One invocation per matrix element (`@workgroup_size(16, 16, 1)`); the
/// kernel writes the dielectrically screened coupling
/// `J_ij = 14.3996 / (ε · r_ij)` (in **eV·Å**, positions in **Å**) for
/// `i ≠ j` and zero on the diagonal. The uniform block carries `n_atoms`, the
/// solvent dielectric constant ε (dimensionless), the Coulomb constant
/// 14.3996 eV·Å and two padding words.
pub const CPM_COULOMB_SHADER: &str = r#"
struct AtomPos {
    x: f32, y: f32, z: f32, _pad: f32,
};

struct Params {
    n_atoms: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
    dielectric: f32,
    ev_per_angstrom: f32,
    _pad3: f32,
    _pad4: f32,
};

@group(0) @binding(0) var<storage, read> positions: array<AtomPos>;
@group(0) @binding(1) var<uniform> params: Params;
@group(0) @binding(2) var<storage, read_write> output: array<f32>;

@compute @workgroup_size(16, 16, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let n = params.n_atoms;
    if (i >= n || j >= n) { return; }
    if (i == j) {
        output[i * n + j] = 0.0;
        return;
    }

    let pi = positions[i];
    let pj = positions[j];
    let dx = pi.x - pj.x;
    let dy = pi.y - pj.y;
    let dz = pi.z - pj.z;
    let r = sqrt(dx * dx + dy * dy + dz * dz);

    if (r > 1e-10) {
        output[i * n + j] = params.ev_per_angstrom / (params.dielectric * r);
    } else {
        output[i * n + j] = 0.0;
    }
}
"#;
