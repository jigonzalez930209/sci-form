//! GPU-accelerated ALPB Born-radii evaluation.

use super::context::{
    bytes_to_f64_vec_from_f32, ceil_div_u32, f32_slice_to_bytes, pack_uniform_values,
    pack_vec3_positions_f32, ComputeBindingDescriptor, ComputeBindingKind,
    ComputeDispatchDescriptor, GpuContext, UniformValue,
};
use crate::solvation_alpb::{intrinsic_radius, AlpbBornRadii};

const GPU_DISPATCH_THRESHOLD: usize = 3;

pub fn compute_born_radii_gpu(
    ctx: &GpuContext,
    elements: &[u8],
    positions: &[[f64; 3]],
    probe_radius: f64,
) -> Result<AlpbBornRadii, String> {
    let n = elements.len();
    if n < GPU_DISPATCH_THRESHOLD {
        return Err("System too small for GPU dispatch".to_string());
    }
    if positions.len() != n {
        return Err("elements/position length mismatch".to_string());
    }

    let intrinsic: Vec<f64> = elements.iter().map(|&z| intrinsic_radius(z)).collect();
    let rho: Vec<f32> = intrinsic
        .iter()
        .map(|radius| (*radius + probe_radius * 0.1) as f32)
        .collect();

    let params = pack_uniform_values(&[
        UniformValue::U32(n as u32),
        UniformValue::U32(0),
        UniformValue::U32(0),
        UniformValue::U32(0),
        UniformValue::F32(1.0),
        UniformValue::F32(0.8),
        UniformValue::F32(4.85),
        UniformValue::F32(0.0),
    ]);

    let descriptor = ComputeDispatchDescriptor {
        label: "alpb born radii".to_string(),
        shader_source: ALPB_BORN_RADII_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [ceil_div_u32(n, 64), 1, 1],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "positions".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: pack_vec3_positions_f32(positions),
            },
            ComputeBindingDescriptor {
                label: "rho".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&rho),
            },
            ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: params,
            },
            ComputeBindingDescriptor {
                label: "output".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: f32_slice_to_bytes(&vec![0.0f32; n]),
            },
        ],
    };

    let mut result = ctx.run_compute(&descriptor)?;
    let bytes = result
        .outputs
        .pop()
        .ok_or("No output from ALPB Born kernel")?;
    if bytes.len() != n * 4 {
        return Err(format!(
            "ALPB Born output size mismatch: expected {}, got {}",
            n * 4,
            bytes.len()
        ));
    }

    let radii = bytes_to_f64_vec_from_f32(&bytes);

    Ok(AlpbBornRadii { radii, intrinsic })
}

pub const ALPB_BORN_RADII_SHADER: &str = r#"
struct AtomPos {
    x: f32, y: f32, z: f32, _pad: f32,
};

struct Params {
    n_atoms: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
    alpha: f32,
    beta: f32,
    gamma: f32,
    _pad3: f32,
};

@group(0) @binding(0) var<storage, read> positions: array<AtomPos>;
@group(0) @binding(1) var<storage, read> rho: array<f32>;
@group(0) @binding(2) var<uniform> params: Params;
@group(0) @binding(3) var<storage, read_write> output: array<f32>;

@compute @workgroup_size(64, 1, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let n = params.n_atoms;
    if (i >= n) {
        return;
    }

    let rho_i = rho[i];
    let pos_i = positions[i];
    var psi: f32 = 0.0;

    for (var j: u32 = 0u; j < n; j = j + 1u) {
        if (i == j) {
            continue;
        }

        let pos_j = positions[j];
        let dx = pos_i.x - pos_j.x;
        let dy = pos_i.y - pos_j.y;
        let dz = pos_i.z - pos_j.z;
        let r_ij = sqrt(dx * dx + dy * dy + dz * dz);
        let rho_j = rho[j];

        if (r_ij > rho_j) {
            let l_ij = max(rho_i, r_ij - rho_j);
            let u_ij = r_ij + rho_j;
            if (u_ij > l_ij) {
                psi += 0.5 * (
                    (1.0 / l_ij) - (1.0 / u_ij)
                    + 0.25 * ((1.0 / u_ij) - (1.0 / l_ij)) * (r_ij * r_ij - rho_j * rho_j)
                    + 0.5 * ((1.0 / (u_ij * u_ij)) - (1.0 / (l_ij * l_ij))) * r_ij
                );
            }
        }
    }

    let psi_scaled = psi * rho_i;
    let tanh_val = tanh(
        params.alpha * psi_scaled
        - params.beta * psi_scaled * psi_scaled
        + params.gamma * psi_scaled * psi_scaled * psi_scaled
    );

    let inv_r_eff = 1.0 / rho_i - tanh_val / rho_i;
    output[i] = select(100.0, 1.0 / inv_r_eff, inv_r_eff > 1e-10);
}
"#;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alpb_gpu_threshold() {
        let ctx = GpuContext::cpu_fallback();
        let result =
            compute_born_radii_gpu(&ctx, &[8, 1], &[[0.0, 0.0, 0.0], [0.7, 0.0, 0.0]], 1.4);
        assert!(result.is_err());
    }
}
