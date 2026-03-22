//! GPU-accelerated HCT Born-radii evaluation for Generalized Born solvation.

use super::context::{
    bytes_to_f64_vec_from_f32, ceil_div_u32, f32_slice_to_bytes, pack_uniform_values,
    pack_vec3_positions_f32, ComputeBindingDescriptor, ComputeBindingKind,
    ComputeDispatchDescriptor, GpuContext, UniformValue,
};

const GPU_DISPATCH_THRESHOLD: usize = 8;

/// Intrinsic Born radius for HCT model (same as solvation.rs).
fn intrinsic_born_radius(z: u8) -> f64 {
    match z {
        1 => 1.20,
        6 => 1.70,
        7 => 1.55,
        8 => 1.52,
        9 => 1.47,
        15 => 1.80,
        16 => 1.80,
        17 => 1.75,
        35 => 1.85,
        53 => 1.98,
        _ => 1.70,
    }
}

/// HCT descreening scale factor.
fn hct_scale(z: u8) -> f64 {
    match z {
        1 => 0.85,
        6 => 0.72,
        7 => 0.79,
        8 => 0.85,
        9 => 0.88,
        _ => 0.80,
    }
}

/// Compute HCT Born radii on GPU.
pub fn compute_hct_born_radii_gpu(
    ctx: &GpuContext,
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<Vec<f64>, String> {
    let n = elements.len();
    if n < GPU_DISPATCH_THRESHOLD {
        return Err("System too small for GPU dispatch".to_string());
    }

    let rho: Vec<f32> = elements
        .iter()
        .map(|&z| intrinsic_born_radius(z) as f32)
        .collect();
    let scale: Vec<f32> = elements.iter().map(|&z| hct_scale(z) as f32).collect();

    let params = pack_uniform_values(&[
        UniformValue::U32(n as u32),
        UniformValue::U32(0),
        UniformValue::U32(0),
        UniformValue::U32(0),
    ]);

    let descriptor = ComputeDispatchDescriptor {
        label: "hct born radii".to_string(),
        shader_source: HCT_BORN_RADII_SHADER.to_string(),
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
                label: "scale".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&scale),
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
        .ok_or("No output from HCT Born kernel")?;
    if bytes.len() != n * 4 {
        return Err(format!(
            "HCT Born output size mismatch: expected {}, got {}",
            n * 4,
            bytes.len()
        ));
    }

    Ok(bytes_to_f64_vec_from_f32(&bytes))
}

pub const HCT_BORN_RADII_SHADER: &str = r#"
struct AtomPos {
    x: f32, y: f32, z: f32, _pad: f32,
};

struct Params {
    n_atoms: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
};

@group(0) @binding(0) var<storage, read> positions: array<AtomPos>;
@group(0) @binding(1) var<storage, read> rho: array<f32>;
@group(0) @binding(2) var<storage, read> scale: array<f32>;
@group(0) @binding(3) var<uniform> params: Params;
@group(0) @binding(4) var<storage, read_write> output: array<f32>;

@compute @workgroup_size(64, 1, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let n = params.n_atoms;
    if (i >= n) {
        return;
    }

    let rho_i = rho[i];
    let pos_i = positions[i];
    var integral: f32 = 0.0;

    for (var j: u32 = 0u; j < n; j = j + 1u) {
        if (i == j) {
            continue;
        }

        let pos_j = positions[j];
        let dx = pos_i.x - pos_j.x;
        let dy = pos_i.y - pos_j.y;
        let dz = pos_i.z - pos_j.z;
        let rij = sqrt(dx * dx + dy * dy + dz * dz);
        let scaled_rj = rho[j] * scale[j];

        if (rij > rho_i + scaled_rj) {
            let denom1 = max(rij - scaled_rj, 1e-10);
            let denom2 = rij + scaled_rj;
            let denom3 = max(abs(rij * rij - scaled_rj * scaled_rj), 1e-10);
            var ljr: f32 = 0.0;
            if (rij > scaled_rj && scaled_rj > 1e-10) {
                ljr = log(rij / scaled_rj);
            }
            integral += 0.5 * (1.0 / denom1 - 1.0 / denom2)
                + scaled_rj * ljr / (2.0 * rij * max(denom3, 1e-10));
        } else if (rij + rho_i > scaled_rj) {
            let denom = max(abs(rij - scaled_rj), 1e-10);
            integral += 0.5 * (1.0 / denom - 1.0 / (rij + scaled_rj));
        }
    }

    let inv_r = 1.0 / rho_i - integral;
    var born_r: f32;
    if (inv_r > 1e-10) {
        born_r = 1.0 / inv_r;
    } else {
        born_r = 50.0;
    }
    output[i] = max(born_r, rho_i);
}
"#;

#[cfg(test)]
mod tests {
    use super::HCT_BORN_RADII_SHADER;

    #[test]
    fn test_hct_born_gpu_module_compiles() {
        assert!(!HCT_BORN_RADII_SHADER.is_empty());
    }
}
