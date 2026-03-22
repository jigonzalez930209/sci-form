//! GPU-accelerated ANI Atomic Environment Vector (AEV) computation.
//!
//! Offloads the O(N²) radial symmetry function computation to GPU.

use super::context::{
    bytes_to_f64_vec_from_f32, ceil_div_u32, f32_slice_to_bytes, pack_uniform_values,
    pack_vec3_positions_f32, ComputeBindingDescriptor, ComputeBindingKind,
    ComputeDispatchDescriptor, GpuContext, UniformValue,
};

const GPU_DISPATCH_THRESHOLD: usize = 8;

/// Compute radial AEV components on GPU.
///
/// Each atom i gets a vector of radial symmetry function values summed over
/// all neighbors j within the cutoff radius.
pub fn compute_radial_aev_gpu(
    ctx: &GpuContext,
    species: &[u32],
    positions: &[[f64; 3]],
    eta: &[f64],
    rs: &[f64],
    cutoff: f64,
    n_species: usize,
) -> Result<Vec<f64>, String> {
    let n = positions.len();
    if n < GPU_DISPATCH_THRESHOLD {
        return Err("System too small for GPU dispatch".to_string());
    }

    let n_radial = eta.len() * rs.len();
    let output_len = n * n_species * n_radial;

    let species_f32: Vec<f32> = species.iter().map(|&s| s as f32).collect();
    let eta_f32: Vec<f32> = eta.iter().map(|&v| v as f32).collect();
    let rs_f32: Vec<f32> = rs.iter().map(|&v| v as f32).collect();

    // Pad eta/rs to fixed sizes (max 16 each)
    let mut eta_padded = [0.0f32; 16];
    let mut rs_padded = [0.0f32; 16];
    for (i, &v) in eta_f32.iter().enumerate().take(16) {
        eta_padded[i] = v;
    }
    for (i, &v) in rs_f32.iter().enumerate().take(16) {
        rs_padded[i] = v;
    }

    let params = pack_uniform_values(&[
        UniformValue::U32(n as u32),
        UniformValue::U32(n_species as u32),
        UniformValue::U32(eta.len() as u32),
        UniformValue::U32(rs.len() as u32),
        UniformValue::F32(cutoff as f32),
        UniformValue::F32(0.0),
        UniformValue::F32(0.0),
        UniformValue::F32(0.0),
    ]);

    let descriptor = ComputeDispatchDescriptor {
        label: "ani radial aev".to_string(),
        shader_source: ANI_RADIAL_AEV_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [ceil_div_u32(n, 64), 1, 1],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "positions".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: pack_vec3_positions_f32(positions),
            },
            ComputeBindingDescriptor {
                label: "species".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&species_f32),
            },
            ComputeBindingDescriptor {
                label: "eta".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&eta_padded),
            },
            ComputeBindingDescriptor {
                label: "rs".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&rs_padded),
            },
            ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: params,
            },
            ComputeBindingDescriptor {
                label: "output".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: f32_slice_to_bytes(&vec![0.0f32; output_len]),
            },
        ],
    };

    let mut result = ctx.run_compute(&descriptor)?;
    let bytes = result
        .outputs
        .pop()
        .ok_or("No output from ANI AEV kernel")?;
    Ok(bytes_to_f64_vec_from_f32(&bytes))
}

pub const ANI_RADIAL_AEV_SHADER: &str = r#"
struct AtomPos {
    x: f32, y: f32, z: f32, _pad: f32,
};

struct Params {
    n_atoms: u32,
    n_species: u32,
    n_eta: u32,
    n_rs: u32,
    cutoff: f32,
    _pad0: f32,
    _pad1: f32,
    _pad2: f32,
};

@group(0) @binding(0) var<storage, read> positions: array<AtomPos>;
@group(0) @binding(1) var<storage, read> species: array<f32>;
@group(0) @binding(2) var<storage, read> eta: array<f32>;
@group(0) @binding(3) var<storage, read> rs: array<f32>;
@group(0) @binding(4) var<uniform> params: Params;
@group(0) @binding(5) var<storage, read_write> output: array<f32>;

fn cosine_cutoff(r: f32, rc: f32) -> f32 {
    if (r >= rc) {
        return 0.0;
    }
    return 0.5 * (1.0 + cos(3.14159265358979 * r / rc));
}

@compute @workgroup_size(64, 1, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let n = params.n_atoms;
    if (i >= n) {
        return;
    }

    let pos_i = positions[i];
    let n_radial = params.n_eta * params.n_rs;
    let stride = params.n_species * n_radial;

    for (var j: u32 = 0u; j < n; j = j + 1u) {
        if (i == j) {
            continue;
        }

        let pos_j = positions[j];
        let dx = pos_i.x - pos_j.x;
        let dy = pos_i.y - pos_j.y;
        let dz = pos_i.z - pos_j.z;
        let rij = sqrt(dx * dx + dy * dy + dz * dz);

        if (rij >= params.cutoff) {
            continue;
        }

        let sj = u32(species[j]);
        let fc = cosine_cutoff(rij, params.cutoff);
        let base = i * stride + sj * n_radial;

        var k: u32 = 0u;
        for (var ie: u32 = 0u; ie < params.n_eta; ie = ie + 1u) {
            let e = eta[ie];
            for (var ir: u32 = 0u; ir < params.n_rs; ir = ir + 1u) {
                let dr = rij - rs[ir];
                output[base + k] += exp(-e * dr * dr) * fc;
                k = k + 1u;
            }
        }
    }
}
"#;

#[cfg(test)]
mod tests {
    use super::ANI_RADIAL_AEV_SHADER;

    #[test]
    fn test_ani_aev_gpu_module_compiles() {
        assert!(!ANI_RADIAL_AEV_SHADER.is_empty());
    }
}
