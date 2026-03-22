//! GPU point-charge ESP grid evaluation.

use super::backend_report::OrbitalGridReport;
use super::context::{
    bytes_to_f64_vec_from_f32, ceil_div_u32, f32_slice_to_bytes, pack_uniform_values,
    pack_vec3_positions_f32, ComputeBindingDescriptor, ComputeBindingKind,
    ComputeDispatchDescriptor, GpuContext, UniformValue,
};
use crate::esp::{compute_esp_grid, EspGrid};

pub fn compute_esp_grid_with_report(
    elements: &[u8],
    positions: &[[f64; 3]],
    mulliken_charges: &[f64],
    spacing: f64,
    padding: f64,
) -> (EspGrid, OrbitalGridReport) {
    let ctx = GpuContext::best_available();
    if ctx.is_gpu_available() {
        match compute_esp_grid_gpu(&ctx, positions, mulliken_charges, spacing, padding) {
            Ok(grid) => {
                let n_points = grid.dims[0] * grid.dims[1] * grid.dims[2];
                return (
                    grid,
                    OrbitalGridReport {
                        backend: ctx.capabilities.backend.clone(),
                        used_gpu: true,
                        attempted_gpu: true,
                        n_points,
                        note: format!("GPU ESP-grid dispatch on {}", ctx.capabilities.backend),
                    },
                );
            }
            Err(_err) => {}
        }
    }

    let grid = compute_esp_grid(elements, positions, mulliken_charges, spacing, padding);
    let n_points = grid.dims[0] * grid.dims[1] * grid.dims[2];
    (
        grid,
        OrbitalGridReport {
            backend: "CPU".to_string(),
            used_gpu: false,
            attempted_gpu: ctx.is_gpu_available(),
            n_points,
            note: if ctx.is_gpu_available() {
                "GPU available but ESP-grid dispatch failed; CPU fallback used".to_string()
            } else {
                "CPU ESP-grid evaluation (GPU not available)".to_string()
            },
        },
    )
}

pub fn compute_esp_grid_gpu(
    ctx: &GpuContext,
    positions: &[[f64; 3]],
    mulliken_charges: &[f64],
    spacing: f64,
    padding: f64,
) -> Result<EspGrid, String> {
    if positions.len() != mulliken_charges.len() {
        return Err("positions/charges length mismatch".to_string());
    }

    let mut min = [f64::MAX; 3];
    let mut max = [f64::MIN; 3];
    for pos in positions {
        for axis in 0..3 {
            min[axis] = min[axis].min(pos[axis]);
            max[axis] = max[axis].max(pos[axis]);
        }
    }

    let origin = [min[0] - padding, min[1] - padding, min[2] - padding];
    let dims = [
        ((max[0] - min[0] + 2.0 * padding) / spacing).ceil() as usize + 1,
        ((max[1] - min[1] + 2.0 * padding) / spacing).ceil() as usize + 1,
        ((max[2] - min[2] + 2.0 * padding) / spacing).ceil() as usize + 1,
    ];
    let total = dims[0] * dims[1] * dims[2];

    let params_bytes = pack_uniform_values(&[
        UniformValue::F32(origin[0] as f32),
        UniformValue::F32(origin[1] as f32),
        UniformValue::F32(origin[2] as f32),
        UniformValue::F32(spacing as f32),
        UniformValue::U32(dims[0] as u32),
        UniformValue::U32(dims[1] as u32),
        UniformValue::U32(dims[2] as u32),
        UniformValue::U32(positions.len() as u32),
    ]);

    let descriptor = ComputeDispatchDescriptor {
        label: "esp grid".to_string(),
        shader_source: ESP_GRID_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [
            ceil_div_u32(dims[0], 8),
            ceil_div_u32(dims[1], 8),
            ceil_div_u32(dims[2], 4),
        ],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "positions".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: pack_vec3_positions_f32(positions),
            },
            ComputeBindingDescriptor {
                label: "charges".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(
                    &mulliken_charges
                        .iter()
                        .map(|value| *value as f32)
                        .collect::<Vec<_>>(),
                ),
            },
            ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: params_bytes,
            },
            ComputeBindingDescriptor {
                label: "output".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: f32_slice_to_bytes(&vec![0.0f32; total]),
            },
        ],
    };

    let mut outputs = ctx.run_compute(&descriptor)?.outputs;
    let bytes = outputs.pop().ok_or("No output from ESP grid kernel")?;
    let values = bytes_to_f64_vec_from_f32(&bytes);
    if values.len() != total {
        return Err(format!(
            "Output size mismatch: expected {}, got {}",
            total,
            values.len()
        ));
    }

    Ok(EspGrid {
        origin,
        spacing,
        dims,
        values,
    })
}

pub const ESP_GRID_SHADER: &str = r#"
struct AtomPos {
    x: f32, y: f32, z: f32, _pad: f32,
};

struct GridParams {
    origin_x: f32, origin_y: f32, origin_z: f32,
    spacing: f32,
    dims_x: u32, dims_y: u32, dims_z: u32,
    n_atoms: u32,
};

@group(0) @binding(0) var<storage, read> positions: array<AtomPos>;
@group(0) @binding(1) var<storage, read> charges: array<f32>;
@group(0) @binding(2) var<uniform> params: GridParams;
@group(0) @binding(3) var<storage, read_write> output: array<f32>;

@compute @workgroup_size(8, 8, 4)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let ix = gid.x;
    let iy = gid.y;
    let iz = gid.z;

    if (ix >= params.dims_x || iy >= params.dims_y || iz >= params.dims_z) {
        return;
    }

    let rx = params.origin_x + f32(ix) * params.spacing;
    let ry = params.origin_y + f32(iy) * params.spacing;
    let rz = params.origin_z + f32(iz) * params.spacing;
    let flat_idx = ix * params.dims_y * params.dims_z + iy * params.dims_z + iz;

    var phi: f32 = 0.0;
    for (var atom: u32 = 0u; atom < params.n_atoms; atom = atom + 1u) {
        let pos = positions[atom];
        let dx = rx - pos.x;
        let dy = ry - pos.y;
        let dz = rz - pos.z;
        let dist = sqrt(dx * dx + dy * dy + dz * dz);
        if (dist < 0.01) { continue; }
        phi += charges[atom] / (dist * 1.88972599);
    }
    output[flat_idx] = phi;
}
"#;
