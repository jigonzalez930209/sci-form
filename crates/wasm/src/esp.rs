//! ESP (electrostatic potential) WASM bindings.

use crate::helpers::*;
use wasm_bindgen::prelude::*;

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
use crate::webgpu::{self, WasmExecutionMode};

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
use sci_form::gpu::context::{
    ceil_div_u32, f32_slice_to_bytes, pack_uniform_values, pack_vec3_positions_f32,
    ComputeBindingDescriptor, ComputeBindingKind, ComputeDispatchDescriptor, UniformValue,
};

#[cfg(all(
    feature = "experimental-gpu",
    target_arch = "wasm32",
    feature = "parallel"
))]
use rayon::prelude::*;

fn esp_grid_geometry(positions: &[[f64; 3]], spacing: f64, padding: f64) -> ([f64; 3], [usize; 3]) {
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

    (origin, dims)
}

fn grid_info_json(origin: [f64; 3], spacing: f64, dims: [usize; 3]) -> String {
    serde_json::json!({
        "origin": origin,
        "spacing": spacing,
        "dims": dims,
    })
    .to_string()
}

/// Compute the full electrostatic potential grid as JSON.
#[wasm_bindgen]
pub fn compute_esp(elements: &str, coords_flat: &str, spacing: f64, padding: f64) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_esp(&elems, &positions, spacing, padding) {
        Ok(grid) => serialize_or_error(&grid),
        Err(e) => json_error(&e),
    }
}

/// Compute ESP grid values as Float64Array (typed-array transfer).
#[wasm_bindgen]
pub fn compute_esp_grid_typed(
    elements: &str,
    coords_flat: &str,
    spacing: f64,
    padding: f64,
) -> js_sys::Float64Array {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(_) => return js_sys::Float64Array::new_with_length(0),
    };
    match sci_form::compute_esp(&elems, &positions, spacing, padding) {
        Ok(grid) => {
            let arr = js_sys::Float64Array::new_with_length(grid.values.len() as u32);
            arr.copy_from(&grid.values);
            arr
        }
        Err(_) => js_sys::Float64Array::new_with_length(0),
    }
}

/// Get ESP grid metadata (origin, spacing, dimensions) as JSON.
#[wasm_bindgen]
pub fn compute_esp_grid_info(
    elements: &str,
    coords_flat: &str,
    spacing: f64,
    padding: f64,
) -> String {
    let (_, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let (origin, dims) = esp_grid_geometry(&positions, spacing, padding);
    grid_info_json(origin, spacing, dims)
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
#[wasm_bindgen]
pub async fn compute_esp_grid_accelerated(
    elements: &str,
    coords_flat: &str,
    spacing: f64,
    padding: f64,
    execution_mode: &str,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };

    let mode = match webgpu::parse_execution_mode(execution_mode) {
        Ok(mode) => mode,
        Err(err) => return json_error(&err),
    };

    match compute_esp_grid_wasm_accelerated(&elems, &positions, spacing, padding, mode).await {
        Ok((grid, backend, used_gpu, note)) => serde_json::json!({
            "origin": grid.origin,
            "spacing": grid.spacing,
            "dims": grid.dims,
            "values": grid.values,
            "backend": backend,
            "used_gpu": used_gpu,
            "mode": execution_mode,
            "note": note,
        })
        .to_string(),
        Err(err) => json_error(&err),
    }
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
#[wasm_bindgen]
pub async fn compute_esp_grid_accelerated_typed(
    elements: &str,
    coords_flat: &str,
    spacing: f64,
    padding: f64,
    execution_mode: &str,
) -> js_sys::Float64Array {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(_) => return js_sys::Float64Array::new_with_length(0),
    };

    let mode = match webgpu::parse_execution_mode(execution_mode) {
        Ok(mode) => mode,
        Err(_) => return js_sys::Float64Array::new_with_length(0),
    };

    match compute_esp_grid_wasm_accelerated(&elems, &positions, spacing, padding, mode).await {
        Ok((grid, _, _, _)) => {
            let arr = js_sys::Float64Array::new_with_length(grid.values.len() as u32);
            arr.copy_from(&grid.values);
            arr
        }
        Err(_) => js_sys::Float64Array::new_with_length(0),
    }
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
async fn compute_esp_grid_wasm_accelerated(
    elements: &[u8],
    positions: &[[f64; 3]],
    spacing: f64,
    padding: f64,
    mode: WasmExecutionMode,
) -> Result<(sci_form::esp::EspGrid, String, bool, String), String> {
    let pop = sci_form::compute_population(elements, positions)?;

    if mode == WasmExecutionMode::Cpu {
        let grid = sci_form::compute_esp(elements, positions, spacing, padding)?;
        return Ok((
            grid,
            "CPU".to_string(),
            false,
            cpu_grid_note("CPU-only ESP execution requested"),
        ));
    }

    let Some(runtime) = webgpu::try_runtime().await else {
        if mode == WasmExecutionMode::Gpu {
            return Err("WebGPU not available and mode=gpu requires GPU".to_string());
        }
        let grid = sci_form::compute_esp(elements, positions, spacing, padding)?;
        return Ok((
            grid,
            "CPU".to_string(),
            false,
            cpu_grid_note("WebGPU unavailable; CPU ESP fallback used"),
        ));
    };

    let (origin, dims) = esp_grid_geometry(positions, spacing, padding);
    if mode == WasmExecutionMode::Hybrid && dims[0] >= 2 {
        let gpu_x = ((dims[0] * 2) / 3).max(1).min(dims[0] - 1);
        let cpu_origin = [origin[0] + gpu_x as f64 * spacing, origin[1], origin[2]];
        let gpu_dims = [gpu_x, dims[1], dims[2]];
        let cpu_dims = [dims[0] - gpu_x, dims[1], dims[2]];

        let mut values = compute_esp_grid_webgpu_slice(
            positions,
            &pop.mulliken_charges,
            origin,
            spacing,
            gpu_dims,
        )
        .await?;
        values.extend(compute_esp_grid_cpu_slice(
            positions,
            &pop.mulliken_charges,
            cpu_origin,
            spacing,
            cpu_dims,
        ));

        return Ok((
            sci_form::esp::EspGrid {
                origin,
                spacing,
                dims,
                values,
            },
            format!("{}+CPU", runtime.backend),
            true,
            cpu_grid_note(
                "Hybrid ESP grid: WebGPU evaluates the leading x-slabs while CPU computes the remainder",
            ),
        ));
    }

    let values =
        compute_esp_grid_webgpu_slice(positions, &pop.mulliken_charges, origin, spacing, dims)
            .await?;
    Ok((
        sci_form::esp::EspGrid {
            origin,
            spacing,
            dims,
            values,
        },
        runtime.backend,
        true,
        "WebGPU ESP grid dispatch".to_string(),
    ))
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
async fn compute_esp_grid_webgpu_slice(
    positions: &[[f64; 3]],
    charges: &[f64],
    origin: [f64; 3],
    spacing: f64,
    dims: [usize; 3],
) -> Result<Vec<f64>, String> {
    let total = dims[0] * dims[1] * dims[2];
    let descriptor = ComputeDispatchDescriptor {
        label: "esp grid wasm".to_string(),
        shader_source: sci_form::gpu::esp_grid_gpu::ESP_GRID_SHADER.to_string(),
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
                    &charges
                        .iter()
                        .map(|value| *value as f32)
                        .collect::<Vec<_>>(),
                ),
            },
            ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: pack_uniform_values(&[
                    UniformValue::F32(origin[0] as f32),
                    UniformValue::F32(origin[1] as f32),
                    UniformValue::F32(origin[2] as f32),
                    UniformValue::F32(spacing as f32),
                    UniformValue::U32(dims[0] as u32),
                    UniformValue::U32(dims[1] as u32),
                    UniformValue::U32(dims[2] as u32),
                    UniformValue::U32(positions.len() as u32),
                ]),
            },
            ComputeBindingDescriptor {
                label: "output".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: f32_slice_to_bytes(&vec![0.0f32; total]),
            },
        ],
    };

    let mut outputs = webgpu::run_compute_async(&descriptor).await?.outputs;
    let bytes = outputs.pop().ok_or("No output from ESP WebGPU kernel")?;
    Ok(webgpu::bytes_to_f64_vec(&bytes))
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
fn compute_esp_grid_cpu_slice(
    positions: &[[f64; 3]],
    charges: &[f64],
    origin: [f64; 3],
    spacing: f64,
    dims: [usize; 3],
) -> Vec<f64> {
    let total = dims[0] * dims[1] * dims[2];

    #[cfg(feature = "parallel")]
    {
        return (0..total)
            .into_par_iter()
            .map(|flat_idx| esp_value_at_index(flat_idx, positions, charges, origin, spacing, dims))
            .collect();
    }

    #[cfg(not(feature = "parallel"))]
    {
        (0..total)
            .map(|flat_idx| esp_value_at_index(flat_idx, positions, charges, origin, spacing, dims))
            .collect()
    }
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
fn esp_value_at_index(
    flat_idx: usize,
    positions: &[[f64; 3]],
    charges: &[f64],
    origin: [f64; 3],
    spacing: f64,
    dims: [usize; 3],
) -> f64 {
    let ny = dims[1];
    let nz = dims[2];
    let ix = flat_idx / (ny * nz);
    let iy = (flat_idx / nz) % ny;
    let iz = flat_idx % nz;

    let rx = origin[0] + ix as f64 * spacing;
    let ry = origin[1] + iy as f64 * spacing;
    let rz = origin[2] + iz as f64 * spacing;

    let mut phi = 0.0;
    for (atom_index, position) in positions.iter().enumerate() {
        let dx = rx - position[0];
        let dy = ry - position[1];
        let dz = rz - position[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        if dist < 0.01 {
            continue;
        }
        phi += charges[atom_index] / (dist * (1.0 / 0.529177));
    }

    phi
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
fn cpu_grid_note(prefix: &str) -> String {
    if cfg!(feature = "parallel") {
        format!(
            "{prefix}; this browser build can use wasm-bindgen-rayon worker threads after initThreadPool()"
        )
    } else {
        prefix.to_string()
    }
}
