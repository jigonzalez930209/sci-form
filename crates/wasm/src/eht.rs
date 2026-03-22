//! EHT (Extended Hückel Theory) WASM bindings.

use crate::helpers::*;
use wasm_bindgen::prelude::*;

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
use crate::webgpu::{self, WasmExecutionMode};

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
use sci_form::gpu::context::{
    ceil_div_u32, f32_slice_to_bytes, pack_uniform_values, ComputeBindingDescriptor,
    ComputeBindingKind, ComputeDispatchDescriptor, UniformValue,
};

#[cfg(all(
    feature = "experimental-gpu",
    target_arch = "wasm32",
    feature = "parallel"
))]
use rayon::prelude::*;

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
const ANG_TO_BOHR: f64 = 1.0 / 0.529177249;

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
#[derive(Clone)]
struct ExpandedEhtBasis {
    basis_bytes: Vec<u8>,
    primitive_bytes: Vec<u8>,
    coefficients: Vec<Vec<f32>>,
}

fn grid_info_json(origin: [f64; 3], spacing: f64, dims: [usize; 3]) -> String {
    serde_json::json!({
        "origin": origin,
        "spacing": spacing,
        "dims": dims,
    })
    .to_string()
}

fn orbital_grid_from_coefficients(
    elements: &str,
    coords_flat: &str,
    coefficients_json: &str,
    mo_index: usize,
    spacing: f64,
) -> Result<sci_form::eht::VolumetricGrid, String> {
    let (elems, positions) = parse_elements_and_positions(elements, coords_flat)?;
    let coefficients: Vec<Vec<f64>> =
        serde_json::from_str(coefficients_json).map_err(|e| format!("bad coefficients: {}", e))?;
    let basis = sci_form::eht::basis::build_basis(&elems, &positions);

    if basis.is_empty() {
        return Err("No basis functions found for orbital evaluation".to_string());
    }
    if mo_index >= coefficients.len() {
        return Err(format!(
            "orbital index {} out of range for {} orbitals",
            mo_index,
            coefficients.len()
        ));
    }
    if coefficients.len() != basis.len() {
        return Err(format!(
            "coefficient row count {} does not match basis size {}",
            coefficients.len(),
            basis.len()
        ));
    }
    if coefficients.iter().any(|row| mo_index >= row.len()) {
        return Err(format!(
            "orbital index {} exceeds coefficient columns",
            mo_index
        ));
    }

    #[cfg(feature = "parallel")]
    {
        Ok(sci_form::eht::evaluate_orbital_on_grid_parallel(
            &basis,
            &coefficients,
            mo_index,
            &positions,
            spacing,
            3.0,
        ))
    }

    #[cfg(not(feature = "parallel"))]
    {
        Ok(sci_form::eht::evaluate_orbital_on_grid(
            &basis,
            &coefficients,
            mo_index,
            &positions,
            spacing,
            3.0,
        ))
    }
}

/// Run an EHT calculation on a molecule.
///
/// - `elements`: JSON array of atomic numbers
/// - `coords_flat`: JSON array of flat xyz coords
/// - `k`: Wolfsberg-Helmholtz constant (0.0 for default 1.75)
#[wasm_bindgen]
pub fn eht_calculate(elements: &str, coords_flat: &str, k: f64) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let k_opt = if k <= 0.0 { None } else { Some(k) };
    match sci_form::eht::solve_eht(&elems, &positions, k_opt) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Run EHT or route to UFF fallback based on support/confidence.
#[wasm_bindgen]
pub fn eht_or_uff_fallback(
    smiles: &str,
    elements: &str,
    coords_flat: &str,
    allow_experimental_eht: bool,
) -> String {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    match sci_form::compute_eht_or_uff_fallback(smiles, &elems, &positions, allow_experimental_eht)
    {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}

/// Generate an orbital isosurface mesh from EHT.
#[wasm_bindgen]
pub fn eht_orbital_mesh(
    elements: &str,
    coords_flat: &str,
    mo_index: usize,
    spacing: f64,
    isovalue: f32,
) -> String {
    let result_json = eht_calculate(elements, coords_flat, 0.0);
    let result: sci_form::eht::EhtResult = match serde_json::from_str(&result_json) {
        Ok(result) => result,
        Err(_) => return result_json,
    };
    let coeff_json = match serde_json::to_string(&result.coefficients) {
        Ok(json) => json,
        Err(e) => return json_error(&e.to_string()),
    };
    let grid =
        match orbital_grid_from_coefficients(elements, coords_flat, &coeff_json, mo_index, spacing)
        {
            Ok(grid) => grid,
            Err(e) => return json_error(&e),
        };
    let mesh = sci_form::eht::marching_cubes(&grid, isovalue);
    serialize_or_error(&mesh)
}

/// Query EHT support metadata.
#[wasm_bindgen]
pub fn eht_support(elements: &str) -> String {
    let elems: Vec<u8> = match parse_elements(elements) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    serialize_or_error(&sci_form::get_eht_support(&elems))
}

/// Compute orbital grid as Float32Array (typed-array transfer).
#[wasm_bindgen]
pub fn eht_orbital_grid_typed(
    elements: &str,
    coords_flat: &str,
    mo_index: usize,
    spacing: f64,
) -> js_sys::Float32Array {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(_) => return js_sys::Float32Array::new_with_length(0),
    };
    let result = match sci_form::eht::solve_eht(&elems, &positions, None) {
        Ok(r) => r,
        Err(_) => return js_sys::Float32Array::new_with_length(0),
    };
    let coeff_json = match serde_json::to_string(&result.coefficients) {
        Ok(json) => json,
        Err(_) => return js_sys::Float32Array::new_with_length(0),
    };
    let grid =
        match orbital_grid_from_coefficients(elements, coords_flat, &coeff_json, mo_index, spacing)
        {
            Ok(grid) => grid,
            Err(_) => return js_sys::Float32Array::new_with_length(0),
        };
    let arr = js_sys::Float32Array::new_with_length(grid.values.len() as u32);
    arr.copy_from(&grid.values);
    arr
}

/// Compute orbital grid from precomputed EHT coefficients as Float32Array.
#[wasm_bindgen]
pub fn eht_orbital_grid_from_coefficients_typed(
    elements: &str,
    coords_flat: &str,
    coefficients_json: &str,
    mo_index: usize,
    spacing: f64,
) -> js_sys::Float32Array {
    let grid = match orbital_grid_from_coefficients(
        elements,
        coords_flat,
        coefficients_json,
        mo_index,
        spacing,
    ) {
        Ok(grid) => grid,
        Err(_) => return js_sys::Float32Array::new_with_length(0),
    };
    let arr = js_sys::Float32Array::new_with_length(grid.values.len() as u32);
    arr.copy_from(&grid.values);
    arr
}

#[wasm_bindgen]
pub fn eht_volumetric_grid_info(
    elements: &str,
    coords_flat: &str,
    spacing: f64,
    padding: f64,
) -> String {
    let (_, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(e) => return json_error(&e),
    };
    let (origin, dims) = sci_form::eht::volume::compute_grid_extents(&positions, padding, spacing);
    grid_info_json(origin, spacing, dims)
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
#[wasm_bindgen]
pub async fn eht_orbital_grid_accelerated(
    elements: &str,
    coords_flat: &str,
    mo_index: usize,
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

    match compute_eht_orbital_grid_wasm_accelerated(
        &elems, &positions, mo_index, spacing, padding, mode,
    )
    .await
    {
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
pub async fn eht_density_grid_accelerated(
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

    match compute_eht_density_grid_wasm_accelerated(&elems, &positions, spacing, padding, mode)
        .await
    {
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
pub async fn eht_orbital_grid_accelerated_typed(
    elements: &str,
    coords_flat: &str,
    mo_index: usize,
    spacing: f64,
    padding: f64,
    execution_mode: &str,
) -> js_sys::Float32Array {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(_) => return js_sys::Float32Array::new_with_length(0),
    };
    let mode = match webgpu::parse_execution_mode(execution_mode) {
        Ok(mode) => mode,
        Err(_) => return js_sys::Float32Array::new_with_length(0),
    };

    match compute_eht_orbital_grid_wasm_accelerated(
        &elems, &positions, mo_index, spacing, padding, mode,
    )
    .await
    {
        Ok((grid, _, _, _)) => {
            let arr = js_sys::Float32Array::new_with_length(grid.values.len() as u32);
            arr.copy_from(&grid.values);
            arr
        }
        Err(_) => js_sys::Float32Array::new_with_length(0),
    }
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
#[wasm_bindgen]
pub async fn eht_density_grid_accelerated_typed(
    elements: &str,
    coords_flat: &str,
    spacing: f64,
    padding: f64,
    execution_mode: &str,
) -> js_sys::Float32Array {
    let (elems, positions) = match parse_elements_and_positions(elements, coords_flat) {
        Ok(v) => v,
        Err(_) => return js_sys::Float32Array::new_with_length(0),
    };
    let mode = match webgpu::parse_execution_mode(execution_mode) {
        Ok(mode) => mode,
        Err(_) => return js_sys::Float32Array::new_with_length(0),
    };

    match compute_eht_density_grid_wasm_accelerated(&elems, &positions, spacing, padding, mode)
        .await
    {
        Ok((grid, _, _, _)) => {
            let arr = js_sys::Float32Array::new_with_length(grid.values.len() as u32);
            arr.copy_from(&grid.values);
            arr
        }
        Err(_) => js_sys::Float32Array::new_with_length(0),
    }
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
pub(crate) async fn compute_eht_orbital_grid_wasm_accelerated(
    elements: &[u8],
    positions: &[[f64; 3]],
    mo_index: usize,
    spacing: f64,
    padding: f64,
    mode: WasmExecutionMode,
) -> Result<(sci_form::eht::VolumetricGrid, String, bool, String), String> {
    let eht_result = sci_form::eht::solve_eht(elements, positions, None)?;
    if mo_index >= eht_result.coefficients.len() {
        return Err(format!(
            "orbital index {} out of range for {} orbitals",
            mo_index,
            eht_result.coefficients.len()
        ));
    }

    if mode == WasmExecutionMode::Cpu {
        let basis = sci_form::eht::basis::build_basis(elements, positions);
        let grid = sci_form::eht::evaluate_orbital_on_grid(
            &basis,
            &eht_result.coefficients,
            mo_index,
            positions,
            spacing,
            padding,
        );
        return Ok((
            grid,
            "CPU".to_string(),
            false,
            cpu_grid_note("CPU-only orbital-grid execution requested"),
        ));
    }

    let Some(runtime) = webgpu::try_runtime().await else {
        if mode == WasmExecutionMode::Gpu {
            return Err("WebGPU not available and mode=gpu requires GPU".to_string());
        }
        let basis = sci_form::eht::basis::build_basis(elements, positions);
        let grid = sci_form::eht::evaluate_orbital_on_grid(
            &basis,
            &eht_result.coefficients,
            mo_index,
            positions,
            spacing,
            padding,
        );
        return Ok((
            grid,
            "CPU".to_string(),
            false,
            cpu_grid_note("WebGPU unavailable; CPU orbital-grid fallback used"),
        ));
    };

    let expanded = expand_eht_basis(
        &sci_form::eht::basis::build_basis(elements, positions),
        &eht_result.coefficients,
    );
    let (origin, dims) = sci_form::eht::volume::compute_grid_extents(positions, padding, spacing);

    if mode == WasmExecutionMode::Hybrid && dims[0] >= 2 {
        let gpu_x = ((dims[0] * 2) / 3).max(1).min(dims[0] - 1);
        let cpu_origin = [origin[0] + gpu_x as f64 * spacing, origin[1], origin[2]];
        let gpu_dims = [gpu_x, dims[1], dims[2]];
        let cpu_dims = [dims[0] - gpu_x, dims[1], dims[2]];

        let mut values =
            compute_eht_orbital_grid_webgpu_slice(&expanded, mo_index, origin, spacing, gpu_dims)
                .await?;
        values.extend(compute_eht_orbital_grid_cpu_slice(
            &expanded, mo_index, cpu_origin, spacing, cpu_dims,
        ));

        return Ok((
            sci_form::eht::VolumetricGrid { origin, spacing, dims, values },
            format!("{}+CPU", runtime.backend),
            true,
            cpu_grid_note(
                "Hybrid orbital grid: WebGPU evaluates the leading x-slabs while CPU computes the remainder",
            ),
        ));
    }

    let values =
        compute_eht_orbital_grid_webgpu_slice(&expanded, mo_index, origin, spacing, dims).await?;
    Ok((
        sci_form::eht::VolumetricGrid {
            origin,
            spacing,
            dims,
            values,
        },
        runtime.backend,
        true,
        "WebGPU orbital-grid dispatch".to_string(),
    ))
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
pub(crate) async fn compute_eht_density_grid_wasm_accelerated(
    elements: &[u8],
    positions: &[[f64; 3]],
    spacing: f64,
    padding: f64,
    mode: WasmExecutionMode,
) -> Result<(sci_form::eht::VolumetricGrid, String, bool, String), String> {
    let eht_result = sci_form::eht::solve_eht(elements, positions, None)?;
    let expanded = expand_eht_basis(
        &sci_form::eht::basis::build_basis(elements, positions),
        &eht_result.coefficients,
    );
    let density = build_density_matrix_from_expanded_coeffs(
        &expanded.coefficients,
        eht_result.n_electrons / 2,
    );
    let (origin, dims) = sci_form::eht::volume::compute_grid_extents(positions, padding, spacing);

    if mode == WasmExecutionMode::Cpu {
        let values = compute_eht_density_grid_cpu_slice(&expanded, &density, origin, spacing, dims);
        return Ok((
            sci_form::eht::VolumetricGrid {
                origin,
                spacing,
                dims,
                values,
            },
            "CPU".to_string(),
            false,
            cpu_grid_note("CPU-only density-grid execution requested"),
        ));
    }

    let Some(runtime) = webgpu::try_runtime().await else {
        if mode == WasmExecutionMode::Gpu {
            return Err("WebGPU not available and mode=gpu requires GPU".to_string());
        }
        let values = compute_eht_density_grid_cpu_slice(&expanded, &density, origin, spacing, dims);
        return Ok((
            sci_form::eht::VolumetricGrid {
                origin,
                spacing,
                dims,
                values,
            },
            "CPU".to_string(),
            false,
            cpu_grid_note("WebGPU unavailable; CPU density-grid fallback used"),
        ));
    };

    if mode == WasmExecutionMode::Hybrid && dims[0] >= 2 {
        let gpu_x = ((dims[0] * 2) / 3).max(1).min(dims[0] - 1);
        let cpu_origin = [origin[0] + gpu_x as f64 * spacing, origin[1], origin[2]];
        let gpu_dims = [gpu_x, dims[1], dims[2]];
        let cpu_dims = [dims[0] - gpu_x, dims[1], dims[2]];

        let mut values =
            compute_eht_density_grid_webgpu_slice(&expanded, &density, origin, spacing, gpu_dims)
                .await?;
        values.extend(compute_eht_density_grid_cpu_slice(
            &expanded, &density, cpu_origin, spacing, cpu_dims,
        ));

        return Ok((
            sci_form::eht::VolumetricGrid { origin, spacing, dims, values },
            format!("{}+CPU", runtime.backend),
            true,
            cpu_grid_note(
                "Hybrid density grid: WebGPU evaluates the leading x-slabs while CPU computes the remainder",
            ),
        ));
    }

    let values =
        compute_eht_density_grid_webgpu_slice(&expanded, &density, origin, spacing, dims).await?;
    Ok((
        sci_form::eht::VolumetricGrid {
            origin,
            spacing,
            dims,
            values,
        },
        runtime.backend,
        true,
        "WebGPU density-grid dispatch".to_string(),
    ))
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
fn expand_eht_basis(
    basis: &[sci_form::eht::basis::AtomicOrbital],
    coefficients: &[Vec<f64>],
) -> ExpandedEhtBasis {
    let mut basis_bytes = Vec::new();
    let mut primitive_bytes = Vec::new();
    let mut expanded_coefficients = Vec::new();

    for (ao_index, orbital) in basis.iter().enumerate() {
        for (term_coefficient, angular) in
            sci_form::eht::basis::orbital_cartesian_terms(orbital.l, orbital.m)
        {
            basis_bytes.extend_from_slice(&(orbital.center[0] as f32).to_ne_bytes());
            basis_bytes.extend_from_slice(&(orbital.center[1] as f32).to_ne_bytes());
            basis_bytes.extend_from_slice(&(orbital.center[2] as f32).to_ne_bytes());
            basis_bytes.extend_from_slice(&(angular[0] as u32).to_ne_bytes());
            basis_bytes.extend_from_slice(&(angular[1] as u32).to_ne_bytes());
            basis_bytes.extend_from_slice(&(angular[2] as u32).to_ne_bytes());
            basis_bytes.extend_from_slice(&(orbital.gaussians.len() as u32).to_ne_bytes());
            basis_bytes.extend_from_slice(&(term_coefficient as f32).to_ne_bytes());

            for primitive_index in 0..3 {
                if primitive_index < orbital.gaussians.len() {
                    let primitive = &orbital.gaussians[primitive_index];
                    let normalized_coeff = primitive.coeff
                        * sci_form::eht::basis::gaussian_cartesian_norm(
                            primitive.alpha,
                            angular[0],
                            angular[1],
                            angular[2],
                        );
                    primitive_bytes.extend_from_slice(&(primitive.alpha as f32).to_ne_bytes());
                    primitive_bytes.extend_from_slice(&(normalized_coeff as f32).to_ne_bytes());
                } else {
                    primitive_bytes.extend_from_slice(&0.0f32.to_ne_bytes());
                    primitive_bytes.extend_from_slice(&0.0f32.to_ne_bytes());
                }
            }

            expanded_coefficients.push(
                coefficients[ao_index]
                    .iter()
                    .map(|value| *value as f32)
                    .collect::<Vec<_>>(),
            );
        }
    }

    ExpandedEhtBasis {
        basis_bytes,
        primitive_bytes,
        coefficients: expanded_coefficients,
    }
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
fn build_density_matrix_from_expanded_coeffs(
    coefficients: &[Vec<f32>],
    n_occupied: usize,
) -> Vec<f32> {
    let n = coefficients.len();
    let mut density = vec![0.0f32; n * n];

    for i in 0..n {
        for j in 0..=i {
            let mut value = 0.0f32;
            for orbital in 0..n_occupied {
                value += coefficients[i][orbital] * coefficients[j][orbital];
            }
            value *= 2.0;
            density[i * n + j] = value;
            density[j * n + i] = value;
        }
    }

    density
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
async fn compute_eht_orbital_grid_webgpu_slice(
    expanded: &ExpandedEhtBasis,
    mo_index: usize,
    origin: [f64; 3],
    spacing: f64,
    dims: [usize; 3],
) -> Result<Vec<f32>, String> {
    let total = dims[0] * dims[1] * dims[2];
    let mo_coefficients: Vec<f32> = expanded
        .coefficients
        .iter()
        .map(|row| row[mo_index])
        .collect();
    let descriptor = ComputeDispatchDescriptor {
        label: "orbital grid wasm".to_string(),
        shader_source: sci_form::gpu::orbital_grid::ORBITAL_GRID_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [
            ceil_div_u32(dims[0], 8),
            ceil_div_u32(dims[1], 8),
            ceil_div_u32(dims[2], 4),
        ],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "basis".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: expanded.basis_bytes.clone(),
            },
            ComputeBindingDescriptor {
                label: "mo_coeffs".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&mo_coefficients),
            },
            ComputeBindingDescriptor {
                label: "primitives".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: expanded.primitive_bytes.clone(),
            },
            ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: pack_uniform_values(&[
                    UniformValue::F32((origin[0] * ANG_TO_BOHR) as f32),
                    UniformValue::F32((origin[1] * ANG_TO_BOHR) as f32),
                    UniformValue::F32((origin[2] * ANG_TO_BOHR) as f32),
                    UniformValue::F32((spacing * ANG_TO_BOHR) as f32),
                    UniformValue::U32(dims[0] as u32),
                    UniformValue::U32(dims[1] as u32),
                    UniformValue::U32(dims[2] as u32),
                    UniformValue::U32(mo_index as u32),
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
    let bytes = outputs
        .pop()
        .ok_or("No output from orbital-grid WebGPU kernel")?;
    Ok(sci_form::gpu::context::bytes_to_f32_vec(&bytes))
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
async fn compute_eht_density_grid_webgpu_slice(
    expanded: &ExpandedEhtBasis,
    density: &[f32],
    origin: [f64; 3],
    spacing: f64,
    dims: [usize; 3],
) -> Result<Vec<f32>, String> {
    let total = dims[0] * dims[1] * dims[2];
    let descriptor = ComputeDispatchDescriptor {
        label: "density grid wasm".to_string(),
        shader_source: sci_form::gpu::density_grid_gpu::DENSITY_GRID_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [
            ceil_div_u32(dims[0], 8),
            ceil_div_u32(dims[1], 8),
            ceil_div_u32(dims[2], 4),
        ],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "basis".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: expanded.basis_bytes.clone(),
            },
            ComputeBindingDescriptor {
                label: "density".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(density),
            },
            ComputeBindingDescriptor {
                label: "primitives".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: expanded.primitive_bytes.clone(),
            },
            ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: pack_uniform_values(&[
                    UniformValue::F32((origin[0] * ANG_TO_BOHR) as f32),
                    UniformValue::F32((origin[1] * ANG_TO_BOHR) as f32),
                    UniformValue::F32((origin[2] * ANG_TO_BOHR) as f32),
                    UniformValue::F32((spacing * ANG_TO_BOHR) as f32),
                    UniformValue::U32(dims[0] as u32),
                    UniformValue::U32(dims[1] as u32),
                    UniformValue::U32(dims[2] as u32),
                    UniformValue::U32(expanded.coefficients.len() as u32),
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
    let bytes = outputs
        .pop()
        .ok_or("No output from density-grid WebGPU kernel")?;
    Ok(sci_form::gpu::context::bytes_to_f32_vec(&bytes))
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
fn compute_eht_orbital_grid_cpu_slice(
    expanded: &ExpandedEhtBasis,
    mo_index: usize,
    origin: [f64; 3],
    spacing: f64,
    dims: [usize; 3],
) -> Vec<f32> {
    let total = dims[0] * dims[1] * dims[2];

    #[cfg(feature = "parallel")]
    {
        return (0..total)
            .into_par_iter()
            .map(|flat_idx| {
                orbital_value_at_index(flat_idx, expanded, mo_index, origin, spacing, dims)
            })
            .collect();
    }

    #[cfg(not(feature = "parallel"))]
    {
        (0..total)
            .map(|flat_idx| {
                orbital_value_at_index(flat_idx, expanded, mo_index, origin, spacing, dims)
            })
            .collect()
    }
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
fn compute_eht_density_grid_cpu_slice(
    expanded: &ExpandedEhtBasis,
    density: &[f32],
    origin: [f64; 3],
    spacing: f64,
    dims: [usize; 3],
) -> Vec<f32> {
    let total = dims[0] * dims[1] * dims[2];

    #[cfg(feature = "parallel")]
    {
        return (0..total)
            .into_par_iter()
            .map(|flat_idx| {
                density_value_at_index(flat_idx, expanded, density, origin, spacing, dims)
            })
            .collect();
    }

    #[cfg(not(feature = "parallel"))]
    {
        (0..total)
            .map(|flat_idx| {
                density_value_at_index(flat_idx, expanded, density, origin, spacing, dims)
            })
            .collect()
    }
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
fn orbital_value_at_index(
    flat_idx: usize,
    expanded: &ExpandedEhtBasis,
    mo_index: usize,
    origin: [f64; 3],
    spacing: f64,
    dims: [usize; 3],
) -> f32 {
    let point = point_from_flat_index(flat_idx, origin, spacing, dims);
    let point_bohr = [
        point[0] * ANG_TO_BOHR,
        point[1] * ANG_TO_BOHR,
        point[2] * ANG_TO_BOHR,
    ];
    let mut psi = 0.0f64;
    let n_basis = expanded.coefficients.len();

    for mu in 0..n_basis {
        let coefficient = expanded.coefficients[mu][mo_index] as f64;
        if coefficient.abs() < 1e-7 {
            continue;
        }
        psi += coefficient * expanded_basis_value(expanded, mu, &point_bohr);
    }

    psi as f32
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
fn density_value_at_index(
    flat_idx: usize,
    expanded: &ExpandedEhtBasis,
    density: &[f32],
    origin: [f64; 3],
    spacing: f64,
    dims: [usize; 3],
) -> f32 {
    let point = point_from_flat_index(flat_idx, origin, spacing, dims);
    let point_bohr = [
        point[0] * ANG_TO_BOHR,
        point[1] * ANG_TO_BOHR,
        point[2] * ANG_TO_BOHR,
    ];
    let n_basis = expanded.coefficients.len();
    let mut phi = vec![0.0f64; n_basis];
    for mu in 0..n_basis {
        phi[mu] = expanded_basis_value(expanded, mu, &point_bohr);
    }

    let mut rho = 0.0f64;
    for mu in 0..n_basis {
        if phi[mu].abs() < 1e-10 {
            continue;
        }
        for nu in 0..n_basis {
            rho += density[mu * n_basis + nu] as f64 * phi[mu] * phi[nu];
        }
    }

    rho as f32
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
fn expanded_basis_value(
    expanded: &ExpandedEhtBasis,
    basis_index: usize,
    point_bohr: &[f64; 3],
) -> f64 {
    let basis_offset = basis_index * 32;
    let primitive_offset = basis_index * 24;

    let center_x = f32::from_ne_bytes(
        expanded.basis_bytes[basis_offset..basis_offset + 4]
            .try_into()
            .unwrap(),
    ) as f64;
    let center_y = f32::from_ne_bytes(
        expanded.basis_bytes[basis_offset + 4..basis_offset + 8]
            .try_into()
            .unwrap(),
    ) as f64;
    let center_z = f32::from_ne_bytes(
        expanded.basis_bytes[basis_offset + 8..basis_offset + 12]
            .try_into()
            .unwrap(),
    ) as f64;
    let lx = u32::from_ne_bytes(
        expanded.basis_bytes[basis_offset + 12..basis_offset + 16]
            .try_into()
            .unwrap(),
    ) as i32;
    let ly = u32::from_ne_bytes(
        expanded.basis_bytes[basis_offset + 16..basis_offset + 20]
            .try_into()
            .unwrap(),
    ) as i32;
    let lz = u32::from_ne_bytes(
        expanded.basis_bytes[basis_offset + 20..basis_offset + 24]
            .try_into()
            .unwrap(),
    ) as i32;
    let n_primitives = u32::from_ne_bytes(
        expanded.basis_bytes[basis_offset + 24..basis_offset + 28]
            .try_into()
            .unwrap(),
    ) as usize;
    let term_coefficient = f32::from_ne_bytes(
        expanded.basis_bytes[basis_offset + 28..basis_offset + 32]
            .try_into()
            .unwrap(),
    ) as f64;

    let dx = point_bohr[0] - center_x;
    let dy = point_bohr[1] - center_y;
    let dz = point_bohr[2] - center_z;
    let r2 = dx * dx + dy * dy + dz * dz;
    let angular = dx.powi(lx) * dy.powi(ly) * dz.powi(lz);

    let mut radial = 0.0f64;
    for primitive_index in 0..n_primitives {
        let offset = primitive_offset + primitive_index * 8;
        let alpha = f32::from_ne_bytes(
            expanded.primitive_bytes[offset..offset + 4]
                .try_into()
                .unwrap(),
        ) as f64;
        let coefficient = f32::from_ne_bytes(
            expanded.primitive_bytes[offset + 4..offset + 8]
                .try_into()
                .unwrap(),
        ) as f64;
        radial += coefficient * (-alpha * r2).exp();
    }

    term_coefficient * angular * radial
}

#[cfg(all(feature = "experimental-gpu", target_arch = "wasm32"))]
fn point_from_flat_index(
    flat_idx: usize,
    origin: [f64; 3],
    spacing: f64,
    dims: [usize; 3],
) -> [f64; 3] {
    let ny = dims[1];
    let nz = dims[2];
    let ix = flat_idx / (ny * nz);
    let iy = (flat_idx / nz) % ny;
    let iz = flat_idx % nz;
    [
        origin[0] + ix as f64 * spacing,
        origin[1] + iy as f64 * spacing,
        origin[2] + iz as f64 * spacing,
    ]
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
