//! GPU electron-density grid evaluation.

use nalgebra::DMatrix;

use super::backend_report::OrbitalGridReport;
use super::context::{
    bytes_to_f32_vec, f32_slice_to_bytes, ComputeBindingDescriptor, ComputeBindingKind,
    ComputeDispatchDescriptor, GpuContext,
};
use crate::gpu::orbital_grid::{evaluate_density_cpu, GridParams};
use crate::scf::basis::{BasisFunction, BasisSet};

fn pack_basis_for_gpu(basis: &BasisSet) -> (Vec<u8>, Vec<u8>) {
    let mut basis_bytes = Vec::new();
    let mut prim_bytes = Vec::new();

    for bf in &basis.functions {
        basis_bytes.extend_from_slice(&(bf.center[0] as f32).to_ne_bytes());
        basis_bytes.extend_from_slice(&(bf.center[1] as f32).to_ne_bytes());
        basis_bytes.extend_from_slice(&(bf.center[2] as f32).to_ne_bytes());
        basis_bytes.extend_from_slice(&bf.angular[0].to_ne_bytes());
        basis_bytes.extend_from_slice(&bf.angular[1].to_ne_bytes());
        basis_bytes.extend_from_slice(&bf.angular[2].to_ne_bytes());
        basis_bytes.extend_from_slice(&(bf.primitives.len() as u32).to_ne_bytes());
        let norm = BasisFunction::normalization(
            bf.primitives.first().map(|p| p.alpha).unwrap_or(1.0),
            bf.angular[0],
            bf.angular[1],
            bf.angular[2],
        );
        basis_bytes.extend_from_slice(&(norm as f32).to_ne_bytes());

        for index in 0..3 {
            if index < bf.primitives.len() {
                prim_bytes.extend_from_slice(&(bf.primitives[index].alpha as f32).to_ne_bytes());
                prim_bytes
                    .extend_from_slice(&(bf.primitives[index].coefficient as f32).to_ne_bytes());
            } else {
                prim_bytes.extend_from_slice(&0.0f32.to_ne_bytes());
                prim_bytes.extend_from_slice(&0.0f32.to_ne_bytes());
            }
        }
    }

    (basis_bytes, prim_bytes)
}

pub fn evaluate_density_with_report(
    basis: &BasisSet,
    density: &DMatrix<f64>,
    params: &GridParams,
) -> (Vec<f64>, OrbitalGridReport) {
    let ctx = GpuContext::best_available();
    if ctx.is_gpu_available() {
        match evaluate_density_gpu(&ctx, basis, density, params) {
            Ok(grid) => {
                return (
                    grid,
                    OrbitalGridReport {
                        backend: ctx.capabilities.backend.clone(),
                        used_gpu: true,
                        attempted_gpu: true,
                        n_points: params.n_points(),
                        note: format!("GPU density-grid dispatch on {}", ctx.capabilities.backend),
                    },
                );
            }
            Err(_err) => {}
        }
    }

    let grid = evaluate_density_cpu(basis, density, params);
    (
        grid,
        OrbitalGridReport {
            backend: "CPU".to_string(),
            used_gpu: false,
            attempted_gpu: ctx.is_gpu_available(),
            n_points: params.n_points(),
            note: if ctx.is_gpu_available() {
                "GPU available but density-grid dispatch failed; CPU fallback used".to_string()
            } else {
                "CPU density-grid evaluation (GPU not available)".to_string()
            },
        },
    )
}

pub fn evaluate_density_gpu(
    ctx: &GpuContext,
    basis: &BasisSet,
    density: &DMatrix<f64>,
    params: &GridParams,
) -> Result<Vec<f64>, String> {
    let n_basis = basis.n_basis;
    let n_points = params.n_points();
    let (basis_bytes, prim_bytes) = pack_basis_for_gpu(basis);
    let density_flat: Vec<f32> = (0..n_basis)
        .flat_map(|mu| (0..n_basis).map(move |nu| density[(mu, nu)] as f32))
        .collect();

    let mut params_bytes = Vec::with_capacity(32);
    for value in &params.origin {
        params_bytes.extend_from_slice(&(*value as f32).to_ne_bytes());
    }
    params_bytes.extend_from_slice(&(params.spacing as f32).to_ne_bytes());
    for dim in &params.dimensions {
        params_bytes.extend_from_slice(&(*dim as u32).to_ne_bytes());
    }
    params_bytes.extend_from_slice(&(n_basis as u32).to_ne_bytes());

    let [nx, ny, nz] = params.dimensions;
    let descriptor = ComputeDispatchDescriptor {
        label: "density grid".to_string(),
        shader_source: DENSITY_GRID_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [
            (nx as u32).div_ceil(8),
            (ny as u32).div_ceil(8),
            (nz as u32).div_ceil(4),
        ],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "basis".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: basis_bytes,
            },
            ComputeBindingDescriptor {
                label: "density".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&density_flat),
            },
            ComputeBindingDescriptor {
                label: "primitives".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: prim_bytes,
            },
            ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: params_bytes,
            },
            ComputeBindingDescriptor {
                label: "output".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: f32_slice_to_bytes(&vec![0.0f32; n_points]),
            },
        ],
    };

    let mut outputs = ctx.run_compute(&descriptor)?.outputs;
    let bytes = outputs.pop().ok_or("No output from density grid kernel")?;
    let values = bytes_to_f32_vec(&bytes);
    if values.len() != n_points {
        return Err(format!(
            "Output size mismatch: expected {}, got {}",
            n_points,
            values.len()
        ));
    }
    Ok(values.into_iter().map(|value| value as f64).collect())
}

pub const DENSITY_GRID_SHADER: &str = r#"
struct BasisFunc {
    center_x: f32, center_y: f32, center_z: f32,
    lx: u32, ly: u32, lz: u32,
    n_primitives: u32,
    norm_coeff: f32,
};

struct GridParams {
    origin_x: f32, origin_y: f32, origin_z: f32,
    spacing: f32,
    dims_x: u32, dims_y: u32, dims_z: u32,
    n_basis: u32,
};

@group(0) @binding(0) var<storage, read> basis: array<BasisFunc>;
@group(0) @binding(1) var<storage, read> density: array<f32>;
@group(0) @binding(2) var<storage, read> primitives: array<vec2<f32>>;
@group(0) @binding(3) var<uniform> params: GridParams;
@group(0) @binding(4) var<storage, read_write> output: array<f32>;

fn phi_at(mu: u32, rx: f32, ry: f32, rz: f32) -> f32 {
    let bf = basis[mu];
    let dx = rx - bf.center_x;
    let dy = ry - bf.center_y;
    let dz = rz - bf.center_z;
    let r2 = dx * dx + dy * dy + dz * dz;

    var angular: f32 = 1.0;
    for (var i: u32 = 0u; i < bf.lx; i = i + 1u) { angular *= dx; }
    for (var i: u32 = 0u; i < bf.ly; i = i + 1u) { angular *= dy; }
    for (var i: u32 = 0u; i < bf.lz; i = i + 1u) { angular *= dz; }

    var radial: f32 = 0.0;
    for (var p: u32 = 0u; p < bf.n_primitives; p = p + 1u) {
        let prim = primitives[mu * 3u + p];
        radial += prim.y * exp(-prim.x * r2);
    }
    return bf.norm_coeff * angular * radial;
}

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

    var rho: f32 = 0.0;
    for (var mu: u32 = 0u; mu < params.n_basis; mu = mu + 1u) {
        let phi_mu = phi_at(mu, rx, ry, rz);
        if (abs(phi_mu) < 1e-7) { continue; }
        for (var nu: u32 = 0u; nu < params.n_basis; nu = nu + 1u) {
            let phi_nu = phi_at(nu, rx, ry, rz);
            rho += density[mu * params.n_basis + nu] * phi_mu * phi_nu;
        }
    }

    output[flat_idx] = rho;
}
"#;
