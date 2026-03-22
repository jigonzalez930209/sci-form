//! Molecular orbital evaluation on 3D grids (GPU + CPU).
//!
//! Computes ψ_i(r) = Σ_μ C_{μi} φ_μ(r) on a regular 3D grid.
//! Each grid point can be evaluated independently — ideal for GPU.
//!
//! GPU path: dispatches ORBITAL_GRID_SHADER via wgpu.
//! CPU path: direct triple-nested loop (always available as fallback).

use super::backend_report::OrbitalGridReport;
use super::context::{
    bytes_to_f32_vec, f32_slice_to_bytes, ComputeBindingDescriptor, ComputeBindingKind,
    ComputeDispatchDescriptor, GpuContext,
};
use crate::scf::basis::{BasisFunction, BasisSet};
use nalgebra::DMatrix;

/// 3D grid parameters.
#[derive(Debug, Clone)]
pub struct GridParams {
    /// Grid origin (x, y, z) in Bohr.
    pub origin: [f64; 3],
    /// Grid spacing in Bohr.
    pub spacing: f64,
    /// Number of grid points [nx, ny, nz].
    pub dimensions: [usize; 3],
}

impl GridParams {
    /// Create grid params enclosing the molecule with padding.
    pub fn from_molecule(positions: &[[f64; 3]], spacing: f64, padding: f64) -> Self {
        let mut min = [f64::MAX; 3];
        let mut max = [f64::MIN; 3];

        for pos in positions {
            for k in 0..3 {
                min[k] = min[k].min(pos[k]);
                max[k] = max[k].max(pos[k]);
            }
        }

        let origin = [min[0] - padding, min[1] - padding, min[2] - padding];
        let dimensions = [
            ((max[0] - min[0] + 2.0 * padding) / spacing).ceil() as usize + 1,
            ((max[1] - min[1] + 2.0 * padding) / spacing).ceil() as usize + 1,
            ((max[2] - min[2] + 2.0 * padding) / spacing).ceil() as usize + 1,
        ];

        Self {
            origin,
            spacing,
            dimensions,
        }
    }

    /// Total number of grid points.
    pub fn n_points(&self) -> usize {
        self.dimensions[0] * self.dimensions[1] * self.dimensions[2]
    }

    /// 3D coordinate of grid point (ix, iy, iz).
    pub fn point(&self, ix: usize, iy: usize, iz: usize) -> [f64; 3] {
        [
            self.origin[0] + ix as f64 * self.spacing,
            self.origin[1] + iy as f64 * self.spacing,
            self.origin[2] + iz as f64 * self.spacing,
        ]
    }

    /// Flat index from 3D indices.
    pub fn flat_index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix * self.dimensions[1] * self.dimensions[2] + iy * self.dimensions[2] + iz
    }
}

/// Result of orbital grid evaluation.
#[derive(Debug, Clone)]
pub struct OrbitalGrid {
    /// Grid values (flat, row-major: x varies slowest).
    pub values: Vec<f64>,
    pub params: GridParams,
    pub orbital_index: usize,
}

/// Evaluate a molecular orbital on a 3D grid with explicit backend reporting.
///
/// Attempts GPU dispatch when available; falls back to CPU otherwise.
pub fn evaluate_orbital_with_report(
    basis: &BasisSet,
    mo_coefficients: &DMatrix<f64>,
    orbital_index: usize,
    params: &GridParams,
) -> (OrbitalGrid, OrbitalGridReport) {
    let ctx = GpuContext::best_available();

    if ctx.is_gpu_available() {
        match evaluate_orbital_gpu(&ctx, basis, mo_coefficients, orbital_index, params) {
            Ok(grid) => {
                let report = OrbitalGridReport {
                    backend: ctx.capabilities.backend.clone(),
                    used_gpu: true,
                    attempted_gpu: true,
                    n_points: params.n_points(),
                    note: format!("GPU dispatch on {}", ctx.capabilities.backend),
                };
                return (grid, report);
            }
            Err(_err) => {
                // Fall through to CPU
            }
        }
    }

    let grid = evaluate_orbital_cpu(basis, mo_coefficients, orbital_index, params);
    let report = OrbitalGridReport {
        backend: "CPU".to_string(),
        used_gpu: false,
        attempted_gpu: ctx.is_gpu_available(),
        n_points: params.n_points(),
        note: if ctx.is_gpu_available() {
            "GPU available but dispatch failed; CPU fallback used".to_string()
        } else {
            "CPU evaluation (GPU not available)".to_string()
        },
    };
    (grid, report)
}

/// Evaluate orbital on CPU (always available).
pub fn evaluate_orbital_cpu(
    basis: &BasisSet,
    mo_coefficients: &DMatrix<f64>,
    orbital_index: usize,
    params: &GridParams,
) -> OrbitalGrid {
    let n_points = params.n_points();
    let mut values = vec![0.0; n_points];
    let n_basis = basis.n_basis;
    let [nx, ny, nz] = params.dimensions;

    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let r = params.point(ix, iy, iz);
                let idx = params.flat_index(ix, iy, iz);

                let mut psi = 0.0;
                for mu in 0..n_basis {
                    let c_mu = mo_coefficients[(mu, orbital_index)];
                    if c_mu.abs() < 1e-15 {
                        continue;
                    }
                    let phi_mu = evaluate_basis_function(&basis.functions[mu], &r);
                    psi += c_mu * phi_mu;
                }
                values[idx] = psi;
            }
        }
    }

    OrbitalGrid {
        values,
        params: params.clone(),
        orbital_index,
    }
}

/// Evaluate electron density ρ(r) = Σ_{μν} P_{μν} φ_μ(r) φ_ν(r) on a 3D grid.
pub fn evaluate_density_cpu(
    basis: &BasisSet,
    density: &DMatrix<f64>,
    params: &GridParams,
) -> Vec<f64> {
    let n_points = params.n_points();
    let mut values = vec![0.0; n_points];
    let n_basis = basis.n_basis;
    let [nx, ny, nz] = params.dimensions;

    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let r = params.point(ix, iy, iz);
                let idx = params.flat_index(ix, iy, iz);

                let phi: Vec<f64> = (0..n_basis)
                    .map(|mu| evaluate_basis_function(&basis.functions[mu], &r))
                    .collect();

                let mut rho = 0.0;
                for mu in 0..n_basis {
                    if phi[mu].abs() < 1e-15 {
                        continue;
                    }
                    for nu in 0..n_basis {
                        rho += density[(mu, nu)] * phi[mu] * phi[nu];
                    }
                }
                values[idx] = rho;
            }
        }
    }
    values
}

/// Evaluate a single contracted Gaussian basis function at point r.
fn evaluate_basis_function(bf: &BasisFunction, r: &[f64; 3]) -> f64 {
    let dx = r[0] - bf.center[0];
    let dy = r[1] - bf.center[1];
    let dz = r[2] - bf.center[2];
    let r2 = dx * dx + dy * dy + dz * dz;

    let angular = dx.powi(bf.angular[0] as i32)
        * dy.powi(bf.angular[1] as i32)
        * dz.powi(bf.angular[2] as i32);

    let mut radial = 0.0;
    for prim in &bf.primitives {
        radial += prim.coefficient * (-prim.alpha * r2).exp();
    }

    BasisFunction::normalization(
        bf.primitives.first().map(|p| p.alpha).unwrap_or(1.0),
        bf.angular[0],
        bf.angular[1],
        bf.angular[2],
    ) * angular
        * radial
}

// ─── GPU dispatch ────────────────────────────────────────────────────────────

/// Pack basis function data for the GPU shader.
///
/// Each basis function → GpuBasisFunc (32 bytes):
///   center: vec3<f32>, lx: u32, ly: u32, lz: u32, n_primitives: u32, coefficient: f32
///
/// Primitives → (alpha: f32, coeff: f32) pairs, max 3 per basis function (STO-3G).
fn pack_basis_for_gpu(basis: &BasisSet) -> (Vec<u8>, Vec<u8>) {
    let mut basis_bytes = Vec::new();
    let mut prim_bytes = Vec::new();

    for bf in &basis.functions {
        // center xyz
        basis_bytes.extend_from_slice(&(bf.center[0] as f32).to_ne_bytes());
        basis_bytes.extend_from_slice(&(bf.center[1] as f32).to_ne_bytes());
        basis_bytes.extend_from_slice(&(bf.center[2] as f32).to_ne_bytes());
        // lx, ly, lz
        basis_bytes.extend_from_slice(&bf.angular[0].to_ne_bytes());
        basis_bytes.extend_from_slice(&bf.angular[1].to_ne_bytes());
        basis_bytes.extend_from_slice(&bf.angular[2].to_ne_bytes());
        // n_primitives
        basis_bytes.extend_from_slice(&(bf.primitives.len() as u32).to_ne_bytes());
        // normalization coefficient
        let norm = BasisFunction::normalization(
            bf.primitives.first().map(|p| p.alpha).unwrap_or(1.0),
            bf.angular[0],
            bf.angular[1],
            bf.angular[2],
        );
        basis_bytes.extend_from_slice(&(norm as f32).to_ne_bytes());

        // Pack primitives (max 3 for STO-3G)
        for i in 0..3 {
            if i < bf.primitives.len() {
                prim_bytes.extend_from_slice(&(bf.primitives[i].alpha as f32).to_ne_bytes());
                prim_bytes.extend_from_slice(&(bf.primitives[i].coefficient as f32).to_ne_bytes());
            } else {
                prim_bytes.extend_from_slice(&0.0f32.to_ne_bytes());
                prim_bytes.extend_from_slice(&0.0f32.to_ne_bytes());
            }
        }
    }

    (basis_bytes, prim_bytes)
}

/// GPU-accelerated orbital grid evaluation.
fn evaluate_orbital_gpu(
    ctx: &GpuContext,
    basis: &BasisSet,
    mo_coefficients: &DMatrix<f64>,
    orbital_index: usize,
    params: &GridParams,
) -> Result<OrbitalGrid, String> {
    let n_basis = basis.n_basis;
    let n_points = params.n_points();

    // Pack basis functions and primitives
    let (basis_bytes, prim_bytes) = pack_basis_for_gpu(basis);

    // Pack MO coefficients for this orbital
    let mo_coeffs: Vec<f32> = (0..n_basis)
        .map(|mu| mo_coefficients[(mu, orbital_index)] as f32)
        .collect();

    // Pack grid params: origin (3×f32) + spacing (f32) + dims (3×u32) + orbital_index (u32) = 32 bytes
    let mut params_bytes = Vec::with_capacity(32);
    for v in &params.origin {
        params_bytes.extend_from_slice(&(*v as f32).to_ne_bytes());
    }
    params_bytes.extend_from_slice(&(params.spacing as f32).to_ne_bytes());
    for d in &params.dimensions {
        params_bytes.extend_from_slice(&(*d as u32).to_ne_bytes());
    }
    params_bytes.extend_from_slice(&(orbital_index as u32).to_ne_bytes());

    // Output buffer
    let output_seed = vec![0.0f32; n_points];

    let [nx, ny, nz] = params.dimensions;
    let wg = [
        (nx as u32).div_ceil(8),
        (ny as u32).div_ceil(8),
        (nz as u32).div_ceil(4),
    ];

    let descriptor = ComputeDispatchDescriptor {
        label: "orbital grid".to_string(),
        shader_source: ORBITAL_GRID_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: wg,
        bindings: vec![
            ComputeBindingDescriptor {
                label: "basis".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: basis_bytes,
            },
            ComputeBindingDescriptor {
                label: "mo_coeffs".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&mo_coeffs),
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
                bytes: f32_slice_to_bytes(&output_seed),
            },
        ],
    };

    let mut result = ctx.run_compute(&descriptor)?.outputs;
    let bytes = result.pop().ok_or("No output from orbital grid kernel")?;
    let f32_values = bytes_to_f32_vec(&bytes);

    if f32_values.len() != n_points {
        return Err(format!(
            "Output size mismatch: expected {}, got {}",
            n_points,
            f32_values.len()
        ));
    }

    let values: Vec<f64> = f32_values.iter().map(|v| *v as f64).collect();

    Ok(OrbitalGrid {
        values,
        params: params.clone(),
        orbital_index,
    })
}

/// WGSL compute shader for orbital grid evaluation.
///
/// Evaluates ψ_i(r) = Σ_μ C_{μi} φ_μ(r) at each grid point.
/// Workgroup size: (8, 8, 4) = 256 threads.
pub const ORBITAL_GRID_SHADER: &str = r#"
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
    orbital_index: u32,
};

@group(0) @binding(0) var<storage, read> basis: array<BasisFunc>;
@group(0) @binding(1) var<storage, read> mo_coeffs: array<f32>;
@group(0) @binding(2) var<storage, read> primitives: array<vec2<f32>>;
@group(0) @binding(3) var<uniform> params: GridParams;
@group(0) @binding(4) var<storage, read_write> output: array<f32>;

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
    let n_basis = arrayLength(&mo_coeffs);

    var psi: f32 = 0.0;

    for (var mu: u32 = 0u; mu < n_basis; mu = mu + 1u) {
        let c_mu = mo_coeffs[mu];
        if (abs(c_mu) < 1e-7) {
            continue;
        }

        let bf = basis[mu];
        let dx = rx - bf.center_x;
        let dy = ry - bf.center_y;
        let dz = rz - bf.center_z;
        let r2 = dx * dx + dy * dy + dz * dz;

        // Angular part
        var angular: f32 = 1.0;
        for (var i: u32 = 0u; i < bf.lx; i = i + 1u) { angular *= dx; }
        for (var i: u32 = 0u; i < bf.ly; i = i + 1u) { angular *= dy; }
        for (var i: u32 = 0u; i < bf.lz; i = i + 1u) { angular *= dz; }

        // Radial part (contracted, max 3 primitives for STO-3G)
        var radial: f32 = 0.0;
        for (var p: u32 = 0u; p < bf.n_primitives; p = p + 1u) {
            let prim = primitives[mu * 3u + p];
            radial += prim.y * exp(-prim.x * r2);
        }

        psi += c_mu * bf.norm_coeff * angular * radial;
    }

    output[flat_idx] = psi;
}
"#;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_grid_params_from_molecule() {
        let positions = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let params = GridParams::from_molecule(&positions, 0.5, 3.0);
        assert!(params.dimensions[0] > 0);
        assert!(params.n_points() > 0);
        assert!(params.origin[0] < -2.0);
    }

    #[test]
    fn test_grid_point_coordinates() {
        let params = GridParams {
            origin: [0.0, 0.0, 0.0],
            spacing: 1.0,
            dimensions: [3, 3, 3],
        };
        let p = params.point(1, 2, 0);
        assert!((p[0] - 1.0).abs() < 1e-12);
        assert!((p[1] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_flat_index() {
        let params = GridParams {
            origin: [0.0, 0.0, 0.0],
            spacing: 1.0,
            dimensions: [3, 4, 5],
        };
        assert_eq!(params.flat_index(0, 0, 0), 0);
        assert_eq!(params.flat_index(0, 0, 1), 1);
        assert_eq!(params.flat_index(0, 1, 0), 5);
        assert_eq!(params.flat_index(1, 0, 0), 20);
    }

    #[test]
    fn test_evaluate_orbital_cpu_h2() {
        // Build H2 basis
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]]; // ~0.74 Å in Bohr
        let basis = BasisSet::sto3g(&elements, &positions);

        // Simple MO coefficients (bonding orbital: equal contribution)
        let n = basis.n_basis;
        let mut coeffs = DMatrix::zeros(n, n);
        let c = 1.0 / (2.0f64).sqrt();
        coeffs[(0, 0)] = c;
        if n > 1 {
            coeffs[(1, 0)] = c;
        }

        let params = GridParams {
            origin: [-2.0, -2.0, -2.0],
            spacing: 0.5,
            dimensions: [5, 5, 13],
        };

        let grid = evaluate_orbital_cpu(&basis, &coeffs, 0, &params);
        assert_eq!(grid.values.len(), params.n_points());

        // ψ should be non-zero near the bond axis
        let center_idx = params.flat_index(2, 2, 5); // near midpoint
        assert!(grid.values[center_idx].abs() > 1e-6);
    }

    #[test]
    fn test_evaluate_orbital_with_report() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]];
        let basis = BasisSet::sto3g(&elements, &positions);

        let n = basis.n_basis;
        let mut coeffs = DMatrix::zeros(n, n);
        coeffs[(0, 0)] = 1.0 / (2.0f64).sqrt();
        if n > 1 {
            coeffs[(1, 0)] = 1.0 / (2.0f64).sqrt();
        }

        let params = GridParams {
            origin: [-1.0, -1.0, -1.0],
            spacing: 1.0,
            dimensions: [3, 3, 5],
        };

        let (grid, report) = evaluate_orbital_with_report(&basis, &coeffs, 0, &params);
        assert_eq!(grid.values.len(), params.n_points());
        assert!(!report.backend.is_empty());
        assert_eq!(report.n_points, params.n_points());
    }

    #[test]
    fn test_evaluate_density_cpu() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]];
        let basis = BasisSet::sto3g(&elements, &positions);

        let n = basis.n_basis;
        // Simple density matrix
        let density = DMatrix::from_fn(n, n, |i, j| if i == j { 1.0 } else { 0.3 });

        let params = GridParams {
            origin: [-1.0, -1.0, -1.0],
            spacing: 1.0,
            dimensions: [3, 3, 4],
        };

        let values = evaluate_density_cpu(&basis, &density, &params);
        assert_eq!(values.len(), params.n_points());
        // Density should be non-negative at most points (positive definite P)
    }

    #[test]
    fn test_pack_basis_for_gpu() {
        let elements = [1u8];
        let positions = [[0.0, 0.0, 0.0]];
        let basis = BasisSet::sto3g(&elements, &positions);

        let (basis_bytes, prim_bytes) = pack_basis_for_gpu(&basis);
        // Each basis function: 8 × f32/u32 = 32 bytes
        assert_eq!(basis_bytes.len(), basis.n_basis * 32);
        // Each basis function: 3 primitives × 2 × f32 = 24 bytes
        assert_eq!(prim_bytes.len(), basis.n_basis * 24);
    }
}
