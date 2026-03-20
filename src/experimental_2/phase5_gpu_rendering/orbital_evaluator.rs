//! Molecular orbital evaluation on 3D grids.
//!
//! Computes ψ_i(r) = Σ_μ C_{μi} φ_μ(r) on a regular 3D grid,
//! where φ_μ are contracted Gaussian basis functions.
//!
//! Each Gaussian primitive contributes:
//!   g(r) = N · x^{lx} y^{ly} z^{lz} · exp(-α|r - R|²)
//!
//! The grid evaluation is designed for GPU parallelization:
//! each grid point can be evaluated independently.

use crate::experimental_2::phase2_quantum_engine::basis_set::BasisSet;


/// 3D grid parameters.
#[derive(Debug, Clone)]
pub struct GridParams {
    /// Grid origin (x, y, z) in Bohr.
    pub origin: [f64; 3],
    /// Grid spacing in Bohr.
    pub spacing: f64,
    /// Number of grid points in each dimension.
    pub dimensions: [usize; 3],
}

impl GridParams {
    /// Create grid params that enclose the molecule with padding.
    pub fn from_molecule(positions: &[[f64; 3]], spacing: f64, padding: f64) -> Self {
        let mut min = [f64::MAX; 3];
        let mut max = [f64::MIN; 3];

        for pos in positions {
            for k in 0..3 {
                min[k] = min[k].min(pos[k]);
                max[k] = max[k].max(pos[k]);
            }
        }

        let origin = [
            min[0] - padding,
            min[1] - padding,
            min[2] - padding,
        ];

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

    /// Get the 3D coordinate of grid point (ix, iy, iz).
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
    /// Grid values (flat array, row-major: x varies slowest).
    pub values: Vec<f64>,
    /// Grid parameters.
    pub params: GridParams,
    /// Index of the evaluated orbital.
    pub orbital_index: usize,
}

/// Evaluate a molecular orbital on a 3D grid.
///
///   ψ_i(r) = Σ_μ C_{μi} φ_μ(r)
///
/// where φ_μ(r) = Σ_p c_p · N · (x-Rx)^lx · (y-Ry)^ly · (z-Rz)^lz · exp(-α_p |r-R|²)
pub fn evaluate_orbital(
    basis: &BasisSet,
    mo_coefficients: &nalgebra::DMatrix<f64>,
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

/// Evaluate electron density on a 3D grid.
///
///   ρ(r) = Σ_μν P_μν φ_μ(r) φ_ν(r)
///
/// This is more expensive than orbital evaluation due to the
/// double sum over basis functions.
pub fn evaluate_density(
    basis: &BasisSet,
    density: &nalgebra::DMatrix<f64>,
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

                // Pre-evaluate all basis functions at this point
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
///
///   φ_μ(r) = N · Σ_p c_p · (x-Rx)^lx · (y-Ry)^ly · (z-Rz)^lz · exp(-α_p |r-R|²)
fn evaluate_basis_function(
    bf: &crate::experimental_2::phase2_quantum_engine::basis_set::BasisFunction,
    r: &[f64; 3],
) -> f64 {
    let dx = r[0] - bf.center[0];
    let dy = r[1] - bf.center[1];
    let dz = r[2] - bf.center[2];
    let r2 = dx * dx + dy * dy + dz * dz;

    // Angular part: x^lx * y^ly * z^lz
    let angular = dx.powi(bf.angular[0] as i32)
        * dy.powi(bf.angular[1] as i32)
        * dz.powi(bf.angular[2] as i32);

    // Radial part: Σ c_p · exp(-α_p · r²)
    let mut radial = 0.0;
    for prim in &bf.primitives {
        radial += prim.coefficient * (-prim.alpha * r2).exp();
    }

    use crate::experimental_2::phase2_quantum_engine::basis_set::BasisFunction;
    BasisFunction::normalization(
        bf.primitives.first().map(|p| p.alpha).unwrap_or(1.0),
        bf.angular[0], bf.angular[1], bf.angular[2],
    ) * angular * radial
}

/// WGSL compute shader source for orbital grid evaluation.
///
/// This shader evaluates a single orbital on a 3D grid in parallel.
/// Each workgroup processes a slice of the grid.
pub const ORBITAL_GRID_SHADER: &str = r#"
struct BasisFunc {
    center: vec3<f32>,
    lx: u32,
    ly: u32,
    lz: u32,
    n_primitives: u32,
    coefficient: f32,
    // Primitives follow in a separate buffer
};

struct GridParams {
    origin: vec3<f32>,
    spacing: f32,
    dims: vec3<u32>,
    orbital_index: u32,
};

@group(0) @binding(0) var<storage, read> basis: array<BasisFunc>;
@group(0) @binding(1) var<storage, read> mo_coeffs: array<f32>;
@group(0) @binding(2) var<storage, read> primitives: array<vec2<f32>>; // (alpha, coeff)
@group(0) @binding(3) var<uniform> params: GridParams;
@group(0) @binding(4) var<storage, read_write> output: array<f32>;

@compute @workgroup_size(8, 8, 4)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let ix = gid.x;
    let iy = gid.y;
    let iz = gid.z;

    if (ix >= params.dims.x || iy >= params.dims.y || iz >= params.dims.z) {
        return;
    }

    let r = params.origin + vec3<f32>(f32(ix), f32(iy), f32(iz)) * params.spacing;
    let flat_idx = ix * params.dims.y * params.dims.z + iy * params.dims.z + iz;

    var psi: f32 = 0.0;
    let n_basis = arrayLength(&basis);

    for (var mu: u32 = 0u; mu < n_basis; mu = mu + 1u) {
        let c_mu = mo_coeffs[mu * params.dims.x + params.orbital_index]; // simplified indexing
        if (abs(c_mu) < 1e-7) {
            continue;
        }

        let bf = basis[mu];
        let dr = r - bf.center;
        let r2 = dot(dr, dr);

        // Angular part
        var angular: f32 = 1.0;
        for (var i: u32 = 0u; i < bf.lx; i = i + 1u) { angular *= dr.x; }
        for (var i: u32 = 0u; i < bf.ly; i = i + 1u) { angular *= dr.y; }
        for (var i: u32 = 0u; i < bf.lz; i = i + 1u) { angular *= dr.z; }

        // Radial part (contracted)
        var radial: f32 = 0.0;
        for (var p: u32 = 0u; p < bf.n_primitives; p = p + 1u) {
            let prim = primitives[mu * 3u + p]; // max 3 primitives for STO-3G
            radial += prim.y * exp(-prim.x * r2);
        }

        psi += c_mu * bf.coefficient * angular * radial;
    }

    output[flat_idx] = psi;
}
"#;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_grid_params_basic() {
        let positions = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let params = GridParams::from_molecule(&positions, 0.5, 3.0);

        assert!(params.dimensions[0] > 0);
        assert!(params.dimensions[1] > 0);
        assert!(params.dimensions[2] > 0);
        assert!(params.n_points() > 0);
        assert!(params.origin[0] < -2.0); // padding should extend beyond molecule
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
        assert!((p[2] - 0.0).abs() < 1e-12);
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
}
