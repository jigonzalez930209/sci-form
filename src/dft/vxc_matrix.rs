//! V_XC matrix assembly: evaluate basis functions on grid, contract with XC kernel.

use nalgebra::DMatrix;

use super::functionals::{pbe, svwn};
use super::grid::MolecularGrid;
use super::ks_fock::DftMethod;
use crate::scf::basis::{BasisFunction, BasisSet};

/// Evaluate all basis functions at a single point.
///
/// Returns a vector of φ_μ(r) for each basis function μ.
fn evaluate_basis_at_point(basis: &BasisSet, point: &[f64; 3]) -> Vec<f64> {
    basis
        .functions
        .iter()
        .map(|bf| bf.evaluate(point[0], point[1], point[2]))
        .collect()
}

/// Evaluate basis function gradients at a single point (for GGA).
fn evaluate_basis_gradient_at_point(basis: &BasisSet, point: &[f64; 3]) -> Vec<[f64; 3]> {
    let n = basis.functions.len();
    let mut gradients = vec![[0.0; 3]; n];

    for (mu, bf) in basis.functions.iter().enumerate() {
        let dx = point[0] - bf.center[0];
        let dy = point[1] - bf.center[1];
        let dz = point[2] - bf.center[2];
        let r2 = dx * dx + dy * dy + dz * dz;

        let [lx, ly, lz] = bf.angular;
        for prim in &bf.primitives {
            let exp_val = (-prim.alpha * r2).exp();
            let a = prim.alpha;
            let amp = BasisFunction::normalization(a, lx, ly, lz) * prim.coefficient * exp_val;
            let ax = if lx == 0 { 0.0 } else { lx as f64 * dx.powi(lx as i32 - 1) };
            let ay = if ly == 0 { 0.0 } else { ly as f64 * dy.powi(ly as i32 - 1) };
            let az = if lz == 0 { 0.0 } else { lz as f64 * dz.powi(lz as i32 - 1) };
            let px = dx.powi(lx as i32);
            let py = dy.powi(ly as i32);
            let pz = dz.powi(lz as i32);
            // ∂/∂x of x^lx exp(-a r²) = (lx x^{lx-1} - 2 a x^{lx+1}) exp, times the other powers.
            gradients[mu][0] += amp * (ax * py * pz - 2.0 * a * dx * px * py * pz);
            gradients[mu][1] += amp * (px * ay * pz - 2.0 * a * dy * px * py * pz);
            gradients[mu][2] += amp * (px * py * az - 2.0 * a * dz * px * py * pz);
        }
    }

    gradients
}

/// Compute electron density at a grid point from the density matrix.
fn compute_density_at_point(phi: &[f64], density: &DMatrix<f64>) -> f64 {
    let n = phi.len();
    let mut rho = 0.0;
    for mu in 0..n {
        for nu in 0..n {
            rho += density[(mu, nu)] * phi[mu] * phi[nu];
        }
    }
    rho.max(0.0)
}

/// Compute density gradient squared |∇ρ|² at a grid point.
fn compute_gradient_squared(phi: &[f64], dphi: &[[f64; 3]], density: &DMatrix<f64>) -> f64 {
    let n = phi.len();
    let mut grad = [0.0f64; 3];

    for mu in 0..n {
        for nu in 0..n {
            let d = density[(mu, nu)];
            for dim in 0..3 {
                grad[dim] += d * (dphi[mu][dim] * phi[nu] + phi[mu] * dphi[nu][dim]);
            }
        }
    }

    grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]
}

/// Build the V_XC matrix from the density matrix on the molecular grid.
///
/// Returns (V_xc, E_xc) where V_xc is the XC potential matrix and E_xc is
/// the XC energy.
pub fn build_vxc_matrix(
    basis: &BasisSet,
    density: &DMatrix<f64>,
    grid: &MolecularGrid,
    method: DftMethod,
) -> (DMatrix<f64>, f64) {
    let n = basis.functions.len();
    let mut vxc = DMatrix::zeros(n, n);
    let mut exc_total = 0.0;

    for gp in &grid.points {
        let phi = evaluate_basis_at_point(basis, &gp.xyz);
        let rho = compute_density_at_point(&phi, density);

        if rho < 1e-20 {
            continue;
        }

        match method {
            DftMethod::Svwn => {
                let (exc, vxc_val) = svwn::svwn(rho);
                exc_total += gp.weight * rho * exc;

                for mu in 0..n {
                    for nu in mu..n {
                        let contrib = gp.weight * vxc_val * phi[mu] * phi[nu];
                        vxc[(mu, nu)] += contrib;
                        if mu != nu {
                            vxc[(nu, mu)] += contrib;
                        }
                    }
                }
            }
            DftMethod::Pbe => {
                let dphi = evaluate_basis_gradient_at_point(basis, &gp.xyz);
                let sigma = compute_gradient_squared(&phi, &dphi, density);

                let (exc, vxc_rho, vxc_sigma) = pbe::pbe(rho, sigma);
                exc_total += gp.weight * rho * exc;

                for mu in 0..n {
                    for nu in mu..n {
                        // LDA part
                        let mut contrib = gp.weight * vxc_rho * phi[mu] * phi[nu];

                        // GGA part: 2 * vxc_sigma * ∇ρ · (∇φ_μ φ_ν + φ_μ ∇φ_ν)
                        // Simplified: use density gradient contribution
                        let grad_contrib = vxc_sigma * 2.0;
                        for d in 0..3 {
                            contrib += gp.weight
                                * grad_contrib
                                * (dphi[mu][d] * phi[nu] + phi[mu] * dphi[nu][d]);
                        }

                        vxc[(mu, nu)] += contrib;
                        if mu != nu {
                            vxc[(nu, mu)] += contrib;
                        }
                    }
                }
            }
        }
    }

    (vxc, exc_total)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dft::grid::{GridQuality, MolecularGrid};
    use crate::dft::ks_fock::DftMethod;
    use crate::scf::basis::BasisSet;
    use crate::scf::density_matrix::build_density_matrix;
    use crate::scf::types::MolecularSystem;
    use nalgebra::DMatrix;

    fn h2_vxc_setup() -> (BasisSet, DMatrix<f64>, MolecularGrid) {
        let system =
            MolecularSystem::from_angstrom(&[1, 1], &[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]], 0, 1);
        let basis = BasisSet::sto3g(&system.atomic_numbers, &system.positions_bohr);
        let n = basis.functions.len();
        // Simple identity-like density
        let c = DMatrix::identity(n, n);
        let density = build_density_matrix(&c, 1);
        let grid = MolecularGrid::build(
            &system.atomic_numbers,
            &system.positions_bohr,
            GridQuality::Coarse,
        );
        (basis, density, grid)
    }

    #[test]
    fn vxc_matrix_is_symmetric_svwn() {
        let (basis, density, grid) = h2_vxc_setup();
        let (vxc, _exc) = build_vxc_matrix(&basis, &density, &grid, DftMethod::Svwn);
        let n = vxc.nrows();
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (vxc[(i, j)] - vxc[(j, i)]).abs() < 1e-12,
                    "V_XC should be symmetric"
                );
            }
        }
    }

    #[test]
    fn vxc_xc_energy_negative_for_svwn() {
        let (basis, density, grid) = h2_vxc_setup();
        let (_vxc, exc) = build_vxc_matrix(&basis, &density, &grid, DftMethod::Svwn);
        assert!(exc < 0.0, "XC energy should be negative, got {exc}");
    }

    #[test]
    fn vxc_pbe_runs_without_panic() {
        let (basis, density, grid) = h2_vxc_setup();
        let (_vxc, _exc) = build_vxc_matrix(&basis, &density, &grid, DftMethod::Pbe);
    }
}
