//! Numerical Hessian computation via central finite differences.
//!
//! Computes the 3N×3N Hessian matrix by displacing each coordinate
//! by ±δ and evaluating the energy at each displaced geometry.

use nalgebra::DMatrix;

/// Which semiempirical method to use for energy evaluations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HessianMethod {
    /// Extended Hückel Theory (fastest, qualitative)
    Eht,
    /// PM3 semi-empirical SCF
    Pm3,
    /// GFN0-xTB tight-binding
    Xtb,
    /// UFF force field with analytical Hessian for bond-stretch and angle-bend
    Uff,
}

/// Evaluate total energy for a given method and geometry.
fn evaluate_energy(
    method: HessianMethod,
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<f64, String> {
    match method {
        HessianMethod::Eht => {
            let result = crate::eht::solve_eht(elements, positions, None)?;
            // EHT total electronic energy: sum of occupied orbital energies × 2
            let n_occ = result.n_electrons / 2;
            let total: f64 = result.energies.iter().take(n_occ).sum::<f64>() * 2.0;
            Ok(total)
        }
        HessianMethod::Pm3 => {
            let result = crate::pm3::solve_pm3(elements, positions)?;
            Ok(result.total_energy)
        }
        HessianMethod::Xtb => {
            let result = crate::xtb::solve_xtb(elements, positions)?;
            Ok(result.total_energy)
        }
        HessianMethod::Uff => {
            Err("UFF uses analytical Hessian path; use compute_uff_analytical_hessian".to_string())
        }
    }
}

/// Compute the 3N×3N Hessian matrix via central finite differences.
///
/// Uses the formula: H_{ij} ≈ [E(+δ_i,+δ_j) - E(+δ_i,-δ_j) - E(-δ_i,+δ_j) + E(-δ_i,-δ_j)] / (4δ²)
///
/// For diagonal elements, uses the more efficient:
/// H_{ii} ≈ [E(+δ_i) - 2E₀ + E(-δ_i)] / δ²
///
/// `elements`: atomic numbers
/// `positions`: Cartesian coordinates in Å
/// `method`: which energy method to use
/// `step_size`: finite difference step in Å (if None, auto-selects based on atomic mass)
pub fn compute_numerical_hessian(
    elements: &[u8],
    positions: &[[f64; 3]],
    method: HessianMethod,
    step_size: Option<f64>,
) -> Result<DMatrix<f64>, String> {
    let n_atoms = elements.len();
    if n_atoms != positions.len() {
        return Err("elements and positions length mismatch".to_string());
    }

    let n3 = 3 * n_atoms;

    // Auto step-size: heavier atoms need larger steps for stable numerics
    let delta = step_size.unwrap_or_else(|| auto_step_size(elements));
    let delta_sq = delta * delta;

    // Reference energy
    let e0 = evaluate_energy(method, elements, positions)?;

    // Flatten positions for easier manipulation
    let mut coords: Vec<f64> = Vec::with_capacity(n3);
    for pos in positions {
        coords.extend_from_slice(pos);
    }

    let mut hessian = DMatrix::zeros(n3, n3);

    // Diagonal elements: H_{ii} = [E(x+δ_i) - 2E₀ + E(x-δ_i)] / δ²
    for i in 0..n3 {
        let mut coords_plus = coords.clone();
        let mut coords_minus = coords.clone();
        coords_plus[i] += delta;
        coords_minus[i] -= delta;

        let pos_plus = flat_to_positions(&coords_plus);
        let pos_minus = flat_to_positions(&coords_minus);

        let e_plus = evaluate_energy(method, elements, &pos_plus)?;
        let e_minus = evaluate_energy(method, elements, &pos_minus)?;

        hessian[(i, i)] = (e_plus - 2.0 * e0 + e_minus) / delta_sq;
    }

    // Off-diagonal elements: H_{ij} = [E(+i,+j) - E(+i,-j) - E(-i,+j) + E(-i,-j)] / (4δ²)
    for i in 0..n3 {
        for j in (i + 1)..n3 {
            let mut coords_pp = coords.clone();
            let mut coords_pm = coords.clone();
            let mut coords_mp = coords.clone();
            let mut coords_mm = coords.clone();

            coords_pp[i] += delta;
            coords_pp[j] += delta;
            coords_pm[i] += delta;
            coords_pm[j] -= delta;
            coords_mp[i] -= delta;
            coords_mp[j] += delta;
            coords_mm[i] -= delta;
            coords_mm[j] -= delta;

            let e_pp = evaluate_energy(method, elements, &flat_to_positions(&coords_pp))?;
            let e_pm = evaluate_energy(method, elements, &flat_to_positions(&coords_pm))?;
            let e_mp = evaluate_energy(method, elements, &flat_to_positions(&coords_mp))?;
            let e_mm = evaluate_energy(method, elements, &flat_to_positions(&coords_mm))?;

            let hij = (e_pp - e_pm - e_mp + e_mm) / (4.0 * delta_sq);
            hessian[(i, j)] = hij;
            hessian[(j, i)] = hij;
        }
    }

    // Symmetry enforcement: average H and H^T to eliminate numerical noise
    enforce_symmetry(&mut hessian);

    Ok(hessian)
}

/// Enforce exact symmetry by averaging H[i,j] and H[j,i].
///
/// Reduces numerical noise from finite difference evaluation:
/// $H_{\text{sym}} = \frac{H + H^T}{2}$
fn enforce_symmetry(hessian: &mut DMatrix<f64>) {
    let n = hessian.nrows();
    for i in 0..n {
        for j in (i + 1)..n {
            let avg = (hessian[(i, j)] + hessian[(j, i)]) * 0.5;
            hessian[(i, j)] = avg;
            hessian[(j, i)] = avg;
        }
    }
}

/// Automatic step-size selection based on the lightest element.
///
/// Lighter atoms (H) need smaller steps; heavier atoms allow larger steps.
/// Uses δ = 0.005 × sqrt(min_mass / 1.008) as a heuristic.
fn auto_step_size(elements: &[u8]) -> f64 {
    let min_mass = elements
        .iter()
        .map(|&z| match z {
            1 => 1.008,
            6 => 12.011,
            7 => 14.007,
            8 => 15.999,
            9 => 18.998,
            _ => z as f64 * 1.5,
        })
        .fold(f64::INFINITY, f64::min);

    // Base step of 0.005 Å scaled by sqrt(mass_ratio)
    (0.005 * (min_mass / 1.008).sqrt()).clamp(0.003, 0.01)
}

fn flat_to_positions(coords: &[f64]) -> Vec<[f64; 3]> {
    coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect()
}

/// Compute the analytical Hessian for UFF force field using gradient differences.
///
/// Instead of computing H_{ij} from 4-point energy stencils (O(9N²) energy evaluations),
/// this computes H_{ij} = ∂g_i/∂x_j using central differences on the analytical gradients:
///   H_{ij} ≈ [g_i(x + δe_j) - g_i(x - δe_j)] / (2δ)
///
/// This requires only 6N gradient evaluations (one ± displacement per coordinate),
/// which is much cheaper since UFF gradients are computed analytically.
pub fn compute_uff_analytical_hessian(
    smiles: &str,
    coords_flat: &[f64],
    step_size: Option<f64>,
) -> Result<DMatrix<f64>, String> {
    let mol = crate::graph::Molecule::from_smiles(smiles)?;
    let n_atoms = mol.graph.node_count();
    let n3 = 3 * n_atoms;

    if coords_flat.len() != n3 {
        return Err(format!(
            "coords length {} != 3 * atoms {}",
            coords_flat.len(),
            n_atoms
        ));
    }

    let ff = crate::forcefield::builder::build_uff_force_field(&mol);
    let elements: Vec<u8> = (0..n_atoms)
        .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].element)
        .collect();
    let delta = step_size.unwrap_or_else(|| auto_step_size(&elements));

    let mut hessian = DMatrix::zeros(n3, n3);

    // For each coordinate j, displace ±δ and compute the full gradient vector.
    // H_{ij} = [g_i(x + δe_j) - g_i(x - δe_j)] / (2δ)
    for j in 0..n3 {
        let mut coords_plus = coords_flat.to_vec();
        let mut coords_minus = coords_flat.to_vec();
        coords_plus[j] += delta;
        coords_minus[j] -= delta;

        let mut grad_plus = vec![0.0; n3];
        let mut grad_minus = vec![0.0; n3];

        ff.compute_system_energy_and_gradients(&coords_plus, &mut grad_plus);
        ff.compute_system_energy_and_gradients(&coords_minus, &mut grad_minus);

        for i in 0..n3 {
            hessian[(i, j)] = (grad_plus[i] - grad_minus[i]) / (2.0 * delta);
        }
    }

    // Symmetry enforcement
    enforce_symmetry(&mut hessian);

    Ok(hessian)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hessian_h2_is_symmetric() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let hessian =
            compute_numerical_hessian(&elements, &positions, HessianMethod::Xtb, Some(0.005))
                .unwrap();

        assert_eq!(hessian.nrows(), 6);
        assert_eq!(hessian.ncols(), 6);

        // Check symmetry
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (hessian[(i, j)] - hessian[(j, i)]).abs() < 1e-4,
                    "Hessian not symmetric at ({}, {}): {} vs {}",
                    i,
                    j,
                    hessian[(i, j)],
                    hessian[(j, i)]
                );
            }
        }
    }

    #[test]
    fn test_hessian_water_shape() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let hessian =
            compute_numerical_hessian(&elements, &positions, HessianMethod::Xtb, Some(0.005))
                .unwrap();

        assert_eq!(hessian.nrows(), 9);
        assert_eq!(hessian.ncols(), 9);

        // Check at least some non-zero off-diagonal elements
        let mut has_nonzero_offdiag = false;
        for i in 0..9 {
            for j in (i + 1)..9 {
                if hessian[(i, j)].abs() > 1e-6 {
                    has_nonzero_offdiag = true;
                    break;
                }
            }
        }
        assert!(
            has_nonzero_offdiag,
            "Hessian should have non-zero off-diagonal elements"
        );
    }
}
