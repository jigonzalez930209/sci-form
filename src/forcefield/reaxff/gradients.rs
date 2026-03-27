//! Analytical gradients for ReaxFF energy terms.
//!
//! Uses numerical differentiation (central differences) for robustness
//! during alpha development. Will be replaced with analytical gradients
//! once the energy terms are validated.

use super::bond_order::compute_bond_orders;
use super::eem::{default_eem_params, solve_eem};
use super::energy::compute_bonded_energy;
use super::nonbonded::compute_nonbonded_energy;
use super::params::ReaxffParams;

/// Compute ReaxFF energy and gradient (numerical, central differences).
///
/// Returns (total_energy_kcal_mol, gradient_flat) where gradient is
/// ∂E/∂r for each coordinate [gx0,gy0,gz0,...].
pub fn compute_reaxff_gradient(
    positions_flat: &[f64],
    elements: &[u8],
    params: &ReaxffParams,
) -> Result<(f64, Vec<f64>), String> {
    let n3 = positions_flat.len();
    let delta = 1e-5; // Å

    // Map elements to parameter indices
    let atom_params: Vec<_> = elements
        .iter()
        .map(|&z| {
            params
                .element_index(z)
                .map(|idx| params.atom_params[idx].clone())
                .unwrap_or_else(|| params.atom_params.last().cloned().unwrap())
        })
        .collect();

    // Reference energy
    let e0 = total_reaxff_energy(positions_flat, elements, &atom_params, params)?;

    // Central-difference gradient — each coordinate displacement is independent.
    // With the `parallel` feature each displacement pair runs on its own thread.
    #[cfg(feature = "parallel")]
    let grad: Vec<f64> = {
        use rayon::prelude::*;
        let result: Result<Vec<f64>, String> = (0..n3)
            .into_par_iter()
            .map(|i| {
                let mut pp = positions_flat.to_vec();
                let mut pm = positions_flat.to_vec();
                pp[i] += delta;
                pm[i] -= delta;
                let ep = total_reaxff_energy(&pp, elements, &atom_params, params)?;
                let em = total_reaxff_energy(&pm, elements, &atom_params, params)?;
                Ok((ep - em) / (2.0 * delta))
            })
            .collect();
        result?
    };

    #[cfg(not(feature = "parallel"))]
    let grad: Vec<f64> = {
        let mut g = vec![0.0; n3];
        let mut pos_plus = positions_flat.to_vec();
        let mut pos_minus = positions_flat.to_vec();
        for i in 0..n3 {
            pos_plus[i] = positions_flat[i] + delta;
            pos_minus[i] = positions_flat[i] - delta;
            let ep = total_reaxff_energy(&pos_plus, elements, &atom_params, params)?;
            let em = total_reaxff_energy(&pos_minus, elements, &atom_params, params)?;
            g[i] = (ep - em) / (2.0 * delta);
            pos_plus[i] = positions_flat[i];
            pos_minus[i] = positions_flat[i];
        }
        g
    };

    Ok((e0, grad))
}

/// Compute total ReaxFF energy (bonded + non-bonded + EEM charges).
fn total_reaxff_energy(
    positions_flat: &[f64],
    elements: &[u8],
    atom_params: &[super::params::ReaxffAtomParams],
    params: &ReaxffParams,
) -> Result<f64, String> {
    // Compute bond orders using per-atom params (sized to match atom count)
    let bo = compute_bond_orders(positions_flat, atom_params, params.cutoff);

    // Build a local params copy whose atom_params matches the actual atom count
    // so that compute_bonded_energy iterates over the correct range.
    let n = atom_params.len();
    let local_params = ReaxffParams {
        atom_params: atom_params.to_vec(),
        bond_de: (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        // Try to reuse global bond_de via element indices when possible
                        let gi = params.element_index(atom_params[i].element);
                        let gj = params.element_index(atom_params[j].element);
                        match (gi, gj) {
                            (Some(a), Some(b)) => params.get_bond_de(a, b),
                            _ => 100.0,
                        }
                    })
                    .collect()
            })
            .collect(),
        ..params.clone()
    };

    // Compute bonded energy with consistent sizing
    let bonded = compute_bonded_energy(positions_flat, &bo, &local_params);

    // Compute EEM charges
    let eem_params: Vec<_> = elements.iter().map(|&z| default_eem_params(z)).collect();
    let charges = solve_eem(positions_flat, &eem_params, 0.0)?;

    // Compute non-bonded energy
    let (e_vdw, e_coul) =
        compute_nonbonded_energy(positions_flat, &charges, elements, params.cutoff);

    Ok(bonded.total + e_vdw + e_coul)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::reaxff::params::ReaxffParams;

    #[test]
    fn gradient_h2_is_finite() {
        let positions = [0.0, 0.0, 0.0, 0.74, 0.0, 0.0];
        let elements = [1u8, 1];
        let params = ReaxffParams::default_chon();
        let (energy, grad) = compute_reaxff_gradient(&positions, &elements, &params).unwrap();
        assert!(energy.is_finite(), "energy must be finite");
        assert_eq!(grad.len(), 6);
        for g in &grad {
            assert!(g.is_finite(), "gradient component must be finite");
        }
    }

    #[test]
    fn gradient_symmetric_h2() {
        // By symmetry, forces on H atoms should be equal and opposite
        let positions = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let elements = [1u8, 1];
        let params = ReaxffParams::default_chon();
        let (_, grad) = compute_reaxff_gradient(&positions, &elements, &params).unwrap();
        // grad[0..3] = force on atom 0, grad[3..6] = force on atom 1
        // By Newton's third law: f0 ≈ -f1
        for d in 0..3 {
            assert!(
                (grad[d] + grad[3 + d]).abs() < 0.1,
                "forces should be approximately opposite in dim {d}: {} vs {}",
                grad[d],
                grad[3 + d]
            );
        }
    }
}
