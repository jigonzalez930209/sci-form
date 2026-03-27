//! Electronegativity Equalization Method (EEM) for ReaxFF charge equilibration.
//!
//! Solves for atomic charges that equalize the electrochemical potential
//! across all atoms, subject to total charge neutrality.

use nalgebra::{DMatrix, DVector};

/// EEM parameters for an atom type.
#[derive(Debug, Clone)]
pub struct EemAtomParams {
    /// Electronegativity (eV).
    pub chi: f64,
    /// Chemical hardness (eV).
    pub eta: f64,
    /// Shielding parameter (Å).
    pub gamma: f64,
}

/// Solve EEM equations to equilibrate charges.
///
/// Minimizes E_charge = Σ_i (χ_i q_i + ½ η_i q_i²) + Σ_{i<j} q_i q_j / r_ij(shielded)
/// subject to Σ_i q_i = Q_total.
///
/// Uses Lagrange multiplier method → (n+1) × (n+1) linear system.
pub fn solve_eem(
    positions_flat: &[f64],
    atom_params: &[EemAtomParams],
    total_charge: f64,
) -> Result<Vec<f64>, String> {
    let n = atom_params.len();
    if n == 0 {
        return Ok(vec![]);
    }

    // Build the (n+1) × (n+1) system: [J λ; 1 0] [q; μ] = [-χ; Q]
    let dim = n + 1;
    let mut a = DMatrix::zeros(dim, dim);
    let mut b = DVector::zeros(dim);

    for i in 0..n {
        // Diagonal: η_i (hardness)
        a[(i, i)] = atom_params[i].eta;

        // Off-diagonal: shielded Coulomb
        for j in (i + 1)..n {
            let dx = positions_flat[3 * i] - positions_flat[3 * j];
            let dy = positions_flat[3 * i + 1] - positions_flat[3 * j + 1];
            let dz = positions_flat[3 * i + 2] - positions_flat[3 * j + 2];
            let r2 = dx * dx + dy * dy + dz * dz;

            // Shielded Coulomb: 1 / (r³ + γ³)^{1/3}
            let gamma_ij = (atom_params[i].gamma + atom_params[j].gamma) / 2.0;
            let shielded = 1.0 / (r2.powf(1.5) + gamma_ij.powi(3)).powf(1.0 / 3.0);

            // Convert to eV·Å units (Coulomb constant in eV·Å)
            let coulomb = 14.4 * shielded; // e²/(4πε₀) ≈ 14.4 eV·Å

            a[(i, j)] = coulomb;
            a[(j, i)] = coulomb;
        }

        // Lagrange multiplier column/row
        a[(i, n)] = 1.0;
        a[(n, i)] = 1.0;

        // RHS: -χ_i
        b[i] = -atom_params[i].chi;
    }

    // Charge neutrality constraint
    b[n] = total_charge;

    // Solve via LU decomposition
    let solution = a.lu().solve(&b).ok_or("EEM linear system is singular")?;

    let charges: Vec<f64> = solution.rows(0, n).iter().copied().collect();
    Ok(charges)
}

/// Default EEM parameters for common elements.
pub fn default_eem_params(element: u8) -> EemAtomParams {
    match element {
        1 => EemAtomParams {
            chi: 3.7248,
            eta: 9.6093,
            gamma: 0.8203,
        }, // H
        6 => EemAtomParams {
            chi: 5.9666,
            eta: 7.0000,
            gamma: 0.9000,
        }, // C
        7 => EemAtomParams {
            chi: 6.8418,
            eta: 8.0555,
            gamma: 0.7820,
        }, // N
        8 => EemAtomParams {
            chi: 8.5000,
            eta: 8.3122,
            gamma: 1.0898,
        }, // O
        _ => EemAtomParams {
            chi: 5.0,
            eta: 7.0,
            gamma: 0.9,
        }, // generic
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eem_h2_charges_sum_to_zero() {
        let positions = [0.0, 0.0, 0.0, 0.74, 0.0, 0.0];
        let params = vec![default_eem_params(1), default_eem_params(1)];
        let charges = solve_eem(&positions, &params, 0.0).unwrap();
        let sum: f64 = charges.iter().sum();
        assert!(
            sum.abs() < 1e-10,
            "charges should sum to total charge (0), got {sum}"
        );
    }

    #[test]
    fn eem_symmetric_molecule_equal_charges() {
        // H-H symmetric — charges should be equal
        let positions = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let params = vec![default_eem_params(1), default_eem_params(1)];
        let charges = solve_eem(&positions, &params, 0.0).unwrap();
        assert!(
            (charges[0] - charges[1]).abs() < 1e-10,
            "symmetric molecule should have equal charges"
        );
    }

    #[test]
    fn eem_nonzero_total_charge() {
        let positions = [0.0, 0.0, 0.0, 1.5, 0.0, 0.0];
        let params = vec![default_eem_params(6), default_eem_params(8)];
        let charges = solve_eem(&positions, &params, 1.0).unwrap();
        let sum: f64 = charges.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-8,
            "charges should sum to total charge (+1), got {sum}"
        );
    }

    #[test]
    fn default_eem_params_known_elements() {
        let h = default_eem_params(1);
        let c = default_eem_params(6);
        assert!(h.eta > 0.0);
        assert!(
            c.chi != h.chi,
            "C and H should have different electronegativities"
        );
    }
}
