//! ReaxFF non-bonded interactions: van der Waals + shielded Coulomb.

use super::taper::taper_function;

/// Compute non-bonded energy (van der Waals + Coulomb with shielding).
///
/// Returns (vdw_energy, coulomb_energy) in kcal/mol.
pub fn compute_nonbonded_energy(
    positions_flat: &[f64],
    charges: &[f64],
    elements: &[u8],
    cutoff: f64,
) -> (f64, f64) {
    let n = elements.len();
    let mut e_vdw = 0.0;
    let mut e_coul = 0.0;

    for i in 0..n {
        let xi = positions_flat[3 * i];
        let yi = positions_flat[3 * i + 1];
        let zi = positions_flat[3 * i + 2];

        for j in (i + 1)..n {
            let dx = xi - positions_flat[3 * j];
            let dy = yi - positions_flat[3 * j + 1];
            let dz = zi - positions_flat[3 * j + 2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();

            if r > cutoff || r < 0.01 {
                continue;
            }

            let tap = taper_function(r, cutoff);

            // van der Waals (Morse-type with shielding)
            let (eps, r0, gamma_w) = vdw_params(elements[i], elements[j]);
            let r_shielded = (r.powi(7) + gamma_w.powi(7)).powf(1.0 / 7.0);
            let rho = r_shielded / r0;
            let _e_morse = eps * (rho.powf(-2.0 * 13.0) - 2.0 * rho.powf(-13.0));
            // Simplified: use standard LJ form
            let sigma = r0 / 2.0f64.powf(1.0 / 6.0);
            let sr6 = (sigma / r_shielded).powi(6);
            let e_lj = 4.0 * eps * (sr6 * sr6 - sr6);
            e_vdw += tap * e_lj;

            // Shielded Coulomb
            let gamma_q: f64 = 0.5; // Å (shielding parameter)
            let r_shield_coul = (r.powi(3) + gamma_q.powi(3)).powf(1.0 / 3.0);
            // 332.0638 converts e²/Å to kcal/mol
            let e_coulomb = 332.0638 * charges[i] * charges[j] / r_shield_coul;
            e_coul += tap * e_coulomb;
        }
    }

    (e_vdw, e_coul)
}

/// Get van der Waals parameters for a pair of element types.
/// Returns (epsilon kcal/mol, r0 Å, gamma shielding Å).
fn vdw_params(z1: u8, z2: u8) -> (f64, f64, f64) {
    let (e1, r1) = element_vdw(z1);
    let (e2, r2) = element_vdw(z2);
    let eps = (e1 * e2).sqrt();
    let r0 = (r1 + r2) / 2.0;
    let gamma = 1.0; // default shielding
    (eps, r0, gamma)
}

/// Single-element vdW parameters (epsilon kcal/mol, r0 Å).
fn element_vdw(z: u8) -> (f64, f64) {
    match z {
        1 => (0.0200, 3.195),  // H
        6 => (0.0951, 3.851),  // C
        7 => (0.0774, 3.660),  // N
        8 => (0.0957, 3.500),  // O
        9 => (0.0725, 3.364),  // F
        16 => (0.2500, 4.035), // S
        17 => (0.2270, 3.947), // Cl
        _ => (0.0500, 3.500),  // generic
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nonbonded_energy_h2_is_finite() {
        let positions = [0.0, 0.0, 0.0, 0.74, 0.0, 0.0];
        let charges = [0.0, 0.0];
        let elements = [1u8, 1];
        let (e_vdw, e_coul) = compute_nonbonded_energy(&positions, &charges, &elements, 10.0);
        assert!(e_vdw.is_finite(), "vdW energy must be finite");
        assert!(e_coul.is_finite(), "Coulomb energy must be finite");
    }

    #[test]
    fn coulomb_zero_for_neutral_atoms() {
        let positions = [0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        let charges = [0.0, 0.0];
        let elements = [1u8, 1];
        let (_e_vdw, e_coul) = compute_nonbonded_energy(&positions, &charges, &elements, 10.0);
        assert!(
            e_coul.abs() < 1e-10,
            "Coulomb should be zero for uncharged atoms"
        );
    }

    #[test]
    fn vdw_repulsive_at_short_range() {
        // Very close — should be repulsive (positive)
        let positions = [0.0, 0.0, 0.0, 0.5, 0.0, 0.0];
        let charges = [0.0, 0.0];
        let elements = [1u8, 1];
        let (e_vdw, _) = compute_nonbonded_energy(&positions, &charges, &elements, 10.0);
        assert!(e_vdw > 0.0, "vdW should be repulsive at very short range");
    }
}
