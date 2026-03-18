//! Short-range basis (SRB) correction for basis set incompleteness.
//!
//! Corrects systematic errors in bond lengths when using minimal basis sets
//! by adding a short-range pair potential:
//! $$E_{SRB} = s_{SRB} \sum_{A>B} \sqrt{Z_A Z_B}
//!   \exp\left(-\gamma_{SRB} R_{AB}^2\right)
//!   [1 - f_{damp}(R_{AB})]$$
//!
//! Reference: Sure, R.; Grimme, S. J. Comput. Chem. 34 (2013): 1672.

/// Compute SRB correction energy (Hartree).
pub fn compute_srb(elements: &[u8], positions_bohr: &[[f64; 3]]) -> f64 {
    let n = elements.len();
    let mut e_srb = 0.0;

    // SRB parameters for HF-3c
    let s_srb = 0.03;
    let gamma_srb = 0.7;

    for a in 0..n {
        for b in (a + 1)..n {
            let r = distance(positions_bohr, a, b);
            if r < 1e-10 {
                continue;
            }

            let za = elements[a] as f64;
            let zb = elements[b] as f64;
            let r_cov = covalent_radius_bohr(elements[a])
                + covalent_radius_bohr(elements[b]);

            let gaussian = (-gamma_srb * (r - r_cov).powi(2)).exp();
            let prefactor = s_srb * (za * zb).sqrt();

            e_srb += prefactor * gaussian;
        }
    }

    e_srb
}

/// Covalent radii in Bohr.
fn covalent_radius_bohr(z: u8) -> f64 {
    let ang = match z {
        1 => 0.31,
        6 => 0.76,
        7 => 0.71,
        8 => 0.66,
        9 => 0.57,
        14 => 1.11,
        15 => 1.07,
        16 => 1.05,
        17 => 1.02,
        _ => 1.0,
    };
    ang * super::basis::ANG_TO_BOHR
}

fn distance(positions: &[[f64; 3]], a: usize, b: usize) -> f64 {
    let dx = positions[a][0] - positions[b][0];
    let dy = positions[a][1] - positions[b][1];
    let dz = positions[a][2] - positions[b][2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::basis::ANG_TO_BOHR;

    #[test]
    fn test_srb_at_equilibrium() {
        // C-C bond at ~1.54 Å
        let r = 1.54 * ANG_TO_BOHR;
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, r]];
        let e = compute_srb(&[6, 6], &positions);
        assert!(e > 0.0, "SRB should be positive");
        assert!(e < 0.5, "SRB should be small: {e}");
    }

    #[test]
    fn test_srb_decays_away_from_covalent() {
        let r1 = 1.54 * ANG_TO_BOHR;
        let r2 = 5.0 * ANG_TO_BOHR;
        let e1 = compute_srb(&[6, 6], &[[0.0, 0.0, 0.0], [0.0, 0.0, r1]]);
        let e2 = compute_srb(&[6, 6], &[[0.0, 0.0, 0.0], [0.0, 0.0, r2]]);
        assert!(e1 > e2, "SRB should be larger near covalent radius");
    }
}
