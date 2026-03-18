//! Geometric counterpoise (gCP) correction for basis set superposition error.
//!
//! Corrects the artificial lowering of energy when small basis sets borrow
//! functions from neighboring atoms. Uses an empirical atom-pair potential.
//!
//! $$E_{gCP} = \sigma \sum_{A>B} \sum_a^{A} e_a^{miss}
//!   \exp\left(-\alpha (R_{AB})^\beta\right) S_{AB}^\gamma N_B^{virt}$$
//!
//! Reference: Kruse, H.; Grimme, S. J. Chem. Phys. 136 (2012): 154101.

/// Compute gCP correction energy (Hartree).
pub fn compute_gcp(elements: &[u8], positions_bohr: &[[f64; 3]]) -> f64 {
    let n = elements.len();
    let mut e_gcp = 0.0;

    // gCP parameters for HF-3c/MINIX
    let sigma = 0.1825;
    let alpha = 1.0265;
    let beta = 1.2148;

    for a in 0..n {
        for b in (a + 1)..n {
            let r = distance(positions_bohr, a, b);
            if r < 1e-10 {
                continue;
            }

            let e_miss_a = emiss(elements[a]);
            let e_miss_b = emiss(elements[b]);
            let nvirt_a = n_virtual(elements[a]);
            let nvirt_b = n_virtual(elements[b]);

            let decay = (-alpha * r.powf(beta)).exp();

            // A borrows from B
            e_gcp += sigma * e_miss_a * decay * nvirt_b as f64;
            // B borrows from A
            e_gcp += sigma * e_miss_b * decay * nvirt_a as f64;
        }
    }

    e_gcp
}

/// Missing basis energy per atom (element-specific empirical values).
fn emiss(z: u8) -> f64 {
    match z {
        1 => 0.0030,
        6 => 0.0450,
        7 => 0.0380,
        8 => 0.0320,
        9 => 0.0280,
        15 => 0.0600,
        16 => 0.0550,
        17 => 0.0480,
        _ => 0.0400,
    }
}

/// Number of virtual orbitals available (proxy for basis borrowing capacity).
fn n_virtual(z: u8) -> usize {
    match z {
        1 => 1,
        6..=9 => 4,
        14..=17 => 9,
        _ => 4,
    }
}

fn distance(positions: &[[f64; 3]], a: usize, b: usize) -> f64 {
    let dx = positions[a][0] - positions[b][0];
    let dy = positions[a][1] - positions[b][1];
    let dz = positions[a][2] - positions[b][2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[cfg(test)]
mod tests {
    use super::super::basis::ANG_TO_BOHR;
    use super::*;

    #[test]
    fn test_gcp_positive() {
        let r = 1.54 * ANG_TO_BOHR;
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, r]];
        let e = compute_gcp(&[6, 6], &positions);
        assert!(e > 0.0, "gCP correction should be positive (repulsive)");
    }

    #[test]
    fn test_gcp_decays_with_distance() {
        let r1 = 1.5 * ANG_TO_BOHR;
        let r2 = 5.0 * ANG_TO_BOHR;
        let e1 = compute_gcp(&[6, 6], &[[0.0, 0.0, 0.0], [0.0, 0.0, r1]]);
        let e2 = compute_gcp(&[6, 6], &[[0.0, 0.0, 0.0], [0.0, 0.0, r2]]);
        assert!(e1 > e2, "gCP should decay with distance");
    }
}
