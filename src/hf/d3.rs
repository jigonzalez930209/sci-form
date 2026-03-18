//! Grimme DFT-D3 dispersion correction with Becke-Johnson damping.
//!
//! Adds long-range van der Waals interactions missing from HF:
//! $$E_{disp} = -\sum_{n=6,8} s_n \sum_{A>B} \frac{C_n^{AB}}{R_{AB}^n + f_{damp}^n}$$
//!
//! Reference: Grimme, S. et al. J. Chem. Phys. 132 (2010): 154104.

use serde::{Deserialize, Serialize};

/// D3 correction result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct D3Result {
    /// Total D3 dispersion energy (Hartree).
    pub energy: f64,
    /// C6 contribution.
    pub e6: f64,
    /// C8 contribution.
    pub e8: f64,
}

/// D3 parameters for HF-3c (BJ damping).
const S6: f64 = 1.0;
const S8: f64 = 0.8777;
const A1: f64 = 0.4171;
const A2: f64 = 2.9149;

/// Compute D3-BJ dispersion correction.
pub fn compute_d3_energy(elements: &[u8], positions_bohr: &[[f64; 3]]) -> D3Result {
    let n = elements.len();
    let mut e6 = 0.0;
    let mut e8 = 0.0;

    for a in 0..n {
        for b in (a + 1)..n {
            let r = distance(positions_bohr, a, b);
            if r < 1e-10 {
                continue;
            }

            let c6 = get_c6(elements[a], elements[b]);
            let c8 = c6_to_c8(c6, elements[a], elements[b]);

            let r0 = (c8 / c6).sqrt();
            let f6 = bj_damping(r, r0, 6);
            let f8 = bj_damping(r, r0, 8);

            e6 -= S6 * c6 / (r.powi(6) + f6);
            e8 -= S8 * c8 / (r.powi(8) + f8);
        }
    }

    D3Result {
        energy: e6 + e8,
        e6,
        e8,
    }
}

/// Becke-Johnson damping function.
fn bj_damping(_r: f64, r0: f64, n: i32) -> f64 {
    let r_cut = A1 * r0 + A2;
    r_cut.powi(n)
}

/// Empirical C6 coefficients (Hartree·Bohr⁶).
fn get_c6(z1: u8, z2: u8) -> f64 {
    let c6_1 = atomic_c6(z1);
    let c6_2 = atomic_c6(z2);
    // Geometric combination rule
    (c6_1 * c6_2).sqrt()
}

/// Approximate atomic C6 values.
fn atomic_c6(z: u8) -> f64 {
    match z {
        1 => 6.50,
        6 => 46.60,
        7 => 24.20,
        8 => 15.60,
        9 => 9.52,
        15 => 185.0,
        16 => 134.0,
        17 => 94.60,
        _ => 30.0, // rough fallback
    }
}

/// Estimate C8 from C6 using Casimir-Polder relation.
fn c6_to_c8(c6: f64, z1: u8, z2: u8) -> f64 {
    let q1 = atomic_q(z1);
    let q2 = atomic_q(z2);
    3.0 * c6 * (q1 * q2).sqrt()
}

/// Atomic quadrupole-like parameter for C6→C8.
fn atomic_q(z: u8) -> f64 {
    match z {
        1 => 2.23,
        6 => 7.40,
        7 => 6.07,
        8 => 5.14,
        9 => 4.09,
        _ => 5.0,
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
    use super::*;
    use super::super::basis::ANG_TO_BOHR;

    #[test]
    fn test_d3_h2() {
        let positions = [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.74 * ANG_TO_BOHR],
        ];
        let result = compute_d3_energy(&[1, 1], &positions);
        // D3 energy should be small and negative
        assert!(result.energy < 0.0, "D3 should be attractive");
        assert!(result.energy.abs() < 0.01, "D3 for H2 should be tiny");
    }

    #[test]
    fn test_d3_ethane_larger() {
        // Two carbons + H atoms → more D3
        let r = 1.54 * ANG_TO_BOHR;
        let positions = [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, r],
        ];
        let result = compute_d3_energy(&[6, 6], &positions);
        assert!(result.energy < 0.0);
    }
}
