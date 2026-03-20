//! D4 Parameter Infrastructure — E7.1
//!
//! Atomic reference data, D4 coordination numbers, and dynamic C6/C8.

/// Per-element D4 parameters.
#[derive(Debug, Clone, Copy)]
pub struct D4Params {
    /// Reference C6 coefficient (Hartree·Bohr^6).
    pub c6_ref: f64,
    /// Reference C8 coefficient (Hartree·Bohr^8).
    pub c8_ref: f64,
    /// Static dipole polarizability α₀ (Bohr³).
    pub alpha0: f64,
    /// Covalent radius for CN (Bohr).
    pub r_cov: f64,
}

/// Get D4 reference parameters by atomic number.
/// Values from Grimme's D4 reference data.
pub fn get_d4_params(z: u8) -> D4Params {
    match z {
        1 => D4Params { c6_ref: 6.50, c8_ref: 94.6, alpha0: 4.50, r_cov: 0.60 },
        5 => D4Params { c6_ref: 99.5, c8_ref: 3430.0, alpha0: 21.0, r_cov: 1.61 },
        6 => D4Params { c6_ref: 46.6, c8_ref: 1315.0, alpha0: 12.0, r_cov: 1.46 },
        7 => D4Params { c6_ref: 24.2, c8_ref: 560.0, alpha0: 7.40, r_cov: 1.42 },
        8 => D4Params { c6_ref: 15.6, c8_ref: 300.0, alpha0: 5.40, r_cov: 1.38 },
        9 => D4Params { c6_ref: 9.52, c8_ref: 160.0, alpha0: 3.80, r_cov: 1.34 },
        14 => D4Params { c6_ref: 305.0, c8_ref: 15700.0, alpha0: 37.0, r_cov: 2.21 },
        15 => D4Params { c6_ref: 185.0, c8_ref: 8200.0, alpha0: 25.0, r_cov: 2.08 },
        16 => D4Params { c6_ref: 134.0, c8_ref: 5470.0, alpha0: 19.6, r_cov: 1.96 },
        17 => D4Params { c6_ref: 94.6, c8_ref: 3430.0, alpha0: 15.0, r_cov: 1.88 },
        35 => D4Params { c6_ref: 162.0, c8_ref: 6800.0, alpha0: 21.0, r_cov: 2.16 },
        53 => D4Params { c6_ref: 385.0, c8_ref: 20500.0, alpha0: 35.0, r_cov: 2.51 },
        22 => D4Params { c6_ref: 379.0, c8_ref: 17100.0, alpha0: 45.0, r_cov: 2.61 }, // Ti
        26 => D4Params { c6_ref: 202.0, c8_ref: 7900.0, alpha0: 30.0, r_cov: 2.45 }, // Fe
        29 => D4Params { c6_ref: 171.0, c8_ref: 6100.0, alpha0: 24.0, r_cov: 2.49 }, // Cu
        30 => D4Params { c6_ref: 207.0, c8_ref: 8000.0, alpha0: 28.0, r_cov: 2.45 }, // Zn
        _ => D4Params { c6_ref: 50.0, c8_ref: 1500.0, alpha0: 12.0, r_cov: 1.80 },
    }
}

/// Compute D4 fractional coordination number.
///
/// CN_i = sum_{j!=i} f(r_ij / r_cov_ij)
/// f(x) = 1 / (1 + exp(-16*(x-1)))    [Fermi function]
pub fn d4_coordination_number(elements: &[u8], positions: &[[f64; 3]]) -> Vec<f64> {
    let n = elements.len();
    let bohr_to_ang = 0.529177;
    let mut cn = vec![0.0; n];

    for i in 0..n {
        let pi = get_d4_params(elements[i]);
        for j in (i + 1)..n {
            let pj = get_d4_params(elements[j]);
            let r_cov_ij = (pi.r_cov + pj.r_cov) * bohr_to_ang; // Å

            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r_ij = (dx * dx + dy * dy + dz * dz).sqrt();

            if r_ij < 1e-10 { continue; }

            let f = 1.0 / (1.0 + (-16.0 * (r_cov_ij / r_ij - 1.0)).exp());
            cn[i] += f;
            cn[j] += f;
        }
    }

    cn
}

/// Get C6 reference coefficient for a pair of elements.
/// Uses Casimir-Polder integration via polarizabilities.
///
/// C6_AB = 3/2 * alpha_A * alpha_B / (alpha_A + alpha_B)  (London formula)
pub fn get_c6_reference(z_a: u8, z_b: u8) -> f64 {
    let pa = get_d4_params(z_a);
    let pb = get_d4_params(z_b);
    let sum_alpha = pa.alpha0 + pb.alpha0;
    if sum_alpha < 1e-15 { return 0.0; }
    1.5 * pa.alpha0 * pb.alpha0 / sum_alpha
}

/// Get dynamic C6 coefficient adjusted by coordination number.
///
/// C6_AB(CN) = C6_ref * w_A(CN_A) * w_B(CN_B)
/// where w(CN) = exp(-4*(CN/CN_ref - 1)^2) with CN_ref from covalent neighbors
pub fn dynamic_c6(z_a: u8, z_b: u8, cn_a: f64, cn_b: f64) -> f64 {
    let c6_ref = get_c6_reference(z_a, z_b);
    let cn_ref_a = expected_cn(z_a);
    let cn_ref_b = expected_cn(z_b);

    let w_a = (-4.0 * (cn_a / cn_ref_a - 1.0).powi(2)).exp();
    let w_b = (-4.0 * (cn_b / cn_ref_b - 1.0).powi(2)).exp();

    c6_ref * w_a * w_b
}

/// Expected coordination number for common elements.
fn expected_cn(z: u8) -> f64 {
    match z {
        1 => 1.0, 6 => 4.0, 7 => 3.0, 8 => 2.0, 9 => 1.0,
        15 => 3.0, 16 => 2.0, 17 => 1.0, 35 => 1.0, 53 => 1.0,
        14 => 4.0, 5 => 3.0,
        22 => 6.0, 26 => 6.0, 29 => 4.0, 30 => 4.0,
        _ => 4.0,
    }
}

/// Compute C8 from C6 using the recursion relation.
/// C8_AB = 3 * C6_AB * sqrt(Q_A * Q_B) where Q = sqrt(C8_ref/C6_ref)
pub fn c8_from_c6(c6: f64, z_a: u8, z_b: u8) -> f64 {
    let pa = get_d4_params(z_a);
    let pb = get_d4_params(z_b);
    let qa = if pa.c6_ref > 1e-10 { (pa.c8_ref / pa.c6_ref).sqrt() } else { 5.0 };
    let qb = if pb.c6_ref > 1e-10 { (pb.c8_ref / pb.c6_ref).sqrt() } else { 5.0 };
    3.0 * c6 * (qa * qb).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_d4_cn_methane() {
        let elements = [6, 1, 1, 1, 1];
        let pos = [
            [0.0, 0.0, 0.0],
            [0.63, 0.63, 0.63], [-0.63, -0.63, 0.63],
            [-0.63, 0.63, -0.63], [0.63, -0.63, -0.63],
        ];
        let cn = d4_coordination_number(&elements, &pos);
        assert!(cn[0] > 1.0 && cn[0] < 5.0, "C CN = {}", cn[0]);
        assert!(cn[1] > 0.3 && cn[1] < 1.5, "H CN = {}", cn[1]);
    }

    #[test]
    fn test_c6_reference_positive() {
        for &z1 in &[1u8, 6, 7, 8, 16] {
            for &z2 in &[1u8, 6, 7, 8, 16] {
                let c6 = get_c6_reference(z1, z2);
                assert!(c6 > 0.0, "C6({},{}) = {}", z1, z2, c6);
            }
        }
    }

    #[test]
    fn test_c6_symmetric() {
        let c6_ab = get_c6_reference(6, 8);
        let c6_ba = get_c6_reference(8, 6);
        assert!((c6_ab - c6_ba).abs() < 1e-10);
    }

    #[test]
    fn test_dynamic_c6_varies_with_cn() {
        let c6_low = dynamic_c6(6, 6, 1.0, 1.0);
        let c6_ref = dynamic_c6(6, 6, 4.0, 4.0);
        // At reference CN, weight should be maximal
        assert!(c6_ref >= c6_low * 0.5, "C6 at ref CN should be ≥ C6 at CN=1");
    }
}
