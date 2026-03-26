//! D4 Parameter Infrastructure — Core
//!
//! Atomic reference data, coordination numbers, and dynamic C6/C8.

/// Per-element D4 parameters.
#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
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
pub fn get_d4_params(z: u8) -> D4Params {
    match z {
        1 => D4Params {
            c6_ref: 6.50,
            c8_ref: 94.6,
            alpha0: 4.50,
            r_cov: 0.60,
        },
        5 => D4Params {
            c6_ref: 99.5,
            c8_ref: 3430.0,
            alpha0: 21.0,
            r_cov: 1.61,
        },
        6 => D4Params {
            c6_ref: 46.6,
            c8_ref: 1315.0,
            alpha0: 12.0,
            r_cov: 1.46,
        },
        7 => D4Params {
            c6_ref: 24.2,
            c8_ref: 560.0,
            alpha0: 7.40,
            r_cov: 1.42,
        },
        8 => D4Params {
            c6_ref: 15.6,
            c8_ref: 300.0,
            alpha0: 5.40,
            r_cov: 1.38,
        },
        9 => D4Params {
            c6_ref: 9.52,
            c8_ref: 160.0,
            alpha0: 3.80,
            r_cov: 1.34,
        },
        14 => D4Params {
            c6_ref: 305.0,
            c8_ref: 15700.0,
            alpha0: 37.0,
            r_cov: 2.21,
        },
        15 => D4Params {
            c6_ref: 185.0,
            c8_ref: 8200.0,
            alpha0: 25.0,
            r_cov: 2.08,
        },
        16 => D4Params {
            c6_ref: 134.0,
            c8_ref: 5470.0,
            alpha0: 19.6,
            r_cov: 1.96,
        },
        17 => D4Params {
            c6_ref: 94.6,
            c8_ref: 3430.0,
            alpha0: 15.0,
            r_cov: 1.88,
        },
        35 => D4Params {
            c6_ref: 162.0,
            c8_ref: 6800.0,
            alpha0: 21.0,
            r_cov: 2.16,
        },
        53 => D4Params {
            c6_ref: 385.0,
            c8_ref: 20500.0,
            alpha0: 35.0,
            r_cov: 2.51,
        },
        22 => D4Params {
            c6_ref: 379.0,
            c8_ref: 17100.0,
            alpha0: 45.0,
            r_cov: 2.61,
        },
        26 => D4Params {
            c6_ref: 202.0,
            c8_ref: 7900.0,
            alpha0: 30.0,
            r_cov: 2.45,
        },
        29 => D4Params {
            c6_ref: 171.0,
            c8_ref: 6100.0,
            alpha0: 24.0,
            r_cov: 2.49,
        },
        30 => D4Params {
            c6_ref: 207.0,
            c8_ref: 8000.0,
            alpha0: 28.0,
            r_cov: 2.45,
        },
        // --- Additional main-group elements ---
        2 => D4Params {
            c6_ref: 1.46,
            c8_ref: 14.1,
            alpha0: 1.38,
            r_cov: 0.46,
        }, // He
        3 => D4Params {
            c6_ref: 1387.0,
            c8_ref: 127800.0,
            alpha0: 164.0,
            r_cov: 2.49,
        }, // Li
        4 => D4Params {
            c6_ref: 214.0,
            c8_ref: 12100.0,
            alpha0: 38.0,
            r_cov: 1.89,
        }, // Be
        10 => D4Params {
            c6_ref: 6.38,
            c8_ref: 76.0,
            alpha0: 2.67,
            r_cov: 1.10,
        }, // Ne
        11 => D4Params {
            c6_ref: 1556.0,
            c8_ref: 165500.0,
            alpha0: 163.0,
            r_cov: 2.99,
        }, // Na
        12 => D4Params {
            c6_ref: 626.0,
            c8_ref: 46400.0,
            alpha0: 71.0,
            r_cov: 2.74,
        }, // Mg
        13 => D4Params {
            c6_ref: 528.0,
            c8_ref: 35300.0,
            alpha0: 60.0,
            r_cov: 2.38,
        }, // Al
        18 => D4Params {
            c6_ref: 64.3,
            c8_ref: 2450.0,
            alpha0: 11.1,
            r_cov: 1.74,
        }, // Ar
        19 => D4Params {
            c6_ref: 3897.0,
            c8_ref: 536300.0,
            alpha0: 290.0,
            r_cov: 3.59,
        }, // K
        20 => D4Params {
            c6_ref: 2190.0,
            c8_ref: 246200.0,
            alpha0: 160.0,
            r_cov: 3.31,
        }, // Ca
        // --- 3d transition metals (remaining) ---
        21 => D4Params {
            c6_ref: 571.0,
            c8_ref: 30500.0,
            alpha0: 55.0,
            r_cov: 2.76,
        }, // Sc
        23 => D4Params {
            c6_ref: 335.0,
            c8_ref: 14800.0,
            alpha0: 40.0,
            r_cov: 2.53,
        }, // V
        24 => D4Params {
            c6_ref: 276.0,
            c8_ref: 11400.0,
            alpha0: 35.0,
            r_cov: 2.49,
        }, // Cr
        25 => D4Params {
            c6_ref: 245.0,
            c8_ref: 9700.0,
            alpha0: 33.0,
            r_cov: 2.49,
        }, // Mn
        27 => D4Params {
            c6_ref: 175.0,
            c8_ref: 6500.0,
            alpha0: 27.0,
            r_cov: 2.38,
        }, // Co
        28 => D4Params {
            c6_ref: 155.0,
            c8_ref: 5500.0,
            alpha0: 24.0,
            r_cov: 2.34,
        }, // Ni
        // --- 4d transition metals ---
        39 => D4Params {
            c6_ref: 700.0,
            c8_ref: 42000.0,
            alpha0: 65.0,
            r_cov: 2.95,
        }, // Y
        40 => D4Params {
            c6_ref: 540.0,
            c8_ref: 29000.0,
            alpha0: 55.0,
            r_cov: 2.80,
        }, // Zr
        42 => D4Params {
            c6_ref: 340.0,
            c8_ref: 15000.0,
            alpha0: 40.0,
            r_cov: 2.61,
        }, // Mo
        44 => D4Params {
            c6_ref: 252.0,
            c8_ref: 10200.0,
            alpha0: 33.0,
            r_cov: 2.53,
        }, // Ru
        45 => D4Params {
            c6_ref: 220.0,
            c8_ref: 8500.0,
            alpha0: 29.0,
            r_cov: 2.49,
        }, // Rh
        46 => D4Params {
            c6_ref: 205.0,
            c8_ref: 7600.0,
            alpha0: 26.0,
            r_cov: 2.49,
        }, // Pd
        47 => D4Params {
            c6_ref: 253.0,
            c8_ref: 10200.0,
            alpha0: 33.0,
            r_cov: 2.72,
        }, // Ag
        48 => D4Params {
            c6_ref: 323.0,
            c8_ref: 14500.0,
            alpha0: 46.0,
            r_cov: 2.76,
        }, // Cd
        // --- 5d transition metals ---
        72 => D4Params {
            c6_ref: 485.0,
            c8_ref: 24000.0,
            alpha0: 50.0,
            r_cov: 2.80,
        }, // Hf
        74 => D4Params {
            c6_ref: 375.0,
            c8_ref: 16500.0,
            alpha0: 42.0,
            r_cov: 2.68,
        }, // W
        76 => D4Params {
            c6_ref: 280.0,
            c8_ref: 11500.0,
            alpha0: 34.0,
            r_cov: 2.53,
        }, // Os
        77 => D4Params {
            c6_ref: 245.0,
            c8_ref: 9500.0,
            alpha0: 30.0,
            r_cov: 2.53,
        }, // Ir
        78 => D4Params {
            c6_ref: 225.0,
            c8_ref: 8500.0,
            alpha0: 28.0,
            r_cov: 2.53,
        }, // Pt
        79 => D4Params {
            c6_ref: 255.0,
            c8_ref: 10500.0,
            alpha0: 36.0,
            r_cov: 2.57,
        }, // Au
        80 => D4Params {
            c6_ref: 305.0,
            c8_ref: 13400.0,
            alpha0: 34.0,
            r_cov: 2.53,
        }, // Hg
        // --- Period 4-5 main group ---
        31 => D4Params {
            c6_ref: 498.0,
            c8_ref: 30600.0,
            alpha0: 50.0,
            r_cov: 2.34,
        }, // Ga
        32 => D4Params {
            c6_ref: 354.0,
            c8_ref: 18600.0,
            alpha0: 40.0,
            r_cov: 2.30,
        }, // Ge
        33 => D4Params {
            c6_ref: 246.0,
            c8_ref: 11400.0,
            alpha0: 30.0,
            r_cov: 2.23,
        }, // As
        34 => D4Params {
            c6_ref: 210.0,
            c8_ref: 9100.0,
            alpha0: 26.0,
            r_cov: 2.19,
        }, // Se
        49 => D4Params {
            c6_ref: 700.0,
            c8_ref: 47000.0,
            alpha0: 65.0,
            r_cov: 2.53,
        }, // In
        50 => D4Params {
            c6_ref: 530.0,
            c8_ref: 32000.0,
            alpha0: 53.0,
            r_cov: 2.53,
        }, // Sn
        51 => D4Params {
            c6_ref: 405.0,
            c8_ref: 21000.0,
            alpha0: 43.0,
            r_cov: 2.49,
        }, // Sb
        52 => D4Params {
            c6_ref: 345.0,
            c8_ref: 16500.0,
            alpha0: 38.0,
            r_cov: 2.49,
        }, // Te
        _ => D4Params {
            c6_ref: 50.0,
            c8_ref: 1500.0,
            alpha0: 12.0,
            r_cov: 1.80,
        },
    }
}

/// Compute D4 fractional coordination number.
pub fn d4_coordination_number(elements: &[u8], positions: &[[f64; 3]]) -> Vec<f64> {
    let n = elements.len();
    let bohr_to_ang = 0.529177;
    let mut cn = vec![0.0; n];

    for i in 0..n {
        let pi = get_d4_params(elements[i]);
        for j in (i + 1)..n {
            let pj = get_d4_params(elements[j]);
            let r_cov_ij = (pi.r_cov + pj.r_cov) * bohr_to_ang;

            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r_ij = (dx * dx + dy * dy + dz * dz).sqrt();

            if r_ij < 1e-10 {
                continue;
            }

            let f = 1.0 / (1.0 + (-16.0 * (r_cov_ij / r_ij - 1.0)).exp());
            cn[i] += f;
            cn[j] += f;
        }
    }

    cn
}

/// Get C6 reference for a pair using the London formula.
pub fn get_c6_reference(z_a: u8, z_b: u8) -> f64 {
    let pa = get_d4_params(z_a);
    let pb = get_d4_params(z_b);
    let sum_alpha = pa.alpha0 + pb.alpha0;
    if sum_alpha < 1e-15 {
        return 0.0;
    }
    1.5 * pa.alpha0 * pb.alpha0 / sum_alpha
}

/// Get dynamic C6 adjusted by coordination number.
pub fn dynamic_c6(z_a: u8, z_b: u8, cn_a: f64, cn_b: f64) -> f64 {
    let c6_ref = get_c6_reference(z_a, z_b);
    let cn_ref_a = expected_cn(z_a);
    let cn_ref_b = expected_cn(z_b);

    let w_a = (-4.0 * (cn_a / cn_ref_a - 1.0).powi(2)).exp();
    let w_b = (-4.0 * (cn_b / cn_ref_b - 1.0).powi(2)).exp();

    c6_ref * w_a * w_b
}

fn expected_cn(z: u8) -> f64 {
    match z {
        1 => 1.0,
        2 => 0.0,                // He
        3 | 11 | 19 => 1.0,      // Li, Na, K
        4 | 12 | 20 => 2.0,      // Be, Mg, Ca
        5 | 13 => 3.0,           // B, Al
        6 | 14 | 32 => 4.0,      // C, Si, Ge
        7 | 15 | 33 => 3.0,      // N, P, As
        8 | 16 | 34 | 52 => 2.0, // O, S, Se, Te
        9 | 17 | 35 | 53 => 1.0, // F, Cl, Br, I
        10 | 18 => 0.0,          // Ne, Ar
        // 3d TMs
        21 => 6.0, // Sc
        22 => 6.0, // Ti
        23 => 6.0, // V
        24 => 6.0, // Cr
        25 => 6.0, // Mn
        26 => 6.0, // Fe
        27 => 6.0, // Co
        28 => 4.0, // Ni
        29 => 4.0, // Cu
        30 => 4.0, // Zn
        // 4d TMs
        39 => 6.0,
        40 => 6.0,
        42 => 6.0,
        44 => 6.0,
        45 => 6.0,
        46 => 4.0,
        47 => 4.0,
        48 => 4.0,
        // 5d TMs
        72 => 6.0,
        74 => 6.0,
        76 => 6.0,
        77 => 6.0,
        78 => 4.0,
        79 => 4.0,
        80 => 4.0,
        // Main group 4-5
        31 | 49 => 3.0, // Ga, In
        50 | 51 => 4.0, // Sn, Sb
        _ => 4.0,
    }
}

/// Compute C8 from C6 using the recursion relation.
pub fn c8_from_c6(c6: f64, z_a: u8, z_b: u8) -> f64 {
    let pa = get_d4_params(z_a);
    let pb = get_d4_params(z_b);
    let qa = if pa.c6_ref > 1e-10 {
        (pa.c8_ref / pa.c6_ref).sqrt()
    } else {
        5.0
    };
    let qb = if pb.c6_ref > 1e-10 {
        (pb.c8_ref / pb.c6_ref).sqrt()
    } else {
        5.0
    };
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
            [0.63, 0.63, 0.63],
            [-0.63, -0.63, 0.63],
            [-0.63, 0.63, -0.63],
            [0.63, -0.63, -0.63],
        ];
        let cn = d4_coordination_number(&elements, &pos);
        assert!(cn[0] > 1.0 && cn[0] < 5.0, "C CN = {}", cn[0]);
        assert!(cn[1] > 0.3 && cn[1] < 1.5, "H CN = {}", cn[1]);
    }

    #[test]
    fn test_c6_reference_positive() {
        for &z in &[1u8, 6, 7, 8, 16] {
            let c6 = get_c6_reference(z, z);
            assert!(c6 > 0.0, "C6({},{}) = {}", z, z, c6);
        }
    }

    #[test]
    fn test_c8_from_c6_positive() {
        let c6 = get_c6_reference(6, 6);
        let c8 = c8_from_c6(c6, 6, 6);
        assert!(c8 > c6, "C8 should be larger than C6");
    }
}
