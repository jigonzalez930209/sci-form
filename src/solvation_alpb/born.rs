//! Born Radii and GB Kernel — Core
//!
//! Hawkins-Cramer-Truhlar (HCT) Born radii with Onufriev-Bashford-Case
//! corrections, and the generalized Born f_GB kernel.

/// Per-atom intrinsic radii (Å) by atomic number (Bondi radii).
pub(crate) fn intrinsic_radius(z: u8) -> f64 {
    match z {
        1 => 1.20,
        6 => 1.70,
        7 => 1.55,
        8 => 1.52,
        9 => 1.47,
        14 => 2.10,
        15 => 1.80,
        16 => 1.80,
        17 => 1.75,
        35 => 1.85,
        53 => 1.98,
        _ => 1.70,
    }
}

/// Result of Born radii calculation.
#[derive(Debug, Clone)]
pub struct AlpbBornRadii {
    /// Effective Born radii per atom (Å).
    pub radii: Vec<f64>,
    /// Intrinsic radii used (Å).
    pub intrinsic: Vec<f64>,
}

/// Compute effective Born radii using the HCT model with OBC correction.
pub fn compute_born_radii(
    elements: &[u8],
    positions: &[[f64; 3]],
    probe_radius: f64,
) -> AlpbBornRadii {
    let n = elements.len();
    let rho: Vec<f64> = elements
        .iter()
        .map(|&z| intrinsic_radius(z) + probe_radius * 0.1)
        .collect();

    let alpha = 1.0;
    let beta = 0.8;
    let gamma = 4.85;

    let mut born_radii = vec![0.0; n];

    for i in 0..n {
        let mut psi = 0.0;
        for j in 0..n {
            if i == j {
                continue;
            }
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r_ij = (dx * dx + dy * dy + dz * dz).sqrt();
            let rho_j = rho[j];

            if r_ij > rho_j {
                let l_ij = rho[i].max(r_ij - rho_j);
                let u_ij = r_ij + rho_j;
                if u_ij > l_ij {
                    psi += 0.5
                        * (1.0 / l_ij - 1.0 / u_ij
                            + 0.25 * (1.0 / u_ij - 1.0 / l_ij) * (r_ij * r_ij - rho_j * rho_j)
                            + 0.5 * (1.0 / u_ij.powi(2) - 1.0 / l_ij.powi(2)) * r_ij);
                }
            }
        }

        let psi_scaled = psi * rho[i];
        let tanh_val =
            (alpha * psi_scaled - beta * psi_scaled.powi(2) + gamma * psi_scaled.powi(3)).tanh();

        let inv_r_eff = 1.0 / rho[i] - tanh_val / rho[i];
        born_radii[i] = if inv_r_eff > 1e-10 {
            1.0 / inv_r_eff
        } else {
            100.0
        };
    }

    AlpbBornRadii {
        radii: born_radii,
        intrinsic: elements.iter().map(|&z| intrinsic_radius(z)).collect(),
    }
}

/// Generalized Born kernel f_GB.
///
/// f_GB = sqrt(r_ij² + R_i R_j exp(-r_ij²/(4 R_i R_j)))
pub fn gb_kernel(r_ij: f64, r_i: f64, r_j: f64) -> f64 {
    let ri_rj = r_i * r_j;
    if ri_rj < 1e-15 {
        return r_ij;
    }
    (r_ij * r_ij + ri_rj * (-r_ij * r_ij / (4.0 * ri_rj)).exp()).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_born_radii_positive() {
        let elements = [8, 1, 1];
        let pos = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let result = compute_born_radii(&elements, &pos, 1.4);
        for &r in &result.radii {
            assert!(r > 0.0 && r < 100.0, "Born radius = {}", r);
        }
    }

    #[test]
    fn test_gb_kernel_asymptotic() {
        let f = gb_kernel(10.0, 1.5, 1.5);
        assert!((f - 10.0).abs() < 0.1, "f_GB at r=10: {}", f);
    }

    #[test]
    fn test_gb_kernel_self() {
        let f = gb_kernel(0.0, 1.5, 1.5);
        assert!((f - 1.5).abs() < 0.01, "f_GB self: {}", f);
    }
}
