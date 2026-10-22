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
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct AlpbBornRadii {
    /// Effective Born radii per atom (Å).
    pub radii: Vec<f64>,
    /// Intrinsic radii used (Å).
    pub intrinsic: Vec<f64>,
}

/// Compute effective Born radii using the HCT model with OBC corrections.
///
/// The pairwise descreening integral is the analytical ("ACE") expression of
/// Onufriev, Bashford and Case (J. Phys. Chem. B 2004, 108, 2648) as used
/// by OpenMM's OBC-I generalised Born model:
///
/// ```text
/// or_i  = radius_i - 0.009
/// I_ij  = 0.5 * ( 1/L - 1/U
///                + 0.25*(1/U^2 - 1/L^2)*(r - sr2^2/r)
///                + 0.5*ln(L/U)/r + C )
/// U     = r + sr2            sr2 = or_j
/// L     = max(or_i, |r - sr2|)
/// C     = 2*(1/or_i - 1/L)   if sr2 >= r + or_i else 0
/// psi_i = sum_j I_ij
/// 1/R_i = 1/or_i - tanh(a*psi - b*psi^2 + c*psi^3)/radius_i   with psi = psi_i*or_i
/// ```
///
/// The intrinsic radii `radius_i` include the probe offset (`+0.1*probe`),
/// as in the rest of the ALPB module.
///
/// # Errors
///
/// Returns `Err` when the geometry is invalid (mismatched lengths, empty
/// system, non-finite coordinates) or when the probe radius is negative or
/// non-finite.
pub fn compute_born_radii(
    elements: &[u8],
    positions: &[[f64; 3]],
    probe_radius: f64,
) -> Result<AlpbBornRadii, String> {
    crate::validation::validate_geometry(elements, positions)?;
    crate::validation::validate_probe_radius(probe_radius)?;
    let n = elements.len();
    let radius: Vec<f64> = elements
        .iter()
        .map(|&z| intrinsic_radius(z) + probe_radius * 0.1)
        .collect();
    // OBC offset ("alpha" of the original OBC formulation).
    const OFFSET: f64 = 0.009;
    let orad: Vec<f64> = radius.iter().map(|&r| (r - OFFSET).max(0.1)).collect();

    let alpha = 1.0;
    let beta = 0.8;
    let gamma = 4.85;

    // Adaptive Born radius cap: largest intrinsic radius + 2× max pairwise
    // distance, with a floor of 50 Å. This prevents non-physical clamping
    // for large molecular systems.
    let mut max_dist: f64 = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            if d > max_dist {
                max_dist = d;
            }
        }
    }
    let born_cap = (max_dist + 10.0).max(50.0);

    let mut born_radii = vec![0.0; n];

    for i in 0..n {
        let or_i = orad[i];
        let mut psi = 0.0f64;

        for j in 0..n {
            if i == j {
                continue;
            }
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < 1e-10 {
                continue;
            }
            let sr2 = orad[j];
            if r + sr2 < or_i {
                // The pair does not reach the descreening sphere of atom i.
                continue;
            }
            let u = r + sr2;
            let d = (r - sr2).abs();
            let l = or_i.max(d);
            let inv_l = 1.0 / l;
            let inv_u = 1.0 / u;
            // Atom j completely inside the descreening sphere of atom i.
            let c = if sr2 >= r + or_i {
                2.0 * (1.0 / or_i - inv_l)
            } else {
                0.0
            };
            let integral = 0.5
                * (inv_l - inv_u
                    + 0.25 * (inv_u * inv_u - inv_l * inv_l) * (r - sr2 * sr2 / r)
                    + 0.5 * (l / u).ln() / r
                    + c);
            psi += integral;
        }

        let psi_scaled = psi * or_i;
        let tanh_val =
            (alpha * psi_scaled - beta * psi_scaled.powi(2) + gamma * psi_scaled.powi(3)).tanh();

        let inv_r_eff = 1.0 / or_i - tanh_val / radius[i];
        born_radii[i] = if inv_r_eff > 1e-10 {
            (1.0 / inv_r_eff).min(born_cap)
        } else {
            born_cap
        };
    }

    Ok(AlpbBornRadii {
        radii: born_radii,
        intrinsic: elements.iter().map(|&z| intrinsic_radius(z)).collect(),
    })
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
        let result = compute_born_radii(&elements, &pos, 1.4).expect("valid input");
        for &r in &result.radii {
            assert!(r > 0.0 && r < 100.0, "Born radius = {}", r);
        }
    }

    #[test]
    fn c_m1_born_radii_reject_a_short_position_list() {
        let elements = [8, 1, 1];
        let pos = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0]];
        let err = compute_born_radii(&elements, &pos, 1.4)
            .expect_err("a short position list must be rejected, not panic");
        assert!(err.contains("positions has 2 entries"), "{}", err);
    }

    #[test]
    fn c_h4_born_radii_reject_a_negative_probe_radius() {
        let elements = [8, 1, 1];
        let pos = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let err = compute_born_radii(&elements, &pos, -1.0)
            .expect_err("a negative probe radius must be rejected");
        assert!(err.contains("probe_radius"), "{}", err);
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
