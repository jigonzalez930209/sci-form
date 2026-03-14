//! Gasteiger-Marsili iterative electronegativity equalization.
//!
//! Computes partial atomic charges by iteratively transferring charge along
//! bonds according to orbital electronegativity differences.
//!
//! Reference: J. Gasteiger & M. Marsili, *Tetrahedron* **36**, 3219–3228, 1980.

use serde::{Deserialize, Serialize};

/// Gasteiger parametrization for one element: a, b, c coefficients for the
/// orbital-electronegativity polynomial  χ(q) = a + b·q + c·q².
#[derive(Debug, Clone, Copy)]
pub struct GasteigerParams {
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

impl GasteigerParams {
    /// Orbital electronegativity at partial charge q.
    #[inline]
    pub fn chi(&self, q: f64) -> f64 {
        self.a + self.b * q + self.c * q * q
    }
}

/// Look up Gasteiger a/b/c coefficients for an element by atomic number.
///
/// Returns `None` for unsupported elements. Parameters sourced from
/// Gasteiger & Marsili (1980) and OpenBabel implementation.
pub fn get_gasteiger_params(z: u8) -> Option<GasteigerParams> {
    match z {
        1 => Some(GasteigerParams {
            a: 7.17,
            b: 6.24,
            c: -0.56,
        }),
        6 => Some(GasteigerParams {
            a: 7.98,
            b: 9.18,
            c: 1.88,
        }),
        7 => Some(GasteigerParams {
            a: 11.54,
            b: 10.82,
            c: 1.36,
        }),
        8 => Some(GasteigerParams {
            a: 14.18,
            b: 12.92,
            c: 1.39,
        }),
        9 => Some(GasteigerParams {
            a: 14.66,
            b: 13.85,
            c: 2.31,
        }),
        15 => Some(GasteigerParams {
            a: 8.90,
            b: 8.24,
            c: 0.96,
        }),
        16 => Some(GasteigerParams {
            a: 10.14,
            b: 9.13,
            c: 1.38,
        }),
        17 => Some(GasteigerParams {
            a: 11.00,
            b: 9.69,
            c: 1.35,
        }),
        35 => Some(GasteigerParams {
            a: 10.08,
            b: 8.47,
            c: 1.16,
        }),
        53 => Some(GasteigerParams {
            a: 9.90,
            b: 7.96,
            c: 0.96,
        }),
        5 => Some(GasteigerParams {
            a: 5.98,
            b: 8.46,
            c: 1.70,
        }),
        14 => Some(GasteigerParams {
            a: 7.30,
            b: 6.56,
            c: 0.68,
        }),
        _ => None,
    }
}

/// Result of a Gasteiger-Marsili charge calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChargeResult {
    /// Partial charges per atom (same order as input).
    pub charges: Vec<f64>,
    /// Number of iterations actually performed.
    pub iterations: usize,
    /// Total charge (should be close to the formal charge sum).
    pub total_charge: f64,
}

/// Compute Gasteiger-Marsili partial charges for a molecule.
///
/// # Arguments
/// - `elements`: atomic numbers for each atom
/// - `bonds`: list of (atom_i, atom_j) index pairs (0-based)
/// - `formal_charges`: formal charge on each atom (usually 0)
/// - `max_iter`: maximum iterations (typically 6–8)
///
/// # Returns
/// `Ok(ChargeResult)` with partial charges, or `Err` if an element is unsupported.
pub fn gasteiger_marsili_charges(
    elements: &[u8],
    bonds: &[(usize, usize)],
    formal_charges: &[i8],
    max_iter: usize,
) -> Result<ChargeResult, String> {
    let n = elements.len();
    if formal_charges.len() != n {
        return Err(format!(
            "formal_charges length {} != elements length {}",
            formal_charges.len(),
            n
        ));
    }

    // Look up parameters for each atom
    let params: Vec<GasteigerParams> = elements
        .iter()
        .map(|&z| {
            get_gasteiger_params(z)
                .ok_or_else(|| format!("No Gasteiger parameters for element Z={}", z))
        })
        .collect::<Result<Vec<_>, _>>()?;

    // Initialize charges from formal charges
    let mut charges: Vec<f64> = formal_charges.iter().map(|&fc| fc as f64).collect();

    // Damping factor starts at 0.5 and halves each iteration
    let mut damping = 0.5;

    let mut actual_iters = 0;
    for _iter in 0..max_iter {
        actual_iters += 1;

        // Compute electronegativity for each atom at current charge
        let chi: Vec<f64> = (0..n).map(|i| params[i].chi(charges[i])).collect();

        // Compute charge transfer along each bond
        let mut delta_q = vec![0.0f64; n];
        for &(i, j) in bonds {
            if i >= n || j >= n {
                continue;
            }
            let chi_diff = chi[j] - chi[i];
            // Transfer is proportional to electronegativity difference
            // Normalize by the "hardness" of the receiving atom
            let dq = if chi_diff > 0.0 {
                // Charge flows from i to j (j is more electronegative)
                damping * chi_diff / params[j].chi(1.0)
            } else {
                // Charge flows from j to i
                damping * chi_diff / params[i].chi(1.0)
            };
            delta_q[i] += dq;
            delta_q[j] -= dq;
        }

        // Apply charge transfers
        let mut max_delta = 0.0f64;
        for i in 0..n {
            charges[i] += delta_q[i];
            max_delta = max_delta.max(delta_q[i].abs());
        }

        damping *= 0.5;

        // Convergence check
        if max_delta < 1e-10 {
            break;
        }
    }

    let total_charge: f64 = charges.iter().sum();

    Ok(ChargeResult {
        charges,
        iterations: actual_iters,
        total_charge,
    })
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_h2_symmetric_charges() {
        let elems = vec![1, 1];
        let bonds = vec![(0, 1)];
        let fc = vec![0, 0];
        let result = gasteiger_marsili_charges(&elems, &bonds, &fc, 6).unwrap();
        assert_eq!(result.charges.len(), 2);
        // Symmetric molecule → charges should be equal (near zero)
        assert!((result.charges[0] - result.charges[1]).abs() < 1e-6);
        assert!(result.charges[0].abs() < 0.01);
    }

    #[test]
    fn test_water_oxygen_negative() {
        // O—H  O—H
        let elems = vec![8, 1, 1];
        let bonds = vec![(0, 1), (0, 2)];
        let fc = vec![0, 0, 0];
        let result = gasteiger_marsili_charges(&elems, &bonds, &fc, 6).unwrap();
        // Oxygen more electronegative → negative charge
        assert!(
            result.charges[0] < -0.1,
            "O charge should be negative: {}",
            result.charges[0]
        );
        // Hydrogens should be positive
        assert!(result.charges[1] > 0.0);
        assert!(result.charges[2] > 0.0);
        // Total should be near zero
        assert!(result.total_charge.abs() < 1e-6);
    }

    #[test]
    fn test_methane_symmetric() {
        // CH₄: C bonded to 4 H
        let elems = vec![6, 1, 1, 1, 1];
        let bonds = vec![(0, 1), (0, 2), (0, 3), (0, 4)];
        let fc = vec![0, 0, 0, 0, 0];
        let result = gasteiger_marsili_charges(&elems, &bonds, &fc, 6).unwrap();
        // All H should have same charge (within tolerance)
        let h_charges: Vec<f64> = result.charges[1..].to_vec();
        for c in &h_charges {
            assert!(
                (c - h_charges[0]).abs() < 1e-10,
                "H charges should be equal"
            );
        }
        assert!(result.total_charge.abs() < 1e-6);
    }

    #[test]
    fn test_co2_carbon_positive() {
        // O=C=O
        let elems = vec![6, 8, 8];
        let bonds = vec![(0, 1), (0, 2)];
        let fc = vec![0, 0, 0];
        let result = gasteiger_marsili_charges(&elems, &bonds, &fc, 6).unwrap();
        // C surrounded by two O → positive
        assert!(
            result.charges[0] > 0.1,
            "C in CO₂ should be positive: {}",
            result.charges[0]
        );
        // Oxygens should be negative and equal
        assert!(result.charges[1] < -0.05);
        assert!((result.charges[1] - result.charges[2]).abs() < 1e-10);
    }

    #[test]
    fn test_hf_fluorine_negative() {
        let elems = vec![1, 9];
        let bonds = vec![(0, 1)];
        let fc = vec![0, 0];
        let result = gasteiger_marsili_charges(&elems, &bonds, &fc, 6).unwrap();
        assert!(
            result.charges[1] < -0.1,
            "F should be very negative: {}",
            result.charges[1]
        );
        assert!(result.charges[0] > 0.1);
        assert!(result.total_charge.abs() < 1e-6);
    }

    #[test]
    fn test_unsupported_element() {
        let elems = vec![2]; // Helium
        let bonds = vec![];
        let fc = vec![0];
        let result = gasteiger_marsili_charges(&elems, &bonds, &fc, 6);
        assert!(result.is_err());
    }

    #[test]
    fn test_formal_charge_preserved() {
        // NH₄⁺ (nitrogen with formal charge +1)
        let elems = vec![7, 1, 1, 1, 1];
        let bonds = vec![(0, 1), (0, 2), (0, 3), (0, 4)];
        let fc = vec![1, 0, 0, 0, 0];
        let result = gasteiger_marsili_charges(&elems, &bonds, &fc, 6).unwrap();
        // Total charge should be ~+1
        assert!(
            (result.total_charge - 1.0).abs() < 0.1,
            "Total charge should be ~+1: {}",
            result.total_charge
        );
    }

    #[test]
    fn test_electronegativity_order() {
        // Get χ at q=0 for different elements
        let h = get_gasteiger_params(1).unwrap().chi(0.0);
        let c = get_gasteiger_params(6).unwrap().chi(0.0);
        let n = get_gasteiger_params(7).unwrap().chi(0.0);
        let o = get_gasteiger_params(8).unwrap().chi(0.0);
        let f = get_gasteiger_params(9).unwrap().chi(0.0);
        // F > O > N > C > H (Pauling-like order)
        assert!(f > o);
        assert!(o > n);
        assert!(n > c);
        assert!(c > h);
    }
}
