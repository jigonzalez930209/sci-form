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
/// Gasteiger & Marsili (1980), OpenBabel implementation, and extended coverage
/// for main-group elements up to period 4.
pub fn get_gasteiger_params(z: u8) -> Option<GasteigerParams> {
    match z {
        // Period 1
        1 => Some(GasteigerParams {
            a: 7.17,
            b: 6.24,
            c: -0.56,
        }),
        // Period 2
        3 => Some(GasteigerParams {
            // Li
            a: 3.01,
            b: 2.78,
            c: 0.30,
        }),
        4 => Some(GasteigerParams {
            // Be
            a: 4.90,
            b: 4.52,
            c: 0.55,
        }),
        5 => Some(GasteigerParams {
            a: 5.98,
            b: 8.46,
            c: 1.70,
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
        // Period 3
        11 => Some(GasteigerParams {
            // Na
            a: 2.84,
            b: 3.00,
            c: 0.40,
        }),
        12 => Some(GasteigerParams {
            // Mg
            a: 3.75,
            b: 3.90,
            c: 0.50,
        }),
        13 => Some(GasteigerParams {
            // Al
            a: 5.47,
            b: 5.10,
            c: 0.65,
        }),
        14 => Some(GasteigerParams {
            a: 7.30,
            b: 6.56,
            c: 0.68,
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
        // Period 4 (selected)
        19 => Some(GasteigerParams {
            // K
            a: 2.42,
            b: 2.60,
            c: 0.35,
        }),
        20 => Some(GasteigerParams {
            // Ca
            a: 3.23,
            b: 3.40,
            c: 0.45,
        }),
        31 => Some(GasteigerParams {
            // Ga
            a: 5.20,
            b: 5.50,
            c: 0.70,
        }),
        32 => Some(GasteigerParams {
            // Ge
            a: 6.90,
            b: 6.80,
            c: 0.80,
        }),
        33 => Some(GasteigerParams {
            // As
            a: 8.30,
            b: 7.80,
            c: 0.90,
        }),
        34 => Some(GasteigerParams {
            // Se
            a: 9.50,
            b: 8.70,
            c: 1.10,
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
        _ => None,
    }
}

/// Hybridization/bonding mode used to select PEOE parameters.
///
/// Mirrors the mode keys of the canonical Gasteiger parameter table
/// (`sp`, `sp2`, `sp3`, sulfur `so`/`so2`, hydrogen `*`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GasteigerMode {
    /// Unknown hybridization: use the element default.
    Default,
    /// Hydrogen.
    Hydrogen,
    /// sp hybridization.
    Sp,
    /// sp2 hybridization.
    Sp2,
    /// sp3 hybridization.
    Sp3,
    /// Sulfur with a single S=O (used when the hybridization is unknown).
    So,
    /// Sulfur with two S=O (used when the hybridization is unknown).
    So2,
}

/// Mode-aware Gasteiger parameter lookup.
///
/// Values follow the canonical Gasteiger-Marsili table (as also published in
/// RDKit's `GasteigerParams.cpp`) so the iterative PEOE implementation can be
/// cross-validated against reference implementations. Elements absent from the
/// canonical table fall back to the extended table used by
/// [`get_gasteiger_params`].
pub fn get_gasteiger_params_mode(z: u8, mode: GasteigerMode) -> Option<GasteigerParams> {
    match (z, mode) {
        (1, _) | (_, GasteigerMode::Hydrogen) => Some(GasteigerParams {
            a: 7.17,
            b: 6.24,
            c: -0.56,
        }),
        (6, GasteigerMode::Sp) => Some(GasteigerParams {
            a: 10.39,
            b: 9.45,
            c: 0.73,
        }),
        (6, GasteigerMode::Sp2) => Some(GasteigerParams {
            a: 8.79,
            b: 9.32,
            c: 1.51,
        }),
        (6, GasteigerMode::Sp3 | GasteigerMode::Default) => Some(GasteigerParams {
            a: 7.98,
            b: 9.18,
            c: 1.88,
        }),
        (7, GasteigerMode::Sp) => Some(GasteigerParams {
            a: 15.68,
            b: 11.7,
            c: -0.27,
        }),
        (7, GasteigerMode::Sp2) => Some(GasteigerParams {
            a: 12.87,
            b: 11.15,
            c: 0.85,
        }),
        (7, GasteigerMode::Sp3 | GasteigerMode::Default) => Some(GasteigerParams {
            a: 11.54,
            b: 10.82,
            c: 1.36,
        }),
        (8, GasteigerMode::Sp2) => Some(GasteigerParams {
            a: 17.07,
            b: 13.79,
            c: 0.47,
        }),
        (8, GasteigerMode::Sp3 | GasteigerMode::Default) => Some(GasteigerParams {
            a: 14.18,
            b: 12.92,
            c: 1.39,
        }),
        (15, GasteigerMode::Sp2) => Some(GasteigerParams {
            a: 9.665,
            b: 8.530,
            c: 0.735,
        }),
        (5, GasteigerMode::Sp2) => Some(GasteigerParams {
            a: 6.420,
            b: 6.807,
            c: 1.322,
        }),
        (5, GasteigerMode::Sp3 | GasteigerMode::Default) => Some(GasteigerParams {
            a: 5.980,
            b: 6.820,
            c: 1.605,
        }),
        (14, GasteigerMode::Sp) => Some(GasteigerParams {
            a: 9.065,
            b: 7.027,
            c: -0.002,
        }),
        (14, GasteigerMode::Sp2) => Some(GasteigerParams {
            a: 7.905,
            b: 6.748,
            c: 0.443,
        }),
        (14, GasteigerMode::Sp3 | GasteigerMode::Default) => Some(GasteigerParams {
            a: 7.300,
            b: 6.567,
            c: 0.657,
        }),
        (4, GasteigerMode::Sp2) => Some(GasteigerParams {
            a: 4.005,
            b: 6.725,
            c: 3.035,
        }),
        (4, GasteigerMode::Sp3 | GasteigerMode::Default) => Some(GasteigerParams {
            a: 3.845,
            b: 6.755,
            c: 3.165,
        }),
        (12, GasteigerMode::Sp) => Some(GasteigerParams {
            a: 4.040,
            b: 5.472,
            c: 1.823,
        }),
        (12, GasteigerMode::Sp2) => Some(GasteigerParams {
            a: 3.565,
            b: 5.572,
            c: 2.197,
        }),
        (12, GasteigerMode::Sp3 | GasteigerMode::Default) => Some(GasteigerParams {
            a: 3.300,
            b: 5.587,
            c: 2.447,
        }),
        (13, GasteigerMode::Sp2) => Some(GasteigerParams {
            a: 5.795,
            b: 5.020,
            c: 0.695,
        }),
        (13, GasteigerMode::Sp3 | GasteigerMode::Default) => Some(GasteigerParams {
            a: 5.375,
            b: 4.953,
            c: 0.867,
        }),
        (16, GasteigerMode::Sp2) => Some(GasteigerParams {
            a: 10.88,
            b: 9.49,
            c: 1.33,
        }),
        (16, GasteigerMode::Sp3 | GasteigerMode::Default | GasteigerMode::So) => {
            Some(GasteigerParams {
                a: 10.14,
                b: 9.13,
                c: 1.38,
            })
        }
        (16, GasteigerMode::So2) => Some(GasteigerParams {
            a: 12.00,
            b: 10.81,
            c: 1.20,
        }),
        _ => get_gasteiger_params(z),
    }
}

/// Constant `ionX` for hydrogen used by the canonical PEOE iteration.
const IONX_HYDROGEN: f64 = 20.02;

/// True when the bond `(i, j, order)` participates in a conjugated system.
///
/// A multiple/aromatic bond is always conjugated; a single bond is conjugated
/// when either endpoint carries a multiple bond to a heavy atom (the classic
/// amide/ester/carboxylate case). Multiple bonds on hypervalent centres
/// (P, S) do not propagate conjugation, matching RDKit's `setConjugation`.
fn bond_is_conjugated(
    elements: &[u8],
    adjacency: &[Vec<(usize, u8)>],
    i: usize,
    j: usize,
    order: u8,
) -> bool {
    if order >= 2 {
        return true;
    }
    let has_multiple = |atom: usize| {
        if elements.get(atom).copied().unwrap_or(0) == 15
            || elements.get(atom).copied().unwrap_or(0) == 16
        {
            return false;
        }
        adjacency[atom]
            .iter()
            .any(|&(k, o)| o >= 2 && elements.get(k).copied().unwrap_or(1) != 1)
    };
    has_multiple(i) || has_multiple(j)
}

/// Split formal charges evenly across same-element conjugated partners.
///
/// Mirrors `Gasteiger::splitChargeConjugated` from the canonical
/// implementation: for each atom carrying a formal charge, walk two conjugated
/// bonds; every same-element atom reached receives an equal share of the total
/// formal charge of the group. This is what turns a carboxylate `C(=O)O⁻` into
/// two equivalent oxygens with charge −0.5.
fn split_charge_conjugated(
    elements: &[u8],
    bonds: &[(usize, usize, u8)],
    formal_charges: &[i8],
    charges: &mut [f64],
) {
    let n = elements.len();
    let mut adjacency: Vec<Vec<(usize, u8)>> = vec![Vec::new(); n];
    for &(i, j, order) in bonds {
        if i < n && j < n {
            adjacency[i].push((j, order));
            adjacency[j].push((i, order));
        }
    }

    for i in 0..n {
        let formal = formal_charges[i] as f64;
        if formal.abs() <= 1e-12 || charges[i].abs() > 1e-12 {
            continue;
        }
        let mut marker = vec![i];
        let mut total = formal;
        for &(mid, order1) in &adjacency[i] {
            if !bond_is_conjugated(elements, &adjacency, i, mid, order1) {
                continue;
            }
            for &(other, order2) in &adjacency[mid] {
                if other == i {
                    continue;
                }
                if !bond_is_conjugated(elements, &adjacency, mid, other, order2) {
                    continue;
                }
                if elements[other] == elements[i] {
                    total += formal_charges[other] as f64;
                    marker.push(other);
                }
            }
        }
        let share = total / marker.len() as f64;
        for &m in &marker {
            charges[m] = share;
        }
    }
}

/// Full Gasteiger-Marsili PEOE implementation with per-atom modes and bond orders.
///
/// Reproduces the canonical iteration:
///
/// - orbital electronegativity `χ_i = a + q·(b + c·q)`
/// - charge transfer along each bond `dq = (χ_j − χ_i) / ionX`, where `ionX` is
///   the electronegativity at unit charge of the *less* electronegative atom
///   (hydrogen uses the constant 20.02)
/// - damping starts at `config.initial_damping` (0.5) and halves each iteration
/// - formal charges are split evenly across same-element conjugated partners
///   before the first iteration
///
/// `bonds` carries the bond order (`1` single, `2` double, `3` triple,
/// `4` aromatic). `modes` optionally gives the per-atom hybridization; when it
/// is `None` every atom uses its element default (sp3-like), which preserves
/// the behavior of the simpler `gasteiger_marsili_charges` API.
pub fn gasteiger_marsili_charges_full(
    elements: &[u8],
    bonds: &[(usize, usize, u8)],
    formal_charges: &[i8],
    modes: Option<&[GasteigerMode]>,
    config: &GasteigerConfig,
) -> Result<ChargeResult, String> {
    config.validate()?;
    let n = elements.len();
    if formal_charges.len() != n {
        return Err(format!(
            "formal_charges length {} != elements length {}",
            formal_charges.len(),
            n
        ));
    }
    if let Some(modes) = modes {
        if modes.len() != n {
            return Err(format!(
                "modes length {} != elements length {}",
                modes.len(),
                n
            ));
        }
    }
    for &(i, j, _) in bonds {
        if i >= n || j >= n {
            return Err(format!(
                "Bond ({}, {}) references atom outside range 0..{}",
                i, j, n
            ));
        }
    }

    // Per-atom parameters, mirroring the reference mode selection.
    let mut params = Vec::with_capacity(n);
    for i in 0..n {
        let requested = modes.map(|m| m[i]).unwrap_or(GasteigerMode::Default);
        let mode = if elements[i] == 1 {
            GasteigerMode::Hydrogen
        } else if matches!(
            requested,
            GasteigerMode::Sp | GasteigerMode::Sp2 | GasteigerMode::Sp3
        ) {
            requested
        } else if elements[i] == 16 && requested == GasteigerMode::Default {
            let n_oxo = bonds
                .iter()
                .filter(|&&(a, b, _)| (a == i && elements[b] == 8) || (b == i && elements[a] == 8))
                .count();
            match n_oxo {
                2 => GasteigerMode::So2,
                1 => GasteigerMode::So,
                _ => GasteigerMode::Sp3,
            }
        } else {
            requested
        };
        params.push(
            get_gasteiger_params_mode(elements[i], mode)
                .ok_or_else(|| format!("No Gasteiger parameters for element Z={}", elements[i]))?,
        );
    }

    // ionX: electronegativity at unit charge; hydrogen uses the constant.
    let ion_x: Vec<f64> = elements
        .iter()
        .zip(params.iter())
        .map(|(&z, p)| {
            if z == 1 {
                IONX_HYDROGEN
            } else {
                p.a + p.b + p.c
            }
        })
        .collect();

    // Initial charges: formal charges, split across conjugated partners.
    let mut charges = vec![0.0f64; n];
    split_charge_conjugated(elements, bonds, formal_charges, &mut charges);

    let mut damping = config.initial_damping;
    let mut max_delta = 0.0f64;
    for _ in 0..config.max_iter {
        let energ: Vec<f64> = (0..n)
            .map(|i| {
                let p = params[i];
                p.a + charges[i] * (p.b + p.c * charges[i])
            })
            .collect();

        let mut delta_q = vec![0.0f64; n];
        for &(i, j, _) in bonds {
            let dx = energ[j] - energ[i];
            let divisor = if dx >= 0.0 { ion_x[i] } else { ion_x[j] };
            if divisor.abs() < 1e-12 {
                continue;
            }
            let transfer = dx / divisor;
            delta_q[i] += transfer;
            delta_q[j] -= transfer;
        }

        max_delta = 0.0;
        for i in 0..n {
            charges[i] += damping * delta_q[i];
            max_delta = max_delta.max((damping * delta_q[i]).abs());
        }
        damping *= 0.5;
    }

    let total_charge: f64 = charges.iter().sum();

    Ok(ChargeResult {
        charges,
        iterations: config.max_iter,
        total_charge,
        converged: max_delta < config.convergence_threshold,
    })
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
    /// Whether the iterative process converged within the tolerance.
    pub converged: bool,
}

/// Configuration for Gasteiger-Marsili charge calculation.
#[derive(Debug, Clone)]
pub struct GasteigerConfig {
    /// Maximum number of iterations (default: 6).
    pub max_iter: usize,
    /// Initial damping factor (default: 0.5, halved each iteration).
    pub initial_damping: f64,
    /// Convergence threshold for max charge transfer (default: 1e-10).
    pub convergence_threshold: f64,
}

impl Default for GasteigerConfig {
    fn default() -> Self {
        GasteigerConfig {
            max_iter: 6,
            initial_damping: 0.5,
            convergence_threshold: 1e-10,
        }
    }
}

impl GasteigerConfig {
    /// Validate the configuration.
    ///
    /// A degenerate configuration used to produce a *silently wrong* result:
    /// `max_iter = 0` skipped the iteration loop entirely and returned a
    /// vector of zeros with `converged = true`, and `initial_damping = 0.0`
    /// froze the charges at their formal values while still claiming
    /// convergence.  Both are rejected instead.
    pub fn validate(&self) -> Result<(), String> {
        if self.max_iter == 0 {
            return Err(format!("max_iter must be >= 1 (got {})", self.max_iter));
        }
        if !(self.initial_damping > 0.0 && self.initial_damping <= 1.0) {
            return Err(format!(
                "initial_damping must be in (0, 1] (got {})",
                self.initial_damping
            ));
        }
        if !self.convergence_threshold.is_finite() || self.convergence_threshold < 0.0 {
            return Err(format!(
                "convergence_threshold must be a finite, non-negative number (got {})",
                self.convergence_threshold
            ));
        }
        Ok(())
    }
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
    let config = GasteigerConfig {
        max_iter,
        ..Default::default()
    };
    gasteiger_marsili_charges_configured(elements, bonds, formal_charges, &config)
}

/// Compute Gasteiger-Marsili partial charges with full configuration.
pub fn gasteiger_marsili_charges_configured(
    elements: &[u8],
    bonds: &[(usize, usize)],
    formal_charges: &[i8],
    config: &GasteigerConfig,
) -> Result<ChargeResult, String> {
    let bonds_with_order: Vec<(usize, usize, u8)> =
        bonds.iter().map(|&(i, j)| (i, j, 1u8)).collect();
    gasteiger_marsili_charges_full(elements, &bonds_with_order, formal_charges, None, config)
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
