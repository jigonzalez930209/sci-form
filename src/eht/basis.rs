//! Atomic orbital basis representations: STO and contracted Gaussian (STO-nG).
//!
//! Provides Slater-type orbitals (STOs) and their expansion into Gaussian primitives
//! for tractable overlap-integral evaluation.

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// A single Gaussian primitive: coefficient × exp(-alpha × r²).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaussianPrimitive {
    /// Contraction coefficient (pre-normalized).
    pub coeff: f64,
    /// Gaussian exponent (bohr⁻²).
    pub alpha: f64,
}

/// A Slater-type orbital represented conceptually.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SlaterOrbital {
    /// Principal quantum number.
    pub n: u8,
    /// Angular momentum l.
    pub l: u8,
    /// Magnetic quantum number m.
    pub m: i8,
    /// Slater exponent ζ.
    pub zeta: f64,
}

/// A single atomic orbital in the molecular basis, located on a specific atom.
#[derive(Debug, Clone)]
pub struct AtomicOrbital {
    /// Index of the atom this orbital belongs to.
    pub atom_index: usize,
    /// Center position in 3D space (bohr).
    pub center: [f64; 3],
    /// Principal quantum number.
    pub n: u8,
    /// Angular momentum l.
    pub l: u8,
    /// Magnetic quantum number m.
    pub m: i8,
    /// VSIP value (eV) for diagonal H_ii.
    pub vsip: f64,
    /// Slater exponent.
    pub zeta: f64,
    /// Label, e.g. "C_2s", "O_2px".
    pub label: String,
    /// Contracted Gaussian primitives (STO-nG expansion).
    pub gaussians: Vec<GaussianPrimitive>,
}

// ─── STO-3G contraction coefficients ─────────────────────────────────────────
// From Hehre, Stewart, Pople, J. Chem. Phys. 51, 2657 (1969).
// These are for the expansion of a Slater orbital with ζ=1.0.
// When ζ ≠ 1.0, scale alpha by ζ².

/// STO-3G coefficients for 1s orbital (n=1, l=0).
static STO3G_1S: [(f64, f64); 3] = [
    (0.154329, 2.227660),
    (0.535328, 0.405771),
    (0.444635, 0.109818),
];

/// STO-3G coefficients for 2s orbital (n=2, l=0).
/// Tabulated for ζ=1.0, two-term linear combination (2s = outer - inner).
static STO3G_2SP_INNER: [(f64, f64); 3] = [
    (-0.099967, 2.581058),
    (0.399513, 0.532013),
    (0.700115, 0.168856),
];

/// STO-3G coefficients for 2p orbital (n=2, l=1).
static STO3G_2SP_OUTER: [(f64, f64); 3] = [
    (0.155916, 2.581058),
    (0.607684, 0.532013),
    (0.391957, 0.168856),
];

/// STO-3G coefficients for 3sp orbital (n=3, l=0 or 1).
static STO3G_3SP_INNER: [(f64, f64); 3] = [
    (-0.219620, 0.994203),
    (0.225595, 0.231031),
    (0.900398, 0.075142),
];

static STO3G_3SP_OUTER: [(f64, f64); 3] = [
    (0.010588, 0.994203),
    (0.595167, 0.231031),
    (0.462001, 0.075142),
];

/// STO-3G for 4sp.
static STO3G_4SP_INNER: [(f64, f64); 3] = [
    (-0.310438, 0.494718),
    (0.015515, 0.120193),
    (1.023468, 0.040301),
];

static STO3G_4SP_OUTER: [(f64, f64); 3] = [
    (-0.063311, 0.494718),
    (0.574459, 0.120193),
    (0.499768, 0.040301),
];

/// STO-3G for 5sp.
static STO3G_5SP_INNER: [(f64, f64); 3] = [
    (-0.383204, 0.289515),
    (-0.159438, 0.073780),
    (1.143456, 0.025335),
];

static STO3G_5SP_OUTER: [(f64, f64); 3] = [
    (-0.094810, 0.289515),
    (0.570520, 0.073780),
    (0.497340, 0.025335),
];

/// Build STO-3G Gaussian primitives for an orbital with given n, l, and ζ.
pub fn sto3g_expansion(n: u8, l: u8, zeta: f64) -> Vec<GaussianPrimitive> {
    let zeta_sq = zeta * zeta;

    let table: &[(f64, f64); 3] = match (n, l) {
        (1, 0) => &STO3G_1S,
        (2, 0) => &STO3G_2SP_INNER,
        (2, 1) => &STO3G_2SP_OUTER,
        (3, 0) => &STO3G_3SP_INNER,
        (3, 1) => &STO3G_3SP_OUTER,
        (4, 0) => &STO3G_4SP_INNER,
        (4, 1) => &STO3G_4SP_OUTER,
        (5, 0) => &STO3G_5SP_INNER,
        (5, 1) => &STO3G_5SP_OUTER,
        _ => return vec![],
    };

    table
        .iter()
        .map(|&(coeff, alpha)| GaussianPrimitive {
            coeff,
            alpha: alpha * zeta_sq,
        })
        .collect()
}

/// Evaluate a normalized s-type Gaussian (2α/π)^{3/4} exp(-α r²) at distance² r2.
pub fn gaussian_s_value(alpha: f64, r2: f64) -> f64 {
    let norm = (2.0 * alpha / PI).powf(0.75);
    norm * (-alpha * r2).exp()
}

/// Evaluate a normalized p-type Gaussian component: N × x × exp(-α r²).
/// `component` is the Cartesian displacement (x, y, or z) from the center.
pub fn gaussian_p_value(alpha: f64, r2: f64, component: f64) -> f64 {
    let norm = (128.0 * alpha.powi(5) / (PI * PI * PI)).powf(0.25);
    norm * component * (-alpha * r2).exp()
}

/// Build the full molecular basis set from atom positions (in Angstrom) and elements.
/// Positions are converted internally to bohr for consistency.
pub fn build_basis(elements: &[u8], positions: &[[f64; 3]]) -> Vec<AtomicOrbital> {
    use super::params::get_params;

    let ang_to_bohr = 1.0 / 0.529177249;
    let mut basis = Vec::new();

    for (atom_idx, (&z, pos)) in elements.iter().zip(positions.iter()).enumerate() {
        let params = match get_params(z) {
            Some(p) => p,
            None => continue,
        };

        let center = [
            pos[0] * ang_to_bohr,
            pos[1] * ang_to_bohr,
            pos[2] * ang_to_bohr,
        ];

        let symbol = params.symbol;

        for orb_def in params.orbitals {
            if orb_def.l == 0 {
                // s orbital: single function
                basis.push(AtomicOrbital {
                    atom_index: atom_idx,
                    center,
                    n: orb_def.n,
                    l: 0,
                    m: 0,
                    vsip: orb_def.vsip,
                    zeta: orb_def.zeta,
                    label: format!("{}_{}", symbol, orb_def.label),
                    gaussians: sto3g_expansion(orb_def.n, 0, orb_def.zeta),
                });
            } else if orb_def.l == 1 {
                // p orbital: px, py, pz (m = -1, 0, +1)
                let p_labels = ["x", "y", "z"];
                let m_values: [i8; 3] = [-1, 0, 1];
                for (idx, &m) in m_values.iter().enumerate() {
                    basis.push(AtomicOrbital {
                        atom_index: atom_idx,
                        center,
                        n: orb_def.n,
                        l: 1,
                        m,
                        vsip: orb_def.vsip,
                        zeta: orb_def.zeta,
                        label: format!("{}_{}{}", symbol, orb_def.label, p_labels[idx]),
                        gaussians: sto3g_expansion(orb_def.n, 1, orb_def.zeta),
                    });
                }
            }
        }
    }

    basis
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sto3g_1s_expansion() {
        let gs = sto3g_expansion(1, 0, 1.0);
        assert_eq!(gs.len(), 3);
        // Coefficients sum roughly to 1 (not exact due to normalization)
        let sum: f64 = gs.iter().map(|g| g.coeff).sum();
        assert!((sum - 1.134292).abs() < 0.01);
    }

    #[test]
    fn test_sto3g_scaling() {
        // With ζ=2.0, alphas should be 4× the ζ=1.0 values
        let gs1 = sto3g_expansion(1, 0, 1.0);
        let gs2 = sto3g_expansion(1, 0, 2.0);
        for (a, b) in gs1.iter().zip(gs2.iter()) {
            assert!((b.alpha - a.alpha * 4.0).abs() < 1e-10);
            assert!((b.coeff - a.coeff).abs() < 1e-10);
        }
    }

    #[test]
    fn test_build_basis_h2() {
        // H₂ at 0.74 Å separation along x-axis
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let basis = build_basis(&elements, &positions);
        // H has 1 basis function (1s) each → 2 total
        assert_eq!(basis.len(), 2);
        assert_eq!(basis[0].label, "H_1s");
        assert_eq!(basis[1].label, "H_1s");
        assert_eq!(basis[0].atom_index, 0);
        assert_eq!(basis[1].atom_index, 1);
    }

    #[test]
    fn test_build_basis_h2o() {
        // H₂O: O at origin, two H's
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let basis = build_basis(&elements, &positions);
        // O: 2s + 2px + 2py + 2pz = 4, H: 1s each = 2.  Total = 6
        assert_eq!(basis.len(), 6);
    }

    #[test]
    fn test_gaussian_s_normalization() {
        // Integral of |g(r)|² over all space should be 1.0 for a normalized Gaussian
        // For (2α/π)^{3/2} exp(-2α r²), the integral is 1.
        // We test at the origin: value should be (2α/π)^{3/4}
        let alpha = 1.0;
        let val = gaussian_s_value(alpha, 0.0);
        let expected = (2.0 * alpha / PI).powf(0.75);
        assert!((val - expected).abs() < 1e-12);
    }

    #[test]
    fn test_gaussian_s_decay() {
        let alpha = 1.0;
        let v0 = gaussian_s_value(alpha, 0.0);
        let v1 = gaussian_s_value(alpha, 1.0);
        let v5 = gaussian_s_value(alpha, 5.0);
        assert!(v0 > v1);
        assert!(v1 > v5);
        assert!(v5 > 0.0);
    }
}
