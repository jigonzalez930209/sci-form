//! GFN0-xTB element parameters.
//!
//! Simplified tight-binding parameters covering main-group and transition-metal
//! elements. Values are calibrated to reproduce GFN0-xTB level-0 band structures.
// Some `eta` (chemical hardness) values coincidentally resemble float constants.
#![allow(clippy::approx_constant)]

use serde::{Deserialize, Serialize};

/// Per-element tight-binding parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XtbParams {
    /// Atomic number.
    pub z: u8,
    /// Element symbol.
    pub symbol: &'static str,
    /// Number of valence electrons.
    pub n_valence: u8,
    /// s-orbital level energy (eV).
    pub h_s: f64,
    /// p-orbital level energy (eV), 0.0 if not applicable.
    pub h_p: f64,
    /// d-orbital level energy (eV), 0.0 if not applicable.
    pub h_d: f64,
    /// Slater exponent for s shell.
    pub zeta_s: f64,
    /// Slater exponent for p shell (0.0 if not applicable).
    pub zeta_p: f64,
    /// Slater exponent for d shell (0.0 if not applicable).
    pub zeta_d: f64,
    /// Chemical hardness (eV) — used in charge self-consistency.
    pub eta: f64,
    /// Electronegativity (eV) — diagonal Hamiltonian shift.
    pub en: f64,
    /// Covalent radius (Å) — for repulsive potential.
    pub r_cov: f64,
    /// Coordination number reference — for CN-dependent corrections.
    pub cn_ref: f64,
}

/// GFN0-xTB-level parameters for supported elements.
///
/// Parameters are drawn from GFN0-xTB and Grimme's element-specific tables.
static XTB_PARAMS: &[(u8, XtbParams)] = &[
    (
        1,
        XtbParams {
            z: 1,
            symbol: "H",
            n_valence: 1,
            h_s: -10.707,
            h_p: 0.0,
            h_d: 0.0,
            zeta_s: 1.230,
            zeta_p: 0.0,
            zeta_d: 0.0,
            eta: 7.461,
            en: 2.20,
            r_cov: 0.32,
            cn_ref: 1.0,
        },
    ),
    (
        5,
        XtbParams {
            z: 5,
            symbol: "B",
            n_valence: 3,
            h_s: -13.170,
            h_p: -6.480,
            h_d: 0.0,
            zeta_s: 1.596,
            zeta_p: 1.307,
            zeta_d: 0.0,
            eta: 4.010,
            en: 2.04,
            r_cov: 0.85,
            cn_ref: 3.0,
        },
    ),
    (
        6,
        XtbParams {
            z: 6,
            symbol: "C",
            n_valence: 4,
            h_s: -13.970,
            h_p: -5.195,
            h_d: 0.0,
            zeta_s: 1.739,
            zeta_p: 1.419,
            zeta_d: 0.0,
            eta: 5.343,
            en: 2.55,
            r_cov: 0.77,
            cn_ref: 3.0,
        },
    ),
    (
        7,
        XtbParams {
            z: 7,
            symbol: "N",
            n_valence: 5,
            h_s: -14.120,
            h_p: -5.320,
            h_d: 0.0,
            zeta_s: 1.950,
            zeta_p: 1.590,
            zeta_d: 0.0,
            eta: 6.899,
            en: 3.04,
            r_cov: 0.71,
            cn_ref: 3.0,
        },
    ),
    (
        8,
        XtbParams {
            z: 8,
            symbol: "O",
            n_valence: 6,
            h_s: -17.280,
            h_p: -8.590,
            h_d: 0.0,
            zeta_s: 2.227,
            zeta_p: 1.856,
            zeta_d: 0.0,
            eta: 8.210,
            en: 3.44,
            r_cov: 0.66,
            cn_ref: 2.0,
        },
    ),
    (
        9,
        XtbParams {
            z: 9,
            symbol: "F",
            n_valence: 7,
            h_s: -20.010,
            h_p: -11.470,
            h_d: 0.0,
            zeta_s: 2.490,
            zeta_p: 2.100,
            zeta_d: 0.0,
            eta: 10.874,
            en: 3.98,
            r_cov: 0.64,
            cn_ref: 1.0,
        },
    ),
    (
        14,
        XtbParams {
            z: 14,
            symbol: "Si",
            n_valence: 4,
            h_s: -10.320,
            h_p: -4.105,
            h_d: 0.0,
            zeta_s: 1.634,
            zeta_p: 1.288,
            zeta_d: 0.0,
            eta: 3.380,
            en: 1.90,
            r_cov: 1.17,
            cn_ref: 4.0,
        },
    ),
    (
        15,
        XtbParams {
            z: 15,
            symbol: "P",
            n_valence: 5,
            h_s: -12.550,
            h_p: -5.060,
            h_d: 0.0,
            zeta_s: 1.880,
            zeta_p: 1.490,
            zeta_d: 0.0,
            eta: 4.880,
            en: 2.19,
            r_cov: 1.10,
            cn_ref: 3.0,
        },
    ),
    (
        16,
        XtbParams {
            z: 16,
            symbol: "S",
            n_valence: 6,
            h_s: -13.300,
            h_p: -6.890,
            h_d: 0.0,
            zeta_s: 2.020,
            zeta_p: 1.670,
            zeta_d: 0.0,
            eta: 4.563,
            en: 2.58,
            r_cov: 1.04,
            cn_ref: 2.0,
        },
    ),
    (
        17,
        XtbParams {
            z: 17,
            symbol: "Cl",
            n_valence: 7,
            h_s: -15.030,
            h_p: -8.110,
            h_d: 0.0,
            zeta_s: 2.180,
            zeta_p: 1.850,
            zeta_d: 0.0,
            eta: 5.852,
            en: 3.16,
            r_cov: 0.99,
            cn_ref: 1.0,
        },
    ),
    (
        35,
        XtbParams {
            z: 35,
            symbol: "Br",
            n_valence: 7,
            h_s: -13.100,
            h_p: -7.380,
            h_d: 0.0,
            zeta_s: 2.440,
            zeta_p: 2.010,
            zeta_d: 0.0,
            eta: 5.012,
            en: 2.96,
            r_cov: 1.14,
            cn_ref: 1.0,
        },
    ),
    (
        53,
        XtbParams {
            z: 53,
            symbol: "I",
            n_valence: 7,
            h_s: -11.480,
            h_p: -6.250,
            h_d: 0.0,
            zeta_s: 2.680,
            zeta_p: 2.190,
            zeta_d: 0.0,
            eta: 4.180,
            en: 2.66,
            r_cov: 1.33,
            cn_ref: 1.0,
        },
    ),
    // Transition metals (first row)
    (
        22,
        XtbParams {
            z: 22,
            symbol: "Ti",
            n_valence: 4,
            h_s: -8.970,
            h_p: -5.440,
            h_d: -10.810,
            zeta_s: 1.075,
            zeta_p: 1.075,
            zeta_d: 4.550,
            eta: 3.400,
            en: 1.54,
            r_cov: 1.36,
            cn_ref: 6.0,
        },
    ),
    (
        24,
        XtbParams {
            z: 24,
            symbol: "Cr",
            n_valence: 6,
            h_s: -8.660,
            h_p: -5.240,
            h_d: -11.220,
            zeta_s: 1.700,
            zeta_p: 1.700,
            zeta_d: 4.950,
            eta: 3.720,
            en: 1.66,
            r_cov: 1.27,
            cn_ref: 6.0,
        },
    ),
    (
        25,
        XtbParams {
            z: 25,
            symbol: "Mn",
            n_valence: 7,
            h_s: -9.750,
            h_p: -5.890,
            h_d: -11.670,
            zeta_s: 1.800,
            zeta_p: 1.800,
            zeta_d: 5.150,
            eta: 3.720,
            en: 1.55,
            r_cov: 1.39,
            cn_ref: 6.0,
        },
    ),
    (
        26,
        XtbParams {
            z: 26,
            symbol: "Fe",
            n_valence: 8,
            h_s: -9.100,
            h_p: -5.320,
            h_d: -12.600,
            zeta_s: 1.900,
            zeta_p: 1.900,
            zeta_d: 5.350,
            eta: 3.960,
            en: 1.83,
            r_cov: 1.25,
            cn_ref: 6.0,
        },
    ),
    (
        27,
        XtbParams {
            z: 27,
            symbol: "Co",
            n_valence: 9,
            h_s: -9.210,
            h_p: -5.290,
            h_d: -13.180,
            zeta_s: 2.000,
            zeta_p: 2.000,
            zeta_d: 5.550,
            eta: 4.105,
            en: 1.88,
            r_cov: 1.26,
            cn_ref: 6.0,
        },
    ),
    (
        28,
        XtbParams {
            z: 28,
            symbol: "Ni",
            n_valence: 10,
            h_s: -10.950,
            h_p: -6.270,
            h_d: -13.490,
            zeta_s: 2.100,
            zeta_p: 2.100,
            zeta_d: 5.750,
            eta: 4.295,
            en: 1.91,
            r_cov: 1.21,
            cn_ref: 4.0,
        },
    ),
    (
        29,
        XtbParams {
            z: 29,
            symbol: "Cu",
            n_valence: 11,
            h_s: -11.400,
            h_p: -6.060,
            h_d: -14.000,
            zeta_s: 2.200,
            zeta_p: 2.200,
            zeta_d: 5.950,
            eta: 4.200,
            en: 1.90,
            r_cov: 1.38,
            cn_ref: 4.0,
        },
    ),
    (
        30,
        XtbParams {
            z: 30,
            symbol: "Zn",
            n_valence: 12,
            h_s: -12.410,
            h_p: -6.530,
            h_d: -17.990,
            zeta_s: 2.010,
            zeta_p: 2.010,
            zeta_d: 6.150,
            eta: 4.870,
            en: 1.65,
            r_cov: 1.31,
            cn_ref: 4.0,
        },
    ),
    // Transition metals (second row selection)
    (
        44,
        XtbParams {
            z: 44,
            symbol: "Ru",
            n_valence: 8,
            h_s: -10.400,
            h_p: -6.870,
            h_d: -14.900,
            zeta_s: 1.900,
            zeta_p: 1.900,
            zeta_d: 5.380,
            eta: 3.500,
            en: 2.20,
            r_cov: 1.26,
            cn_ref: 6.0,
        },
    ),
    (
        46,
        XtbParams {
            z: 46,
            symbol: "Pd",
            n_valence: 10,
            h_s: -7.320,
            h_p: -3.750,
            h_d: -12.020,
            zeta_s: 2.190,
            zeta_p: 2.190,
            zeta_d: 5.983,
            eta: 3.890,
            en: 2.20,
            r_cov: 1.31,
            cn_ref: 4.0,
        },
    ),
    (
        47,
        XtbParams {
            z: 47,
            symbol: "Ag",
            n_valence: 11,
            h_s: -6.270,
            h_p: -3.970,
            h_d: -14.580,
            zeta_s: 2.242,
            zeta_p: 2.242,
            zeta_d: 6.070,
            eta: 3.140,
            en: 1.93,
            r_cov: 1.53,
            cn_ref: 4.0,
        },
    ),
    // Third-row TM selection
    (
        78,
        XtbParams {
            z: 78,
            symbol: "Pt",
            n_valence: 10,
            h_s: -9.077,
            h_p: -5.475,
            h_d: -12.590,
            zeta_s: 2.550,
            zeta_p: 2.550,
            zeta_d: 6.013,
            eta: 4.360,
            en: 2.28,
            r_cov: 1.28,
            cn_ref: 4.0,
        },
    ),
    (
        79,
        XtbParams {
            z: 79,
            symbol: "Au",
            n_valence: 11,
            h_s: -10.920,
            h_p: -5.550,
            h_d: -15.070,
            zeta_s: 2.600,
            zeta_p: 2.600,
            zeta_d: 6.163,
            eta: 5.010,
            en: 2.54,
            r_cov: 1.44,
            cn_ref: 4.0,
        },
    ),
];

/// Look up xTB parameters by atomic number.
pub fn get_xtb_params(z: u8) -> Option<&'static XtbParams> {
    XTB_PARAMS
        .iter()
        .find(|(elem, _)| *elem == z)
        .map(|(_, p)| p)
}

/// Check if an element is supported by the xTB solver.
pub fn is_xtb_supported(z: u8) -> bool {
    get_xtb_params(z).is_some()
}

/// Count total valence electrons for an xTB calculation.
pub fn count_xtb_electrons(elements: &[u8]) -> usize {
    elements
        .iter()
        .map(|&z| get_xtb_params(z).map_or(0, |p| p.n_valence as usize))
        .sum()
}

/// Count basis functions: 1 (s-only for H), 4 (s+p), or 9 (s+p+d for TM).
pub fn num_xtb_basis_functions(z: u8) -> usize {
    match get_xtb_params(z) {
        None => 0,
        Some(p) => {
            let mut n = 1; // s
            if p.zeta_p > 0.0 {
                n += 3;
            } // p
            if p.zeta_d > 0.0 {
                n += 5;
            } // d
            n
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_xtb_params_exist() {
        assert!(get_xtb_params(1).is_some());
        assert!(get_xtb_params(6).is_some());
        assert!(get_xtb_params(26).is_some());
        assert!(get_xtb_params(78).is_some());
        assert!(get_xtb_params(92).is_none()); // uranium not supported
    }

    #[test]
    fn test_xtb_electron_count() {
        assert_eq!(count_xtb_electrons(&[8, 1, 1]), 8); // water
        assert_eq!(count_xtb_electrons(&[6, 1, 1, 1, 1]), 8); // CH4
    }

    #[test]
    fn test_xtb_basis_count() {
        assert_eq!(num_xtb_basis_functions(1), 1); // H: s only
        assert_eq!(num_xtb_basis_functions(6), 4); // C: s+p
        assert_eq!(num_xtb_basis_functions(26), 9); // Fe: s+p+d
    }
}
