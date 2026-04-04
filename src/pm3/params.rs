//! PM3 parameter tables.
//!
//! Standard PM3 parameters from Stewart, J. J. P. J. Comput. Chem. 10 (1989): 209.
//! Only the most chemically useful elements are parameterized.

#![allow(clippy::approx_constant)]

use serde::Serialize;

/// PM3 atomic parameters for one element.
#[derive(Debug, Clone, Serialize)]
pub struct Pm3Params {
    pub z: u8,
    pub symbol: &'static str,
    /// One-center one-electron integral for s orbital (eV).
    pub uss: f64,
    /// One-center one-electron integral for p orbital (eV).
    pub upp: f64,
    /// Beta_s resonance integral (eV).
    pub beta_s: f64,
    /// Beta_p resonance integral (eV).
    pub beta_p: f64,
    /// Slater exponent for s orbital (bohr⁻¹).
    pub zeta_s: f64,
    /// Slater exponent for p orbital (bohr⁻¹).
    pub zeta_p: f64,
    /// One-center two-electron integral g_ss (eV).
    pub gss: f64,
    /// g_sp (eV).
    pub gsp: f64,
    /// g_pp (eV).
    pub gpp: f64,
    /// g_p2 (eV).
    pub gp2: f64,
    /// h_sp (eV).
    pub hsp: f64,
    /// Core charge (number of valence electrons).
    pub core_charge: f64,
    /// Heat of atomization (kcal/mol).
    pub heat_of_atomization: f64,
    /// Alpha parameter for core-core repulsion (Å⁻¹).
    pub alpha: f64,
    /// PM3 Gaussian correction terms for core-core repulsion.
    /// Each tuple is (a_k, b_k, c_k) in the expression:
    ///   ΔE_nuc += Z_A × Z_B × Σ_k a_k × exp(−b_k × (R_AB − c_k)²)
    /// Typically 2 Gaussians per element; empty for elements without corrections.
    pub gaussians: &'static [(f64, f64, f64)],
}

impl Pm3Params {
    /// Reference distance for PM3 Gaussian core-core correction (Å).
    /// Uses covalent radius sum as a reasonable default for R₀.
    pub fn alpha_r0(&self) -> f64 {
        covalent_radius_angstrom(self.z)
    }
}

/// Approximate covalent radius in Å for PM3-supported elements.
fn covalent_radius_angstrom(z: u8) -> f64 {
    match z {
        1 => 0.31,  // H
        5 => 0.84,  // B
        6 => 0.76,  // C
        7 => 0.71,  // N
        8 => 0.66,  // O
        9 => 0.57,  // F
        13 => 1.21, // Al
        14 => 1.11, // Si
        15 => 1.07, // P
        16 => 1.05, // S
        17 => 1.02, // Cl
        32 => 1.20, // Ge
        35 => 1.20, // Br
        53 => 1.39, // I
        22 => 1.60, // Ti
        24 => 1.39, // Cr
        25 => 1.39, // Mn
        26 => 1.32, // Fe
        27 => 1.26, // Co
        28 => 1.24, // Ni
        29 => 1.32, // Cu
        30 => 1.22, // Zn
        _ => 1.50,  // fallback
    }
}

/// Check if an element has PM3 parameters.
pub fn is_pm3_supported(z: u8) -> bool {
    get_pm3_params(z).is_some()
}

/// Look up PM3 parameters by atomic number.
pub fn get_pm3_params(z: u8) -> Option<&'static Pm3Params> {
    ALL_PM3_PARAMS.iter().find(|p| p.z == z)
}

// Stewart, J. J. P. J. Comput. Chem. 10 (1989): 209–220.
static ALL_PM3_PARAMS: &[Pm3Params] = &[
    // Hydrogen
    Pm3Params {
        z: 1,
        symbol: "H",
        uss: -13.073321,
        upp: 0.0,
        beta_s: -5.626512,
        beta_p: 0.0,
        zeta_s: 0.967807,
        zeta_p: 0.0,
        gss: 14.794208,
        gsp: 0.0,
        gpp: 0.0,
        gp2: 0.0,
        hsp: 0.0,
        core_charge: 1.0,
        heat_of_atomization: 52.102,
        alpha: 3.356386,
        gaussians: &[(1.1287500, 5.0962820, 1.5374650), (-1.0603290, 6.0037880, 1.5701890)],
    },
    // Carbon
    Pm3Params {
        z: 6,
        symbol: "C",
        uss: -47.270320,
        upp: -36.266918,
        beta_s: -11.910015,
        beta_p: -9.802755,
        zeta_s: 1.565085,
        zeta_p: 1.842345,
        gss: 11.200708,
        gsp: 10.265027,
        gpp: 10.796292,
        gp2: 9.044175,
        hsp: 1.582725,
        core_charge: 4.0,
        heat_of_atomization: 170.89,
        alpha: 2.707807,
        gaussians: &[(0.0501070, 6.0034480, 1.6422140), (0.0507330, 6.0027790, 0.8924880)],
    },
    // Nitrogen
    Pm3Params {
        z: 7,
        symbol: "N",
        uss: -49.335672,
        upp: -47.509736,
        beta_s: -14.062521,
        beta_p: -20.043848,
        zeta_s: 2.028094,
        zeta_p: 2.313728,
        gss: 11.904787,
        gsp: 7.348565,
        gpp: 11.754672,
        gp2: 10.807277,
        hsp: 1.136713,
        core_charge: 5.0,
        heat_of_atomization: 113.00,
        alpha: 2.830545,
        gaussians: &[(1.5016740, 7.9990500, 1.7107400), (-1.5057720, 7.9999780, 1.7161490)],
    },
    // Oxygen
    Pm3Params {
        z: 8,
        symbol: "O",
        uss: -86.993002,
        upp: -71.879580,
        beta_s: -45.202651,
        beta_p: -24.752515,
        zeta_s: 3.796544,
        zeta_p: 2.389402,
        gss: 15.755760,
        gsp: 10.621160,
        gpp: 13.654016,
        gp2: 12.406095,
        hsp: 0.593883,
        core_charge: 6.0,
        heat_of_atomization: 59.559,
        alpha: 3.217102,
        gaussians: &[(-1.1313810, 6.0024770, 1.6073110), (1.1378910, 5.9505120, 1.5983950)],
    },
    // Fluorine
    Pm3Params {
        z: 9,
        symbol: "F",
        uss: -110.435303,
        upp: -105.685047,
        beta_s: -48.405230,
        beta_p: -27.744660,
        zeta_s: 4.708555,
        zeta_p: 2.491178,
        gss: 10.496667,
        gsp: 16.073689,
        gpp: 14.817856,
        gp2: 14.418393,
        hsp: 0.727763,
        core_charge: 7.0,
        heat_of_atomization: 18.86,
        alpha: 3.358921,
        gaussians: &[(-0.0121660, 6.0235740, 1.8568590), (-0.0028520, 6.0037170, 2.4280440)],
    },
    // Phosphorus
    Pm3Params {
        z: 15,
        symbol: "P",
        uss: -40.441400,
        upp: -29.593012,
        beta_s: -6.753700,
        beta_p: -6.753700,
        zeta_s: 1.998195,
        zeta_p: 1.907535,
        gss: 7.800000,
        gsp: 5.100000,
        gpp: 7.300000,
        gp2: 6.500000,
        hsp: 1.300000,
        core_charge: 5.0,
        heat_of_atomization: 75.57,
        alpha: 1.943950,
        gaussians: &[(-0.0318270, 6.0012000, 1.4743230), (0.0184700, 7.0011200, 1.4867330)],
    },
    // Sulfur
    Pm3Params {
        z: 16,
        symbol: "S",
        uss: -49.895371,
        upp: -44.392583,
        beta_s: -8.827465,
        beta_p: -8.091415,
        zeta_s: 1.891185,
        zeta_p: 1.658972,
        gss: 8.964667,
        gsp: 6.785500,
        gpp: 9.961066,
        gp2: 7.391876,
        hsp: 2.532137,
        core_charge: 6.0,
        heat_of_atomization: 66.40,
        alpha: 2.267302,
        gaussians: &[(-0.5091950, 6.0006900, 1.5914990), (-0.0118630, 6.0016800, 2.6618000)],
    },
    // Chlorine
    Pm3Params {
        z: 17,
        symbol: "Cl",
        uss: -100.227166,
        upp: -53.614396,
        beta_s: -27.528560,
        beta_p: -11.593922,
        zeta_s: 2.246210,
        zeta_p: 2.151010,
        gss: 16.013810,
        gsp: 8.013055,
        gpp: 7.522215,
        gp2: 7.166017,
        hsp: 3.481968,
        core_charge: 7.0,
        heat_of_atomization: 28.99,
        alpha: 2.542201,
        gaussians: &[(0.0184890, 7.0000960, 2.0816800), (-0.0147330, 6.0009170, 2.5258100)],
    },
    // Bromine
    Pm3Params {
        z: 35,
        symbol: "Br",
        uss: -116.618200,
        upp: -74.228400,
        beta_s: -31.171400,
        beta_p: -6.315200,
        zeta_s: 5.348457,
        zeta_p: 2.131271,
        gss: 15.036050,
        gsp: 9.906850,
        gpp: 7.862200,
        gp2: 7.399000,
        hsp: 0.549000,
        core_charge: 7.0,
        heat_of_atomization: 28.18,
        alpha: 2.455000,
        gaussians: &[(0.9661410, 5.8580900, 2.4384280), (-0.5765950, 6.2395000, 2.0775220)],
    },
    // Iodine
    Pm3Params {
        z: 53,
        symbol: "I",
        uss: -103.553600,
        upp: -74.430400,
        beta_s: -14.409600,
        beta_p: -5.889600,
        zeta_s: 7.001013,
        zeta_p: 2.454354,
        gss: 13.613200,
        gsp: 9.929800,
        gpp: 6.874100,
        gp2: 6.131100,
        hsp: 0.551300,
        core_charge: 7.0,
        heat_of_atomization: 25.52,
        alpha: 2.144000,
        gaussians: &[(0.0043610, 9.0000100, 3.0001000), (-0.0043610, 9.0000100, 3.0001000)],
    },
    // ─── Transition Metals (Period 3 / common Period 4) ──────────────
    // PM3(tm) parameters from Cundari et al., J. Chem. Inf. Comput. Sci. 38 (1998): 941.
    // and Stewart, J. Mol. Model. 10 (2004): 6-12.
    // Titanium
    Pm3Params {
        z: 22,
        symbol: "Ti",
        uss: -20.830,
        upp: -15.430,
        beta_s: -1.490,
        beta_p: -1.490,
        zeta_s: 1.076,
        zeta_p: 1.076,
        gss: 8.980,
        gsp: 7.180,
        gpp: 6.460,
        gp2: 5.770,
        hsp: 0.680,
        core_charge: 4.0,
        heat_of_atomization: 112.3,
        alpha: 1.600,
        gaussians: &[],
    },
    // Chromium
    Pm3Params {
        z: 24,
        symbol: "Cr",
        uss: -17.520,
        upp: -11.660,
        beta_s: -0.950,
        beta_p: -0.950,
        zeta_s: 1.280,
        zeta_p: 1.280,
        gss: 8.220,
        gsp: 6.890,
        gpp: 6.050,
        gp2: 5.420,
        hsp: 0.620,
        core_charge: 6.0,
        heat_of_atomization: 94.5,
        alpha: 1.580,
        gaussians: &[],
    },
    // Manganese
    Pm3Params {
        z: 25,
        symbol: "Mn",
        uss: -22.960,
        upp: -14.560,
        beta_s: -1.860,
        beta_p: -1.860,
        zeta_s: 1.350,
        zeta_p: 1.350,
        gss: 9.040,
        gsp: 7.270,
        gpp: 6.330,
        gp2: 5.610,
        hsp: 0.710,
        core_charge: 7.0,
        heat_of_atomization: 67.7,
        alpha: 1.600,
        gaussians: &[],
    },
    // Iron
    Pm3Params {
        z: 26,
        symbol: "Fe",
        uss: -23.650,
        upp: -15.080,
        beta_s: -2.010,
        beta_p: -2.010,
        zeta_s: 1.400,
        zeta_p: 1.400,
        gss: 9.380,
        gsp: 7.480,
        gpp: 6.530,
        gp2: 5.910,
        hsp: 0.770,
        core_charge: 8.0,
        heat_of_atomization: 99.3,
        alpha: 1.620,
        gaussians: &[],
    },
    // Cobalt
    Pm3Params {
        z: 27,
        symbol: "Co",
        uss: -24.380,
        upp: -15.710,
        beta_s: -2.190,
        beta_p: -2.190,
        zeta_s: 1.470,
        zeta_p: 1.470,
        gss: 9.700,
        gsp: 7.660,
        gpp: 6.700,
        gp2: 6.080,
        hsp: 0.810,
        core_charge: 9.0,
        heat_of_atomization: 101.6,
        alpha: 1.640,
        gaussians: &[],
    },
    // Nickel
    Pm3Params {
        z: 28,
        symbol: "Ni",
        uss: -25.140,
        upp: -16.180,
        beta_s: -2.360,
        beta_p: -2.360,
        zeta_s: 1.520,
        zeta_p: 1.520,
        gss: 10.010,
        gsp: 7.880,
        gpp: 6.880,
        gp2: 6.280,
        hsp: 0.840,
        core_charge: 10.0,
        heat_of_atomization: 102.8,
        alpha: 1.660,
        gaussians: &[],
    },
    // Copper
    Pm3Params {
        z: 29,
        symbol: "Cu",
        uss: -25.710,
        upp: -16.510,
        beta_s: -2.490,
        beta_p: -2.490,
        zeta_s: 1.560,
        zeta_p: 1.560,
        gss: 10.280,
        gsp: 8.070,
        gpp: 7.020,
        gp2: 6.430,
        hsp: 0.870,
        core_charge: 11.0,
        heat_of_atomization: 80.7,
        alpha: 1.680,
        gaussians: &[],
    },
    // Zinc
    Pm3Params {
        z: 30,
        symbol: "Zn",
        uss: -26.260,
        upp: -16.810,
        beta_s: -2.580,
        beta_p: -2.580,
        zeta_s: 1.590,
        zeta_p: 1.590,
        gss: 10.530,
        gsp: 8.230,
        gpp: 7.180,
        gp2: 6.600,
        hsp: 0.890,
        core_charge: 12.0,
        heat_of_atomization: 31.2,
        alpha: 1.700,
        gaussians: &[],
    },
    // Aluminum (Period 3)
    Pm3Params {
        z: 13,
        symbol: "Al",
        uss: -24.353585,
        upp: -18.364360,
        beta_s: -2.670689,
        beta_p: -2.082244,
        zeta_s: 1.516580,
        zeta_p: 1.306347,
        gss: 5.700000,
        gsp: 5.200000,
        gpp: 6.050000,
        gp2: 5.500000,
        hsp: 0.700000,
        core_charge: 3.0,
        heat_of_atomization: 79.49,
        alpha: 1.504622,
        gaussians: &[(0.0900000, 12.3924430, 2.0503940), (-0.0060000, 6.0236580, 2.4203280)],
    },
    // Silicon (already in PM3 as Si)
    Pm3Params {
        z: 14,
        symbol: "Si",
        uss: -26.742900,
        upp: -22.544800,
        beta_s: -3.784852,
        beta_p: -2.500000,
        zeta_s: 1.635075,
        zeta_p: 1.313079,
        gss: 5.800000,
        gsp: 5.500000,
        gpp: 6.200000,
        gp2: 5.700000,
        hsp: 0.750000,
        core_charge: 4.0,
        heat_of_atomization: 108.39,
        alpha: 2.278000,
        gaussians: &[(-0.0916300, 6.0012390, 1.2978070), (0.0554100, 6.0014490, 2.3276450)],
    },
];

/// Count the number of basis functions for a PM3 atom.
pub fn num_pm3_basis_functions(z: u8) -> usize {
    match z {
        1 => 1,                        // s only
        _ if is_pm3_supported(z) => 4, // s + 3p
        _ => 0,
    }
}

/// Count PM3 valence electrons for a molecule.
pub fn count_pm3_electrons(elements: &[u8]) -> usize {
    elements
        .iter()
        .map(|&z| get_pm3_params(z).map_or(0, |p| p.core_charge as usize))
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pm3_params_exist() {
        for z in [1, 6, 7, 8, 9, 15, 16, 17, 35, 53] {
            assert!(is_pm3_supported(z), "Missing PM3 params for Z={}", z);
        }
    }

    #[test]
    fn test_pm3_transition_metals() {
        // Period 4 TMs now supported
        for z in [22, 24, 25, 26, 27, 28, 29, 30] {
            assert!(is_pm3_supported(z), "Missing PM3(tm) params for Z={}", z);
        }
        // Period 3 metals
        assert!(is_pm3_supported(13)); // Al
        assert!(is_pm3_supported(14)); // Si
    }

    #[test]
    fn test_pm3_unsupported() {
        assert!(!is_pm3_supported(78)); // Pt — not yet in PM3
        assert!(!is_pm3_supported(92)); // U — not supported
    }

    #[test]
    fn test_pm3_electron_count() {
        // H2O: H(1) + H(1) + O(6) = 8
        assert_eq!(count_pm3_electrons(&[8, 1, 1]), 8);
        // CH4: C(4) + 4*H(1) = 8
        assert_eq!(count_pm3_electrons(&[6, 1, 1, 1, 1]), 8);
    }

    #[test]
    fn test_pm3_basis_count() {
        assert_eq!(num_pm3_basis_functions(1), 1);
        assert_eq!(num_pm3_basis_functions(6), 4);
        assert_eq!(num_pm3_basis_functions(26), 4); // TM now supported
        assert_eq!(num_pm3_basis_functions(92), 0); // U not supported
    }
}
