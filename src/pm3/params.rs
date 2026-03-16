//! PM3 parameter tables.
//!
//! Standard PM3 parameters from Stewart, J. J. P. J. Comput. Chem. 10 (1989): 209.
//! Only the most chemically useful elements are parameterized.

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
    fn test_pm3_unsupported() {
        assert!(!is_pm3_supported(26)); // Fe
        assert!(!is_pm3_supported(78)); // Pt
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
        assert_eq!(num_pm3_basis_functions(26), 0);
    }
}
