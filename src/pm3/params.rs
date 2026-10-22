//! PM3 parameter tables.
//!
//! Standard PM3 parameters from Stewart, J. J. P. J. Comput. Chem. 10 (1989): 209.
//! Only the most chemically useful elements are parameterized.

#![allow(clippy::approx_constant)]

use serde::Serialize;

use super::params_table::PM3_MAIN_GROUP;

/// PM3 atomic parameters for one element.
#[derive(Debug, Clone, Serialize)]
pub struct Pm3Params {
    /// Atomic number.
    pub z: u8,
    /// Element symbol.
    pub symbol: &'static str,
    /// Whether the entry comes from the PM3(tm) transition-metal
    /// parameterization (Cundari et al.) instead of Stewart's PM3.
    pub tm: bool,
    /// Principal quantum number of the valence sp shell (MOPAC `nspqn`).
    pub nspqn: u8,
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
    /// Atomic heat of formation (kcal/mol, MOPAC `eheat`): the reference that
    /// turns the SCF energy into a heat of formation.
    pub heat_of_atomization: f64,
    /// Total energy of the isolated atom (eV, MOPAC `eisol`).
    pub eisol: f64,
    /// Alpha parameter for core-core repulsion (Å⁻¹).
    pub alpha: f64,
    /// PM3 Gaussian correction terms for core-core repulsion.
    /// Each tuple is (a_k, b_k, c_k) in MOPAC's expression:
    ///   E_nuc += Z_A Z_B / R × Σ_k a_k × exp(−b_k × (R − c_k)²)
    pub gaussians: &'static [(f64, f64, f64)],
}

/// Check if an element has PM3 parameters.
pub fn is_pm3_supported(z: u8) -> bool {
    get_pm3_params(z).is_some()
}

/// Look up PM3 parameters by atomic number.
pub fn get_pm3_params(z: u8) -> Option<&'static Pm3Params> {
    PM3_MAIN_GROUP
        .iter()
        .find(|p| p.z == z)
        .or_else(|| PM3_TM.iter().find(|p| p.z == z))
}

/// Which PM3 parameterization an element uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum Pm3Method {
    /// Stewart's PM3 parameterization (main-group elements).
    Pm3,
    /// PM3(tm) transition-metal parameterization (Cundari et al., 1998).
    Pm3Tm,
}

/// Parameterization used for `z`, or `None` when the element is unsupported.
///
/// A molecule mixing main-group and transition-metal elements is evaluated
/// with both sets at once: the NDDO core is shared, but the transition metals
/// take their integrals from PM3(tm).
pub fn pm3_method(z: u8) -> Option<Pm3Method> {
    get_pm3_params(z).map(|p| {
        if p.tm {
            Pm3Method::Pm3Tm
        } else {
            Pm3Method::Pm3
        }
    })
}

/// PM3(tm) transition-metal entries (Cundari et al., J. Chem. Inf. Comput.
/// Sci. 38 (1998): 941), preserved from the original table: MOPAC's PM3 file
/// does not parameterize these elements.
static PM3_TM: &[Pm3Params] = &[
    Pm3Params {
        z: 22,
        tm: true,
        nspqn: 4,
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
        eisol: -32.680000,
        heat_of_atomization: 112.3,
        alpha: 1.600,
        gaussians: &[],
    },
    Pm3Params {
        z: 24,
        tm: true,
        nspqn: 4,
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
        eisol: -17.520000,
        heat_of_atomization: 94.5,
        alpha: 1.580,
        gaussians: &[],
    },
    Pm3Params {
        z: 25,
        tm: true,
        nspqn: 4,
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
        eisol: -36.880000,
        heat_of_atomization: 67.7,
        alpha: 1.600,
        gaussians: &[],
    },
    Pm3Params {
        z: 26,
        tm: true,
        nspqn: 4,
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
        eisol: -37.920000,
        heat_of_atomization: 99.3,
        alpha: 1.620,
        gaussians: &[],
    },
    Pm3Params {
        z: 27,
        tm: true,
        nspqn: 4,
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
        eisol: -39.060000,
        heat_of_atomization: 101.6,
        alpha: 1.640,
        gaussians: &[],
    },
    Pm3Params {
        z: 28,
        tm: true,
        nspqn: 4,
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
        eisol: -40.270000,
        heat_of_atomization: 102.8,
        alpha: 1.660,
        gaussians: &[],
    },
    Pm3Params {
        z: 29,
        tm: true,
        nspqn: 4,
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
        eisol: -25.710000,
        heat_of_atomization: 80.7,
        alpha: 1.680,
        gaussians: &[],
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
    fn test_pm3_method_delimits_tm_and_main_group() {
        assert_eq!(pm3_method(6), Some(Pm3Method::Pm3));
        assert_eq!(pm3_method(17), Some(Pm3Method::Pm3));
        assert_eq!(pm3_method(26), Some(Pm3Method::Pm3Tm));
        // Zn is a main-group PM3 element in MOPAC's table.
        assert_eq!(pm3_method(30), Some(Pm3Method::Pm3));
        // Unsupported elements are reported as such.
        assert_eq!(pm3_method(21), None);
        assert_eq!(pm3_method(92), None);
        // The TM flag matches the declared method.
        for p in PM3_MAIN_GROUP.iter().chain(PM3_TM.iter()) {
            let expected = if p.tm {
                Pm3Method::Pm3Tm
            } else {
                Pm3Method::Pm3
            };
            // Main-group entries shadow PM3(tm) duplicates (e.g. Zn).
            if p.tm && PM3_MAIN_GROUP.iter().any(|m| m.z == p.z) {
                continue;
            }
            assert_eq!(pm3_method(p.z), Some(expected), "Z={}", p.z);
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
