//! EHT parameter tables: VSIP values and Slater exponents per element/orbital.
//!
//! Hoffmann parameters for valence orbitals of common elements.

use serde::Serialize;

/// Definition of a single atomic orbital for EHT.
#[derive(Debug, Clone, Serialize)]
pub struct OrbitalDef {
    /// Principal quantum number n.
    pub n: u8,
    /// Angular momentum quantum number l (0=s, 1=p, 2=d).
    pub l: u8,
    /// Label, e.g. "2s", "2p".
    pub label: &'static str,
    /// Valence State Ionization Potential (eV), used as H_ii.
    pub vsip: f64,
    /// Slater exponent ζ (bohr⁻¹).
    pub zeta: f64,
}

/// Full EHT parameter set for one element.
#[derive(Debug, Clone, Serialize)]
pub struct EhtParams {
    /// Atomic number.
    pub z: u8,
    /// Element symbol.
    pub symbol: &'static str,
    /// Valence orbital definitions.
    pub orbitals: &'static [OrbitalDef],
}

/// Wolfsberg-Helmholtz constant (dimensionless).
pub const DEFAULT_K: f64 = 1.75;

// ─── Parameter tables ────────────────────────────────────────────────────────

static H_ORBITALS: [OrbitalDef; 1] = [OrbitalDef {
    n: 1,
    l: 0,
    label: "1s",
    vsip: -13.6,
    zeta: 1.3,
}];

static C_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef { n: 2, l: 0, label: "2s", vsip: -21.4, zeta: 1.625 },
    OrbitalDef { n: 2, l: 1, label: "2p", vsip: -11.4, zeta: 1.625 },
];

static N_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef { n: 2, l: 0, label: "2s", vsip: -26.0, zeta: 1.950 },
    OrbitalDef { n: 2, l: 1, label: "2p", vsip: -13.4, zeta: 1.950 },
];

static O_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef { n: 2, l: 0, label: "2s", vsip: -32.3, zeta: 2.275 },
    OrbitalDef { n: 2, l: 1, label: "2p", vsip: -14.8, zeta: 2.275 },
];

static F_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef { n: 2, l: 0, label: "2s", vsip: -40.0, zeta: 2.425 },
    OrbitalDef { n: 2, l: 1, label: "2p", vsip: -18.1, zeta: 2.425 },
];

static CL_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef { n: 3, l: 0, label: "3s", vsip: -26.3, zeta: 2.183 },
    OrbitalDef { n: 3, l: 1, label: "3p", vsip: -14.2, zeta: 1.733 },
];

static BR_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef { n: 4, l: 0, label: "4s", vsip: -22.07, zeta: 2.588 },
    OrbitalDef { n: 4, l: 1, label: "4p", vsip: -13.1, zeta: 2.131 },
];

static I_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef { n: 5, l: 0, label: "5s", vsip: -18.0, zeta: 2.679 },
    OrbitalDef { n: 5, l: 1, label: "5p", vsip: -12.7, zeta: 2.322 },
];

static S_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef { n: 3, l: 0, label: "3s", vsip: -20.0, zeta: 1.817 },
    OrbitalDef { n: 3, l: 1, label: "3p", vsip: -11.0, zeta: 1.817 },
];

static P_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef { n: 3, l: 0, label: "3s", vsip: -18.6, zeta: 1.600 },
    OrbitalDef { n: 3, l: 1, label: "3p", vsip: -14.0, zeta: 1.600 },
];

static SI_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef { n: 3, l: 0, label: "3s", vsip: -17.3, zeta: 1.383 },
    OrbitalDef { n: 3, l: 1, label: "3p", vsip: -9.2, zeta: 1.383 },
];

static B_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef { n: 2, l: 0, label: "2s", vsip: -15.2, zeta: 1.300 },
    OrbitalDef { n: 2, l: 1, label: "2p", vsip: -8.5, zeta: 1.300 },
];

/// All supported element parameter sets.
static ALL_PARAMS: &[EhtParams] = &[
    EhtParams { z: 1,  symbol: "H",  orbitals: &H_ORBITALS },
    EhtParams { z: 5,  symbol: "B",  orbitals: &B_ORBITALS },
    EhtParams { z: 6,  symbol: "C",  orbitals: &C_ORBITALS },
    EhtParams { z: 7,  symbol: "N",  orbitals: &N_ORBITALS },
    EhtParams { z: 8,  symbol: "O",  orbitals: &O_ORBITALS },
    EhtParams { z: 9,  symbol: "F",  orbitals: &F_ORBITALS },
    EhtParams { z: 14, symbol: "Si", orbitals: &SI_ORBITALS },
    EhtParams { z: 15, symbol: "P",  orbitals: &P_ORBITALS },
    EhtParams { z: 16, symbol: "S",  orbitals: &S_ORBITALS },
    EhtParams { z: 17, symbol: "Cl", orbitals: &CL_ORBITALS },
    EhtParams { z: 35, symbol: "Br", orbitals: &BR_ORBITALS },
    EhtParams { z: 53, symbol: "I",  orbitals: &I_ORBITALS },
];

/// Look up EHT parameters by atomic number.
pub fn get_params(z: u8) -> Option<&'static EhtParams> {
    ALL_PARAMS.iter().find(|p| p.z == z)
}

/// Count the total number of valence basis functions for an element.
/// s → 1, p → 3, d → 5.
pub fn num_basis_functions(z: u8) -> usize {
    match get_params(z) {
        Some(p) => p
            .orbitals
            .iter()
            .map(|o| match o.l {
                0 => 1,
                1 => 3,
                2 => 5,
                _ => 0,
            })
            .sum(),
        None => 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hydrogen_params() {
        let p = get_params(1).expect("H params");
        assert_eq!(p.symbol, "H");
        assert_eq!(p.orbitals.len(), 1);
        assert!((p.orbitals[0].vsip - (-13.6)).abs() < 1e-10);
        assert!((p.orbitals[0].zeta - 1.3).abs() < 1e-10);
    }

    #[test]
    fn test_carbon_params() {
        let p = get_params(6).expect("C params");
        assert_eq!(p.symbol, "C");
        assert_eq!(p.orbitals.len(), 2);
        // 2s + 2px,2py,2pz = 4 basis functions
        assert_eq!(num_basis_functions(6), 4);
    }

    #[test]
    fn test_oxygen_params() {
        let p = get_params(8).expect("O params");
        assert_eq!(p.symbol, "O");
        assert!((p.orbitals[0].vsip - (-32.3)).abs() < 1e-10);
        assert!((p.orbitals[1].vsip - (-14.8)).abs() < 1e-10);
    }

    #[test]
    fn test_all_elements_have_params() {
        for z in [1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53] {
            assert!(get_params(z).is_some(), "Missing params for Z={}", z);
        }
    }

    #[test]
    fn test_unknown_element_returns_none() {
        assert!(get_params(118).is_none());
    }

    #[test]
    fn test_basis_function_count() {
        assert_eq!(num_basis_functions(1), 1);  // H: 1s
        assert_eq!(num_basis_functions(6), 4);  // C: 2s + 2p(3)
        assert_eq!(num_basis_functions(8), 4);  // O: 2s + 2p(3)
    }
}
