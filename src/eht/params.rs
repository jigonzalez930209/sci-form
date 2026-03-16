//! EHT parameter tables: VSIP values and Slater exponents per element/orbital.
//!
//! Hoffmann parameters for valence orbitals of common elements.

use serde::{Deserialize, Serialize};

/// Confidence level for the current EHT parameterization.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SupportLevel {
    Unsupported,
    Experimental,
    Supported,
}

/// Summary of EHT support for a specific element set.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EhtSupport {
    /// Overall support level across the full molecule.
    pub level: SupportLevel,
    /// Whether the element set contains transition metals.
    pub has_transition_metals: bool,
    /// Supported, validated elements for the current EHT implementation.
    pub supported_elements: Vec<u8>,
    /// Supported but still provisional / uncalibrated elements.
    pub provisional_elements: Vec<u8>,
    /// Elements with no EHT parameter set.
    pub unsupported_elements: Vec<u8>,
    /// Human-readable warnings for callers that want to surface caveats.
    pub warnings: Vec<String>,
}

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

/// Return true if the element belongs to the d-block series currently covered by EHT.
pub fn is_transition_metal(z: u8) -> bool {
    matches!(z, 21..=30 | 39..=48 | 72..=80)
}

/// Return the support level for a single element.
pub fn support_level_for_element(z: u8) -> SupportLevel {
    if get_params(z).is_none() {
        SupportLevel::Unsupported
    } else if is_transition_metal(z) {
        // Calibrated against Alvarez/Ammeter/Hoffmann literature tables.
        SupportLevel::Experimental
    } else {
        SupportLevel::Supported
    }
}

/// Analyze a molecule-level EHT support profile from atomic numbers.
pub fn analyze_eht_support(elements: &[u8]) -> EhtSupport {
    let mut supported_elements = Vec::new();
    let mut provisional_elements = Vec::new();
    let mut unsupported_elements = Vec::new();

    for &z in elements {
        match support_level_for_element(z) {
            SupportLevel::Supported => {
                if !supported_elements.contains(&z) {
                    supported_elements.push(z);
                }
            }
            SupportLevel::Experimental => {
                if !provisional_elements.contains(&z) {
                    provisional_elements.push(z);
                }
            }
            SupportLevel::Unsupported => {
                if !unsupported_elements.contains(&z) {
                    unsupported_elements.push(z);
                }
            }
        }
    }

    supported_elements.sort_unstable();
    provisional_elements.sort_unstable();
    unsupported_elements.sort_unstable();

    let mut warnings = Vec::new();
    if !provisional_elements.is_empty() {
        warnings.push(format!(
            "Transition-metal EHT parameters are provisional for elements {:?}; results should be treated as experimental until benchmark calibration is completed.",
            provisional_elements
        ));
    }
    if !unsupported_elements.is_empty() {
        warnings.push(format!(
            "No EHT parameters are available for elements {:?}.",
            unsupported_elements
        ));
    }

    let level = if !unsupported_elements.is_empty() {
        SupportLevel::Unsupported
    } else if !provisional_elements.is_empty() {
        SupportLevel::Experimental
    } else {
        SupportLevel::Supported
    };

    EhtSupport {
        level,
        has_transition_metals: !provisional_elements.is_empty(),
        supported_elements,
        provisional_elements,
        unsupported_elements,
        warnings,
    }
}

// ─── Parameter tables ────────────────────────────────────────────────────────

static H_ORBITALS: [OrbitalDef; 1] = [OrbitalDef {
    n: 1,
    l: 0,
    label: "1s",
    vsip: -13.6,
    zeta: 1.3,
}];

static C_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef {
        n: 2,
        l: 0,
        label: "2s",
        vsip: -21.4,
        zeta: 1.625,
    },
    OrbitalDef {
        n: 2,
        l: 1,
        label: "2p",
        vsip: -11.4,
        zeta: 1.625,
    },
];

static N_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef {
        n: 2,
        l: 0,
        label: "2s",
        vsip: -26.0,
        zeta: 1.950,
    },
    OrbitalDef {
        n: 2,
        l: 1,
        label: "2p",
        vsip: -13.4,
        zeta: 1.950,
    },
];

static O_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef {
        n: 2,
        l: 0,
        label: "2s",
        vsip: -32.3,
        zeta: 2.275,
    },
    OrbitalDef {
        n: 2,
        l: 1,
        label: "2p",
        vsip: -14.8,
        zeta: 2.275,
    },
];

static F_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef {
        n: 2,
        l: 0,
        label: "2s",
        vsip: -40.0,
        zeta: 2.425,
    },
    OrbitalDef {
        n: 2,
        l: 1,
        label: "2p",
        vsip: -18.1,
        zeta: 2.425,
    },
];

static CL_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef {
        n: 3,
        l: 0,
        label: "3s",
        vsip: -26.3,
        zeta: 2.183,
    },
    OrbitalDef {
        n: 3,
        l: 1,
        label: "3p",
        vsip: -14.2,
        zeta: 1.733,
    },
];

static BR_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef {
        n: 4,
        l: 0,
        label: "4s",
        vsip: -22.07,
        zeta: 2.588,
    },
    OrbitalDef {
        n: 4,
        l: 1,
        label: "4p",
        vsip: -13.1,
        zeta: 2.131,
    },
];

static I_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef {
        n: 5,
        l: 0,
        label: "5s",
        vsip: -18.0,
        zeta: 2.679,
    },
    OrbitalDef {
        n: 5,
        l: 1,
        label: "5p",
        vsip: -12.7,
        zeta: 2.322,
    },
];

static S_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef {
        n: 3,
        l: 0,
        label: "3s",
        vsip: -20.0,
        zeta: 1.817,
    },
    OrbitalDef {
        n: 3,
        l: 1,
        label: "3p",
        vsip: -11.0,
        zeta: 1.817,
    },
];

static P_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef {
        n: 3,
        l: 0,
        label: "3s",
        vsip: -18.6,
        zeta: 1.600,
    },
    OrbitalDef {
        n: 3,
        l: 1,
        label: "3p",
        vsip: -14.0,
        zeta: 1.600,
    },
];

static SI_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef {
        n: 3,
        l: 0,
        label: "3s",
        vsip: -17.3,
        zeta: 1.383,
    },
    OrbitalDef {
        n: 3,
        l: 1,
        label: "3p",
        vsip: -9.2,
        zeta: 1.383,
    },
];

static B_ORBITALS: [OrbitalDef; 2] = [
    OrbitalDef {
        n: 2,
        l: 0,
        label: "2s",
        vsip: -15.2,
        zeta: 1.300,
    },
    OrbitalDef {
        n: 2,
        l: 1,
        label: "2p",
        vsip: -8.5,
        zeta: 1.300,
    },
];

// ─── First-row transition metals (3d series): 4s + 4p + 3d ──────────────────
//
// Literature-calibrated parameters from:
//   Alvarez, S. "Tables of Parameters for Extended Hückel Calculations" (1993)
//   Ammeter, J. H. et al. Helv. Chim. Acta 61 (1978): 1.
//   Hoffmann, R. J. Chem. Phys. 39 (1963): 1397.
//
// VSIP values (eV) correspond to valence state ionization potentials.
// Slater exponents (ζ, bohr⁻¹) are single-zeta values from Alvarez tables.
// The outer np shell is included per Ammeter convention for proper σ/π metal-ligand bonding.
static SC_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 4, l: 0, label: "4s", vsip: -8.87, zeta: 1.300 },
    OrbitalDef { n: 4, l: 1, label: "4p", vsip: -2.75, zeta: 1.300 },
    OrbitalDef { n: 3, l: 2, label: "3d", vsip: -8.51, zeta: 4.350 },
];
static TI_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 4, l: 0, label: "4s", vsip: -8.97, zeta: 1.075 },
    OrbitalDef { n: 4, l: 1, label: "4p", vsip: -5.44, zeta: 1.075 },
    OrbitalDef { n: 3, l: 2, label: "3d", vsip: -10.81, zeta: 4.550 },
];
static V_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 4, l: 0, label: "4s", vsip: -8.81, zeta: 1.300 },
    OrbitalDef { n: 4, l: 1, label: "4p", vsip: -5.52, zeta: 1.300 },
    OrbitalDef { n: 3, l: 2, label: "3d", vsip: -11.00, zeta: 4.750 },
];
static CR_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 4, l: 0, label: "4s", vsip: -8.66, zeta: 1.700 },
    OrbitalDef { n: 4, l: 1, label: "4p", vsip: -5.24, zeta: 1.700 },
    OrbitalDef { n: 3, l: 2, label: "3d", vsip: -11.22, zeta: 4.950 },
];
static MN_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 4, l: 0, label: "4s", vsip: -9.75, zeta: 1.800 },
    OrbitalDef { n: 4, l: 1, label: "4p", vsip: -5.89, zeta: 1.800 },
    OrbitalDef { n: 3, l: 2, label: "3d", vsip: -11.67, zeta: 5.150 },
];
static FE_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 4, l: 0, label: "4s", vsip: -9.10, zeta: 1.900 },
    OrbitalDef { n: 4, l: 1, label: "4p", vsip: -5.32, zeta: 1.900 },
    OrbitalDef { n: 3, l: 2, label: "3d", vsip: -12.60, zeta: 5.350 },
];
static CO_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 4, l: 0, label: "4s", vsip: -9.21, zeta: 2.000 },
    OrbitalDef { n: 4, l: 1, label: "4p", vsip: -5.29, zeta: 2.000 },
    OrbitalDef { n: 3, l: 2, label: "3d", vsip: -13.18, zeta: 5.550 },
];
static NI_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 4, l: 0, label: "4s", vsip: -10.95, zeta: 2.100 },
    OrbitalDef { n: 4, l: 1, label: "4p", vsip: -6.27, zeta: 2.100 },
    OrbitalDef { n: 3, l: 2, label: "3d", vsip: -13.49, zeta: 5.750 },
];
static CU_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 4, l: 0, label: "4s", vsip: -11.40, zeta: 2.200 },
    OrbitalDef { n: 4, l: 1, label: "4p", vsip: -6.06, zeta: 2.200 },
    OrbitalDef { n: 3, l: 2, label: "3d", vsip: -14.00, zeta: 5.950 },
];
static ZN_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 4, l: 0, label: "4s", vsip: -12.41, zeta: 2.010 },
    OrbitalDef { n: 4, l: 1, label: "4p", vsip: -6.53, zeta: 2.010 },
    OrbitalDef { n: 3, l: 2, label: "3d", vsip: -17.10, zeta: 6.150 },
];

// ─── Second-row transition metals (4d series): 5s + 5p + 4d ─────────────────
//
// Alvarez, S. "Tables of Parameters for Extended Hückel Calculations" (1993).
static Y_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 5, l: 0, label: "5s", vsip: -7.29, zeta: 1.390 },
    OrbitalDef { n: 5, l: 1, label: "5p", vsip: -4.37, zeta: 1.390 },
    OrbitalDef { n: 4, l: 2, label: "4d", vsip: -8.46, zeta: 3.310 },
];
static ZR_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 5, l: 0, label: "5s", vsip: -8.12, zeta: 1.520 },
    OrbitalDef { n: 5, l: 1, label: "5p", vsip: -5.12, zeta: 1.520 },
    OrbitalDef { n: 4, l: 2, label: "4d", vsip: -10.14, zeta: 3.840 },
];
static NB_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 5, l: 0, label: "5s", vsip: -10.10, zeta: 1.640 },
    OrbitalDef { n: 5, l: 1, label: "5p", vsip: -6.86, zeta: 1.640 },
    OrbitalDef { n: 4, l: 2, label: "4d", vsip: -12.10, zeta: 4.080 },
];
static MO_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 5, l: 0, label: "5s", vsip: -8.34, zeta: 1.730 },
    OrbitalDef { n: 5, l: 1, label: "5p", vsip: -5.24, zeta: 1.730 },
    OrbitalDef { n: 4, l: 2, label: "4d", vsip: -10.50, zeta: 4.540 },
];
static TC_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 5, l: 0, label: "5s", vsip: -9.00, zeta: 1.820 },
    OrbitalDef { n: 5, l: 1, label: "5p", vsip: -5.60, zeta: 1.820 },
    OrbitalDef { n: 4, l: 2, label: "4d", vsip: -11.20, zeta: 4.900 },
];
static RU_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 5, l: 0, label: "5s", vsip: -10.40, zeta: 1.900 },
    OrbitalDef { n: 5, l: 1, label: "5p", vsip: -6.87, zeta: 1.900 },
    OrbitalDef { n: 4, l: 2, label: "4d", vsip: -14.90, zeta: 5.380 },
];
static RH_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 5, l: 0, label: "5s", vsip: -8.09, zeta: 2.135 },
    OrbitalDef { n: 5, l: 1, label: "5p", vsip: -4.57, zeta: 2.135 },
    OrbitalDef { n: 4, l: 2, label: "4d", vsip: -12.50, zeta: 4.290 },
];
static PD_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 5, l: 0, label: "5s", vsip: -7.32, zeta: 2.190 },
    OrbitalDef { n: 5, l: 1, label: "5p", vsip: -3.75, zeta: 2.190 },
    OrbitalDef { n: 4, l: 2, label: "4d", vsip: -12.02, zeta: 5.983 },
];
static AG_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 5, l: 0, label: "5s", vsip: -6.27, zeta: 2.242 },
    OrbitalDef { n: 5, l: 1, label: "5p", vsip: -3.97, zeta: 2.242 },
    OrbitalDef { n: 4, l: 2, label: "4d", vsip: -14.58, zeta: 6.070 },
];
static CD_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 5, l: 0, label: "5s", vsip: -11.79, zeta: 2.300 },
    OrbitalDef { n: 5, l: 1, label: "5p", vsip: -6.10, zeta: 2.300 },
    OrbitalDef { n: 4, l: 2, label: "4d", vsip: -17.84, zeta: 6.330 },
];

// ─── Third-row transition metals (5d series): 6s + 6p + 5d ──────────────────
//
// Alvarez, S. "Tables of Parameters for Extended Hückel Calculations" (1993).
// Pt values from Thorn, D. L.; Hoffmann, R. Inorg. Chem. 17 (1978): 126.
static HF_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 6, l: 0, label: "6s", vsip: -8.20, zeta: 1.720 },
    OrbitalDef { n: 6, l: 1, label: "6p", vsip: -4.65, zeta: 1.720 },
    OrbitalDef { n: 5, l: 2, label: "5d", vsip: -11.18, zeta: 4.360 },
];
static TA_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 6, l: 0, label: "6s", vsip: -10.79, zeta: 1.830 },
    OrbitalDef { n: 6, l: 1, label: "6p", vsip: -6.86, zeta: 1.830 },
    OrbitalDef { n: 5, l: 2, label: "5d", vsip: -12.10, zeta: 4.762 },
];
static W_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 6, l: 0, label: "6s", vsip: -8.26, zeta: 1.890 },
    OrbitalDef { n: 6, l: 1, label: "6p", vsip: -5.17, zeta: 1.890 },
    OrbitalDef { n: 5, l: 2, label: "5d", vsip: -10.37, zeta: 4.982 },
];
static RE_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 6, l: 0, label: "6s", vsip: -9.36, zeta: 1.980 },
    OrbitalDef { n: 6, l: 1, label: "6p", vsip: -5.96, zeta: 1.980 },
    OrbitalDef { n: 5, l: 2, label: "5d", vsip: -12.66, zeta: 5.343 },
];
static OS_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 6, l: 0, label: "6s", vsip: -8.17, zeta: 2.070 },
    OrbitalDef { n: 6, l: 1, label: "6p", vsip: -4.81, zeta: 2.070 },
    OrbitalDef { n: 5, l: 2, label: "5d", vsip: -11.84, zeta: 5.571 },
];
static IR_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 6, l: 0, label: "6s", vsip: -11.36, zeta: 2.200 },
    OrbitalDef { n: 6, l: 1, label: "6p", vsip: -4.50, zeta: 2.200 },
    OrbitalDef { n: 5, l: 2, label: "5d", vsip: -12.17, zeta: 5.796 },
];
static PT_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 6, l: 0, label: "6s", vsip: -9.077, zeta: 2.554 },
    OrbitalDef { n: 6, l: 1, label: "6p", vsip: -5.475, zeta: 2.554 },
    OrbitalDef { n: 5, l: 2, label: "5d", vsip: -12.59, zeta: 6.013 },
];
static AU_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 6, l: 0, label: "6s", vsip: -10.92, zeta: 2.602 },
    OrbitalDef { n: 6, l: 1, label: "6p", vsip: -5.55, zeta: 2.602 },
    OrbitalDef { n: 5, l: 2, label: "5d", vsip: -15.07, zeta: 6.163 },
];
static HG_ORBITALS: [OrbitalDef; 3] = [
    OrbitalDef { n: 6, l: 0, label: "6s", vsip: -13.68, zeta: 2.649 },
    OrbitalDef { n: 6, l: 1, label: "6p", vsip: -8.47, zeta: 2.649 },
    OrbitalDef { n: 5, l: 2, label: "5d", vsip: -17.50, zeta: 6.350 },
];

/// All supported element parameter sets.
static ALL_PARAMS: &[EhtParams] = &[
    EhtParams {
        z: 1,
        symbol: "H",
        orbitals: &H_ORBITALS,
    },
    EhtParams {
        z: 5,
        symbol: "B",
        orbitals: &B_ORBITALS,
    },
    EhtParams {
        z: 6,
        symbol: "C",
        orbitals: &C_ORBITALS,
    },
    EhtParams {
        z: 7,
        symbol: "N",
        orbitals: &N_ORBITALS,
    },
    EhtParams {
        z: 8,
        symbol: "O",
        orbitals: &O_ORBITALS,
    },
    EhtParams {
        z: 9,
        symbol: "F",
        orbitals: &F_ORBITALS,
    },
    EhtParams {
        z: 14,
        symbol: "Si",
        orbitals: &SI_ORBITALS,
    },
    EhtParams {
        z: 15,
        symbol: "P",
        orbitals: &P_ORBITALS,
    },
    EhtParams {
        z: 16,
        symbol: "S",
        orbitals: &S_ORBITALS,
    },
    EhtParams {
        z: 17,
        symbol: "Cl",
        orbitals: &CL_ORBITALS,
    },
    EhtParams {
        z: 21,
        symbol: "Sc",
        orbitals: &SC_ORBITALS,
    },
    EhtParams {
        z: 22,
        symbol: "Ti",
        orbitals: &TI_ORBITALS,
    },
    EhtParams {
        z: 23,
        symbol: "V",
        orbitals: &V_ORBITALS,
    },
    EhtParams {
        z: 24,
        symbol: "Cr",
        orbitals: &CR_ORBITALS,
    },
    EhtParams {
        z: 25,
        symbol: "Mn",
        orbitals: &MN_ORBITALS,
    },
    EhtParams {
        z: 26,
        symbol: "Fe",
        orbitals: &FE_ORBITALS,
    },
    EhtParams {
        z: 27,
        symbol: "Co",
        orbitals: &CO_ORBITALS,
    },
    EhtParams {
        z: 28,
        symbol: "Ni",
        orbitals: &NI_ORBITALS,
    },
    EhtParams {
        z: 29,
        symbol: "Cu",
        orbitals: &CU_ORBITALS,
    },
    EhtParams {
        z: 30,
        symbol: "Zn",
        orbitals: &ZN_ORBITALS,
    },
    EhtParams {
        z: 35,
        symbol: "Br",
        orbitals: &BR_ORBITALS,
    },
    EhtParams {
        z: 39,
        symbol: "Y",
        orbitals: &Y_ORBITALS,
    },
    EhtParams {
        z: 40,
        symbol: "Zr",
        orbitals: &ZR_ORBITALS,
    },
    EhtParams {
        z: 41,
        symbol: "Nb",
        orbitals: &NB_ORBITALS,
    },
    EhtParams {
        z: 42,
        symbol: "Mo",
        orbitals: &MO_ORBITALS,
    },
    EhtParams {
        z: 43,
        symbol: "Tc",
        orbitals: &TC_ORBITALS,
    },
    EhtParams {
        z: 44,
        symbol: "Ru",
        orbitals: &RU_ORBITALS,
    },
    EhtParams {
        z: 45,
        symbol: "Rh",
        orbitals: &RH_ORBITALS,
    },
    EhtParams {
        z: 46,
        symbol: "Pd",
        orbitals: &PD_ORBITALS,
    },
    EhtParams {
        z: 47,
        symbol: "Ag",
        orbitals: &AG_ORBITALS,
    },
    EhtParams {
        z: 48,
        symbol: "Cd",
        orbitals: &CD_ORBITALS,
    },
    EhtParams {
        z: 53,
        symbol: "I",
        orbitals: &I_ORBITALS,
    },
    EhtParams {
        z: 72,
        symbol: "Hf",
        orbitals: &HF_ORBITALS,
    },
    EhtParams {
        z: 73,
        symbol: "Ta",
        orbitals: &TA_ORBITALS,
    },
    EhtParams {
        z: 74,
        symbol: "W",
        orbitals: &W_ORBITALS,
    },
    EhtParams {
        z: 75,
        symbol: "Re",
        orbitals: &RE_ORBITALS,
    },
    EhtParams {
        z: 76,
        symbol: "Os",
        orbitals: &OS_ORBITALS,
    },
    EhtParams {
        z: 77,
        symbol: "Ir",
        orbitals: &IR_ORBITALS,
    },
    EhtParams {
        z: 78,
        symbol: "Pt",
        orbitals: &PT_ORBITALS,
    },
    EhtParams {
        z: 79,
        symbol: "Au",
        orbitals: &AU_ORBITALS,
    },
    EhtParams {
        z: 80,
        symbol: "Hg",
        orbitals: &HG_ORBITALS,
    },
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
        assert_eq!(num_basis_functions(1), 1); // H: 1s
        assert_eq!(num_basis_functions(6), 4); // C: 2s + 2p(3)
        assert_eq!(num_basis_functions(8), 4); // O: 2s + 2p(3)
    }

    #[test]
    fn test_transition_metals_are_experimental() {
        assert_eq!(support_level_for_element(26), SupportLevel::Experimental);
        assert!(is_transition_metal(26));
    }

    #[test]
    fn test_support_summary_collects_warnings() {
        let support = analyze_eht_support(&[8, 26, 118]);
        assert_eq!(support.level, SupportLevel::Unsupported);
        assert!(support.has_transition_metals);
        assert_eq!(support.supported_elements, vec![8]);
        assert_eq!(support.provisional_elements, vec![26]);
        assert_eq!(support.unsupported_elements, vec![118]);
        assert_eq!(support.warnings.len(), 2);
    }
}
