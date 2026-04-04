//! GFN2-xTB parameters from the official parametrization.
//!
//! Values extracted from tblite (Bannwarth, Ehlert, Grimme, JCTC 2019, 15, 1652).
//! All self-energies are in eV, Slater exponents in bohr⁻¹,
//! Hubbard parameters in Hartree, repulsion parameters dimensionless.

/// Per-element GFN2-xTB parameter set.
#[derive(Debug, Clone, Copy)]
pub struct Gfn2ElementParams {
    pub z: u8,
    /// Number of shells in the basis.
    pub nshell: usize,
    /// Angular momentum per shell: 0=s, 1=p, 2=d.
    pub ang_shell: [u8; 3],
    /// Self-energy / atomic level per shell (eV).
    pub selfenergy: [f64; 3],
    /// Slater exponent per shell (bohr⁻¹).
    pub slater: [f64; 3],
    /// Reference occupation per shell.
    pub ref_occ: [f64; 3],
    /// Hubbard parameter (Hartree).
    pub hubbard: f64,
    /// Shell-resolved Hubbard scaling: 1 + shell_hubbard[l].
    pub shell_hubbard: [f64; 3],
    /// Repulsion exponent α.
    pub rep_alpha: f64,
    /// Repulsion effective nuclear charge Z_eff.
    pub rep_zeff: f64,
    /// Number of valence electrons.
    pub n_valence: u8,
    /// Principal quantum number per shell.
    pub pqn: [u8; 3],
    /// Number of Gaussian primitives per shell for STO-nG.
    pub ngauss: [u8; 3],
    /// Shell polynomial coefficients (per shell, already ×0.01 scaled).
    /// Indexed by shell index (matching ang_shell ordering).
    pub shpoly: [f64; 3],
    /// CN-dependent self-energy shift per shell (eV).
    /// SE_eff = SE - kCN * CN.
    pub kcn: [f64; 3],
    /// Third-order Hubbard derivative (already ×0.1 scaled, Hartree).
    pub gam3: f64,
    /// Atomic radius (bohr) for shell polynomial distance scaling.
    pub atomic_rad: f64,
    /// Pauling electronegativity.
    pub pauling_en: f64,
    /// Element kind (1=main group, 2=transition metal) for gam3shell lookup.
    pub kind: u8,
    /// Multipole dipole XC kernel (already scaled by 0.01).
    pub dkernel: f64,
    /// Multipole quadrupole XC kernel (already scaled by 0.01).
    pub qkernel: f64,
    /// Multipole damping radius (bohr).
    pub mp_rad: f64,
    /// Valence CN for multipole radii.
    pub mp_vcn: f64,
}

/// Conversion factor: 1 eV → Hartree.
pub const EV_TO_HARTREE: f64 = 1.0 / 27.21138505;
/// Repulsion exponent for heavy-atom pairs.
pub const REP_KEXP: f64 = 1.5;
/// Repulsion exponent for light-atom pairs (involving H or He).
pub const REP_KEXP_LIGHT: f64 = 1.0;
/// Repulsion R exponent.
pub const REP_REXP: f64 = 1.0;

/// Global GFN2-xTB Hamiltonian scaling parameters.
/// Off-diagonal K factor by angular momentum pair: kdiag[l].
pub const KDIAG: [f64; 5] = [1.85, 2.23, 2.23, 2.23, 2.23];
/// Electronegativity scaling of off-diagonal elements.
pub const ENSCALE: f64 = 2.0e-2;
/// Exponent for Slater-ratio scaling of off-diagonal elements.
pub const WEXP: f64 = 0.5;

/// Ångström to bohr conversion factor.
pub const AATOAU: f64 = 1.8897259886;

/// Third-order shell-resolved scaling factors.
/// gam3shell(kind, l): kind=1 (main group), kind=2 (TM).
/// Indexed as [s, p, d, f] = [[1.0, 1.0], [0.5, 0.5], [0.25, 0.25], [0.25, 0.25]].
pub const GAM3_SHELL: [[f64; 2]; 4] = [
    [1.0, 1.0],     // s
    [0.5, 0.5],     // p
    [0.25, 0.25],   // d
    [0.25, 0.25],   // f
];

/// Gamma function exponent (alphaj in xtb).
pub const GAMMA_EXP: f64 = 2.0;

/// Multipole damping exponents.
pub const MP_DMP3: f64 = 3.0;
pub const MP_DMP5: f64 = 4.0;
/// Multipole CN-dependent radius shift.
pub const MP_SHIFT: f64 = 1.2;
/// Multipole radius sigmoid exponent.
pub const MP_KEXP_MULTI: f64 = 4.0;
/// Multipole maximum radius (bohr).
pub const MP_RMAX: f64 = 5.0;

/// D4 dispersion parameters for GFN2-xTB.
pub const D4_S6: f64 = 1.0;
pub const D4_S8: f64 = 2.7;
pub const D4_A1: f64 = 0.52;
pub const D4_A2: f64 = 5.0;
pub const D4_S9: f64 = 5.0;

// ---- GFN2-xTB element parameters ----
// Source: tblite/tblite gfn2.f90 (official GFN2-xTB parametrization)

static GFN2_PARAMS: &[Gfn2ElementParams] = &[
    // Z=1, H
    Gfn2ElementParams {
        z: 1,
        nshell: 1,
        ang_shell: [0, 0, 0],
        selfenergy: [-10.707211, 0.0, 0.0],
        slater: [1.230000, 0.0, 0.0],
        ref_occ: [1.0, 0.0, 0.0],
        hubbard: 0.405771,
        shell_hubbard: [1.0, 1.0, 1.0],
        rep_alpha: 2.213717,
        rep_zeff: 1.105388,
        n_valence: 1,
        pqn: [1, 0, 0],
        ngauss: [3, 0, 0],
        shpoly: [-0.00953618, 0.0, 0.0],
        kcn: [-0.0500000, 0.0, 0.0],
        gam3: 0.0800000,
        atomic_rad: 0.32 * AATOAU,
        pauling_en: 2.20,
        kind: 1,
        dkernel: 0.05563889,
        qkernel: 0.00027431,
        mp_rad: 1.4,
        mp_vcn: 1.0,
    },
    // Z=5, B
    Gfn2ElementParams {
        z: 5,
        nshell: 2,
        ang_shell: [0, 1, 0],
        selfenergy: [-9.224376, -7.419002, 0.0],
        slater: [1.479444, 1.479805, 0.0],
        ref_occ: [2.0, 1.0, 0.0],
        hubbard: 0.513556,
        shell_hubbard: [1.0, 1.3994080, 1.0],
        rep_alpha: 1.373856,
        rep_zeff: 7.192431,
        n_valence: 3,
        pqn: [2, 2, 0],
        ngauss: [4, 4, 0],
        shpoly: [-0.05183150, -0.02453322, 0.0],
        kcn: [0.0120462, -0.0141086, 0.0],
        gam3: 0.0946104,
        atomic_rad: 0.84 * AATOAU,
        pauling_en: 2.04,
        kind: 1,
        dkernel: -0.00481186,
        qkernel: -0.00058228,
        mp_rad: 5.0,
        mp_vcn: 3.0,
    },
    // Z=6, C
    Gfn2ElementParams {
        z: 6,
        nshell: 2,
        ang_shell: [0, 1, 0],
        selfenergy: [-13.970922, -10.063292, 0.0],
        slater: [2.096432, 1.800000, 0.0],
        ref_occ: [1.0, 3.0, 0.0],
        hubbard: 0.538015,
        shell_hubbard: [1.0, 1.1056358, 1.0],
        rep_alpha: 1.247655,
        rep_zeff: 4.231078,
        n_valence: 4,
        pqn: [2, 2, 0],
        ngauss: [4, 4, 0],
        shpoly: [-0.02294321, -0.00271102, 0.0],
        kcn: [-0.0102144, 0.0161657, 0.0],
        gam3: 0.1500000,
        atomic_rad: 0.75 * AATOAU,
        pauling_en: 2.55,
        kind: 1,
        dkernel: -0.00411674,
        qkernel: 0.00213583,
        mp_rad: 3.0,
        mp_vcn: 3.0,
    },
    // Z=7, N
    Gfn2ElementParams {
        z: 7,
        nshell: 2,
        ang_shell: [0, 1, 0],
        selfenergy: [-16.686243, -12.523956, 0.0],
        slater: [2.339881, 2.014332, 0.0],
        ref_occ: [1.5, 3.5, 0.0],
        hubbard: 0.461493,
        shell_hubbard: [1.0, 1.1164892, 1.0],
        rep_alpha: 1.682689,
        rep_zeff: 5.242592,
        n_valence: 5,
        pqn: [2, 2, 0],
        ngauss: [4, 4, 0],
        shpoly: [-0.08506003, -0.02504201, 0.0],
        kcn: [-0.1955336, 0.0561076, 0.0],
        gam3: -0.0639780,
        atomic_rad: 0.71 * AATOAU,
        pauling_en: 3.04,
        kind: 1,
        dkernel: 0.03521273,
        qkernel: 0.02026786,
        mp_rad: 1.9,
        mp_vcn: 3.0,
    },
    // Z=8, O
    Gfn2ElementParams {
        z: 8,
        nshell: 2,
        ang_shell: [0, 1, 0],
        selfenergy: [-20.229985, -15.503117, 0.0],
        slater: [2.439742, 2.137023, 0.0],
        ref_occ: [2.0, 4.0, 0.0],
        hubbard: 0.451896,
        shell_hubbard: [1.0, 1.1497020, 1.0],
        rep_alpha: 2.165712,
        rep_zeff: 5.784415,
        n_valence: 6,
        pqn: [2, 2, 0],
        ngauss: [4, 4, 0],
        shpoly: [-0.14955291, -0.03350819, 0.0],
        kcn: [0.0117826, -0.0145102, 0.0],
        gam3: -0.0517134,
        atomic_rad: 0.64 * AATOAU,
        pauling_en: 3.44,
        kind: 1,
        dkernel: -0.04935670,
        qkernel: -0.00310828,
        mp_rad: 1.8,
        mp_vcn: 2.0,
    },
    // Z=9, F
    Gfn2ElementParams {
        z: 9,
        nshell: 2,
        ang_shell: [0, 1, 0],
        selfenergy: [-23.458179, -15.746583, 0.0],
        slater: [2.416361, 2.308399, 0.0],
        ref_occ: [2.0, 5.0, 0.0],
        hubbard: 0.531518,
        shell_hubbard: [1.0, 1.1677376, 1.0],
        rep_alpha: 2.421394,
        rep_zeff: 7.021486,
        n_valence: 7,
        pqn: [2, 2, 0],
        ngauss: [4, 4, 0],
        shpoly: [-0.13011924, -0.12300828, 0.0],
        kcn: [0.0394362, -0.0538373, 0.0],
        gam3: 0.1426212,
        atomic_rad: 0.60 * AATOAU,
        pauling_en: 3.98,
        kind: 1,
        dkernel: -0.08339183,
        qkernel: -0.00245955,
        mp_rad: 2.4,
        mp_vcn: 1.0,
    },
    // Z=14, Si
    Gfn2ElementParams {
        z: 14,
        nshell: 3,
        ang_shell: [0, 1, 2],
        selfenergy: [-14.360932, -6.915131, -1.825036],
        slater: [1.773917, 1.718996, 1.250000],
        ref_occ: [1.5, 2.5, 0.0],
        hubbard: 0.720000,
        shell_hubbard: [1.0000000, 0.4419958, 0.7700000],
        rep_alpha: 1.187323,
        rep_zeff: 40.001111,
        n_valence: 4,
        pqn: [3, 3, 3],
        ngauss: [4, 4, 3],
        shpoly: [0.02358522, -0.07900406, 0.11366185],
        kcn: [0.1858479, -0.1383073, -0.1935494],
        gam3: 0.1936289,
        atomic_rad: 1.14 * AATOAU,
        pauling_en: 1.90,
        kind: 1,
        dkernel: -0.00025750,
        qkernel: -0.00080000,
        mp_rad: 3.9,
        mp_vcn: 3.0,
    },
    // Z=15, P
    Gfn2ElementParams {
        z: 15,
        nshell: 3,
        ang_shell: [0, 1, 2],
        selfenergy: [-17.518756, -9.842286, -0.444893],
        slater: [1.816945, 1.903247, 1.167533],
        ref_occ: [1.5, 3.5, 0.0],
        hubbard: 0.297739,
        shell_hubbard: [1.0000000, 0.8441940, 0.6500000],
        rep_alpha: 1.143343,
        rep_zeff: 19.683502,
        n_valence: 5,
        pqn: [3, 3, 3],
        ngauss: [4, 4, 3],
        shpoly: [-0.19831771, -0.05515577, 0.26397535],
        kcn: [0.0547610, -0.0489930, 0.2429507],
        gam3: 0.0711291,
        atomic_rad: 1.09 * AATOAU,
        pauling_en: 2.19,
        kind: 1,
        dkernel: 0.02110225,
        qkernel: 0.00028679,
        mp_rad: 2.1,
        mp_vcn: 3.0,
    },
    // Z=16, S
    Gfn2ElementParams {
        z: 16,
        nshell: 3,
        ang_shell: [0, 1, 2],
        selfenergy: [-20.029654, -11.377694, -0.420282],
        slater: [1.981333, 2.025643, 1.702555],
        ref_occ: [2.0, 4.0, 0.0],
        hubbard: 0.339971,
        shell_hubbard: [1.0000000, 0.8914134, 0.7500000],
        rep_alpha: 1.214553,
        rep_zeff: 14.995090,
        n_valence: 6,
        pqn: [3, 3, 3],
        ngauss: [4, 4, 3],
        shpoly: [-0.25855520, -0.08048064, 0.25993857],
        kcn: [-0.0256951, -0.0098465, 0.2007690],
        gam3: -0.0501722,
        atomic_rad: 1.04 * AATOAU,
        pauling_en: 2.58,
        kind: 1,
        dkernel: -0.00151117,
        qkernel: 0.00442859,
        mp_rad: 3.1,
        mp_vcn: 3.0,
    },
    // Z=17, Cl
    Gfn2ElementParams {
        z: 17,
        nshell: 3,
        ang_shell: [0, 1, 2],
        selfenergy: [-29.278781, -12.673758, -0.240338],
        slater: [2.485265, 2.199650, 2.476089],
        ref_occ: [2.0, 5.0, 0.0],
        hubbard: 0.248514,
        shell_hubbard: [1.0000000, 1.4989400, 1.5000000],
        rep_alpha: 1.577144,
        rep_zeff: 17.353134,
        n_valence: 7,
        pqn: [3, 3, 3],
        ngauss: [4, 4, 3],
        shpoly: [-0.16562004, -0.06986430, 0.38045622],
        kcn: [0.0617972, -0.0181618, 0.1672768],
        gam3: 0.1495483,
        atomic_rad: 1.00 * AATOAU,
        pauling_en: 3.16,
        kind: 1,
        dkernel: -0.02536958,
        qkernel: 0.00122783,
        mp_rad: 2.5,
        mp_vcn: 1.0,
    },
    // Z=35, Br
    Gfn2ElementParams {
        z: 35,
        nshell: 3,
        ang_shell: [0, 1, 2],
        selfenergy: [-23.583718, -12.588824, 0.047980],
        slater: [2.077587, 2.263120, 1.845038],
        ref_occ: [2.0, 5.0, 0.0],
        hubbard: 0.261253,
        shell_hubbard: [1.0, 1.5203002, 1.4],
        rep_alpha: 1.296174,
        rep_zeff: 32.845361,
        n_valence: 7,
        pqn: [4, 4, 4],
        ngauss: [4, 4, 3],
        shpoly: [-0.25005079, -0.14520078, 0.36614038],
        kcn: [0.0006150, -0.0058347, 0.2250180],
        gam3: 0.1300000,
        atomic_rad: 1.17 * AATOAU,
        pauling_en: 2.96,
        kind: 1,
        dkernel: -0.01088586,
        qkernel: 0.00216935,
        mp_rad: 3.9,
        mp_vcn: 1.0,
    },
    // Z=53, I
    Gfn2ElementParams {
        z: 53,
        nshell: 3,
        ang_shell: [0, 1, 2],
        selfenergy: [-20.949407, -12.180159, -0.266596],
        slater: [2.159500, 2.308379, 1.691185],
        ref_occ: [2.0, 5.0, 0.0],
        hubbard: 0.383124,
        shell_hubbard: [1.0000000, 0.9661265, 1.3000000],
        rep_alpha: 1.017946,
        rep_zeff: 63.319176,
        n_valence: 7,
        pqn: [5, 5, 5],
        ngauss: [4, 4, 3],
        shpoly: [-0.26957547, -0.14183312, 0.28211905],
        kcn: [-0.0506150, 0.0084766, 0.3077127],
        gam3: 0.1200000,
        atomic_rad: 1.33 * AATOAU,
        pauling_en: 2.66,
        kind: 1,
        dkernel: -0.00603648,
        qkernel: 0.00069660,
        mp_rad: 5.0,
        mp_vcn: 1.0,
    },
    // ---- Transition metals (3d) ----
    // Z=22, Ti  (shell order: d=0, s=1, p=2)
    Gfn2ElementParams {
        z: 22,
        nshell: 3,
        ang_shell: [2, 0, 1],
        selfenergy: [-7.234331, -10.900000, -1.928783],
        slater: [1.849994, 1.469983, 0.957410],
        ref_occ: [1.0, 1.0, 1.0],
        hubbard: 0.513586,
        shell_hubbard: [0.5078886, 1.0000000, 0.6200000],
        rep_alpha: 0.723104,
        rep_zeff: 12.036336,
        n_valence: 4,
        pqn: [3, 4, 4],
        ngauss: [3, 4, 4],
        shpoly: [-0.27724389, 0.04561234, 0.51801626],
        kcn: [0.1028188, 0.1007021, -0.0237074],
        gam3: 0.1767268,
        atomic_rad: 1.36 * AATOAU,
        pauling_en: 1.54,
        kind: 2,
        dkernel: -0.00434506,
        qkernel: 0.00059660,
        mp_rad: 5.0,
        mp_vcn: 4.0,
    },
    // Z=24, Cr
    Gfn2ElementParams {
        z: 24,
        nshell: 3,
        ang_shell: [2, 0, 1],
        selfenergy: [-7.209794, -9.201304, -0.696957],
        slater: [1.568211, 1.395427, 1.080270],
        ref_occ: [1.0, 1.0, 4.0],
        hubbard: 0.396299,
        shell_hubbard: [1.7405872, 1.0000000, 0.5300000],
        rep_alpha: 0.966993,
        rep_zeff: 19.517914,
        n_valence: 6,
        pqn: [3, 4, 4],
        ngauss: [3, 4, 4],
        shpoly: [-0.27971622, 0.13376234, 0.48092152],
        kcn: [0.0289291, -0.0232087, -0.0188919],
        gam3: 0.0300000,
        atomic_rad: 1.22 * AATOAU,
        pauling_en: 1.66,
        kind: 2,
        dkernel: 0.00149669,
        qkernel: 0.00137744,
        mp_rad: 5.0,
        mp_vcn: 6.0,
    },
    // Z=25, Mn
    Gfn2ElementParams {
        z: 25,
        nshell: 3,
        ang_shell: [2, 0, 1],
        selfenergy: [-10.120933, -5.617346, -4.198724],
        slater: [1.839250, 1.222190, 1.240215],
        ref_occ: [1.0, 1.0, 5.0],
        hubbard: 0.346651,
        shell_hubbard: [1.0545811, 1.0000000, 0.4000000],
        rep_alpha: 1.071100,
        rep_zeff: 18.760605,
        n_valence: 7,
        pqn: [3, 4, 4],
        ngauss: [3, 4, 4],
        shpoly: [-0.31255885, 0.28519691, 0.26346555],
        kcn: [-0.0195827, -0.0275000, -0.0015839],
        gam3: 0.0600000,
        atomic_rad: 1.19 * AATOAU,
        pauling_en: 1.55,
        kind: 2,
        dkernel: -0.00759168,
        qkernel: 0.00229903,
        mp_rad: 5.0,
        mp_vcn: 6.0,
    },
    // Z=26, Fe
    Gfn2ElementParams {
        z: 26,
        nshell: 3,
        ang_shell: [2, 0, 1],
        selfenergy: [-10.035473, -5.402911, -3.308988],
        slater: [1.911049, 1.022393, 1.294467],
        ref_occ: [1.0, 1.0, 6.0],
        hubbard: 0.271594,
        shell_hubbard: [1.4046615, 1.0000000, 0.3500000],
        rep_alpha: 1.113422,
        rep_zeff: 20.360089,
        n_valence: 8,
        pqn: [3, 4, 4],
        ngauss: [3, 4, 4],
        shpoly: [-0.28614961, 0.11527794, 0.39459890],
        kcn: [-0.0274654, -0.4049876, -0.0756480],
        gam3: -0.0500000,
        atomic_rad: 1.16 * AATOAU,
        pauling_en: 1.83,
        kind: 2,
        dkernel: 0.00412929,
        qkernel: 0.00267734,
        mp_rad: 5.0,
        mp_vcn: 6.0,
    },
    // Z=27, Co
    Gfn2ElementParams {
        z: 27,
        nshell: 3,
        ang_shell: [2, 0, 1],
        selfenergy: [-10.580430, -8.596723, -2.585753],
        slater: [2.326507, 1.464221, 1.298678],
        ref_occ: [1.0, 1.0, 7.0],
        hubbard: 0.477760,
        shell_hubbard: [0.7581507, 1.0000000, 0.3500000],
        rep_alpha: 1.241717,
        rep_zeff: 27.127744,
        n_valence: 9,
        pqn: [3, 4, 4],
        ngauss: [3, 4, 4],
        shpoly: [-0.22355636, 0.09168460, 0.25424719],
        kcn: [0.0121980, -0.0227872, 0.0076513],
        gam3: 0.0300000,
        atomic_rad: 1.11 * AATOAU,
        pauling_en: 1.88,
        kind: 2,
        dkernel: -0.00247938,
        qkernel: 0.00048237,
        mp_rad: 5.0,
        mp_vcn: 6.0,
    },
    // Z=28, Ni
    Gfn2ElementParams {
        z: 28,
        nshell: 3,
        ang_shell: [2, 0, 1],
        selfenergy: [-12.712236, -8.524281, -2.878873],
        slater: [2.430756, 1.469945, 1.317046],
        ref_occ: [1.0, 1.0, 8.0],
        hubbard: 0.344970,
        shell_hubbard: [0.9388812, 1.0000000, 0.4000000],
        rep_alpha: 1.077516,
        rep_zeff: 10.533269,
        n_valence: 10,
        pqn: [3, 4, 4],
        ngauss: [3, 4, 4],
        shpoly: [-0.25385640, 0.20839550, 0.30886445],
        kcn: [-0.0066417, 0.0310301, 0.0226796],
        gam3: -0.0200000,
        atomic_rad: 1.10 * AATOAU,
        pauling_en: 1.91,
        kind: 2,
        dkernel: -0.01261887,
        qkernel: -0.00080000,
        mp_rad: 5.0,
        mp_vcn: 4.0,
    },
    // Z=29, Cu
    Gfn2ElementParams {
        z: 29,
        nshell: 3,
        ang_shell: [2, 0, 1],
        selfenergy: [-9.506548, -6.922958, -2.267723],
        slater: [2.375425, 1.550837, 1.984703],
        ref_occ: [1.0, 0.0, 10.0],
        hubbard: 0.202969,
        shell_hubbard: [2.3333066, 1.0000000, 1.0700000],
        rep_alpha: 0.998768,
        rep_zeff: 9.913846,
        n_valence: 11,
        pqn: [3, 4, 4],
        ngauss: [3, 4, 4],
        shpoly: [-0.26508943, 0.17798264, 0.14977818],
        kcn: [-0.0173684, 0.3349047, -0.2619446],
        gam3: 0.0500000,
        atomic_rad: 1.12 * AATOAU,
        pauling_en: 1.90,
        kind: 2,
        dkernel: -0.00700000,
        qkernel: -0.00345631,
        mp_rad: 5.0,
        mp_vcn: 4.0,
    },
    // Z=30, Zn
    Gfn2ElementParams {
        z: 30,
        nshell: 3,
        ang_shell: [2, 0, 1],
        selfenergy: [-7.177294, -0.991895, 0.000000],
        slater: [1.664847, 1.176434, 0.000000],
        ref_occ: [2.0, 0.0, 10.0],
        hubbard: 0.564152,
        shell_hubbard: [1.0000000, 1.0000000, 1.0684343],
        rep_alpha: 1.160262,
        rep_zeff: 22.099503,
        n_valence: 12,
        pqn: [3, 4, 4],
        ngauss: [3, 4, 4],
        shpoly: [0.00000000, -0.09240315, 0.22271839],
        kcn: [0.0000000, 0.2011910, -0.0055135],
        gam3: 0.2312896,
        atomic_rad: 1.18 * AATOAU,
        pauling_en: 1.65,
        kind: 2,
        dkernel: -0.00100000,
        qkernel: 0.00007658,
        mp_rad: 5.0,
        mp_vcn: 2.0,
    },
    // ---- Second-row TM ----
    // Z=44, Ru
    Gfn2ElementParams {
        z: 44,
        nshell: 3,
        ang_shell: [2, 0, 1],
        selfenergy: [-10.285405, -5.332608, -3.307153],
        slater: [2.102697, 1.749643, 1.348322],
        ref_occ: [1.0, 1.0, 6.0],
        hubbard: 0.711205,
        shell_hubbard: [0.6947695, 1.0000000, 0.3500000],
        rep_alpha: 1.031669,
        rep_zeff: 28.070247,
        n_valence: 8,
        pqn: [4, 5, 5],
        ngauss: [4, 4, 3],
        shpoly: [-0.27583213, 0.10106270, 0.34028722],
        kcn: [-0.0263914, 0.4471162, -0.0034723],
        gam3: -0.0500000,
        atomic_rad: 1.25 * AATOAU,
        pauling_en: 2.20,
        kind: 2,
        dkernel: -0.00081922,
        qkernel: 0.00359228,
        mp_rad: 5.0,
        mp_vcn: 6.0,
    },
    // Z=46, Pd
    Gfn2ElementParams {
        z: 46,
        nshell: 3,
        ang_shell: [2, 0, 1],
        selfenergy: [-11.963518, -9.714059, -2.035281],
        slater: [2.353691, 1.828354, 1.333352],
        ref_occ: [1.0, 0.0, 10.0],
        hubbard: 0.273310,
        shell_hubbard: [1.0931707, 1.0000000, 0.4000000],
        rep_alpha: 1.092745,
        rep_zeff: 28.674700,
        n_valence: 10,
        pqn: [4, 5, 5],
        ngauss: [4, 4, 3],
        shpoly: [-0.27173113, 0.06200145, 0.45341322],
        kcn: [0.0060285, 0.0266820, 0.0503075],
        gam3: 0.0800000,
        atomic_rad: 1.20 * AATOAU,
        pauling_en: 2.20,
        kind: 2,
        dkernel: -0.00310361,
        qkernel: -0.00040485,
        mp_rad: 5.0,
        mp_vcn: 4.0,
    },
    // Z=47, Ag
    Gfn2ElementParams {
        z: 47,
        nshell: 3,
        ang_shell: [2, 0, 1],
        selfenergy: [-9.591083, -8.083960, -2.934333],
        slater: [2.843549, 1.798462, 1.266649],
        ref_occ: [2.0, 0.0, 10.0],
        hubbard: 0.263740,
        shell_hubbard: [1.8024848, 1.0000000, 0.9700000],
        rep_alpha: 0.678344,
        rep_zeff: 6.493286,
        n_valence: 11,
        pqn: [4, 5, 5],
        ngauss: [4, 4, 3],
        shpoly: [-0.16490742, 0.01091490, 0.11561444],
        kcn: [-0.0062719, -0.0065794, 0.1677171],
        gam3: 0.0200000,
        atomic_rad: 1.28 * AATOAU,
        pauling_en: 1.93,
        kind: 2,
        dkernel: -0.00800314,
        qkernel: -0.00020810,
        mp_rad: 5.0,
        mp_vcn: 4.0,
    },
    // ---- Third-row TM ----
    // Z=78, Pt
    Gfn2ElementParams {
        z: 78,
        nshell: 3,
        ang_shell: [2, 0, 1],
        selfenergy: [-12.047728, -10.482306, -3.778297],
        slater: [2.704492, 2.329136, 1.623286],
        ref_occ: [1.0, 1.0, 8.0],
        hubbard: 0.353574,
        shell_hubbard: [0.9796788, 1.0000000, 0.4000000],
        rep_alpha: 1.204081,
        rep_zeff: 53.881425,
        n_valence: 10,
        pqn: [5, 6, 6],
        ngauss: [3, 6, 6],
        shpoly: [-0.27243898, 0.06737556, 0.19259455],
        kcn: [-0.0204828, 0.1139530, 0.1408029],
        gam3: 0.0600000,
        atomic_rad: 1.23 * AATOAU,
        pauling_en: 2.28,
        kind: 2,
        dkernel: -0.00423684,
        qkernel: 0.00098019,
        mp_rad: 5.0,
        mp_vcn: 4.0,
    },
    // Z=79, Au
    Gfn2ElementParams {
        z: 79,
        nshell: 3,
        ang_shell: [2, 0, 1],
        selfenergy: [-9.578599, -7.688552, 0.883399],
        slater: [3.241287, 2.183171, 2.084484],
        ref_occ: [2.0, 0.0, 10.0],
        hubbard: 0.438997,
        shell_hubbard: [1.0614126, 1.0000000, 0.4000000],
        rep_alpha: 0.919210,
        rep_zeff: 14.711475,
        n_valence: 11,
        pqn: [5, 6, 6],
        ngauss: [3, 6, 6],
        shpoly: [-0.06410815, 0.04691539, 0.25250274],
        kcn: [-0.0154462, 0.1479337, 0.1048065],
        gam3: 0.0850000,
        atomic_rad: 1.24 * AATOAU,
        pauling_en: 2.54,
        kind: 2,
        dkernel: 0.00393418,
        qkernel: -0.00020320,
        mp_rad: 5.0,
        mp_vcn: 4.0,
    },
];

/// Look up GFN2-xTB parameters by atomic number.
pub fn get_gfn2_params(z: u8) -> Option<&'static Gfn2ElementParams> {
    GFN2_PARAMS.iter().find(|p| p.z == z)
}

/// Check if an element is supported by GFN2-xTB.
pub fn is_gfn2_supported(z: u8) -> bool {
    get_gfn2_params(z).is_some()
}

/// Compute effective Hubbard parameter (Hartree) as occupation-weighted average
/// across shells: η_eff = Σ_sh(n_sh × η × shell_hubbard_sh) / Σ_sh(n_sh).
/// This approximates shell-resolved SCC with an atom-level gamma.
pub fn effective_hubbard(p: &Gfn2ElementParams) -> f64 {
    let mut num = 0.0;
    let mut den = 0.0;
    for sh in 0..p.nshell {
        let n = p.ref_occ[sh];
        num += n * p.hubbard * p.shell_hubbard[sh];
        den += n;
    }
    if den > 0.0 {
        num / den
    } else {
        p.hubbard
    }
}

/// Count basis functions for GFN2-xTB: s=1, p=3, d=5.
pub fn gfn2_basis_size(z: u8) -> usize {
    match get_gfn2_params(z) {
        None => 0,
        Some(p) => {
            let mut n = 0;
            for i in 0..p.nshell {
                n += match p.ang_shell[i] {
                    0 => 1,
                    1 => 3,
                    2 => 5,
                    _ => 0,
                };
            }
            n
        }
    }
}

/// Count total valence electrons.
pub fn gfn2_electron_count(elements: &[u8]) -> usize {
    elements
        .iter()
        .map(|&z| get_gfn2_params(z).map_or(0, |p| p.n_valence as usize))
        .sum()
}

/// Compute kshell factor for off-diagonal H0 scaling.
/// GFN2 uses: kshell(k, l) = merge(2.0, (kdiag[l]+kdiag[k])/2, ...)
pub fn kshell_factor(l_a: u8, l_b: u8) -> f64 {
    let la = l_a as usize;
    let lb = l_b as usize;
    // sp/ps and sd/ds pairs use k=2.0 (from tblite kshell function)
    if (la == 2 && (lb == 0 || lb == 1)) || (lb == 2 && (la == 0 || la == 1)) {
        2.0
    } else {
        (KDIAG[la] + KDIAG[lb]) / 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gfn2_params_hydrogen() {
        let p = get_gfn2_params(1).unwrap();
        assert_eq!(p.nshell, 1);
        assert!((p.selfenergy[0] - (-10.707211)).abs() < 1e-4);
        assert!((p.slater[0] - 1.23).abs() < 1e-4);
        assert_eq!(p.n_valence, 1);
    }

    #[test]
    fn test_gfn2_params_carbon() {
        let p = get_gfn2_params(6).unwrap();
        assert_eq!(p.nshell, 2);
        assert!((p.selfenergy[0] - (-13.970922)).abs() < 1e-4);
        assert!((p.selfenergy[1] - (-10.063292)).abs() < 1e-4);
        assert_eq!(p.n_valence, 4);
    }

    #[test]
    fn test_gfn2_basis_size() {
        assert_eq!(gfn2_basis_size(1), 1); // H: s
        assert_eq!(gfn2_basis_size(6), 4); // C: s+p
        assert_eq!(gfn2_basis_size(26), 9); // Fe: d+s+p
    }

    #[test]
    fn test_gfn2_electron_count() {
        assert_eq!(gfn2_electron_count(&[8, 1, 1]), 8); // H2O
        assert_eq!(gfn2_electron_count(&[6, 1, 1, 1, 1]), 8); // CH4
    }
}
