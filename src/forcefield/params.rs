/// UFF atom type parameters.
/// Source: Rappé, A. K. et al. "UFF, a full periodic table force field for molecular mechanics
/// and molecular dynamics simulations." J. Am. Chem. Soc. 114 (1992): 10024–10035. Table 1.
#[derive(Debug, Clone)]
pub struct UffAtomParams {
    pub label: &'static str,
    pub element: u8, // Atomic number
    pub r1: f64,     // Bond radius (Angstroms) — used in bond length formula
    pub theta0: f64, // Equilibrium valence angle (degrees, convert to radians in builder)
    pub x1: f64,     // VDW distance parameter (Angstroms)
    pub d1: f64,     // VDW well depth (kcal/mol)
    pub zeta: f64,   // GMP electronegativity scale factor
    pub z_star: f64, // Effective nuclear charge for bond force constant
    pub chi: f64,    // GMP electronegativity (for bond length correction)
    pub v_tors: f64, // Torsion barrier Vi (kcal/mol), used as sqrt(Vi*Vj)/2 (eq. 16)
}

/// Retrieve UFF parameters by atom type label.
pub fn get_uff_params(label: &str) -> Option<UffAtomParams> {
    match label {
        // Hydrogen
        "H_" => Some(UffAtomParams {
            label: "H_",
            element: 1,
            r1: 0.354,
            theta0: 180.0,
            x1: 2.886,
            d1: 0.044,
            zeta: 0.712,
            z_star: 0.712,
            chi: 4.528,
            v_tors: 0.0,
        }),
        // Carbon
        "C_3" => Some(UffAtomParams {
            label: "C_3",
            element: 6,
            r1: 0.757,
            theta0: 109.47,
            x1: 3.851,
            d1: 0.105,
            zeta: 1.912,
            z_star: 1.912,
            chi: 5.343,
            v_tors: 2.119,
        }),
        "C_2" => Some(UffAtomParams {
            label: "C_2",
            element: 6,
            r1: 0.732,
            theta0: 120.0,
            x1: 3.851,
            d1: 0.105,
            zeta: 1.912,
            z_star: 1.912,
            chi: 5.343,
            v_tors: 2.119,
        }),
        "C_R" => Some(UffAtomParams {
            label: "C_R",
            element: 6,
            r1: 0.729,
            theta0: 120.0,
            x1: 3.851,
            d1: 0.105,
            zeta: 1.912,
            z_star: 1.912,
            chi: 5.343,
            v_tors: 2.119,
        }),
        "C_1" => Some(UffAtomParams {
            label: "C_1",
            element: 6,
            r1: 0.706,
            theta0: 180.0,
            x1: 3.851,
            d1: 0.105,
            zeta: 1.912,
            z_star: 1.912,
            chi: 5.343,
            v_tors: 2.119,
        }),
        // Nitrogen
        "N_3" => Some(UffAtomParams {
            label: "N_3",
            element: 7,
            r1: 0.700,
            theta0: 106.7,
            x1: 3.660,
            d1: 0.069,
            zeta: 2.544,
            z_star: 2.544,
            chi: 6.899,
            v_tors: 0.450,
        }),
        "N_2" => Some(UffAtomParams {
            label: "N_2",
            element: 7,
            r1: 0.685,
            theta0: 111.2,
            x1: 3.660,
            d1: 0.069,
            zeta: 2.544,
            z_star: 2.544,
            chi: 6.899,
            v_tors: 0.450,
        }),
        "N_R" => Some(UffAtomParams {
            label: "N_R",
            element: 7,
            r1: 0.699,
            theta0: 120.0,
            x1: 3.660,
            d1: 0.069,
            zeta: 2.544,
            z_star: 2.544,
            chi: 6.899,
            v_tors: 0.450,
        }),
        "N_1" => Some(UffAtomParams {
            label: "N_1",
            element: 7,
            r1: 0.656,
            theta0: 180.0,
            x1: 3.660,
            d1: 0.069,
            zeta: 2.544,
            z_star: 2.544,
            chi: 6.899,
            v_tors: 0.450,
        }),
        // Oxygen
        "O_3" => Some(UffAtomParams {
            label: "O_3",
            element: 8,
            r1: 0.658,
            theta0: 104.51,
            x1: 3.500,
            d1: 0.060,
            zeta: 2.300,
            z_star: 2.300,
            chi: 8.741,
            v_tors: 2.100,
        }),
        "O_2" => Some(UffAtomParams {
            label: "O_2",
            element: 8,
            r1: 0.634,
            theta0: 120.0,
            x1: 3.500,
            d1: 0.060,
            zeta: 2.300,
            z_star: 2.300,
            chi: 8.741,
            v_tors: 2.100,
        }),
        "O_R" => Some(UffAtomParams {
            label: "O_R",
            element: 8,
            r1: 0.680,
            theta0: 110.0,
            x1: 3.500,
            d1: 0.060,
            zeta: 2.300,
            z_star: 2.300,
            chi: 8.741,
            v_tors: 2.100,
        }),
        "O_1" => Some(UffAtomParams {
            label: "O_1",
            element: 8,
            r1: 0.639,
            theta0: 180.0,
            x1: 3.500,
            d1: 0.060,
            zeta: 2.300,
            z_star: 2.300,
            chi: 8.741,
            v_tors: 2.100,
        }),
        // Halogens
        "F_" => Some(UffAtomParams {
            label: "F_",
            element: 9,
            r1: 0.668,
            theta0: 180.0,
            x1: 3.364,
            d1: 0.050,
            zeta: 1.735,
            z_star: 1.735,
            chi: 10.874,
            v_tors: 1.082,
        }),
        "Cl" => Some(UffAtomParams {
            label: "Cl",
            element: 17,
            r1: 1.029,
            theta0: 109.47,
            x1: 3.947,
            d1: 0.227,
            zeta: 2.315,
            z_star: 2.315,
            chi: 8.564,
            v_tors: 1.673,
        }),
        "Br" => Some(UffAtomParams {
            label: "Br",
            element: 35,
            r1: 1.137,
            theta0: 109.47,
            x1: 4.189,
            d1: 0.251,
            zeta: 2.283,
            z_star: 2.283,
            chi: 7.790,
            v_tors: 1.682,
        }),
        "I_" => Some(UffAtomParams {
            label: "I_",
            element: 53,
            r1: 1.353,
            theta0: 109.47,
            x1: 4.500,
            d1: 0.339,
            zeta: 2.650,
            z_star: 2.650,
            chi: 6.822,
            v_tors: 1.530,
        }),
        // Sulfur
        "S_3" => Some(UffAtomParams {
            label: "S_3",
            element: 16,
            r1: 1.020,
            theta0: 92.1,
            x1: 4.035,
            d1: 0.274,
            zeta: 2.703,
            z_star: 2.703,
            chi: 10.000,
            v_tors: 1.948,
        }),
        "S_2" => Some(UffAtomParams {
            label: "S_2",
            element: 16,
            r1: 0.987,
            theta0: 120.0,
            x1: 4.035,
            d1: 0.274,
            zeta: 2.703,
            z_star: 2.703,
            chi: 10.000,
            v_tors: 1.948,
        }),
        "S_R" => Some(UffAtomParams {
            label: "S_R",
            element: 16,
            r1: 0.990,
            theta0: 120.0,
            x1: 4.035,
            d1: 0.274,
            zeta: 2.703,
            z_star: 2.703,
            chi: 10.000,
            v_tors: 1.948,
        }),
        // Phosphorus
        "P_3" => Some(UffAtomParams {
            label: "P_3",
            element: 15,
            r1: 0.930,
            theta0: 93.8,
            x1: 4.147,
            d1: 0.305,
            zeta: 2.862,
            z_star: 2.862,
            chi: 8.000,
            v_tors: 2.421,
        }),
        // Boron
        "B_2" => Some(UffAtomParams {
            label: "B_2",
            element: 5,
            r1: 0.828,
            theta0: 120.0,
            x1: 3.637,
            d1: 0.180,
            zeta: 1.755,
            z_star: 1.755,
            chi: 4.961,
            v_tors: 0.0,
        }),
        // Silicon
        "Si3" => Some(UffAtomParams {
            label: "Si3",
            element: 14,
            r1: 1.117,
            theta0: 109.47,
            x1: 4.295,
            d1: 0.402,
            zeta: 2.323,
            z_star: 2.323,
            chi: 4.168,
            v_tors: 1.225,
        }),
        // ─── Transition metals (Rappé 1992, Table 1) ─────────────────────
        "Ti3+4" => Some(UffAtomParams {
            label: "Ti3+4",
            element: 22,
            r1: 1.412,
            theta0: 109.47,
            x1: 3.175,
            d1: 0.017,
            zeta: 2.659,
            z_star: 2.659,
            chi: 3.470,
            v_tors: 0.0,
        }),
        "V_3+5" => Some(UffAtomParams {
            label: "V_3+5",
            element: 23,
            r1: 1.402,
            theta0: 109.47,
            x1: 3.144,
            d1: 0.016,
            zeta: 2.679,
            z_star: 2.679,
            chi: 3.650,
            v_tors: 0.0,
        }),
        "Cr6+3" => Some(UffAtomParams {
            label: "Cr6+3",
            element: 24,
            r1: 1.345,
            theta0: 90.0,
            x1: 3.023,
            d1: 0.015,
            zeta: 2.463,
            z_star: 2.463,
            chi: 3.415,
            v_tors: 0.0,
        }),
        "Mn6+2" => Some(UffAtomParams {
            label: "Mn6+2",
            element: 25,
            r1: 1.382,
            theta0: 90.0,
            x1: 2.961,
            d1: 0.013,
            zeta: 2.430,
            z_star: 2.430,
            chi: 3.325,
            v_tors: 0.0,
        }),
        "Fe3+2" | "Fe6+2" => Some(UffAtomParams {
            label: "Fe3+2",
            element: 26,
            r1: 1.335,
            theta0: 109.47,
            x1: 2.912,
            d1: 0.013,
            zeta: 2.430,
            z_star: 2.430,
            chi: 3.760,
            v_tors: 0.0,
        }),
        "Co6+3" => Some(UffAtomParams {
            label: "Co6+3",
            element: 27,
            r1: 1.241,
            theta0: 90.0,
            x1: 2.872,
            d1: 0.014,
            zeta: 2.430,
            z_star: 2.430,
            chi: 4.105,
            v_tors: 0.0,
        }),
        "Ni4+2" => Some(UffAtomParams {
            label: "Ni4+2",
            element: 28,
            r1: 1.164,
            theta0: 90.0,
            x1: 2.834,
            d1: 0.015,
            zeta: 2.430,
            z_star: 2.430,
            chi: 4.465,
            v_tors: 0.0,
        }),
        "Cu3+1" => Some(UffAtomParams {
            label: "Cu3+1",
            element: 29,
            r1: 1.302,
            theta0: 109.47,
            x1: 3.495,
            d1: 0.005,
            zeta: 1.756,
            z_star: 1.756,
            chi: 4.200,
            v_tors: 0.0,
        }),
        "Zn3+2" => Some(UffAtomParams {
            label: "Zn3+2",
            element: 30,
            r1: 1.193,
            theta0: 109.47,
            x1: 2.763,
            d1: 0.124,
            zeta: 1.308,
            z_star: 1.308,
            chi: 5.106,
            v_tors: 0.0,
        }),
        // ─── Additional main-group elements ──────────────────────────────
        "Ge3" => Some(UffAtomParams {
            label: "Ge3",
            element: 32,
            r1: 1.197,
            theta0: 109.47,
            x1: 4.280,
            d1: 0.379,
            zeta: 2.789,
            z_star: 2.789,
            chi: 4.051,
            v_tors: 0.701,
        }),
        "As3" => Some(UffAtomParams {
            label: "As3",
            element: 33,
            r1: 1.211,
            theta0: 92.1,
            x1: 4.230,
            d1: 0.309,
            zeta: 2.864,
            z_star: 2.864,
            chi: 5.188,
            v_tors: 1.476,
        }),
        "Se3" => Some(UffAtomParams {
            label: "Se3",
            element: 34,
            r1: 1.190,
            theta0: 90.6,
            x1: 4.205,
            d1: 0.291,
            zeta: 2.764,
            z_star: 2.764,
            chi: 6.428,
            v_tors: 0.335,
        }),
        "Sn3" => Some(UffAtomParams {
            label: "Sn3",
            element: 50,
            r1: 1.399,
            theta0: 109.47,
            x1: 4.392,
            d1: 0.567,
            zeta: 3.198,
            z_star: 3.198,
            chi: 4.000,
            v_tors: 0.199,
        }),
        "Sb3" => Some(UffAtomParams {
            label: "Sb3",
            element: 51,
            r1: 1.404,
            theta0: 91.6,
            x1: 4.420,
            d1: 0.449,
            zeta: 3.510,
            z_star: 3.510,
            chi: 4.899,
            v_tors: 1.100,
        }),
        "Te3" => Some(UffAtomParams {
            label: "Te3",
            element: 52,
            r1: 1.386,
            theta0: 90.25,
            x1: 4.470,
            d1: 0.398,
            zeta: 3.526,
            z_star: 3.526,
            chi: 5.816,
            v_tors: 0.300,
        }),
        "Mo3+6" => Some(UffAtomParams {
            label: "Mo3+6",
            element: 42,
            r1: 1.386,
            theta0: 109.47,
            x1: 3.052,
            d1: 0.056,
            zeta: 3.056,
            z_star: 3.056,
            chi: 3.465,
            v_tors: 0.0,
        }),
        "Pd4+2" => Some(UffAtomParams {
            label: "Pd4+2",
            element: 46,
            r1: 1.338,
            theta0: 90.0,
            x1: 2.899,
            d1: 0.048,
            zeta: 2.899,
            z_star: 2.899,
            chi: 4.320,
            v_tors: 0.0,
        }),
        "Ag1+1" => Some(UffAtomParams {
            label: "Ag1+1",
            element: 47,
            r1: 1.386,
            theta0: 180.0,
            x1: 3.148,
            d1: 0.036,
            zeta: 1.956,
            z_star: 1.956,
            chi: 4.436,
            v_tors: 0.0,
        }),
        "Pt4+2" => Some(UffAtomParams {
            label: "Pt4+2",
            element: 78,
            r1: 1.364,
            theta0: 90.0,
            x1: 2.754,
            d1: 0.080,
            zeta: 2.754,
            z_star: 2.754,
            chi: 4.790,
            v_tors: 0.0,
        }),
        "Au4+3" => Some(UffAtomParams {
            label: "Au4+3",
            element: 79,
            r1: 1.255,
            theta0: 90.0,
            x1: 3.293,
            d1: 0.039,
            zeta: 2.309,
            z_star: 2.309,
            chi: 5.770,
            v_tors: 0.0,
        }),
        _ => None,
    }
}

/// Compute the UFF equilibrium bond length between two atom types.
///
/// Rappé eq. 2–4:
/// r_ij = r_i + r_j + r_BO - r_EN
/// r_BO = -λ · (r_i + r_j) · ln(BO)          [λ = 0.1332]
/// r_EN = r_i · r_j · (√χ_i - √χ_j)² / (χ_i·r_i + χ_j·r_j)
pub fn get_uff_bond_length(
    pi: &UffAtomParams,
    pj: &UffAtomParams,
    bond_order: &crate::graph::BondOrder,
) -> f64 {
    const LAMBDA: f64 = 0.1332;

    let bo_value: f64 = match bond_order {
        crate::graph::BondOrder::Single => 1.0,
        crate::graph::BondOrder::Aromatic => 1.5,
        crate::graph::BondOrder::Double => 2.0,
        crate::graph::BondOrder::Triple => 3.0,
        _ => 1.0,
    };

    let r_bo = if (bo_value - 1.0).abs() < 1e-10 {
        0.0
    } else {
        -LAMBDA * (pi.r1 + pj.r1) * bo_value.ln()
    };

    let r_en = if (pi.chi - pj.chi).abs() > 1e-6 {
        let sqrt_diff = pi.chi.sqrt() - pj.chi.sqrt();
        pi.r1 * pj.r1 * sqrt_diff * sqrt_diff / (pi.chi * pi.r1 + pj.chi * pj.r1)
    } else {
        0.0
    };

    (pi.r1 + pj.r1 + r_bo - r_en).max(0.5)
}

/// UFF bond stretching force constant K_ij (Rappé eq. 6).
/// E_bond = ½ K_ij (r - r₀)²    →   K_ij = 664.12 · Z*_i · Z*_j / r₀³
/// 664.12 = 2 × 332.06 where 332.06 ≈ e²Nₐ/(4πε₀) in kcal·Å/e² (Coulomb constant)
pub fn get_uff_bond_force_constant(pi: &UffAtomParams, pj: &UffAtomParams, r0: f64) -> f64 {
    664.12 * pi.z_star * pj.z_star / r0.powi(3)
}

/// UFF angle bending force constant K_ijk (Rappé eq. 13, simplified).
/// K_ijk = 332.06 · Z*_i · Z*_k / (r_ij · r_jk · r_ik²)
/// 332.06 ≈ e²Nₐ/(4πε₀) in kcal·Å/e² (Coulomb constant for UFF energy units)
/// where r_ik² = r_ij² + r_jk² - 2·r_ij·r_jk·cos(θ₀)
pub fn get_uff_angle_force_constant(
    pi: &UffAtomParams,
    pk: &UffAtomParams,
    r_ij: f64,
    r_jk: f64,
    theta0_rad: f64,
) -> f64 {
    let r_ik_sq = r_ij * r_ij + r_jk * r_jk - 2.0 * r_ij * r_jk * theta0_rad.cos();
    if r_ik_sq < 1e-10 || r_ij < 1e-10 || r_jk < 1e-10 {
        return 0.0;
    }
    332.06 * pi.z_star * pk.z_star / (r_ij * r_jk * r_ik_sq)
}

/// UFF torsion barriers for bond j–k (Rappé Section 3D, eq. 16).
/// Returns (V_barrier, periodicity_n, cos_phi0).
/// V_barrier = ½ · √(V_i · V_j)  for sp3–sp3 bonds.
pub fn get_uff_torsion_params(
    hyb_j: &crate::graph::Hybridization,
    hyb_k: &crate::graph::Hybridization,
    vj: f64,
    vk: f64,
    bond_order: &crate::graph::BondOrder,
) -> (f64, f64, f64) {
    use crate::graph::BondOrder::*;
    use crate::graph::Hybridization::*;

    match (hyb_j, hyb_k, bond_order) {
        // sp3–sp3: three-fold (staggered → anti preferred)
        (SP3, SP3, _) => (0.5 * (vj * vk).sqrt(), 3.0, 1.0),
        // sp3–sp2 mixed: six-fold (nearly flat)
        (SP3, SP2, _) | (SP2, SP3, _) => (1.0, 6.0, 1.0),
        // sp2–sp2 conjugated (aromatic or double bond): strong planarity
        (SP2, SP2, Aromatic) | (SP2, SP2, Double) => (5.0, 2.0, -1.0),
        // sp2–sp2 single (biaryl/amide): moderate planarity
        (SP2, SP2, Single) => (10.0, 2.0, -1.0),
        // sp (linear): no torsion
        (SP, _, _) | (_, SP, _) => (0.0, 1.0, 1.0),
        _ => (1.0, 3.0, 1.0),
    }
}
