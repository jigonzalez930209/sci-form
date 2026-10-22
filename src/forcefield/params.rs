//! UFF parameter tables and the rules that turn them into force constants.
//!
//! The tabulated values are the canonical UFF numbers (RDKit `Params.cpp`,
//! which reproduce Rappé et al., *J. Am. Chem. Soc.* **1992**, 114, 10024).
//! The lookup helpers implement the paper's formulas for equilibrium bond
//! lengths (eq. 2–4), bond force constants (eq. 6), angle force constants
//! (eq. 13), torsion barriers (eq. 16/17) and inversion coefficients.
//!
//! Lengths are in Å, angles in degrees in the tables (converted to radians by
//! the builder), and energies in kcal/mol.

/// UFF atom type parameters.
/// Source: Rappé, A. K. et al. "UFF, a full periodic table force field for molecular mechanics
/// and molecular dynamics simulations." J. Am. Chem. Soc. 114 (1992): 10024–10035. Table 1.
///
/// One instance describes one UFF atom type (e.g. `C_3`, `N_2`, `Pt4+2`).
#[derive(Debug, Clone)]
pub struct UffAtomParams {
    /// UFF atom-type label, e.g. `"C_3"` or `"H_"`.
    pub label: &'static str,
    /// Atomic number of the element this type belongs to.
    pub element: u8,
    /// Bond radius `r_1` in Å — used in the equilibrium bond-length formula.
    pub r1: f64,
    /// Equilibrium valence angle `θ_0` in **degrees** (convert to radians before use).
    pub theta0: f64,
    /// van der Waals distance parameter `x_1` in Å.
    pub x1: f64,
    /// van der Waals well depth `D_1` in kcal/mol.
    pub d1: f64,
    /// Geometric-mean electronegativity scale factor `ζ` (dimensionless).
    pub zeta: f64,
    /// Effective nuclear charge `Z*` used in the bond force-constant formula (dimensionless).
    pub z_star: f64,
    /// Geometric-mean electronegativity `χ` used in the bond-length correction (dimensionless).
    pub chi: f64,
    /// Torsion barrier `V_i` in kcal/mol (UFF eq. 16).
    pub v_tors: f64,
    /// Torsion contribution `U_j` in kcal/mol for sp²–sp³ bonds (UFF eq. 17).
    pub u1: f64,
}

/// Retrieve UFF parameters by atom type label.
/// Canonical UFF atom parameters (RDKit `Params.cpp` values).
pub fn get_uff_params(label: &str) -> Option<UffAtomParams> {
    macro_rules! uff_atom {
        ($label:expr, $z:expr, $r1:expr, $theta0:expr, $x1:expr, $d1:expr, $zeta:expr, $z_star:expr, $v1:expr, $u1:expr, $xi:expr) => {
            Some(UffAtomParams {
                label: $label,
                element: $z,
                r1: $r1,
                theta0: $theta0,
                x1: $x1,
                d1: $d1,
                zeta: $zeta,
                z_star: $z_star,
                chi: $xi,
                v_tors: $v1,
                u1: $u1,
            })
        };
    }
    match label {
        "H_" => {
            uff_atom!("H_", 1, 0.354, 180.0f64, 2.886, 0.044, 12.0f64, 0.712, 0.0f64, 0.0f64, 4.528)
        }
        "C_3" => {
            uff_atom!("C_3", 6, 0.757, 109.47, 3.851, 0.105, 12.73, 1.912, 2.119, 2.0f64, 5.343)
        }
        "C_2" => {
            uff_atom!("C_2", 6, 0.732, 120.0f64, 3.851, 0.105, 12.73, 1.912, 0.0f64, 2.0f64, 5.343)
        }
        "C_R" => {
            uff_atom!("C_R", 6, 0.729, 120.0f64, 3.851, 0.105, 12.73, 1.912, 0.0f64, 2.0f64, 5.343)
        }
        "C_1" => {
            uff_atom!("C_1", 6, 0.706, 180.0f64, 3.851, 0.105, 12.73, 1.912, 0.0f64, 2.0f64, 5.343)
        }
        "N_3" => uff_atom!("N_3", 7, 0.7, 106.7, 3.66, 0.069, 13.407, 2.544, 0.45, 2.0f64, 6.899),
        "N_2" => {
            uff_atom!("N_2", 7, 0.685, 111.2, 3.66, 0.069, 13.407, 2.544, 0.0f64, 2.0f64, 6.899)
        }
        "N_R" => {
            uff_atom!("N_R", 7, 0.699, 120.0f64, 3.66, 0.069, 13.407, 2.544, 0.0f64, 2.0f64, 6.899)
        }
        "N_1" => {
            uff_atom!("N_1", 7, 0.656, 180.0f64, 3.66, 0.069, 13.407, 2.544, 0.0f64, 2.0f64, 6.899)
        }
        "O_3" => uff_atom!("O_3", 8, 0.658, 104.51, 3.5, 0.06, 14.085, 2.3, 0.018, 2.0f64, 8.741),
        "O_2" => {
            uff_atom!("O_2", 8, 0.634, 120.0f64, 3.5, 0.06, 14.085, 2.3, 0.0f64, 2.0f64, 8.741)
        }
        "O_R" => uff_atom!("O_R", 8, 0.68, 110.0f64, 3.5, 0.06, 14.085, 2.3, 0.0f64, 2.0f64, 8.741),
        "O_1" => {
            uff_atom!("O_1", 8, 0.639, 180.0f64, 3.5, 0.06, 14.085, 2.3, 0.0f64, 2.0f64, 8.741)
        }
        "F_" => {
            uff_atom!("F_", 9, 0.668, 180.0f64, 3.364, 0.05, 14.762, 1.735, 0.0f64, 2.0f64, 10.874)
        }
        "Cl" => {
            uff_atom!("Cl", 17, 1.044, 180.0f64, 3.947, 0.227, 14.866, 2.348, 0.0f64, 1.25, 8.564)
        }
        "Br" => {
            uff_atom!("Br", 35, 1.192, 180.0f64, 4.189, 0.251, 15.0f64, 2.519, 0.0f64, 0.7, 7.79)
        }
        "I_" => uff_atom!("I_", 53, 1.382, 180.0f64, 4.5, 0.339, 15.0f64, 2.65, 0.0f64, 0.2, 6.822),
        "S_3" => {
            uff_atom!("S_3", 16, 1.02, 92.1, 4.035, 0.274, 2.703, 2.703, 1.948, 0.0f64, 10.0f64)
        }
        "S_2" => {
            uff_atom!("S_2", 16, 0.854, 120.0f64, 4.035, 0.274, 13.969, 2.703, 0.0f64, 1.25, 6.928)
        }
        "S_R" => {
            uff_atom!("S_R", 16, 1.077, 92.2, 4.035, 0.274, 13.969, 2.703, 0.0f64, 1.25, 6.928)
        }
        "P_3" => {
            uff_atom!("P_3", 15, 0.93, 93.8, 4.147, 0.305, 2.862, 2.862, 2.421, 0.0f64, 8.0f64)
        }
        "B_2" => {
            uff_atom!("B_2", 5, 0.828, 120.0f64, 4.083, 0.18, 12.052, 1.755, 0.0f64, 2.0f64, 5.11)
        }
        "Si3" => {
            uff_atom!("Si3", 14, 1.117, 109.47, 4.295, 0.402, 12.175, 2.323, 1.225, 1.25, 4.168)
        }
        "Ti3+4" => {
            uff_atom!("Ti3+4", 22, 1.412, 109.47, 3.175, 0.017, 2.659, 2.659, 0.0f64, 0.0f64, 3.47)
        }
        "V_3+5" => {
            uff_atom!("V_3+5", 23, 1.402, 109.47, 3.144, 0.016, 2.679, 2.679, 0.0f64, 0.0f64, 3.65)
        }
        "Cr6+3" => uff_atom!(
            "Cr6+3", 24, 1.345, 90.0f64, 3.023, 0.015, 2.463, 2.463, 0.0f64, 0.0f64, 3.415
        ),
        "Mn6+2" => {
            uff_atom!("Mn6+2", 25, 1.382, 90.0f64, 2.961, 0.013, 2.43, 2.43, 0.0f64, 0.0f64, 3.325)
        }
        "Fe3+2" | "Fe6+2" => {
            uff_atom!("Fe3+2", 26, 1.335, 109.47, 2.912, 0.013, 2.43, 2.43, 0.0f64, 0.0f64, 3.76)
        }
        "Co6+3" => {
            uff_atom!("Co6+3", 27, 1.241, 90.0f64, 2.872, 0.014, 2.43, 2.43, 0.0f64, 0.0f64, 4.105)
        }
        "Ni4+2" => {
            uff_atom!("Ni4+2", 28, 1.164, 90.0f64, 2.834, 0.015, 2.43, 2.43, 0.0f64, 0.0f64, 4.465)
        }
        "Cu3+1" => {
            uff_atom!("Cu3+1", 29, 1.302, 109.47, 3.495, 0.005, 1.756, 1.756, 0.0f64, 0.0f64, 4.2)
        }
        "Zn3+2" => {
            uff_atom!("Zn3+2", 30, 1.193, 109.47, 2.763, 0.124, 1.308, 1.308, 0.0f64, 0.0f64, 5.106)
        }
        "Ge3" => uff_atom!("Ge3", 32, 1.197, 109.47, 4.28, 0.379, 2.789, 2.789, 0.701, 0.7, 4.051),
        "As3" => uff_atom!("As3", 33, 1.211, 92.1, 4.23, 0.309, 2.864, 2.864, 1.476, 0.7, 5.188),
        "Se3" => uff_atom!("Se3", 34, 1.19, 90.6, 4.205, 0.291, 2.764, 2.764, 0.335, 0.7, 6.428),
        "Sn3" => {
            uff_atom!("Sn3", 50, 1.399, 109.47, 4.392, 0.567, 3.198, 3.198, 0.199, 0.2, 4.0f64)
        }
        "Sb3" => uff_atom!("Sb3", 51, 1.404, 91.6, 4.42, 0.449, 3.51, 3.51, 1.1, 0.2, 4.899),
        "Te3" => uff_atom!("Te3", 52, 1.386, 90.25, 4.47, 0.398, 3.526, 3.526, 0.3, 0.2, 5.816),
        "Mo3+6" => {
            uff_atom!("Mo3+6", 42, 1.386, 109.47, 3.052, 0.056, 3.056, 3.056, 0.0f64, 0.0f64, 3.465)
        }
        "Pd4+2" => {
            uff_atom!("Pd4+2", 46, 1.338, 90.0f64, 2.899, 0.048, 2.899, 2.899, 0.0f64, 0.0f64, 4.32)
        }
        "Ag1+1" => uff_atom!(
            "Ag1+1", 47, 1.386, 180.0f64, 3.148, 0.036, 1.956, 1.956, 0.0f64, 0.0f64, 4.436
        ),
        "Pt4+2" => {
            uff_atom!("Pt4+2", 78, 1.364, 90.0f64, 2.754, 0.08, 2.754, 2.754, 0.0f64, 0.0f64, 4.79)
        }
        "Au4+3" => {
            uff_atom!("Au4+3", 79, 1.255, 90.0f64, 3.293, 0.039, 2.309, 2.309, 0.0f64, 0.0f64, 5.77)
        }
        "B_3" => uff_atom!("B_3", 5, 0.838, 109.47, 4.083, 0.18, 12.052, 1.755, 0.0, 2.0, 5.11),
        "P_3+3" => {
            uff_atom!("P_3+3", 15, 1.101, 93.8, 4.147, 0.305, 13.072, 2.863, 2.4, 1.25, 5.463)
        }
        "P_3+5" => {
            uff_atom!("P_3+5", 15, 1.056, 109.47, 4.147, 0.305, 13.072, 2.863, 2.4, 1.25, 5.463)
        }
        "S_3+2" => {
            uff_atom!("S_3+2", 16, 1.064, 92.1, 4.035, 0.274, 13.969, 2.703, 0.484, 1.25, 6.928)
        }
        "S_3+4" => {
            uff_atom!("S_3+4", 16, 1.049, 103.2, 4.035, 0.274, 13.969, 2.703, 0.484, 1.25, 6.928)
        }
        "S_3+6" => {
            uff_atom!("S_3+6", 16, 1.027, 109.47, 4.035, 0.274, 13.969, 2.703, 0.484, 1.25, 6.928)
        }
        _ => None,
    }
}

/// Compute the UFF equilibrium bond length between two atom types.
///
/// Rappé eq. 2–4:
/// r_ij = r_i + r_j + r_BO - r_EN
/// r_BO = -λ · (r_i + r_j) · ln(BO)          [λ = 0.1332]
/// r_EN = r_i · r_j · (√χ_i - √χ_j)² / (χ_i·r_i + χ_j·r_j)
///
/// The returned length is in **Å** and is clamped to a minimum of 0.5 Å.
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
///
/// `r0` is the equilibrium bond length in Å; the returned force constant is in
/// **kcal/mol/Å²**.
pub fn get_uff_bond_force_constant(pi: &UffAtomParams, pj: &UffAtomParams, r0: f64) -> f64 {
    664.12 * pi.z_star * pj.z_star / r0.powi(3)
}

/// UFF angle bending force constant K_ijk (Rappé eq. 13, simplified).
/// K_ijk = 332.06 · Z*_i · Z*_k / (r_ij · r_jk · r_ik²)
/// 332.06 ≈ e²Nₐ/(4πε₀) in kcal·Å/e² (Coulomb constant for UFF energy units)
/// where r_ik² = r_ij² + r_jk² - 2·r_ij·r_jk·cos(θ₀)
///
/// `r_ij` and `r_jk` are bond lengths in Å and `theta0_rad` the equilibrium
/// angle in radians; the returned force constant is in **kcal/mol/rad²**.
pub fn get_uff_angle_force_constant(
    pi: &UffAtomParams,
    pk: &UffAtomParams,
    r_ij: f64,
    r_jk: f64,
    theta0_rad: f64,
) -> f64 {
    // Canonical RDKit formula (`AngleBend.cpp::calcAngleForceConstant`):
    //   r13    = sqrt(r12^2 + r23^2 - 2 r12 r23 cos(theta0))
    //   beta   = 2G / (r12 r23)
    //   ka     = beta * Z1_i * Z1_k / r13^5 * r12 r23 * inner
    //   inner  = 3 r12 r23 (1 - cos^2) - r13^2 cos
    const G: f64 = 332.06;
    if r_ij < 1e-10 || r_jk < 1e-10 {
        return 0.0;
    }
    let cos0 = theta0_rad.cos();
    let r13_sq = r_ij * r_ij + r_jk * r_jk - 2.0 * r_ij * r_jk * cos0;
    if r13_sq < 1e-12 {
        return 0.0;
    }
    let r13 = r13_sq.sqrt();
    let beta = 2.0 * G / (r_ij * r_jk);
    let pre = beta * pi.z_star * pk.z_star / r13.powi(5);
    let inner = 3.0 * r_ij * r_jk * (1.0 - cos0 * cos0) - r13_sq * cos0;
    pre * r_ij * r_jk * inner
}

/// True for the UFF "group 6" elements used by the torsion special cases.
///
/// Group 6 in the UFF paper is O, S, Se, Te and Po; the special case applies
/// to single bonds between two of them.
pub fn is_uff_group6(z: u8) -> bool {
    matches!(z, 8 | 16 | 34 | 52 | 84)
}

/// UFF equation 17: sp2-sp2 torsion barrier.
///
/// `bond_order` is the numeric bond order of the central bond and `uj`/`uk`
/// the `U_j` torsion parameters of the two atom types (kcal/mol). The
/// returned barrier is in **kcal/mol**.
pub fn uff_equation17(bond_order: f64, uj: f64, uk: f64) -> f64 {
    5.0 * (uj * uk).sqrt() * (1.0 + 4.18 * bond_order.ln())
}

/// UFF torsion parameters for the central bond j-k.
///
/// Returns `(V, periodicity_n, cos_term)` such that
/// `E = V/2 * (1 - cos_term * cos(n φ))`, exactly as the canonical
/// implementation (`TorsionAngle.cpp::calcTorsionParams`).
///
/// `bond_order` is the numeric order of the j-k bond (1, 2, 3 or 1.5 for
/// aromatic); `end_atom_is_sp2` is true when either terminal atom of the
/// four-atom torsion is sp2.
///
/// Returns `(V, periodicity_n, cos_term)` where `V` is the barrier in
/// kcal/mol, `n` the (dimensionless) periodicity and `cos_term` the
/// (dimensionless) phase factor `cos(φ₀)`.
pub fn get_uff_torsion_params(
    hyb_j: &crate::graph::Hybridization,
    hyb_k: &crate::graph::Hybridization,
    at_num_j: u8,
    at_num_k: u8,
    pj: &UffAtomParams,
    pk: &UffAtomParams,
    bond_order: f64,
    end_atom_is_sp2: bool,
) -> (f64, f64, f64) {
    use crate::graph::Hybridization::*;

    let both_sp3 = matches!(hyb_j, SP3) && matches!(hyb_k, SP3);
    let both_sp2 = matches!(hyb_j, SP2) && matches!(hyb_k, SP2);
    if both_sp3 {
        let mut v = (pj.v_tors * pk.v_tors).sqrt();
        let mut n = 3.0;
        let mut cos_term = -1.0;
        // Special case: single bond between two group-6 elements.
        if (bond_order - 1.0).abs() < 1e-9 && is_uff_group6(at_num_j) && is_uff_group6(at_num_k) {
            let v2: f64 = if at_num_j == 8 { 2.0 } else { 6.8 };
            let v3: f64 = if at_num_k == 8 { 2.0 } else { 6.8 };
            v = (v2 * v3).sqrt();
            n = 2.0;
            cos_term = -1.0;
        }
        (v, n, cos_term)
    } else if both_sp2 {
        (uff_equation17(bond_order, pj.u1, pk.u1), 2.0, 1.0)
    } else {
        // sp2-sp3: default V = 1, periodicity 6, phi0 = 0
        let mut v = 1.0;
        let mut n = 6.0;
        let mut cos_term = 1.0;
        if (bond_order - 1.0).abs() < 1e-9 {
            let j_is_sp3_g6 = matches!(hyb_j, SP3) && is_uff_group6(at_num_j);
            let k_is_sp3_g6 = matches!(hyb_k, SP3) && is_uff_group6(at_num_k);
            if (j_is_sp3_g6 && !is_uff_group6(at_num_k))
                || (k_is_sp3_g6 && !is_uff_group6(at_num_j))
            {
                v = uff_equation17(bond_order, pj.u1, pk.u1);
                n = 2.0;
                cos_term = -1.0;
            } else if end_atom_is_sp2 {
                v = 2.0;
                n = 3.0;
                cos_term = -1.0;
            }
        }
        (v, n, cos_term)
    }
}

/// UFF inversion coefficients and per-contribution force constant.
///
/// Returns `(K_per_contrib, C0, C1, C2)` with `K` already divided by three,
/// matching `Utils::calcInversionCoefficientsAndForceConstant` used by the
/// three inversion contributions per sp2 centre.
///
/// `K_per_contrib` is in kcal/mol/rad²; `C0`, `C1` and `C2` are dimensionless.
/// `(0.0, 0.0, 0.0, 0.0)` is returned for elements without a UFF inversion
/// parameterisation.
pub fn get_uff_inversion_params(at2_atomic_num: u8, is_c_bound_to_o: bool) -> (f64, f64, f64, f64) {
    use std::f64::consts::PI;
    if matches!(at2_atomic_num, 6..=8) {
        let k = if is_c_bound_to_o { 50.0 } else { 6.0 };
        (k / 3.0, 1.0, -1.0, 0.0)
    } else {
        let w0 = match at2_atomic_num {
            15 => 84.4339_f64,
            33 => 86.9735,
            51 => 87.7047,
            83 => 90.0,
            _ => return (0.0, 0.0, 0.0, 0.0),
        } * PI
            / 180.0;
        let c2 = 1.0;
        let c1 = -4.0 * w0.cos();
        let c0 = -(c1 * w0.cos() + c2 * (2.0 * w0).cos());
        let k = 22.0 / (c0 + c1 + c2);
        (k / 3.0, c0, c1, c2)
    }
}
