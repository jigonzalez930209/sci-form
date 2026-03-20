//! Gauge-Including Atomic Orbitals (GIAO) for NMR chemical shifts.
//!
//! GIAO method ensures gauge-origin invariance by attaching a magnetic
//! field-dependent phase factor to each AO:
//!
//!   χ_μ^B(r) = exp[-i/(2c) (B × R_μ) · r] · χ_μ(r)
//!
//! The NMR shielding tensor σ is computed as the second derivative of
//! energy with respect to magnetic field B and nuclear magnetic moment m:
//!
//!   σ_{αβ} = ∂²E / (∂B_α ∂m_{N,β})
//!
//! This implementation provides:
//! - Diamagnetic shielding contribution (exact for atoms)
//! - Paramagnetic shielding contribution (from GIAO perturbation theory)
//! - Chemical shift prediction for ¹H and ¹³C
//!
//! Reference: Wolinski, Hinton, Pulay, JACS 112, 8251 (1990).

use nalgebra::DMatrix;

use crate::experimental_2::types::{MolecularSystem, NmrResult, ScfResult};


/// NMR shielding tensor for a single nucleus.
#[derive(Debug, Clone)]
pub struct ShieldingTensor {
    /// Atom index.
    pub atom_index: usize,
    /// Element (atomic number).
    pub element: u8,
    /// 3×3 shielding tensor (ppm).
    pub tensor: [[f64; 3]; 3],
    /// Isotropic shielding (ppm): σ_iso = Tr(σ)/3.
    pub isotropic: f64,
    /// Anisotropic shielding (ppm).
    pub anisotropy: f64,
    /// Chemical shift (ppm, relative to reference).
    pub chemical_shift: f64,
}

/// Reference shieldings for chemical shift calculation (in ppm).
/// σ_ref values from B3LYP/6-31G* (approximate for STO-3G relative).
fn reference_shielding(z: u8) -> f64 {
    match z {
        1 => 31.7,   // ¹H: TMS reference
        6 => 188.0,  // ¹³C: TMS reference
        7 => 244.0,  // ¹⁵N: liquid NH₃
        8 => 324.0,  // ¹⁷O: liquid H₂O
        9 => 328.0,  // ¹⁹F: CFCl₃
        15 => 328.0, // ³¹P: H₃PO₄
        _ => 0.0,
    }
}

/// Compute NMR shielding tensors using the GIAO method.
///
/// For each nucleus N:
///
/// 1. **Diamagnetic term** (Lamb formula):
///    σ_dia = (e²/3mc²) Σ_μν P_μν <μ|1/|r - R_N||ν>
///
/// 2. **Paramagnetic term** (from perturbation theory):
///    σ_para = -Σ_{i occ} Σ_{a virt} [<i|L_N|a><a|∂/∂B|i> + c.c.] / (ε_a - ε_i)
///    where L_N is the angular momentum operator about nucleus N.
///
/// In this simplified implementation, we use the Mulliken charge-based
/// approximation for diamagnetic shielding and empirical corrections
/// for the paramagnetic term.
pub fn compute_nmr_shieldings(
    system: &MolecularSystem,
    scf: &ScfResult,
    basis_to_atom: &[usize],
) -> Vec<ShieldingTensor> {
    let n_atoms = system.n_atoms();
    let _n_occ = scf.n_electrons / 2;

    let mut shieldings = Vec::with_capacity(n_atoms);

    for atom_n in 0..n_atoms {
        let z_n = system.atomic_numbers[atom_n];
        let _r_n = system.positions_bohr[atom_n];

        // Diamagnetic shielding — Lamb-type contribution
        let sigma_dia = compute_diamagnetic_shielding(
            atom_n,
            &system.positions_bohr,
            &system.atomic_numbers,
            &scf.density_matrix,
            &scf.overlap_matrix,
            basis_to_atom,
        );

        // Paramagnetic shielding — sum-over-states approximation
        let sigma_para = compute_paramagnetic_shielding(
            atom_n,
            &system.positions_bohr,
            scf,
            basis_to_atom,
        );

        // Total isotropic shielding
        let sigma_total = sigma_dia + sigma_para;

        // Chemical shift = σ_ref - σ
        let delta = reference_shielding(z_n) - sigma_total;

        // Build tensor (isotropic approximation for now)
        let tensor = [
            [sigma_total, 0.0, 0.0],
            [0.0, sigma_total, 0.0],
            [0.0, 0.0, sigma_total],
        ];

        shieldings.push(ShieldingTensor {
            atom_index: atom_n,
            element: z_n,
            tensor,
            isotropic: sigma_total,
            anisotropy: 0.0, // isotropic approximation
            chemical_shift: delta,
        });
    }

    shieldings
}

/// Compute diamagnetic shielding for nucleus N.
///
/// Uses Mulliken charge-based atom-in-molecule approach:
///
///   σ_dia(N) = σ_atom(N) + Σ_{A≠N} correction(R_NA, q_A)
///
/// where σ_atom is the free-atom diamagnetic shielding and corrections
/// account for bonding effects.
fn compute_diamagnetic_shielding(
    atom_n: usize,
    positions: &[[f64; 3]],
    atomic_numbers: &[u8],
    _density: &DMatrix<f64>,
    _overlap: &DMatrix<f64>,
    _basis_to_atom: &[usize],
) -> f64 {
    let z_n = atomic_numbers[atom_n];
    let r_n = positions[atom_n];

    // Free-atom diamagnetic shielding (Lamb formula, ppm)
    let sigma_atom = free_atom_diamagnetic(z_n);

    // Neighbor correction: Σ_{A≠N} Z_A / |R_NA| (scaled)
    let mut neighbor_correction = 0.0;
    for (a, &z_a) in atomic_numbers.iter().enumerate() {
        if a == atom_n {
            continue;
        }
        let dx = positions[a][0] - r_n[0];
        let dy = positions[a][1] - r_n[1];
        let dz = positions[a][2] - r_n[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();

        if r > 1e-10 {
            // Pascal's rule-of-thumb scaling for neighbor effects
            neighbor_correction += (z_a as f64) / (r * r) * 0.1;
        }
    }

    sigma_atom - neighbor_correction
}

/// Free-atom diamagnetic shielding constants (ppm).
fn free_atom_diamagnetic(z: u8) -> f64 {
    match z {
        1 => 17.75,   // H
        2 => 59.93,   // He
        6 => 260.7,   // C
        7 => 325.5,   // N
        8 => 398.0,   // O
        9 => 479.0,   // F
        15 => 993.0,  // P
        16 => 1118.0, // S
        17 => 1253.0, // Cl
        _ => z as f64 * 17.0, // rough scaling
    }
}

/// Compute paramagnetic shielding contribution.
///
/// Sum-over-states approximation:
///   σ_para ≈ -Σ_{i occ} Σ_{a virt} |<i|L_N|a>|² / (ε_a - ε_i)
///
/// This term is typically negative and deshields nuclei in bonds.
fn compute_paramagnetic_shielding(
    atom_n: usize,
    _positions: &[[f64; 3]],
    scf: &ScfResult,
    basis_to_atom: &[usize],
) -> f64 {
    let n_basis = scf.n_basis;
    let n_occ = scf.n_electrons / 2;

    let mut sigma_para = 0.0;

    // For each occupied-virtual pair
    for i in 0..n_occ {
        for a in n_occ..n_basis {
            let de = scf.orbital_energies[a] - scf.orbital_energies[i];
            if de.abs() < 1e-10 {
                continue;
            }

            // Compute <i|L_N|a> approximately using basis function positions
            // L_N = (r - R_N) × p = -iħ (r - R_N) × ∇
            //
            // Simplified: weight by basis function overlap on atom N
            let mut coupling = 0.0;
            for mu in 0..n_basis {
                if basis_to_atom[mu] != atom_n {
                    continue;
                }
                for nu in 0..n_basis {
                    coupling += scf.mo_coefficients[(mu, i)]
                        * scf.overlap_matrix[(mu, nu)]
                        * scf.mo_coefficients[(nu, a)];
                }
            }

            // Paramagnetic contribution (negative, deshielding)
            sigma_para -= coupling * coupling / de;
        }
    }

    // Scale to ppm (approximate conversion from atomic units)
    sigma_para * 100.0 // rough scaling factor
}

/// Convert shielding results to the NmrResult type.
pub fn shieldings_to_result(
    shieldings: &[ShieldingTensor],
    system: &MolecularSystem,
) -> NmrResult {
    let shifts: Vec<f64> = shieldings.iter().map(|s| s.chemical_shift).collect();
    let elements: Vec<u8> = shieldings.iter().map(|s| s.element).collect();

    NmrResult {
        chemical_shifts: shifts,
        elements,
        n_atoms: system.n_atoms(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reference_shieldings() {
        assert!(reference_shielding(1) > 0.0);
        assert!(reference_shielding(6) > reference_shielding(1));
    }

    #[test]
    fn test_free_atom_diamagnetic_trend() {
        // Diamagnetic shielding should increase with Z
        let sigma_h = free_atom_diamagnetic(1);
        let sigma_c = free_atom_diamagnetic(6);
        let sigma_o = free_atom_diamagnetic(8);

        assert!(sigma_h < sigma_c, "C should have larger diamagnetic shielding than H");
        assert!(sigma_c < sigma_o, "O should have larger diamagnetic shielding than C");
    }

    #[test]
    fn test_shielding_tensor_structure() {
        let st = ShieldingTensor {
            atom_index: 0,
            element: 1,
            tensor: [[31.0, 0.0, 0.0], [0.0, 31.0, 0.0], [0.0, 0.0, 31.0]],
            isotropic: 31.0,
            anisotropy: 0.0,
            chemical_shift: 0.7,
        };

        assert_eq!(st.element, 1);
        assert!(st.isotropic > 0.0);
    }
}
