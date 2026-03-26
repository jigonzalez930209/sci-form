//! Molecular dipole moment calculation.
//!
//! μ = μ_nuclear + μ_electronic
//!
//! μ_nuclear = Σ_A Z_A * R_A   (valence charges × positions)
//! μ_electronic = -Σ_{μν} P_{μν} <μ|r|ν>
//!
//! For EHT with localized Gaussian basis, the electronic dipole is
//! approximated using Mulliken-partitioned charges:
//!   μ_electronic ≈ -Σ_{μ} (PS)_{μμ} * R_{atom(μ)}
//!
//! This gives:  μ_A = q_A * R_A  where q_A = Z_A - Σ_{μ∈A} (PS)_{μμ}
//! i.e. μ = Σ_A q_Mulliken_A * R_A  (with sign convention)

use serde::{Deserialize, Serialize};

/// Conversion factor: e·Å → Debye (1 e·Å = 4.8032 D)
const EANG_TO_DEBYE: f64 = 4.80321;

/// Result of a dipole moment calculation.
///
/// The total dipole μ = μ_nuclear + μ_electronic.
/// `nuclear_dipole` and `electronic_dipole` are `Some` when computed from
/// an electronic-structure method (EHT); they are `None` when computed
/// directly from pre-built Mulliken charges.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DipoleResult {
    /// Dipole vector components (Debye).
    pub vector: [f64; 3],
    /// Dipole magnitude (Debye).
    pub magnitude: f64,
    /// Unit: always "Debye".
    pub unit: String,
    /// Nuclear contribution to the dipole (Debye), if decomposition is available.
    pub nuclear_dipole: Option<[f64; 3]>,
    /// Electronic contribution to the dipole (Debye), if decomposition is available.
    pub electronic_dipole: Option<[f64; 3]>,
}

/// Compute the molecular dipole moment from Mulliken charges and positions.
///
/// μ = Σ_A q_A * R_A  (in e·Å, then converted to Debye)
///
/// - `mulliken_charges`: per-atom partial charges from population analysis
/// - `positions`: atom positions in Ångström
pub fn compute_dipole(mulliken_charges: &[f64], positions: &[[f64; 3]]) -> DipoleResult {
    let n = mulliken_charges.len();
    let mut mu = [0.0, 0.0, 0.0];

    // Dipole = Σ q_A * R_A  (charge × position)
    // Note: positive charge at position R contributes +q*R
    for i in 0..n {
        let q = mulliken_charges[i];
        mu[0] += q * positions[i][0];
        mu[1] += q * positions[i][1];
        mu[2] += q * positions[i][2];
    }

    // Convert e·Å to Debye
    mu[0] *= EANG_TO_DEBYE;
    mu[1] *= EANG_TO_DEBYE;
    mu[2] *= EANG_TO_DEBYE;

    let magnitude = (mu[0] * mu[0] + mu[1] * mu[1] + mu[2] * mu[2]).sqrt();

    DipoleResult {
        vector: mu,
        magnitude,
        unit: "Debye".to_string(),
        nuclear_dipole: None,
        electronic_dipole: None,
    }
}

/// Compute dipole directly from EHT results + molecule geometry.
///
/// Runs Mulliken population analysis internally.
pub fn compute_dipole_from_eht(
    elements: &[u8],
    positions: &[[f64; 3]],
    coefficients: &[Vec<f64>],
    n_electrons: usize,
) -> DipoleResult {
    let basis = crate::eht::basis::build_basis(elements, positions);
    let overlap = crate::eht::build_overlap_matrix(&basis);
    let charges =
        crate::population::mulliken_charges(elements, &basis, &overlap, coefficients, n_electrons);

    // Decompose into nuclear and electronic contributions.
    // Nuclear: μ_nuc = Σ_A Z_A * R_A  (valence Z — nuclear charges from EHT basis).
    // Electronic: μ_el = -Σ_A pop_A * R_A  (Mulliken populations).
    let mut mu_nuc = [0.0; 3];
    let mut mu_el = [0.0; 3];
    for (i, &z) in elements.iter().enumerate() {
        let z_val = crate::population::population::valence_electrons(z);
        let pop = z_val - charges[i]; // population = Z_val - q
        for k in 0..3 {
            mu_nuc[k] += z_val * positions[i][k];
            mu_el[k] -= pop * positions[i][k];
        }
    }
    for k in 0..3 {
        mu_nuc[k] *= EANG_TO_DEBYE;
        mu_el[k] *= EANG_TO_DEBYE;
    }

    let mut result = compute_dipole(&charges, positions);
    result.nuclear_dipole = Some(mu_nuc);
    result.electronic_dipole = Some(mu_el);
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eht::solve_eht;

    #[test]
    fn test_h2_zero_dipole() {
        // Symmetric H₂ → dipole ≈ 0
        let elems = vec![1u8, 1];
        let pos = vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let result = solve_eht(&elems, &pos, None).unwrap();
        let dipole =
            compute_dipole_from_eht(&elems, &pos, &result.coefficients, result.n_electrons);

        assert!(
            dipole.magnitude < 0.1,
            "H₂ should be nonpolar, got {:.3} D",
            dipole.magnitude
        );
    }

    #[test]
    fn test_water_nonzero_dipole() {
        let elems = vec![8u8, 1, 1];
        let pos = vec![[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let result = solve_eht(&elems, &pos, None).unwrap();
        let dipole =
            compute_dipole_from_eht(&elems, &pos, &result.coefficients, result.n_electrons);

        // Water has a significant dipole (EHT may not match ab initio exactly,
        // but should be clearly nonzero)
        assert!(
            dipole.magnitude > 0.1,
            "Water should be polar, got {:.3} D",
            dipole.magnitude
        );
    }

    #[test]
    fn test_dipole_vector_components() {
        let elems = vec![8u8, 1, 1];
        let pos = vec![[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let result = solve_eht(&elems, &pos, None).unwrap();
        let dipole =
            compute_dipole_from_eht(&elems, &pos, &result.coefficients, result.n_electrons);

        // z component should be ~0 (planar molecule in xy)
        assert!(
            dipole.vector[2].abs() < 0.1,
            "z-component should be ~0 for planar molecule"
        );
        // x component should be ~0 (symmetric about yz)
        assert!(
            dipole.vector[0].abs() < 0.1,
            "x-component should be ~0 for symmetric water"
        );
    }

    #[test]
    fn test_dipole_units_debye() {
        let charges = vec![0.5, -0.5];
        let pos = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let dipole = compute_dipole(&charges, &pos);
        assert_eq!(dipole.unit, "Debye");
        // 0.5 e × 1 Å = 0.5 × 4.803 = ~2.4 D
        assert!(
            (dipole.magnitude - 2.4016).abs() < 0.01,
            "Expected ~2.4 D, got {:.3}",
            dipole.magnitude
        );
    }
}
