//! ReaxFF parameter set and parsing.
// The tables below transcribe published parameter values verbatim; the
// decimal literals are kept as published rather than replaced by symbolic
// constants, so `clippy::approx_constant` does not apply here.
#![allow(clippy::approx_constant)]

use serde::{Deserialize, Serialize};

/// Parameters for a single atom type in ReaxFF.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReaxffAtomParams {
    /// Atomic number.
    pub element: u8,
    /// σ bond equilibrium distance (Å).
    pub r_sigma: f64,
    /// π bond equilibrium distance (Å).
    pub r_pi: f64,
    /// ππ bond equilibrium distance (Å).
    pub r_pipi: f64,
    /// Bond order parameters.
    ///
    /// The six `p_bo` values are the empirical exponents of the ReaxFF
    /// bond-order function (dimensionless); they are published in the
    /// `ffield.reaxff` files of van Duin et al., *J. Phys. Chem. A* **2001**,
    /// 105, 9393.
    pub p_bo1: f64,
    /// Second bond-order exponent (dimensionless).
    pub p_bo2: f64,
    /// Third bond-order exponent (dimensionless).
    pub p_bo3: f64,
    /// Fourth bond-order exponent (dimensionless).
    pub p_bo4: f64,
    /// Fifth bond-order exponent (dimensionless).
    pub p_bo5: f64,
    /// Sixth bond-order exponent (dimensionless).
    pub p_bo6: f64,
    /// Valence (nominal number of bonds).
    pub valence: f64,
}

/// Full ReaxFF parameter set.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReaxffParams {
    /// Per-atom parameters.
    pub atom_params: Vec<ReaxffAtomParams>,
    /// Bond dissociation energy parameters.
    pub bond_de: Vec<Vec<f64>>,
    /// Bond energy parameters.
    pub p_be1: f64,
    /// Bond energy exponent (global fallback; per-pair values in
    /// [`ReaxffParams::pair_table`] take precedence).
    pub p_be2: f64,
    /// Over-coordination penalty (kcal/mol).
    pub p_ovun1: f64,
    /// Over-coordination threshold (dimensionless bond-order sum).
    pub p_val3: f64,
    /// Under-coordination penalty (kcal/mol).
    pub p_ovun5: f64,
    /// Equilibrium angles (radians), indexed as `[i][j][k]`.
    pub equilibrium_angles: Vec<f64>,
    /// Angle force constants.
    pub angle_force_constants: Vec<f64>,
    /// Non-bonded cutoff (Å).
    pub cutoff: f64,
    /// Per-element-pair σ-bond data transcribed from `ffield.reax.cho`
    /// (Chenoweth/van Duin/Goddard 2008, C/H/O combustion). Takes precedence
    /// over the atom-averaged rules for the covered pairs.
    pub pair_table: Vec<ReaxffPairParams>,
}

/// Per-pair σ-bond parameters from a ReaxFF force field file.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ReaxffPairParams {
    /// First element (atomic number).
    pub z1: u8,
    /// Second element (atomic number).
    pub z2: u8,
    /// Pair σ equilibrium distance (Å; off-diagonal value for hetero pairs).
    pub r_sigma: f64,
    /// σ bond-order exponent prefactor (`p_bo1`).
    pub p_bo1: f64,
    /// σ bond-order distance exponent (`p_bo2`).
    pub p_bo2: f64,
    /// σ dissociation energy (kcal/mol, `De_sigma`).
    pub de_sigma: f64,
    /// Bond-energy screening (`p_be1`).
    pub p_be1: f64,
    /// Bond-energy screening exponent (`p_be2`).
    pub p_be2: f64,
}

impl ReaxffParams {
    /// Build default C/H/O/N parameters for ReaxFF.
    pub fn default_chon() -> Self {
        let atom_params = vec![
            // H
            ReaxffAtomParams {
                element: 1,
                r_sigma: 0.656,
                r_pi: 0.0,
                r_pipi: 0.0,
                p_bo1: -0.0500,
                p_bo2: 6.9136,
                p_bo3: 0.0,
                p_bo4: 6.0,
                p_bo5: 0.0,
                p_bo6: 6.0,
                valence: 1.0,
            },
            // C
            ReaxffAtomParams {
                element: 6,
                r_sigma: 1.3825,
                r_pi: 1.1359,
                r_pipi: 1.2104,
                p_bo1: -0.0777,
                p_bo2: 6.7268,
                p_bo3: -0.1000,
                p_bo4: 9.1628,
                p_bo5: -0.1418,
                p_bo6: 13.3056,
                valence: 4.0,
            },
            // N
            ReaxffAtomParams {
                element: 7,
                r_sigma: 1.2380,
                r_pi: 1.1748,
                r_pipi: 1.0630,
                p_bo1: -0.1000,
                p_bo2: 6.8773,
                p_bo3: -0.1193,
                p_bo4: 7.8431,
                p_bo5: -0.1418,
                p_bo6: 13.1260,
                valence: 3.0,
            },
            // O
            ReaxffAtomParams {
                element: 8,
                r_sigma: 1.2477,
                r_pi: 1.0863,
                r_pipi: 0.0,
                p_bo1: -0.1363,
                p_bo2: 5.6346,
                p_bo3: -0.1743,
                p_bo4: 7.6279,
                p_bo5: 0.0,
                p_bo6: 6.0,
                valence: 2.0,
            },
        ];

        // Build bond DE matrix (simplified)
        let n = atom_params.len();
        let mut bond_de = vec![vec![100.0; n]; n]; // kcal/mol
                                                   // H-H, C-C, etc.
        bond_de[0][0] = 104.2; // H-H
        bond_de[1][1] = 145.0; // C-C
        bond_de[0][1] = 99.0; // H-C
        bond_de[1][0] = 99.0;
        bond_de[0][3] = 111.0; // H-O
        bond_de[3][0] = 111.0;
        bond_de[1][3] = 85.0; // C-O
        bond_de[3][1] = 85.0;
        bond_de[1][2] = 73.0; // C-N
        bond_de[2][1] = 73.0;

        Self {
            atom_params,
            bond_de,
            p_be1: -0.2,
            p_be2: 6.25,
            p_ovun1: 50.0,
            p_val3: 3.0,
            p_ovun5: 10.0,
            equilibrium_angles: vec![
                std::f64::consts::FRAC_PI_2 * 2.0 / 3.0 * std::f64::consts::PI;
                100
            ],
            angle_force_constants: vec![50.0; 100],
            cutoff: 10.0,
            pair_table: ffield_cho_pairs(),
        }
    }

    /// Look up per-pair σ data for elements (`z1`, `z2`), order-independent.
    pub fn pair_params(&self, z1: u8, z2: u8) -> Option<&ReaxffPairParams> {
        self.pair_table
            .iter()
            .find(|p| (p.z1 == z1 && p.z2 == z2) || (p.z1 == z2 && p.z2 == z1))
    }

    /// Get bond dissociation energy for pair (i, j) by atom index.
    pub fn get_bond_de(&self, i: usize, j: usize) -> f64 {
        if i < self.bond_de.len() && j < self.bond_de[i].len() {
            self.bond_de[i][j]
        } else {
            100.0 // default
        }
    }

    /// Get equilibrium angle for triplet i-j-k.
    pub fn get_equilibrium_angle(&self, _i: usize, _j: usize, _k: usize) -> f64 {
        109.47_f64.to_radians() // tetrahedral default
    }

    /// Get angle force constant for triplet i-j-k.
    pub fn get_angle_force_constant(&self, _i: usize, _j: usize, _k: usize) -> f64 {
        50.0 // kcal/mol/rad²
    }

    /// Map element to parameter index.
    pub fn element_index(&self, z: u8) -> Option<usize> {
        self.atom_params.iter().position(|a| a.element == z)
    }
}

/// C/H/O σ-pair data transcribed from `ffield.reax.cho`
/// (Chenoweth/van Duin/Goddard, J. Phys. Chem. A 112, 1040–1053 (2008)):
/// bond section (`De_sigma`, `p_be1`, `p_be2`, `p_bo1`, `p_bo2`) plus the
/// off-diagonal `ro(sigma)` for hetero pairs (homo pairs use atomic `ro`).
fn ffield_cho_pairs() -> Vec<ReaxffPairParams> {
    vec![
        ReaxffPairParams {
            z1: 6,
            z2: 6,
            r_sigma: 1.3825,
            p_bo1: -0.0750,
            p_bo2: 6.8316,
            de_sigma: 156.5953,
            p_be1: -0.8157,
            p_be2: 0.4527,
        },
        ReaxffPairParams {
            z1: 6,
            z2: 1,
            r_sigma: 1.1203,
            p_bo1: -0.0500,
            p_bo2: 6.8315,
            de_sigma: 170.2316,
            p_be1: -0.5931,
            p_be2: 5.2267,
        },
        ReaxffPairParams {
            z1: 1,
            z2: 1,
            r_sigma: 0.7853,
            p_bo1: -0.0593,
            p_bo2: 4.8358,
            de_sigma: 156.0973,
            p_be1: -0.1377,
            p_be2: 2.9907,
        },
        ReaxffPairParams {
            z1: 6,
            z2: 8,
            r_sigma: 1.2775,
            p_bo1: -0.1463,
            p_bo2: 5.2913,
            de_sigma: 160.4802,
            p_be1: -0.3873,
            p_be2: 0.5341,
        },
        ReaxffPairParams {
            z1: 8,
            z2: 8,
            r_sigma: 1.2477,
            p_bo1: -0.1302,
            p_bo2: 6.2919,
            de_sigma: 60.1463,
            p_be1: -0.2802,
            p_be2: 0.2441,
        },
        ReaxffPairParams {
            z1: 1,
            z2: 8,
            r_sigma: 0.9013,
            p_bo1: -0.0657,
            p_bo2: 5.0451,
            de_sigma: 180.4373,
            p_be1: -0.8074,
            p_be2: 1.2490,
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_chon_has_four_elements() {
        let params = ReaxffParams::default_chon();
        assert_eq!(params.atom_params.len(), 4);
    }

    #[test]
    fn element_index_finds_carbon() {
        let params = ReaxffParams::default_chon();
        assert!(params.element_index(6).is_some(), "carbon should be found");
    }

    #[test]
    fn element_index_none_for_missing() {
        let params = ReaxffParams::default_chon();
        assert!(
            params.element_index(79).is_none(),
            "gold should not be in CHON"
        );
    }

    #[test]
    fn bond_de_is_positive() {
        let params = ReaxffParams::default_chon();
        let de = params.get_bond_de(0, 1);
        assert!(de > 0.0, "bond dissociation energy should be positive");
    }

    #[test]
    fn equilibrium_angle_is_tetrahedral() {
        let params = ReaxffParams::default_chon();
        let angle = params.get_equilibrium_angle(0, 1, 2);
        let tetra = 109.47_f64.to_radians();
        assert!((angle - tetra).abs() < 0.01);
    }
}
