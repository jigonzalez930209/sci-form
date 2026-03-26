//! STO-3G minimal basis set with contracted Gaussian primitives.
//!
//! Implements the Hehre-Stewart-Pople STO-3G basis (1969) where each
//! Slater-type orbital is approximated by 3 Gaussian primitives:
//!
//!   χ_STO(r) ≈ Σ_{i=1}^{3} c_i · g(α_i, r)
//!
//! where g(α, r) = N · r^l · exp(-α·r²) is a Cartesian Gaussian primitive.
//!
//! # Supported Elements
//!
//! H(1), He(2), Li(3), Be(4), B(5), C(6), N(7), O(8), F(9), Ne(10),
//! Si(14), P(15), S(16), Cl(17), Br(35), I(53)

use std::f64::consts::PI;

/// A single Gaussian primitive: g(r) = coeff · exp(-alpha · r²)
#[derive(Debug, Clone, Copy)]
pub struct GaussianPrimitive {
    /// Gaussian exponent α (Bohr⁻²).
    pub alpha: f64,
    /// Contraction coefficient (pre-normalized).
    pub coefficient: f64,
}

/// A contracted Gaussian shell (a group of primitives sharing a center).
#[derive(Debug, Clone)]
pub struct ContractedShell {
    /// Atom index this shell belongs to.
    pub atom_index: usize,
    /// Center coordinates in Bohr [x, y, z].
    pub center: [f64; 3],
    /// Angular momentum quantum number (0=s, 1=p, 2=d).
    pub l: u32,
    /// Gaussian primitives forming the contraction.
    pub primitives: Vec<GaussianPrimitive>,
}

impl ContractedShell {
    /// Number of Cartesian components for this shell.
    /// s=1, p=3, d=6 (Cartesian)
    pub fn n_cartesian(&self) -> usize {
        match self.l {
            0 => 1,
            1 => 3,
            2 => 6,
            _ => ((self.l + 1) * (self.l + 2) / 2) as usize,
        }
    }

    /// Evaluate the radial part at distance r² from center (s-type only).
    pub fn evaluate_s(&self, r_sq: f64) -> f64 {
        self.primitives
            .iter()
            .map(|p| p.coefficient * (-p.alpha * r_sq).exp())
            .sum()
    }
}

/// A basis function: one component of a contracted shell.
#[derive(Debug, Clone)]
pub struct BasisFunction {
    /// Index of the atom this function is centered on.
    pub atom_index: usize,
    /// Center in Bohr.
    pub center: [f64; 3],
    /// Angular momentum: (l_x, l_y, l_z) such that l_x + l_y + l_z = l.
    pub angular: [u32; 3],
    /// Total angular momentum l = l_x + l_y + l_z.
    pub l_total: u32,
    /// Contracted Gaussian primitives.
    pub primitives: Vec<GaussianPrimitive>,
}

impl BasisFunction {
    /// Normalization constant for a Cartesian Gaussian.
    ///
    /// N = (2α/π)^{3/4} · (4α)^{l/2} · sqrt(1 / ((2lx-1)!! (2ly-1)!! (2lz-1)!!))
    pub fn normalization(alpha: f64, lx: u32, ly: u32, lz: u32) -> f64 {
        let l = lx + ly + lz;
        let prefactor = (2.0 * alpha / PI).powf(0.75) * (4.0 * alpha).powf(l as f64 / 2.0);
        let denom = odd_double_factorial(2 * lx as i32 - 1)
            * odd_double_factorial(2 * ly as i32 - 1)
            * odd_double_factorial(2 * lz as i32 - 1);
        prefactor / denom.sqrt()
    }

    /// Evaluate this basis function at a point (x, y, z) in Bohr.
    pub fn evaluate(&self, x: f64, y: f64, z: f64) -> f64 {
        let dx = x - self.center[0];
        let dy = y - self.center[1];
        let dz = z - self.center[2];
        let r_sq = dx * dx + dy * dy + dz * dz;

        let angular_part = dx.powi(self.angular[0] as i32)
            * dy.powi(self.angular[1] as i32)
            * dz.powi(self.angular[2] as i32);

        let radial: f64 = self
            .primitives
            .iter()
            .map(|p| {
                let norm =
                    Self::normalization(p.alpha, self.angular[0], self.angular[1], self.angular[2]);
                norm * p.coefficient * (-p.alpha * r_sq).exp()
            })
            .sum();

        angular_part * radial
    }
}

/// Double factorial: n!! = n · (n-2) · ... · 1
/// Convention: 0!! = 1, (-1)!! = 1
#[cfg(test)]
fn double_factorial(n: u32) -> u64 {
    if n <= 1 {
        return 1;
    }
    let mut result = 1u64;
    let mut k = n;
    while k > 1 {
        result *= k as u64;
        k -= 2;
    }
    result
}

/// Odd double factorial: (2n-1)!! for normalization.
/// Accepts negative inputs: (-1)!! = 1, (-3)!! = 1.
fn odd_double_factorial(n: i32) -> f64 {
    if n <= 0 {
        return 1.0;
    }
    let mut acc = 1.0;
    let mut k = n;
    while k > 0 {
        acc *= k as f64;
        k -= 2;
    }
    acc
}

/// Complete basis set for a molecular system.
#[derive(Debug, Clone)]
pub struct BasisSet {
    /// All basis functions in canonical order.
    pub functions: Vec<BasisFunction>,
    /// Shell structure (groups of functions sharing exponents).
    pub shells: Vec<ContractedShell>,
    /// Number of basis functions.
    pub n_basis: usize,
    /// Mapping: basis function index → atom index.
    pub function_to_atom: Vec<usize>,
}

impl BasisSet {
    /// Build STO-3G basis set for a molecular system.
    pub fn sto3g(elements: &[u8], positions_bohr: &[[f64; 3]]) -> Self {
        let mut functions = Vec::new();
        let mut shells = Vec::new();
        let mut function_to_atom = Vec::new();

        for (atom_idx, (&z, &center)) in elements.iter().zip(positions_bohr.iter()).enumerate() {
            let atom_shells = get_sto3g_shells(z, atom_idx, center);
            for shell in &atom_shells {
                match shell.l {
                    0 => {
                        functions.push(BasisFunction {
                            atom_index: atom_idx,
                            center,
                            angular: [0, 0, 0],
                            l_total: 0,
                            primitives: shell.primitives.clone(),
                        });
                        function_to_atom.push(atom_idx);
                    }
                    1 => {
                        for (lx, ly, lz) in [(1, 0, 0), (0, 1, 0), (0, 0, 1)] {
                            functions.push(BasisFunction {
                                atom_index: atom_idx,
                                center,
                                angular: [lx, ly, lz],
                                l_total: 1,
                                primitives: shell.primitives.clone(),
                            });
                            function_to_atom.push(atom_idx);
                        }
                    }
                    2 => {
                        for (lx, ly, lz) in [
                            (2, 0, 0),
                            (1, 1, 0),
                            (1, 0, 1),
                            (0, 2, 0),
                            (0, 1, 1),
                            (0, 0, 2),
                        ] {
                            functions.push(BasisFunction {
                                atom_index: atom_idx,
                                center,
                                angular: [lx, ly, lz],
                                l_total: 2,
                                primitives: shell.primitives.clone(),
                            });
                            function_to_atom.push(atom_idx);
                        }
                    }
                    _ => {}
                }
            }
            shells.extend(atom_shells);
        }

        let n_basis = functions.len();
        BasisSet {
            functions,
            shells,
            n_basis,
            function_to_atom,
        }
    }
}

/// Return the STO-3G shells for a given element.
///
/// Parameters from Hehre, Stewart, Pople, J. Chem. Phys. 51, 2657 (1969).
fn get_sto3g_shells(z: u8, atom_index: usize, center: [f64; 3]) -> Vec<ContractedShell> {
    match z {
        // Hydrogen: 1s
        1 => vec![ContractedShell {
            atom_index,
            center,
            l: 0,
            primitives: vec![
                GaussianPrimitive {
                    alpha: 3.42525091,
                    coefficient: 0.15432897,
                },
                GaussianPrimitive {
                    alpha: 0.62353014,
                    coefficient: 0.53532814,
                },
                GaussianPrimitive {
                    alpha: 0.16885540,
                    coefficient: 0.44463454,
                },
            ],
        }],

        // Helium: 1s
        2 => vec![ContractedShell {
            atom_index,
            center,
            l: 0,
            primitives: vec![
                GaussianPrimitive {
                    alpha: 6.36242139,
                    coefficient: 0.15432897,
                },
                GaussianPrimitive {
                    alpha: 1.15892300,
                    coefficient: 0.53532814,
                },
                GaussianPrimitive {
                    alpha: 0.31364979,
                    coefficient: 0.44463454,
                },
            ],
        }],

        // Carbon: 1s, 2s, 2p
        6 => vec![
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 71.6168370,
                        coefficient: 0.15432897,
                    },
                    GaussianPrimitive {
                        alpha: 13.0450960,
                        coefficient: 0.53532814,
                    },
                    GaussianPrimitive {
                        alpha: 3.53051220,
                        coefficient: 0.44463454,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 2.94124940,
                        coefficient: -0.09996723,
                    },
                    GaussianPrimitive {
                        alpha: 0.68348310,
                        coefficient: 0.39951283,
                    },
                    GaussianPrimitive {
                        alpha: 0.22228990,
                        coefficient: 0.70011547,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 2.94124940,
                        coefficient: 0.15591627,
                    },
                    GaussianPrimitive {
                        alpha: 0.68348310,
                        coefficient: 0.60768372,
                    },
                    GaussianPrimitive {
                        alpha: 0.22228990,
                        coefficient: 0.39195739,
                    },
                ],
            },
        ],

        // Nitrogen: 1s, 2s, 2p
        7 => vec![
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 99.1061690,
                        coefficient: 0.15432897,
                    },
                    GaussianPrimitive {
                        alpha: 18.0523120,
                        coefficient: 0.53532814,
                    },
                    GaussianPrimitive {
                        alpha: 4.88566020,
                        coefficient: 0.44463454,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 3.78045590,
                        coefficient: -0.09996723,
                    },
                    GaussianPrimitive {
                        alpha: 0.87849660,
                        coefficient: 0.39951283,
                    },
                    GaussianPrimitive {
                        alpha: 0.28571440,
                        coefficient: 0.70011547,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 3.78045590,
                        coefficient: 0.15591627,
                    },
                    GaussianPrimitive {
                        alpha: 0.87849660,
                        coefficient: 0.60768372,
                    },
                    GaussianPrimitive {
                        alpha: 0.28571440,
                        coefficient: 0.39195739,
                    },
                ],
            },
        ],

        // Oxygen: 1s, 2s, 2p
        8 => vec![
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 130.709320,
                        coefficient: 0.15432897,
                    },
                    GaussianPrimitive {
                        alpha: 23.8088610,
                        coefficient: 0.53532814,
                    },
                    GaussianPrimitive {
                        alpha: 6.44360830,
                        coefficient: 0.44463454,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 5.03315130,
                        coefficient: -0.09996723,
                    },
                    GaussianPrimitive {
                        alpha: 1.16959610,
                        coefficient: 0.39951283,
                    },
                    GaussianPrimitive {
                        alpha: 0.38038900,
                        coefficient: 0.70011547,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 5.03315130,
                        coefficient: 0.15591627,
                    },
                    GaussianPrimitive {
                        alpha: 1.16959610,
                        coefficient: 0.60768372,
                    },
                    GaussianPrimitive {
                        alpha: 0.38038900,
                        coefficient: 0.39195739,
                    },
                ],
            },
        ],

        // Fluorine: 1s, 2s, 2p
        9 => vec![
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 166.679130,
                        coefficient: 0.15432897,
                    },
                    GaussianPrimitive {
                        alpha: 30.3608120,
                        coefficient: 0.53532814,
                    },
                    GaussianPrimitive {
                        alpha: 8.21682070,
                        coefficient: 0.44463454,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 6.46480320,
                        coefficient: -0.09996723,
                    },
                    GaussianPrimitive {
                        alpha: 1.50228120,
                        coefficient: 0.39951283,
                    },
                    GaussianPrimitive {
                        alpha: 0.48858850,
                        coefficient: 0.70011547,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 6.46480320,
                        coefficient: 0.15591627,
                    },
                    GaussianPrimitive {
                        alpha: 1.50228120,
                        coefficient: 0.60768372,
                    },
                    GaussianPrimitive {
                        alpha: 0.48858850,
                        coefficient: 0.39195739,
                    },
                ],
            },
        ],

        // Phosphorus: 1s, 2s, 2p, 3s, 3p
        15 => vec![
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 508.291310,
                        coefficient: 0.15432897,
                    },
                    GaussianPrimitive {
                        alpha: 92.5891370,
                        coefficient: 0.53532814,
                    },
                    GaussianPrimitive {
                        alpha: 25.0571730,
                        coefficient: 0.44463454,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 18.5172490,
                        coefficient: -0.09996723,
                    },
                    GaussianPrimitive {
                        alpha: 4.30422160,
                        coefficient: 0.39951283,
                    },
                    GaussianPrimitive {
                        alpha: 1.39999670,
                        coefficient: 0.70011547,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 18.5172490,
                        coefficient: 0.15591627,
                    },
                    GaussianPrimitive {
                        alpha: 4.30422160,
                        coefficient: 0.60768372,
                    },
                    GaussianPrimitive {
                        alpha: 1.39999670,
                        coefficient: 0.39195739,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 1.56367880,
                        coefficient: -0.09996723,
                    },
                    GaussianPrimitive {
                        alpha: 0.36368650,
                        coefficient: 0.39951283,
                    },
                    GaussianPrimitive {
                        alpha: 0.11828520,
                        coefficient: 0.70011547,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 1.56367880,
                        coefficient: 0.15591627,
                    },
                    GaussianPrimitive {
                        alpha: 0.36368650,
                        coefficient: 0.60768372,
                    },
                    GaussianPrimitive {
                        alpha: 0.11828520,
                        coefficient: 0.39195739,
                    },
                ],
            },
        ],

        // Sulfur: 1s, 2s, 2p, 3s, 3p
        16 => vec![
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 598.642680,
                        coefficient: 0.15432897,
                    },
                    GaussianPrimitive {
                        alpha: 109.046680,
                        coefficient: 0.53532814,
                    },
                    GaussianPrimitive {
                        alpha: 29.5121090,
                        coefficient: 0.44463454,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 22.1564680,
                        coefficient: -0.09996723,
                    },
                    GaussianPrimitive {
                        alpha: 5.14855900,
                        coefficient: 0.39951283,
                    },
                    GaussianPrimitive {
                        alpha: 1.67441430,
                        coefficient: 0.70011547,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 22.1564680,
                        coefficient: 0.15591627,
                    },
                    GaussianPrimitive {
                        alpha: 5.14855900,
                        coefficient: 0.60768372,
                    },
                    GaussianPrimitive {
                        alpha: 1.67441430,
                        coefficient: 0.39195739,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 1.80579080,
                        coefficient: -0.09996723,
                    },
                    GaussianPrimitive {
                        alpha: 0.41988400,
                        coefficient: 0.39951283,
                    },
                    GaussianPrimitive {
                        alpha: 0.13655330,
                        coefficient: 0.70011547,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 1.80579080,
                        coefficient: 0.15591627,
                    },
                    GaussianPrimitive {
                        alpha: 0.41988400,
                        coefficient: 0.60768372,
                    },
                    GaussianPrimitive {
                        alpha: 0.13655330,
                        coefficient: 0.39195739,
                    },
                ],
            },
        ],

        // Chlorine: 1s, 2s, 2p, 3s, 3p
        17 => vec![
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 696.408520,
                        coefficient: 0.15432897,
                    },
                    GaussianPrimitive {
                        alpha: 126.888800,
                        coefficient: 0.53532814,
                    },
                    GaussianPrimitive {
                        alpha: 34.3399080,
                        coefficient: 0.44463454,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 25.9670530,
                        coefficient: -0.09996723,
                    },
                    GaussianPrimitive {
                        alpha: 6.03406090,
                        coefficient: 0.39951283,
                    },
                    GaussianPrimitive {
                        alpha: 1.96235810,
                        coefficient: 0.70011547,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 25.9670530,
                        coefficient: 0.15591627,
                    },
                    GaussianPrimitive {
                        alpha: 6.03406090,
                        coefficient: 0.60768372,
                    },
                    GaussianPrimitive {
                        alpha: 1.96235810,
                        coefficient: 0.39195739,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 2.14407210,
                        coefficient: -0.09996723,
                    },
                    GaussianPrimitive {
                        alpha: 0.49841410,
                        coefficient: 0.39951283,
                    },
                    GaussianPrimitive {
                        alpha: 0.16208590,
                        coefficient: 0.70011547,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 2.14407210,
                        coefficient: 0.15591627,
                    },
                    GaussianPrimitive {
                        alpha: 0.49841410,
                        coefficient: 0.60768372,
                    },
                    GaussianPrimitive {
                        alpha: 0.16208590,
                        coefficient: 0.39195739,
                    },
                ],
            },
        ],

        // Default fallback: use scaled hydrogen-like 1s
        _ => {
            let zeta = 0.5 * (z as f64).sqrt();
            vec![ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 3.42525091 * zeta * zeta,
                        coefficient: 0.15432897,
                    },
                    GaussianPrimitive {
                        alpha: 0.62353014 * zeta * zeta,
                        coefficient: 0.53532814,
                    },
                    GaussianPrimitive {
                        alpha: 0.16885540 * zeta * zeta,
                        coefficient: 0.44463454,
                    },
                ],
            }]
        }
    }
}

/// Number of valence electrons for common elements.
pub fn valence_electrons(z: u8) -> usize {
    match z {
        1 => 1,
        2 => 2,
        3 => 1,
        4 => 2,
        5 => 3,
        6 => 4,
        7 => 5,
        8 => 6,
        9 => 7,
        10 => 8,
        14 => 4,
        15 => 5,
        16 => 6,
        17 => 7,
        35 => 7,
        53 => 7,
        _ => z as usize,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hydrogen_basis_size() {
        let basis = BasisSet::sto3g(&[1], &[[0.0, 0.0, 0.0]]);
        assert_eq!(basis.n_basis, 1);
        assert_eq!(basis.functions[0].l_total, 0);
    }

    #[test]
    fn test_carbon_basis_size() {
        let basis = BasisSet::sto3g(&[6], &[[0.0, 0.0, 0.0]]);
        assert_eq!(basis.n_basis, 5); // 1s + 2s + 2px + 2py + 2pz
    }

    #[test]
    fn test_water_basis_size() {
        let basis = BasisSet::sto3g(
            &[8, 1, 1],
            &[[0.0, 0.0, 0.0], [1.43, 1.11, 0.0], [-1.43, 1.11, 0.0]],
        );
        assert_eq!(basis.n_basis, 7); // O(5) + H(1) + H(1)
    }

    #[test]
    fn test_double_factorial() {
        assert_eq!(double_factorial(0), 1);
        assert_eq!(double_factorial(1), 1);
        assert_eq!(double_factorial(5), 15);
        assert_eq!(double_factorial(6), 48);
    }

    #[test]
    fn test_s_function_evaluation() {
        let basis = BasisSet::sto3g(&[1], &[[0.0, 0.0, 0.0]]);
        let val = basis.functions[0].evaluate(0.0, 0.0, 0.0);
        assert!(val > 0.0);
    }
}
