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
        let mut shells = Vec::new();
        for (atom_idx, (&z, &center)) in elements.iter().zip(positions_bohr.iter()).enumerate() {
            shells.extend(get_sto3g_shells(z, atom_idx, center));
        }
        Self::from_shells(shells)
    }

    /// Assemble a basis from contracted shells already assigned to atoms.
    pub fn from_shells(shells: Vec<ContractedShell>) -> Self {
        let mut functions = Vec::new();
        let mut function_to_atom = Vec::new();
        for shell in &shells {
            match shell.l {
                0 => {
                    functions.push(BasisFunction {
                        atom_index: shell.atom_index,
                        center: shell.center,
                        angular: [0, 0, 0],
                        l_total: 0,
                        primitives: shell.primitives.clone(),
                    });
                    function_to_atom.push(shell.atom_index);
                }
                1 => {
                    for (lx, ly, lz) in [(1, 0, 0), (0, 1, 0), (0, 0, 1)] {
                        functions.push(BasisFunction {
                            atom_index: shell.atom_index,
                            center: shell.center,
                            angular: [lx, ly, lz],
                            l_total: 1,
                            primitives: shell.primitives.clone(),
                        });
                        function_to_atom.push(shell.atom_index);
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
                            atom_index: shell.atom_index,
                            center: shell.center,
                            angular: [lx, ly, lz],
                            l_total: 2,
                            primitives: shell.primitives.clone(),
                        });
                        function_to_atom.push(shell.atom_index);
                    }
                }
                _ => {}
            }
        }
        let n_basis = functions.len();
        // Normalize contracted functions (PySCF convention): rescale the
        // coefficients so that the self-overlap of every AO is exactly 1.
        for func in functions.iter_mut() {
            let mut self_overlap = 0.0;
            for prim_a in &func.primitives {
                let norm_a = BasisFunction::normalization(
                    prim_a.alpha,
                    func.angular[0],
                    func.angular[1],
                    func.angular[2],
                );
                for prim_b in &func.primitives {
                    let norm_b = BasisFunction::normalization(
                        prim_b.alpha,
                        func.angular[0],
                        func.angular[1],
                        func.angular[2],
                    );
                    self_overlap += norm_a
                        * prim_a.coefficient
                        * norm_b
                        * prim_b.coefficient
                        * super::gaussian_integrals::overlap_primitive(
                            prim_a.alpha,
                            &func.center,
                            func.angular,
                            prim_b.alpha,
                            &func.center,
                            func.angular,
                        );
                }
            }
            if self_overlap > 0.0 {
                let scale = 1.0 / self_overlap.sqrt();
                for prim in func.primitives.iter_mut() {
                    prim.coefficient *= scale;
                }
            }
        }
        BasisSet {
            functions,
            shells,
            n_basis,
            function_to_atom,
        }
    }
}

/// Hehre–Stewart–Pople STO-3G parameters exist for these atomic numbers.
/// Every other element currently receives a one-function hydrogenic fallback.
pub fn sto3g_tabulated(z: u8) -> bool {
    matches!(z, 1 | 2 | 6 | 7 | 8 | 9 | 15 | 16 | 17)
}

/// Return the STO-3G shells for a given element.
///
/// Parameters from Hehre, Stewart, Pople, J. Chem. Phys. 51, 2657 (1969).
pub(crate) fn get_sto3g_shells(z: u8, atom_index: usize, center: [f64; 3]) -> Vec<ContractedShell> {
    match z {
        // Hydrogen: 1s
        1 => vec![ContractedShell {
            atom_index,
            center,
            l: 0,
            primitives: vec![
                GaussianPrimitive {
                    alpha: 3.425_250_91,
                    coefficient: 0.154_328_97,
                },
                GaussianPrimitive {
                    alpha: 0.623_913_73,
                    coefficient: 0.535_328_14,
                },
                GaussianPrimitive {
                    alpha: 0.168_855_4,
                    coefficient: 0.444_634_54,
                },
            ],
        }],

        2 => vec![ContractedShell {
            atom_index,
            center,
            l: 0,
            primitives: vec![
                GaussianPrimitive {
                    alpha: 6.362_421_39,
                    coefficient: 0.154_328_97,
                },
                GaussianPrimitive {
                    alpha: 1.158_923,
                    coefficient: 0.535_328_14,
                },
                GaussianPrimitive {
                    alpha: 0.313_649_79,
                    coefficient: 0.444_634_54,
                },
            ],
        }],

        6 => vec![
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 71.616_837,
                        coefficient: 0.154_328_97,
                    },
                    GaussianPrimitive {
                        alpha: 13.045_096,
                        coefficient: 0.535_328_14,
                    },
                    GaussianPrimitive {
                        alpha: 3.5305122,
                        coefficient: 0.444_634_54,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 2.941_249_4,
                        coefficient: -0.099_967_23,
                    },
                    GaussianPrimitive {
                        alpha: 0.683_483_1,
                        coefficient: 0.399_512_83,
                    },
                    GaussianPrimitive {
                        alpha: 0.222_289_9,
                        coefficient: 0.700_115_47,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 2.941_249_4,
                        coefficient: 0.15591627,
                    },
                    GaussianPrimitive {
                        alpha: 0.683_483_1,
                        coefficient: 0.607_683_72,
                    },
                    GaussianPrimitive {
                        alpha: 0.222_289_9,
                        coefficient: 0.391_957_39,
                    },
                ],
            },
        ],

        7 => vec![
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 99.106_169,
                        coefficient: 0.154_328_97,
                    },
                    GaussianPrimitive {
                        alpha: 18.052_312,
                        coefficient: 0.535_328_14,
                    },
                    GaussianPrimitive {
                        alpha: 4.885_660_2,
                        coefficient: 0.444_634_54,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 3.780_455_9,
                        coefficient: -0.099_967_23,
                    },
                    GaussianPrimitive {
                        alpha: 0.878_496_6,
                        coefficient: 0.399_512_83,
                    },
                    GaussianPrimitive {
                        alpha: 0.285_714_4,
                        coefficient: 0.700_115_47,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 3.780_455_9,
                        coefficient: 0.15591627,
                    },
                    GaussianPrimitive {
                        alpha: 0.878_496_6,
                        coefficient: 0.607_683_72,
                    },
                    GaussianPrimitive {
                        alpha: 0.285_714_4,
                        coefficient: 0.391_957_39,
                    },
                ],
            },
        ],

        8 => vec![
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 130.709_32,
                        coefficient: 0.154_328_97,
                    },
                    GaussianPrimitive {
                        alpha: 23.808861,
                        coefficient: 0.535_328_14,
                    },
                    GaussianPrimitive {
                        alpha: 6.443_608_3,
                        coefficient: 0.444_634_54,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 5.033_151_3,
                        coefficient: -0.099_967_23,
                    },
                    GaussianPrimitive {
                        alpha: 1.169_596_1,
                        coefficient: 0.399_512_83,
                    },
                    GaussianPrimitive {
                        alpha: 0.380_389,
                        coefficient: 0.700_115_47,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 5.033_151_3,
                        coefficient: 0.15591627,
                    },
                    GaussianPrimitive {
                        alpha: 1.169_596_1,
                        coefficient: 0.607_683_72,
                    },
                    GaussianPrimitive {
                        alpha: 0.380_389,
                        coefficient: 0.391_957_39,
                    },
                ],
            },
        ],

        9 => vec![
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 166.679_13,
                        coefficient: 0.154_328_97,
                    },
                    GaussianPrimitive {
                        alpha: 30.360_812,
                        coefficient: 0.535_328_14,
                    },
                    GaussianPrimitive {
                        alpha: 8.216_820_7,
                        coefficient: 0.444_634_54,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 6.464_803_2,
                        coefficient: -0.099_967_23,
                    },
                    GaussianPrimitive {
                        alpha: 1.502_281_2,
                        coefficient: 0.399_512_83,
                    },
                    GaussianPrimitive {
                        alpha: 0.488_588_5,
                        coefficient: 0.700_115_47,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 6.464_803_2,
                        coefficient: 0.15591627,
                    },
                    GaussianPrimitive {
                        alpha: 1.502_281_2,
                        coefficient: 0.607_683_72,
                    },
                    GaussianPrimitive {
                        alpha: 0.488_588_5,
                        coefficient: 0.391_957_39,
                    },
                ],
            },
        ],

        15 => vec![
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 468.3656378,
                        coefficient: 0.154_328_967_3,
                    },
                    GaussianPrimitive {
                        alpha: 85.313_385_59,
                        coefficient: 0.535_328_142_3,
                    },
                    GaussianPrimitive {
                        alpha: 23.089_131_56,
                        coefficient: 0.444_634_542_2,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 28.032_639_58,
                        coefficient: -0.099_967_229_19,
                    },
                    GaussianPrimitive {
                        alpha: 6.514_182_577,
                        coefficient: 0.3995128261,
                    },
                    GaussianPrimitive {
                        alpha: 2.118_614_352,
                        coefficient: 0.700_115_468_9,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 1.743_103_231,
                        coefficient: -0.219_620_369,
                    },
                    GaussianPrimitive {
                        alpha: 0.486_321_377_1,
                        coefficient: 0.225_595_433_6,
                    },
                    GaussianPrimitive {
                        alpha: 0.190_342_890_9,
                        coefficient: 0.900_398_426,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 28.032_639_58,
                        coefficient: 0.155_916_275,
                    },
                    GaussianPrimitive {
                        alpha: 6.514_182_577,
                        coefficient: 0.607_683_718_6,
                    },
                    GaussianPrimitive {
                        alpha: 2.118_614_352,
                        coefficient: 0.3919573931,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 1.743_103_231,
                        coefficient: 0.01058760429,
                    },
                    GaussianPrimitive {
                        alpha: 0.486_321_377_1,
                        coefficient: 0.5951670053,
                    },
                    GaussianPrimitive {
                        alpha: 0.190_342_890_9,
                        coefficient: 0.462_001_012,
                    },
                ],
            },
        ],

        16 => vec![
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 533.1257359,
                        coefficient: 0.154_328_967_3,
                    },
                    GaussianPrimitive {
                        alpha: 97.109_518_3,
                        coefficient: 0.535_328_142_3,
                    },
                    GaussianPrimitive {
                        alpha: 26.281_625_42,
                        coefficient: 0.444_634_542_2,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 33.329_751_73,
                        coefficient: -0.099_967_229_19,
                    },
                    GaussianPrimitive {
                        alpha: 7.745_117_521,
                        coefficient: 0.3995128261,
                    },
                    GaussianPrimitive {
                        alpha: 2.518_952_599,
                        coefficient: 0.700_115_468_9,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 2.029194274,
                        coefficient: -0.219_620_369,
                    },
                    GaussianPrimitive {
                        alpha: 0.566_140_051_8,
                        coefficient: 0.225_595_433_6,
                    },
                    GaussianPrimitive {
                        alpha: 0.221_583_379_2,
                        coefficient: 0.900_398_426,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 33.329_751_73,
                        coefficient: 0.155_916_275,
                    },
                    GaussianPrimitive {
                        alpha: 7.745_117_521,
                        coefficient: 0.607_683_718_6,
                    },
                    GaussianPrimitive {
                        alpha: 2.518_952_599,
                        coefficient: 0.3919573931,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 2.029194274,
                        coefficient: 0.01058760429,
                    },
                    GaussianPrimitive {
                        alpha: 0.566_140_051_8,
                        coefficient: 0.5951670053,
                    },
                    GaussianPrimitive {
                        alpha: 0.221_583_379_2,
                        coefficient: 0.462_001_012,
                    },
                ],
            },
        ],

        17 => vec![
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 601.345_613_6,
                        coefficient: 0.154_328_967_3,
                    },
                    GaussianPrimitive {
                        alpha: 109.5358542,
                        coefficient: 0.535_328_142_3,
                    },
                    GaussianPrimitive {
                        alpha: 29.644_676_86,
                        coefficient: 0.444_634_542_2,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 38.96041889,
                        coefficient: -0.099_967_229_19,
                    },
                    GaussianPrimitive {
                        alpha: 9.053_563_477,
                        coefficient: 0.3995128261,
                    },
                    GaussianPrimitive {
                        alpha: 2.944_499_834,
                        coefficient: 0.700_115_468_9,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 0,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 2.129_386_495,
                        coefficient: -0.219_620_369,
                    },
                    GaussianPrimitive {
                        alpha: 0.594_093_427_4,
                        coefficient: 0.225_595_433_6,
                    },
                    GaussianPrimitive {
                        alpha: 0.232_524_141,
                        coefficient: 0.900_398_426,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 38.96041889,
                        coefficient: 0.155_916_275,
                    },
                    GaussianPrimitive {
                        alpha: 9.053_563_477,
                        coefficient: 0.607_683_718_6,
                    },
                    GaussianPrimitive {
                        alpha: 2.944_499_834,
                        coefficient: 0.3919573931,
                    },
                ],
            },
            ContractedShell {
                atom_index,
                center,
                l: 1,
                primitives: vec![
                    GaussianPrimitive {
                        alpha: 2.129_386_495,
                        coefficient: 0.01058760429,
                    },
                    GaussianPrimitive {
                        alpha: 0.594_093_427_4,
                        coefficient: 0.5951670053,
                    },
                    GaussianPrimitive {
                        alpha: 0.232_524_141,
                        coefficient: 0.462_001_012,
                    },
                ],
            },
        ],
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
                        alpha: 0.62391373 * zeta * zeta,
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
