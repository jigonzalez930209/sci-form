//! Gaussian basis set definitions for HF-3c.
//!
//! Minimal basis sets for Hartree-Fock: STO-3G contractions for s and p shells.
//! Each primitive Gaussian has the form:
//! $$g(\vec{r}) = N x^l y^m z^n \exp(-\alpha |\vec{r} - \vec{R}|^2)$$

/// Angular momentum type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShellType {
    /// s-type (l=0): 1 function.
    S,
    /// p-type (l=1): 3 functions (px, py, pz).
    P,
}

/// A contracted Gaussian shell centered on an atom.
#[derive(Debug, Clone)]
pub struct Shell {
    /// Center index (atom index).
    pub center_idx: usize,
    /// Center position [x, y, z] in Å (converted to Bohr internally).
    pub center: [f64; 3],
    /// Shell type (S or P).
    pub shell_type: ShellType,
    /// Primitive exponents (α values in Bohr^-2).
    pub exponents: Vec<f64>,
    /// Contraction coefficients (normalized).
    pub coefficients: Vec<f64>,
}

impl Shell {
    /// Number of basis functions in this shell.
    pub fn n_functions(&self) -> usize {
        match self.shell_type {
            ShellType::S => 1,
            ShellType::P => 3,
        }
    }

    /// Normalize contraction coefficients by multiplying by primitive norms
    /// and then normalizing the entire contracted shell.
    pub fn normalize(&mut self) {
        // 1. Multiply by primitive norms
        for (i, &alpha) in self.exponents.iter().enumerate() {
            let n_prim = match self.shell_type {
                ShellType::S => (2.0 * alpha / std::f64::consts::PI).powf(0.75),
                ShellType::P => (128.0 * alpha.powi(5) / std::f64::consts::PI.powi(3)).powf(0.25),
            };
            self.coefficients[i] *= n_prim;
        }

        // 2. Compute self overlap of the unnormalized contracted shell
        let mut sum_s = 0.0;
        for (i, &a) in self.exponents.iter().enumerate() {
            for (j, &b) in self.exponents.iter().enumerate() {
                let gamma = a + b;
                let overlap = match self.shell_type {
                    ShellType::S => (std::f64::consts::PI / gamma).powf(1.5),
                    ShellType::P => {
                        let pre = (std::f64::consts::PI / gamma).powf(1.5);
                        pre * (1.0 / (2.0 * gamma)) // <px | px>
                    }
                };
                sum_s += self.coefficients[i] * self.coefficients[j] * overlap;
            }
        }

        // 3. Normalize the final contracted shell
        let norm_factor = 1.0 / sum_s.sqrt();
        for c in &mut self.coefficients {
            *c *= norm_factor;
        }
    }
}

/// Complete basis set for a molecular system.
#[derive(Debug, Clone)]
pub struct BasisSet {
    /// All shells in the basis.
    pub shells: Vec<Shell>,
}

impl BasisSet {
    /// Total number of basis functions.
    pub fn n_basis(&self) -> usize {
        self.shells.iter().map(|s| s.n_functions()).sum()
    }
}

/// Convert an HF-3c basis to the equivalent `scf` basis set.
///
/// The `scf` integral engine is cross-validated against PySCF
/// (docs/validation/cross-library-validation.es.md, fichas F36/F37), so the
/// HF-3c paths reuse it instead of the legacy local Obara-Saika code.
///
/// HF-3c contraction coefficients already include the primitive
/// normalization factors; the `scf` convention expects bare contraction
/// coefficients (the primitive normalization is applied when evaluating
/// integrals), so those factors are divided out here.
pub(crate) fn to_scf_basis(basis: &BasisSet) -> crate::scf::basis::BasisSet {
    use crate::scf::basis::{
        BasisFunction, BasisSet as ScfBasisSet, ContractedShell, GaussianPrimitive,
    };

    let mut shells = Vec::with_capacity(basis.shells.len());
    for shell in &basis.shells {
        let (l, lx, ly, lz) = match shell.shell_type {
            ShellType::S => (0u32, 0u32, 0u32, 0u32),
            ShellType::P => (1u32, 1u32, 0u32, 0u32),
        };
        let primitives = shell
            .exponents
            .iter()
            .zip(shell.coefficients.iter())
            .map(|(&alpha, &coefficient)| {
                let norm = BasisFunction::normalization(alpha, lx, ly, lz);
                GaussianPrimitive {
                    alpha,
                    coefficient: coefficient / norm,
                }
            })
            .collect();
        shells.push(ContractedShell {
            atom_index: shell.center_idx,
            center: shell.center,
            l,
            primitives,
        });
    }
    ScfBasisSet::from_shells(shells)
}

/// Map each AO index to its parent atom index.
pub fn ao_to_atom_map(basis: &BasisSet) -> Vec<usize> {
    let mut map = Vec::with_capacity(basis.n_basis());
    for shell in &basis.shells {
        for _ in 0..shell.n_functions() {
            map.push(shell.center_idx);
        }
    }
    map
}

/// Ångström to Bohr conversion factor.
pub const ANG_TO_BOHR: f64 = 1.8897259886;

/// Build an STO-3G basis set for the given atoms.
pub fn build_sto3g_basis(elements: &[u8], positions: &[[f64; 3]]) -> BasisSet {
    let mut shells = Vec::new();
    for (idx, (&z, pos)) in elements.iter().zip(positions.iter()).enumerate() {
        let center = [
            pos[0] * ANG_TO_BOHR,
            pos[1] * ANG_TO_BOHR,
            pos[2] * ANG_TO_BOHR,
        ];
        let atom_shells = sto3g_shells(z, idx, center);
        shells.extend(atom_shells);
    }
    BasisSet { shells }
}

/// STO-3G shell definitions per element.
///
/// Exponents and coefficients from Hehre, Stewart, Pople (1969).
fn sto3g_shells(z: u8, center_idx: usize, center: [f64; 3]) -> Vec<Shell> {
    let mut shells = match z {
        1 => vec![Shell {
            center_idx,
            center,
            shell_type: ShellType::S,
            exponents: vec![3.425_250_91, 0.623_913_73, 0.168_855_4],
            coefficients: vec![0.154_328_97, 0.535_328_14, 0.444_634_54],
        }],

        6 => vec![
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![71.616_837, 13.045_096, 3.5305122],
                coefficients: vec![0.154_328_97, 0.535_328_14, 0.444_634_54],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![2.941_249_4, 0.683_483_1, 0.222_289_9],
                coefficients: vec![-0.099_967_23, 0.399_512_83, 0.700_115_47],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![2.941_249_4, 0.683_483_1, 0.222_289_9],
                coefficients: vec![0.15591627, 0.607_683_72, 0.391_957_39],
            },
        ],

        7 => vec![
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![99.106_169, 18.052_312, 4.885_660_2],
                coefficients: vec![0.154_328_97, 0.535_328_14, 0.444_634_54],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![3.780_455_9, 0.878_496_6, 0.285_714_4],
                coefficients: vec![-0.099_967_23, 0.399_512_83, 0.700_115_47],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![3.780_455_9, 0.878_496_6, 0.285_714_4],
                coefficients: vec![0.15591627, 0.607_683_72, 0.391_957_39],
            },
        ],

        8 => vec![
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![130.709_32, 23.808861, 6.443_608_3],
                coefficients: vec![0.154_328_97, 0.535_328_14, 0.444_634_54],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![5.033_151_3, 1.169_596_1, 0.380_389],
                coefficients: vec![-0.099_967_23, 0.399_512_83, 0.700_115_47],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![5.033_151_3, 1.169_596_1, 0.380_389],
                coefficients: vec![0.15591627, 0.607_683_72, 0.391_957_39],
            },
        ],

        9 => vec![
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![166.679_13, 30.360_812, 8.216_820_7],
                coefficients: vec![0.154_328_97, 0.535_328_14, 0.444_634_54],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![6.464_803_2, 1.502_281_2, 0.488_588_5],
                coefficients: vec![-0.099_967_23, 0.399_512_83, 0.700_115_47],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![6.464_803_2, 1.502_281_2, 0.488_588_5],
                coefficients: vec![0.15591627, 0.607_683_72, 0.391_957_39],
            },
        ],

        15 => vec![
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![468.3656378, 85.313_385_59, 23.089_131_56],
                coefficients: vec![0.154_328_967_3, 0.535_328_142_3, 0.444_634_542_2],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![28.032_639_58, 6.514_182_577, 2.118_614_352],
                coefficients: vec![-0.099_967_229_19, 0.3995128261, 0.700_115_468_9],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![1.743_103_231, 0.486_321_377_1, 0.190_342_890_9],
                coefficients: vec![-0.219_620_369, 0.225_595_433_6, 0.900_398_426],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![28.032_639_58, 6.514_182_577, 2.118_614_352],
                coefficients: vec![0.155_916_275, 0.607_683_718_6, 0.3919573931],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![1.743_103_231, 0.486_321_377_1, 0.190_342_890_9],
                coefficients: vec![0.01058760429, 0.5951670053, 0.462_001_012],
            },
        ],

        16 => vec![
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![533.1257359, 97.109_518_3, 26.281_625_42],
                coefficients: vec![0.154_328_967_3, 0.535_328_142_3, 0.444_634_542_2],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![33.329_751_73, 7.745_117_521, 2.518_952_599],
                coefficients: vec![-0.099_967_229_19, 0.3995128261, 0.700_115_468_9],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![2.029194274, 0.566_140_051_8, 0.221_583_379_2],
                coefficients: vec![-0.219_620_369, 0.225_595_433_6, 0.900_398_426],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![33.329_751_73, 7.745_117_521, 2.518_952_599],
                coefficients: vec![0.155_916_275, 0.607_683_718_6, 0.3919573931],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![2.029194274, 0.566_140_051_8, 0.221_583_379_2],
                coefficients: vec![0.01058760429, 0.5951670053, 0.462_001_012],
            },
        ],

        17 => vec![
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![601.345_613_6, 109.5358542, 29.644_676_86],
                coefficients: vec![0.154_328_967_3, 0.535_328_142_3, 0.444_634_542_2],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![38.96041889, 9.053_563_477, 2.944_499_834],
                coefficients: vec![-0.099_967_229_19, 0.3995128261, 0.700_115_468_9],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![2.129_386_495, 0.594_093_427_4, 0.232_524_141],
                coefficients: vec![-0.219_620_369, 0.225_595_433_6, 0.900_398_426],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![38.96041889, 9.053_563_477, 2.944_499_834],
                coefficients: vec![0.155_916_275, 0.607_683_718_6, 0.3919573931],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![2.129_386_495, 0.594_093_427_4, 0.232_524_141],
                coefficients: vec![0.01058760429, 0.5951670053, 0.462_001_012],
            },
        ],

        _ => {
            // Fallback: hydrogen-like 1s with scaled exponent
            vec![Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![
                    3.42525091 * (z as f64 / 1.0).powi(2),
                    0.62391373 * (z as f64 / 1.0).powi(2),
                    0.16885540 * (z as f64 / 1.0).powi(2),
                ],
                coefficients: vec![0.15432897, 0.53532814, 0.44463454],
            }]
        }
    };

    for shell in &mut shells {
        shell.normalize();
    }

    shells
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_h2_basis() {
        let basis = build_sto3g_basis(&[1, 1], &[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]]);
        assert_eq!(basis.n_basis(), 2); // Two 1s orbitals
        assert_eq!(basis.shells.len(), 2);
    }

    #[test]
    fn test_water_basis() {
        let elements = [8u8, 1, 1];
        let positions = [
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469],
        ];
        let basis = build_sto3g_basis(&elements, &positions);
        // O: 1s + 2s + 2p(×3) = 5, 2×H: 1s each = 2 → total = 7
        assert_eq!(basis.n_basis(), 7);
    }
}
