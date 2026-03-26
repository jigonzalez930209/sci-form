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
            exponents: vec![3.42525091, 0.62391373, 0.16885540],
            coefficients: vec![0.15432897, 0.53532814, 0.44463454],
        }],
        6 => vec![
            // 1s core
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![71.6168370, 13.0450960, 3.5305122],
                coefficients: vec![0.15432897, 0.53532814, 0.44463454],
            },
            // 2s valence
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![2.9412494, 0.6834831, 0.2222899],
                coefficients: vec![-0.09996723, 0.39951283, 0.70011547],
            },
            // 2p valence
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![2.9412494, 0.6834831, 0.2222899],
                coefficients: vec![0.15591627, 0.60768372, 0.39195739],
            },
        ],
        7 => vec![
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![99.1061690, 18.0523120, 4.8856602],
                coefficients: vec![0.15432897, 0.53532814, 0.44463454],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![3.7804559, 0.8784966, 0.2857144],
                coefficients: vec![-0.09996723, 0.39951283, 0.70011547],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![3.7804559, 0.8784966, 0.2857144],
                coefficients: vec![0.15591627, 0.60768372, 0.39195739],
            },
        ],
        8 => vec![
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![130.7093200, 23.8088610, 6.4436083],
                coefficients: vec![0.15432897, 0.53532814, 0.44463454],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![5.0331513, 1.1695961, 0.3803890],
                coefficients: vec![-0.09996723, 0.39951283, 0.70011547],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![5.0331513, 1.1695961, 0.3803890],
                coefficients: vec![0.15591627, 0.60768372, 0.39195739],
            },
        ],
        9 => vec![
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![166.6791300, 30.3608120, 8.2168207],
                coefficients: vec![0.15432897, 0.53532814, 0.44463454],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![6.4648032, 1.5022812, 0.4885885],
                coefficients: vec![-0.09996723, 0.39951283, 0.70011547],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![6.4648032, 1.5022812, 0.4885885],
                coefficients: vec![0.15591627, 0.60768372, 0.39195739],
            },
        ],
        15 => vec![
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![508.291310, 92.5891370, 25.0571730],
                coefficients: vec![0.15432897, 0.53532814, 0.44463454],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![18.5172490, 4.30422160, 1.39999670],
                coefficients: vec![-0.09996723, 0.39951283, 0.70011547],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![18.5172490, 4.30422160, 1.39999670],
                coefficients: vec![0.15591627, 0.60768372, 0.39195739],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![1.56367880, 0.36368650, 0.11828520],
                coefficients: vec![-0.09996723, 0.39951283, 0.70011547],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![1.56367880, 0.36368650, 0.11828520],
                coefficients: vec![0.15591627, 0.60768372, 0.39195739],
            },
        ],
        16 => vec![
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![598.642680, 109.046680, 29.5121090],
                coefficients: vec![0.15432897, 0.53532814, 0.44463454],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![22.1564680, 5.14855900, 1.67441430],
                coefficients: vec![-0.09996723, 0.39951283, 0.70011547],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![22.1564680, 5.14855900, 1.67441430],
                coefficients: vec![0.15591627, 0.60768372, 0.39195739],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![1.80579080, 0.41988400, 0.13655330],
                coefficients: vec![-0.09996723, 0.39951283, 0.70011547],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![1.80579080, 0.41988400, 0.13655330],
                coefficients: vec![0.15591627, 0.60768372, 0.39195739],
            },
        ],
        17 => vec![
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![696.408520, 126.888800, 34.3399080],
                coefficients: vec![0.15432897, 0.53532814, 0.44463454],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![25.9670530, 6.03406090, 1.96235810],
                coefficients: vec![-0.09996723, 0.39951283, 0.70011547],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![25.9670530, 6.03406090, 1.96235810],
                coefficients: vec![0.15591627, 0.60768372, 0.39195739],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::S,
                exponents: vec![2.14407210, 0.49841410, 0.16208590],
                coefficients: vec![-0.09996723, 0.39951283, 0.70011547],
            },
            Shell {
                center_idx,
                center,
                shell_type: ShellType::P,
                exponents: vec![2.14407210, 0.49841410, 0.16208590],
                coefficients: vec![0.15591627, 0.60768372, 0.39195739],
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
