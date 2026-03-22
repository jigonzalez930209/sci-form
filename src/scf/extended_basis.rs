//! Extended basis set library: 3-21G and 6-31G.
//!
//! Provides larger basis sets beyond the minimal STO-3G for improved
//! accuracy in SCF and post-SCF calculations.

use super::basis::{BasisFunction, BasisSet, GaussianPrimitive};

/// Available basis set choices.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BasisSetType {
    /// Minimal basis: 1 function per AO shell.
    Sto3g,
    /// Split-valence: 2 functions per valence shell.
    Basis321g,
    /// Double-zeta split-valence: better polarization.
    Basis631g,
}

/// Build a basis set of the specified type.
pub fn build_basis_set(
    basis_type: BasisSetType,
    elements: &[u8],
    positions_bohr: &[[f64; 3]],
) -> BasisSet {
    match basis_type {
        BasisSetType::Sto3g => BasisSet::sto3g(elements, positions_bohr),
        BasisSetType::Basis321g => build_321g(elements, positions_bohr),
        BasisSetType::Basis631g => build_631g(elements, positions_bohr),
    }
}

/// Build 3-21G basis set.
///
/// Split-valence basis: inner shell contracted, valence split into
/// 2-GTO and 1-GTO components.
fn build_321g(elements: &[u8], positions_bohr: &[[f64; 3]]) -> BasisSet {
    let mut functions = Vec::new();

    for (idx, (&z, pos)) in elements.iter().zip(positions_bohr.iter()).enumerate() {
        let params = get_321g_params(z);
        for shell in params {
            let lx = shell.angular[0];
            let ly = shell.angular[1];
            let lz = shell.angular[2];
            let primitives: Vec<GaussianPrimitive> = shell
                .exponents
                .iter()
                .zip(shell.coefficients.iter())
                .map(|(&alpha, &coeff)| GaussianPrimitive {
                    alpha,
                    coefficient: coeff,
                })
                .collect();
            functions.push(BasisFunction {
                center: *pos,
                angular: [lx, ly, lz],
                primitives,
                l_total: lx + ly + lz,
                atom_index: idx,
            });
        }
    }

    let function_to_atom: Vec<usize> = functions.iter().map(|f| f.atom_index).collect();
    let n_basis = functions.len();
    BasisSet {
        functions,
        shells: vec![],
        n_basis,
        function_to_atom,
    }
}

/// Build 6-31G basis set.
fn build_631g(elements: &[u8], positions_bohr: &[[f64; 3]]) -> BasisSet {
    let mut functions = Vec::new();

    for (idx, (&z, pos)) in elements.iter().zip(positions_bohr.iter()).enumerate() {
        let params = get_631g_params(z);
        for shell in params {
            let lx = shell.angular[0];
            let ly = shell.angular[1];
            let lz = shell.angular[2];
            let primitives: Vec<GaussianPrimitive> = shell
                .exponents
                .iter()
                .zip(shell.coefficients.iter())
                .map(|(&alpha, &coeff)| GaussianPrimitive {
                    alpha,
                    coefficient: coeff,
                })
                .collect();
            functions.push(BasisFunction {
                center: *pos,
                angular: [lx, ly, lz],
                primitives,
                l_total: lx + ly + lz,
                atom_index: idx,
            });
        }
    }

    let function_to_atom: Vec<usize> = functions.iter().map(|f| f.atom_index).collect();
    let n_basis = functions.len();
    BasisSet {
        functions,
        shells: vec![],
        n_basis,
        function_to_atom,
    }
}

struct ShellParams {
    angular: [u32; 3],
    exponents: Vec<f64>,
    coefficients: Vec<f64>,
}

/// 3-21G parameters for common elements.
fn get_321g_params(z: u8) -> Vec<ShellParams> {
    match z {
        1 => vec![
            // H: 1s (2 GTOs contracted)
            ShellParams {
                angular: [0, 0, 0],
                exponents: vec![5.4471780, 0.8245470],
                coefficients: vec![0.1562850, 0.9046910],
            },
            // H: 1s' (1 GTO)
            ShellParams {
                angular: [0, 0, 0],
                exponents: vec![0.1831920],
                coefficients: vec![1.0000000],
            },
        ],
        6 => {
            let mut shells = vec![
                // C: 1s inner (3 GTOs)
                ShellParams {
                    angular: [0, 0, 0],
                    exponents: vec![172.2560000, 25.9109000, 5.5333500],
                    coefficients: vec![0.0617669, 0.3587940, 0.7007130],
                },
                // C: 2s (2 GTOs)
                ShellParams {
                    angular: [0, 0, 0],
                    exponents: vec![3.6649800, 0.7705450],
                    coefficients: vec![-0.3958970, 1.2158400],
                },
                // C: 2s' (1 GTO)
                ShellParams {
                    angular: [0, 0, 0],
                    exponents: vec![0.1958570],
                    coefficients: vec![1.0000000],
                },
            ];
            // C: 2p (2 GTOs) - px, py, pz
            for angular in &[[1, 0, 0], [0, 1, 0], [0, 0, 1]] {
                shells.push(ShellParams {
                    angular: [angular[0] as u32, angular[1] as u32, angular[2] as u32],
                    exponents: vec![3.6649800, 0.7705450],
                    coefficients: vec![0.2364600, 0.8606190],
                });
            }
            // C: 2p' (1 GTO) - px, py, pz
            for angular in &[[1, 0, 0], [0, 1, 0], [0, 0, 1]] {
                shells.push(ShellParams {
                    angular: [angular[0] as u32, angular[1] as u32, angular[2] as u32],
                    exponents: vec![0.1958570],
                    coefficients: vec![1.0000000],
                });
            }
            shells
        }
        7 => {
            let mut shells = vec![
                ShellParams {
                    angular: [0, 0, 0],
                    exponents: vec![242.7660000, 36.4851000, 7.8144900],
                    coefficients: vec![0.0598657, 0.3529550, 0.7065130],
                },
                ShellParams {
                    angular: [0, 0, 0],
                    exponents: vec![5.4252200, 1.1491500],
                    coefficients: vec![-0.4133010, 1.2244200],
                },
                ShellParams {
                    angular: [0, 0, 0],
                    exponents: vec![0.2832050],
                    coefficients: vec![1.0000000],
                },
            ];
            for angular in &[[1u32, 0, 0], [0, 1, 0], [0, 0, 1]] {
                shells.push(ShellParams {
                    angular: *angular,
                    exponents: vec![5.4252200, 1.1491500],
                    coefficients: vec![0.2379720, 0.8589530],
                });
                shells.push(ShellParams {
                    angular: *angular,
                    exponents: vec![0.2832050],
                    coefficients: vec![1.0000000],
                });
            }
            shells
        }
        8 => {
            let mut shells = vec![
                ShellParams {
                    angular: [0, 0, 0],
                    exponents: vec![322.0370000, 48.4308000, 10.4206000],
                    coefficients: vec![0.0592394, 0.3515000, 0.7076580],
                },
                ShellParams {
                    angular: [0, 0, 0],
                    exponents: vec![7.4029400, 1.5762000],
                    coefficients: vec![-0.4044530, 1.2215600],
                },
                ShellParams {
                    angular: [0, 0, 0],
                    exponents: vec![0.3736840],
                    coefficients: vec![1.0000000],
                },
            ];
            for angular in &[[1u32, 0, 0], [0, 1, 0], [0, 0, 1]] {
                shells.push(ShellParams {
                    angular: *angular,
                    exponents: vec![7.4029400, 1.5762000],
                    coefficients: vec![0.2445860, 0.8539550],
                });
                shells.push(ShellParams {
                    angular: *angular,
                    exponents: vec![0.3736840],
                    coefficients: vec![1.0000000],
                });
            }
            shells
        }
        // For unsupported: fall back to STO-3G-like minimal
        _ => get_321g_fallback(z),
    }
}

fn get_321g_fallback(z: u8) -> Vec<ShellParams> {
    // Minimal s-type for unsupported elements
    let alpha = 0.3 * z as f64;
    vec![ShellParams {
        angular: [0, 0, 0],
        exponents: vec![alpha],
        coefficients: vec![1.0],
    }]
}

/// 6-31G parameters for common elements.
fn get_631g_params(z: u8) -> Vec<ShellParams> {
    match z {
        1 => vec![
            // H: 1s (3 GTOs)
            ShellParams {
                angular: [0, 0, 0],
                exponents: vec![18.7311370, 2.8253937, 0.6401217],
                coefficients: vec![0.03349460, 0.23472695, 0.81375733],
            },
            // H: 1s' (1 GTO)
            ShellParams {
                angular: [0, 0, 0],
                exponents: vec![0.1612778],
                coefficients: vec![1.0000000],
            },
        ],
        6 => {
            let mut shells = vec![
                // C: 1s inner (6 GTOs)
                ShellParams {
                    angular: [0, 0, 0],
                    exponents: vec![
                        3047.5249000,
                        457.3695100,
                        103.9486900,
                        29.2101550,
                        9.2866630,
                        3.1639270,
                    ],
                    coefficients: vec![
                        0.0018347, 0.0140373, 0.0688426, 0.2321844, 0.4679413, 0.3623120,
                    ],
                },
                // C: 2s (3 GTOs)
                ShellParams {
                    angular: [0, 0, 0],
                    exponents: vec![7.8682724, 1.8812885, 0.5442493],
                    coefficients: vec![-0.1193324, -0.1608542, 1.1434564],
                },
                // C: 2s' (1 GTO)
                ShellParams {
                    angular: [0, 0, 0],
                    exponents: vec![0.1687144],
                    coefficients: vec![1.0000000],
                },
            ];
            for angular in &[[1u32, 0, 0], [0, 1, 0], [0, 0, 1]] {
                shells.push(ShellParams {
                    angular: *angular,
                    exponents: vec![7.8682724, 1.8812885, 0.5442493],
                    coefficients: vec![0.0689991, 0.3164240, 0.7443083],
                });
                shells.push(ShellParams {
                    angular: *angular,
                    exponents: vec![0.1687144],
                    coefficients: vec![1.0000000],
                });
            }
            shells
        }
        8 => {
            let mut shells = vec![
                ShellParams {
                    angular: [0, 0, 0],
                    exponents: vec![
                        5484.6717000,
                        825.2349500,
                        188.0469600,
                        52.9645000,
                        16.8975700,
                        5.7996353,
                    ],
                    coefficients: vec![
                        0.0018311, 0.0139501, 0.0684451, 0.2327143, 0.4701929, 0.3585209,
                    ],
                },
                ShellParams {
                    angular: [0, 0, 0],
                    exponents: vec![15.5396160, 3.5999336, 1.0137618],
                    coefficients: vec![-0.1107775, -0.1480263, 1.1307670],
                },
                ShellParams {
                    angular: [0, 0, 0],
                    exponents: vec![0.2700058],
                    coefficients: vec![1.0000000],
                },
            ];
            for angular in &[[1u32, 0, 0], [0, 1, 0], [0, 0, 1]] {
                shells.push(ShellParams {
                    angular: *angular,
                    exponents: vec![15.5396160, 3.5999336, 1.0137618],
                    coefficients: vec![0.0708743, 0.3397528, 0.7271586],
                });
                shells.push(ShellParams {
                    angular: *angular,
                    exponents: vec![0.2700058],
                    coefficients: vec![1.0000000],
                });
            }
            shells
        }
        _ => get_321g_fallback(z),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_321g_h2() {
        let bs = build_basis_set(
            BasisSetType::Basis321g,
            &[1, 1],
            &[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]],
        );
        // H: 2 functions each in 3-21G (contracted + diffuse)
        assert_eq!(bs.n_basis, 4);
    }

    #[test]
    fn test_631g_larger_than_sto3g() {
        let sto3g = BasisSet::sto3g(&[6, 8], &[[0.0, 0.0, 0.0], [2.2, 0.0, 0.0]]);
        let b631g = build_basis_set(
            BasisSetType::Basis631g,
            &[6, 8],
            &[[0.0, 0.0, 0.0], [2.2, 0.0, 0.0]],
        );
        assert!(b631g.n_basis > sto3g.n_basis);
    }
}
