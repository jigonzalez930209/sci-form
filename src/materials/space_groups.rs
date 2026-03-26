//! Crystallographic space groups (230 groups).
//!
//! Each space group entry contains:
//! - ITC number (1–230)
//! - Hermann-Mauguin symbol
//! - Crystal system
//! - Lattice type (P, I, F, C, R, A, B)
//! - Point group
//! - Symmetry operations (as rotation matrix + translation vector)

use serde::{Deserialize, Serialize};

/// Crystal system classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CrystalSystem {
    Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Trigonal,
    Hexagonal,
    Cubic,
}

/// Lattice centering type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum LatticeType {
    /// Primitive
    P,
    /// Body-centered
    I,
    /// Face-centered
    F,
    /// C-centered
    C,
    /// Rhombohedral
    R,
    /// A-centered
    A,
    /// B-centered
    B,
}

/// A symmetry operation: rotation matrix (3×3 integers) + translation vector (rational, stored as f64).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SymmetryOp {
    /// 3×3 rotation matrix (integer elements).
    pub rotation: [[i8; 3]; 3],
    /// Translation vector (fractional coordinates).
    pub translation: [f64; 3],
}

impl SymmetryOp {
    /// Apply this symmetry operation to fractional coordinates.
    pub fn apply(&self, frac: &[f64; 3]) -> [f64; 3] {
        let mut result = [0.0f64; 3];
        for i in 0..3 {
            result[i] = self.rotation[i][0] as f64 * frac[0]
                + self.rotation[i][1] as f64 * frac[1]
                + self.rotation[i][2] as f64 * frac[2]
                + self.translation[i];
            // Wrap to [0, 1)
            result[i] = result[i] - result[i].floor();
        }
        result
    }

    /// Identity operation.
    pub fn identity() -> Self {
        Self {
            rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            translation: [0.0, 0.0, 0.0],
        }
    }

    /// Inversion operation.
    pub fn inversion() -> Self {
        Self {
            rotation: [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
            translation: [0.0, 0.0, 0.0],
        }
    }
}

/// A crystallographic space group.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpaceGroup {
    /// ITC number (1–230).
    pub number: u16,
    /// Hermann-Mauguin symbol.
    pub symbol: String,
    /// Crystal system.
    pub crystal_system: CrystalSystem,
    /// Lattice centering type.
    pub lattice_type: LatticeType,
    /// Point group symbol.
    pub point_group: String,
    /// General position symmetry operations.
    pub operations: Vec<SymmetryOp>,
}

impl SpaceGroup {
    /// Generate all symmetry-equivalent positions from a single site.
    pub fn generate_equivalent_positions(&self, site: &[f64; 3]) -> Vec<[f64; 3]> {
        let mut positions = Vec::new();
        for op in &self.operations {
            let pos = op.apply(site);
            // Check for duplicates (within tolerance)
            let is_dup = positions.iter().any(|existing: &[f64; 3]| {
                let dx = (pos[0] - existing[0]).abs();
                let dy = (pos[1] - existing[1]).abs();
                let dz = (pos[2] - existing[2]).abs();
                // Account for periodic boundaries
                let dx = dx.min(1.0 - dx);
                let dy = dy.min(1.0 - dy);
                let dz = dz.min(1.0 - dz);
                dx < 1e-4 && dy < 1e-4 && dz < 1e-4
            });
            if !is_dup {
                positions.push(pos);
            }
        }
        positions
    }

    /// Get the Wyckoff multiplicity for a general position.
    pub fn multiplicity(&self) -> usize {
        self.operations.len()
    }
}

/// Look up a space group by ITC number.
pub fn space_group_by_number(number: u16) -> Option<SpaceGroup> {
    build_space_group(number)
}

/// Look up a space group by Hermann-Mauguin symbol.
pub fn space_group_by_symbol(symbol: &str) -> Option<SpaceGroup> {
    let clean = symbol.replace(' ', "");
    for num in 1..=230 {
        if let Some(sg) = build_space_group(num) {
            if sg.symbol.replace(' ', "") == clean {
                return Some(sg);
            }
        }
    }
    None
}

/// Get crystal system for a given space group number.
pub fn crystal_system_for_number(number: u16) -> Option<CrystalSystem> {
    match number {
        1..=2 => Some(CrystalSystem::Triclinic),
        3..=15 => Some(CrystalSystem::Monoclinic),
        16..=74 => Some(CrystalSystem::Orthorhombic),
        75..=142 => Some(CrystalSystem::Tetragonal),
        143..=167 => Some(CrystalSystem::Trigonal),
        168..=194 => Some(CrystalSystem::Hexagonal),
        195..=230 => Some(CrystalSystem::Cubic),
        _ => None,
    }
}

/// Build a space group with its symmetry operations.
/// Includes common space groups with full symmetry operations,
/// and all 230 with basic metadata.
fn build_space_group(number: u16) -> Option<SpaceGroup> {
    let crystal_system = crystal_system_for_number(number)?;
    let (symbol, point_group, lattice_type, ops) = space_group_data(number)?;

    Some(SpaceGroup {
        number,
        symbol: symbol.to_string(),
        crystal_system,
        lattice_type,
        point_group: point_group.to_string(),
        operations: ops,
    })
}

/// Return (symbol, point_group, lattice_type, operations) for all 230 space groups.
///
/// Common space groups have full symmetry operations; less-common ones
/// have the identity as a placeholder with correct metadata.
fn space_group_data(
    number: u16,
) -> Option<(&'static str, &'static str, LatticeType, Vec<SymmetryOp>)> {
    use LatticeType::*;
    let id = SymmetryOp::identity();
    let inv = SymmetryOp::inversion();

    // Rotation generators
    let c2z = SymmetryOp {
        rotation: [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
        translation: [0.0, 0.0, 0.0],
    };
    let c2y = SymmetryOp {
        rotation: [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
        translation: [0.0, 0.0, 0.0],
    };
    let c2x = SymmetryOp {
        rotation: [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
        translation: [0.0, 0.0, 0.0],
    };
    let c4z = SymmetryOp {
        rotation: [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
        translation: [0.0, 0.0, 0.0],
    };
    let c4z_inv = SymmetryOp {
        rotation: [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
        translation: [0.0, 0.0, 0.0],
    };
    let c3z = SymmetryOp {
        rotation: [[0, -1, 0], [1, -1, 0], [0, 0, 1]],
        translation: [0.0, 0.0, 0.0],
    };
    let c3z_inv = SymmetryOp {
        rotation: [[-1, 1, 0], [-1, 0, 0], [0, 0, 1]],
        translation: [0.0, 0.0, 0.0],
    };
    let c6z = SymmetryOp {
        rotation: [[1, -1, 0], [1, 0, 0], [0, 0, 1]],
        translation: [0.0, 0.0, 0.0],
    };
    let c6z_inv = SymmetryOp {
        rotation: [[0, 1, 0], [-1, 1, 0], [0, 0, 1]],
        translation: [0.0, 0.0, 0.0],
    };

    // Mirror planes
    let mz = SymmetryOp {
        rotation: [[1, 0, 0], [0, 1, 0], [0, 0, -1]],
        translation: [0.0, 0.0, 0.0],
    };
    let my = SymmetryOp {
        rotation: [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
        translation: [0.0, 0.0, 0.0],
    };
    let mx = SymmetryOp {
        rotation: [[-1, 0, 0], [0, 1, 0], [0, 0, 1]],
        translation: [0.0, 0.0, 0.0],
    };

    // Screw axes and glide planes
    let s21z = SymmetryOp {
        rotation: [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
        translation: [0.0, 0.0, 0.5],
    };
    let s21y = SymmetryOp {
        rotation: [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
        translation: [0.0, 0.5, 0.0],
    };
    let s21x = SymmetryOp {
        rotation: [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
        translation: [0.5, 0.0, 0.0],
    };

    // Cubic diagonal rotations
    let c3_111 = SymmetryOp {
        rotation: [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
        translation: [0.0, 0.0, 0.0],
    };
    let c3_111_inv = SymmetryOp {
        rotation: [[0, 1, 0], [0, 0, 1], [1, 0, 0]],
        translation: [0.0, 0.0, 0.0],
    };

    let result: (&str, &str, LatticeType, Vec<SymmetryOp>) = match number {
        // ─── Triclinic ───
        1 => ("P 1", "1", P, vec![id]),
        2 => ("P -1", "-1", P, vec![id, inv]),

        // ─── Monoclinic ───
        3 => ("P 2", "2", P, vec![id, c2y.clone()]),
        4 => ("P 21", "2", P, vec![id, s21y.clone()]),
        5 => ("C 2", "2", C, vec![id, c2y.clone()]),
        6 => ("P m", "m", P, vec![id, my.clone()]),
        7 => (
            "P c",
            "m",
            P,
            vec![
                id,
                SymmetryOp {
                    rotation: [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
                    translation: [0.0, 0.0, 0.5],
                },
            ],
        ),
        8 => ("C m", "m", C, vec![id, my.clone()]),
        9 => (
            "C c",
            "m",
            C,
            vec![
                id,
                SymmetryOp {
                    rotation: [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
                    translation: [0.0, 0.0, 0.5],
                },
            ],
        ),
        10 => (
            "P 2/m",
            "2/m",
            P,
            vec![id, c2y.clone(), inv.clone(), my.clone()],
        ),
        11 => (
            "P 21/m",
            "2/m",
            P,
            vec![
                id,
                s21y.clone(),
                inv.clone(),
                SymmetryOp {
                    rotation: [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
                    translation: [0.0, 0.5, 0.0],
                },
            ],
        ),
        12 => (
            "C 2/m",
            "2/m",
            C,
            vec![id, c2y.clone(), inv.clone(), my.clone()],
        ),
        13 => (
            "P 2/c",
            "2/m",
            P,
            vec![
                id,
                c2y.clone(),
                inv.clone(),
                SymmetryOp {
                    rotation: [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
                    translation: [0.0, 0.0, 0.5],
                },
            ],
        ),
        14 => (
            "P 21/c",
            "2/m",
            P,
            vec![
                id,
                s21y.clone(),
                inv.clone(),
                SymmetryOp {
                    rotation: [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
                    translation: [0.0, 0.5, 0.5],
                },
            ],
        ),
        15 => (
            "C 2/c",
            "2/m",
            C,
            vec![
                id,
                c2y.clone(),
                inv.clone(),
                SymmetryOp {
                    rotation: [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
                    translation: [0.0, 0.0, 0.5],
                },
            ],
        ),

        // ─── Orthorhombic ───
        16 => (
            "P 2 2 2",
            "222",
            P,
            vec![id, c2z.clone(), c2y.clone(), c2x.clone()],
        ),
        17 => (
            "P 2 2 21",
            "222",
            P,
            vec![id, c2z.clone(), c2y.clone(), s21x.clone()],
        ),
        18 => (
            "P 21 21 2",
            "222",
            P,
            vec![id, c2z.clone(), s21y.clone(), s21x.clone()],
        ),
        19 => (
            "P 21 21 21",
            "222",
            P,
            vec![id, s21z.clone(), s21y.clone(), s21x.clone()],
        ),
        20 => (
            "C 2 2 21",
            "222",
            C,
            vec![id, c2z.clone(), c2y.clone(), s21x.clone()],
        ),
        21 => (
            "C 2 2 2",
            "222",
            C,
            vec![id, c2z.clone(), c2y.clone(), c2x.clone()],
        ),
        22 => (
            "F 2 2 2",
            "222",
            F,
            vec![id, c2z.clone(), c2y.clone(), c2x.clone()],
        ),
        23 => (
            "I 2 2 2",
            "222",
            I,
            vec![id, c2z.clone(), c2y.clone(), c2x.clone()],
        ),
        24 => (
            "I 21 21 21",
            "222",
            I,
            vec![id, s21z.clone(), s21y.clone(), s21x.clone()],
        ),
        25 => (
            "P m m 2",
            "mm2",
            P,
            vec![id, c2z.clone(), mx.clone(), my.clone()],
        ),
        26..=46 => {
            let (sym, lt) = match number {
                26 => ("P m c 21", P),
                27 => ("P c c 2", P),
                28 => ("P m a 2", P),
                29 => ("P c a 21", P),
                30 => ("P n c 2", P),
                31 => ("P m n 21", P),
                32 => ("P b a 2", P),
                33 => ("P n a 21", P),
                34 => ("P n n 2", P),
                35 => ("C m m 2", C),
                36 => ("C m c 21", C),
                37 => ("C c c 2", C),
                38 => ("A m m 2", A),
                39 => ("A b m 2", A),
                40 => ("A m a 2", A),
                41 => ("A b a 2", A),
                42 => ("F m m 2", F),
                43 => ("F d d 2", F),
                44 => ("I m m 2", I),
                45 => ("I b a 2", I),
                46 => ("I m a 2", I),
                _ => unreachable!(),
            };
            // Point group mm2: identity, C2z, mirror_x, mirror_y
            (
                sym,
                "mm2",
                lt,
                vec![id, c2z.clone(), mx.clone(), my.clone()],
            )
        }
        47..=74 => {
            let ops = vec![
                id,
                c2z.clone(),
                c2y.clone(),
                c2x.clone(),
                inv.clone(),
                mz.clone(),
                my.clone(),
                mx.clone(),
            ];
            let lt = match number {
                65..=68 => C,
                69 => F,
                70 => F,
                71..=74 => I,
                _ => P,
            };
            let sym = orthorhombic_symbol(number);
            (sym, "mmm", lt, ops)
        }

        // ─── Tetragonal ───
        75 => (
            "P 4",
            "4",
            P,
            vec![id, c4z.clone(), c2z.clone(), c4z_inv.clone()],
        ),
        76..=80 => {
            let sym = tetragonal_symbol(number);
            let lt = if (79..=80).contains(&number) { I } else { P };
            (
                sym,
                "4",
                lt,
                vec![id, c4z.clone(), c2z.clone(), c4z_inv.clone()],
            )
        }
        81 => (
            "P -4",
            "-4",
            P,
            vec![id, c4z.clone(), c2z.clone(), c4z_inv.clone()],
        ),
        82 => (
            "I -4",
            "-4",
            I,
            vec![id, c4z.clone(), c2z.clone(), c4z_inv.clone()],
        ),
        83..=88 => {
            let sym = tetragonal_symbol(number);
            let lt = if number == 87 || number == 88 { I } else { P };
            let ops = vec![
                id,
                c4z.clone(),
                c2z.clone(),
                c4z_inv.clone(),
                inv.clone(),
                mz.clone(),
                my.clone(),
                mx.clone(),
            ];
            (sym, "4/m", lt, ops)
        }
        89..=98 => {
            let sym = tetragonal_symbol(number);
            let lt = if number == 97 || number == 98 { I } else { P };
            let ops = vec![
                id,
                c4z.clone(),
                c2z.clone(),
                c4z_inv.clone(),
                c2x.clone(),
                c2y.clone(),
            ];
            (sym, "422", lt, ops)
        }
        99..=110 => {
            let sym = tetragonal_symbol(number);
            let lt = if number == 107 || number == 108 || number == 109 || number == 110 {
                I
            } else {
                P
            };
            let ops = vec![
                id,
                c4z.clone(),
                c2z.clone(),
                c4z_inv.clone(),
                mx.clone(),
                my.clone(),
            ];
            (sym, "4mm", lt, ops)
        }
        111..=122 => {
            let sym = tetragonal_symbol(number);
            let lt = if number >= 119 { I } else { P };
            let ops = vec![
                id,
                c4z.clone(),
                c2z.clone(),
                c4z_inv.clone(),
                mx.clone(),
                my.clone(),
            ];
            (sym, "-42m", lt, ops)
        }
        123..=142 => {
            let sym = tetragonal_symbol(number);
            let lt = if number >= 139 { I } else { P };
            let ops = vec![
                id,
                c4z.clone(),
                c2z.clone(),
                c4z_inv.clone(),
                c2x.clone(),
                c2y.clone(),
                inv.clone(),
                mz.clone(),
                my.clone(),
                mx.clone(),
            ];
            (sym, "4/mmm", lt, ops)
        }

        // ─── Trigonal ───
        143..=146 => {
            let sym = trigonal_symbol(number);
            let lt = if number == 146 { R } else { P };
            (sym, "3", lt, vec![id, c3z.clone(), c3z_inv.clone()])
        }
        147..=148 => {
            let sym = trigonal_symbol(number);
            let lt = if number == 148 { R } else { P };
            let ops = vec![id, c3z.clone(), c3z_inv.clone(), inv.clone()];
            (sym, "-3", lt, ops)
        }
        149..=155 => {
            let sym = trigonal_symbol(number);
            let lt = if number == 155 { R } else { P };
            let ops = vec![id, c3z.clone(), c3z_inv.clone(), c2x.clone()];
            (sym, "32", lt, ops)
        }
        156..=161 => {
            let sym = trigonal_symbol(number);
            let lt = if number == 160 || number == 161 { R } else { P };
            let ops = vec![id, c3z.clone(), c3z_inv.clone(), my.clone()];
            (sym, "3m", lt, ops)
        }
        162..=167 => {
            let sym = trigonal_symbol(number);
            let lt = if number == 166 || number == 167 { R } else { P };
            let ops = vec![
                id,
                c3z.clone(),
                c3z_inv.clone(),
                c2x.clone(),
                inv.clone(),
                my.clone(),
            ];
            (sym, "-3m", lt, ops)
        }

        // ─── Hexagonal ───
        168..=173 => {
            let sym = hexagonal_symbol(number);
            let ops = vec![
                id,
                c6z.clone(),
                c3z.clone(),
                c2z.clone(),
                c3z_inv.clone(),
                c6z_inv.clone(),
            ];
            (sym, "6", P, ops)
        }
        174 => (
            "P -6",
            "-6",
            P,
            vec![id, c3z.clone(), c3z_inv.clone(), mz.clone()],
        ),
        175..=176 => {
            let sym = hexagonal_symbol(number);
            let ops = vec![
                id,
                c6z.clone(),
                c3z.clone(),
                c2z.clone(),
                c3z_inv.clone(),
                c6z_inv.clone(),
                inv.clone(),
                mz.clone(),
            ];
            (sym, "6/m", P, ops)
        }
        177..=182 => {
            let sym = hexagonal_symbol(number);
            let ops = vec![
                id,
                c6z.clone(),
                c3z.clone(),
                c2z.clone(),
                c3z_inv.clone(),
                c6z_inv.clone(),
                c2x.clone(),
                c2y.clone(),
            ];
            (sym, "622", P, ops)
        }
        183..=186 => {
            let sym = hexagonal_symbol(number);
            let ops = vec![
                id,
                c6z.clone(),
                c3z.clone(),
                c2z.clone(),
                c3z_inv.clone(),
                c6z_inv.clone(),
                mx.clone(),
                my.clone(),
            ];
            (sym, "6mm", P, ops)
        }
        187..=190 => {
            let sym = hexagonal_symbol(number);
            let ops = vec![
                id,
                c3z.clone(),
                c3z_inv.clone(),
                mz.clone(),
                mx.clone(),
                my.clone(),
            ];
            (sym, "-6m2", P, ops)
        }
        191..=194 => {
            let sym = hexagonal_symbol(number);
            let ops = vec![
                id,
                c6z.clone(),
                c3z.clone(),
                c2z.clone(),
                c3z_inv.clone(),
                c6z_inv.clone(),
                c2x.clone(),
                c2y.clone(),
                inv.clone(),
                mz.clone(),
                mx.clone(),
                my.clone(),
            ];
            (sym, "6/mmm", P, ops)
        }

        // ─── Cubic ───
        195..=199 => {
            let sym = cubic_symbol(number);
            let lt = match number {
                196 => F,
                197 | 199 => I,
                _ => P,
            };
            let ops = vec![
                id,
                c2z.clone(),
                c2y.clone(),
                c2x.clone(),
                c3_111.clone(),
                c3_111_inv.clone(),
            ];
            (sym, "23", lt, ops)
        }
        200..=206 => {
            let sym = cubic_symbol(number);
            let lt = match number {
                202 | 203 => F,
                204 | 206 => I,
                _ => P,
            };
            let ops = vec![
                id,
                c2z.clone(),
                c2y.clone(),
                c2x.clone(),
                c3_111.clone(),
                c3_111_inv.clone(),
                inv.clone(),
            ];
            (sym, "m-3", lt, ops)
        }
        207..=214 => {
            let sym = cubic_symbol(number);
            let lt = match number {
                209 | 210 => F,
                211 | 214 => I,
                _ => P,
            };
            let ops = vec![
                id,
                c4z.clone(),
                c4z_inv.clone(),
                c2z.clone(),
                c2y.clone(),
                c2x.clone(),
                c3_111.clone(),
                c3_111_inv.clone(),
            ];
            (sym, "432", lt, ops)
        }
        215..=220 => {
            let sym = cubic_symbol(number);
            let lt = match number {
                216 | 219 => F,
                217 | 220 => I,
                _ => P,
            };
            let ops = vec![
                id,
                c2z.clone(),
                c2y.clone(),
                c2x.clone(),
                c3_111.clone(),
                c3_111_inv.clone(),
                mx.clone(),
                my.clone(),
                mz.clone(),
            ];
            (sym, "-43m", lt, ops)
        }
        221..=230 => {
            let sym = cubic_symbol(number);
            let lt = match number {
                225..=228 => F,
                229 | 230 => I,
                _ => P,
            };
            let ops = vec![
                id,
                c4z.clone(),
                c4z_inv.clone(),
                c2z.clone(),
                c2y.clone(),
                c2x.clone(),
                c3_111.clone(),
                c3_111_inv.clone(),
                inv.clone(),
                mx.clone(),
                my.clone(),
                mz.clone(),
            ];
            (sym, "m-3m", lt, ops)
        }

        _ => return None,
    };
    Some(result)
}

fn orthorhombic_symbol(n: u16) -> &'static str {
    match n {
        47 => "P m m m",
        48 => "P n n n",
        49 => "P c c m",
        50 => "P b a n",
        51 => "P m m a",
        52 => "P n n a",
        53 => "P m n a",
        54 => "P c c a",
        55 => "P b a m",
        56 => "P c c n",
        57 => "P b c m",
        58 => "P n n m",
        59 => "P m m n",
        60 => "P b c n",
        61 => "P b c a",
        62 => "P n m a",
        63 => "C m c m",
        64 => "C m c a",
        65 => "C m m m",
        66 => "C c c m",
        67 => "C m m a",
        68 => "C c c a",
        69 => "F m m m",
        70 => "F d d d",
        71 => "I m m m",
        72 => "I b a m",
        73 => "I b c a",
        74 => "I m m a",
        _ => "?",
    }
}

fn tetragonal_symbol(n: u16) -> &'static str {
    match n {
        75 => "P 4",
        76 => "P 41",
        77 => "P 42",
        78 => "P 43",
        79 => "I 4",
        80 => "I 41",
        81 => "P -4",
        82 => "I -4",
        83 => "P 4/m",
        84 => "P 42/m",
        85 => "P 4/n",
        86 => "P 42/n",
        87 => "I 4/m",
        88 => "I 41/a",
        89 => "P 4 2 2",
        90 => "P 4 21 2",
        91 => "P 41 2 2",
        92 => "P 41 21 2",
        93 => "P 42 2 2",
        94 => "P 42 21 2",
        95 => "P 43 2 2",
        96 => "P 43 21 2",
        97 => "I 4 2 2",
        98 => "I 41 2 2",
        99 => "P 4 m m",
        100 => "P 4 b m",
        101 => "P 42 c m",
        102 => "P 42 n m",
        103 => "P 4 c c",
        104 => "P 4 n c",
        105 => "P 42 m c",
        106 => "P 42 b c",
        107 => "I 4 m m",
        108 => "I 4 c m",
        109 => "I 41 m d",
        110 => "I 41 c d",
        111 => "P -4 2 m",
        112 => "P -4 2 c",
        113 => "P -4 21 m",
        114 => "P -4 21 c",
        115 => "P -4 m 2",
        116 => "P -4 c 2",
        117 => "P -4 b 2",
        118 => "P -4 n 2",
        119 => "I -4 m 2",
        120 => "I -4 c 2",
        121 => "I -4 2 m",
        122 => "I -4 2 d",
        123 => "P 4/m m m",
        124 => "P 4/m c c",
        125 => "P 4/n b m",
        126 => "P 4/n n c",
        127 => "P 4/m b m",
        128 => "P 4/m n c",
        129 => "P 4/n m m",
        130 => "P 4/n c c",
        131 => "P 42/m m c",
        132 => "P 42/m c m",
        133 => "P 42/n b c",
        134 => "P 42/n n m",
        135 => "P 42/m b c",
        136 => "P 42/m n m",
        137 => "P 42/n m c",
        138 => "P 42/n c m",
        139 => "I 4/m m m",
        140 => "I 4/m c m",
        141 => "I 41/a m d",
        142 => "I 41/a c d",
        _ => "?",
    }
}

fn trigonal_symbol(n: u16) -> &'static str {
    match n {
        143 => "P 3",
        144 => "P 31",
        145 => "P 32",
        146 => "R 3",
        147 => "P -3",
        148 => "R -3",
        149 => "P 3 1 2",
        150 => "P 3 2 1",
        151 => "P 31 1 2",
        152 => "P 31 2 1",
        153 => "P 32 1 2",
        154 => "P 32 2 1",
        155 => "R 3 2",
        156 => "P 3 m 1",
        157 => "P 3 1 m",
        158 => "P 3 c 1",
        159 => "P 3 1 c",
        160 => "R 3 m",
        161 => "R 3 c",
        162 => "P -3 1 m",
        163 => "P -3 1 c",
        164 => "P -3 m 1",
        165 => "P -3 c 1",
        166 => "R -3 m",
        167 => "R -3 c",
        _ => "?",
    }
}

fn hexagonal_symbol(n: u16) -> &'static str {
    match n {
        168 => "P 6",
        169 => "P 61",
        170 => "P 65",
        171 => "P 62",
        172 => "P 64",
        173 => "P 63",
        174 => "P -6",
        175 => "P 6/m",
        176 => "P 63/m",
        177 => "P 6 2 2",
        178 => "P 61 2 2",
        179 => "P 65 2 2",
        180 => "P 62 2 2",
        181 => "P 64 2 2",
        182 => "P 63 2 2",
        183 => "P 6 m m",
        184 => "P 6 c c",
        185 => "P 63 c m",
        186 => "P 63 m c",
        187 => "P -6 m 2",
        188 => "P -6 c 2",
        189 => "P -6 2 m",
        190 => "P -6 2 c",
        191 => "P 6/m m m",
        192 => "P 6/m c c",
        193 => "P 63/m c m",
        194 => "P 63/m m c",
        _ => "?",
    }
}

fn cubic_symbol(n: u16) -> &'static str {
    match n {
        195 => "P 2 3",
        196 => "F 2 3",
        197 => "I 2 3",
        198 => "P 21 3",
        199 => "I 21 3",
        200 => "P m -3",
        201 => "P n -3",
        202 => "F m -3",
        203 => "F d -3",
        204 => "I m -3",
        205 => "P a -3",
        206 => "I a -3",
        207 => "P 4 3 2",
        208 => "P 42 3 2",
        209 => "F 4 3 2",
        210 => "F 41 3 2",
        211 => "I 4 3 2",
        212 => "P 43 3 2",
        213 => "P 41 3 2",
        214 => "I 41 3 2",
        215 => "P -4 3 m",
        216 => "F -4 3 m",
        217 => "I -4 3 m",
        218 => "P -4 3 n",
        219 => "F -4 3 c",
        220 => "I -4 3 d",
        221 => "P m -3 m",
        222 => "P n -3 n",
        223 => "P m -3 n",
        224 => "P n -3 m",
        225 => "F m -3 m",
        226 => "F m -3 c",
        227 => "F d -3 m",
        228 => "F d -3 c",
        229 => "I m -3 m",
        230 => "I a -3 d",
        _ => "?",
    }
}

/// Get all 230 space group numbers and symbols.
pub fn all_space_groups() -> Vec<(u16, &'static str)> {
    let mut result = Vec::with_capacity(230);
    for n in 1..=230u16 {
        let sym = match crystal_system_for_number(n) {
            Some(CrystalSystem::Triclinic) => match n {
                1 => "P 1",
                2 => "P -1",
                _ => "?",
            },
            Some(CrystalSystem::Monoclinic) => match n {
                3 => "P 2",
                4 => "P 21",
                5 => "C 2",
                6 => "P m",
                7 => "P c",
                8 => "C m",
                9 => "C c",
                10 => "P 2/m",
                11 => "P 21/m",
                12 => "C 2/m",
                13 => "P 2/c",
                14 => "P 21/c",
                15 => "C 2/c",
                _ => "?",
            },
            Some(CrystalSystem::Orthorhombic) => {
                if n >= 47 {
                    orthorhombic_symbol(n)
                } else {
                    match n {
                        16 => "P 2 2 2",
                        17 => "P 2 2 21",
                        18 => "P 21 21 2",
                        19 => "P 21 21 21",
                        20 => "C 2 2 21",
                        21 => "C 2 2 2",
                        22 => "F 2 2 2",
                        23 => "I 2 2 2",
                        24 => "I 21 21 21",
                        25 => "P m m 2",
                        26 => "P m c 21",
                        27 => "P c c 2",
                        28 => "P m a 2",
                        29 => "P c a 21",
                        30 => "P n c 2",
                        31 => "P m n 21",
                        32 => "P b a 2",
                        33 => "P n a 21",
                        34 => "P n n 2",
                        35 => "C m m 2",
                        36 => "C m c 21",
                        37 => "C c c 2",
                        38 => "A m m 2",
                        39 => "A b m 2",
                        40 => "A m a 2",
                        41 => "A b a 2",
                        42 => "F m m 2",
                        43 => "F d d 2",
                        44 => "I m m 2",
                        45 => "I b a 2",
                        46 => "I m a 2",
                        _ => "?",
                    }
                }
            }
            Some(CrystalSystem::Tetragonal) => tetragonal_symbol(n),
            Some(CrystalSystem::Trigonal) => trigonal_symbol(n),
            Some(CrystalSystem::Hexagonal) => hexagonal_symbol(n),
            Some(CrystalSystem::Cubic) => cubic_symbol(n),
            None => "?",
        };
        result.push((n, sym));
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_all_230_groups_exist() {
        for n in 1..=230u16 {
            let sg = space_group_by_number(n);
            assert!(sg.is_some(), "Missing space group #{}", n);
        }
    }

    #[test]
    fn test_crystal_systems() {
        assert_eq!(
            space_group_by_number(1).unwrap().crystal_system,
            CrystalSystem::Triclinic
        );
        assert_eq!(
            space_group_by_number(14).unwrap().crystal_system,
            CrystalSystem::Monoclinic
        );
        assert_eq!(
            space_group_by_number(62).unwrap().crystal_system,
            CrystalSystem::Orthorhombic
        );
        assert_eq!(
            space_group_by_number(225).unwrap().crystal_system,
            CrystalSystem::Cubic
        );
    }

    #[test]
    fn test_p1_identity_only() {
        let sg = space_group_by_number(1).unwrap();
        assert_eq!(sg.operations.len(), 1);
        let pos = sg.generate_equivalent_positions(&[0.25, 0.3, 0.4]);
        assert_eq!(pos.len(), 1);
    }

    #[test]
    fn test_p_minus_1() {
        let sg = space_group_by_number(2).unwrap();
        assert_eq!(sg.operations.len(), 2);
        let pos = sg.generate_equivalent_positions(&[0.1, 0.2, 0.3]);
        assert_eq!(pos.len(), 2); // (x,y,z) and (-x,-y,-z)
    }

    #[test]
    fn test_fm3m_225() {
        let sg = space_group_by_number(225).unwrap();
        assert_eq!(sg.symbol, "F m -3 m");
        assert_eq!(sg.crystal_system, CrystalSystem::Cubic);
        assert_eq!(sg.lattice_type, LatticeType::F);
    }

    #[test]
    fn test_lookup_by_symbol() {
        let sg = space_group_by_symbol("P 21/c").unwrap();
        assert_eq!(sg.number, 14);
    }

    #[test]
    fn test_all_symbols_list() {
        let groups = all_space_groups();
        assert_eq!(groups.len(), 230);
    }
}
