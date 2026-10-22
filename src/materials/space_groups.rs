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
///
/// The seven crystal systems, ordered from lowest to highest symmetry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CrystalSystem {
    /// Triclinic: `a ≠ b ≠ c`, `α ≠ β ≠ γ ≠ 90°`.
    Triclinic,
    /// Monoclinic: two equal lengths, two equal angles, one 90° angle.
    Monoclinic,
    /// Orthorhombic: `a ≠ b ≠ c`, all angles 90°.
    Orthorhombic,
    /// Tetragonal: `a = b ≠ c`, all angles 90°.
    Tetragonal,
    /// Trigonal (rhombohedral): `a = b = c`, `α = β = γ ≠ 90°`.
    Trigonal,
    /// Hexagonal: `a = b ≠ c`, `α = β = 90°`, `γ = 120°`.
    Hexagonal,
    /// Cubic: `a = b = c`, all angles 90°.
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

/// Build a space group with its symmetry operations (all 230, ITC settings).
fn build_space_group(number: u16) -> Option<SpaceGroup> {
    let crystal_system = crystal_system_for_number(number)?;
    let entry = super::space_group_table::group_table(number)?;
    let operations: Vec<SymmetryOp> = entry
        .ops
        .iter()
        .map(|(rotation, translation)| SymmetryOp {
            rotation: *rotation,
            translation: [
                translation[0] as f64 / 12.0,
                translation[1] as f64 / 12.0,
                translation[2] as f64 / 12.0,
            ],
        })
        .collect();
    Some(SpaceGroup {
        number,
        symbol: entry.symbol.to_string(),
        crystal_system,
        lattice_type: entry.lattice,
        point_group: entry.point_group.to_string(),
        operations,
    })
}

/// List every space group as `(ITC number, Hermann-Mauguin symbol)`.
pub fn all_space_groups() -> Vec<(u16, &'static str)> {
    (1..=230u16)
        .filter_map(|n| super::space_group_table::group_table(n).map(|entry| (n, entry.symbol)))
        .collect()
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
