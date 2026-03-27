//! Molecular integration grid for DFT.
//!
//! Combines Euler-Maclaurin radial quadrature with Lebedev angular quadrature
//! and Becke partitioning for multi-center numerical integration.

use super::becke::becke_weights;
use super::lebedev::lebedev_grid;

/// A single quadrature point with position, weight, and atom assignment.
#[derive(Debug, Clone)]
pub struct GridPoint {
    /// Cartesian position (Bohr).
    pub xyz: [f64; 3],
    /// Combined quadrature weight.
    pub weight: f64,
    /// Index of the atom this point belongs to.
    pub atom_index: usize,
}

/// Molecular integration grid.
#[derive(Debug, Clone)]
pub struct MolecularGrid {
    /// All quadrature points.
    pub points: Vec<GridPoint>,
    /// Total number of points.
    pub n_points: usize,
}

/// Grid quality settings.
#[derive(Debug, Clone, Copy)]
pub enum GridQuality {
    /// Coarse: ~50 radial × 6 angular.
    Coarse,
    /// Medium: ~75 radial × 26 angular.
    Medium,
    /// Fine: ~100 radial × 110 angular.
    Fine,
    /// Very fine: ~150 radial × 302 angular.
    VeryFine,
}

impl GridQuality {
    fn radial_points(self) -> usize {
        match self {
            Self::Coarse => 50,
            Self::Medium => 75,
            Self::Fine => 100,
            Self::VeryFine => 150,
        }
    }

    fn angular_points(self) -> usize {
        match self {
            Self::Coarse => 6,
            Self::Medium => 26,
            Self::Fine => 110,
            Self::VeryFine => 302,
        }
    }
}

impl MolecularGrid {
    /// Build a molecular integration grid.
    ///
    /// Uses Treutler-Ahlrichs radial mapping, Lebedev angular quadrature,
    /// and Becke partitioning across atomic centers.
    pub fn build(atomic_numbers: &[u8], positions_bohr: &[[f64; 3]], quality: GridQuality) -> Self {
        let n_atoms = atomic_numbers.len();
        let n_rad = quality.radial_points();
        let n_ang = quality.angular_points();
        let angular = lebedev_grid(n_ang);

        let mut points = Vec::with_capacity(n_atoms * n_rad * angular.len());

        for (atom_idx, &center) in positions_bohr.iter().enumerate() {
            let r_bragg = bragg_slater_radius(atomic_numbers[atom_idx]);

            for i_rad in 0..n_rad {
                // Treutler-Ahlrichs radial mapping: r = r_bragg * (1+x)/(1-x)
                // with Gauss-Chebyshev abscissa x_i
                let x = ((i_rad as f64 + 0.5) / n_rad as f64) * std::f64::consts::PI;
                let cos_x = x.cos();
                let r = r_bragg * (1.0 + cos_x) / (1.0 - cos_x + 1e-15);
                let r = r.max(1e-10);

                // Radial weight: Euler-Maclaurin with Jacobian
                let w_rad = std::f64::consts::PI / n_rad as f64 * x.sin() * r_bragg * 2.0
                    / ((1.0 - cos_x + 1e-15).powi(2))
                    * r
                    * r;

                for &(theta_phi, w_ang) in &angular {
                    let xyz = [
                        center[0] + r * theta_phi[0],
                        center[1] + r * theta_phi[1],
                        center[2] + r * theta_phi[2],
                    ];

                    // Becke partition weight
                    let w_becke = becke_weights(&xyz, atom_idx, positions_bohr, atomic_numbers);

                    let weight = w_rad * w_ang * w_becke;

                    if weight.abs() > 1e-30 {
                        points.push(GridPoint {
                            xyz,
                            weight,
                            atom_index: atom_idx,
                        });
                    }
                }
            }
        }

        let n_points = points.len();
        Self { points, n_points }
    }
}

/// Bragg-Slater atomic radii in Bohr.
fn bragg_slater_radius(z: u8) -> f64 {
    // Radii from J.C. Slater, J. Chem. Phys. 41, 3199 (1964)
    // Converted to Bohr (1 Å = 1.8897259886 Bohr)
    let ang_to_bohr = 1.889_725_988_6;
    let r_angstrom = match z {
        1 => 0.25,
        2 => 0.31,
        3 => 1.45,
        4 => 1.05,
        5 => 0.85,
        6 => 0.70,
        7 => 0.65,
        8 => 0.60,
        9 => 0.50,
        10 => 0.38,
        11 => 1.80,
        12 => 1.50,
        13 => 1.25,
        14 => 1.10,
        15 => 1.00,
        16 => 1.00,
        17 => 1.00,
        18 => 0.71,
        _ => 1.50, // generic fallback
    };
    r_angstrom * ang_to_bohr
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grid_builds_for_single_atom() {
        let grid = MolecularGrid::build(&[1], &[[0.0, 0.0, 0.0]], GridQuality::Coarse);
        assert!(grid.n_points > 0, "grid should have points");
        assert_eq!(grid.points.len(), grid.n_points);
    }

    #[test]
    fn grid_quality_ordering() {
        let g1 = MolecularGrid::build(&[1], &[[0.0, 0.0, 0.0]], GridQuality::Coarse);
        let g2 = MolecularGrid::build(&[1], &[[0.0, 0.0, 0.0]], GridQuality::Fine);
        assert!(
            g2.n_points > g1.n_points,
            "Fine grid should have more points than Coarse"
        );
    }

    #[test]
    fn grid_weights_are_positive() {
        let grid = MolecularGrid::build(&[8], &[[0.0, 0.0, 0.0]], GridQuality::Medium);
        for gp in &grid.points {
            assert!(gp.weight >= 0.0, "grid weights must be non-negative");
        }
    }

    #[test]
    fn two_atom_grid_has_both_centers() {
        let grid = MolecularGrid::build(
            &[1, 1],
            &[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            GridQuality::Coarse,
        );
        let has_atom0 = grid.points.iter().any(|p| p.atom_index == 0);
        let has_atom1 = grid.points.iter().any(|p| p.atom_index == 1);
        assert!(
            has_atom0 && has_atom1,
            "grid should sample around both atoms"
        );
    }
}
