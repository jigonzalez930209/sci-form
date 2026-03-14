use nalgebra::{DMatrix, Vector3};

/// Calculate the dihedral (torsion) angle for 4 ordered points: A-B-C-D.
/// The output is in radians, in the range [-pi, pi].
pub fn calculate_dihedral_angle(
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
    coords: &DMatrix<f32>,
) -> f32 {
    let dim = coords.ncols();
    assert!(dim >= 3, "Coordinates must have at least 3 dimensions");

    // Extract position vectors
    let p1 = Vector3::new(coords[(idx1, 0)], coords[(idx1, 1)], coords[(idx1, 2)]);
    let p2 = Vector3::new(coords[(idx2, 0)], coords[(idx2, 1)], coords[(idx2, 2)]);
    let p3 = Vector3::new(coords[(idx3, 0)], coords[(idx3, 1)], coords[(idx3, 2)]);
    let p4 = Vector3::new(coords[(idx4, 0)], coords[(idx4, 1)], coords[(idx4, 2)]);

    // Bond vectors
    let b1 = p2 - p1;
    let b2 = p3 - p2;
    let b3 = p4 - p3;

    // Normal vectors to the planes
    let n1 = b1.cross(&b2);
    let n2 = b2.cross(&b3);

    // Frame vectors
    let n1_unit = n1.normalize();
    let n2_unit = n2.normalize();

    // Check for degenerate cases (collinear atoms) where normalization fails
    if n1.norm() < 1e-6 || n2.norm() < 1e-6 {
        return 0.0;
    }

    // Vector orthogonal to n1 and b2
    let m = n1_unit.cross(&b2.normalize());

    // Compute the components of n2_unit
    let x = n1_unit.dot(&n2_unit);
    let y = m.dot(&n2_unit);

    // Dihedral angle
    y.atan2(x)
}

/// Calculate the ETKDG torsion energy penalty.
/// This is a simplified mathematical mock representing the preference of a dihedral
/// to be close to an experimental preferred angle `preferred_angle_rad`.
///
/// The function assumes a harmonic well around the generic preference.
/// Real ETKDG in RDKit builds heavily on empirical datasets over hybridization types.
pub fn calculate_torsion_penalty(
    current_angle_rad: f32,
    preferred_angle_rad: f32,
    force_constant: f32,
) -> f32 {
    // We use a cosine penalty function to properly account for angular periodicity
    // E = k * (1 - cos(current - preferred))
    force_constant * (1.0 - (current_angle_rad - preferred_angle_rad).cos())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    #[test]
    fn test_dihedral_angle_trans() {
        // Anti-periplanar (trans) should be ~PI or -PI (180 deg)
        let mut coords = DMatrix::from_element(4, 3, 0.0);
        coords.set_row(0, &nalgebra::RowVector3::new(0.0, 1.0, 0.0));
        coords.set_row(1, &nalgebra::RowVector3::new(0.0, 0.0, 0.0));
        coords.set_row(2, &nalgebra::RowVector3::new(1.0, 0.0, 0.0));
        coords.set_row(3, &nalgebra::RowVector3::new(1.0, -1.0, 0.0));

        let angle = calculate_dihedral_angle(0, 1, 2, 3, &coords);
        let angle_deg = angle.to_degrees();
        assert!((angle_deg.abs() - 180.0).abs() < 1e-4);
    }

    #[test]
    fn test_dihedral_angle_gauche() {
        // Gauche should be ~ 60 or -60 deg
        let mut coords = DMatrix::from_element(4, 3, 0.0);
        coords.set_row(0, &nalgebra::RowVector3::new(0.0, 1.0, 0.0));
        coords.set_row(1, &nalgebra::RowVector3::new(0.0, 0.0, 0.0));
        coords.set_row(2, &nalgebra::RowVector3::new(1.0, 0.0, 0.0));
        coords.set_row(3, &nalgebra::RowVector3::new(1.0, 0.0, 1.0)); // Out of plane

        let angle = calculate_dihedral_angle(0, 1, 2, 3, &coords);
        let angle_deg = angle.to_degrees();
        assert!((angle_deg - -90.0).abs() < 1e-4); // Actually Orthogonal (90 degrees) for this manual structure
    }

    #[test]
    fn test_torsion_penalty() {
        // Perfect angle = 0 penalty
        let penalty = calculate_torsion_penalty(PI, PI, 10.0);
        assert!((penalty - 0.0).abs() < 1e-5);

        // Max penalty at 180 deviation: 1 - cos(180) = 1 - (-1) = 2
        // Penalty = 10 * 2 = 20
        let penalty2 = calculate_torsion_penalty(PI, 0.0, 10.0);
        assert!((penalty2 - 20.0).abs() < 1e-5);
    }
}
