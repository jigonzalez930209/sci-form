//! Kabsch algorithm for optimal rotation alignment and RMSD computation.
//!
//! Provides a clean f64 public API for molecular alignment.

/// Result of a Kabsch alignment.
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    /// RMSD after optimal alignment (Å).
    pub rmsd: f64,
    /// 3×3 rotation matrix (row-major).
    pub rotation: [[f64; 3]; 3],
    /// Translation applied to `coords` centroid.
    pub translation: [f64; 3],
    /// Aligned coordinates: flat [x0,y0,z0, x1,...].
    pub aligned_coords: Vec<f64>,
}

/// Compute RMSD between two conformers after Kabsch alignment.
///
/// `coords`: flat [x0,y0,z0, x1,y1,z1,...] mobile structure.
/// `reference`: flat [x0,y0,z0,...] reference structure.
/// Both must have the same number of atoms.
pub fn compute_rmsd(coords: &[f64], reference: &[f64]) -> f64 {
    align_coordinates(coords, reference).rmsd
}

/// Kabsch alignment: find optimal rotation mapping `coords` onto `reference`.
///
/// `coords`, `reference`: flat [x0,y0,z0, x1,y1,z1,...] in Å.
pub fn align_coordinates(coords: &[f64], reference: &[f64]) -> AlignmentResult {
    if coords.len() != reference.len() || coords.len() % 3 != 0 {
        return AlignmentResult {
            rmsd: f64::NAN,
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            translation: [0.0; 3],
            aligned_coords: coords.to_vec(),
        };
    }
    let n = coords.len() / 3;

    if n == 0 {
        return AlignmentResult {
            rmsd: 0.0,
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            translation: [0.0; 3],
            aligned_coords: Vec::new(),
        };
    }

    // Compute centroids
    let mut c1 = [0.0f64; 3];
    let mut c2 = [0.0f64; 3];
    for i in 0..n {
        for k in 0..3 {
            c1[k] += coords[i * 3 + k];
            c2[k] += reference[i * 3 + k];
        }
    }
    for k in 0..3 {
        c1[k] /= n as f64;
        c2[k] /= n as f64;
    }

    // Build H = Σ (p_i - c1)(q_i - c2)^T  (3×3)
    let mut h = [[0.0f64; 3]; 3];
    for i in 0..n {
        let p = [
            coords[i * 3] - c1[0],
            coords[i * 3 + 1] - c1[1],
            coords[i * 3 + 2] - c1[2],
        ];
        let q = [
            reference[i * 3] - c2[0],
            reference[i * 3 + 1] - c2[1],
            reference[i * 3 + 2] - c2[2],
        ];
        for r in 0..3 {
            for c in 0..3 {
                h[r][c] += p[r] * q[c];
            }
        }
    }

    // SVD via nalgebra
    let h_mat = nalgebra::Matrix3::new(
        h[0][0], h[0][1], h[0][2], h[1][0], h[1][1], h[1][2], h[2][0], h[2][1], h[2][2],
    );
    let svd = h_mat.svd(true, true);
    let (u, v_t) = match (svd.u, svd.v_t) {
        (Some(u), Some(v_t)) => (u, v_t),
        _ => {
            // SVD failed — return identity alignment with raw RMSD
            let mut sum_sq = 0.0;
            for i in 0..coords.len() {
                let diff = coords[i] - reference[i];
                sum_sq += diff * diff;
            }
            return AlignmentResult {
                rmsd: (sum_sq / n as f64).sqrt(),
                rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                translation: [0.0; 3],
                aligned_coords: coords.to_vec(),
            };
        }
    };
    let v = v_t.transpose();

    // Handle reflection: correct the sign of the smallest singular value
    // to ensure a proper rotation (det(R) = +1).
    // For coplanar molecules (one singular value ≈ 0), this prevents
    // an improper rotation (reflection) from being selected.
    let mut d = nalgebra::Matrix3::<f64>::identity();
    if (v * u.transpose()).determinant() < 0.0 {
        d[(2, 2)] = -1.0;
    }
    let r_mat = v * d * u.transpose();

    // Build rotation as row-major array
    let rotation = [
        [r_mat[(0, 0)], r_mat[(0, 1)], r_mat[(0, 2)]],
        [r_mat[(1, 0)], r_mat[(1, 1)], r_mat[(1, 2)]],
        [r_mat[(2, 0)], r_mat[(2, 1)], r_mat[(2, 2)]],
    ];

    let translation = [c2[0] - c1[0], c2[1] - c1[1], c2[2] - c1[2]];

    // Apply rotation and compute RMSD
    let mut aligned = vec![0.0f64; coords.len()];
    let mut sum_sq = 0.0;
    for i in 0..n {
        let p = [
            coords[i * 3] - c1[0],
            coords[i * 3 + 1] - c1[1],
            coords[i * 3 + 2] - c1[2],
        ];
        for k in 0..3 {
            let rotated = r_mat[(k, 0)] * p[0] + r_mat[(k, 1)] * p[1] + r_mat[(k, 2)] * p[2];
            aligned[i * 3 + k] = rotated + c2[k];
        }
        for k in 0..3 {
            let diff = aligned[i * 3 + k] - reference[i * 3 + k];
            sum_sq += diff * diff;
        }
    }
    let rmsd = (sum_sq / n as f64).sqrt();

    AlignmentResult {
        rmsd,
        rotation,
        translation,
        aligned_coords: aligned,
    }
}

/// Quaternion-based optimal rotation alignment (Coutsias et al. 2004).
///
/// Uses the quaternion eigenvector method which is numerically more stable
/// than SVD for near-degenerate cases.  Produces the same result as Kabsch
/// but avoids explicit SVD decomposition.
pub fn align_quaternion(coords: &[f64], reference: &[f64]) -> AlignmentResult {
    if coords.len() != reference.len() || coords.len() % 3 != 0 {
        return AlignmentResult {
            rmsd: f64::NAN,
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            translation: [0.0; 3],
            aligned_coords: coords.to_vec(),
        };
    }
    let n = coords.len() / 3;

    if n == 0 {
        return AlignmentResult {
            rmsd: 0.0,
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            translation: [0.0; 3],
            aligned_coords: Vec::new(),
        };
    }

    // Centroids
    let mut c1 = [0.0f64; 3];
    let mut c2 = [0.0f64; 3];
    for i in 0..n {
        for k in 0..3 {
            c1[k] += coords[i * 3 + k];
            c2[k] += reference[i * 3 + k];
        }
    }
    for k in 0..3 {
        c1[k] /= n as f64;
        c2[k] /= n as f64;
    }

    // Build cross-covariance elements: R_ij = Σ (p_i - c1_i)(q_j - c2_j)
    let mut r = [[0.0f64; 3]; 3];
    for i in 0..n {
        let p = [
            coords[i * 3] - c1[0],
            coords[i * 3 + 1] - c1[1],
            coords[i * 3 + 2] - c1[2],
        ];
        let q = [
            reference[i * 3] - c2[0],
            reference[i * 3 + 1] - c2[1],
            reference[i * 3 + 2] - c2[2],
        ];
        for a in 0..3 {
            for b in 0..3 {
                r[a][b] += p[a] * q[b];
            }
        }
    }

    // Build the 4×4 symmetric key matrix (Davenport/Coutsias):
    //   F = [[ Sxx+Syy+Szz, Syz-Szy, Szx-Sxz, Sxy-Syx ],
    //        [ Syz-Szy, Sxx-Syy-Szz, Sxy+Syx, Szx+Sxz ],
    //        [ Szx-Sxz, Sxy+Syx, -Sxx+Syy-Szz, Syz+Szy ],
    //        [ Sxy-Syx, Szx+Sxz, Syz+Szy, -Sxx-Syy+Szz ]]
    let sxx = r[0][0];
    let sxy = r[0][1];
    let sxz = r[0][2];
    let syx = r[1][0];
    let syy = r[1][1];
    let syz = r[1][2];
    let szx = r[2][0];
    let szy = r[2][1];
    let szz = r[2][2];

    let f = nalgebra::Matrix4::new(
        sxx + syy + szz,
        syz - szy,
        szx - sxz,
        sxy - syx,
        syz - szy,
        sxx - syy - szz,
        sxy + syx,
        szx + sxz,
        szx - sxz,
        sxy + syx,
        -sxx + syy - szz,
        syz + szy,
        sxy - syx,
        szx + sxz,
        syz + szy,
        -sxx - syy + szz,
    );

    // The optimal rotation quaternion is the eigenvector of F with the largest eigenvalue
    let eig = f.symmetric_eigen();
    let mut best_idx = 0;
    let mut best_val = eig.eigenvalues[0];
    for i in 1..4 {
        if eig.eigenvalues[i] > best_val {
            best_val = eig.eigenvalues[i];
            best_idx = i;
        }
    }

    let q0 = eig.eigenvectors[(0, best_idx)];
    let q1 = eig.eigenvectors[(1, best_idx)];
    let q2 = eig.eigenvectors[(2, best_idx)];
    let q3 = eig.eigenvectors[(3, best_idx)];

    // quaternion → rotation matrix
    let rotation = [
        [
            q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3,
            2.0 * (q1 * q2 - q0 * q3),
            2.0 * (q1 * q3 + q0 * q2),
        ],
        [
            2.0 * (q1 * q2 + q0 * q3),
            q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3,
            2.0 * (q2 * q3 - q0 * q1),
        ],
        [
            2.0 * (q1 * q3 - q0 * q2),
            2.0 * (q2 * q3 + q0 * q1),
            q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3,
        ],
    ];

    let translation = [c2[0] - c1[0], c2[1] - c1[1], c2[2] - c1[2]];

    // Apply rotation and compute RMSD
    let mut aligned = vec![0.0f64; coords.len()];
    let mut sum_sq = 0.0;
    for i in 0..n {
        let p = [
            coords[i * 3] - c1[0],
            coords[i * 3 + 1] - c1[1],
            coords[i * 3 + 2] - c1[2],
        ];
        for k in 0..3 {
            let rotated = rotation[k][0] * p[0] + rotation[k][1] * p[1] + rotation[k][2] * p[2];
            aligned[i * 3 + k] = rotated + c2[k];
        }
        for k in 0..3 {
            let diff = aligned[i * 3 + k] - reference[i * 3 + k];
            sum_sq += diff * diff;
        }
    }
    let rmsd = (sum_sq / n as f64).sqrt();

    AlignmentResult {
        rmsd,
        rotation,
        translation,
        aligned_coords: aligned,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identical_zero_rmsd() {
        let coords = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let rmsd = compute_rmsd(&coords, &coords);
        assert!(rmsd < 1e-10);
    }

    #[test]
    fn test_translated_zero_rmsd() {
        // Pure translation → RMSD should be ~0 after alignment.
        let reference = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let coords: Vec<f64> = reference.iter().map(|x| x + 5.0).collect();
        let rmsd = compute_rmsd(&coords, &reference);
        assert!(rmsd < 1e-10, "got rmsd = {rmsd}");
    }

    #[test]
    fn test_rotation_90deg_z() {
        // 90° rotation around Z-axis → RMSD should be ~0 after alignment.
        let reference = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0];
        // Rotate 90° around Z: (x,y) -> (-y,x)
        let rotated = vec![0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, 0.0];
        let rmsd = compute_rmsd(&rotated, &reference);
        assert!(rmsd < 1e-10, "got rmsd = {rmsd}");
    }

    #[test]
    fn test_known_rmsd() {
        // Slightly perturbed structure → nonzero RMSD.
        let reference = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let perturbed = vec![0.1, 0.0, 0.0, 1.0, 0.1, 0.0, 0.0, 1.0, 0.1];
        let rmsd = compute_rmsd(&perturbed, &reference);
        assert!(rmsd > 0.01);
        assert!(rmsd < 1.0);
    }

    #[test]
    fn test_aligned_coords_returned() {
        let reference = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let coords: Vec<f64> = reference.iter().map(|x| x + 10.0).collect();
        let result = align_coordinates(&coords, &reference);
        assert_eq!(result.aligned_coords.len(), 9);
        // Aligned should be close to reference
        for i in 0..9 {
            assert!(
                (result.aligned_coords[i] - reference[i]).abs() < 1e-8,
                "mismatch at index {i}"
            );
        }
    }

    #[test]
    fn test_reflection_handling() {
        // Mirror image (reflection): Kabsch should handle determinant < 0.
        let reference = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let reflected = vec![-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let result = align_coordinates(&reflected, &reference);
        // Should still give a valid RMSD (may not be zero for true reflection)
        assert!(result.rmsd.is_finite());
    }

    #[test]
    fn test_quaternion_identical() {
        let coords = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let result = align_quaternion(&coords, &coords);
        assert!(result.rmsd < 1e-10);
    }

    #[test]
    fn test_quaternion_translated() {
        let reference = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let coords: Vec<f64> = reference.iter().map(|x| x + 5.0).collect();
        let result = align_quaternion(&coords, &reference);
        assert!(result.rmsd < 1e-10, "got rmsd = {}", result.rmsd);
    }

    #[test]
    fn test_quaternion_rotated_90() {
        let reference = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0];
        let rotated = vec![0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, 0.0];
        let result = align_quaternion(&rotated, &reference);
        assert!(result.rmsd < 1e-10, "got rmsd = {}", result.rmsd);
    }

    #[test]
    fn test_quaternion_matches_kabsch() {
        // Both methods should give the same RMSD.
        let reference = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.5, 0.5, 1.0];
        let perturbed = vec![
            0.1, -0.05, 0.02, 1.1, 0.1, -0.05, -0.1, 0.9, 0.1, 0.6, 0.4, 1.1,
        ];

        let kabsch = align_coordinates(&perturbed, &reference);
        let quat = align_quaternion(&perturbed, &reference);

        assert!(
            (kabsch.rmsd - quat.rmsd).abs() < 1e-8,
            "Kabsch RMSD = {}, Quaternion RMSD = {}",
            kabsch.rmsd,
            quat.rmsd,
        );

        // Aligned coords should also match
        for i in 0..reference.len() {
            assert!(
                (kabsch.aligned_coords[i] - quat.aligned_coords[i]).abs() < 1e-6,
                "aligned mismatch at {}: {:.8} vs {:.8}",
                i,
                kabsch.aligned_coords[i],
                quat.aligned_coords[i],
            );
        }
    }
}
