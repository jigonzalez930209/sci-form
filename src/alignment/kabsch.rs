//! Kabsch algorithm for optimal rotation alignment and RMSD computation.
//!
//! Provides a clean f64 public API for molecular alignment.
//!
//! The superposition follows Kabsch, *Acta Cryst.* **1976**, A32, 922: the
//! centroids are removed, the optimal proper rotation is obtained from the
//! singular-value decomposition of the covariance matrix, and the smallest
//! singular value is sign-flipped when needed so that the result is a proper
//! rotation (`det(R) = +1`) rather than a reflection. The alternative
//! quaternion route follows Coutsias, Dill & Church, *Protein Sci.* **2004**,
//! 13, 275, and is numerically more stable for near-degenerate cases.
//!
//! Coordinates are flat `[x0, y0, z0, x1, ...]` arrays in Å; the returned
//! RMSD and translation are in Å.

/// Result of a Kabsch alignment.
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    /// RMSD after optimal alignment (Å).
    ///
    /// `f64::NAN` is never produced.  When the alignment could not be
    /// performed (mismatched input lengths) the RMSD is the sentinel
    /// [`AlignmentResult::INVALID_RMSD`] and [`AlignmentResult::error`] holds
    /// the reason (blocker A-H6).
    pub rmsd: f64,
    /// 3×3 rotation matrix (row-major).
    pub rotation: [[f64; 3]; 3],
    /// Translation applied to `coords` centroid, in Å.
    ///
    /// The aligned coordinates are `rotation · (p - c_mobile) + c_reference`,
    /// i.e. this is the offset `c_reference - c_mobile` expressed in the
    /// reference frame.
    pub translation: [f64; 3],
    /// Aligned coordinates: flat [x0,y0,z0, x1,...] in Å.
    pub aligned_coords: Vec<f64>,
    /// Reason why the alignment failed, if any.
    ///
    /// `None` means the alignment succeeded.
    pub error: Option<String>,
}

impl AlignmentResult {
    /// Sentinel RMSD used when the alignment could not be performed.
    ///
    /// A negative RMSD is impossible for a real alignment, so it can never
    /// be confused with "the structures do not overlap" (RMSD = 0).
    pub const INVALID_RMSD: f64 = -1.0;

    /// Whether the alignment failed.
    pub fn is_valid(&self) -> bool {
        self.error.is_none()
    }
}

/// Compute RMSD between two conformers after Kabsch alignment.
///
/// `coords`: flat [x0,y0,z0, x1,y1,z1,...] mobile structure.
/// `reference`: flat [x0,y0,z0,...] reference structure.
/// Both must have the same number of atoms.
///
/// A length mismatch yields [`AlignmentResult::INVALID_RMSD`]; use
/// [`try_compute_rmsd`] to get an explicit error instead.
pub fn compute_rmsd(coords: &[f64], reference: &[f64]) -> f64 {
    align_coordinates(coords, reference).rmsd
}

/// Fallible variant of [`compute_rmsd`]: returns an explicit error when the
/// two structures cannot be compared (mismatched lengths).
///
/// # Errors
///
/// Returns `Err` if the two coordinate arrays have different lengths, or if
/// their length is not a multiple of three (a partial `x, y, z` triplet).
pub fn try_compute_rmsd(coords: &[f64], reference: &[f64]) -> Result<f64, String> {
    try_align_coordinates(coords, reference).map(|result| result.rmsd)
}

/// Kabsch alignment: find optimal rotation mapping `coords` onto `reference`.
///
/// `coords`, `reference`: flat [x0,y0,z0, x1,y1,z1,...] in Å.
///
/// When the two structures have a different length the alignment cannot be
/// performed.  The result then carries [`AlignmentResult::INVALID_RMSD`] and
/// an [`AlignmentResult::error`] message; use [`try_align_coordinates`] to
/// obtain a `Result` instead (blocker A-H6: the previous implementation
/// returned `f64::NAN`, which the WASM/CLI bindings serialised as invalid
/// JSON).
pub fn align_coordinates(coords: &[f64], reference: &[f64]) -> AlignmentResult {
    match try_align_coordinates(coords, reference) {
        Ok(result) => result,
        Err(reason) => AlignmentResult {
            rmsd: AlignmentResult::INVALID_RMSD,
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            translation: [0.0; 3],
            aligned_coords: coords.to_vec(),
            error: Some(reason),
        },
    }
}

/// Fallible variant of [`align_coordinates`].
///
/// # Errors
///
/// Returns `Err` if the two coordinate arrays have different lengths, or if
/// their length is not a multiple of three (a partial `x, y, z` triplet).
pub fn try_align_coordinates(coords: &[f64], reference: &[f64]) -> Result<AlignmentResult, String> {
    if coords.len() != reference.len() {
        return Err(format!(
            "cannot align coordinate sets of different length: {} vs {} values",
            coords.len(),
            reference.len()
        ));
    }
    if !coords.len().is_multiple_of(3) {
        return Err(format!(
            "coordinate length {} is not a multiple of 3 (x, y, z triplets)",
            coords.len()
        ));
    }
    Ok(align_coordinates_impl(coords, reference))
}

fn align_coordinates_impl(coords: &[f64], reference: &[f64]) -> AlignmentResult {
    let n = coords.len() / 3;

    if n == 0 {
        return AlignmentResult {
            rmsd: 0.0,
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            translation: [0.0; 3],
            aligned_coords: Vec::new(),
            error: None,
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
                error: None,
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
        error: None,
    }
}

/// Quaternion-based optimal rotation alignment (Coutsias et al. 2004).
///
/// Uses the quaternion eigenvector method which is numerically more stable
/// than SVD for near-degenerate cases.  Produces the same result as Kabsch
/// but avoids explicit SVD decomposition.
///
/// A length mismatch yields [`AlignmentResult::INVALID_RMSD`] plus an
/// [`AlignmentResult::error`] message; use [`try_align_quaternion`] for a
/// `Result` (blocker A-H6).
pub fn align_quaternion(coords: &[f64], reference: &[f64]) -> AlignmentResult {
    match try_align_quaternion(coords, reference) {
        Ok(result) => result,
        Err(reason) => AlignmentResult {
            rmsd: AlignmentResult::INVALID_RMSD,
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            translation: [0.0; 3],
            aligned_coords: coords.to_vec(),
            error: Some(reason),
        },
    }
}

/// Fallible variant of [`align_quaternion`].
///
/// # Errors
///
/// Returns `Err` if the two coordinate arrays have different lengths, or if
/// their length is not a multiple of three (a partial `x, y, z` triplet).
pub fn try_align_quaternion(coords: &[f64], reference: &[f64]) -> Result<AlignmentResult, String> {
    if coords.len() != reference.len() {
        return Err(format!(
            "cannot align coordinate sets of different length: {} vs {} values",
            coords.len(),
            reference.len()
        ));
    }
    if !coords.len().is_multiple_of(3) {
        return Err(format!(
            "coordinate length {} is not a multiple of 3 (x, y, z triplets)",
            coords.len()
        ));
    }
    Ok(align_quaternion_impl(coords, reference))
}

fn align_quaternion_impl(coords: &[f64], reference: &[f64]) -> AlignmentResult {
    let n = coords.len() / 3;

    if n == 0 {
        return AlignmentResult {
            rmsd: 0.0,
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            translation: [0.0; 3],
            aligned_coords: Vec::new(),
            error: None,
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
        error: None,
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

    /// A-H6 regression test: a length mismatch must never be reported as a
    /// `NaN` RMSD (which the WASM/CLI bindings serialised as invalid JSON).
    #[test]
    fn test_length_mismatch_is_an_error_not_a_nan() {
        let coords = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let reference = vec![0.0, 0.0, 0.0];

        let result = align_coordinates(&coords, &reference);
        assert!(
            result.rmsd.is_finite(),
            "a length mismatch must not produce NaN or infinity"
        );
        assert_eq!(result.rmsd, AlignmentResult::INVALID_RMSD);
        assert!(!result.is_valid());
        let error = result.error.expect("the reason must be reported");
        assert!(
            error.contains("different length"),
            "the error must explain the mismatch, got: {error}"
        );

        // The fallible API reports the error directly.
        assert!(try_align_coordinates(&coords, &reference).is_err());
        assert!(try_compute_rmsd(&coords, &reference).is_err());

        // A length that is not a multiple of 3 is an error as well.
        let odd = vec![0.0, 0.0, 0.0, 1.0];
        let error = try_align_coordinates(&odd, &odd)
            .expect_err("a partial coordinate triplet must be rejected");
        assert!(
            error.contains("multiple of 3"),
            "the error must explain the mismatch, got: {error}"
        );

        // The quaternion variant behaves the same way.
        let result = align_quaternion(&coords, &reference);
        assert!(result.rmsd.is_finite());
        assert_eq!(result.rmsd, AlignmentResult::INVALID_RMSD);
        assert!(result.error.is_some());

        // A valid alignment carries no error.
        let valid = align_coordinates(&coords, &coords);
        assert!(valid.error.is_none());
        assert!(valid.is_valid());
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
