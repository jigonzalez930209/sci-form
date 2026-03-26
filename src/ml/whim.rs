//! WHIM (Weighted Holistic Invariant Molecular) descriptors.
//!
//! 3D molecular descriptors based on statistical indices of the atomic
//! coordinate distribution projected onto principal axes. Includes:
//! - Directional WHIM (per principal axis): λ, ν, γ, η
//! - Global WHIM: T (total spread), A (anisotropy), D (directional), K (kurtosis)
//!
//! Reference: Todeschini & Gramatica, 3D Molecular Descriptors (1997).

use serde::{Deserialize, Serialize};

/// WHIM descriptor result for a molecule.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WhimDescriptors {
    /// Eigenvalues of the weighted covariance matrix (λ1 ≥ λ2 ≥ λ3).
    pub eigenvalues: [f64; 3],
    /// Directional WHIM: spread along each principal axis.
    pub lambda: [f64; 3],
    /// Directional shape: ν_i = λ_i / (λ1 + λ2 + λ3).
    pub nu: [f64; 3],
    /// Directional symmetry: γ_i (third moment / λ_i^1.5).
    pub gamma: [f64; 3],
    /// Directional kurtosis: η_i (fourth moment / λ_i^2).
    pub eta: [f64; 3],
    /// Total spread: T = λ1 + λ2 + λ3.
    pub total_spread: f64,
    /// Anisotropy: A = 1 − (3 * product / T^2)^(1/3).
    /// 0 = isotropic (sphere), 1 = maximally anisotropic.
    pub anisotropy: f64,
    /// Directional index: D = product of directional symmetries.
    pub directional: f64,
    /// Global kurtosis: K = sum of η_i.
    pub kurtosis: f64,
}

/// Atomic weighting scheme for WHIM descriptors.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum WhimWeighting {
    /// Unit weights (unweighted).
    Unit,
    /// Atomic mass weights.
    Mass,
    /// Van der Waals volume weights.
    Volume,
    /// Pauling electronegativity weights.
    Electronegativity,
    /// Atomic polarizability weights.
    Polarizability,
}

/// Compute WHIM descriptors from 3D coordinates with specified weighting.
pub fn compute_whim(
    elements: &[u8],
    positions: &[[f64; 3]],
    weighting: WhimWeighting,
) -> WhimDescriptors {
    let n = elements.len().min(positions.len());
    if n == 0 {
        return WhimDescriptors {
            eigenvalues: [0.0; 3],
            lambda: [0.0; 3],
            nu: [1.0 / 3.0; 3],
            gamma: [0.0; 3],
            eta: [0.0; 3],
            total_spread: 0.0,
            anisotropy: 0.0,
            directional: 0.0,
            kurtosis: 0.0,
        };
    }

    // Compute weights
    let weights: Vec<f64> = elements[..n]
        .iter()
        .map(|&z| atom_weight(z, weighting))
        .collect();
    let w_sum: f64 = weights.iter().sum();
    if w_sum < 1e-12 {
        return WhimDescriptors {
            eigenvalues: [0.0; 3],
            lambda: [0.0; 3],
            nu: [1.0 / 3.0; 3],
            gamma: [0.0; 3],
            eta: [0.0; 3],
            total_spread: 0.0,
            anisotropy: 0.0,
            directional: 0.0,
            kurtosis: 0.0,
        };
    }

    // Weighted centroid
    let mut centroid = [0.0f64; 3];
    for i in 0..n {
        for d in 0..3 {
            centroid[d] += weights[i] * positions[i][d];
        }
    }
    for d in 0..3 {
        centroid[d] /= w_sum;
    }

    // Weighted covariance matrix (3×3 symmetric)
    let mut cov = [[0.0f64; 3]; 3];
    for i in 0..n {
        let dx = [
            positions[i][0] - centroid[0],
            positions[i][1] - centroid[1],
            positions[i][2] - centroid[2],
        ];
        for r in 0..3 {
            for c in r..3 {
                cov[r][c] += weights[i] * dx[r] * dx[c];
            }
        }
    }
    for r in 0..3 {
        for c in r..3 {
            cov[r][c] /= w_sum;
            if c != r {
                cov[c][r] = cov[r][c];
            }
        }
    }

    // Eigendecomposition of 3×3 symmetric matrix (Jacobi iteration)
    let (eigenvalues, eigenvectors) = eigendecompose_3x3(&cov);

    // Sort eigenvalues descending
    let mut sorted: Vec<(f64, usize)> = eigenvalues
        .iter()
        .enumerate()
        .map(|(i, &v)| (v, i))
        .collect();
    sorted.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

    let lambda = [
        sorted[0].0.max(0.0),
        sorted[1].0.max(0.0),
        sorted[2].0.max(0.0),
    ];

    // Project onto principal axes and compute higher moments
    let mut gamma = [0.0f64; 3];
    let mut eta = [0.0f64; 3];

    for axis in 0..3 {
        let ax_idx = sorted[axis].1;
        let ev = &eigenvectors[ax_idx];
        let lam = lambda[axis];

        if lam < 1e-12 {
            continue;
        }

        let mut m3 = 0.0;
        let mut m4 = 0.0;
        for i in 0..n {
            let dx = [
                positions[i][0] - centroid[0],
                positions[i][1] - centroid[1],
                positions[i][2] - centroid[2],
            ];
            let proj = dx[0] * ev[0] + dx[1] * ev[1] + dx[2] * ev[2];
            let p2 = proj * proj;
            m3 += weights[i] * proj * p2;
            m4 += weights[i] * p2 * p2;
        }
        m3 /= w_sum;
        m4 /= w_sum;

        gamma[axis] = m3 / lam.powf(1.5);
        eta[axis] = m4 / (lam * lam);
    }

    let total_spread = lambda[0] + lambda[1] + lambda[2];
    let nu = if total_spread > 1e-12 {
        [
            lambda[0] / total_spread,
            lambda[1] / total_spread,
            lambda[2] / total_spread,
        ]
    } else {
        [1.0 / 3.0; 3]
    };

    let product = lambda[0] * lambda[1] * lambda[2];
    let anisotropy = if total_spread > 1e-12 {
        let ratio = (3.0 * product) / (total_spread * total_spread);
        1.0 - ratio.cbrt().min(1.0)
    } else {
        0.0
    };

    let directional = nu[0] * nu[1] * nu[2];
    let kurtosis = eta[0] + eta[1] + eta[2];

    WhimDescriptors {
        eigenvalues: lambda,
        lambda,
        nu,
        gamma,
        eta,
        total_spread,
        anisotropy,
        directional,
        kurtosis,
    }
}

/// Eigendecompose a 3×3 symmetric matrix using Jacobi rotations.
/// Returns (eigenvalues, eigenvectors) where eigenvectors[i] is the i-th eigenvector.
fn eigendecompose_3x3(mat: &[[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut a = *mat;
    let mut v = [[0.0f64; 3]; 3];
    for i in 0..3 {
        v[i][i] = 1.0;
    }

    for _ in 0..100 {
        // Find largest off-diagonal element
        let mut p = 0;
        let mut q = 1;
        let mut max_val = a[0][1].abs();
        for i in 0..3 {
            for j in (i + 1)..3 {
                if a[i][j].abs() > max_val {
                    max_val = a[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }

        if max_val < 1e-14 {
            break;
        }

        let theta = if (a[p][p] - a[q][q]).abs() < 1e-14 {
            std::f64::consts::FRAC_PI_4
        } else {
            0.5 * (2.0 * a[p][q] / (a[p][p] - a[q][q])).atan()
        };

        let c = theta.cos();
        let s = theta.sin();

        // Rotate matrix
        let mut new_a = a;
        new_a[p][p] = c * c * a[p][p] + 2.0 * s * c * a[p][q] + s * s * a[q][q];
        new_a[q][q] = s * s * a[p][p] - 2.0 * s * c * a[p][q] + c * c * a[q][q];
        new_a[p][q] = 0.0;
        new_a[q][p] = 0.0;

        for r in 0..3 {
            if r != p && r != q {
                new_a[r][p] = c * a[r][p] + s * a[r][q];
                new_a[p][r] = new_a[r][p];
                new_a[r][q] = -s * a[r][p] + c * a[r][q];
                new_a[q][r] = new_a[r][q];
            }
        }
        a = new_a;

        // Rotate eigenvectors
        let mut new_v = v;
        for r in 0..3 {
            new_v[r][p] = c * v[r][p] + s * v[r][q];
            new_v[r][q] = -s * v[r][p] + c * v[r][q];
        }
        v = new_v;
    }

    let eigenvalues = [a[0][0], a[1][1], a[2][2]];
    // Transpose v so eigenvectors[i] is column i
    let eigenvectors = [
        [v[0][0], v[1][0], v[2][0]],
        [v[0][1], v[1][1], v[2][1]],
        [v[0][2], v[1][2], v[2][2]],
    ];

    (eigenvalues, eigenvectors)
}

/// Get atomic weight for WHIM weighting scheme.
fn atom_weight(z: u8, scheme: WhimWeighting) -> f64 {
    match scheme {
        WhimWeighting::Unit => 1.0,
        WhimWeighting::Mass => atomic_mass(z),
        WhimWeighting::Volume => vdw_volume(z),
        WhimWeighting::Electronegativity => electronegativity(z),
        WhimWeighting::Polarizability => polarizability(z),
    }
}

fn atomic_mass(z: u8) -> f64 {
    match z {
        1 => 1.008,
        5 => 10.81,
        6 => 12.011,
        7 => 14.007,
        8 => 15.999,
        9 => 18.998,
        14 => 28.086,
        15 => 30.974,
        16 => 32.065,
        17 => 35.453,
        35 => 79.904,
        53 => 126.904,
        _ => z as f64 * 1.5, // rough fallback
    }
}

fn vdw_volume(z: u8) -> f64 {
    // Approximate van der Waals volumes in Å³
    match z {
        1 => 7.24,
        6 => 20.58,
        7 => 15.60,
        8 => 14.71,
        9 => 13.31,
        15 => 24.43,
        16 => 24.43,
        17 => 22.45,
        35 => 26.52,
        53 => 32.52,
        _ => 20.0,
    }
}

fn electronegativity(z: u8) -> f64 {
    match z {
        1 => 2.20,
        5 => 2.04,
        6 => 2.55,
        7 => 3.04,
        8 => 3.44,
        9 => 3.98,
        14 => 1.90,
        15 => 2.19,
        16 => 2.58,
        17 => 3.16,
        35 => 2.96,
        53 => 2.66,
        _ => 2.00,
    }
}

fn polarizability(z: u8) -> f64 {
    // Approximate atomic polarizabilities in Å³
    match z {
        1 => 0.387,
        6 => 1.026,
        7 => 0.802,
        8 => 0.637,
        9 => 0.440,
        14 => 3.707,
        15 => 3.222,
        16 => 2.900,
        17 => 2.315,
        35 => 3.013,
        53 => 5.415,
        _ => 1.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_whim_water() {
        let elements = vec![8, 1, 1];
        let positions = vec![
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469],
        ];

        let w = compute_whim(&elements, &positions, WhimWeighting::Unit);
        assert!(w.total_spread > 0.0);
        assert!(w.lambda[0] >= w.lambda[1]);
        assert!(w.lambda[1] >= w.lambda[2]);
        assert!((w.nu[0] + w.nu[1] + w.nu[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_whim_mass_weighted() {
        let elements = vec![8, 1, 1];
        let positions = vec![
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469],
        ];

        let w = compute_whim(&elements, &positions, WhimWeighting::Mass);
        assert!(w.total_spread > 0.0);
    }

    #[test]
    fn test_whim_empty() {
        let w = compute_whim(&[], &[], WhimWeighting::Unit);
        assert_eq!(w.total_spread, 0.0);
    }
}
