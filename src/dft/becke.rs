//! Becke partitioning scheme for multi-center molecular integration.
//!
//! Implements the Becke (1988) atomic partitioning with Bragg-Slater radius adjustments.

/// Compute the Becke partition weight for a grid point assigned to atom `atom_idx`.
///
/// Uses the Becke-Stratmann scheme with smoothed step functions (k=3 iterations).
pub fn becke_weights(
    point: &[f64; 3],
    atom_idx: usize,
    positions: &[[f64; 3]],
    atomic_numbers: &[u8],
) -> f64 {
    let n_atoms = positions.len();
    if n_atoms == 1 {
        return 1.0;
    }

    // Compute raw Becke partition values for each atom
    let mut p = vec![1.0; n_atoms];

    for i in 0..n_atoms {
        for j in (i + 1)..n_atoms {
            let r_ij = distance(&positions[i], &positions[j]);
            if r_ij < 1e-12 {
                continue;
            }

            let r_ip = distance(&positions[i], point);
            let r_jp = distance(&positions[j], point);

            // Elliptical coordinate: μ = (r_ip - r_jp) / r_ij
            let mu = (r_ip - r_jp) / r_ij;

            // Adjust for heteronuclear pairs using Bragg-Slater radii
            let chi = bragg_slater_ratio(atomic_numbers[i], atomic_numbers[j]);
            let u = (chi - 1.0) / (chi + 1.0);
            let a = (u / (u * u - 1.0)).clamp(-0.5, 0.5);
            let nu = mu + a * (1.0 - mu * mu);

            // Apply smoothed step function (3 iterations of Becke's function)
            let s = step_function_3(nu);

            p[i] *= s;
            p[j] *= 1.0 - s;
        }
    }

    let total: f64 = p.iter().sum();
    if total.abs() < 1e-30 {
        return 1.0 / n_atoms as f64;
    }
    p[atom_idx] / total
}

/// Becke's smoothed step function: s(μ) = 1/2 (1 - p₃(μ))
/// where p₃ is the 3-iteration polynomial.
fn step_function_3(mu: f64) -> f64 {
    let mu = mu.clamp(-1.0, 1.0);
    let mut f = mu;
    // Three iterations of p(x) = 3/2 x - 1/2 x³
    for _ in 0..3 {
        f = 1.5 * f - 0.5 * f * f * f;
    }
    0.5 * (1.0 - f)
}

/// Bragg-Slater radius ratio χ = R_i / R_j.
fn bragg_slater_ratio(z_i: u8, z_j: u8) -> f64 {
    let r_i = bragg_radius(z_i);
    let r_j = bragg_radius(z_j);
    r_i / r_j
}

/// Bragg-Slater atomic radii in Bohr.
fn bragg_radius(z: u8) -> f64 {
    let ang_to_bohr = 1.889_725_988_6;
    let r_ang = match z {
        1 => 0.25,
        2 => 0.31,
        6 => 0.70,
        7 => 0.65,
        8 => 0.60,
        9 => 0.50,
        14 => 1.10,
        15 => 1.00,
        16 => 1.00,
        17 => 1.00,
        _ => 1.50,
    };
    r_ang * ang_to_bohr
}

fn distance(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_atom_weight_is_one() {
        let pos = [[0.0, 0.0, 0.0]];
        let z = [6u8];
        let w = becke_weights(&[1.0, 0.0, 0.0], 0, &pos, &z);
        assert!((w - 1.0).abs() < 1e-10);
    }

    #[test]
    fn two_atoms_midpoint() {
        let pos = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let z = [6u8, 6];
        // At midpoint, weights should be ~0.5 for homonuclear pair
        let w0 = becke_weights(&[2.0, 0.0, 0.0], 0, &pos, &z);
        let w1 = becke_weights(&[2.0, 0.0, 0.0], 1, &pos, &z);
        assert!((w0 + w1 - 1.0).abs() < 1e-10);
        assert!((w0 - 0.5).abs() < 0.1);
    }
}
