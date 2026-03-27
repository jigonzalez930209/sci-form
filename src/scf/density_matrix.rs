//! Density matrix construction from MO coefficients.
//!
//! P_μν = 2 Σ_{k ∈ occ} C_μk C_νk
//!
//! The factor of 2 accounts for double occupancy in RHF (closed-shell).

use nalgebra::DMatrix;

/// Build the density matrix from MO coefficients.
///
/// P_μν = 2 Σ_{k=0}^{n_occ-1} C_μk · C_νk
pub fn build_density_matrix(coefficients: &DMatrix<f64>, n_occupied: usize) -> DMatrix<f64> {
    let n = coefficients.nrows();
    let mut p = DMatrix::zeros(n, n);

    for i in 0..n {
        for j in 0..=i {
            let mut p_ij = 0.0;
            for k in 0..n_occupied {
                p_ij += coefficients[(i, k)] * coefficients[(j, k)];
            }
            p_ij *= 2.0;
            p[(i, j)] = p_ij;
            p[(j, i)] = p_ij;
        }
    }

    p
}

/// Build the density matrix with fractional occupation numbers (FON).
///
/// Useful for HOMO-LUMO near-degeneracy where integer occupation causes
/// SCF oscillations. Applies a Fermi smearing around the HOMO-LUMO gap.
///
/// P_μν = Σ_k f_k · C_μk · C_νk  (f_k = Fermi occupation 0..2)
pub fn build_density_matrix_fon(
    coefficients: &DMatrix<f64>,
    orbital_energies: &[f64],
    n_electrons: usize,
    temperature_au: f64,
) -> DMatrix<f64> {
    let n = coefficients.nrows();
    let n_orb = orbital_energies.len();
    let _n_occ = n_electrons / 2;

    // Determine chemical potential (Fermi level) by bisection
    let mu = find_fermi_level(orbital_energies, n_electrons, temperature_au);

    // Compute fractional occupations
    let occupations: Vec<f64> = orbital_energies
        .iter()
        .map(|&e| 2.0 * fermi_dirac(e, mu, temperature_au))
        .collect();

    let mut p = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..=i {
            let mut p_ij = 0.0;
            for k in 0..n_orb.min(n) {
                p_ij += occupations[k] * coefficients[(i, k)] * coefficients[(j, k)];
            }
            p[(i, j)] = p_ij;
            p[(j, i)] = p_ij;
        }
    }
    p
}

fn fermi_dirac(energy: f64, mu: f64, temperature: f64) -> f64 {
    if temperature < 1e-15 {
        return if energy < mu {
            1.0
        } else if energy > mu {
            0.0
        } else {
            0.5
        };
    }
    let x = (energy - mu) / temperature;
    if x > 50.0 {
        0.0
    } else if x < -50.0 {
        1.0
    } else {
        1.0 / (1.0 + x.exp())
    }
}

fn find_fermi_level(orbital_energies: &[f64], n_electrons: usize, temperature: f64) -> f64 {
    let target = n_electrons as f64;
    let mut mu_lo = orbital_energies
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min)
        - 1.0;
    let mut mu_hi = orbital_energies
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max)
        + 1.0;

    for _ in 0..100 {
        let mu = 0.5 * (mu_lo + mu_hi);
        let n: f64 = orbital_energies
            .iter()
            .map(|&e| 2.0 * fermi_dirac(e, mu, temperature))
            .sum();
        if n < target {
            mu_lo = mu;
        } else {
            mu_hi = mu;
        }
        if (mu_hi - mu_lo).abs() < 1e-12 {
            break;
        }
    }
    0.5 * (mu_lo + mu_hi)
}

/// Compute the density matrix change (RMS difference).
pub fn density_rms_change(p_new: &DMatrix<f64>, p_old: &DMatrix<f64>) -> f64 {
    let n = p_new.nrows();
    let diff = p_new - p_old;
    let mut sum_sq = 0.0;
    for i in 0..n {
        for j in 0..n {
            sum_sq += diff[(i, j)] * diff[(i, j)];
        }
    }
    (sum_sq / (n * n) as f64).sqrt()
}

/// Compute the number of electrons from the density matrix.
///
/// N_e = Tr(PS) = Σ_μν P_μν S_μν
pub fn electron_count(p: &DMatrix<f64>, s: &DMatrix<f64>) -> f64 {
    let ps = p * s;
    let mut trace = 0.0;
    for i in 0..ps.nrows() {
        trace += ps[(i, i)];
    }
    trace
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_density_matrix_symmetric() {
        let c = DMatrix::from_row_slice(3, 3, &[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
        let p = build_density_matrix(&c, 1);

        for i in 0..3 {
            for j in 0..3 {
                assert!((p[(i, j)] - p[(j, i)]).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn test_density_trace() {
        let c = DMatrix::identity(3, 3);
        let p = build_density_matrix(&c, 2);
        let s = DMatrix::identity(3, 3);
        let n_e = electron_count(&p, &s);
        assert!((n_e - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_density_change() {
        let p1 = DMatrix::identity(2, 2);
        let mut p2 = DMatrix::identity(2, 2);
        p2[(0, 0)] = 1.1;
        let rms = density_rms_change(&p1, &p2);
        assert!(rms > 0.0);
    }

    #[test]
    fn test_fon_density_conserves_electrons() {
        let c = DMatrix::identity(4, 4);
        let energies = [-1.0, -0.5, 0.5, 1.0];
        let n_electrons = 4; // 2 occupied orbs
        let temp_au = 0.001; // very cold — should be close to integer occupation
        let p = build_density_matrix_fon(&c, &energies, n_electrons, temp_au);
        let s = DMatrix::identity(4, 4);
        let n_e = electron_count(&p, &s);
        assert!(
            (n_e - 4.0).abs() < 0.1,
            "FON should conserve ~4 electrons, got {n_e}"
        );
    }

    #[test]
    fn test_fon_at_zero_temp_matches_integer() {
        let c = DMatrix::identity(3, 3);
        let energies = [-1.0, -0.5, 0.5];
        let p_fon = build_density_matrix_fon(&c, &energies, 2, 1e-20);
        let p_int = build_density_matrix(&c, 1);
        // At zero temperature, FON should match integer occupation
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (p_fon[(i, j)] - p_int[(i, j)]).abs() < 1e-6,
                    "FON(T→0) should match integer density at ({i},{j})"
                );
            }
        }
    }

    #[test]
    fn test_fon_high_temp_smears_occupation() {
        let c = DMatrix::identity(3, 3);
        let energies = [-0.1, 0.0, 0.1]; // near-degenerate
        let p = build_density_matrix_fon(&c, &energies, 2, 1.0);
        // All orbitals should have partial occupation — no P_ii should be exactly 0 or 2
        for i in 0..3 {
            assert!(p[(i, i)] > 0.01, "at high T, all orbitals get partial occ");
            assert!(p[(i, i)] < 1.99, "at high T, no orbital is fully occupied");
        }
    }
}
