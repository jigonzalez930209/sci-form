//! Self-Consistent Field (SCF) solver with DIIS acceleration.
//!
//! Solves the Roothaan-Hall equations iteratively:
//! $$\mathbf{FC} = \mathbf{SC}\boldsymbol{\varepsilon}$$
//!
//! Uses Löwdin orthogonalization ($S^{-1/2}$) to transform the generalized
//! eigenvalue problem into a standard one, then applies DIIS for convergence.

use super::fock::{build_fock, electronic_energy};
use nalgebra::{DMatrix, DVector, SymmetricEigen};

/// SCF convergence result.
#[derive(Debug, Clone)]
pub struct ScfResult {
    /// Converged electronic energy (Hartree).
    pub energy: f64,
    /// Orbital energies (eigenvalues).
    pub orbital_energies: Vec<f64>,
    /// MO coefficient matrix C (columns are MOs).
    pub coefficients: DMatrix<f64>,
    /// Converged density matrix P.
    pub density: DMatrix<f64>,
    /// Number of SCF iterations.
    pub iterations: usize,
    /// Whether SCF converged.
    pub converged: bool,
}

/// SCF solver configuration.
pub struct ScfConfig {
    pub max_iter: usize,
    pub energy_threshold: f64,
    pub density_threshold: f64,
    pub diis_size: usize,
}

impl Default for ScfConfig {
    fn default() -> Self {
        ScfConfig {
            max_iter: 100,
            energy_threshold: 1e-8,
            density_threshold: 1e-6,
            diis_size: 6,
        }
    }
}

/// Run the SCF procedure.
pub fn solve_scf(
    h_core: &DMatrix<f64>,
    s_mat: &DMatrix<f64>,
    eris: &[f64],
    n_electrons: usize,
    config: &ScfConfig,
) -> ScfResult {
    let n = h_core.nrows();
    let n_occ = n_electrons / 2;

    // Löwdin orthogonalization: S^{-1/2}
    let s_half_inv = lowdin_orthogonalization(s_mat);

    // Initial guess: diagonalize H_core
    let (mut energies, mut coeffs) = diagonalize_fock(h_core, &s_half_inv);
    let mut density = build_density(&coeffs, n_occ);

    let mut prev_energy = 0.0;
    let mut converged = false;
    let mut iterations = 0;

    // DIIS storage
    let mut diis_focks: Vec<DMatrix<f64>> = Vec::new();
    let mut diis_errors: Vec<DMatrix<f64>> = Vec::new();

    for iter in 0..config.max_iter {
        iterations = iter + 1;

        let fock = build_fock(h_core, &density, eris, n);
        let energy = electronic_energy(&density, h_core, &fock);

        // DIIS error: FPS - SPF
        let error = &fock * &density * s_mat - s_mat * &density * &fock;
        let error_norm = error.iter().map(|v| v * v).sum::<f64>().sqrt();

        // DIIS extrapolation
        diis_focks.push(fock.clone());
        diis_errors.push(error);
        if diis_focks.len() > config.diis_size {
            diis_focks.remove(0);
            diis_errors.remove(0);
        }

        let fock_diis = if diis_focks.len() >= 2 {
            diis_extrapolate(&diis_focks, &diis_errors)
        } else {
            fock
        };

        let (new_energies, new_coeffs) = diagonalize_fock(&fock_diis, &s_half_inv);
        let new_density = build_density(&new_coeffs, n_occ);

        let de = (energy - prev_energy).abs();

        if de < config.energy_threshold && error_norm < config.density_threshold {
            converged = true;
            energies = new_energies;
            coeffs = new_coeffs;
            density = new_density;
            break;
        }

        prev_energy = energy;
        energies = new_energies;
        coeffs = new_coeffs;
        density = new_density;
    }

    let final_energy = electronic_energy(
        &density,
        h_core,
        &build_fock(h_core, &density, eris, n),
    );

    ScfResult {
        energy: final_energy,
        orbital_energies: energies.as_slice().to_vec(),
        coefficients: coeffs,
        density,
        iterations,
        converged,
    }
}

fn lowdin_orthogonalization(s: &DMatrix<f64>) -> DMatrix<f64> {
    let eigen = SymmetricEigen::new(s.clone());
    let n = s.nrows();
    let mut s_inv_half = DMatrix::zeros(n, n);

    for i in 0..n {
        let val = eigen.eigenvalues[i];
        if val > 1e-10 {
            let factor = 1.0 / val.sqrt();
            let col = eigen.eigenvectors.column(i);
            s_inv_half += factor * &col * col.transpose();
        }
    }
    s_inv_half
}

fn diagonalize_fock(
    fock: &DMatrix<f64>,
    s_half_inv: &DMatrix<f64>,
) -> (DVector<f64>, DMatrix<f64>) {
    let f_prime = s_half_inv.transpose() * fock * s_half_inv;
    let eigen = SymmetricEigen::new(f_prime);

    // Sort eigenvalues/vectors by energy
    let n = eigen.eigenvalues.len();
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| {
        eigen.eigenvalues[a]
            .partial_cmp(&eigen.eigenvalues[b])
            .unwrap()
    });

    let sorted_energies = DVector::from_fn(n, |i, _| eigen.eigenvalues[indices[i]]);
    let sorted_vecs = DMatrix::from_fn(n, n, |r, c| eigen.eigenvectors[(r, indices[c])]);

    // Back-transform: C = S^{-1/2} C'
    let coeffs = s_half_inv * sorted_vecs;
    (sorted_energies, coeffs)
}

fn build_density(coeffs: &DMatrix<f64>, n_occ: usize) -> DMatrix<f64> {
    let n = coeffs.nrows();
    let mut density = DMatrix::zeros(n, n);
    for i in 0..n_occ {
        let col = coeffs.column(i);
        density += 2.0 * &col * col.transpose();
    }
    density
}

fn diis_extrapolate(focks: &[DMatrix<f64>], errors: &[DMatrix<f64>]) -> DMatrix<f64> {
    let m = errors.len();
    let mut b = DMatrix::zeros(m + 1, m + 1);

    for i in 0..m {
        for j in 0..=i {
            let bij: f64 = errors[i].iter().zip(errors[j].iter()).map(|(a, b)| a * b).sum();
            b[(i, j)] = bij;
            b[(j, i)] = bij;
        }
    }
    for i in 0..m {
        b[(m, i)] = -1.0;
        b[(i, m)] = -1.0;
    }

    let mut rhs = DVector::zeros(m + 1);
    rhs[m] = -1.0;

    // Solve B·c = rhs using pseudo-inverse
    let svd = b.svd(true, true);
    let c = match svd.solve(&rhs, 1e-10) {
        Ok(c) => c,
        Err(_) => {
            // Fallback: use latest Fock
            return focks.last().unwrap().clone();
        }
    };

    let mut f_diis = DMatrix::zeros(focks[0].nrows(), focks[0].ncols());
    for i in 0..m {
        f_diis += c[i] * &focks[i];
    }
    f_diis
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lowdin_identity() {
        let s = DMatrix::identity(3, 3);
        let s_inv = lowdin_orthogonalization(&s);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (s_inv[(i, j)] - expected).abs() < 1e-10,
                    "S^{{-1/2}} of identity should be identity"
                );
            }
        }
    }
}
