//! Fock matrix construction.
//!
//! F_μν = H⁰_μν + G_μν
//!
//! where G is the two-electron contribution:
//!   G_μν = Σ_{λσ} P_{λσ} [(μν|λσ) - ½(μλ|νσ)]
//!
//! Also provides DFTB-type Fock matrix construction via SCC γ-matrix.

use nalgebra::DMatrix;

use super::two_electron::TwoElectronIntegrals;

/// Build the Fock matrix from core Hamiltonian, density matrix,
/// and two-electron integrals (Hartree-Fock).
///
/// F_μν = H⁰_μν + Σ_{λσ} P_{λσ} [(μν|λσ) - ½(μλ|νσ)]
pub fn build_fock_matrix(
    h_core: &DMatrix<f64>,
    density: &DMatrix<f64>,
    eris: &TwoElectronIntegrals,
) -> DMatrix<f64> {
    let n = h_core.nrows();
    let mut fock = h_core.clone();

    for mu in 0..n {
        for nu in 0..=mu {
            let mut g_mn = 0.0;
            for lam in 0..n {
                for sig in 0..n {
                    let p_ls = density[(lam, sig)];
                    // Coulomb: (μν|λσ)
                    let j = eris.get(mu, nu, lam, sig);
                    // Exchange: (μλ|νσ)
                    let k = eris.get(mu, lam, nu, sig);
                    g_mn += p_ls * (j - 0.5 * k);
                }
            }
            fock[(mu, nu)] += g_mn;
            if mu != nu {
                fock[(nu, mu)] += g_mn;
            }
        }
    }

    fock
}

/// Build the Fock matrix for SCC-DFTB.
///
/// F_μν = H⁰_μν + ½ S_μν Σ_C (Δq_A + Δq_C) γ_{AC}
///
/// where A is the atom of basis function μ, and Δq are Mulliken charge
/// differences from the reference.
pub fn build_fock_dftb(
    h_core: &DMatrix<f64>,
    overlap: &DMatrix<f64>,
    delta_charges: &[f64],
    gamma_matrix: &DMatrix<f64>,
    basis_to_atom: &[usize],
) -> DMatrix<f64> {
    let n = h_core.nrows();
    let mut fock = h_core.clone();

    for mu in 0..n {
        let atom_a = basis_to_atom[mu];
        for nu in 0..=mu {
            let atom_b = basis_to_atom[nu];

            let mut shift = 0.0;
            let n_atoms = delta_charges.len();
            for c in 0..n_atoms {
                shift += (delta_charges[atom_a] + delta_charges[c]) * gamma_matrix[(atom_a, c)];
                shift += (delta_charges[atom_b] + delta_charges[c]) * gamma_matrix[(atom_b, c)];
            }
            shift *= 0.25; // ½ × ½ from symmetric sum

            let correction = overlap[(mu, nu)] * shift;
            fock[(mu, nu)] += correction;
            if mu != nu {
                fock[(nu, mu)] += correction;
            }
        }
    }

    fock
}

/// Build the γ-matrix for SCC-DFTB Coulomb interactions.
///
/// γ_{AB} = 1 / |R_A - R_B| for A ≠ B
/// γ_{AA} = η_A (chemical hardness / Hubbard parameter)
pub fn build_gamma_matrix(eta: &[f64], positions_bohr: &[[f64; 3]]) -> DMatrix<f64> {
    let n = eta.len();
    let mut gamma = DMatrix::zeros(n, n);

    for a in 0..n {
        gamma[(a, a)] = eta[a];
        for b in (a + 1)..n {
            let dx = positions_bohr[a][0] - positions_bohr[b][0];
            let dy = positions_bohr[a][1] - positions_bohr[b][1];
            let dz = positions_bohr[a][2] - positions_bohr[b][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            let g = 1.0 / r;
            gamma[(a, b)] = g;
            gamma[(b, a)] = g;
        }
    }

    gamma
}

#[cfg(test)]
mod tests {
    use super::super::basis::BasisSet;
    use super::*;

    #[test]
    fn test_fock_symmetric() {
        let basis = BasisSet::sto3g(&[1, 1], &[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]);
        let n = basis.n_basis;
        let h = DMatrix::from_fn(n, n, |i, j| if i == j { -1.0 } else { -0.3 });
        let p = DMatrix::from_fn(n, n, |i, j| if i == j { 1.0 } else { 0.2 });
        let eris = TwoElectronIntegrals::compute(&basis);

        let f = build_fock_matrix(&h, &p, &eris);

        for i in 0..n {
            for j in 0..n {
                assert!(
                    (f[(i, j)] - f[(j, i)]).abs() < 1e-12,
                    "Fock not symmetric at ({}, {})",
                    i,
                    j
                );
            }
        }
    }

    #[test]
    fn test_gamma_matrix_symmetric() {
        let eta = vec![0.5, 0.5];
        let pos = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let gamma = build_gamma_matrix(&eta, &pos);

        assert!((gamma[(0, 0)] - 0.5).abs() < 1e-14);
        assert!((gamma[(1, 1)] - 0.5).abs() < 1e-14);
        assert!((gamma[(0, 1)] - 0.5).abs() < 1e-14); // 1/2.0
        assert!((gamma[(0, 1)] - gamma[(1, 0)]).abs() < 1e-14);
    }

    #[test]
    fn test_fock_dftb_basic() {
        let n = 2;
        let h = DMatrix::from_fn(n, n, |i, j| if i == j { -1.0 } else { -0.2 });
        let s = DMatrix::from_fn(n, n, |i, j| if i == j { 1.0 } else { 0.5 });
        let delta_q = vec![0.1, -0.1];
        let gamma = DMatrix::from_fn(2, 2, |i, j| if i == j { 0.5 } else { 0.3 });
        let basis_to_atom = vec![0, 1];

        let f = build_fock_dftb(&h, &s, &delta_q, &gamma, &basis_to_atom);

        // Should be symmetric
        assert!((f[(0, 1)] - f[(1, 0)]).abs() < 1e-12);
    }
}
