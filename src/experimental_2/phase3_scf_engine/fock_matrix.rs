//! Fock matrix F = H⁰ + G(P) construction.
//!
//! The Fock matrix in the Hartree-Fock method:
//!
//!   F_μν = H⁰_μν + G_μν
//!   G_μν = Σ_{λσ} P_{λσ} [(μν|λσ) - 0.5·(μλ|νσ)]
//!
//! where (μν|λσ) are two-electron repulsion integrals.
//! The first term is Coulomb (J) and the second is exchange (K).

use nalgebra::DMatrix;

use super::super::phase2_quantum_engine::two_electron::TwoElectronIntegrals;

/// Build the Fock matrix from core Hamiltonian, density, and ERIs.
///
/// F_μν = H_core_μν + Σ_{λσ} P_{λσ} [(μν|λσ) - 0.5(μλ|νσ)]
pub fn build_fock_matrix(
    h_core: &DMatrix<f64>,
    density: &DMatrix<f64>,
    eris: &TwoElectronIntegrals,
) -> DMatrix<f64> {
    let n = h_core.nrows();
    let mut f = h_core.clone();

    for mu in 0..n {
        for nu in 0..=mu {
            let mut g_mn = 0.0;

            for lam in 0..n {
                for sig in 0..n {
                    let p_ls = density[(lam, sig)];
                    // Coulomb: (μν|λσ)
                    let j = eris.get(mu, nu, lam, sig);
                    // Exchange: -0.5·(μλ|νσ)
                    let k = eris.get(mu, lam, nu, sig);
                    g_mn += p_ls * (j - 0.5 * k);
                }
            }

            f[(mu, nu)] += g_mn;
            if mu != nu {
                f[(nu, mu)] += g_mn;
            }
        }
    }

    f
}

/// Build the Fock matrix for DFTB/xTB-style Hamiltonians.
///
/// In SCC-DFTB, the Fock matrix includes charge-dependent corrections:
///
///   F_μν = H⁰_μν + 0.5·S_μν · Σ_C (γ_AC + γ_BC)·Δq_C
///
/// where A,B host μ,ν and γ_AB is the Coulomb interaction function.
pub fn build_fock_dftb(
    h_core: &DMatrix<f64>,
    overlap: &DMatrix<f64>,
    delta_charges: &[f64],
    gamma_matrix: &DMatrix<f64>,
    basis_to_atom: &[usize],
) -> DMatrix<f64> {
    let n = h_core.nrows();
    let mut f = h_core.clone();

    for mu in 0..n {
        let atom_a = basis_to_atom[mu];
        for nu in 0..n {
            let atom_b = basis_to_atom[nu];

            // Charge shift: 0.5 · S_μν · Σ_C (γ_AC + γ_BC) · Δq_C
            let mut shift = 0.0;
            for (c, &dq) in delta_charges.iter().enumerate() {
                shift += (gamma_matrix[(atom_a, c)] + gamma_matrix[(atom_b, c)]) * dq;
            }
            f[(mu, nu)] += 0.5 * overlap[(mu, nu)] * shift;
        }
    }

    f
}

/// Compute the Coulomb interaction matrix γ_AB for SCC-DFTB.
///
/// γ_AB = 1 / sqrt((1/η_A + 1/η_B)² / 4 + R_AB²)
///
/// where η is the chemical hardness (Hubbard parameter) and
/// R_AB is the interatomic distance.
pub fn build_gamma_matrix(
    eta: &[f64],
    positions_bohr: &[[f64; 3]],
) -> DMatrix<f64> {
    let n = eta.len();
    let mut gamma = DMatrix::zeros(n, n);

    for a in 0..n {
        // Diagonal: γ_AA = η_A (on-site)
        gamma[(a, a)] = eta[a];

        for b in (a + 1)..n {
            let dx = positions_bohr[a][0] - positions_bohr[b][0];
            let dy = positions_bohr[a][1] - positions_bohr[b][1];
            let dz = positions_bohr[a][2] - positions_bohr[b][2];
            let r_sq = dx * dx + dy * dy + dz * dz;

            let avg_eta_inv = 0.5 * (1.0 / eta[a] + 1.0 / eta[b]);
            let gamma_ab = 1.0 / (avg_eta_inv * avg_eta_inv + r_sq).sqrt();

            gamma[(a, b)] = gamma_ab;
            gamma[(b, a)] = gamma_ab;
        }
    }

    gamma
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fock_equals_h_core_when_no_electrons() {
        let h = DMatrix::from_row_slice(2, 2, &[-1.0, -0.5, -0.5, -1.0]);
        let p = DMatrix::zeros(2, 2);

        // Need ERI storage
        use crate::experimental_2::phase2_quantum_engine::basis_set::BasisSet;
        use crate::experimental_2::phase2_quantum_engine::two_electron::TwoElectronIntegrals;

        let basis = BasisSet::sto3g(&[1, 1], &[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]);
        let eris = TwoElectronIntegrals::compute(&basis);

        let f = build_fock_matrix(&h, &p, &eris);
        // With zero density, Fock should equal core Hamiltonian
        for i in 0..2 {
            for j in 0..2 {
                assert!((f[(i, j)] - h[(i, j)]).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn test_gamma_matrix_symmetric() {
        let eta = vec![0.5, 0.4, 0.3];
        let pos = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let gamma = build_gamma_matrix(&eta, &pos);

        for i in 0..3 {
            for j in 0..3 {
                assert!((gamma[(i, j)] - gamma[(j, i)]).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn test_gamma_diagonal_equals_eta() {
        let eta = vec![0.5, 0.4];
        let pos = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let gamma = build_gamma_matrix(&eta, &pos);

        assert!((gamma[(0, 0)] - 0.5).abs() < 1e-14);
        assert!((gamma[(1, 1)] - 0.4).abs() < 1e-14);
    }
}
