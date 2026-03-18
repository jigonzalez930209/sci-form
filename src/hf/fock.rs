//! Fock matrix construction from one-electron and two-electron integrals.
//!
//! The Fock matrix is:
//! $$F_{\mu\nu} = H^{core}_{\mu\nu} + \sum_{\lambda\sigma} P_{\lambda\sigma}
//!   \left[(\mu\nu|\lambda\sigma) - \frac{1}{2}(\mu\lambda|\nu\sigma)\right]$$

use super::integrals::get_eri;
use nalgebra::DMatrix;

/// Build the Fock matrix from core Hamiltonian and density matrix.
///
/// `h_core`: one-electron integrals (T + V).
/// `density`: density matrix P.
/// `eris`: two-electron integrals (flat, 8-fold symmetry).
/// `n_basis`: number of basis functions.
pub fn build_fock(
    h_core: &DMatrix<f64>,
    density: &DMatrix<f64>,
    eris: &[f64],
    n_basis: usize,
) -> DMatrix<f64> {
    let mut fock = h_core.clone();

    for mu in 0..n_basis {
        for nu in 0..=mu {
            let mut g = 0.0;
            for lam in 0..n_basis {
                for sig in 0..n_basis {
                    let p = density[(lam, sig)];
                    if p.abs() < 1e-15 {
                        continue;
                    }
                    // Coulomb: (μν|λσ)
                    let j = get_eri(eris, mu, nu, lam, sig, n_basis);
                    // Exchange: (μλ|νσ)
                    let k = get_eri(eris, mu, lam, nu, sig, n_basis);
                    g += p * (j - 0.5 * k);
                }
            }
            fock[(mu, nu)] += g;
            if mu != nu {
                fock[(nu, mu)] = fock[(mu, nu)];
            }
        }
    }
    fock
}

/// Compute the electronic energy from density, h_core, and Fock matrix.
///
/// $$E_{elec} = \frac{1}{2} \sum_{\mu\nu} P_{\mu\nu} (H^{core}_{\mu\nu} + F_{\mu\nu})$$
pub fn electronic_energy(
    density: &DMatrix<f64>,
    h_core: &DMatrix<f64>,
    fock: &DMatrix<f64>,
) -> f64 {
    let n = density.nrows();
    let mut energy = 0.0;
    for mu in 0..n {
        for nu in 0..n {
            energy += density[(mu, nu)] * (h_core[(mu, nu)] + fock[(mu, nu)]);
        }
    }
    0.5 * energy
}

/// Compute nuclear repulsion energy.
///
/// $$E_{nuc} = \sum_{A>B} \frac{Z_A Z_B}{R_{AB}}$$
pub fn nuclear_repulsion(elements: &[u8], positions_bohr: &[[f64; 3]]) -> f64 {
    let n = elements.len();
    let mut e_nuc = 0.0;
    for a in 0..n {
        for b in (a + 1)..n {
            let dx = positions_bohr[a][0] - positions_bohr[b][0];
            let dy = positions_bohr[a][1] - positions_bohr[b][1];
            let dz = positions_bohr[a][2] - positions_bohr[b][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r > 1e-10 {
                e_nuc += elements[a] as f64 * elements[b] as f64 / r;
            }
        }
    }
    e_nuc
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_nuclear_repulsion_h2() {
        // H2 at 0.74 Å = 1.398 Bohr → V_nn = 1/1.398 ≈ 0.715
        let elements = [1u8, 1];
        let r = 0.74 * 1.8897259886; // Bohr
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, r]];
        let e_nuc = nuclear_repulsion(&elements, &positions);
        assert!(
            (e_nuc - 1.0 / r).abs() < 1e-10,
            "E_nuc = {e_nuc}, expected {}",
            1.0 / r
        );
    }

    #[test]
    fn test_electronic_energy_trace() {
        let n = 2;
        let p = DMatrix::from_row_slice(n, n, &[1.0, 0.0, 0.0, 1.0]);
        let h = DMatrix::from_row_slice(n, n, &[-1.0, 0.5, 0.5, -1.0]);
        let f = DMatrix::from_row_slice(n, n, &[-0.5, 0.3, 0.3, -0.5]);
        let e = electronic_energy(&p, &h, &f);
        // E = 0.5 * Tr(P(H+F)) = 0.5 * (P[0,0]*(H[0,0]+F[0,0]) + P[1,1]*(H[1,1]+F[1,1]))
        let expected = 0.5 * (1.0 * (-1.0 + -0.5) + 1.0 * (-1.0 + -0.5));
        assert!((e - expected).abs() < 1e-10);
    }
}
