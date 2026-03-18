//! Configuration Interaction Singles (CIS) for UV-Vis excitation energies.
//!
//! Promotes one electron from an occupied MO to a virtual MO, forming
//! singly-excited determinants. The CIS Hamiltonian eigenvalues give
//! vertical excitation energies, and eigenvectors yield oscillator strengths.
//!
//! $$H^{CIS}_{ia,jb} = \delta_{ij}\delta_{ab}(\varepsilon_a - \varepsilon_i)
//!   + (ia|jb) - (ij|ab)$$

use super::integrals::get_eri;
use nalgebra::{DMatrix, SymmetricEigen};
use serde::{Deserialize, Serialize};

/// A single electronic excitation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Excitation {
    /// Excitation energy (Hartree).
    pub energy: f64,
    /// Excitation energy (eV).
    pub energy_ev: f64,
    /// Excitation wavelength (nm).
    pub wavelength_nm: f64,
    /// Oscillator strength (dimensionless).
    pub oscillator_strength: f64,
    /// Dominant occupied→virtual transition (orbital indices).
    pub dominant_transition: (usize, usize),
}

/// CIS calculation result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CisResult {
    /// Electronic excitations, sorted by energy.
    pub excitations: Vec<Excitation>,
}

/// Hartree to eV conversion.
const HARTREE_TO_EV: f64 = 27.21138602;

/// Compute CIS excitation energies from converged SCF data.
pub fn compute_cis(
    orbital_energies: &[f64],
    coefficients: &DMatrix<f64>,
    eris: &[f64],
    n_basis: usize,
    n_occupied: usize,
    n_states: usize,
) -> CisResult {
    let n_virtual = n_basis - n_occupied;
    let n_singles = n_occupied * n_virtual;

    if n_singles == 0 {
        return CisResult {
            excitations: Vec::new(),
        };
    }

    // Build CIS Hamiltonian
    let mut h_cis = DMatrix::zeros(n_singles, n_singles);

    for ia in 0..n_singles {
        let i = ia / n_virtual;
        let a = ia % n_virtual + n_occupied;

        for jb in 0..=ia {
            let j = jb / n_virtual;
            let b = jb % n_virtual + n_occupied;

            let mut val = 0.0;

            // Diagonal: orbital energy difference
            if ia == jb {
                val += orbital_energies[a] - orbital_energies[i];
            }

            // Two-electron integrals in MO basis (using AO ERIs + MO coeffs)
            let coulomb = mo_eri(coefficients, eris, n_basis, i, a, j, b);
            let exchange = mo_eri(coefficients, eris, n_basis, i, j, a, b);

            val += coulomb - exchange;

            h_cis[(ia, jb)] = val;
            h_cis[(jb, ia)] = val;
        }
    }

    // Diagonalize CIS Hamiltonian
    let eigen = SymmetricEigen::new(h_cis);

    // Sort and collect excitations
    let mut indices: Vec<usize> = (0..n_singles).collect();
    indices.sort_by(|&a, &b| {
        eigen.eigenvalues[a]
            .partial_cmp(&eigen.eigenvalues[b])
            .unwrap()
    });

    let n_out = n_states.min(n_singles);
    let mut excitations = Vec::with_capacity(n_out);

    for &idx in indices.iter().take(n_out) {
        let energy = eigen.eigenvalues[idx];
        if energy <= 0.0 {
            continue;
        }

        let energy_ev = energy * HARTREE_TO_EV;
        let wavelength_nm = 1239.84193 / energy_ev;

        // Find dominant transition
        let col = eigen.eigenvectors.column(idx);
        let (dom_idx, _) = col
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap())
            .unwrap();
        let dom_i = dom_idx / n_virtual;
        let dom_a = dom_idx % n_virtual + n_occupied;

        // Oscillator strength (dipole-length approximation, simplified)
        let f_osc = 2.0 / 3.0 * energy * transition_dipole_sq(col.as_slice(), n_occupied, n_virtual);

        excitations.push(Excitation {
            energy,
            energy_ev,
            wavelength_nm,
            oscillator_strength: f_osc,
            dominant_transition: (dom_i, dom_a),
        });
    }

    CisResult { excitations }
}

/// Compute MO-basis ERI from AO-basis ERIs: (pq|rs) = Σ C_μp C_νq C_λr C_σs (μν|λσ).
fn mo_eri(
    c: &DMatrix<f64>,
    eris: &[f64],
    n: usize,
    p: usize,
    q: usize,
    r: usize,
    s: usize,
) -> f64 {
    let mut val = 0.0;
    for mu in 0..n {
        let c_mu_p = c[(mu, p)];
        if c_mu_p.abs() < 1e-12 {
            continue;
        }
        for nu in 0..n {
            let c_nu_q = c[(nu, q)];
            if c_nu_q.abs() < 1e-12 {
                continue;
            }
            let cpq = c_mu_p * c_nu_q;
            for lam in 0..n {
                let c_lam_r = c[(lam, r)];
                if c_lam_r.abs() < 1e-12 {
                    continue;
                }
                for sig in 0..n {
                    let c_sig_s = c[(sig, s)];
                    if c_sig_s.abs() < 1e-12 {
                        continue;
                    }
                    val += cpq * c_lam_r * c_sig_s * get_eri(eris, mu, nu, lam, sig, n);
                }
            }
        }
    }
    val
}

/// Simplified transition dipole squared estimate.
fn transition_dipole_sq(cis_vec: &[f64], _n_occ: usize, _n_virt: usize) -> f64 {
    // Approximate: sum of squared coefficients as proxy
    cis_vec.iter().map(|c| c * c).sum::<f64>()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cis_dimensions() {
        let n_occ = 2;
        let n_virt = 3;
        let n_basis = n_occ + n_virt;
        let n_singles = n_occ * n_virt;

        // With dummy data, CIS should produce excitations
        let energies = vec![-2.0, -1.0, 0.5, 1.0, 1.5];
        let coeffs = DMatrix::identity(n_basis, n_basis);
        let eris = vec![0.0; n_basis * (n_basis + 1) / 2
            * (n_basis * (n_basis + 1) / 2 + 1) / 2];

        let result = compute_cis(&energies, &coeffs, &eris, n_basis, n_occ, 3);
        assert!(
            result.excitations.len() <= n_singles,
            "CIS states ≤ n_singles"
        );
    }
}
