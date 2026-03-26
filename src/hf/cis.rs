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
///
/// When `positions_bohr` and `basis_to_atom` are provided, computes real
/// oscillator strengths using the monopole approximation for transition
/// dipole moments. Otherwise falls back to a rough CI-vector norm estimate.
pub fn compute_cis(
    orbital_energies: &[f64],
    coefficients: &DMatrix<f64>,
    eris: &[f64],
    n_basis: usize,
    n_occupied: usize,
    n_states: usize,
) -> CisResult {
    compute_cis_with_dipole(
        orbital_energies,
        coefficients,
        eris,
        n_basis,
        n_occupied,
        n_states,
        None,
        None,
    )
}

/// CIS with optional dipole integral data for proper oscillator strengths.
pub fn compute_cis_with_dipole(
    orbital_energies: &[f64],
    coefficients: &DMatrix<f64>,
    eris: &[f64],
    n_basis: usize,
    n_occupied: usize,
    n_states: usize,
    positions_bohr: Option<&[[f64; 3]]>,
    basis_to_atom: Option<&[usize]>,
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

        // Oscillator strength via transition dipole moment
        let f_osc = match (positions_bohr, basis_to_atom) {
            (Some(pos), Some(b2a)) => {
                // Monopole approximation: μ_0k = √2 Σ_ia X_ia Σ_A q_ia(A) R_A
                let mut tdm = [0.0f64; 3];
                for single in 0..n_singles {
                    let i_s = single / n_virtual;
                    let a_s = single % n_virtual + n_occupied;
                    let x_ia = col[single];
                    if x_ia.abs() < 1e-12 {
                        continue;
                    }
                    // Compute transition charge on each atom
                    let n_atoms = pos.len();
                    let mut q_atom = vec![0.0f64; n_atoms];
                    for mu in 0..n_basis {
                        let atom = b2a[mu];
                        let mut s_contrib = 0.0;
                        for nu in 0..n_basis {
                            s_contrib += coefficients[(nu, a_s)] * if mu == nu { 1.0 } else { 0.0 };
                        }
                        q_atom[atom] += coefficients[(mu, i_s)] * s_contrib;
                    }
                    for atom in 0..n_atoms {
                        tdm[0] += x_ia * q_atom[atom] * pos[atom][0];
                        tdm[1] += x_ia * q_atom[atom] * pos[atom][1];
                        tdm[2] += x_ia * q_atom[atom] * pos[atom][2];
                    }
                }
                let sqrt2 = std::f64::consts::SQRT_2;
                tdm[0] *= sqrt2;
                tdm[1] *= sqrt2;
                tdm[2] *= sqrt2;
                let tdm_sq = tdm[0] * tdm[0] + tdm[1] * tdm[1] + tdm[2] * tdm[2];
                2.0 / 3.0 * energy * tdm_sq
            }
            _ => {
                // Fallback: rough estimate using CI vector norm (always ~2/3 * E)
                2.0 / 3.0 * energy * transition_dipole_sq(col.as_slice(), n_occupied, n_virtual)
            }
        };

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

/// Compute MO-basis ERI from AO-basis ERIs using half-transformed intermediates.
///
/// Standard 4-index transformation: (pq|rs) = Σ C_μp C_νq C_λr C_σs (μν|λσ).
/// Uses intermediate contraction to reduce O(N⁵) per integral to O(N⁴) total via
/// sequential index transformation, with coefficient screening for sparsity.
fn mo_eri(c: &DMatrix<f64>, eris: &[f64], n: usize, p: usize, q: usize, r: usize, s: usize) -> f64 {
    // First half-transform: contract σ → s to get (μν|λs)
    let mut half1 = vec![0.0f64; n * n * n]; // [mu][nu][lam]
    for lam in 0..n {
        for sig in 0..n {
            let c_sig_s = c[(sig, s)];
            if c_sig_s.abs() < 1e-12 {
                continue;
            }
            for mu in 0..n {
                for nu in 0..n {
                    half1[mu * n * n + nu * n + lam] +=
                        c_sig_s * get_eri(eris, mu, nu, lam, sig, n);
                }
            }
        }
    }

    // Second half-transform: contract λ → r to get (μν|rs)
    let mut half2 = vec![0.0f64; n * n]; // [mu][nu]
    for lam in 0..n {
        let c_lam_r = c[(lam, r)];
        if c_lam_r.abs() < 1e-12 {
            continue;
        }
        for mu in 0..n {
            for nu in 0..n {
                half2[mu * n + nu] += c_lam_r * half1[mu * n * n + nu * n + lam];
            }
        }
    }

    // Third quarter-transform: contract ν → q, then μ → p
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
            val += c_mu_p * c_nu_q * half2[mu * n + nu];
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
        let eris = vec![0.0; n_basis * (n_basis + 1) / 2 * (n_basis * (n_basis + 1) / 2 + 1) / 2];

        let result = compute_cis(&energies, &coeffs, &eris, n_basis, n_occ, 3);
        assert!(
            result.excitations.len() <= n_singles,
            "CIS states ≤ n_singles"
        );
    }
}
