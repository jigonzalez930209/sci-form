//! Configuration Interaction Singles and Doubles (CISD).
//!
//! Extends CIS with doubly-excited determinants for improved correlation
//! energies and more accurate excited-state descriptions.
//!
//! $$H^{CISD}_{ia,jb} = H^{CIS}_{ia,jb}$$
//! $$H^{CISD}_{ijab,klcd} = (\varepsilon_a + \varepsilon_b - \varepsilon_i - \varepsilon_j)\delta_{ik}\delta_{jl}\delta_{ac}\delta_{bd}
//!   + \langle ab||cd \rangle \delta_{ij} + \langle ij||kl \rangle \delta_{ab} - P(ij)P(ab)(ia|jb)$$

use super::integrals::get_eri;
use nalgebra::{DMatrix, SymmetricEigen};
use serde::{Deserialize, Serialize};

/// CISD excitation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CisdExcitation {
    /// Excitation energy (Hartree).
    pub energy: f64,
    /// Excitation energy (eV).
    pub energy_ev: f64,
    /// Wavelength (nm).
    pub wavelength_nm: f64,
    /// Oscillator strength.
    pub oscillator_strength: f64,
    /// Dominant transition: (occ_indices, virt_indices).
    pub dominant_transition: Vec<(usize, usize)>,
    /// Character: "singles" or "doubles".
    pub character: String,
    /// Singles weight (fraction of singles contribution).
    pub singles_weight: f64,
}

/// CISD result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CisdResult {
    /// Excitations sorted by energy.
    pub excitations: Vec<CisdExcitation>,
    /// Ground-state correlation energy (Hartree).
    pub correlation_energy: f64,
    /// Total number of CSFs in the CI space.
    pub n_csfs: usize,
    /// Approximation used: "full" for exact CISD, "perturbative" for CIS+MP2.
    pub approximation: String,
}

const HARTREE_TO_EV: f64 = 27.21138602;

/// Compute CISD excitation energies from converged SCF data.
///
/// This builds the full CIS+D Hamiltonian including single and double
/// excitations. For large systems the matrix can be very large
/// (O(N²_occ × N²_virt)), so this is practical only for small molecules.
pub fn compute_cisd(
    orbital_energies: &[f64],
    coefficients: &DMatrix<f64>,
    eris: &[f64],
    n_basis: usize,
    n_occupied: usize,
    n_states: usize,
) -> CisdResult {
    let n_virtual = n_basis - n_occupied;
    let n_singles = n_occupied * n_virtual;

    // Count unique doubles: i<j, a<b
    let n_doubles = n_occupied * (n_occupied - 1) / 2 * n_virtual * (n_virtual - 1) / 2;

    let n_csfs = n_singles + n_doubles;

    if n_csfs == 0 || n_singles == 0 {
        return CisdResult {
            excitations: Vec::new(),
            correlation_energy: 0.0,
            n_csfs: 0,
            approximation: "full".to_string(),
        };
    }

    // For very large CI spaces, truncate to keep tractable
    let max_csfs = 5000;
    if n_csfs > max_csfs {
        // Fall back to CIS + MP2-style doubles correction
        return compute_cisd_perturbative(
            orbital_energies,
            coefficients,
            eris,
            n_basis,
            n_occupied,
            n_states,
        );
    }

    // Build full CISD Hamiltonian
    let mut h = DMatrix::zeros(n_csfs, n_csfs);

    // Index mappings for doubles
    let doubles_map = build_doubles_map(n_occupied, n_virtual);

    // Singles-Singles block (same as CIS)
    for ia in 0..n_singles {
        let i = ia / n_virtual;
        let a = ia % n_virtual + n_occupied;

        for jb in 0..=ia {
            let j = jb / n_virtual;
            let b = jb % n_virtual + n_occupied;

            let mut val = 0.0;
            if ia == jb {
                val += orbital_energies[a] - orbital_energies[i];
            }

            val += mo_eri_cisd(coefficients, eris, n_basis, i, a, j, b)
                - mo_eri_cisd(coefficients, eris, n_basis, i, j, a, b);

            h[(ia, jb)] = val;
            h[(jb, ia)] = val;
        }
    }

    // Doubles-Doubles block
    for (idx_d1, &(i, j, a, b)) in doubles_map.iter().enumerate() {
        let row = n_singles + idx_d1;

        // Diagonal
        h[(row, row)] = orbital_energies[a + n_occupied] + orbital_energies[b + n_occupied]
            - orbital_energies[i]
            - orbital_energies[j];

        // Off-diagonal doubles-doubles (simplified: direct + exchange)
        for (idx_d2, &(k, l, c, d)) in doubles_map.iter().enumerate() {
            if idx_d2 >= idx_d1 {
                break;
            }
            let col = n_singles + idx_d2;

            let val = if i == k && j == l {
                // Same occupied pair: (ab||cd)
                mo_eri_cisd(
                    coefficients,
                    eris,
                    n_basis,
                    a + n_occupied,
                    b + n_occupied,
                    c + n_occupied,
                    d + n_occupied,
                ) - mo_eri_cisd(
                    coefficients,
                    eris,
                    n_basis,
                    a + n_occupied,
                    d + n_occupied,
                    c + n_occupied,
                    b + n_occupied,
                )
            } else if a == c && b == d {
                // Same virtual pair: (ij||kl)
                mo_eri_cisd(coefficients, eris, n_basis, i, j, k, l)
                    - mo_eri_cisd(coefficients, eris, n_basis, i, l, k, j)
            } else {
                0.0
            };

            h[(row, col)] = val;
            h[(col, row)] = val;
        }
    }

    // Singles-Doubles coupling: <Φ_i^a|H|Φ_{ij}^{ab}>
    // Full coupling includes all permutations of occupied/virtual index matching.
    for ia in 0..n_singles {
        let i_s = ia / n_virtual;
        let a_s = ia % n_virtual + n_occupied;

        for (idx_d, &(i_d, j_d, a_d, b_d)) in doubles_map.iter().enumerate() {
            let col = n_singles + idx_d;
            let a_d_abs = a_d + n_occupied;
            let b_d_abs = b_d + n_occupied;

            // <S|H|D> coupling via Slater-Condon rules for one-electron difference
            // between single and double excitations.
            // The coupling is nonzero when the single shares one occupied and one
            // virtual index with the double. Four cases arise from antisymmetry:
            let mut val = 0.0;

            if i_s == i_d {
                // Shared occupied index i_s == i_d; differs in j_d and virtual pair
                val += mo_eri_cisd(coefficients, eris, n_basis, a_s, j_d, a_d_abs, b_d_abs)
                    - mo_eri_cisd(coefficients, eris, n_basis, a_s, b_d_abs, a_d_abs, j_d);
            }
            if i_s == j_d {
                // Shared occupied index i_s == j_d; antisymmetric permutation
                val -= mo_eri_cisd(coefficients, eris, n_basis, a_s, i_d, a_d_abs, b_d_abs)
                    - mo_eri_cisd(coefficients, eris, n_basis, a_s, b_d_abs, a_d_abs, i_d);
            }
            if a_s == a_d_abs {
                // Shared virtual index a_s == a_d; differs in occupied pair and b_d
                val += mo_eri_cisd(coefficients, eris, n_basis, i_s, j_d, i_d, b_d_abs)
                    - mo_eri_cisd(coefficients, eris, n_basis, i_s, b_d_abs, i_d, j_d);
            }
            if a_s == b_d_abs {
                // Shared virtual index a_s == b_d; antisymmetric permutation
                val -= mo_eri_cisd(coefficients, eris, n_basis, i_s, j_d, i_d, a_d_abs)
                    - mo_eri_cisd(coefficients, eris, n_basis, i_s, a_d_abs, i_d, j_d);
            }

            h[(ia, col)] = val;
            h[(col, ia)] = val;
        }
    }

    // Diagonalize
    let eigen = SymmetricEigen::new(h);

    // Sort eigenvalues
    let mut indices: Vec<usize> = (0..n_csfs).collect();
    indices.sort_by(|&a, &b| {
        eigen.eigenvalues[a]
            .partial_cmp(&eigen.eigenvalues[b])
            .unwrap()
    });

    // Ground state: lowest eigenvalue gives correlation energy
    let correlation_energy = if !indices.is_empty() {
        let ground_idx = indices[0];
        let e0 = eigen.eigenvalues[ground_idx];
        if e0 < 0.0 {
            e0
        } else {
            0.0
        }
    } else {
        0.0
    };

    // Collect excited states
    let n_out = n_states.min(n_csfs.saturating_sub(1));
    let mut excitations = Vec::with_capacity(n_out);

    let start = if correlation_energy < 0.0 { 1 } else { 0 };
    for &idx in indices.iter().skip(start).take(n_out) {
        let energy = eigen.eigenvalues[idx] - correlation_energy;
        if energy <= 0.0 {
            continue;
        }

        let energy_ev = energy * HARTREE_TO_EV;
        let wavelength_nm = if energy_ev > 0.0 {
            1239.84193 / energy_ev
        } else {
            f64::INFINITY
        };

        // Analyze character
        let col = eigen.eigenvectors.column(idx);
        let singles_weight: f64 = col.rows(0, n_singles).iter().map(|c| c * c).sum();

        // Find dominant transitions
        let mut transitions = Vec::new();
        let mut max_coeff = 0.0f64;
        let mut dom_idx = 0;
        for (k, &c) in col.iter().enumerate() {
            if c.abs() > max_coeff {
                max_coeff = c.abs();
                dom_idx = k;
            }
        }

        if dom_idx < n_singles {
            let i = dom_idx / n_virtual;
            let a = dom_idx % n_virtual + n_occupied;
            transitions.push((i, a));
        } else {
            let d_idx = dom_idx - n_singles;
            if d_idx < doubles_map.len() {
                let (i, j, a, b) = doubles_map[d_idx];
                transitions.push((i, a + n_occupied));
                transitions.push((j, b + n_occupied));
            }
        }

        let character = if singles_weight > 0.5 {
            "singles".to_string()
        } else {
            "doubles".to_string()
        };

        let f_osc = 2.0 / 3.0 * energy * singles_weight;

        excitations.push(CisdExcitation {
            energy,
            energy_ev,
            wavelength_nm,
            oscillator_strength: f_osc,
            dominant_transition: transitions,
            character,
            singles_weight,
        });
    }

    CisdResult {
        excitations,
        correlation_energy,
        n_csfs,
        approximation: "full".to_string(),
    }
}

/// Build index mapping for unique doubles: (i<j) → (a<b).
fn build_doubles_map(n_occ: usize, n_virt: usize) -> Vec<(usize, usize, usize, usize)> {
    let mut map = Vec::new();
    for i in 0..n_occ {
        for j in (i + 1)..n_occ {
            for a in 0..n_virt {
                for b in (a + 1)..n_virt {
                    map.push((i, j, a, b));
                }
            }
        }
    }
    map
}

/// Perturbative CISD approximation for large CI spaces.
///
/// Uses CIS for excitations + MP2-like correction for doubles contribution.
fn compute_cisd_perturbative(
    orbital_energies: &[f64],
    coefficients: &DMatrix<f64>,
    eris: &[f64],
    n_basis: usize,
    n_occupied: usize,
    n_states: usize,
) -> CisdResult {
    let n_virtual = n_basis - n_occupied;

    // Get CIS excitations
    let cis = super::cis::compute_cis(
        orbital_energies,
        coefficients,
        eris,
        n_basis,
        n_occupied,
        n_states,
    );

    // MP2 correlation energy from doubles
    // Denominator is always negative: ε_i + ε_j - ε_a - ε_b < 0 (occ < virt)
    let mut e_corr = 0.0;
    for i in 0..n_occupied {
        for j in (i + 1)..n_occupied {
            for a in n_occupied..n_basis {
                for b in (a + 1)..n_basis {
                    let ijab = mo_eri_cisd(coefficients, eris, n_basis, i, a, j, b);
                    let ijba = mo_eri_cisd(coefficients, eris, n_basis, i, b, j, a);
                    let antisym = ijab - ijba;
                    let denom = orbital_energies[i] + orbital_energies[j]
                        - orbital_energies[a]
                        - orbital_energies[b];
                    // Only include if denominator is negative (physical)
                    if denom < -1e-10 {
                        e_corr += antisym * antisym / denom;
                    }
                }
            }
        }
    }

    let excitations = cis
        .excitations
        .into_iter()
        .map(|e| CisdExcitation {
            energy: e.energy,
            energy_ev: e.energy_ev,
            wavelength_nm: e.wavelength_nm,
            oscillator_strength: e.oscillator_strength,
            dominant_transition: vec![e.dominant_transition],
            character: "singles".to_string(),
            singles_weight: 1.0,
        })
        .collect();

    CisdResult {
        excitations,
        correlation_energy: e_corr,
        n_csfs: n_occupied * n_virtual
            + n_occupied * (n_occupied - 1) / 2 * n_virtual * (n_virtual - 1) / 2,
        approximation: "perturbative".to_string(),
    }
}

/// MO-basis ERI.
fn mo_eri_cisd(
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cisd_h2() {
        // H2 minimal basis: 2 basis functions, 1 occupied, 1 virtual
        // 1 single, 0 doubles → should gracefully reduce to CIS
        let orbital_energies = vec![-0.6, 0.7];
        let c = DMatrix::from_row_slice(2, 2, &[0.7, 0.7, 0.7, -0.7]);
        let eris = vec![0.0; 16]; // Simplified
        let result = compute_cisd(&orbital_energies, &c, &eris, 2, 1, 5);
        assert_eq!(result.n_csfs, 1); // Only 1 single, 0 doubles
    }

    #[test]
    fn test_doubles_map() {
        let map = build_doubles_map(3, 4);
        // 3 choose 2 = 3, 4 choose 2 = 6 → 18 doubles
        assert_eq!(map.len(), 18);
        // Check ordering
        for &(i, j, a, b) in &map {
            assert!(i < j);
            assert!(a < b);
        }
    }
}
