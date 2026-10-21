//! Configuration Interaction Singles and Doubles (CISD).
//!
//! Determinant-based (spin-orbital) CISD over the RHF reference: the CI space
//! contains the reference determinant, all singly excited determinants and all
//! doubly excited determinants, including the *same-orbital* doubles
//! (`(iα iβ) → (aα aβ)`, the dominant correlation determinant of minimal
//! basis sets such as H₂/STO-3G) and mixed-spin doubles.
//!
//! Matrix elements use the standard Slater–Condon rules in the canonical
//! determinant basis:
//!
//! ```text
//!   <D|H|D>          = Σ_i h_ii + ½ Σ_ij [(ii|jj) − (ij|ji)]
//!   <D|H|D_i^a>      = σ · [h_ia + Σ_j ((ia|jj) − (ij|ja))]
//!   <D|H|D_{ij}^{ab}> = σ · [(ia|jb) − (ib|ja)]
//! ```
//!
//! with `σ` the Slater–Condon phase: the excitation operator is taken in the
//! standard form `a†_a a†_b a_j a_i` (removals ascending, creations
//! *descending* — the ordering for which the value formula above pairs the
//! ascending removed list with the ascending added list). Getting this
//! ordering wrong is invisible for the H₂O ground state but corrupts the
//! larger spaces (CH₂O/C₂H₄ off by ~1e-2 Ha), see `tests/cisd_brute.rs`.
//! One-electron integrals are recovered from the canonical
//! Fock condition `h_pq = ε_p δ_pq − Σ_j (2(pq|jj) − (pj|jq))`, so the API
//! stays identical to the previous (spatial-orbital) implementation.
//!
//! The ground-state correlation energy is `λ_min − E_ref` (with `E_ref` the
//! reference determinant energy = the converged RHF energy), not the lowest
//! CIS eigenvalue of the previous version. Validated against PySCF
//! `ci.cisd.CISD` (gate F41, `pyscf_tda.json::cisd_correlation_hartree`).

use super::mo_transform::ao_to_mo_transform;
use nalgebra::{DMatrix, SymmetricEigen};
use serde::{Deserialize, Serialize};

/// CISD excitation.
///
/// Energies are the CI eigenvalues relative to the CI ground state; CISD is
/// not size-consistent for excited states (use the CIS path, F41, when that
/// matters).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CisdExcitation {
    /// Excitation energy (Hartree, relative to the CI ground state).
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
    /// Excitations sorted by energy (relative to the CI ground state).
    pub excitations: Vec<CisdExcitation>,
    /// Ground-state correlation energy (Hartree): λ_min − E_ref.
    pub correlation_energy: f64,
    /// Total number of determinants in the CI space.
    pub n_csfs: usize,
    /// Approximation used: "full" for exact CISD, "perturbative" for CIS+MP2.
    pub approximation: String,
}

const HARTREE_TO_EV: f64 = 27.21138602;

/// Largest determinant space built explicitly; above this the perturbative
/// (CIS + MP2) path is used.
const MAX_DETERMINANTS: usize = 5000;

/// Occupation of one determinant: alpha and beta spin-orbital bitmasks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Determinant {
    alpha: u64,
    beta: u64,
}

impl Determinant {
    #[cfg(test)]
    fn popcount(&self) -> usize {
        self.alpha.count_ones() as usize + self.beta.count_ones() as usize
    }
}

/// Parity of the excitation operator mapping the canonical determinant `from`
/// onto the canonical determinant `to` (allocation-free, mask-based).
fn phase_between(from: &Determinant, to: &Determinant) -> f64 {
    let mut alpha = from.alpha;
    let mut beta = from.beta;
    let mut parity = 0u32;
    // Removals in ascending global order (alpha block first).
    let mut rem = from.alpha & !to.alpha;
    while rem != 0 {
        let bit = rem.trailing_zeros();
        parity += (alpha & ((1u64 << bit) - 1)).count_ones();
        alpha &= !(1u64 << bit);
        rem &= rem - 1;
    }
    let mut rem = from.beta & !to.beta;
    while rem != 0 {
        let bit = rem.trailing_zeros();
        parity += (beta & ((1u64 << bit) - 1)).count_ones();
        beta &= !(1u64 << bit);
        rem &= rem - 1;
    }
    // Additions in *descending* global order, matching the standard
    // Slater-Condon operator a†_a a†_b a_j a_i (i<j, a<b) whose value formula
    // (ia|jb) - (ib|ja) pairs the ascending removed list with the ascending
    // added list. Each creation anticommutes past the occupied orbitals
    // *below* its index (canonical ordering). Verified element-by-element
    // against the explicit second-quantized reference (tests/cisd_brute.rs).
    let mut add = to.alpha & !from.alpha;
    while add != 0 {
        let bit = 63 - add.leading_zeros();
        parity += (alpha & ((1u64 << bit) - 1)).count_ones();
        alpha |= 1u64 << bit;
        add &= !(1u64 << bit);
    }
    let mut add = to.beta & !from.beta;
    while add != 0 {
        let bit = 63 - add.leading_zeros();
        parity += (beta & ((1u64 << bit) - 1)).count_ones();
        beta |= 1u64 << bit;
        add &= !(1u64 << bit);
    }
    if parity.is_multiple_of(2) {
        1.0
    } else {
        -1.0
    }
}

/// Enumerate the CISD determinant space of a closed-shell reference.
fn build_space(n_basis: usize, n_occupied: usize) -> Vec<Determinant> {
    let n_virtual = n_basis - n_occupied;
    let ref_alpha: u64 = if n_occupied >= 64 {
        u64::MAX
    } else {
        (1u64 << n_occupied) - 1
    };
    let reference = Determinant {
        alpha: ref_alpha,
        beta: ref_alpha,
    };
    let mut space = vec![reference];

    let mut occ: Vec<usize> = (0..n_occupied).collect();
    occ.extend((0..n_occupied).map(|k| n_basis + k));
    let virt: Vec<usize> = (0..2 * n_virtual)
        .map(|k| {
            let spin = k / n_virtual;
            let spatial = k % n_virtual + n_occupied;
            spin * n_basis + spatial
        })
        .collect();
    let spin_of = |idx: usize| idx / n_basis;

    let set = |det: &mut Determinant, idx: usize, value: bool| {
        let mask = 1u64 << (idx % n_basis);
        if idx / n_basis == 0 {
            if value {
                det.alpha |= mask;
            } else {
                det.alpha &= !mask;
            }
        } else if value {
            det.beta |= mask;
        } else {
            det.beta &= !mask;
        }
    };

    // Singles (spin-conserving).
    for &i in &occ {
        for &a in &virt {
            if spin_of(i) != spin_of(a) {
                continue;
            }
            let mut det = reference;
            set(&mut det, i, false);
            set(&mut det, a, true);
            space.push(det);
        }
    }

    // Doubles (spin-conserving: the spin pattern of the occupied pair must
    // match the spin pattern of the virtual pair).
    for (ii, &i) in occ.iter().enumerate() {
        for &j in occ.iter().skip(ii + 1) {
            let n_alpha_occ = (spin_of(i) == 0) as usize + (spin_of(j) == 0) as usize;
            for (aa, &a) in virt.iter().enumerate() {
                for &b in virt.iter().skip(aa + 1) {
                    let n_alpha_virt = (spin_of(a) == 0) as usize + (spin_of(b) == 0) as usize;
                    if n_alpha_occ != n_alpha_virt {
                        continue;
                    }
                    let mut det = reference;
                    set(&mut det, i, false);
                    set(&mut det, j, false);
                    set(&mut det, a, true);
                    set(&mut det, b, true);
                    space.push(det);
                }
            }
        }
    }
    space
}

/// Compute CISD energies from converged SCF data.
///
/// `orbital_energies`: canonical RHF orbital energies (Hartree, spatial MOs).
/// `coefficients`: MO coefficient matrix (n_ao × n_mo).
/// `eris`: AO ERIs in the packed layout of [`super::integrals::get_eri`].
/// `n_occupied`: number of occupied *spatial* orbitals.
/// `n_states`: number of excited states to report.
pub fn compute_cisd(
    orbital_energies: &[f64],
    coefficients: &DMatrix<f64>,
    eris: &[f64],
    n_basis: usize,
    n_occupied: usize,
    n_states: usize,
) -> CisdResult {
    let n_virtual = n_basis - n_occupied;
    if n_occupied == 0 || n_virtual == 0 || n_basis > 64 {
        return compute_cisd_perturbative(
            orbital_energies,
            coefficients,
            eris,
            n_basis,
            n_occupied,
            n_states,
        );
    }

    let space = build_space(n_basis, n_occupied);
    if space.len() > MAX_DETERMINANTS {
        return compute_cisd_perturbative(
            orbital_energies,
            coefficients,
            eris,
            n_basis,
            n_occupied,
            n_states,
        );
    }

    // Determinant space + Hamiltonian (shared builder).
    let (space, h) = build_hamiltonian(orbital_energies, coefficients, eris, n_basis, n_occupied);
    let n_dets = space.len();
    let e_ref = h[(0, 0)];

    // Diagonalize.
    let eigen = SymmetricEigen::new(h);
    let mut indices: Vec<usize> = (0..n_dets).collect();
    indices.sort_by(|&a, &b| {
        eigen.eigenvalues[a]
            .partial_cmp(&eigen.eigenvalues[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let e_ground = eigen.eigenvalues[indices[0]];
    let correlation_energy = e_ground - e_ref;

    // Excited states relative to the CI ground state.
    let n_out = n_states.min(n_dets.saturating_sub(1));
    let mut excitations = Vec::with_capacity(n_out);
    for &idx in indices.iter().skip(1).take(n_out) {
        let energy = eigen.eigenvalues[idx] - e_ground;
        if energy <= 0.0 {
            continue;
        }
        let energy_ev = energy * HARTREE_TO_EV;
        let wavelength_nm = 1239.84193 / energy_ev;

        let col = eigen.eigenvectors.column(idx);
        // Determinant 0 is the reference, 1..=n_singles the singles.
        let n_singles = 2 * n_occupied * n_virtual;
        let singles_end = (1 + n_singles).min(n_dets);
        let singles_weight: f64 = col.rows(1, singles_end - 1).iter().map(|c| c * c).sum();

        let mut max_coeff = 0.0f64;
        let mut dom_idx = 0;
        for (k, &c) in col.iter().enumerate() {
            if c.abs() > max_coeff {
                max_coeff = c.abs();
                dom_idx = k;
            }
        }

        let mut transitions = Vec::new();
        if dom_idx > 0 {
            let det = space[dom_idx];
            let reference = space[0];
            for spin in 0..2 {
                let (det_mask, ref_mask) = if spin == 0 {
                    (det.alpha, reference.alpha)
                } else {
                    (det.beta, reference.beta)
                };
                let mut removed = det_mask & !ref_mask;
                let mut added = ref_mask & !det_mask;
                while removed != 0 && added != 0 {
                    transitions.push((
                        removed.trailing_zeros() as usize,
                        added.trailing_zeros() as usize,
                    ));
                    removed &= removed - 1;
                    added &= added - 1;
                }
            }
        }

        let character = if singles_weight > 0.5 {
            "singles"
        } else {
            "doubles"
        };

        // Crude oscillator proxy (no dipole integrals in this API): the
        // transition is only bright through its singles component.
        let f_osc = 2.0 / 3.0 * energy * singles_weight;

        excitations.push(CisdExcitation {
            energy,
            energy_ev,
            wavelength_nm,
            oscillator_strength: f_osc,
            dominant_transition: transitions,
            character: character.to_string(),
            singles_weight,
        });
    }

    CisdResult {
        excitations,
        correlation_energy,
        n_csfs: n_dets,
        approximation: "full".to_string(),
    }
}

/// Perturbative CISD approximation for large CI spaces.
///
/// Uses CIS for excitations + the MP2 correlation energy for the doubles
/// contribution.
fn compute_cisd_perturbative(
    orbital_energies: &[f64],
    coefficients: &DMatrix<f64>,
    eris: &[f64],
    n_basis: usize,
    n_occupied: usize,
    n_states: usize,
) -> CisdResult {
    let n_virtual = n_basis - n_occupied;
    if n_virtual == 0 {
        return CisdResult {
            excitations: Vec::new(),
            correlation_energy: 0.0,
            n_csfs: 0,
            approximation: "perturbative".to_string(),
        };
    }

    let mo = ao_to_mo_transform(coefficients, eris, n_basis, n_basis);

    // Closed-shell MP2 correlation energy (spin-summed form, spatial MOs):
    //   E = Σ_ijab (ia|jb)[2(ia|jb) − (ib|ja)] / (ε_i+ε_j−ε_a−ε_b)
    // with the sum over all occupied/virtual spatial orbitals (the spin
    // summation is already contained in the 2(ia|jb) − (ib|ja) factor).
    // Verified against PySCF `mp.MP2`: H2O/STO-3G −0.0355456516 Ha.
    let mut e_corr = 0.0;
    for i in 0..n_occupied {
        for j in 0..n_occupied {
            for a in n_occupied..n_basis {
                for b in n_occupied..n_basis {
                    let ijab = mo.get(i, a, j, b);
                    let ijba = mo.get(i, b, j, a);
                    let denom = orbital_energies[i] + orbital_energies[j]
                        - orbital_energies[a]
                        - orbital_energies[b];
                    if denom >= -1e-12 {
                        continue;
                    }
                    e_corr += ijab * (2.0 * ijab - ijba) / denom;
                }
            }
        }
    }

    let cis = super::cis::compute_cis(
        orbital_energies,
        coefficients,
        eris,
        n_basis,
        n_occupied,
        n_states,
    );

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

    // Determinant count of the CISD space (same convention as the full path):
    // reference + singles + (αα→αα) + (αβ→αβ) + (ββ→ββ) doubles.
    let n_pairs_occ = n_occupied * (n_occupied - 1) / 2;
    let n_pairs_virt = n_virtual * (n_virtual - 1) / 2;
    let n_dets = 1
        + 2 * n_occupied * n_virtual
        + n_occupied * n_occupied * n_virtual * n_virtual
        + 2 * n_pairs_occ * n_pairs_virt;

    CisdResult {
        excitations,
        correlation_energy: e_corr,
        n_csfs: n_dets,
        approximation: "perturbative".to_string(),
    }
}

/// Debug helper: determinant space as (alpha_mask, beta_mask) pairs.
#[doc(hidden)]
pub fn debug_space(n_basis: usize, n_occupied: usize) -> Vec<(u64, u64)> {
    build_space(n_basis, n_occupied)
        .into_iter()
        .map(|d| (d.alpha, d.beta))
        .collect()
}

/// Build the CISD determinant space and Hamiltonian matrix.
///
/// Shared by [`compute_cisd`] and the cross-check test that compares the
/// matrix against an explicit second-quantized reference
/// (`tests/cisd_brute.rs`, exposed through [`debug_matrix`]).
fn build_hamiltonian(
    orbital_energies: &[f64],
    coefficients: &DMatrix<f64>,
    eris: &[f64],
    n_basis: usize,
    n_occupied: usize,
) -> (Vec<Determinant>, DMatrix<f64>) {
    let space = build_space(n_basis, n_occupied);
    let mo = ao_to_mo_transform(coefficients, eris, n_basis, n_basis);
    let spin_of = |idx: usize| idx / n_basis;
    let spatial_of = |idx: usize| idx % n_basis;
    let mut h_mo = DMatrix::zeros(n_basis, n_basis);
    for p in 0..n_basis {
        for q in 0..n_basis {
            let mut fock = 0.0;
            for j in 0..n_occupied {
                fock += 2.0 * mo.get(p, q, j, j) - mo.get(p, j, j, q);
            }
            h_mo[(p, q)] = if p == q {
                orbital_energies[p] - fock
            } else {
                -fock
            };
        }
    }
    let so_h = |p: usize, q: usize| -> f64 {
        if spin_of(p) == spin_of(q) {
            h_mo[(spatial_of(p), spatial_of(q))]
        } else {
            0.0
        }
    };
    let so_eri = |p: usize, q: usize, r: usize, s: usize| -> f64 {
        if spin_of(p) == spin_of(q) && spin_of(r) == spin_of(s) {
            mo.get(spatial_of(p), spatial_of(q), spatial_of(r), spatial_of(s))
        } else {
            0.0
        }
    };
    let n_dets = space.len();
    let mut occ_lists: Vec<u32> = Vec::with_capacity(n_dets * 2 * n_occupied);
    let mut occ_offsets: Vec<usize> = Vec::with_capacity(n_dets + 1);
    for det in &space {
        occ_offsets.push(occ_lists.len());
        let mut a = det.alpha;
        while a != 0 {
            occ_lists.push(a.trailing_zeros());
            a &= a - 1;
        }
        let mut b = det.beta;
        while b != 0 {
            occ_lists.push(n_basis as u32 + b.trailing_zeros());
            b &= b - 1;
        }
    }
    occ_offsets.push(occ_lists.len());
    let mut h = DMatrix::zeros(n_dets, n_dets);
    for row in 0..n_dets {
        for col in 0..=row {
            let mut val = 0.0;
            if row == col {
                let occ = &occ_lists[occ_offsets[row]..occ_offsets[row + 1]];
                for &p in occ {
                    let p = p as usize;
                    val += so_h(p, p);
                    for &q in occ {
                        let q = q as usize;
                        val += 0.5 * (so_eri(p, p, q, q) - so_eri(p, q, q, p));
                    }
                }
            } else {
                let d1 = space[row];
                let d2 = space[col];
                let diff_alpha = (d1.alpha ^ d2.alpha).count_ones() as usize;
                let diff_beta = (d1.beta ^ d2.beta).count_ones() as usize;
                match diff_alpha + diff_beta {
                    2 => {
                        let (from, to) = if diff_alpha == 2 {
                            (
                                (d2.alpha & !d1.alpha).trailing_zeros() as usize,
                                (d1.alpha & !d2.alpha).trailing_zeros() as usize,
                            )
                        } else {
                            (
                                (d2.beta & !d1.beta).trailing_zeros() as usize + n_basis,
                                (d1.beta & !d2.beta).trailing_zeros() as usize + n_basis,
                            )
                        };
                        let phase = phase_between(&d2, &d1);
                        let mut one = so_h(from, to);
                        for &j in &occ_lists[occ_offsets[col]..occ_offsets[col + 1]] {
                            let j = j as usize;
                            one += so_eri(from, to, j, j) - so_eri(from, j, j, to);
                        }
                        val = phase * one;
                    }
                    4 => {
                        let mut removed = [0usize; 2];
                        let mut added = [0usize; 2];
                        let mut n_rem = 0;
                        let mut n_add = 0;
                        for spin in 0..2 {
                            let (mask2, mask1) = if spin == 0 {
                                (d2.alpha, d1.alpha)
                            } else {
                                (d2.beta, d1.beta)
                            };
                            let mut gone = mask2 & !mask1;
                            while gone != 0 {
                                removed[n_rem] = gone.trailing_zeros() as usize + spin * n_basis;
                                n_rem += 1;
                                gone &= gone - 1;
                            }
                            let mut new = mask1 & !mask2;
                            while new != 0 {
                                added[n_add] = new.trailing_zeros() as usize + spin * n_basis;
                                n_add += 1;
                                new &= new - 1;
                            }
                        }
                        let phase = phase_between(&d2, &d1);
                        val = phase
                            * (so_eri(removed[0], added[0], removed[1], added[1])
                                - so_eri(removed[0], added[1], removed[1], added[0]));
                    }
                    _ => {}
                }
            }
            h[(row, col)] = val;
            h[(col, row)] = val;
        }
    }
    (space, h)
}

/// Debug helper: CISD Hamiltonian matrix (same builder as `compute_cisd`).
#[doc(hidden)]
pub fn debug_matrix(
    orbital_energies: &[f64],
    coefficients: &DMatrix<f64>,
    eris: &[f64],
    n_basis: usize,
    n_occupied: usize,
) -> DMatrix<f64> {
    build_hamiltonian(orbital_energies, coefficients, eris, n_basis, n_occupied).1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_space_h2() {
        // H2/STO-3G: 1 occupied + 1 virtual spatial → 2 singles (α, β) and
        // one same-orbital double (0α0β → 1α1β).
        let space = build_space(2, 1);
        assert_eq!(space.len(), 4);
        let double = space.last().unwrap();
        assert_eq!(double.popcount(), 2);
        assert_eq!(double.alpha, 0b10);
        assert_eq!(double.beta, 0b10);
    }

    #[test]
    fn test_build_space_counts() {
        // 2 occupied + 2 virtual spatial (4 + 4 spin orbitals):
        //   singles: 2 (spin) · 2 · 2 = 8
        //   doubles: (αα→αα) 1 + (αβ→αβ) 4·4 + (ββ→ββ) 1 = 18
        let space = build_space(4, 2);
        assert_eq!(space.len(), 1 + 8 + 18);
        // Every determinant keeps N_alpha = N_beta = 2.
        for det in &space {
            assert_eq!(det.alpha.count_ones(), 2);
            assert_eq!(det.beta.count_ones(), 2);
        }
    }

    #[test]
    fn test_cisd_h2_two_determinants() {
        // H2/STO-3G-like model with the same-orbital double present.
        let orbital_energies = vec![-0.5783, 0.6704];
        // Identity coefficients so the "AO" ERIs below are the MO integrals.
        let c = DMatrix::identity(2, 2);
        let n = 2;
        // Minimal ERIs (values in the ballpark of STO-3G H2), packed layout.
        let n_pair = n * (n + 1) / 2;
        let mut eris = vec![0.0; n_pair * (n_pair + 1) / 2];
        let mut set = |i: usize, j: usize, k: usize, l: usize, v: f64| {
            eris[super::super::integrals::eri_index(i, j, k, l, n)] = v;
        };
        set(0, 0, 0, 0, 0.6746);
        set(0, 0, 1, 1, 0.6636);
        set(1, 1, 1, 1, 0.6975);
        set(0, 1, 0, 1, 0.1813);

        let result = compute_cisd(&orbital_energies, &c, &eris, n, 1, 3);
        assert_eq!(result.approximation, "full");
        assert_eq!(result.n_csfs, 4);
        // The same-orbital double must lower the energy below the reference.
        assert!(
            result.correlation_energy < -0.01,
            "correlation = {}",
            result.correlation_energy
        );
    }

    #[test]
    fn test_cisd_perturbative_fallback() {
        // 20 basis / 4 occupied spatial → 5632 determinants > MAX_DETERMINANTS.
        let n = 20;
        let orbital_energies: Vec<f64> = (0..n).map(|i| -0.5 + 0.05 * i as f64).collect();
        let c = DMatrix::identity(n, n);
        let n_pair = n * (n + 1) / 2;
        let eris = vec![0.0; n_pair * (n_pair + 1) / 2];
        let result = compute_cisd(&orbital_energies, &c, &eris, n, 4, 3);
        assert_eq!(result.approximation, "perturbative");
    }
}
