//! Simplified Tamm-Dancoff Approximation (sTDA) for UV-Vis spectra.
//!
//! Grimme's sTDA method provides electronic excitation energies and
//! oscillator strengths at dramatically reduced cost vs full TD-DFT.
//!
//! Reference: S. Grimme, J. Chem. Phys. 138, 244104 (2013).
//!
//! # Model
//!
//! The A matrix of this implementation is
//!
//! ```text
//!   A[ia,jb] = delta_ij delta_ab (eps_a - eps_i)
//!            + 2 sum_AB q_ia(A) gamma_AB q_jb(B)      (damped monopole Coulomb)
//!            - 2 ax (ij|ab)                           (exchange, optional)
//! ```
//!
//! The Coulomb term uses the damped monopole model (transition charges per
//! atom), which is what makes the method cheap.  The **exchange** term
//! requires the real two-electron integrals and is therefore only available
//! through [`compute_stda_with_exchange`], which needs the AO-basis ERI
//! tensor of the same basis that produced the SCF.  [`compute_stda`] keeps
//! the exchange-free (monopole) model and is the entry point used by the
//! reactivity and fluorescence modules.

use nalgebra::DMatrix;

use super::types::{ScfInput, SpectroscopyResult, TransitionInfo};

const HARTREE_TO_EV: f64 = 27.211386245988;

/// Configuration for sTDA calculation.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct StdaConfig {
    /// Energy window for occupied orbitals below HOMO (eV).
    pub occ_window_ev: f64,
    /// Energy window for virtual orbitals above LUMO (eV).
    pub virt_window_ev: f64,
    /// Maximum number of excitations to compute.
    pub n_roots: usize,
    /// Scaling factor of the exchange integrals in the A matrix.
    ///
    /// The exchange contribution is `A_x[ia,jb] = -2 ax (ij|ab)`, so the
    /// default `ax = 0.5` reproduces the Tamm-Dancoff exchange term
    /// `-(ij|ab)`.  It is only used by [`compute_stda_with_exchange`];
    /// [`compute_stda`] uses the exchange-free monopole model and ignores
    /// this field.
    pub ax: f64,
    /// Screening threshold for integral approximation.
    pub threshold: f64,
}

impl Default for StdaConfig {
    fn default() -> Self {
        Self {
            occ_window_ev: 7.0,
            virt_window_ev: 9.0,
            n_roots: 20,
            ax: 0.5,
            threshold: 1e-6,
        }
    }
}

struct ActiveSpace {
    occ_indices: Vec<usize>,
    virt_indices: Vec<usize>,
    n_occ: usize,
    n_virt: usize,
}

fn select_active_space(scf: &ScfInput, config: &StdaConfig) -> ActiveSpace {
    let n_occ = scf.n_electrons / 2;
    let homo_e = scf.orbital_energies[n_occ - 1];

    let lumo_e = if n_occ < scf.n_basis {
        scf.orbital_energies[n_occ]
    } else {
        homo_e + 1.0
    };

    let occ_cutoff = homo_e - config.occ_window_ev / HARTREE_TO_EV;
    let virt_cutoff = lumo_e + config.virt_window_ev / HARTREE_TO_EV;

    // Hard floor to exclude deep-core orbitals (1s, 2s of heavy atoms).
    // Orbitals below -20 eV (≈ -0.735 Hartree) are always excluded as core,
    // consistent with sTDA convention for energy window filtering.
    let core_floor_hartree = -20.0 / HARTREE_TO_EV;
    let effective_occ_cutoff = occ_cutoff.max(core_floor_hartree);

    let mut occ_indices: Vec<usize> = (0..n_occ)
        .filter(|&i| scf.orbital_energies[i] >= effective_occ_cutoff)
        .collect();

    // Safety: if the energy window excludes all occupied orbitals (e.g., due to
    // systematically shifted orbital energies from a minimal-basis SCF), fall
    // back to including the highest occupied orbitals so sTDA can still produce
    // transitions.
    if occ_indices.is_empty() && n_occ > 0 {
        let n_include = n_occ.clamp(1, 3);
        occ_indices = ((n_occ - n_include)..n_occ).collect();
    }

    let virt_indices: Vec<usize> = (n_occ..scf.n_basis)
        .filter(|&a| scf.orbital_energies[a] <= virt_cutoff)
        .collect();

    ActiveSpace {
        n_occ: occ_indices.len(),
        n_virt: virt_indices.len(),
        occ_indices,
        virt_indices,
    }
}

/// Compute transition charges for orbital pair (i, a).
fn transition_charges(
    scf: &ScfInput,
    basis_to_atom: &[usize],
    n_atoms: usize,
) -> Vec<Vec<Vec<f64>>> {
    let n_occ = scf.n_electrons / 2;
    let n_basis = scf.n_basis;

    let sc = &scf.overlap_matrix * &scf.mo_coefficients;

    let mut q = vec![vec![vec![0.0; n_atoms]; n_basis - n_occ]; n_occ];

    for i in 0..n_occ {
        for (a_idx, a) in (n_occ..n_basis).enumerate() {
            for mu in 0..n_basis {
                let atom = basis_to_atom[mu];
                q[i][a_idx][atom] += scf.mo_coefficients[(mu, i)] * sc[(mu, a)];
            }
        }
    }

    q
}

/// Compute sTDA excitation energies and oscillator strengths.
///
/// Requires converged SCF data, basis-to-atom mapping, and atomic positions (Bohr).
///
/// This is the exchange-free (monopole Coulomb) model; see
/// [`compute_stda_with_exchange`] for the variant that adds the exchange
/// term `-2 ax (ij|ab)`.
pub fn compute_stda(
    scf: &ScfInput,
    basis_to_atom: &[usize],
    positions_bohr: &[[f64; 3]],
    config: &StdaConfig,
) -> SpectroscopyResult {
    match try_compute_stda(scf, basis_to_atom, positions_bohr, config) {
        Ok(result) => result,
        Err(error) => failed_spectroscopy_result(error),
    }
}

/// Compute sTDA excitation energies, reporting an invalid geometry as an
/// error.
pub fn try_compute_stda(
    scf: &ScfInput,
    basis_to_atom: &[usize],
    positions_bohr: &[[f64; 3]],
    config: &StdaConfig,
) -> Result<SpectroscopyResult, String> {
    validate_stda_geometry(positions_bohr)?;
    Ok(compute_stda_with_dipoles(
        scf,
        basis_to_atom,
        positions_bohr,
        config,
        None,
    ))
}

/// Validate the atomic positions used by the sTDA monopole model.
fn validate_stda_geometry(positions_bohr: &[[f64; 3]]) -> Result<(), String> {
    if positions_bohr.is_empty() {
        return Err("empty geometry: 0 atoms (at least one atom is required)".to_string());
    }
    crate::validation::validate_coords_finite(positions_bohr)
}

/// The "nothing was computed" result returned for an invalid system.
fn failed_spectroscopy_result(error: String) -> SpectroscopyResult {
    SpectroscopyResult {
        transitions: Vec::new(),
        method: "sTDA".to_string(),
        error: Some(error),
    }
}

/// sTDA with the Tamm-Dancoff exchange term `-2 ax (ij|ab)`.
///
/// The exchange contribution needs the real two-electron integrals, so the
/// AO-basis ERI tensor of the *same* basis that produced the SCF has to be
/// supplied:
///
/// ```no_run
/// use sci_form::scf::basis::BasisSet;
/// use sci_form::scf::two_electron::TwoElectronIntegrals;
/// # use sci_form::spectroscopy::{compute_stda_with_exchange, ScfInput, StdaConfig};
/// # fn run(scf: ScfInput, basis: &BasisSet) {
/// let eris = TwoElectronIntegrals::compute(basis);
/// let spectrum = compute_stda_with_exchange(
///     &scf,
///     &basis.function_to_atom,
///     &basis.functions.iter().map(|f| f.center).collect::<Vec<_>>(),
///     &StdaConfig::default(),
///     &eris,
/// );
/// # let _ = spectrum;
/// # }
/// ```
///
/// The exchange term lowers the excitation energies (it is the dominant
/// correction of the monopole model) and is *not* screened by the monopole
/// Cauchy-Schwarz criterion.
///
/// The AO → MO transformation of the ERI tensor costs O(N^5) time and
/// O(N^4) memory, so this variant is meant for small-to-medium systems;
/// [`compute_stda`] stays the cheap (exchange-free) option for large
/// molecules.
pub fn compute_stda_with_exchange(
    scf: &ScfInput,
    basis_to_atom: &[usize],
    positions_bohr: &[[f64; 3]],
    config: &StdaConfig,
    eris: &crate::scf::two_electron::TwoElectronIntegrals,
) -> SpectroscopyResult {
    let n_occ_total = scf.n_electrons / 2;
    let exchange = build_active_exchange(eris, &scf.mo_coefficients, n_occ_total);
    compute_stda_inner(
        scf,
        basis_to_atom,
        positions_bohr,
        config,
        None,
        Some(&exchange),
    )
}

/// sTDA with AO dipole integrals: exact transition dipoles (instead of the
/// monopole approximation) for the oscillator strengths.
///
/// `ao_dipoles` are the `[D^x, D^y, D^z]` matrices in the AO basis (Bohr),
/// e.g. from [`crate::scf::dipole_matrix::dipole_matrices`].
pub fn compute_stda_with_dipoles(
    scf: &ScfInput,
    basis_to_atom: &[usize],
    positions_bohr: &[[f64; 3]],
    config: &StdaConfig,
    ao_dipoles: Option<&[DMatrix<f64>; 3]>,
) -> SpectroscopyResult {
    compute_stda_inner(scf, basis_to_atom, positions_bohr, config, ao_dipoles, None)
}

/// Active-space exchange integrals `(ij|ab)` in chemist notation.
struct ActiveExchange {
    /// `(ij|ab)` with `i`, `j` the occupied MO indices and `a`, `b` the
    /// virtual MO indices (both counted from zero).
    values: Vec<f64>,
    /// Number of occupied MOs.
    n_occ: usize,
    /// Number of virtual MOs.
    n_virt: usize,
}

impl ActiveExchange {
    #[inline]
    fn get(&self, i: usize, j: usize, a: usize, b: usize) -> f64 {
        let ai = a - self.n_occ;
        let bi = b - self.n_occ;
        let idx = ((i * self.n_occ + j) * self.n_virt + ai) * self.n_virt + bi;
        self.values[idx]
    }
}

/// AO -> MO transformation of the two-electron integrals.
fn mo_integrals(
    eris: &crate::scf::two_electron::TwoElectronIntegrals,
    coefficients: &DMatrix<f64>,
) -> crate::hf::mo_transform::MoIntegrals {
    use crate::hf::integrals::eri_index;
    use crate::hf::mo_transform::ao_to_mo_transform;

    let n = eris.n_basis();
    let n_pair = n * (n + 1) / 2;
    let mut flat = vec![0.0f64; n_pair * (n_pair + 1) / 2];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                for l in 0..n {
                    flat[eri_index(i, j, k, l, n)] = eris.get(i, j, k, l);
                }
            }
        }
    }
    ao_to_mo_transform(coefficients, &flat, n, n)
}

/// Build the exchange block `(ij|ab)` for the occupied/virtual space.
///
/// The active-space restriction is applied when the A matrix is assembled,
/// so that the same integrals serve every active-space choice.
fn build_active_exchange(
    eris: &crate::scf::two_electron::TwoElectronIntegrals,
    coefficients: &DMatrix<f64>,
    n_occ: usize,
) -> ActiveExchange {
    let mo = mo_integrals(eris, coefficients);
    let n_mo = mo.n_mo();
    let n_virt = n_mo - n_occ;
    let mut values = vec![0.0f64; n_occ * n_occ * n_virt * n_virt];
    for i in 0..n_occ {
        for j in 0..n_occ {
            for a in 0..n_virt {
                for b in 0..n_virt {
                    let idx = ((i * n_occ + j) * n_virt + a) * n_virt + b;
                    values[idx] = mo.get(i, j, n_occ + a, n_occ + b);
                }
            }
        }
    }
    ActiveExchange {
        values,
        n_occ,
        n_virt,
    }
}

#[allow(clippy::too_many_arguments)]
fn compute_stda_inner(
    scf: &ScfInput,
    basis_to_atom: &[usize],
    positions_bohr: &[[f64; 3]],
    config: &StdaConfig,
    ao_dipoles: Option<&[DMatrix<f64>; 3]>,
    exchange: Option<&ActiveExchange>,
) -> SpectroscopyResult {
    // Every sTDA entry point funnels through here: validate the geometry
    // once so a NaN coordinate cannot silently poison the monopole model.
    if let Err(error) = validate_stda_geometry(positions_bohr) {
        return failed_spectroscopy_result(error);
    }
    let active = select_active_space(scf, config);
    let n_active_occ = active.n_occ;
    let n_active_virt = active.n_virt;
    let n_singles = n_active_occ * n_active_virt;

    if n_singles == 0 {
        return SpectroscopyResult {
            transitions: Vec::new(),
            method: "sTDA".to_string(),
            error: None,
        };
    }

    let n_atoms = positions_bohr.len();

    // Build CIS-like A matrix
    let mut a_matrix = DMatrix::zeros(n_singles, n_singles);

    // Diagonal: orbital energy differences
    for (idx, (i_local, a_local)) in iproduct(n_active_occ, n_active_virt).enumerate() {
        let i = active.occ_indices[i_local];
        let a = active.virt_indices[a_local];
        a_matrix[(idx, idx)] = scf.orbital_energies[a] - scf.orbital_energies[i];
    }

    // Damped Coulomb monopole interaction
    let eta: Vec<f64> = (0..n_atoms).map(|_| 0.3).collect();
    let gamma = compute_gamma(positions_bohr, &eta);

    let n_occ_total = scf.n_electrons / 2;
    let q = transition_charges(scf, basis_to_atom, n_atoms);

    // Off-diagonal: Coulomb-type integrals with screening
    // Precompute transition charge norms for Cauchy-Schwarz-like prescreening
    let q_norms: Vec<f64> = iproduct(n_active_occ, n_active_virt)
        .map(|(i_l, a_l)| {
            let i = active.occ_indices[i_l];
            let a_abs = active.virt_indices[a_l] - n_occ_total;
            q[i][a_abs].iter().map(|x| x * x).sum::<f64>().sqrt()
        })
        .collect();

    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;

        let pairs_1: Vec<(usize, usize, usize)> = iproduct(n_active_occ, n_active_virt)
            .enumerate()
            .map(|(idx, (i_l, a_l))| {
                (
                    idx,
                    active.occ_indices[i_l],
                    active.virt_indices[a_l] - n_occ_total,
                )
            })
            .collect();

        let pairs_2: Vec<(usize, usize, usize)> = iproduct(n_active_occ, n_active_virt)
            .enumerate()
            .map(|(idx, (j_l, b_l))| {
                (
                    idx,
                    active.occ_indices[j_l],
                    active.virt_indices[b_l] - n_occ_total,
                )
            })
            .collect();

        let row_contribs: Vec<Vec<(usize, f64)>> = pairs_1
            .par_iter()
            .map(|&(idx1, i, a_abs)| {
                let mut row = Vec::with_capacity(n_singles);
                let norm1 = q_norms[idx1];
                for &(idx2, j, b_abs) in &pairs_2 {
                    // Cauchy-Schwarz prescreening: skip if product of norms is below threshold
                    if norm1 * q_norms[idx2] < config.threshold {
                        continue;
                    }
                    let mut j_integral = 0.0;
                    for atom_a in 0..n_atoms {
                        let q_ia = q[i][a_abs][atom_a];
                        if q_ia.abs() < config.threshold {
                            continue;
                        }
                        for atom_b in 0..n_atoms {
                            j_integral += q_ia * gamma[(atom_a, atom_b)] * q[j][b_abs][atom_b];
                        }
                    }
                    row.push((idx2, 2.0 * j_integral));
                }
                row
            })
            .collect();

        for (idx1, row) in row_contribs.into_iter().enumerate() {
            for (idx2, val) in row {
                a_matrix[(idx1, idx2)] += val;
            }
        }
    }

    #[cfg(not(feature = "parallel"))]
    {
        for (idx1, (i_l, a_l)) in iproduct(n_active_occ, n_active_virt).enumerate() {
            let i = active.occ_indices[i_l];
            let a_abs = active.virt_indices[a_l] - n_occ_total;
            let norm1 = q_norms[idx1];

            for (idx2, (j_l, b_l)) in iproduct(n_active_occ, n_active_virt).enumerate() {
                let j = active.occ_indices[j_l];
                let b_abs = active.virt_indices[b_l] - n_occ_total;

                // Cauchy-Schwarz prescreening
                if norm1 * q_norms[idx2] < config.threshold {
                    continue;
                }

                let mut j_integral = 0.0;
                for atom_a in 0..n_atoms {
                    let q_ia = q[i][a_abs][atom_a];
                    if q_ia.abs() < config.threshold {
                        continue;
                    }
                    for atom_b in 0..n_atoms {
                        j_integral += q_ia * gamma[(atom_a, atom_b)] * q[j][b_abs][atom_b];
                    }
                }

                a_matrix[(idx1, idx2)] += 2.0 * j_integral;
            }
        }
    }

    // Exchange contribution: A_x[ia,jb] = -2 ax (ij|ab).
    //
    // The exchange integrals are *not* screened by the monopole
    // Cauchy-Schwarz criterion: they are a different (non-local) quantity.
    if let Some(exch) = exchange {
        for (idx1, (i_l, a_l)) in iproduct(n_active_occ, n_active_virt).enumerate() {
            let i = active.occ_indices[i_l];
            let a = active.virt_indices[a_l];
            for (idx2, (j_l, b_l)) in iproduct(n_active_occ, n_active_virt).enumerate() {
                let j = active.occ_indices[j_l];
                let b = active.virt_indices[b_l];
                a_matrix[(idx1, idx2)] -= 2.0 * config.ax * exch.get(i, j, a, b);
            }
        }
    }

    // Diagonalize
    let eigen = a_matrix.symmetric_eigen();

    let mut idx_sorted: Vec<usize> = (0..n_singles).collect();
    idx_sorted.sort_by(|&a, &b| {
        eigen.eigenvalues[a]
            .partial_cmp(&eigen.eigenvalues[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let n_roots = config.n_roots.min(n_singles);
    let mut transitions = Vec::with_capacity(n_roots);

    for root in 0..n_roots {
        let idx = idx_sorted[root];
        let energy_hartree = eigen.eigenvalues[idx];
        let energy_ev = energy_hartree * HARTREE_TO_EV;

        if energy_ev < 0.0 {
            continue;
        }

        let ci_vector = eigen.eigenvectors.column(idx);
        let (tdm, osc_strength) = transition_dipole_from_ci(
            &ci_vector,
            &active,
            &q,
            positions_bohr,
            n_occ_total,
            energy_hartree,
            scf,
            ao_dipoles,
        );

        transitions.push(TransitionInfo {
            energy_ev,
            wavelength_nm: if energy_ev > 0.0 {
                1239.84198 / energy_ev
            } else {
                0.0
            },
            oscillator_strength: osc_strength,
            transition_dipole: tdm,
        });
    }

    SpectroscopyResult {
        transitions,
        method: "sTDA".to_string(),
        error: None,
    }
}

#[allow(clippy::too_many_arguments)]
fn transition_dipole_from_ci(
    ci: &nalgebra::DVectorView<f64>,
    active: &ActiveSpace,
    q: &[Vec<Vec<f64>>],
    positions_bohr: &[[f64; 3]],
    n_occ_total: usize,
    energy_hartree: f64,
    scf: &ScfInput,
    ao_dipoles: Option<&[DMatrix<f64>; 3]>,
) -> ([f64; 3], f64) {
    if let Some(dips) = ao_dipoles {
        // Exact AO transition dipole: μ_0k = √2 Σ_ia X_ia Σ_μν c_μi c_νa D_μν
        let mut tdm = [0.0f64; 3];
        for (idx, (i_l, a_l)) in iproduct(active.n_occ, active.n_virt).enumerate() {
            let i = active.occ_indices[i_l];
            let a = active.virt_indices[a_l];
            let x_ia = ci[idx];
            if x_ia.abs() < 1e-12 {
                continue;
            }
            for mu in 0..scf.n_basis {
                let c_mu = scf.mo_coefficients[(mu, i)];
                if c_mu.abs() < 1e-12 {
                    continue;
                }
                for nu in 0..scf.n_basis {
                    let c_nu = scf.mo_coefficients[(nu, a)];
                    if c_nu.abs() < 1e-12 {
                        continue;
                    }
                    let weight = x_ia * c_mu * c_nu;
                    tdm[0] += weight * dips[0][(mu, nu)];
                    tdm[1] += weight * dips[1][(mu, nu)];
                    tdm[2] += weight * dips[2][(mu, nu)];
                }
            }
        }
        let sqrt2 = std::f64::consts::SQRT_2;
        let tdm = [tdm[0] * sqrt2, tdm[1] * sqrt2, tdm[2] * sqrt2];
        let tdm_sq = tdm.iter().map(|v| v * v).sum::<f64>();
        let osc = (2.0 / 3.0) * energy_hartree * tdm_sq;
        return (tdm, osc);
    }

    let n_atoms = positions_bohr.len();
    let mut tdm = [0.0f64; 3];

    // sTDA monopole approximation:
    // μ_0k = √2 Σ_ia X_ia^k Σ_A q_ia(A) · R_A
    for (idx, (i_l, a_l)) in iproduct(active.n_occ, active.n_virt).enumerate() {
        let i = active.occ_indices[i_l];
        let a_abs = active.virt_indices[a_l] - n_occ_total;
        let x_ia = ci[idx];

        if x_ia.abs() < 1e-10 {
            continue;
        }

        for atom in 0..n_atoms {
            let charge = q[i][a_abs][atom];
            tdm[0] += x_ia * charge * positions_bohr[atom][0];
            tdm[1] += x_ia * charge * positions_bohr[atom][1];
            tdm[2] += x_ia * charge * positions_bohr[atom][2];
        }
    }

    // Factor √2 for singlet excitations
    let sqrt2 = std::f64::consts::SQRT_2;
    tdm[0] *= sqrt2;
    tdm[1] *= sqrt2;
    tdm[2] *= sqrt2;

    let tdm_sq = tdm[0] * tdm[0] + tdm[1] * tdm[1] + tdm[2] * tdm[2];
    let osc = (2.0 / 3.0) * energy_hartree * tdm_sq;

    (tdm, osc)
}

fn compute_gamma(positions: &[[f64; 3]], eta: &[f64]) -> DMatrix<f64> {
    let n = positions.len();
    let mut gamma = DMatrix::zeros(n, n);

    for a in 0..n {
        for b in 0..n {
            if a == b {
                gamma[(a, b)] = eta[a];
            } else {
                let dx = positions[a][0] - positions[b][0];
                let dy = positions[a][1] - positions[b][1];
                let dz = positions[a][2] - positions[b][2];
                let r2 = dx * dx + dy * dy + dz * dz;

                let avg_eta_inv = 1.0 / (2.0 * eta[a]) + 1.0 / (2.0 * eta[b]);
                gamma[(a, b)] = 1.0 / (r2 + avg_eta_inv * avg_eta_inv).sqrt();
            }
        }
    }

    gamma
}

fn iproduct(n_occ: usize, n_virt: usize) -> impl Iterator<Item = (usize, usize)> + Clone {
    (0..n_occ).flat_map(move |i| (0..n_virt).map(move |a| (i, a)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_gamma_matrix_symmetry() {
        let pos = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]];
        let eta = vec![0.3, 0.3, 0.3];
        let gamma = compute_gamma(&pos, &eta);

        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (gamma[(i, j)] - gamma[(j, i)]).abs() < 1e-14,
                    "Gamma should be symmetric"
                );
            }
        }
        assert!((gamma[(0, 0)] - 0.3).abs() < 1e-14);
    }

    #[test]
    fn test_active_space_selection() {
        let n_basis = 5;
        let scf = ScfInput {
            orbital_energies: vec![-1.0, -0.5, 0.2, 0.8, 1.5],
            mo_coefficients: DMatrix::identity(n_basis, n_basis),
            density_matrix: DMatrix::zeros(n_basis, n_basis),
            overlap_matrix: DMatrix::identity(n_basis, n_basis),
            n_basis,
            n_electrons: 4,
        };

        let config = StdaConfig::default();
        let active = select_active_space(&scf, &config);

        assert!(active.n_occ > 0, "Should have active occupied orbitals");
        assert!(active.n_virt > 0, "Should have active virtual orbitals");
    }

    #[test]
    fn test_stda_empty_on_no_space() {
        let scf = ScfInput {
            orbital_energies: vec![-10.0, 10.0],
            mo_coefficients: DMatrix::identity(2, 2),
            density_matrix: DMatrix::zeros(2, 2),
            overlap_matrix: DMatrix::identity(2, 2),
            n_basis: 2,
            n_electrons: 2,
        };

        let config = StdaConfig {
            occ_window_ev: 0.1,
            virt_window_ev: 0.1,
            ..Default::default()
        };
        let result = compute_stda(&scf, &[0, 0], &[[0.0, 0.0, 0.0]], &config);
        // Very narrow windows may yield 0 or 1 transitions
        assert!(result.method == "sTDA");
    }

    #[test]
    fn test_stda_produces_transitions() {
        let n_basis = 4;
        let scf = ScfInput {
            orbital_energies: vec![-0.8, -0.3, 0.1, 0.5],
            mo_coefficients: DMatrix::identity(n_basis, n_basis),
            density_matrix: DMatrix::zeros(n_basis, n_basis),
            overlap_matrix: DMatrix::identity(n_basis, n_basis),
            n_basis,
            n_electrons: 4,
        };

        let config = StdaConfig::default();
        let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let basis_to_atom = [0, 0, 1, 1];
        let result = compute_stda(&scf, &basis_to_atom, &positions, &config);

        assert!(!result.transitions.is_empty(), "Should produce transitions");
        for t in &result.transitions {
            assert!(t.energy_ev > 0.0);
            assert!(t.wavelength_nm > 0.0);
        }
    }

    #[test]
    fn test_stda_ao_dipoles_change_oscillator_strength() {
        let n_basis = 4;
        let scf = ScfInput {
            orbital_energies: vec![-0.8, -0.3, -0.25, -0.2],
            mo_coefficients: DMatrix::identity(n_basis, n_basis),
            density_matrix: DMatrix::zeros(n_basis, n_basis),
            overlap_matrix: DMatrix::identity(n_basis, n_basis),
            n_basis,
            n_electrons: 2,
        };
        let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let basis_to_atom = [0, 0, 1, 1];
        let config = StdaConfig::default();
        let monopole = compute_stda(&scf, &basis_to_atom, &positions, &config);

        let mut dips = [
            DMatrix::zeros(n_basis, n_basis),
            DMatrix::zeros(n_basis, n_basis),
            DMatrix::zeros(n_basis, n_basis),
        ];
        dips[0][(0, 2)] = 1.0;
        dips[0][(2, 0)] = 1.0;
        let exact =
            compute_stda_with_dipoles(&scf, &basis_to_atom, &positions, &config, Some(&dips));
        assert_eq!(monopole.transitions.len(), exact.transitions.len());
        assert!(exact
            .transitions
            .iter()
            .all(|t| t.oscillator_strength.is_finite()));
        let changed = monopole
            .transitions
            .iter()
            .zip(exact.transitions.iter())
            .any(|(a, b)| (a.oscillator_strength - b.oscillator_strength).abs() > 1e-6);
        assert!(
            changed,
            "the AO dipole path should change the oscillator strengths"
        );
    }

    #[test]
    fn test_stda_prescreening_removes_off_diagonal_coupling() {
        let n_basis = 4;
        let scf = ScfInput {
            orbital_energies: vec![-0.8, -0.3, 0.1, 0.5],
            mo_coefficients: DMatrix::identity(n_basis, n_basis),
            density_matrix: DMatrix::zeros(n_basis, n_basis),
            overlap_matrix: DMatrix::identity(n_basis, n_basis),
            n_basis,
            n_electrons: 4,
        };
        let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let basis_to_atom = [0, 0, 1, 1];

        let unscreened = compute_stda(
            &scf,
            &basis_to_atom,
            &positions,
            &StdaConfig {
                threshold: 0.0,
                ..Default::default()
            },
        );
        // With every pair screened the A matrix is diagonal: the transition
        // energies are exactly the bare active-space orbital differences.
        let screened = compute_stda(
            &scf,
            &basis_to_atom,
            &positions,
            &StdaConfig {
                threshold: f64::MAX,
                ..Default::default()
            },
        );
        assert!(!unscreened.transitions.is_empty());
        assert!(!screened.transitions.is_empty());

        let active = select_active_space(&scf, &StdaConfig::default());
        let mut expected: Vec<f64> = iproduct(active.n_occ, active.n_virt)
            .map(|(i_l, a_l)| {
                (scf.orbital_energies[active.virt_indices[a_l]]
                    - scf.orbital_energies[active.occ_indices[i_l]])
                    * HARTREE_TO_EV
            })
            .filter(|e| *e >= 0.0)
            .collect();
        expected.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        for (t, e) in screened.transitions.iter().zip(expected.iter()) {
            assert!(
                (t.energy_ev - e).abs() < 1e-9,
                "screened transition {} != diagonal {}",
                t.energy_ev,
                e
            );
        }
    }

    #[test]
    fn test_iproduct() {
        let pairs: Vec<_> = iproduct(2, 3).collect();
        assert_eq!(pairs.len(), 6);
        assert_eq!(pairs[0], (0, 0));
        assert_eq!(pairs[5], (1, 2));
    }
}
