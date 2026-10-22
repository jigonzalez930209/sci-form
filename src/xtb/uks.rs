//! Unrestricted (spin-polarized) GFN2-xTB — collinear UKS with α/β densities.
//!
//! Structural mirror of `gfn2::solve_gfn2` (same H0/S/CN/γ/AES/D4 machinery,
//! same Hartree/bohr internals and eV outputs) with a spin-resolved SCC loop:
//! separate α/β densities (`Pa`, `Pb`), per-shell magnetizations
//! `m = pop_α − pop_β`, and the tblite on-site spin interaction
//! (`tblite_spin`, LGPL-3.0-or-later):
//!
//!   E_spin = +½ Σ_A Σ_{l,l'} m_{Al} W^A_{ll'} m_{Al'}
//!   F_α ⊃ +½·S·(V_i + V_j),  F_β ⊃ −½·S·(V_i + V_j),  V = W·m
//!
//! (signs from δE_spin/δPσ; the majority spin is lowered since W < 0).
//! AES multipoles, third order, D4 and repulsion use total quantities only,
//! exactly as in tblite. Spin constants come from `spin_params` (transcribed
//! verbatim from tblite `data/spin.f90`).
//!
//! With `uhf = 0` the spin shifts vanish and the solver reproduces
//! `solve_gfn2` (guarded by `uks_uhf0_matches_restricted`); with `uhf = Nα−Nβ`
//! it reproduces `tblite.interface.Calculator(..., uhf=N)` (gate F45-UKS).

use super::broyden::BroydenMixer;
use super::d4::D4Model;
use super::gfn0::{BOHR_PER_ANGSTROM, EV_PER_HARTREE};
use super::gfn1::{
    build_gfn1_basis, build_overlap_hamiltonian, compute_halogen_bond, compute_repulsion,
    gfn1_valence, nsao,
};
use super::gfn1_params::{get_gfn1_params, GFN1_GLOBALS, GFN1_NVALENCE};
use super::gfn2::{
    build_gfn2_basis_map, compute_coordination_numbers, compute_gfn2_repulsion,
    compute_halogen_bond_correction, get_fermi_filling,
};
use super::gfn2_params::{self, EV_TO_HARTREE, GAM3_SHELL};
use super::solver::ANGSTROM_TO_BOHR;
use super::spin_params::spin_w;
use super::sto_overlap::sto_multipole_with_ng;
use crate::hf::d3::{compute_d3_energy_with_params, coordination_numbers_exp, D3Params};
use nalgebra::DMatrix;

/// Default electronic temperature in Hartree (≈300 K, as in `gfn2`).
const ELECTRONIC_KT: f64 = 9.500e-4;

/// Euclidean distance between two points already in bohr.
fn dist_raw(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Per-shell data (as in `gfn2`, plus the angular momentum for W lookup).
struct ShellInfo {
    atom: usize,
    l: u8,
    eta: f64,
    ref_occ: f64,
    gam3_sh: f64,
}

/// Result of an unrestricted GFN2-xTB calculation (energies in eV).
pub struct UksResult {
    /// Total energy (eV).
    pub total_energy: f64,
    /// Electronic energy incl. spin polarization (eV).
    pub electronic_energy: f64,
    /// Repulsive energy (eV).
    pub repulsive_energy: f64,
    /// Three-body (ATM) D4 dispersion energy (eV).
    pub dispersion_energy: f64,
    /// Self-consistent two-body D4 energy (eV, in `electronic_energy`).
    pub dispersion_2_body_energy: f64,
    /// Halogen bond correction (eV).
    pub halogen_bond_energy: f64,
    /// Spin-polarization energy ½ΣmWm (eV, in `electronic_energy`).
    pub spin_energy: f64,
    /// Mulliken charges per atom (total), in units of `e`.
    pub mulliken_charges: Vec<f64>,
    /// Magnetization per atom (pop_α − pop_β), in units of `e`.
    pub spin_moments: Vec<f64>,
    /// Magnetization per shell, in units of `e`.
    pub shell_magnetization: Vec<f64>,
    /// Number of alpha electrons.
    pub n_alpha: f64,
    /// Number of beta electrons.
    pub n_beta: f64,
    /// Alpha HOMO/LUMO (eV).
    pub homo_alpha: f64,
    /// Alpha LUMO (eV).
    pub lumo_alpha: f64,
    /// Beta HOMO (eV).
    pub homo_beta: f64,
    /// Beta LUMO (eV).
    pub lumo_beta: f64,
    /// Alpha/beta orbital energies (eV).
    pub orbital_energies_alpha: Vec<f64>,
    /// Beta orbital energies (eV).
    pub orbital_energies_beta: Vec<f64>,
    /// SCC iterations.
    pub scc_iterations: usize,
    /// Whether the SCC converged.
    pub converged: bool,
}

/// Build the SCC Fock matrix for one spin channel.
///
/// `shifts` are the per-shell potentials entering with the same convention
/// as the restricted solver (`F = H0 − S·shifts − AES terms`).
#[allow(clippy::too_many_arguments)]
fn build_fock(
    h_mat: &DMatrix<f64>,
    s_mat: &DMatrix<f64>,
    shifts: &[f64],
    vat_mp: &[f64],
    vdp: &[[f64; 3]],
    vqp: &[[f64; 6]],
    dpint: &[DMatrix<f64>; 3],
    qpint: &[DMatrix<f64>; 6],
    basis_to_shell: &[usize],
    basis_map: &[(usize, usize, u8, u8)],
    n_basis: usize,
) -> DMatrix<f64> {
    let mut h_scc = DMatrix::zeros(n_basis, n_basis);
    for i in 0..n_basis {
        let shi = basis_to_shell[i];
        let ai = basis_map[i].0;
        let mut fii = h_mat[(i, i)] - shifts[shi] - vat_mp[ai];
        for k in 0..3 {
            fii -= dpint[k][(i, i)] * vdp[ai][k];
        }
        for c in 0..6 {
            fii -= qpint[c][(i, i)] * vqp[ai][c];
        }
        h_scc[(i, i)] = fii;
        for j in (i + 1)..n_basis {
            let shj = basis_to_shell[j];
            let aj = basis_map[j].0;
            let mut fij = h_mat[(i, j)] - 0.5 * s_mat[(i, j)] * (shifts[shi] + shifts[shj]);
            fij -= 0.5 * s_mat[(i, j)] * (vat_mp[ai] + vat_mp[aj]);
            for k in 0..3 {
                fij -= 0.5 * (dpint[k][(i, j)] * vdp[ai][k] + dpint[k][(j, i)] * vdp[aj][k]);
            }
            for c in 0..6 {
                fij -= 0.5 * (qpint[c][(i, j)] * vqp[ai][c] + qpint[c][(j, i)] * vqp[aj][c]);
            }
            h_scc[(i, j)] = fij;
            h_scc[(j, i)] = fij;
        }
    }
    h_scc
}

/// Löwdin diagonalization returning ascending eigenvalues + AO coefficients.
fn diagonalize(
    s_half_inv: &DMatrix<f64>,
    fock: &DMatrix<f64>,
    n_basis: usize,
) -> (Vec<f64>, DMatrix<f64>) {
    let f_prime = s_half_inv * fock * s_half_inv;
    let eigen = f_prime.symmetric_eigen();
    let mut indices: Vec<usize> = (0..n_basis).collect();
    indices.sort_by(|&a, &b| {
        eigen.eigenvalues[a]
            .partial_cmp(&eigen.eigenvalues[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut evals = vec![0.0f64; n_basis];
    let c_prime = &eigen.eigenvectors;
    let c_full = s_half_inv * c_prime;
    let mut coeff = DMatrix::zeros(n_basis, n_basis);
    for (new_idx, &old_idx) in indices.iter().enumerate() {
        evals[new_idx] = eigen.eigenvalues[old_idx];
        for i in 0..n_basis {
            coeff[(i, new_idx)] = c_full[(i, old_idx)];
        }
    }
    (evals, coeff)
}

/// Density matrix from occupations (no factor-2; caller handles the spin).
fn density_from(coeff: &DMatrix<f64>, focc: &[f64], n_basis: usize) -> DMatrix<f64> {
    let mut density = DMatrix::zeros(n_basis, n_basis);
    for i in 0..n_basis {
        for j in 0..n_basis {
            let mut val = 0.0;
            for k in 0..n_basis {
                val += focc[k] * coeff[(i, k)] * coeff[(j, k)];
            }
            density[(i, j)] = val;
        }
    }
    density
}

/// Solve unrestricted GFN2-xTB for `uhf = Nα − Nβ` unpaired electrons.
///
/// Total charge is neutral (as in `solve_gfn2`); open-shell anions/cations
/// are out of scope. `uhf = 0` reproduces the restricted solver.
///
/// `positions` are Cartesian coordinates in Ångström. The energies in the
/// returned [`UksResult`] are in eV.
///
/// # Errors
///
/// Returns `Err` if `elements` and `positions` have different lengths, if
/// an element has no GFN2 parameters, if `uhf` is incompatible with the
/// electron count (`N ± uhf` must both be even and `0 <= uhf <= N`), or if
/// the basis set ends up empty.
///
/// # Panics
///
/// Panics if an element has no GFN2 parameters. This cannot happen for the
/// molecules this function accepts, because the parameter lookup is
/// validated before the calculation starts.
pub fn solve_gfn2_unrestricted(
    elements: &[u8],
    positions: &[[f64; 3]],
    uhf: u32,
) -> Result<UksResult, String> {
    if elements.len() != positions.len() {
        return Err(format!(
            "elements ({}) and positions ({}) length mismatch",
            elements.len(),
            positions.len()
        ));
    }
    for &z in elements {
        if gfn2_params::get_gfn2_params(z).is_none() {
            return Err(format!("Element Z={} not supported by GFN2-xTB", z));
        }
    }

    let n_atoms = elements.len();
    let basis_map = build_gfn2_basis_map(elements);
    let n_basis = basis_map.len();
    let n_electrons = gfn2_params::gfn2_electron_count(elements);
    let n_el = n_electrons as i64;
    let nuhf = uhf as i64;
    if (n_el + nuhf) % 2 != 0 || (n_el - nuhf) % 2 != 0 || nuhf < 0 || nuhf > n_el {
        return Err(format!(
            "uhf={uhf} incompatible with {n_electrons} electrons (need N±uhf even, 0<=uhf<=N)"
        ));
    }
    let n_alpha = ((n_el + nuhf) / 2) as f64;
    let n_beta = ((n_el - nuhf) / 2) as f64;
    if n_basis == 0 {
        return Err("No basis functions for GFN2-xTB".to_string());
    }

    // ── Setup (structural copy of `gfn2::solve_gfn2`) ──
    let pos_bohr: Vec<[f64; 3]> = positions
        .iter()
        .map(|p| {
            [
                p[0] * ANGSTROM_TO_BOHR,
                p[1] * ANGSTROM_TO_BOHR,
                p[2] * ANGSTROM_TO_BOHR,
            ]
        })
        .collect();

    let mut s_mat = DMatrix::zeros(n_basis, n_basis);
    let mut dpint = [
        DMatrix::zeros(n_basis, n_basis),
        DMatrix::zeros(n_basis, n_basis),
        DMatrix::zeros(n_basis, n_basis),
    ];
    let mut qpint = [
        DMatrix::zeros(n_basis, n_basis),
        DMatrix::zeros(n_basis, n_basis),
        DMatrix::zeros(n_basis, n_basis),
        DMatrix::zeros(n_basis, n_basis),
        DMatrix::zeros(n_basis, n_basis),
        DMatrix::zeros(n_basis, n_basis),
    ];
    for i in 0..n_basis {
        s_mat[(i, i)] = 1.0;
        let (atom_a, sh_a, l_a, m_a) = basis_map[i];
        let pa = gfn2_params::get_gfn2_params(elements[atom_a]).unwrap();
        {
            let za = pa.slater[sh_a];
            if za > 1e-10 {
                let (_, _, qij) = sto_multipole_with_ng(
                    pa.pqn[sh_a],
                    l_a,
                    m_a,
                    za,
                    &pos_bohr[atom_a],
                    pa.ngauss[sh_a],
                    pa.pqn[sh_a],
                    l_a,
                    m_a,
                    za,
                    &pos_bohr[atom_a],
                    pa.ngauss[sh_a],
                );
                for c in 0..6 {
                    qpint[c][(i, i)] = qij[c];
                }
            }
        }
        for j in (i + 1)..n_basis {
            let (atom_b, sh_b, l_b, m_b) = basis_map[j];
            let pb = gfn2_params::get_gfn2_params(elements[atom_b]).unwrap();
            let za = pa.slater[sh_a];
            let zb = pb.slater[sh_b];
            if za < 1e-10 || zb < 1e-10 {
                continue;
            }
            if atom_a == atom_b {
                let (_, dij, qij) = sto_multipole_with_ng(
                    pa.pqn[sh_a],
                    l_a,
                    m_a,
                    za,
                    &pos_bohr[atom_a],
                    pa.ngauss[sh_a],
                    pb.pqn[sh_b],
                    l_b,
                    m_b,
                    zb,
                    &pos_bohr[atom_b],
                    pb.ngauss[sh_b],
                );
                for alpha in 0..3 {
                    dpint[alpha][(i, j)] = dij[alpha];
                    dpint[alpha][(j, i)] = dij[alpha];
                }
                for c in 0..6 {
                    qpint[c][(i, j)] = qij[c];
                    qpint[c][(j, i)] = qij[c];
                }
                continue;
            }
            let za = pa.slater[sh_a];
            let zb = pb.slater[sh_b];
            if za < 1e-10 || zb < 1e-10 {
                continue;
            }
            let (sij, dij, qij) = sto_multipole_with_ng(
                pa.pqn[sh_a],
                l_a,
                m_a,
                za,
                &pos_bohr[atom_a],
                pa.ngauss[sh_a],
                pb.pqn[sh_b],
                l_b,
                m_b,
                zb,
                &pos_bohr[atom_b],
                pb.ngauss[sh_b],
            );
            s_mat[(i, j)] = sij;
            s_mat[(j, i)] = sij;
            let rab = [
                pos_bohr[atom_a][0] - pos_bohr[atom_b][0],
                pos_bohr[atom_a][1] - pos_bohr[atom_b][1],
                pos_bohr[atom_a][2] - pos_bohr[atom_b][2],
            ];
            for alpha in 0..3 {
                dpint[alpha][(i, j)] = dij[alpha];
                dpint[alpha][(j, i)] = dij[alpha] + rab[alpha] * sij;
            }
            for c in 0..6 {
                qpint[c][(i, j)] = qij[c];
            }
            let shift_raw = [
                2.0 * dij[0] * rab[0] + sij * rab[0] * rab[0],
                dij[0] * rab[1] + dij[1] * rab[0] + sij * rab[0] * rab[1],
                2.0 * dij[1] * rab[1] + sij * rab[1] * rab[1],
                dij[0] * rab[2] + dij[2] * rab[0] + sij * rab[0] * rab[2],
                dij[1] * rab[2] + dij[2] * rab[1] + sij * rab[1] * rab[2],
                2.0 * dij[2] * rab[2] + sij * rab[2] * rab[2],
            ];
            let tr_s = 0.5 * (shift_raw[0] + shift_raw[2] + shift_raw[5]);
            qpint[0][(j, i)] = qij[0] + 1.5 * shift_raw[0] - tr_s;
            qpint[1][(j, i)] = qij[1] + 1.5 * shift_raw[1];
            qpint[2][(j, i)] = qij[2] + 1.5 * shift_raw[2] - tr_s;
            qpint[3][(j, i)] = qij[3] + 1.5 * shift_raw[3];
            qpint[4][(j, i)] = qij[4] + 1.5 * shift_raw[4];
            qpint[5][(j, i)] = qij[5] + 1.5 * shift_raw[5] - tr_s;
        }
    }

    let cn = compute_coordination_numbers(elements, &pos_bohr);

    let mut h_mat = DMatrix::zeros(n_basis, n_basis);
    for i in 0..n_basis {
        let (atom_a, sh_a, _, _) = basis_map[i];
        let pa = gfn2_params::get_gfn2_params(elements[atom_a]).unwrap();
        let se_eff = (pa.selfenergy[sh_a] - pa.kcn[sh_a] * cn[atom_a]) * EV_TO_HARTREE;
        h_mat[(i, i)] = se_eff;
    }
    for i in 0..n_basis {
        let (atom_a, sh_a, l_a, _) = basis_map[i];
        let pa = gfn2_params::get_gfn2_params(elements[atom_a]).unwrap();
        for j in (i + 1)..n_basis {
            let (atom_b, sh_b, l_b, _) = basis_map[j];
            if atom_a == atom_b {
                continue;
            }
            let pb = gfn2_params::get_gfn2_params(elements[atom_b]).unwrap();
            if s_mat[(i, j)].abs() < 1e-15 {
                continue;
            }
            let zi = pa.slater[sh_a];
            let zj = pb.slater[sh_b];
            let zij = if zi > 1e-10 && zj > 1e-10 {
                let sum = zi + zj;
                (2.0 * (zi * zj).sqrt() / sum).sqrt()
            } else {
                1.0
            };
            let den = pa.pauling_en - pb.pauling_en;
            let enp = 1.0 + 0.02 * den * den;
            let k = gfn2_params::kshell_factor(l_a, l_b);
            let r_bohr = dist_raw(&pos_bohr[atom_a], &pos_bohr[atom_b]);
            let rad_sum = pa.atomic_rad + pb.atomic_rad;
            let rr = if rad_sum > 1e-10 {
                (r_bohr / rad_sum).sqrt()
            } else {
                0.0
            };
            let shpoly = (1.0 + pa.shpoly[sh_a] * rr) * (1.0 + pb.shpoly[sh_b] * rr);
            let se_a = h_mat[(i, i)];
            let se_b = h_mat[(j, j)];
            let hij = 0.5 * (se_a + se_b) * s_mat[(i, j)] * zij * k * enp * shpoly;
            h_mat[(i, j)] = hij;
            h_mat[(j, i)] = hij;
        }
    }

    let e_rep = compute_gfn2_repulsion(elements, positions);
    let d4_model = D4Model::new(elements, &pos_bohr);

    let mut shells: Vec<ShellInfo> = Vec::new();
    let mut basis_to_shell: Vec<usize> = Vec::with_capacity(n_basis);
    for (i, &z) in elements.iter().enumerate() {
        let p = gfn2_params::get_gfn2_params(z).unwrap();
        for sh in 0..p.nshell {
            let sh_idx = shells.len();
            let l = p.ang_shell[sh];
            let kind_idx = (p.kind as usize).saturating_sub(1);
            let l_idx = (l as usize).min(3);
            let gam3 = p.gam3 * GAM3_SHELL[l_idx][kind_idx];
            shells.push(ShellInfo {
                atom: i,
                l,
                eta: p.hubbard * p.shell_hubbard[sh],
                ref_occ: p.ref_occ[sh],
                gam3_sh: gam3,
            });
            let n_m = 2 * l as usize + 1;
            for _ in 0..n_m {
                basis_to_shell.push(sh_idx);
            }
        }
    }
    let n_shells = shells.len();

    let gexp = 2.0f64;
    let mut gamma_sh = vec![vec![0.0f64; n_shells]; n_shells];
    for i in 0..n_shells {
        let si = &shells[i];
        gamma_sh[i][i] = si.eta;
        for j in (i + 1)..n_shells {
            let sj = &shells[j];
            let eta_avg = 0.5 * (si.eta + sj.eta);
            if si.atom == sj.atom {
                gamma_sh[i][j] = eta_avg;
                gamma_sh[j][i] = eta_avg;
            } else {
                let r_bohr = dist_raw(&pos_bohr[si.atom], &pos_bohr[sj.atom]);
                let gamma = (r_bohr.powf(gexp) + eta_avg.powf(-gexp)).powf(-1.0 / gexp);
                gamma_sh[i][j] = gamma;
                gamma_sh[j][i] = gamma;
            }
        }
    }

    let s_eigen = s_mat.clone().symmetric_eigen();
    let mut s_half_inv = DMatrix::zeros(n_basis, n_basis);
    for k in 0..n_basis {
        let val = s_eigen.eigenvalues[k];
        if val > 1e-8 {
            let inv_sqrt = 1.0 / val.sqrt();
            let col = s_eigen.eigenvectors.column(k);
            for i in 0..n_basis {
                for j in 0..n_basis {
                    s_half_inv[(i, j)] += inv_sqrt * col[i] * col[j];
                }
            }
        }
    }

    let mut mrad = vec![0.0f64; n_atoms];
    for iat in 0..n_atoms {
        let p = gfn2_params::get_gfn2_params(elements[iat]).unwrap();
        let arg = cn[iat] - p.mp_vcn - gfn2_params::MP_SHIFT;
        let t1 = (-gfn2_params::MP_KEXP_MULTI * arg).exp();
        let t2 = (gfn2_params::MP_RMAX - p.mp_rad) / (1.0 + t1);
        mrad[iat] = p.mp_rad + t2;
    }

    let mut amat_sd = vec![vec![[0.0f64; 3]; n_atoms]; n_atoms];
    let mut amat_dd = vec![vec![[[0.0f64; 3]; 3]; n_atoms]; n_atoms];
    let mut amat_sq = vec![vec![[0.0f64; 6]; n_atoms]; n_atoms];
    for iat in 0..n_atoms {
        for jat in 0..n_atoms {
            if iat == jat {
                continue;
            }
            let vec3 = [
                pos_bohr[iat][0] - pos_bohr[jat][0],
                pos_bohr[iat][1] - pos_bohr[jat][1],
                pos_bohr[iat][2] - pos_bohr[jat][2],
            ];
            let r1 = (vec3[0] * vec3[0] + vec3[1] * vec3[1] + vec3[2] * vec3[2]).sqrt();
            let g1 = 1.0 / r1;
            let g3 = g1 * g1 * g1;
            let g5 = g3 * g1 * g1;
            let rr = 0.5 * (mrad[jat] + mrad[iat]) * g1;
            let fdmp3 = 1.0 / (1.0 + 6.0 * rr.powf(gfn2_params::MP_DMP3));
            let fdmp5 = 1.0 / (1.0 + 6.0 * rr.powf(gfn2_params::MP_DMP5));
            for k in 0..3 {
                amat_sd[iat][jat][k] += vec3[k] * g3 * fdmp3;
            }
            for a in 0..3 {
                for b in 0..3 {
                    let mut val = -3.0 * vec3[a] * vec3[b] * g5 * fdmp5;
                    if a == b {
                        val += g3 * fdmp5;
                    }
                    amat_dd[iat][jat][a][b] += val;
                }
            }
            amat_sq[iat][jat][0] += vec3[0] * vec3[0] * g5 * fdmp5;
            amat_sq[iat][jat][1] += 2.0 * vec3[0] * vec3[1] * g5 * fdmp5;
            amat_sq[iat][jat][2] += vec3[1] * vec3[1] * g5 * fdmp5;
            amat_sq[iat][jat][3] += 2.0 * vec3[0] * vec3[2] * g5 * fdmp5;
            amat_sq[iat][jat][4] += 2.0 * vec3[1] * vec3[2] * g5 * fdmp5;
            amat_sq[iat][jat][5] += vec3[2] * vec3[2] * g5 * fdmp5;
        }
    }

    // ── Spin-resolved SCC loop ──
    let max_iter = 250;
    let mut sh_charges = vec![0.0f64; n_shells];
    let mut sh_mag = vec![0.0f64; n_shells];
    let mut atom_charges = vec![0.0f64; n_atoms];
    let mut evals_a = vec![0.0f64; n_basis];
    let mut evals_b = vec![0.0f64; n_basis];
    let mut converged = false;
    let mut scc_iter = 0;
    let mut prev_e_elec = 0.0f64;
    let mut e_d4_2body_ha = 0.0f64;
    let mut sh_shifts = vec![0.0f64; n_shells];
    let mut dpat = vec![[0.0f64; 3]; n_atoms];
    let mut qpat = vec![[0.0f64; 6]; n_atoms];
    let mut vdp = vec![[0.0f64; 3]; n_atoms];
    let mut vqp = vec![[0.0f64; 6]; n_atoms];
    let mut vat_mp = vec![0.0f64; n_atoms];

    let mix_len = 2 * n_shells + n_atoms * 3 + n_atoms * 6;
    let mut mixer = BroydenMixer::new(mix_len, max_iter, 0.4);
    let flatten = |sh: &[f64], mag: &[f64], dp: &[[f64; 3]], qp: &[[f64; 6]]| -> Vec<f64> {
        let mut v = Vec::with_capacity(mix_len);
        v.extend_from_slice(sh);
        v.extend_from_slice(mag);
        for d in dp {
            v.extend_from_slice(d);
        }
        for q in qp {
            v.extend_from_slice(q);
        }
        v
    };

    for iter in 0..max_iter {
        scc_iter = iter + 1;
        if iter > 0 {
            mixer.step().unwrap_or(());
            let mut mixed = vec![0.0f64; mix_len];
            mixer.get(&mut mixed);
            sh_charges.copy_from_slice(&mixed[..n_shells]);
            sh_mag.copy_from_slice(&mixed[n_shells..2 * n_shells]);
            for iat in 0..n_atoms {
                for k in 0..3 {
                    dpat[iat][k] = mixed[2 * n_shells + iat * 3 + k];
                }
                for c in 0..6 {
                    qpat[iat][c] = mixed[2 * n_shells + n_atoms * 3 + iat * 6 + c];
                }
            }
        }

        for i in 0..n_shells {
            sh_shifts[i] = 0.0;
            for j in 0..n_shells {
                sh_shifts[i] += gamma_sh[i][j] * sh_charges[j];
            }
            sh_shifts[i] += shells[i].gam3_sh * sh_charges[i] * sh_charges[i];
        }

        // On-site spin potential V[sh] = Σ_{sh' ∈ same atom} W_{ll'} m_{sh'}.
        let mut v_spin = vec![0.0f64; n_shells];
        for sh in 0..n_shells {
            let at = shells[sh].atom;
            let z = elements[at];
            let mut v = 0.0;
            for shp in 0..n_shells {
                if shells[shp].atom != at {
                    continue;
                }
                v += spin_w(shells[sh].l, shells[shp].l, z) * sh_mag[shp];
            }
            v_spin[sh] = v;
        }
        // F_α carries −(Vq − ½Vspin), F_β −(Vq + ½Vspin): the majority spin
        // is lowered (W < 0), as in tblite (V_α = ½(Vq + Vm)).
        let mut shifts_a = vec![0.0f64; n_shells];
        let mut shifts_b = vec![0.0f64; n_shells];
        for sh in 0..n_shells {
            shifts_a[sh] = sh_shifts[sh] - 0.5 * v_spin[sh];
            shifts_b[sh] = sh_shifts[sh] + 0.5 * v_spin[sh];
        }

        // AES + D4 potentials from TOTAL quantities (as in tblite).
        let mut qat = vec![0.0f64; n_atoms];
        for sh in 0..n_shells {
            qat[shells[sh].atom] += sh_charges[sh];
        }
        for iat in 0..n_atoms {
            vdp[iat] = [0.0; 3];
            vqp[iat] = [0.0; 6];
            vat_mp[iat] = 0.0;
        }
        for iat in 0..n_atoms {
            for jat in 0..n_atoms {
                if iat == jat {
                    continue;
                }
                for k in 0..3 {
                    vdp[iat][k] += amat_sd[jat][iat][k] * qat[jat];
                    for l in 0..3 {
                        vdp[iat][k] += amat_dd[iat][jat][k][l] * dpat[jat][l];
                    }
                }
                for c in 0..6 {
                    vqp[iat][c] += amat_sq[iat][jat][c] * qat[jat];
                }
                for k in 0..3 {
                    vat_mp[iat] += amat_sd[iat][jat][k] * dpat[jat][k];
                }
                for c in 0..6 {
                    vat_mp[iat] += amat_sq[iat][jat][c] * qpat[jat][c];
                }
            }
        }
        let qp_scale: [f64; 6] = [1.0, 2.0, 1.0, 2.0, 2.0, 1.0];
        for iat in 0..n_atoms {
            let p = gfn2_params::get_gfn2_params(elements[iat]).unwrap();
            for k in 0..3 {
                vdp[iat][k] += 2.0 * p.dkernel * dpat[iat][k];
            }
            for c in 0..6 {
                vqp[iat][c] += 2.0 * p.qkernel * qpat[iat][c] * qp_scale[c];
            }
        }
        let d4_weights = d4_model.weight_references(&qat);
        let d4_vat = d4_model.get_potential(&d4_weights);
        for iat in 0..n_atoms {
            vat_mp[iat] += d4_vat[iat];
        }

        let f_a = build_fock(
            &h_mat,
            &s_mat,
            &shifts_a,
            &vat_mp,
            &vdp,
            &vqp,
            &dpint,
            &qpint,
            &basis_to_shell,
            &basis_map,
            n_basis,
        );
        let f_b = build_fock(
            &h_mat,
            &s_mat,
            &shifts_b,
            &vat_mp,
            &vdp,
            &vqp,
            &dpint,
            &qpint,
            &basis_to_shell,
            &basis_map,
            n_basis,
        );
        let (ea, ca) = diagonalize(&s_half_inv, &f_a, n_basis);
        let (eb, cb) = diagonalize(&s_half_inv, &f_b, n_basis);
        evals_a = ea;
        evals_b = eb;

        let (focc_a, ts_a) = get_fermi_filling(&evals_a, n_alpha, ELECTRONIC_KT);
        let (focc_b, ts_b) = get_fermi_filling(&evals_b, n_beta, ELECTRONIC_KT);
        let pa = density_from(&ca, &focc_a, n_basis);
        let pb = density_from(&cb, &focc_b, n_basis);

        // Shell charges (total) and magnetizations.
        let mut new_sh_charges = vec![0.0f64; n_shells];
        let mut new_sh_mag = vec![0.0f64; n_shells];
        for mu in 0..n_basis {
            let sh = basis_to_shell[mu];
            let mut pop_a = 0.0;
            let mut pop_b = 0.0;
            for nu in 0..n_basis {
                pop_a += pa[(mu, nu)] * s_mat[(nu, mu)];
                pop_b += pb[(mu, nu)] * s_mat[(nu, mu)];
            }
            new_sh_charges[sh] -= pop_a + pop_b;
            new_sh_mag[sh] += pop_a - pop_b;
        }
        for sh in 0..n_shells {
            new_sh_charges[sh] += shells[sh].ref_occ;
        }

        // AES multipoles from the TOTAL density.
        let mut pt = pa.clone();
        pt += &pb;
        let mut new_dpat = vec![[0.0f64; 3]; n_atoms];
        let mut new_qpat = vec![[0.0f64; 6]; n_atoms];
        for mu in 0..n_basis {
            let iat = basis_map[mu].0;
            for nu in 0..n_basis {
                let p_nu_mu = pt[(nu, mu)];
                if p_nu_mu.abs() < 1e-15 {
                    continue;
                }
                for k in 0..3 {
                    new_dpat[iat][k] -= p_nu_mu * dpint[k][(mu, nu)];
                }
                for c in 0..6 {
                    new_qpat[iat][c] -= p_nu_mu * qpint[c][(mu, nu)];
                }
            }
        }

        // Energies (Hartree).
        let mut e_band = 0.0;
        for i in 0..n_basis {
            for j in 0..n_basis {
                e_band += pa[(i, j)] * h_mat[(i, j)] + pb[(i, j)] * h_mat[(i, j)];
            }
        }
        let mut e_scc = 0.0;
        for si in 0..n_shells {
            for sj in 0..n_shells {
                e_scc += gamma_sh[si][sj] * new_sh_charges[si] * new_sh_charges[sj];
            }
        }
        e_scc *= 0.5;
        let mut e_3rd = 0.0;
        for sh in 0..n_shells {
            e_3rd += shells[sh].gam3_sh * new_sh_charges[sh].powi(3);
        }
        e_3rd /= 3.0;
        // Spin energy ½Σ_A Σ_{ll'} m W m (NEW magnetizations, as for charges).
        let mut e_spin = 0.0;
        for sh in 0..n_shells {
            let at = shells[sh].atom;
            let z = elements[at];
            for shp in 0..n_shells {
                if shells[shp].atom != at {
                    continue;
                }
                e_spin +=
                    0.5 * new_sh_mag[sh] * spin_w(shells[sh].l, shells[shp].l, z) * new_sh_mag[shp];
            }
        }

        let mut new_qat = vec![0.0f64; n_atoms];
        for sh in 0..n_shells {
            new_qat[shells[sh].atom] += new_sh_charges[sh];
        }
        let mut e_aes = 0.0;
        for iat in 0..n_atoms {
            let mut vd_iat = [0.0f64; 3];
            let mut vq_iat = [0.0f64; 6];
            for jat in 0..n_atoms {
                if iat == jat {
                    continue;
                }
                for k in 0..3 {
                    vd_iat[k] += amat_sd[jat][iat][k] * new_qat[jat];
                    for l in 0..3 {
                        vd_iat[k] += 0.5 * amat_dd[iat][jat][k][l] * new_dpat[jat][l];
                    }
                }
                for c in 0..6 {
                    vq_iat[c] += amat_sq[iat][jat][c] * new_qat[jat];
                }
            }
            for k in 0..3 {
                e_aes += new_dpat[iat][k] * vd_iat[k];
            }
            for c in 0..6 {
                e_aes += new_qpat[iat][c] * vq_iat[c];
            }
        }
        for iat in 0..n_atoms {
            let p = gfn2_params::get_gfn2_params(elements[iat]).unwrap();
            let mut dk_e = 0.0;
            for k in 0..3 {
                dk_e += new_dpat[iat][k] * new_dpat[iat][k];
            }
            e_aes += p.dkernel * dk_e;
            let mut qk_e = 0.0;
            for c in 0..6 {
                qk_e += new_qpat[iat][c] * new_qpat[iat][c] * qp_scale[c];
            }
            e_aes += p.qkernel * qk_e;
        }

        let e_elec = e_band + e_scc + e_3rd + e_aes + e_spin + ts_a + ts_b;
        let d4_w_energy = d4_model.weight_references(&new_qat);
        let e_d4_sc = d4_model.get_energy(&d4_w_energy);
        let e_elec = e_elec + e_d4_sc;
        e_d4_2body_ha = e_d4_sc;

        let x_in_tmp = flatten(&sh_charges, &sh_mag, &dpat, &qpat);
        let x_out_tmp = flatten(&new_sh_charges, &new_sh_mag, &new_dpat, &new_qpat);
        let rms_dq: f64 = {
            let sum_sq: f64 = x_in_tmp
                .iter()
                .zip(x_out_tmp.iter())
                .map(|(i, o)| (o - i).powi(2))
                .sum();
            (sum_sq / mix_len as f64).sqrt()
        };
        let energy_conv = (e_elec - prev_e_elec).abs() < 1e-6 && iter > 0;
        let charge_conv = rms_dq < 2e-5;

        if energy_conv && charge_conv {
            converged = true;
            prev_e_elec = e_elec;
            sh_charges = new_sh_charges;
            sh_mag = new_sh_mag;
            break;
        }
        prev_e_elec = e_elec;

        let x_in = flatten(&sh_charges, &sh_mag, &dpat, &qpat);
        mixer.set(&x_in);
        let x_out = flatten(&new_sh_charges, &new_sh_mag, &new_dpat, &new_qpat);
        mixer.diff(&x_out);
    }

    for sh in 0..n_shells {
        atom_charges[shells[sh].atom] += sh_charges[sh];
    }
    let mut spin_moments = vec![0.0f64; n_atoms];
    for sh in 0..n_shells {
        spin_moments[shells[sh].atom] += sh_mag[sh];
    }

    let e_elec_ha = prev_e_elec;
    let disp_ha = d4_model.get_atm_energy(&pos_bohr);
    let xb_ha = compute_halogen_bond_correction(elements, positions);
    let total_ha = e_elec_ha + e_rep + disp_ha + xb_ha;

    // Recompute the spin energy from converged moments for reporting.
    let mut e_spin_rep = 0.0;
    for sh in 0..n_shells {
        let at = shells[sh].atom;
        let z = elements[at];
        for shp in 0..n_shells {
            if shells[shp].atom != at {
                continue;
            }
            e_spin_rep += 0.5 * sh_mag[sh] * spin_w(shells[sh].l, shells[shp].l, z) * sh_mag[shp];
        }
    }

    let orb_a_ev: Vec<f64> = evals_a.iter().map(|e| e * EV_PER_HARTREE).collect();
    let orb_b_ev: Vec<f64> = evals_b.iter().map(|e| e * EV_PER_HARTREE).collect();
    let na_idx = (n_alpha.round() as usize)
        .saturating_sub(1)
        .min(n_basis - 1);
    let nb_idx = (n_beta.round() as usize).saturating_sub(1).min(n_basis - 1);
    let (homo_alpha, lumo_alpha) = (
        orb_a_ev[na_idx],
        if na_idx + 1 < n_basis {
            orb_a_ev[na_idx + 1]
        } else {
            orb_a_ev[na_idx]
        },
    );
    let (homo_beta, lumo_beta) = (
        orb_b_ev[nb_idx],
        if nb_idx + 1 < n_basis {
            orb_b_ev[nb_idx + 1]
        } else {
            orb_b_ev[nb_idx]
        },
    );

    Ok(UksResult {
        total_energy: total_ha * EV_PER_HARTREE,
        electronic_energy: e_elec_ha * EV_PER_HARTREE,
        repulsive_energy: e_rep * EV_PER_HARTREE,
        dispersion_energy: disp_ha * EV_PER_HARTREE,
        dispersion_2_body_energy: e_d4_2body_ha * EV_PER_HARTREE,
        halogen_bond_energy: xb_ha * EV_PER_HARTREE,
        spin_energy: e_spin_rep * EV_PER_HARTREE,
        mulliken_charges: atom_charges,
        spin_moments,
        shell_magnetization: sh_mag,
        n_alpha,
        n_beta,
        homo_alpha,
        lumo_alpha,
        homo_beta,
        lumo_beta,
        orbital_energies_alpha: orb_a_ev,
        orbital_energies_beta: orb_b_ev,
        scc_iterations: scc_iter,
        converged,
    })
}

/// Result of an unrestricted GFN1-xTB calculation (energies in eV).
pub struct Gfn1UksResult {
    /// Total energy (eV).
    pub total_energy: f64,
    /// Electronic energy incl. spin polarization (eV).
    pub electronic_energy: f64,
    /// Repulsive energy (eV).
    pub repulsive_energy: f64,
    /// D3(BJ) dispersion energy (eV).
    pub dispersion_energy: f64,
    /// Halogen bond correction (eV).
    pub halogen_bond_energy: f64,
    /// Spin-polarization energy ½ΣmWm (eV, in `electronic_energy`).
    pub spin_energy: f64,
    /// Mulliken charges per atom (total), in units of `e`.
    pub mulliken_charges: Vec<f64>,
    /// Magnetization per atom (pop_α − pop_β), in units of `e`.
    pub spin_moments: Vec<f64>,
    /// Number of alpha electrons.
    pub n_alpha: f64,
    /// Number of beta electrons.
    pub n_beta: f64,
    /// Alpha HOMO/LUMO (eV).
    pub homo_alpha: f64,
    /// Alpha LUMO (eV).
    pub lumo_alpha: f64,
    /// Beta HOMO (eV).
    pub homo_beta: f64,
    /// Beta LUMO (eV).
    pub lumo_beta: f64,
    /// SCC iterations.
    pub scc_iterations: usize,
    /// Whether the SCC converged.
    pub converged: bool,
}

/// Solve unrestricted GFN1-xTB for `uhf = Nα − Nβ` unpaired electrons.
///
/// Same spin model as the GFN2 solver above (on-site W from `spin_params`,
/// AES-free GFN1 electrostatics on total quantities). `uhf = 0` reproduces
/// `solve_gfn1`.
///
/// `positions` are Cartesian coordinates in Ångström. The energies in the
/// returned [`Gfn1UksResult`] are in eV.
///
/// # Errors
///
/// Returns `Err` if `elements` and `positions` have different lengths, if
/// an element has no GFN1 parameters, or if `uhf` is incompatible with the
/// electron count.
///
/// # Panics
///
/// Panics if an element has no GFN1 parameters. This cannot happen for the
/// molecules this function accepts, because the parameter lookup is
/// validated before the calculation starts.
pub fn solve_gfn1_unrestricted(
    elements: &[u8],
    positions: &[[f64; 3]],
    uhf: u32,
) -> Result<Gfn1UksResult, String> {
    if elements.len() != positions.len() {
        return Err(format!(
            "elements ({}) and positions ({}) length mismatch",
            elements.len(),
            positions.len()
        ));
    }
    for &z in elements {
        if get_gfn1_params(z).is_none() {
            return Err(format!("Element Z={} not supported by GFN1-xTB", z));
        }
    }

    let n_atoms = elements.len();
    let pos_bohr: Vec<[f64; 3]> = positions
        .iter()
        .map(|p| {
            [
                p[0] * BOHR_PER_ANGSTROM,
                p[1] * BOHR_PER_ANGSTROM,
                p[2] * BOHR_PER_ANGSTROM,
            ]
        })
        .collect();

    let basis = build_gfn1_basis(elements);
    let n_basis = basis.n_basis;
    if n_basis == 0 {
        return Err("No basis functions for GFN1-xTB".to_string());
    }

    let n_electrons: usize = elements
        .iter()
        .map(|&z| GFN1_NVALENCE[z as usize - 1] as usize)
        .sum();
    let n_el = n_electrons as i64;
    let nuhf = uhf as i64;
    if (n_el + nuhf) % 2 != 0 || (n_el - nuhf) % 2 != 0 || nuhf < 0 || nuhf > n_el {
        return Err(format!(
            "uhf={uhf} incompatible with {n_electrons} electrons (need N±uhf even, 0<=uhf<=N)"
        ));
    }
    let n_alpha = ((n_el + nuhf) / 2) as f64;
    let n_beta = ((n_el - nuhf) / 2) as f64;

    let cn = coordination_numbers_exp(elements, &pos_bohr);
    let (s_mat, h_mat) = build_overlap_hamiltonian(&basis, elements, &pos_bohr, &cn);

    struct ShellMeta {
        atom: usize,
        l: u8,
        hardness: f64,
        ref_occ: f64,
    }
    let mut shells: Vec<ShellMeta> = Vec::new();
    let mut basis_to_shell = vec![0usize; n_basis];
    for (i, &z) in elements.iter().enumerate() {
        let p = get_gfn1_params(z).expect("GFN1 params");
        let val = gfn1_valence(&p.ang_shell[..p.nshell]);
        for sh in 0..p.nshell {
            let sh_idx = shells.len();
            let l = p.ang_shell[sh];
            shells.push(ShellMeta {
                atom: i,
                l,
                hardness: p.hubbard * p.shell_hubbard[sh],
                ref_occ: if val[sh] { p.refocc[sh] } else { 0.0 },
            });
            for ii in 0..nsao(l) {
                basis_to_shell[basis.shell_ao[i][sh] + ii] = sh_idx;
            }
        }
    }
    let n_shells = shells.len();

    let hd: Vec<f64> = elements
        .iter()
        .map(|&z| get_gfn1_params(z).expect("GFN1 params").hubbard_derivs)
        .collect();

    let gexp = GFN1_GLOBALS.gexp;
    let mut gamma = vec![vec![0.0f64; n_shells]; n_shells];
    for i in 0..n_shells {
        gamma[i][i] = shells[i].hardness;
    }
    for i in 0..n_shells {
        for j in (i + 1)..n_shells {
            let gam = 2.0 / (1.0 / shells[i].hardness + 1.0 / shells[j].hardness);
            let val = if shells[i].atom == shells[j].atom {
                gam
            } else {
                let a = shells[i].atom;
                let b = shells[j].atom;
                let dx = pos_bohr[a][0] - pos_bohr[b][0];
                let dy = pos_bohr[a][1] - pos_bohr[b][1];
                let dz = pos_bohr[a][2] - pos_bohr[b][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                1.0 / (r.powf(gexp) + gam.powf(-gexp)).powf(1.0 / gexp)
            };
            gamma[i][j] = val;
            gamma[j][i] = val;
        }
    }

    let s_eigen = s_mat.clone().symmetric_eigen();
    let mut s_half_inv = DMatrix::zeros(n_basis, n_basis);
    for k in 0..n_basis {
        let val = s_eigen.eigenvalues[k];
        if val > 1e-8 {
            let inv_sqrt = 1.0 / val.sqrt();
            let col = s_eigen.eigenvectors.column(k);
            for i in 0..n_basis {
                for j in 0..n_basis {
                    s_half_inv[(i, j)] += inv_sqrt * col[i] * col[j];
                }
            }
        }
    }

    // ── Spin-resolved SCC loop (mirrors `solve_gfn1`) ──
    let max_iter = 250;
    let mut qsh = vec![0.0f64; n_shells];
    let mut msh = vec![0.0f64; n_shells];
    let mut mixer = BroydenMixer::new(2 * n_shells, max_iter, 0.4);
    let mut evals_a = vec![0.0f64; n_basis];
    let mut evals_b = vec![0.0f64; n_basis];
    let mut converged = false;
    let mut scc_iter = 0;
    let mut e_elec = 0.0f64;
    let mut prev_e = 0.0f64;

    for iter in 0..max_iter {
        scc_iter = iter + 1;
        if iter > 0 {
            mixer.step().unwrap_or(());
            let mut mixed = vec![0.0f64; 2 * n_shells];
            mixer.get(&mut mixed);
            qsh.copy_from_slice(&mixed[..n_shells]);
            msh.copy_from_slice(&mixed[n_shells..]);
        }

        let mut qat = vec![0.0f64; n_atoms];
        for sh in 0..n_shells {
            qat[shells[sh].atom] += qsh[sh];
        }
        let mut vsh = vec![0.0f64; n_shells];
        for i in 0..n_shells {
            let mut acc = 0.0;
            for j in 0..n_shells {
                acc += gamma[i][j] * qsh[j];
            }
            vsh[i] = acc;
        }
        let vat: Vec<f64> = (0..n_atoms).map(|a| qat[a] * qat[a] * hd[a]).collect();

        // On-site spin potential (same W model as GFN2).
        let mut v_spin = vec![0.0f64; n_shells];
        for sh in 0..n_shells {
            let at = shells[sh].atom;
            let z = elements[at];
            let mut v = 0.0;
            for shp in 0..n_shells {
                if shells[shp].atom != at {
                    continue;
                }
                v += spin_w(shells[sh].l, shells[shp].l, z) * msh[shp];
            }
            v_spin[sh] = v;
        }

        // F_α/β = H0 − ½S((vsh∓½Vspin + vat) + (…)).
        let mut f_a = DMatrix::zeros(n_basis, n_basis);
        let mut f_b = DMatrix::zeros(n_basis, n_basis);
        for i in 0..n_basis {
            let sh_i = basis_to_shell[i];
            let at_i = basis.funcs[i].0;
            let via = vsh[sh_i] - 0.5 * v_spin[sh_i] + vat[at_i];
            let vib = vsh[sh_i] + 0.5 * v_spin[sh_i] + vat[at_i];
            for j in i..n_basis {
                let sh_j = basis_to_shell[j];
                let at_j = basis.funcs[j].0;
                let vja = vsh[sh_j] - 0.5 * v_spin[sh_j] + vat[at_j];
                let vjb = vsh[sh_j] + 0.5 * v_spin[sh_j] + vat[at_j];
                let fa = h_mat[(i, j)] - 0.5 * s_mat[(i, j)] * (via + vja);
                let fb = h_mat[(i, j)] - 0.5 * s_mat[(i, j)] * (vib + vjb);
                f_a[(i, j)] = fa;
                f_b[(i, j)] = fb;
                if i != j {
                    f_a[(j, i)] = fa;
                    f_b[(j, i)] = fb;
                }
            }
        }

        let (ea, ca) = diagonalize(&s_half_inv, &f_a, n_basis);
        let (eb, cb) = diagonalize(&s_half_inv, &f_b, n_basis);
        evals_a = ea;
        evals_b = eb;

        let (focc_a, ts_a) = get_fermi_filling(&evals_a, n_alpha, ELECTRONIC_KT);
        let (focc_b, ts_b) = get_fermi_filling(&evals_b, n_beta, ELECTRONIC_KT);
        let pa = density_from(&ca, &focc_a, n_basis);
        let pb = density_from(&cb, &focc_b, n_basis);

        let mut new_qsh = vec![0.0f64; n_shells];
        let mut new_msh = vec![0.0f64; n_shells];
        for sh in 0..n_shells {
            new_qsh[sh] = shells[sh].ref_occ;
        }
        for mu in 0..n_basis {
            let sh = basis_to_shell[mu];
            let mut pop_a = 0.0;
            let mut pop_b = 0.0;
            for nu in 0..n_basis {
                pop_a += pa[(mu, nu)] * s_mat[(nu, mu)];
                pop_b += pb[(mu, nu)] * s_mat[(nu, mu)];
            }
            new_qsh[sh] -= pop_a + pop_b;
            new_msh[sh] += pop_a - pop_b;
        }
        let mut new_qat = vec![0.0f64; n_atoms];
        for sh in 0..n_shells {
            new_qat[shells[sh].atom] += new_qsh[sh];
        }

        let mut e_band = 0.0;
        for i in 0..n_basis {
            for j in 0..n_basis {
                e_band += pa[(i, j)] * h_mat[(i, j)] + pb[(i, j)] * h_mat[(i, j)];
            }
        }
        let mut e2 = 0.0;
        for i in 0..n_shells {
            for j in 0..n_shells {
                e2 += new_qsh[i] * gamma[i][j] * new_qsh[j];
            }
        }
        e2 *= 0.5;
        let mut e3 = 0.0;
        for a in 0..n_atoms {
            e3 += new_qat[a].powi(3) * hd[a];
        }
        e3 /= 3.0;
        let mut e_spin = 0.0;
        for sh in 0..n_shells {
            let at = shells[sh].atom;
            let z = elements[at];
            for shp in 0..n_shells {
                if shells[shp].atom != at {
                    continue;
                }
                e_spin += 0.5 * new_msh[sh] * spin_w(shells[sh].l, shells[shp].l, z) * new_msh[shp];
            }
        }
        e_elec = e_band + e2 + e3 + e_spin + ts_a + ts_b;

        let rms_dq: f64 = {
            let sum_sq: f64 = qsh
                .iter()
                .chain(msh.iter())
                .zip(new_qsh.iter().chain(new_msh.iter()))
                .map(|(o, n)| (n - o).powi(2))
                .sum();
            (sum_sq / (2 * n_shells) as f64).sqrt()
        };
        let e_conv = iter > 0 && (e_elec - prev_e).abs() < 1e-6;
        if e_conv && rms_dq < 2e-5 {
            converged = true;
            qsh = new_qsh;
            msh = new_msh;
            break;
        }
        prev_e = e_elec;

        let mut x_in = Vec::with_capacity(2 * n_shells);
        x_in.extend_from_slice(&qsh);
        x_in.extend_from_slice(&msh);
        mixer.set(&x_in);
        let mut x_out = Vec::with_capacity(2 * n_shells);
        x_out.extend_from_slice(&new_qsh);
        x_out.extend_from_slice(&new_msh);
        mixer.diff(&x_out);
    }

    let mut atom_charges = vec![0.0f64; n_atoms];
    let mut spin_moments = vec![0.0f64; n_atoms];
    for sh in 0..n_shells {
        atom_charges[shells[sh].atom] += qsh[sh];
        spin_moments[shells[sh].atom] += msh[sh];
    }

    let e_rep = compute_repulsion(elements, &pos_bohr);
    let d3 = compute_d3_energy_with_params(elements, &pos_bohr, D3Params::gfn1());
    let e_xb = compute_halogen_bond(elements, &pos_bohr);
    let total_ha = e_elec + e_rep + d3.energy + e_xb;

    let mut e_spin_rep = 0.0;
    for sh in 0..n_shells {
        let at = shells[sh].atom;
        let z = elements[at];
        for shp in 0..n_shells {
            if shells[shp].atom != at {
                continue;
            }
            e_spin_rep += 0.5 * msh[sh] * spin_w(shells[sh].l, shells[shp].l, z) * msh[shp];
        }
    }

    let orb_a_ev: Vec<f64> = evals_a.iter().map(|e| e * EV_PER_HARTREE).collect();
    let orb_b_ev: Vec<f64> = evals_b.iter().map(|e| e * EV_PER_HARTREE).collect();
    let na_idx = (n_alpha.round() as usize)
        .saturating_sub(1)
        .min(n_basis - 1);
    let nb_idx = (n_beta.round() as usize).saturating_sub(1).min(n_basis - 1);
    let (homo_alpha, lumo_alpha) = (
        orb_a_ev[na_idx],
        if na_idx + 1 < n_basis {
            orb_a_ev[na_idx + 1]
        } else {
            orb_a_ev[na_idx]
        },
    );
    let (homo_beta, lumo_beta) = (
        orb_b_ev[nb_idx],
        if nb_idx + 1 < n_basis {
            orb_b_ev[nb_idx + 1]
        } else {
            orb_b_ev[nb_idx]
        },
    );

    Ok(Gfn1UksResult {
        total_energy: total_ha * EV_PER_HARTREE,
        electronic_energy: e_elec * EV_PER_HARTREE,
        repulsive_energy: e_rep * EV_PER_HARTREE,
        dispersion_energy: d3.energy * EV_PER_HARTREE,
        halogen_bond_energy: e_xb * EV_PER_HARTREE,
        spin_energy: e_spin_rep * EV_PER_HARTREE,
        mulliken_charges: atom_charges,
        spin_moments,
        n_alpha,
        n_beta,
        homo_alpha,
        lumo_alpha,
        homo_beta,
        lumo_beta,
        scc_iterations: scc_iter,
        converged,
    })
}

#[cfg(test)]
mod tests {
    #[test]
    fn spin_constants_match_tblite_table() {
        // tblite data/spin.f90 (Hartree): H(Z=1) ss, B(Z=5) sp, C(Z=6) sp/pp, O(Z=8) pp.
        assert!((super::super::spin_params::SPIN_W[0][0] - (-0.0716250)).abs() < 1e-9);
        assert!((super::super::spin_params::SPIN_W[4][1] - (-0.0219370)).abs() < 1e-9);
        assert!((super::super::spin_params::SPIN_W[5][1] - (-0.0250250)).abs() < 1e-9);
        assert!((super::super::spin_params::SPIN_W[5][2] - (-0.0226750)).abs() < 1e-9);
        assert!((super::super::spin_params::SPIN_W[7][2] - (-0.0278500)).abs() < 1e-9);
        // lidx symmetry: W(s,p) == W(p,s).
        assert_eq!(
            super::super::spin_params::spin_w(0, 1, 8),
            super::super::spin_params::spin_w(1, 0, 8)
        );
    }
}
