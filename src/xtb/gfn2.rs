//! GFN2-xTB solver with the official parametrization.
//!
//! Implements the GFN2-xTB tight-binding method using the exact parameters
//! from the original publication (Bannwarth, Ehlert, Grimme, JCTC 2019, 15, 1652).
//!
//! All internal calculations use Hartree/bohr units; final results are in eV.

use super::gfn2_params::{self, EV_TO_HARTREE, GAM3_SHELL};
use super::solver::{distance_bohr, ANGSTROM_TO_BOHR, EV_PER_HARTREE};
use super::sto_overlap::sto_multipole_with_ng;
use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

/// Euclidean distance between two points already in bohr.
fn dist_raw(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Default electronic temperature in Hartree (≈300 K).
/// kB = 3.166808578545117e-6 Ha/K, T = 300 K → kT ≈ 9.500e-4 Ha.
const ELECTRONIC_KT: f64 = 9.500e-4;

/// Compute Fermi occupation numbers for a restricted (spin-paired) system.
///
/// Given orbital energies `emo` (sorted ascending) and total number of
/// electrons `nel`, finds the Fermi level and returns fractional occupations
/// in [0, 1] (per spin-channel; multiply by 2 for total). Also returns the
/// electronic entropy `ts` in Hartree.
fn get_fermi_filling(emo: &[f64], nel: f64, kt: f64) -> (Vec<f64>, f64) {
    let nao = emo.len();
    let mut focc = vec![0.0f64; nao];

    // Aufbau filling to find HOMO
    let homo_float = nel; // nel is already per-spin-channel for restricted
    let homo = homo_float.floor() as usize;
    for i in 0..homo.min(nao) {
        focc[i] = 1.0;
    }
    if homo < nao {
        let frac = nel - homo as f64;
        if frac > 0.0 {
            focc[homo] = frac;
        }
    }

    if kt < 1e-15 || nao == 0 {
        // Zero temperature: just return Aufbau filling + zero entropy
        return (focc, 0.0);
    }

    // Find initial Fermi energy as midpoint of HOMO/LUMO
    let homo_idx = if homo > 0 { homo - 1 } else { 0 };
    let lumo_idx = homo.min(nao - 1);
    let mut e_fermi = 0.5 * (emo[homo_idx] + emo[lumo_idx]);

    // Newton iteration to find Fermi level
    let max_cycle = 200;
    let thr = f64::EPSILON.sqrt();
    let sqrttiny = f64::MIN_POSITIVE.sqrt();

    for _ in 0..max_cycle {
        let mut total_number = 0.0;
        let mut total_dfermi = 0.0;
        for i in 0..nao {
            let x = (emo[i] - e_fermi) / kt;
            if x < 50.0 {
                let ex = x.exp();
                let fermifunct = 1.0 / (ex + 1.0);
                let dfermifunct = ex / (kt * (ex + 1.0).powi(2));
                focc[i] = fermifunct;
                total_number += fermifunct;
                total_dfermi += dfermifunct;
            } else {
                focc[i] = 0.0;
            }
        }
        let change = if total_dfermi > sqrttiny {
            ((nel - total_number) / total_dfermi).clamp(-1.0, 1.0)
        } else {
            0.0
        };
        e_fermi += change;
        if (nel - total_number).abs() <= thr {
            break;
        }
    }

    // Electronic entropy: ts = Σ log(f^f * (1-f)^(1-f)) * kt
    let mut ts = 0.0;
    for i in 0..nao {
        let f = focc[i];
        if f > 1e-15 && (1.0 - f) > 1e-15 {
            ts += f * f.ln() + (1.0 - f) * (1.0 - f).ln();
        }
    }
    ts *= kt;
    // For restricted calculation, entropy from both spin channels
    ts *= 2.0;

    (focc, ts)
}

/// GFN2-xTB calculation result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gfn2Result {
    /// Orbital energies (eV).
    pub orbital_energies: Vec<f64>,
    /// Electronic energy (eV).
    pub electronic_energy: f64,
    /// Repulsive energy (eV).
    pub repulsive_energy: f64,
    /// Dispersion energy (eV) — D4.
    pub dispersion_energy: f64,
    /// Halogen bond correction energy (eV).
    pub halogen_bond_energy: f64,
    /// Total energy (eV).
    pub total_energy: f64,
    /// Number of basis functions.
    pub n_basis: usize,
    /// Number of electrons.
    pub n_electrons: usize,
    /// HOMO energy (eV).
    pub homo_energy: f64,
    /// LUMO energy (eV).
    pub lumo_energy: f64,
    /// HOMO-LUMO gap (eV).
    pub gap: f64,
    /// Mulliken charges per atom.
    pub mulliken_charges: Vec<f64>,
    /// Atomic dipole moments (Debye, per atom).
    pub atomic_dipoles: Vec<[f64; 3]>,
    /// Atomic quadrupole moment traces.
    pub atomic_quadrupoles: Vec<f64>,
    /// SCC iterations.
    pub scc_iterations: usize,
    /// Whether the SCC converged.
    pub converged: bool,
}

/// Build basis map using GFN2 shell structure: (atom_idx, shell_idx, l, m).
fn build_gfn2_basis_map(elements: &[u8]) -> Vec<(usize, usize, u8, u8)> {
    let mut basis = Vec::new();
    for (i, &z) in elements.iter().enumerate() {
        if let Some(p) = gfn2_params::get_gfn2_params(z) {
            for sh in 0..p.nshell {
                let l = p.ang_shell[sh];
                let n_m = 2 * l as usize + 1; // s=1, p=3, d=5
                for m in 0..n_m {
                    basis.push((i, sh, l, m as u8));
                }
            }
        }
    }
    basis
}

/// Solve GFN2-xTB for a molecular system.
///
/// Uses the official GFN2-xTB parametrization with consistent Hartree/bohr
/// internal units. Returns all energies in eV.
pub fn solve_gfn2(elements: &[u8], positions: &[[f64; 3]]) -> Result<Gfn2Result, String> {
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
    let n_occ = n_electrons / 2;
    let debug_diag = std::env::var("GFN2_DEBUG").is_ok();

    if n_basis == 0 {
        return Err("No basis functions for GFN2-xTB".to_string());
    }

    // ── Overlap matrix (STO-nG with GFN2 Slater exponents) ──
    // Positions in bohr for the overlap integral computation
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
    // Dipole integrals: dpint[alpha][(i,j)] for alpha=0,1,2 (x,y,z)
    let mut dpint = [
        DMatrix::zeros(n_basis, n_basis),
        DMatrix::zeros(n_basis, n_basis),
        DMatrix::zeros(n_basis, n_basis),
    ];
    // Quadrupole integrals: qpint[c][(i,j)] for c=0..5 (xx,xy,yy,xz,yz,zz)
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
        // Compute diagonal quadrupole integrals <μ|Θ|μ> (nonzero for p-orbitals).
        // Dipole diagonal is zero by parity; overlap diagonal is 1 (already set).
        {
            let za = pa.slater[sh_a];
            if za > 1e-10 {
                let (_, _, qij) = sto_multipole_with_ng(
                    pa.pqn[sh_a], l_a, m_a, za, &pos_bohr[atom_a], pa.ngauss[sh_a],
                    pa.pqn[sh_a], l_a, m_a, za, &pos_bohr[atom_a], pa.ngauss[sh_a],
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
                // Same-atom: overlap is identity (already set), but dipole/quadrupole
                // integrals <s|(r-R_A)|p> etc. can be nonzero. Compute them.
                let (_, dij, qij) = sto_multipole_with_ng(
                    pa.pqn[sh_a], l_a, m_a, za, &pos_bohr[atom_a], pa.ngauss[sh_a],
                    pb.pqn[sh_b], l_b, m_b, zb, &pos_bohr[atom_b], pb.ngauss[sh_b],
                );
                // Same atom: both about atom_a = atom_b, no shift needed
                for alpha in 0..3 {
                    dpint[alpha][(i, j)] = dij[alpha];
                    dpint[alpha][(j, i)] = dij[alpha]; // symmetric since same center
                }
                for c in 0..6 {
                    qpint[c][(i, j)] = qij[c];
                    qpint[c][(j, i)] = qij[c]; // symmetric since same center
                }
                continue;
            }
            let za = pa.slater[sh_a];
            let zb = pb.slater[sh_b];
            if za < 1e-10 || zb < 1e-10 {
                continue;
            }
            // Use multipole routine to get overlap + dipole + quadrupole in one call
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
            // Dipole integrals: dpint[k][(i,j)] about center of atom_a (atom of i)
            //                   dpint[k][(j,i)] about center of atom_b (atom of j)
            // Shift: d_B = d_A + (R_A - R_B) * S
            let rab = [
                pos_bohr[atom_a][0] - pos_bohr[atom_b][0],
                pos_bohr[atom_a][1] - pos_bohr[atom_b][1],
                pos_bohr[atom_a][2] - pos_bohr[atom_b][2],
            ];
            for alpha in 0..3 {
                dpint[alpha][(i, j)] = dij[alpha];
                dpint[alpha][(j, i)] = dij[alpha] + rab[alpha] * sij;
            }
            // Quadrupole integrals: qpint[c][(i,j)] about atom_a, qpint[c][(j,i)] about atom_b
            // Shift: q_B(a,b) = q_A(a,b) + d_A(a)*ΔR(b) + d_A(b)*ΔR(a) + S*ΔR(a)*ΔR(b)
            // Then apply traceless correction to get the shifted version
            // Indices: 0=xx, 1=xy, 2=yy, 3=xz, 4=yz, 5=zz
            for c in 0..6 {
                qpint[c][(i, j)] = qij[c];
            }
            // Compute shifted traceless quadrupole about B
            // q_B = q_A + shift_terms, where shift terms come from origin change
            let shift_raw = [
                2.0 * dij[0] * rab[0] + sij * rab[0] * rab[0], // xx
                dij[0] * rab[1] + dij[1] * rab[0] + sij * rab[0] * rab[1], // xy
                2.0 * dij[1] * rab[1] + sij * rab[1] * rab[1], // yy
                dij[0] * rab[2] + dij[2] * rab[0] + sij * rab[0] * rab[2], // xz
                dij[1] * rab[2] + dij[2] * rab[1] + sij * rab[1] * rab[2], // yz
                2.0 * dij[2] * rab[2] + sij * rab[2] * rab[2], // zz
            ];
            // Apply traceless projection to the shift
            let tr_s = 0.5 * (shift_raw[0] + shift_raw[2] + shift_raw[5]);
            qpint[0][(j, i)] = qij[0] + 1.5 * shift_raw[0] - tr_s; // xx
            qpint[1][(j, i)] = qij[1] + 1.5 * shift_raw[1];        // xy
            qpint[2][(j, i)] = qij[2] + 1.5 * shift_raw[2] - tr_s; // yy
            qpint[3][(j, i)] = qij[3] + 1.5 * shift_raw[3];        // xz
            qpint[4][(j, i)] = qij[4] + 1.5 * shift_raw[4];        // yz
            qpint[5][(j, i)] = qij[5] + 1.5 * shift_raw[5] - tr_s; // zz
        }
    }

    if debug_diag {
        eprintln!("=== Dipole integral matrix dpint[2] (z-component) ===");
        for i in 0..n_basis {
            let mut row = String::new();
            for j in 0..n_basis {
                row.push_str(&format!("{:10.6} ", dpint[2][(i, j)]));
            }
            eprintln!("dpint_z[{}] = [{}]", i, row.trim());
        }
        eprintln!("=== Quadrupole qpint[5] (zz-component) ===");
        for i in 0..n_basis {
            let mut row = String::new();
            for j in 0..n_basis {
                row.push_str(&format!("{:10.6} ", qpint[5][(i, j)]));
            }
            eprintln!("qpint_zz[{}] = [{}]", i, row.trim());
        }
    }

    // ── Coordination numbers (D3-type, for CN-dependent self-energies) ──
    let cn = compute_coordination_numbers(elements, &pos_bohr);

    // ── Hamiltonian H0 (Hartree) — full tblite formula ──
    let mut h_mat = DMatrix::zeros(n_basis, n_basis);
    // Diagonal: GFN2 self-energies with CN shift (eV → Hartree)
    // SE_eff(sh) = SE(sh) - kCN(sh) * CN(atom)
    for i in 0..n_basis {
        let (atom_a, sh_a, _, _) = basis_map[i];
        let pa = gfn2_params::get_gfn2_params(elements[atom_a]).unwrap();
        let se_eff = (pa.selfenergy[sh_a] - pa.kcn[sh_a] * cn[atom_a]) * EV_TO_HARTREE;
        h_mat[(i, i)] = se_eff;
    }
    // Off-diagonal: full tblite H0 with Slater ratio, EN scaling, shell polynomial
    //   hij = 0.5 * (SE_eff_a + SE_eff_b) * S_μν * zij * kpair * kshell * enp * shpoly
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

            // Slater ratio: zij = (2*sqrt(zi*zj)/(zi+zj))^wexp, wexp=0.5
            let zi = pa.slater[sh_a];
            let zj = pb.slater[sh_b];
            let zij = if zi > 1e-10 && zj > 1e-10 {
                let sum = zi + zj;
                (2.0 * (zi * zj).sqrt() / sum).sqrt()
            } else {
                1.0
            };

            // Electronegativity scaling: enp = 1 + 0.02*(EN_A - EN_B)²
            let den = pa.pauling_en - pb.pauling_en;
            let enp = 1.0 + 0.02 * den * den;

            // kshell factor
            let k = gfn2_params::kshell_factor(l_a, l_b);

            // Shell polynomial: rr = sqrt(R / (radA + radB)) where R is in bohr
            let r_bohr = dist_raw(&pos_bohr[atom_a], &pos_bohr[atom_b]);
            let rad_sum = pa.atomic_rad + pb.atomic_rad;
            let rr = if rad_sum > 1e-10 {
                (r_bohr / rad_sum).sqrt()
            } else {
                0.0
            };
            let shpoly = (1.0 + pa.shpoly[sh_a] * rr) * (1.0 + pb.shpoly[sh_b] * rr);

            // H0 off-diagonal
            let se_a = h_mat[(i, i)]; // already in Hartree with CN shift
            let se_b = h_mat[(j, j)];
            let hij = 0.5 * (se_a + se_b) * s_mat[(i, j)] * zij * k * enp * shpoly;
            h_mat[(i, j)] = hij;
            h_mat[(j, i)] = hij;
        }
    }

    // ── DEBUG: dump intermediate values ──
    if debug_diag {
        eprintln!("=== GFN2 DEBUG ===");
        eprintln!("n_atoms={}, n_basis={}, n_electrons={}, n_occ={}", n_atoms, n_basis, n_electrons, n_occ);
        eprintln!("CN: {:?}", cn);
        eprintln!("--- Overlap matrix ---");
        for i in 0..n_basis {
            let row: Vec<f64> = (0..n_basis).map(|j| s_mat[(i,j)]).collect();
            eprintln!("S[{}]: {:?}", i, row);
        }
        eprintln!("--- H0 diagonal (Hartree) ---");
        for i in 0..n_basis {
            let (atom, sh, l, m) = basis_map[i];
            eprintln!("H0[{},{}] = {:.10} (atom={}, sh={}, l={}, m={})", i, i, h_mat[(i,i)], atom, sh, l, m);
        }
        eprintln!("--- H0 off-diagonal (Hartree, nonzero) ---");
        for i in 0..n_basis {
            for j in (i+1)..n_basis {
                if h_mat[(i,j)].abs() > 1e-15 {
                    eprintln!("H0[{},{}] = {:.10}", i, j, h_mat[(i,j)]);
                }
            }
        }
        // H0-only eigenvalues
        let s_eig_dbg = s_mat.clone().symmetric_eigen();
        let mut s_half_inv_dbg = DMatrix::zeros(n_basis, n_basis);
        for k in 0..n_basis {
            let val = s_eig_dbg.eigenvalues[k];
            if val > 1e-8 {
                let inv_sqrt = 1.0 / val.sqrt();
                let col = s_eig_dbg.eigenvectors.column(k);
                for ii in 0..n_basis { for jj in 0..n_basis { s_half_inv_dbg[(ii,jj)] += inv_sqrt * col[ii] * col[jj]; } }
            }
        }
        let f_dbg = &s_half_inv_dbg * &h_mat * &s_half_inv_dbg;
        let eig_dbg = f_dbg.symmetric_eigen();
        let mut h0_evals: Vec<f64> = eig_dbg.eigenvalues.iter().cloned().collect();
        h0_evals.sort_by(|a,b| a.partial_cmp(b).unwrap());
        eprintln!("--- H0 eigenvalues (Hartree) ---");
        for (i, e) in h0_evals.iter().enumerate() {
            eprintln!("  ε_H0[{}] = {:.10} ({:.6} eV)", i, e, e * EV_PER_HARTREE);
        }
        eprintln!("=== END GFN2 DEBUG ===");
    }

    let e_rep = compute_gfn2_repulsion(elements, positions);

    // ── Self-consistent D4 dispersion model (precompute before SCC) ──
    let d4_model = super::d4::D4Model::new(elements, &pos_bohr);

    // ── Shell-resolved SCC setup ──
    struct ShellInfo {
        atom: usize,
        eta: f64,      // shell Hubbard = hubbard * shell_hubbard[sh] (Hartree)
        ref_occ: f64,  // reference occupation
        gam3_sh: f64,  // third-order shell Hubbard derivative (Hartree)
    }
    let mut shells: Vec<ShellInfo> = Vec::new();
    let mut basis_to_shell: Vec<usize> = Vec::with_capacity(n_basis);
    for (i, &z) in elements.iter().enumerate() {
        let p = gfn2_params::get_gfn2_params(z).unwrap();
        for sh in 0..p.nshell {
            let sh_idx = shells.len();
            let l = p.ang_shell[sh];
            // Third-order: gam3_shell = gam3_atom * GAM3_SHELL_FACTOR[kind-1][l]
            let kind_idx = (p.kind as usize).saturating_sub(1);
            let l_idx = (l as usize).min(3);
            let gam3 = p.gam3 * GAM3_SHELL[l_idx][kind_idx];
            shells.push(ShellInfo {
                atom: i,
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


    // Shell-pair gamma matrix (n_shells × n_shells)
    // GFN2 uses the arithmetic-average gamma kernel with gexp=2.0:
    //   γ(A,B) = 1/√(R² + η_avg^{-gexp})^{1/gexp}
    // which for gexp=2 simplifies to the standard Klopman-Ohno form.
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
                // γ = (R^gexp + η_avg^{-gexp})^{-1/gexp}
                let gamma = (r_bohr.powf(gexp) + eta_avg.powf(-gexp)).powf(-1.0 / gexp);
                gamma_sh[i][j] = gamma;
                gamma_sh[j][i] = gamma;
            }
        }
    }

    // ── Löwdin S^{-1/2} ──
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

    // ── SCC loop (all in Hartree, shell-resolved) ──
    // ── AES: Multipole damping radii (CN-dependent) ──
    let mut mrad = vec![0.0f64; n_atoms];
    for iat in 0..n_atoms {
        let p = gfn2_params::get_gfn2_params(elements[iat]).unwrap();
        let arg = cn[iat] - p.mp_vcn - gfn2_params::MP_SHIFT;
        let t1 = (-gfn2_params::MP_KEXP_MULTI * arg).exp();
        let t2 = (gfn2_params::MP_RMAX - p.mp_rad) / (1.0 + t1);
        mrad[iat] = p.mp_rad + t2;
    }

    // ── AES: Atom-pair interaction matrices (0d, finite molecule) ──
    // amat_sd[iat][jat] = [3] (charge-dipole)
    // amat_dd[iat][jat] = [3][3] (dipole-dipole)
    // amat_sq[iat][jat] = [6] (charge-quadrupole)
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

            // charge-dipole: amat_sd[:, jat, iat] += vec * g3 * fdmp3
            for k in 0..3 {
                amat_sd[iat][jat][k] += vec3[k] * g3 * fdmp3;
            }
            // dipole-dipole: amat_dd[:, jat, :, iat] += I*g3*fdmp5 - 3*vec⊗vec*g5*fdmp5
            for a in 0..3 {
                for b in 0..3 {
                    let mut val = -3.0 * vec3[a] * vec3[b] * g5 * fdmp5;
                    if a == b {
                        val += g3 * fdmp5;
                    }
                    amat_dd[iat][jat][a][b] += val;
                }
            }
            // charge-quadrupole: tc components
            // tc = [v1v1, 2*v1v2, v2v2, 2*v1v3, 2*v2v3, v3v3] * g5 * fdmp5
            amat_sq[iat][jat][0] += vec3[0] * vec3[0] * g5 * fdmp5; // xx
            amat_sq[iat][jat][1] += 2.0 * vec3[0] * vec3[1] * g5 * fdmp5; // 2xy
            amat_sq[iat][jat][2] += vec3[1] * vec3[1] * g5 * fdmp5; // yy
            amat_sq[iat][jat][3] += 2.0 * vec3[0] * vec3[2] * g5 * fdmp5; // 2xz
            amat_sq[iat][jat][4] += 2.0 * vec3[1] * vec3[2] * g5 * fdmp5; // 2yz
            amat_sq[iat][jat][5] += vec3[2] * vec3[2] * g5 * fdmp5; // zz
        }
    }

    let max_iter = 250;
    let mut sh_charges = vec![0.0f64; n_shells]; // shell-level Δq
    let mut atom_charges = vec![0.0f64; n_atoms]; // atom-level (for output)
    let mut orbital_energies = vec![0.0f64; n_basis];
    let mut coefficients = DMatrix::zeros(n_basis, n_basis);
    let mut converged = false;
    let mut scc_iter = 0;
    let mut prev_e_elec = 0.0f64;
    let mut sh_shifts = vec![0.0f64; n_shells];
    // AES multipole moments (per atom)
    let mut dpat = vec![[0.0f64; 3]; n_atoms]; // atomic dipole moments
    let mut qpat = vec![[0.0f64; 6]; n_atoms]; // atomic traceless quadrupoles
    // AES multipole potentials (per atom)
    let mut vdp = vec![[0.0f64; 3]; n_atoms]; // dipole potential
    let mut vqp = vec![[0.0f64; 6]; n_atoms]; // quadrupole potential
    let mut vat_mp = vec![0.0f64; n_atoms]; // charge potential from multipoles

    // Broyden mixer (matching tblite exactly)
    // Mix all variables: [sh_charges, dpat_flat, qpat_flat]
    let mix_len = n_shells + n_atoms * 3 + n_atoms * 6;
    let mut mixer = super::broyden::BroydenMixer::new(mix_len, max_iter, 0.4);

    // Helper: flatten dpat/qpat into a contiguous Vec for mixer
    let flatten_multipoles = |sh: &[f64], dp: &[[f64; 3]], qp: &[[f64; 6]]| -> Vec<f64> {
        let mut v = Vec::with_capacity(mix_len);
        v.extend_from_slice(sh);
        for d in dp { v.extend_from_slice(d); }
        for q in qp { v.extend_from_slice(q); }
        v
    };

    for iter in 0..max_iter {
        scc_iter = iter + 1;

        // Broyden mixing: apply mixing step and retrieve mixed charges (for iter > 0)
        if iter > 0 {
            mixer.next().unwrap_or(());
            // Retrieve mixed values back into sh_charges, dpat, qpat
            let mut mixed = vec![0.0f64; mix_len];
            mixer.get(&mut mixed);
            sh_charges.copy_from_slice(&mixed[..n_shells]);
            for iat in 0..n_atoms {
                for k in 0..3 {
                    dpat[iat][k] = mixed[n_shells + iat * 3 + k];
                }
                for c in 0..6 {
                    qpat[iat][c] = mixed[n_shells + n_atoms * 3 + iat * 6 + c];
                }
            }
        }

        // Compute per-shell SCC potential: V_sh = Σ_{sh'} γ(sh,sh') * Δq_{sh'} + Γ_sh * Δq_sh²
        // (second-order + third-order derivative contribution to the Fock matrix)
        for i in 0..n_shells {
            sh_shifts[i] = 0.0;
            for j in 0..n_shells {
                sh_shifts[i] += gamma_sh[i][j] * sh_charges[j];
            }
            // Third-order contribution to the potential: d/dq(1/3 * Γ * Δq³) = Γ * Δq²
            sh_shifts[i] += shells[i].gam3_sh * sh_charges[i] * sh_charges[i];
        }

        // ── AES: Compute multipole potentials from current multipole moments ──
        // Aggregate atom charges from shell charges for AES
        let mut qat = vec![0.0f64; n_atoms];
        for sh in 0..n_shells {
            qat[shells[sh].atom] += sh_charges[sh];
        }
        // vdp = amat_sd @ qat + amat_dd @ dpat  (dipole potential)
        // vqp = amat_sq @ qat                    (quadrupole potential)
        // vat_mp = amat_sd^T @ dpat + amat_sq^T @ qpat  (charge potential from multipoles)
        for iat in 0..n_atoms {
            vdp[iat] = [0.0; 3];
            vqp[iat] = [0.0; 6];
            vat_mp[iat] = 0.0;
        }
        // My convention: amat_sd[iat][jat] = tblite amat_sd(:, jat, iat)
        // tblite: vdp(:, iat) = Σ_jat amat_sd(:, jat, iat) * qat(jat)
        //       = Σ_jat amat_sd[iat][jat][:] * qat[jat]
        // tblite: vat(iat) += Σ_jat amat_sd(:, jat, iat)^T . dpat(:, jat) [trans='T']
        //       = Σ_jat Σ_k amat_sd[iat][jat][k] * dpat[jat][k]
        for iat in 0..n_atoms {
            for jat in 0..n_atoms {
                if iat == jat {
                    continue;
                }
                for k in 0..3 {
                    // vdp[iat] += amat_sd[jat][iat] * qat[jat]  (no-transpose gemv contracts 3rd index)
                    vdp[iat][k] += amat_sd[jat][iat][k] * qat[jat];
                    // vdp[iat] += amat_dd[iat][jat] @ dpat[jat]
                    for l in 0..3 {
                        vdp[iat][k] += amat_dd[iat][jat][k][l] * dpat[jat][l];
                    }
                }
                for c in 0..6 {
                    // vqp[iat] += amat_sq[iat][jat] * qat[jat]
                    vqp[iat][c] += amat_sq[iat][jat][c] * qat[jat];
                }
                // vat_mp: gemv with trans='T' on amat_sd
                // vat(iat) += Σ_jat Σ_k amat_sd[iat][jat][k] * dpat[jat][k]
                for k in 0..3 {
                    vat_mp[iat] += amat_sd[iat][jat][k] * dpat[jat][k];
                }
                // vat(iat) += Σ_jat Σ_c amat_sq[iat][jat][c] * qpat[jat][c]
                for c in 0..6 {
                    vat_mp[iat] += amat_sq[iat][jat][c] * qpat[jat][c];
                }
            }
        }

        // Add kernel potentials: vdp += 2*dkernel*dpat, vqp += 2*qkernel*qpat*mpscale
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

        // ── SC-D4 potential: add dispersion potential to atom-level charge potential ──
        let d4_weights = d4_model.weight_references(&qat);
        let d4_vat = d4_model.get_potential(&d4_weights);
        for iat in 0..n_atoms {
            vat_mp[iat] += d4_vat[iat];
        }

        // Build SCC Fock matrix with AES contribution:
        //   F(i,j) = H0(i,j) - 0.5*S(i,j)*(Vsh(i)+Vsh(j))
        //          - 0.5*S(i,j)*(vat_mp(ai)+vat_mp(aj))
        //          - 0.5*(Σ_k dpint_aboutA(k,i,j)*vdp(k,ai) + dpint_aboutB(k,i,j)*vdp(k,aj))
        //          - 0.5*(Σ_c qpint_aboutA(c,i,j)*vqp(c,ai) + qpint_aboutB(c,i,j)*vqp(c,aj))
        // where dpint[k][(i,j)] = about atom_of_i, dpint[k][(j,i)] = about atom_of_j
        let mut h_scc = DMatrix::zeros(n_basis, n_basis);
        for i in 0..n_basis {
            let shi = basis_to_shell[i];
            let ai = basis_map[i].0;
            // Diagonal: S_μμ = 1, both centers are the same atom
            let mut fii = h_mat[(i, i)] - sh_shifts[shi] - vat_mp[ai];
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
                let mut fij =
                    h_mat[(i, j)] - 0.5 * s_mat[(i, j)] * (sh_shifts[shi] + sh_shifts[shj]);
                // AES charge potential from multipoles
                fij -= 0.5 * s_mat[(i, j)] * (vat_mp[ai] + vat_mp[aj]);
                // AES dipole integral contribution (non-symmetric integrals!)
                // dpint[k][(i,j)] = about ai, dpint[k][(j,i)] = about aj
                for k in 0..3 {
                    fij -= 0.5 * (dpint[k][(i, j)] * vdp[ai][k] + dpint[k][(j, i)] * vdp[aj][k]);
                }
                // AES quadrupole integral contribution (non-symmetric!)
                for c in 0..6 {
                    fij -= 0.5 * (qpint[c][(i, j)] * vqp[ai][c] + qpint[c][(j, i)] * vqp[aj][c]);
                }
                h_scc[(i, j)] = fij;
                h_scc[(j, i)] = fij;
            }
        }

        // Solve HC = SCε via Löwdin orthogonalization
        let f_prime = &s_half_inv * &h_scc * &s_half_inv;
        let eigen = f_prime.symmetric_eigen();

        let mut indices: Vec<usize> = (0..n_basis).collect();
        indices.sort_by(|&a, &b| {
            eigen.eigenvalues[a]
                .partial_cmp(&eigen.eigenvalues[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        for (new_idx, &old_idx) in indices.iter().enumerate() {
            orbital_energies[new_idx] = eigen.eigenvalues[old_idx];
        }

        let c_prime = &eigen.eigenvectors;
        let c_full = &s_half_inv * c_prime;
        for new_idx in 0..n_basis {
            let old_idx = indices[new_idx];
            for i in 0..n_basis {
                coefficients[(i, new_idx)] = c_full[(i, old_idx)];
            }
        }

        // Build density matrix using Fermi occupation numbers
        // nel_per_spin = n_electrons / 2 (restricted case)
        let nel_per_spin = n_electrons as f64 / 2.0;
        let (focc, ts) = get_fermi_filling(&orbital_energies, nel_per_spin, ELECTRONIC_KT);
        let mut density = DMatrix::zeros(n_basis, n_basis);
        for i in 0..n_basis {
            for j in 0..n_basis {
                let mut val = 0.0;
                for k in 0..n_basis {
                    val += focc[k] * coefficients[(i, k)] * coefficients[(j, k)];
                }
                density[(i, j)] = 2.0 * val; // factor 2 for spin
            }
        }

        // Shell-resolved Mulliken charges
        let ps = &density * &s_mat;
        let mut new_sh_charges = vec![0.0f64; n_shells];
        for mu in 0..n_basis {
            let sh = basis_to_shell[mu];
            new_sh_charges[sh] -= ps[(mu, mu)]; // subtract population
        }
        for sh in 0..n_shells {
            new_sh_charges[sh] += shells[sh].ref_occ; // Δq = ref_occ - pop
        }

        // ── AES: Compute Mulliken atomic multipoles from density ──
        // dpat(α, iat) = -Σ_{μ on iat} Σ_ν P(ν,μ) * dpint(α, μ, ν)
        //   where dpint(α, μ, ν) is about atom of μ (= iat)
        // qpat(c, iat) = -Σ_{μ on iat} Σ_ν P(ν,μ) * qpint(c, μ, ν)
        let mut new_dpat = vec![[0.0f64; 3]; n_atoms];
        let mut new_qpat = vec![[0.0f64; 6]; n_atoms];
        for mu in 0..n_basis {
            let iat = basis_map[mu].0;
            for nu in 0..n_basis {
                let p_nu_mu = density[(nu, mu)];
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

        // Electronic energy (Hartree):
        //   E_elec = Tr(P * H0) + E_SCC_shell + E_3rd + E_AES
        let mut e_band = 0.0; // Tr(P * H0)
        for i in 0..n_basis {
            for j in 0..n_basis {
                e_band += density[(i, j)] * h_mat[(i, j)];
            }
        }
        let mut e_scc = 0.0;
        for si in 0..n_shells {
            for sj in 0..n_shells {
                e_scc += gamma_sh[si][sj] * new_sh_charges[si] * new_sh_charges[sj];
            }
        }
        e_scc *= 0.5;
        // Third-order energy: E_3rd = (1/3) * Σ_sh Γ_sh * Δq_sh³
        let mut e_3rd = 0.0;
        for sh in 0..n_shells {
            e_3rd += shells[sh].gam3_sh * new_sh_charges[sh].powi(3);
        }
        e_3rd /= 3.0;

        // AES energy: E_AES = Σ_iat (dpat · vd + qpat · vq) + E_kernel
        // where vd = amat_sd @ qat + 0.5 * amat_dd @ dpat
        // and vq = amat_sq @ qat
        let mut new_qat = vec![0.0f64; n_atoms];
        for sh in 0..n_shells {
            new_qat[shells[sh].atom] += new_sh_charges[sh];
        }
        let mut e_aes = 0.0;
        // Compute vd and vq for energy (without kernel)
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
        // Kernel energy: E_kernel_d = Σ dkernel * (dpat · dpat*mpscale)
        //                E_kernel_q = Σ qkernel * (qpat · qpat*mpscale)
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

        // Split AES energy into components for debugging
        let mut e_aes_dd = 0.0;
        let mut e_aes_sd = 0.0;
        let mut e_aes_sq = 0.0;
        let mut e_aes_dk = 0.0;
        let mut e_aes_qk = 0.0;
        if debug_diag {
            for iat in 0..n_atoms {
                let mut vd_sd = [0.0f64; 3];
                let mut vd_dd = [0.0f64; 3];
                let mut vq_sq = [0.0f64; 6];
                for jat in 0..n_atoms {
                    if iat == jat { continue; }
                    for k in 0..3 {
                        vd_sd[k] += amat_sd[jat][iat][k] * new_qat[jat];
                        for l in 0..3 {
                            vd_dd[k] += 0.5 * amat_dd[iat][jat][k][l] * new_dpat[jat][l];
                        }
                    }
                    for c in 0..6 {
                        vq_sq[c] += amat_sq[iat][jat][c] * new_qat[jat];
                    }
                }
                for k in 0..3 {
                    e_aes_sd += new_dpat[iat][k] * vd_sd[k];
                    e_aes_dd += new_dpat[iat][k] * vd_dd[k];
                }
                for c in 0..6 {
                    e_aes_sq += new_qpat[iat][c] * vq_sq[c];
                }
                let p = gfn2_params::get_gfn2_params(elements[iat]).unwrap();
                let mut dk_e = 0.0;
                for k in 0..3 { dk_e += new_dpat[iat][k] * new_dpat[iat][k]; }
                e_aes_dk += p.dkernel * dk_e;
                let mut qk_e = 0.0;
                for c in 0..6 { qk_e += new_qpat[iat][c] * new_qpat[iat][c] * qp_scale[c]; }
                e_aes_qk += p.qkernel * qk_e;
            }
        }

        let e_elec = e_band + e_scc + e_3rd + e_aes + ts;

        // ── SC-D4 energy: add self-consistent pairwise dispersion to electronic energy ──
        let d4_w_energy = d4_model.weight_references(&new_qat);
        let e_d4_sc = d4_model.get_energy(&d4_w_energy);
        let e_elec = e_elec + e_d4_sc;

        if debug_diag {
            let e_total_iter = e_elec + e_rep;
            eprintln!("SCC iter {}: e_total={:.12} e_band={:.12} e_scc={:.10} e_3rd={:.10} e_aes={:.10} e_d4_sc={:.10} ts={:.10}", 
                      scc_iter, e_total_iter, e_band, e_scc, e_3rd, e_aes, e_d4_sc, ts);
            eprintln!("  AES breakdown: e_sd={:.10} e_dd={:.10} e_sq={:.10} e_dk={:.10} e_qk={:.10}",
                      e_aes_sd, e_aes_dd, e_aes_sq, e_aes_dk, e_aes_qk);
            for iat in 0..n_atoms {
                eprintln!("  atom {}: dpat=[{:.8},{:.8},{:.8}] qpat=[{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}] q={:.12}",
                    iat, new_dpat[iat][0], new_dpat[iat][1], new_dpat[iat][2],
                    new_qpat[iat][0], new_qpat[iat][1], new_qpat[iat][2],
                    new_qpat[iat][3], new_qpat[iat][4], new_qpat[iat][5],
                    new_qat[iat]);
            }
        }

        // Convergence check (matching tblite: econv=1e-6, pconv=2e-5)
        let econv = 1e-6; // Hartree — tblite uses 1e-6 * accuracy (accuracy=1.0)
        let pconv = 2e-5; // tblite uses 2e-5 * accuracy
        // Compute RMS residual for density convergence (same metric as tblite's mixer%get_error)
        let x_in_tmp = flatten_multipoles(&sh_charges, &dpat, &qpat);
        let x_out_tmp = flatten_multipoles(&new_sh_charges, &new_dpat, &new_qpat);
        let rms_dq: f64 = {
            let sum_sq: f64 = x_in_tmp.iter().zip(x_out_tmp.iter())
                .map(|(i, o)| (o - i).powi(2))
                .sum();
            (sum_sq / mix_len as f64).sqrt()
        };
        let energy_conv = (e_elec - prev_e_elec).abs() < econv && iter > 0;
        let charge_conv = rms_dq < pconv;

        if energy_conv && charge_conv {
            converged = true;
            prev_e_elec = e_elec;
            sh_charges = new_sh_charges;
            dpat = new_dpat;
            qpat = new_qpat;
            break;
        }
        prev_e_elec = e_elec;

        // Broyden mixer: set current input, then diff with new output
        // tblite flow: set_mixer(q_in) BEFORE diag, diff_mixer(q_out) AFTER Mulliken
        // Here we do it at end of loop since q_in = sh_charges/dpat/qpat and q_out = new_*
        let x_in = flatten_multipoles(&sh_charges, &dpat, &qpat);
        mixer.set(&x_in);
        let x_out = flatten_multipoles(&new_sh_charges, &new_dpat, &new_qpat);
        mixer.diff(&x_out);
    }

    // Aggregate atom charges from shell charges
    for sh in 0..n_shells {
        atom_charges[shells[sh].atom] += sh_charges[sh];
    }

    let e_elec_ha = prev_e_elec;

    // ── D4 ATM three-body dispersion (non-SC, post-SCF, Hartree) ──
    let disp_ha = d4_model.get_atm_energy(&pos_bohr);

    // ── Halogen bond (Hartree) ──
    let xb_ha = compute_halogen_bond_correction(elements, positions);

    // ── Total energy (Hartree → eV) ──
    let total_ha = e_elec_ha + e_rep + disp_ha + xb_ha;

    // DEBUG: energy decomposition
    let orb_ev: Vec<f64> = orbital_energies.iter().map(|e| e * EV_PER_HARTREE).collect();
    let homo_idx = if n_occ > 0 { n_occ - 1 } else { 0 };
    let lumo_idx = n_occ.min(n_basis - 1);
    let homo_energy = orb_ev[homo_idx];
    let lumo_energy = if n_occ < n_basis {
        orb_ev[lumo_idx]
    } else {
        homo_energy
    };
    let gap = if n_occ < n_basis {
        lumo_energy - homo_energy
    } else {
        0.0
    };

    Ok(Gfn2Result {
        orbital_energies: orb_ev,
        electronic_energy: e_elec_ha * EV_PER_HARTREE,
        repulsive_energy: e_rep * EV_PER_HARTREE,
        dispersion_energy: disp_ha * EV_PER_HARTREE,
        halogen_bond_energy: xb_ha * EV_PER_HARTREE,
        total_energy: total_ha * EV_PER_HARTREE,
        n_basis,
        n_electrons,
        homo_energy,
        lumo_energy,
        gap,
        mulliken_charges: atom_charges,
        atomic_dipoles: dpat,
        atomic_quadrupoles: qpat.iter().map(|q| q[0] + q[2] + q[5]).collect(),
        scc_iterations: scc_iter,
        converged,
    })
}

/// Compute coordination numbers using the double-exponential counting function
/// as used in GFN2-xTB (mctc-lib `dexp` CN type).
///
///   count(R) = exp_count(ka, R, rc) * exp_count(kb, R, rc + r_shift)
///   exp_count(k, R, r0) = 1 / (1 + exp(-k * (r0/R - 1)))
///
/// with ka=10, kb=20, r_shift=2.0 bohr, cutoff=25.0 bohr.
/// Covalent radii: Pyykko & Atsumi 2009 scaled by 4/3 and converted to bohr
/// (same as mctc-lib `covalent_rad_d3`).
fn compute_coordination_numbers(elements: &[u8], pos_bohr: &[[f64; 3]]) -> Vec<f64> {
    let n = elements.len();
    const KA: f64 = 10.0; // steepness of first counting function
    const KB: f64 = 20.0; // steepness of second counting function
    const R_SHIFT: f64 = 2.0; // offset of second counting function (bohr)
    const CUTOFF: f64 = 25.0; // real-space cutoff (bohr)

    #[inline]
    fn exp_count(k: f64, r: f64, r0: f64) -> f64 {
        1.0 / (1.0 + (-k * (r0 / r - 1.0)).exp())
    }

    let mut cn = vec![0.0f64; n];
    for i in 0..n {
        let ri_cov = covalent_radius_d3_bohr(elements[i]);
        for j in (i + 1)..n {
            let rj_cov = covalent_radius_d3_bohr(elements[j]);
            let r_ij = dist_raw(&pos_bohr[i], &pos_bohr[j]);
            if r_ij < 0.1 || r_ij > CUTOFF {
                continue;
            }
            let rc = ri_cov + rj_cov;
            let count = exp_count(KA, r_ij, rc) * exp_count(KB, r_ij, rc + R_SHIFT);
            cn[i] += count;
            cn[j] += count;
        }
    }
    cn
}

/// D3 covalent radii in bohr for CN computation.
/// Pyykko & Atsumi 2009 (metals decreased 10%) scaled by 4/3
/// (matching mctc-lib `covalent_rad_d3` from `covrad.f90`).
fn covalent_radius_d3_bohr(z: u8) -> f64 {
    // Values from mctc-lib covrad.f90 (covalent_rad_2009), NOT atomicrad.f90
    let r_angstrom = match z {
        1 => 0.32,   // H
        2 => 0.46,   // He
        3 => 1.20,   // Li
        4 => 0.94,   // Be
        5 => 0.77,   // B
        6 => 0.75,   // C
        7 => 0.71,   // N
        8 => 0.63,   // O
        9 => 0.64,   // F
        10 => 0.67,  // Ne
        11 => 1.40,  // Na
        12 => 1.25,  // Mg
        13 => 1.13,  // Al
        14 => 1.04,  // Si
        15 => 1.10,  // P
        16 => 1.02,  // S
        17 => 0.99,  // Cl
        18 => 0.96,  // Ar
        22 => 1.22,  // Ti
        24 => 1.10,  // Cr
        25 => 1.07,  // Mn
        26 => 1.04,  // Fe
        27 => 1.00,  // Co
        28 => 0.99,  // Ni
        29 => 1.01,  // Cu
        30 => 1.09,  // Zn
        35 => 1.14,  // Br
        44 => 1.13,  // Ru
        46 => 1.08,  // Pd
        47 => 1.15,  // Ag
        53 => 1.32,  // I
        78 => 1.12,  // Pt
        79 => 1.13,  // Au
        _ => 1.20,
    };
    (4.0 / 3.0) * r_angstrom * ANGSTROM_TO_BOHR
}

/// GFN2-xTB gamma matrix (atom-level, kept for reference).
#[allow(dead_code)]
fn build_gfn2_gamma(elements: &[u8], positions: &[[f64; 3]]) -> Vec<Vec<f64>> {
    let n = elements.len();
    let mut gm = vec![vec![0.0f64; n]; n];
    for a in 0..n {
        let pa = gfn2_params::get_gfn2_params(elements[a]).unwrap();
        // Self-interaction: γ_AA = η_A (Hartree)
        gm[a][a] = pa.hubbard;
        for b in (a + 1)..n {
            let pb = gfn2_params::get_gfn2_params(elements[b]).unwrap();
            let r_bohr = distance_bohr(&positions[a], &positions[b]);
            // γ_AB = 1/√(R² + (1/η_avg)²), η_avg = (η_A + η_B)/2
            let eta_avg = 0.5 * (pa.hubbard + pb.hubbard);
            let gamma = 1.0 / (r_bohr.powi(2) + (1.0 / eta_avg).powi(2)).sqrt();
            gm[a][b] = gamma;
            gm[b][a] = gamma;
        }
    }
    gm
}

/// GFN2-xTB repulsion energy (Hartree).
///
/// V = Σ_{A>B} Z_eff_A × Z_eff_B / R × exp(−(α_A·α_B·R)^(kexp/2))
fn compute_gfn2_repulsion(elements: &[u8], positions: &[[f64; 3]]) -> f64 {
    let n = elements.len();
    let mut e_rep = 0.0;
    for i in 0..n {
        let pi = gfn2_params::get_gfn2_params(elements[i]).unwrap();
        for j in (i + 1)..n {
            let pj = gfn2_params::get_gfn2_params(elements[j]).unwrap();
            let r_bohr = distance_bohr(&positions[i], &positions[j]);
            if r_bohr < 0.2 {
                continue;
            }
            // kexp_light only when BOTH atoms are light (H/He, Z ≤ 2)
            // tblite: merge(kexp, kexp_light, izp > 2 .or. jzp > 2)
            let either_heavy = elements[i] > 2 || elements[j] > 2;
            let kexp = if either_heavy {
                gfn2_params::REP_KEXP
            } else {
                gfn2_params::REP_KEXP_LIGHT
            };
            let rexp = gfn2_params::REP_REXP;
            // Geometric mean of alpha, product of zeff (tblite convention)
            let alpha_pair = (pi.rep_alpha * pj.rep_alpha).sqrt();
            let zeff_pair = pi.rep_zeff * pj.rep_zeff;
            // E = zeff * exp(-alpha * R^kexp) / R^rexp
            let r1k = r_bohr.powf(kexp);
            let exa = (-alpha_pair * r1k).exp();
            let r1r = r_bohr.powf(rexp);
            let e_pair = zeff_pair * exa / r1r;
            e_rep += e_pair;
        }
    }
    e_rep
}

/// D4 dispersion: charge-dependent C6 coefficients (returns Hartree).
fn compute_d4_dispersion(elements: &[u8], positions: &[[f64; 3]], charges: &[f64]) -> f64 {
    let n = elements.len();
    let mut e_disp = 0.0;

    let s6 = gfn2_params::D4_S6;
    let s8 = gfn2_params::D4_S8;
    let a1 = gfn2_params::D4_A1;
    let a2_bohr = gfn2_params::D4_A2 * ANGSTROM_TO_BOHR;

    for i in 0..n {
        for j in (i + 1)..n {
            let r_bohr = distance_bohr(&positions[i], &positions[j]);
            if r_bohr < 0.2 || r_bohr > 95.0 {
                continue;
            }

            let q_scale_i = (-0.08 * charges[i].powi(2)).exp();
            let q_scale_j = (-0.08 * charges[j].powi(2)).exp();

            let c6_base = get_c6_d4(elements[i], elements[j]);
            let c6 = c6_base * q_scale_i * q_scale_j;
            let q_ij = get_r2r4_d4(elements[i]) * get_r2r4_d4(elements[j]);
            let c8 = 3.0 * c6 * q_ij * q_ij;

            let r0 = (c8 / (c6 + 1e-30)).sqrt();
            let f6 = 1.0 / (r_bohr.powi(6) + (a1 * r0 + a2_bohr).powi(6));
            let f8 = 1.0 / (r_bohr.powi(8) + (a1 * r0 + a2_bohr).powi(8));

            e_disp -= s6 * c6 * f6 + s8 * c8 * f8;
        }
    }

    e_disp // Hartree
}

/// Halogen bond correction (returns Hartree).
fn compute_halogen_bond_correction(elements: &[u8], positions: &[[f64; 3]]) -> f64 {
    let n = elements.len();
    let mut e_xb = 0.0;

    let halogens = [17u8, 35, 53];
    let bases = [7u8, 8, 16];

    for i in 0..n {
        if !halogens.contains(&elements[i]) {
            continue;
        }
        let bonded_atom = (0..n)
            .filter(|&k| k != i)
            .min_by(|&a, &b| {
                let da = dist(&positions[i], &positions[a]);
                let db = dist(&positions[i], &positions[b]);
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            });

        for j in 0..n {
            if i == j || !bases.contains(&elements[j]) {
                continue;
            }
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r > 4.0 {
                continue;
            }
            let r0 = get_xb_r0(elements[i], elements[j]);
            let strength = get_xb_strength(elements[i]);

            let cos2_theta = if let Some(r_atom) = bonded_atom {
                let rx = positions[i][0] - positions[r_atom][0];
                let ry = positions[i][1] - positions[r_atom][1];
                let rz = positions[i][2] - positions[r_atom][2];
                let r_ri = (rx * rx + ry * ry + rz * rz).sqrt();
                if r_ri > 0.01 && r > 0.01 {
                    let cos_theta = -(rx * dx + ry * dy + rz * dz) / (r_ri * r);
                    cos_theta.powi(2)
                } else {
                    1.0
                }
            } else {
                1.0
            };

            let x = r / r0;
            let damp = 1.0 / (1.0 + (-20.0 * (x - 1.0)).exp());
            e_xb -= strength * damp * cos2_theta * (-2.0 * (r - r0).powi(2)).exp();
        }
    }

    e_xb // Hartree
}

fn dist(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn get_xb_r0(z_hal: u8, z_base: u8) -> f64 {
    let r_hal = match z_hal {
        17 => 1.75,
        35 => 1.85,
        53 => 1.98,
        _ => 1.80,
    };
    let r_base = match z_base {
        7 => 1.55,
        8 => 1.52,
        16 => 1.80,
        _ => 1.60,
    };
    r_hal + r_base
}

fn get_xb_strength(z_hal: u8) -> f64 {
    match z_hal {
        17 => 0.005,
        35 => 0.010,
        53 => 0.015,
        _ => 0.005,
    }
}

fn get_c6_d4(z1: u8, z2: u8) -> f64 {
    let c6_1 = atomic_c6_d4(z1);
    let c6_2 = atomic_c6_d4(z2);
    (2.0 * c6_1 * c6_2) / (c6_1 + c6_2 + 1e-30)
}

fn atomic_c6_d4(z: u8) -> f64 {
    match z {
        1 => 6.5,
        6 => 46.6,
        7 => 24.2,
        8 => 15.6,
        9 => 9.52,
        14 => 305.0,
        15 => 185.0,
        16 => 134.0,
        17 => 94.6,
        35 => 162.0,
        53 => 408.0,
        22 => 1044.0,
        24 => 602.0,
        25 => 552.0,
        26 => 482.0,
        27 => 408.0,
        28 => 373.0,
        29 => 253.0,
        30 => 284.0,
        _ => 50.0,
    }
}

fn get_r2r4_d4(z: u8) -> f64 {
    match z {
        1 => 2.00,
        6 => 3.09,
        7 => 2.71,
        8 => 2.44,
        9 => 1.91,
        14 => 4.17,
        15 => 3.63,
        16 => 3.49,
        17 => 3.01,
        35 => 3.47,
        53 => 4.38,
        _ => 3.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gfn2_water() {
        let elements = vec![8u8, 1, 1];
        let positions = vec![
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469],
        ];
        let result = solve_gfn2(&elements, &positions);
        assert!(result.is_ok());
        let r = result.unwrap();
        assert!(r.total_energy.is_finite());
        assert!(r.dispersion_energy.is_finite());
        assert!(r.gap > 0.0);
    }

    #[test]
    fn test_gfn2_energy_decomposition() {
        let elements = vec![8u8, 1, 1];
        let positions = vec![
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469],
        ];
        let r = solve_gfn2(&elements, &positions).unwrap();
        eprintln!("=== H2O GFN2 Decomposition ===");
        eprintln!("  E_elec:  {:.6} Ha", r.electronic_energy / 27.2114);
        eprintln!("  E_rep:   {:.6} Ha", r.repulsive_energy / 27.2114);
        eprintln!("  E_disp:  {:.6} Ha", r.dispersion_energy / 27.2114);
        eprintln!("  E_xb:    {:.6} Ha", r.halogen_bond_energy / 27.2114);
        eprintln!("  E_total: {:.6} Ha  (tblite: -5.070370 Ha)", r.total_energy / 27.2114);
        eprintln!("  n_basis={}, n_elec={}, SCC={}, conv={}", r.n_basis, r.n_electrons, r.scc_iterations, r.converged);
        eprintln!("  charges: {:?}", r.mulliken_charges);

        // H2 test
        let h2_r = solve_gfn2(&[1u8, 1], &[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]).unwrap();
        eprintln!("=== H2 GFN2 Decomposition ===");
        eprintln!("  E_elec:  {:.6} Ha", h2_r.electronic_energy / 27.2114);
        eprintln!("  E_rep:   {:.6} Ha", h2_r.repulsive_energy / 27.2114);
        eprintln!("  E_disp:  {:.6} Ha", h2_r.dispersion_energy / 27.2114);
        eprintln!("  E_total: {:.6} Ha  (tblite: -0.981984 Ha)", h2_r.total_energy / 27.2114);
    }

    #[test]
    fn test_gfn2_dispersion() {
        // Test D4 charge-dependent dispersion
        let elements = vec![6u8, 6];
        let positions = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = vec![0.0, 0.0];
        let e_neutral = compute_d4_dispersion(&elements, &positions, &charges);

        let charges_polar = vec![0.5, -0.5];
        let e_polar = compute_d4_dispersion(&elements, &positions, &charges_polar);

        // Charged atoms should have reduced dispersion
        assert!(e_polar.abs() < e_neutral.abs());
    }

    #[test]
    fn test_gfn2_h2_detailed_debug() {
        let elements = vec![1u8, 1];
        let positions = vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let r = solve_gfn2(&elements, &positions).unwrap();
        let e_ha = r.total_energy / EV_PER_HARTREE;
        eprintln!("=== H2 Detailed ===");
        eprintln!("  E_elec:  {:.10} Ha", r.electronic_energy / EV_PER_HARTREE);
        eprintln!("  E_rep:   {:.10} Ha", r.repulsive_energy / EV_PER_HARTREE);
        eprintln!("  E_disp:  {:.10} Ha", r.dispersion_energy / EV_PER_HARTREE);
        eprintln!("  E_total: {:.10} Ha  (tblite: -0.9819836925 Ha)", e_ha);
        eprintln!("  error:   {:.6}%", (e_ha - (-0.9819836925)).abs() / 0.9819836925 * 100.0);

        // Now manually compute overlap and H0 for H2 to check
        use crate::xtb::solver::ANGSTROM_TO_BOHR;
        let pos_bohr: Vec<[f64; 3]> = positions.iter().map(|p| [p[0] * ANGSTROM_TO_BOHR, p[1] * ANGSTROM_TO_BOHR, p[2] * ANGSTROM_TO_BOHR]).collect();
        let r_bohr = dist_raw(&pos_bohr[0], &pos_bohr[1]);
        eprintln!("  R(H-H): {:.10} bohr ({:.6} Å)", r_bohr, 0.74);

        let p = gfn2_params::get_gfn2_params(1).unwrap();
        eprintln!("  H slater: {:.6}, pqn: {}, ngauss: {}", p.slater[0], p.pqn[0], p.ngauss[0]);
        eprintln!("  H selfenergy: {:.6} eV = {:.10} Ha", p.selfenergy[0], p.selfenergy[0] * EV_TO_HARTREE);
        eprintln!("  H kcn: {:.6} eV", p.kcn[0]);
        eprintln!("  H pauling_en: {:.4}", p.pauling_en);
        eprintln!("  H shpoly[0]: {:.10}", p.shpoly[0]);
        eprintln!("  H atomic_rad: {:.10} bohr", p.atomic_rad);
        eprintln!("  H hubbard: {:.10} Ha", p.hubbard);
        eprintln!("  H shell_hubbard[0]: {:.10}", p.shell_hubbard[0]);

        // Overlap S12
        let s12 = crate::xtb::sto_overlap::sto_overlap_with_ng(
            p.pqn[0], 0, 0, p.slater[0], &pos_bohr[0], p.ngauss[0],
            p.pqn[0], 0, 0, p.slater[0], &pos_bohr[1], p.ngauss[0],
        );
        eprintln!("  S12 overlap: {:.10}", s12);

        // zij factor
        let zij = (2.0_f64 * (p.slater[0] * p.slater[0]).sqrt() / (p.slater[0] + p.slater[0])).sqrt();
        eprintln!("  zij: {:.10} (should be 1.0 for same element)", zij);

        // enp factor
        let den = p.pauling_en - p.pauling_en;
        let enp = 1.0 + 0.02 * den * den;
        eprintln!("  enp: {:.10} (should be 1.0 for same element)", enp);

        // kshell
        let k = gfn2_params::kshell_factor(0, 0);
        eprintln!("  kshell(s,s): {:.10}", k);

        // Compute CN
        let cn = compute_coordination_numbers(&elements, &pos_bohr);
        eprintln!("  CN: [{:.10}, {:.10}]", cn[0], cn[1]);

        // SE_eff with CN shift
        let se_eff = (p.selfenergy[0] - p.kcn[0] * cn[0]) * EV_TO_HARTREE;
        eprintln!("  SE_eff: {:.10} Ha = {:.6} eV", se_eff, se_eff / EV_TO_HARTREE);

        // shpoly
        let rad_sum = p.atomic_rad + p.atomic_rad;
        let rr = (r_bohr / rad_sum).sqrt();
        let shpoly = (1.0 + p.shpoly[0] * rr) * (1.0 + p.shpoly[0] * rr);
        eprintln!("  rr: {:.10}", rr);
        eprintln!("  shpoly: {:.10}", shpoly);

        // H0_12 = 0.5 * (SE_eff_1 + SE_eff_2) * S12 * zij * k * enp * shpoly
        let h12 = 0.5 * (se_eff + se_eff) * s12 * zij * k * enp * shpoly;
        eprintln!("  H0_12: {:.10} Ha", h12);
        eprintln!("  H0_11 = SE_eff = {:.10} Ha", se_eff);

        // Exact eigenvalues of 2x2 system
        let h11 = se_eff;
        let s_diag = 1.0;
        eprintln!("  --- 2x2 generalized eigenvalue ---");
        eprintln!("  H = [[{:.8}, {:.8}], [{:.8}, {:.8}]]", h11, h12, h12, h11);
        eprintln!("  S = [[{:.8}, {:.8}], [{:.8}, {:.8}]]", s_diag, s12, s12, s_diag);

        // Eigenvalues of H*c = e*S*c for 2x2 symmetric:
        // e_bond = (H11 + H12) / (1 + S12)
        // e_anti = (H11 - H12) / (1 - S12)
        let e_bond = (h11 + h12) / (1.0 + s12);
        let e_anti = (h11 - h12) / (1.0 - s12);
        eprintln!("  E_bond = {:.10} Ha, E_anti = {:.10} Ha", e_bond, e_anti);
        eprintln!("  E_elec(H0) = 2 * E_bond = {:.10} Ha", 2.0 * e_bond);
    }

    #[test]
    fn test_gfn2_acetylene_debug() {
        // C2H2 (acetylene): H-C≡C-H along z-axis
        // C≡C = 1.203 Å, C-H = 1.060 Å
        let elements = vec![6u8, 6, 1, 1];
        let positions = vec![
            [0.0, 0.0, -1.203 / 2.0],
            [0.0, 0.0, 1.203 / 2.0],
            [0.0, 0.0, -(1.203 / 2.0 + 1.060)],
            [0.0, 0.0, 1.203 / 2.0 + 1.060],
        ];
        let r = solve_gfn2(&elements, &positions).unwrap();
        let e_ha = r.total_energy / EV_PER_HARTREE;
        eprintln!("=== C2H2 (Acetylene) GFN2 ===");
        eprintln!("  E_elec:  {:.10} Ha", r.electronic_energy / EV_PER_HARTREE);
        eprintln!("  E_rep:   {:.10} Ha", r.repulsive_energy / EV_PER_HARTREE);
        eprintln!("  E_disp:  {:.10} Ha", r.dispersion_energy / EV_PER_HARTREE);
        eprintln!("  E_total: {:.10} Ha  (tblite: -5.2064375858 Ha)", e_ha);
        eprintln!("  error:   {:.6}%", (e_ha - (-5.2064375858)).abs() / 5.2064375858 * 100.0);
        eprintln!("  n_basis={}, n_elec={}, SCC={}, conv={}", r.n_basis, r.n_electrons, r.scc_iterations, r.converged);
        eprintln!("  charges: {:?}", r.mulliken_charges);
        eprintln!("  HOMO: {:.6} eV, LUMO: {:.6} eV, gap: {:.6} eV", r.homo_energy, r.lumo_energy, r.gap);
    }
}
