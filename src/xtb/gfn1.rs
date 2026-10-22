//! GFN1-xTB: faithful port of tblite's native `GFN1-xTB` calculator.
//!
//! Term-by-term port of tblite v0.7.0: Slater-basis H0 with the GFN1
//! kpair/kshell/enscale enhancement, shell-resolved second-order
//! electrostatics (harmonic-averaged Hubbard hardness, Klopman-Ohno kernel
//! with `gexp = 2`), atom-resolved onsite third order, effective screened
//! repulsion, D3(BJ) dispersion with the GFN1 damping parameters and the
//! classical halogen-bond correction.
//!
//! Units: Hartree/bohr internally, eV in the result. Reference formulas:
//! `h0.f90` (Hamiltonian), `coulomb/charge/{effective,type}.f90` +
//! `coulomb/thirdorder.f90` (electrostatics), `repulsion/effective.f90`
//! (repulsion), `disp/d3.f90` (dispersion), `classical/halogen.f90`
//! (halogen bond), `basis/ortho.f90` (Schmidt orthogonalization) and
//! `scf/{iterator,potential}.f90` (SCF driver). Parameter source:
//! `src/tblite/xtb/gfn1.f90` (LGPL-3.0-or-later), parsed into
//! `tests/fixtures/validation/gfn1_params.json` and generated into
//! [`super::gfn1_params`] by `scripts/validation/{parse,gen}_gfn1_params.py`.
//!
//! Cross-validation: F44 in `tests/cross_validation.rs` against tblite.

use super::broyden::BroydenMixer;
use super::gfn0::{BOHR_PER_ANGSTROM, EV_PER_HARTREE, GFN0_ATOMIC_RAD};
use super::gfn1_params::{get_gfn1_params, GFN1_GLOBALS, GFN1_NVALENCE};
use super::gfn2::get_fermi_filling;
use super::sto_overlap::{cart_prim_overlap, sto_table_prims};
use crate::hf::d3::{compute_d3_energy_with_params, coordination_numbers_exp, D3Params};
use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

/// Electronic temperature in Hartree (tblite default, ≈300 K).
const ELECTRONIC_KT: f64 = 9.500e-4;

/// GFN1-xTB calculation result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gfn1Result {
    /// Orbital energies (eV).
    pub orbital_energies: Vec<f64>,
    /// Electronic energy (eV).
    pub electronic_energy: f64,
    /// Repulsive energy (eV).
    pub repulsive_energy: f64,
    /// Dispersion energy (eV) — D3-BJ total.
    pub dispersion_energy: f64,
    /// Two-body C6 contribution of the D3-BJ dispersion (eV).
    pub dispersion_e6: f64,
    /// Two-body C8 contribution of the D3-BJ dispersion (eV).
    pub dispersion_e8: f64,
    /// Halogen-bond correction energy (eV).
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
    /// Shell-resolved charges (Δq per atom, indexed s/p/d).
    pub shell_charges: Vec<Vec<f64>>,
    /// SCC iterations.
    pub scc_iterations: usize,
    /// Whether the SCC converged.
    pub converged: bool,
}

/// kshell matrix from tblite's rule: 2.08 for s-p pairs,
/// (kdiag[k]+kdiag[l])/2 otherwise with kdiag = [1.85, 2.25, 2, 2, 2].
fn gfn1_kshell(k: u8, l: u8) -> f64 {
    if (k == 0 && l == 1) || (k == 1 && l == 0) {
        2.08
    } else {
        let kd = GFN1_GLOBALS.kdiag;
        0.5 * (kd[k as usize] + kd[l as usize])
    }
}

/// Transition-metal d-block row (tblite `get_dblock_row`): 1 (3d), 2 (4d),
/// 3 (5d), 0 otherwise.
fn dblock_row(z: u8) -> u8 {
    match z {
        21..=29 => 1,
        39..=47 => 2,
        57..=79 => 3,
        _ => 0,
    }
}

/// Pair parameter (tblite `get_pair_param`): special H/metal pairs,
/// d-block average, default 1.
fn gfn1_kpair(zi: u8, zj: u8) -> f64 {
    if zi == 1 && zj == 1 {
        return 0.96;
    }
    for (x, y, v) in [
        (5u8, 1u8, 0.95),
        (7, 1, 1.04),
        (28, 1, 0.90),
        (75, 1, 0.80),
        (78, 1, 0.80),
        (15, 5, 0.97),
        (14, 7, 1.01),
    ] {
        if (zi == x && zj == y) || (zi == y && zj == x) {
            return v;
        }
    }
    let kp = [1.1, 1.2, 1.2];
    let (ti, tj) = (dblock_row(zi), dblock_row(zj));
    if ti > 0 && tj > 0 {
        0.5 * (kp[ti as usize - 1] + kp[tj as usize - 1])
    } else {
        1.0
    }
}

/// First shell of each angular momentum is valence (tblite valence rule).
pub(crate) fn gfn1_valence(ang: &[u8]) -> Vec<bool> {
    let mut seen = [false; 5];
    ang.iter()
        .map(|&l| {
            if seen[l as usize] {
                false
            } else {
                seen[l as usize] = true;
                true
            }
        })
        .collect()
}

/// Number of SAOs per shell (s=1, p=3, d=5).
pub(crate) fn nsao(l: u8) -> usize {
    match l {
        0 => 1,
        1 => 3,
        2 => 5,
        _ => panic!("GFN1: unsupported shell angular momentum"),
    }
}

/// One contracted shell with absolute primitives.
pub(crate) struct Gfn1Shell {
    l: u8,
    prims: Vec<(f64, f64)>,
}

/// tblite's `slater_to_gauss(..., norm=.true.)` primitive normalization:
/// `(2α/π)^{3/4} (4α)^{l/2} / sqrt((2l-1)!!)`.
fn slater_prim_norm(alpha: f64, l: u8) -> f64 {
    use std::f64::consts::PI;
    let dfact: f64 = match l {
        0 | 1 => 1.0,
        2 => 3.0,
        _ => panic!("GFN1: unsupported shell angular momentum"),
    };
    (2.0 * alpha / PI).powf(0.75) * (4.0 * alpha).powf(l as f64 / 2.0) / dfact.sqrt()
}

/// Build shell contractions for an atom (tblite `add_basis` / `orthogonalize`).
///
/// Repeated-l shells are Schmidt-orthogonalized against the earlier same-l
/// shell exactly as tblite's `orthogonalize()` does: in the normalized
/// primitive-coefficient convention of `slater_to_gauss` with the bare
/// `(π/(a+b))^{3/2}` kernel, followed by a bare renormalization. The result is
/// converted back to bare coefficients (our overlap machinery applies the
/// per-primitive normalization itself).
fn gfn1_atom_shells(z: u8) -> Vec<Gfn1Shell> {
    use std::f64::consts::PI;
    let p = get_gfn1_params(z).expect("GFN1 params");
    let mut seen: Vec<(u8, Vec<(f64, f64)>)> = Vec::new();
    let mut out = Vec::new();
    for sh in 0..p.nshell {
        let l = p.ang_shell[sh];
        let zeta = p.slater[sh];
        let table = sto_table_prims(p.pqn[sh], l, p.ngauss[sh]);
        let mut cprim: Vec<(f64, f64)> = table
            .iter()
            .map(|(ar, ca)| {
                let a = ar * zeta * zeta;
                (a, ca * slater_prim_norm(a, l))
            })
            .collect();
        if let Some((_, ref prev)) = seen.iter().find(|(ll, _)| *ll == l) {
            // tblite orthogonalize(): bare-kernel Schmidt + bare renormalization.
            let mut ov = 0.0;
            for (ar, ca) in cprim.iter() {
                for (br, cb) in prev.iter() {
                    ov += ca * cb * (PI / (ar + br)).powf(1.5);
                }
            }
            let mut combined = cprim.clone();
            for (br, cb) in prev.iter() {
                combined.push((*br, -ov * cb));
            }
            let mut norm2 = 0.0;
            for (ar, ca) in combined.iter() {
                for (br, cb) in combined.iter() {
                    norm2 += ca * cb * (PI / (ar + br)).powf(1.5);
                }
            }
            let nrm = norm2.sqrt();
            for (_, c) in combined.iter_mut() {
                *c /= nrm;
            }
            cprim = combined;
        } else {
            seen.push((l, cprim.clone()));
        }
        let prims: Vec<(f64, f64)> = cprim
            .iter()
            .map(|&(a, c)| (a, c / slater_prim_norm(a, l)))
            .collect();
        out.push(Gfn1Shell { l, prims });
    }
    out
}

/// Cartesian component for an s/p SAO index.
fn ca_comps(l: u8, m: usize) -> ((u8, u8, u8), f64) {
    match l {
        0 => ((0, 0, 0), 1.0),
        1 => match m {
            0 => ((1, 0, 0), 1.0),
            1 => ((0, 1, 0), 1.0),
            _ => ((0, 0, 1), 1.0),
        },
        _ => panic!("ca_comps: d needs sao_pair"),
    }
}

/// Real spherical d orbital as a Cartesian combination (rotation-invariant
/// energies make the SAO convention irrelevant as long as S and H share it):
/// 0=dxy, 1=dxz, 2=dyz, 3=dx2-y2, 4=dz2.
fn sph_d_comp(m: usize) -> Vec<((u8, u8, u8), f64)> {
    let s3_4 = (3.0f64 / 4.0).sqrt();
    match m {
        0 => vec![((1, 1, 0), 1.0)],
        1 => vec![((1, 0, 1), 1.0)],
        2 => vec![((0, 1, 1), 1.0)],
        3 => vec![((2, 0, 0), s3_4), ((0, 2, 0), -s3_4)],
        _ => vec![((2, 0, 0), -0.5), ((0, 2, 0), -0.5), ((0, 0, 2), 1.0)],
    }
}

/// Cartesian expansion lists for an SAO pair.
type SaoExpansion = Vec<((u8, u8, u8), f64)>;

fn sao_pair(l_a: u8, l_b: u8, ii: usize, jj: usize) -> (SaoExpansion, SaoExpansion) {
    let la = if l_a == 2 {
        sph_d_comp(ii)
    } else {
        vec![ca_comps(l_a, ii)]
    };
    let lb = if l_b == 2 {
        sph_d_comp(jj)
    } else {
        vec![ca_comps(l_b, jj)]
    };
    (la, lb)
}

/// Basis bookkeeping for the GFN1 SAO basis.
pub(crate) struct Gfn1Basis {
    /// Contracted shells per atom.
    pub(crate) atom_shells: Vec<Vec<Gfn1Shell>>,
    /// (atom, shell, sao index) per basis function.
    pub(crate) funcs: Vec<(usize, usize, usize)>,
    /// First basis-function index of each (atom, shell).
    pub(crate) shell_ao: Vec<Vec<usize>>,
    /// Number of basis functions.
    pub(crate) n_basis: usize,
}

pub(crate) fn build_gfn1_basis(elements: &[u8]) -> Gfn1Basis {
    let atom_shells: Vec<Vec<Gfn1Shell>> = elements.iter().map(|&z| gfn1_atom_shells(z)).collect();
    let mut funcs: Vec<(usize, usize, usize)> = Vec::new();
    let mut shell_ao: Vec<Vec<usize>> = Vec::new();
    for (at, shells) in atom_shells.iter().enumerate() {
        let mut offs = Vec::new();
        for (sh, shell) in shells.iter().enumerate() {
            offs.push(funcs.len());
            for m in 0..nsao(shell.l) {
                funcs.push((at, sh, m));
            }
        }
        shell_ao.push(offs);
    }
    let n_basis = funcs.len();
    Gfn1Basis {
        atom_shells,
        funcs,
        shell_ao,
        n_basis,
    }
}

/// SAO overlap between shell `sh_a` on atom `a` and shell `sh_b` on atom `b`.
fn shell_overlap(
    basis: &Gfn1Basis,
    pos_bohr: &[[f64; 3]],
    a: usize,
    sh_a: usize,
    b: usize,
    sh_b: usize,
) -> DMatrix<f64> {
    let ca = &basis.atom_shells[a][sh_a];
    let cb = &basis.atom_shells[b][sh_b];
    let ma = nsao(ca.l);
    let mb = nsao(cb.l);
    let mut block = DMatrix::zeros(ma, mb);
    for ii in 0..ma {
        for jj in 0..mb {
            let (la_list, lb_list) = sao_pair(ca.l, cb.l, ii, jj);
            let mut sij = 0.0;
            for ((lxa, lya, lza), wa) in &la_list {
                for ((lxb, lyb, lzb), wb) in &lb_list {
                    sij += wa
                        * wb
                        * cart_prim_overlap(
                            *lxa,
                            *lya,
                            *lza,
                            &ca.prims,
                            &pos_bohr[a],
                            *lxb,
                            *lyb,
                            *lzb,
                            &cb.prims,
                            &pos_bohr[b],
                        );
                }
            }
            block[(ii, jj)] = sij;
        }
    }
    block
}

/// GFN1 overlap and H0 matrices (Hartree), following tblite `get_hamiltonian`.
pub(crate) fn build_overlap_hamiltonian(
    basis: &Gfn1Basis,
    elements: &[u8],
    pos_bohr: &[[f64; 3]],
    cn: &[f64],
) -> (DMatrix<f64>, DMatrix<f64>) {
    let n = basis.n_basis;
    let nat = elements.len();
    let mut s_mat = DMatrix::zeros(n, n);
    let mut h_mat = DMatrix::zeros(n, n);

    // CN-shifted self energies (Hartree).
    let se: Vec<Vec<f64>> = elements
        .iter()
        .enumerate()
        .map(|(i, &z)| {
            let p = get_gfn1_params(z).expect("GFN1 params");
            (0..p.nshell)
                .map(|sh| p.selfenergy[sh] - p.kcn[sh] * cn[i])
                .collect()
        })
        .collect();
    let valence: Vec<Vec<bool>> = elements
        .iter()
        .map(|&z| {
            let p = get_gfn1_params(z).expect("GFN1 params");
            gfn1_valence(&p.ang_shell[..p.nshell])
        })
        .collect();

    // Helper to stamp a shell block (diagonal blocks on the same atom,
    // off-diagonal blocks between two atoms) into S and H.
    let stamp = |s_mat: &mut DMatrix<f64>,
                 h_mat: &mut DMatrix<f64>,
                 a: usize,
                 sh_a: usize,
                 b: usize,
                 sh_b: usize,
                 block_s: DMatrix<f64>,
                 hij: f64| {
        let ia = basis.shell_ao[a][sh_a];
        let ib = basis.shell_ao[b][sh_b];
        let (ma, mb) = (block_s.nrows(), block_s.ncols());
        for ii in 0..ma {
            for jj in 0..mb {
                let v = block_s[(ii, jj)];
                s_mat[(ia + ii, ib + jj)] = v;
                s_mat[(ib + jj, ia + ii)] = v;
                let h = hij * v;
                h_mat[(ia + ii, ib + jj)] = h;
                h_mat[(ib + jj, ia + ii)] = h;
            }
        }
    };

    // Onsite blocks: overlap x 0.5*(se_i + se_j) with no shpoly (rr = 0).
    for a in 0..nat {
        let nsh = basis.atom_shells[a].len();
        for sh_a in 0..nsh {
            for sh_b in 0..nsh {
                let hij = 0.5 * (se[a][sh_a] + se[a][sh_b]);
                let block = if sh_a == sh_b {
                    let mut id = DMatrix::zeros(
                        nsao(basis.atom_shells[a][sh_a].l),
                        nsao(basis.atom_shells[a][sh_b].l),
                    );
                    for i in 0..id.nrows() {
                        id[(i, i)] = 1.0;
                    }
                    id
                } else {
                    shell_overlap(basis, pos_bohr, a, sh_a, a, sh_b)
                };
                stamp(&mut s_mat, &mut h_mat, a, sh_a, a, sh_b, block, hij);
            }
        }
    }

    // Interatomic blocks with the full GFN1 enhancement factor.
    for a in 0..nat {
        for b in (a + 1)..nat {
            let pa = get_gfn1_params(elements[a]).expect("GFN1 params");
            let pb = get_gfn1_params(elements[b]).expect("GFN1 params");
            let den = crate::xtb::gfn0::pauling_en(elements[a])
                - crate::xtb::gfn0::pauling_en(elements[b]);
            let den2 = den * den;
            let dx = pos_bohr[a][0] - pos_bohr[b][0];
            let dy = pos_bohr[a][1] - pos_bohr[b][1];
            let dz = pos_bohr[a][2] - pos_bohr[b][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            let rad_sum = GFN0_ATOMIC_RAD[elements[a] as usize - 1]
                + GFN0_ATOMIC_RAD[elements[b] as usize - 1];
            let rr = (r / rad_sum).sqrt();
            for sh_a in 0..basis.atom_shells[a].len() {
                let l_a = basis.atom_shells[a][sh_a].l;
                for sh_b in 0..basis.atom_shells[b].len() {
                    let l_b = basis.atom_shells[b][sh_b].l;
                    let km = match (valence[a][sh_a], valence[b][sh_b]) {
                        (true, true) => {
                            gfn1_kpair(elements[a], elements[b])
                                * gfn1_kshell(l_b, l_a)
                                * (1.0 + GFN1_GLOBALS.enscale * den2)
                        }
                        (true, false) => 0.5 * (gfn1_kshell(l_a, l_a) + GFN1_GLOBALS.kdiff),
                        (false, true) => 0.5 * (gfn1_kshell(l_b, l_b) + GFN1_GLOBALS.kdiff),
                        (false, false) => GFN1_GLOBALS.kdiff,
                    };
                    let shpoly =
                        (1.0 + pa.shpoly[l_a as usize] * rr) * (1.0 + pb.shpoly[l_b as usize] * rr);
                    let hij = 0.5 * (se[a][sh_a] + se[b][sh_b]) * km * shpoly;
                    let block = shell_overlap(basis, pos_bohr, a, sh_a, b, sh_b);
                    stamp(&mut s_mat, &mut h_mat, a, sh_a, b, sh_b, block, hij);
                }
            }
        }
    }

    (s_mat, h_mat)
}

/// GFN1 effective screened repulsion (Hartree), tblite `repulsion/effective.f90`.
pub(crate) fn compute_repulsion(elements: &[u8], pos_bohr: &[[f64; 3]]) -> f64 {
    let n = elements.len();
    let kexp = GFN1_GLOBALS.rep_kexp;
    let rexp = GFN1_GLOBALS.rep_rexp;
    let cutoff2 = 25.0f64.powi(2);
    let mut e_rep = 0.0;
    for i in 0..n {
        let pi = get_gfn1_params(elements[i]).expect("GFN1 params");
        for j in (i + 1)..n {
            let pj = get_gfn1_params(elements[j]).expect("GFN1 params");
            let dx = pos_bohr[i][0] - pos_bohr[j][0];
            let dy = pos_bohr[i][1] - pos_bohr[j][1];
            let dz = pos_bohr[i][2] - pos_bohr[j][2];
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 > cutoff2 || r2 < 1e-12 {
                continue;
            }
            let r = r2.sqrt();
            let alpha = (pi.rep_alpha * pj.rep_alpha).sqrt();
            let zeff = pi.rep_zeff * pj.rep_zeff;
            e_rep += zeff * (-alpha * r.powf(kexp)).exp() / r.powf(rexp);
        }
    }
    e_rep
}

/// Classical halogen-bond correction (Hartree), tblite `classical/halogen.f90`.
pub(crate) fn compute_halogen_bond(elements: &[u8], pos_bohr: &[[f64; 3]]) -> f64 {
    const HALOGENS: [u8; 4] = [17, 35, 53, 85];
    const ACCEPTORS: [u8; 4] = [7, 8, 15, 16];
    const CUTOFF: f64 = 20.0;
    const ALP: f64 = 6.0;
    const LJ: f64 = 12.0;

    let n = elements.len();
    if !elements.iter().any(|z| HALOGENS.contains(z))
        || !elements.iter().any(|z| ACCEPTORS.contains(z))
    {
        return 0.0;
    }
    let damping = GFN1_GLOBALS.halogen_damping;
    let rad_scale = GFN1_GLOBALS.halogen_radscale;
    let rad = |z: u8| GFN0_ATOMIC_RAD[(z - 1) as usize] * rad_scale;
    let dist2 = |a: usize, b: usize| {
        let dx = pos_bohr[a][0] - pos_bohr[b][0];
        let dy = pos_bohr[a][1] - pos_bohr[b][1];
        let dz = pos_bohr[a][2] - pos_bohr[b][2];
        dx * dx + dy * dy + dz * dz
    };

    let mut e_xb = 0.0;
    for x in 0..n {
        if !HALOGENS.contains(&elements[x]) {
            continue;
        }
        let cc = get_gfn1_params(elements[x])
            .expect("GFN1 params")
            .halogen_bond;
        // Nearest neighbour of the halogen (tblite `get_xbond_list`).
        let mut nearest: Option<(usize, f64)> = None;
        for k in 0..n {
            if k == x {
                continue;
            }
            let r = dist2(x, k).sqrt();
            if r > 0.0 && nearest.is_none_or(|(_, d)| r < d) {
                nearest = Some((k, r));
            }
        }
        let Some((k, rkx)) = nearest else { continue };
        let d2kx = rkx * rkx;
        for j in 0..n {
            if !ACCEPTORS.contains(&elements[j]) {
                continue;
            }
            let d2jx = dist2(x, j);
            let d2jk = dist2(j, k);
            let rjx = d2jx.sqrt();
            if !(1e-18..=CUTOFF).contains(&rjx) {
                continue;
            }
            let term = (d2kx + d2jx - d2jk) / (rkx * rjx);
            let aterm = (0.5 - 0.25 * term).powf(ALP);
            let r0jx = rad(elements[x]) + rad(elements[j]);
            let t13 = r0jx / rjx;
            let t14 = t13.powf(LJ);
            e_xb += aterm * cc * (t14 - damping * t13.powf(LJ * 0.5)) / (1.0 + t14);
        }
    }
    e_xb
}

/// Solve GFN1-xTB for a molecular system (neutral, closed shell).
///
/// `positions` are Cartesian coordinates in Ångström. The returned
/// [`Gfn1Result`] carries energies in eV.
///
/// # Errors
///
/// Returns `Err` if the geometry is invalid, if `elements` and `positions`
/// have different lengths, or if an element has no GFN1 parameters.
///
/// # Panics
///
/// Panics if an element has no GFN1 parameters. This cannot happen for the
/// molecules this function accepts, because the parameter lookup is
/// validated before the calculation starts.
pub fn solve_gfn1(elements: &[u8], positions: &[[f64; 3]]) -> Result<Gfn1Result, String> {
    crate::validation::validate_geometry(elements, positions)?;
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
    let n_occ = n_electrons / 2;

    // Coordination numbers: tblite `cn_count%exp` (kcn = 16, D3 covalent radii).
    let cn = coordination_numbers_exp(elements, &pos_bohr);

    // Overlap and H0 (Hartree).
    let (s_mat, h_mat) = build_overlap_hamiltonian(&basis, elements, &pos_bohr, &cn);

    // ── Shell bookkeeping ──
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

    // Atom-resolved Hubbard derivatives (third order is atom resolved in GFN1).
    let hd: Vec<f64> = elements
        .iter()
        .map(|&z| get_gfn1_params(z).expect("GFN1 params").hubbard_derivs)
        .collect();

    // ── Second-order electrostatics kernel ──
    // gamma_ij = (R_ij^gexp + gam^-gexp)^(-1/gexp) with harmonic-averaged
    // shell hardness `gam`; same-atom pairs use `gam` directly.
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

    // ── SCC loop ──
    let max_iter = 250;
    let mut qsh = vec![0.0f64; n_shells];
    let mut mixer = BroydenMixer::new(n_shells, max_iter, 0.4);
    let mut orbital_energies = vec![0.0f64; n_basis];
    let mut coefficients = DMatrix::zeros(n_basis, n_basis);
    let mut converged = false;
    let mut scc_iter = 0;
    let mut e_elec = 0.0f64;
    let mut prev_e = 0.0f64;

    for iter in 0..max_iter {
        scc_iter = iter + 1;

        if iter > 0 {
            mixer.step().unwrap_or(());
            let mut mixed = vec![0.0f64; n_shells];
            mixer.get(&mut mixed);
            qsh.copy_from_slice(&mixed);
        }

        // Atom charges and potential shifts.
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

        // F = H0 - 0.5*S*(v_i + v_j)
        let mut f_mat = DMatrix::zeros(n_basis, n_basis);
        for i in 0..n_basis {
            let vi = vsh[basis_to_shell[i]] + vat[basis.funcs[i].0];
            for j in i..n_basis {
                let vj = vsh[basis_to_shell[j]] + vat[basis.funcs[j].0];
                let fij = h_mat[(i, j)] - 0.5 * s_mat[(i, j)] * (vi + vj);
                f_mat[(i, j)] = fij;
                if i != j {
                    f_mat[(j, i)] = fij;
                }
            }
        }

        // Solve FC = SCε via Löwdin orthogonalization.
        let f_prime = &s_half_inv * &f_mat * &s_half_inv;
        let eigen = f_prime.symmetric_eigen();
        let mut indices: Vec<usize> = (0..n_basis).collect();
        indices.sort_by(|&a, &b| {
            eigen.eigenvalues[a]
                .partial_cmp(&eigen.eigenvalues[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let c_full = &s_half_inv * &eigen.eigenvectors;
        for (new_idx, &old_idx) in indices.iter().enumerate() {
            orbital_energies[new_idx] = eigen.eigenvalues[old_idx];
            for i in 0..n_basis {
                coefficients[(i, new_idx)] = c_full[(i, old_idx)];
            }
        }

        // Fermi filling (per spin channel) and density matrix.
        let nel_per_spin = n_electrons as f64 / 2.0;
        let (focc, ts) = get_fermi_filling(&orbital_energies, nel_per_spin, ELECTRONIC_KT);
        let mut cw = coefficients.clone();
        for k in 0..n_basis {
            let f = focc[k];
            for i in 0..n_basis {
                cw[(i, k)] *= f;
            }
        }
        let density = (&cw * coefficients.transpose()) * 2.0;

        // Shell- and atom-resolved Mulliken charges.
        let ps = &density * &s_mat;
        let mut new_qsh = vec![0.0f64; n_shells];
        for sh in 0..n_shells {
            new_qsh[sh] = shells[sh].ref_occ;
        }
        for mu in 0..n_basis {
            new_qsh[basis_to_shell[mu]] -= ps[(mu, mu)];
        }
        let mut new_qat = vec![0.0f64; n_atoms];
        for sh in 0..n_shells {
            new_qat[shells[sh].atom] += new_qsh[sh];
        }

        // Electronic energy: Tr(P H0) + E2 + E3 + T·S.
        let mut e_band = 0.0;
        for i in 0..n_basis {
            for j in 0..n_basis {
                e_band += density[(i, j)] * h_mat[(i, j)];
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
        e_elec = e_band + e2 + e3 + ts;

        let rms_dq: f64 = {
            let sum_sq: f64 = qsh
                .iter()
                .zip(new_qsh.iter())
                .map(|(o, n)| (n - o).powi(2))
                .sum();
            (sum_sq / n_shells as f64).sqrt()
        };
        let e_conv = iter > 0 && (e_elec - prev_e).abs() < 1e-6;
        if e_conv && rms_dq < 2e-5 {
            converged = true;
            qsh = new_qsh;
            break;
        }
        prev_e = e_elec;

        mixer.set(&qsh);
        mixer.diff(&new_qsh);
    }

    // Final charges and dispersion/repulsion/halogen terms.
    let mut atom_charges = vec![0.0f64; n_atoms];
    let mut atom_shell_charges = vec![vec![0.0f64; 3]; n_atoms];
    for sh in 0..n_shells {
        atom_charges[shells[sh].atom] += qsh[sh];
        atom_shell_charges[shells[sh].atom][shells[sh].l as usize] = qsh[sh];
    }

    let e_rep = compute_repulsion(elements, &pos_bohr);
    let d3 = compute_d3_energy_with_params(elements, &pos_bohr, D3Params::gfn1());
    let e_xb = compute_halogen_bond(elements, &pos_bohr);
    let total_ha = e_elec + e_rep + d3.energy + e_xb;

    let orb_ev: Vec<f64> = orbital_energies
        .iter()
        .map(|e| e * EV_PER_HARTREE)
        .collect();
    let homo_idx = if n_occ > 0 { n_occ - 1 } else { 0 };
    let homo_energy = orb_ev[homo_idx];
    let lumo_energy = if n_occ < n_basis {
        orb_ev[n_occ]
    } else {
        homo_energy
    };

    Ok(Gfn1Result {
        orbital_energies: orb_ev,
        electronic_energy: e_elec * EV_PER_HARTREE,
        repulsive_energy: e_rep * EV_PER_HARTREE,
        dispersion_energy: d3.energy * EV_PER_HARTREE,
        dispersion_e6: d3.e6 * EV_PER_HARTREE,
        dispersion_e8: d3.e8 * EV_PER_HARTREE,
        halogen_bond_energy: e_xb * EV_PER_HARTREE,
        total_energy: total_ha * EV_PER_HARTREE,
        n_basis,
        n_electrons,
        homo_energy,
        lumo_energy,
        gap: lumo_energy - homo_energy,
        mulliken_charges: atom_charges,
        shell_charges: atom_shell_charges,
        scc_iterations: scc_iter,
        converged,
    })
}

#[cfg(test)]
mod decomposition_tests {
    use super::*;

    fn toluene() -> (Vec<u8>, Vec<[f64; 3]>) {
        let conf = crate::embed("Cc1ccccc1", 42);
        let pos: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        (conf.elements, pos)
    }

    #[test]
    fn dispersion_components_sum_to_total() {
        let (elements, positions) = toluene();
        let result = solve_gfn1(&elements, &positions).unwrap();
        let sum = result.dispersion_e6 + result.dispersion_e8;
        assert!(
            (sum - result.dispersion_energy).abs() < 1e-9,
            "{} vs {}",
            sum,
            result.dispersion_energy
        );
        assert!(result.dispersion_e6 <= 0.0 && result.dispersion_e8 <= 0.0);
    }
}
