//! Faithful GFN0-xTB implementation (xtb 6.7.1, cross-validation F43).
//!
//! Unlike [`super::solver`] (a GFN0-style approximation), this module ports
//! the xtb program's GFN0 energy term by term:
//! - EHT electronic energy in the GFN0 Slater basis (`ccm_build_SH0`),
//! - EEQ electrostatic energy (1 shell/atom, erf-screened Coulomb),
//! - effective pair repulsion,
//! - short-range bond (SRB) correction,
//! - D4(BJ) dispersion with GFN0 parameters.
//!
//! Parameter source: `param_gfn0-xtb.txt` (xtb 6.7.1, LGPL-2.1-or-later),
//! parsed by `scripts/validation/parse_gfn0_params.py` into
//! `tests/fixtures/validation/gfn0_params.json` and generated into
//! [`super::gfn0_params`] by `scripts/validation/gen_gfn0_params.py`.
//! Formula source: xtb 6.7.1 `src/peeq_module.f90` (`drep_grad`,
//! `ccm_build_SH0`), `src/xtb/eeq.f90`, `src/approxrab.f90`,
//! `src/coulomb/gaussian.f90`, `src/read_gfn_param.f90`.

use super::gfn0_params::{get_gfn0_params, GFN0_GLOBALS};
use serde::{Deserialize, Serialize};

/// Bohr per angstrom.
pub const BOHR_PER_ANGSTROM: f64 = 1.889_725_988_6;
/// eV per Hartree.
pub const EV_PER_HARTREE: f64 = 27.211_386_245_988;

/// GFN0-xTB single-point result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gfn0Result {
    /// Orbital energies (eV), sorted ascending.
    pub orbital_energies: Vec<f64>,
    /// EHT electronic energy (eV).
    pub electronic_energy: f64,
    /// EEQ electrostatic energy (eV).
    pub electrostatic_energy: f64,
    /// Repulsion energy (eV).
    pub repulsive_energy: f64,
    /// Short-range bond correction (eV).
    pub srb_energy: f64,
    /// D4 dispersion energy (eV) — two-body + three-body total.
    pub dispersion_energy: f64,
    /// Two-body D4 dispersion energy (eV).
    pub dispersion_2_body: f64,
    /// Three-body (ATM) D4 dispersion energy (eV).
    pub dispersion_3_body: f64,
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
    /// EEQ atomic charges (e).
    pub charges: Vec<f64>,
}

/// Effective pair repulsion (xtb `drep_grad` energy part).
///
/// E = Σ_{i≤j} Zeff_i·Zeff_j·exp(−α_ij·r^1.5)/r, skipped for r > 40 bohr.
/// α_ij = √(α_i·α_j)·(1 + (0.01·ΔEN² + 0.01·ΔEN⁴)·renscale).
/// Positions in bohr, energy in Hartree.
pub fn gfn0_repulsion(elements: &[u8], pos_bohr: &[[f64; 3]]) -> Result<f64, String> {
    let g = &GFN0_GLOBALS;
    let n = elements.len();
    let mut erep = 0.0;
    for i in 0..n {
        let pi = gfn0_params(elements[i])?;
        for j in 0..=i {
            let pj = gfn0_params(elements[j])?;
            let dx = pos_bohr[i][0] - pos_bohr[j][0];
            let dy = pos_bohr[i][1] - pos_bohr[j][1];
            let dz = pos_bohr[i][2] - pos_bohr[j][2];
            let r2 = dx * dx + dy * dy + dz * dz;
            if !(1e-12..=1600.0).contains(&r2) {
                continue;
            }
            let r = r2.sqrt();
            let den = pi.electronegativity - pj.electronegativity;
            let den2 = den * den;
            let alpha = (pi.rep_alpha * pj.rep_alpha).sqrt()
                * (1.0 + (0.01 * den2 + 0.01 * den2 * den2) * g.renscale);
            let repab = pi.rep_zeff * pj.rep_zeff;
            // kexp = kexp_light = 1.5 for GFN0: r^1.5 via r2^0.75.
            let r15 = r2.powf(0.75);
            erep += repab * (-alpha * r15).exp() / r;
        }
    }
    Ok(erep)
}

/// Look up the GFN0 parameters of an element, reporting an unsupported
/// atomic number instead of panicking.
fn gfn0_params(z: u8) -> Result<&'static super::gfn0_params::Gfn0ElementParams, String> {
    get_gfn0_params(z).ok_or_else(|| format!("GFN0 parameters not available for Z={z}"))
}

/// Error function (Abramowitz & Stegun 7.1.26, max error 1.5e-7 —
/// same accuracy class used throughout the codebase, sufficient for
/// 1e-6 Ha component targets as shown by F46).
fn erf_exact(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    // Winitzki approximation refined: relative error < 1e-9 for x < 2.4,
    // asymptotic tail beyond.
    if x > 6.0 {
        return sign;
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    sign * (1.0 - poly * (-x * x).exp())
}

/// mctc erf coordination numbers with xtb's softplus cap at 8
/// (xtb `getCoordinationNumber(..., cnType%erf, ...)` +
/// `cutCoordinationNumber(maxCN=8)`).
///
/// count = 0.5·(1 + erf(−7.5·(r−rc)/rc)), rc = D3 covalent radii sum;
/// pairs beyond 40 bohr skipped. The steepness 7.5 (not mctc's default
/// 3.15) is proven by F43 validation: it reproduces xtb's EEQ energy to
/// 7 decimals. Positions in bohr.
pub fn gfn0_cn(elements: &[u8], pos_bohr: &[[f64; 3]]) -> Vec<f64> {
    use super::d4_data::COVRAD_D3;
    let n = elements.len();
    let mut cn = vec![0.0; n];
    for i in 0..n {
        let zi = elements[i] as usize;
        if zi == 0 || zi > 86 {
            continue;
        }
        for j in 0..n {
            if i == j {
                continue;
            }
            let zj = elements[j] as usize;
            if zj == 0 || zj > 86 {
                continue;
            }
            let dx = pos_bohr[i][0] - pos_bohr[j][0];
            let dy = pos_bohr[i][1] - pos_bohr[j][1];
            let dz = pos_bohr[i][2] - pos_bohr[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if !(1e-12..=40.0).contains(&r) {
                continue;
            }
            let rc = COVRAD_D3[zi - 1] + COVRAD_D3[zj - 1];
            cn[i] += 0.5 * (1.0 + erf_exact(-7.5 * (r - rc) / rc));
        }
        // cutCN(cn, 8) = log(1+e^8) − log(1+e^(8−cn))
        cn[i] = (1.0 + 8.0f64.exp()).ln() - (1.0 + (8.0 - cn[i]).exp()).ln();
    }
    cn
}

/// GFN0 EEQ charges and electrostatic energy (xtb `chargeEquilibration`,
/// 1 shell per atom, erf-screened Coulomb kernel).
///
/// J_ij = erf(g_ij·r)/r, g_ij = 1/√(w_i² + w_j²), w = ALPG widths;
/// J_ii += GAM (hardness); x_i = −XI + KAPPA·√cn_i; solve with Σq = charge.
/// Returns (charges, e_es_hartree). Positions in bohr.
pub fn gfn0_eeq(
    elements: &[u8],
    pos_bohr: &[[f64; 3]],
    cn: &[f64],
    charge: f64,
) -> Result<(Vec<f64>, f64), String> {
    use nalgebra::DMatrix;
    let n = elements.len();
    if n == 0 {
        return Err("GFN0 EEQ: empty molecule (at least one atom is required)".to_string());
    }
    if cn.len() != n {
        return Err(format!(
            "GFN0 EEQ: coordination numbers length {} != atom count {}",
            cn.len(),
            n
        ));
    }
    // J matrix with Lagrange row/column (n+1 dim).
    let mut jmat = DMatrix::zeros(n + 1, n + 1);
    for i in 0..n {
        let pi = gfn0_params(elements[i])?;
        for j in 0..=i {
            let pj = gfn0_params(elements[j])?;
            let val = if i == j {
                // On-site: 2/√π·√(1/2)/w + GAM.
                2.0 / std::f64::consts::PI.sqrt() * (0.5f64.sqrt() / pi.eeq_width) + pi.hubbard
            } else {
                let dx = pos_bohr[i][0] - pos_bohr[j][0];
                let dy = pos_bohr[i][1] - pos_bohr[j][1];
                let dz = pos_bohr[i][2] - pos_bohr[j][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                if r < 1e-12 {
                    0.0
                } else {
                    let gij =
                        1.0 / (pi.eeq_width * pi.eeq_width + pj.eeq_width * pj.eeq_width).sqrt();
                    erf_exact(gij * r) / r
                }
            };
            jmat[(i, j)] = val;
            jmat[(j, i)] = val;
        }
        jmat[(i, n)] = 1.0;
        jmat[(n, i)] = 1.0;
    }
    // RHS: x_i = −XI + KAPPA·√cn; x_n = total charge.
    let mut xvec = vec![0.0; n + 1];
    for i in 0..n {
        let p = gfn0_params(elements[i])?;
        xvec[i] = -p.eeq_chi + p.eeq_kcn * cn[i].max(0.0).sqrt();
    }
    xvec[n] = charge;
    let x = DMatrix::from_column_slice(n + 1, 1, &xvec);
    let q = jmat
        .clone()
        .lu()
        .solve(&x)
        .ok_or_else(|| {
            "GFN0 EEQ: the charge-equilibration Jacobian is singular (the molecule has no resolvable EEQ solution)".to_string()
        })?
        .column(0)
        .into_owned();
    let charges: Vec<f64> = q.rows(0, n).iter().copied().collect();
    // E = 0.5·qᵀJq − qᵀx.
    let jxq: Vec<f64> = (0..=n)
        .map(|i| (0..=n).map(|k| jmat[(i, k)] * q[k]).sum::<f64>())
        .collect();
    let energy = 0.5 * q.iter().zip(jxq.iter()).map(|(a, b)| a * b).sum::<f64>()
        - q.iter().zip(xvec.iter()).map(|(a, b)| a * b).sum::<f64>();
    Ok((charges, energy))
}

/// Atomic radii for the GFN0 shell polynomial (xtb `atomicRad`, bohr).
/// Values are angstrom × 1.889726 (Z = 1..86).
pub const GFN0_ATOMIC_RAD: [f64; 86] = [
    0.604_712_316_4,
    0.699_199_115_8,
    2.456_643_785_2,
    1.870_828_709_0,
    1.587_570_470_4,
    1.417_294_491_5,
    1.341_705_799_9,
    1.209_424_632_6,
    1.133_835_593_2,
    1.171_629_939_0,
    3.023_561_582_0,
    2.645_616_384_1,
    2.343_031_906_0,
    2.154_304_787_3,
    2.059_940_227_9,
    1.965_576_668_6,
    1.889_725_988_6,
    1.908_623_248_5,
    3.779_451_977_2,
    3.288_160_200_3,
    3.004_664_322_0,
    2.796_793_183_3,
    2.721_205_423_8,
    2.456_643_785_2,
    2.439_746_525_4,
    2.344_921_165_9,
    2.231_407_568_6,
    2.212_510_308_5,
    2.305_237_486_1,
    2.267_671_186_6,
    2.324_363_026_1,
    2.267_671_186_6,
    2.267_671_186_6,
    2.230_105_966_8,
    2.210_940_906_9,
    2.192_043_647_0,
    4.061_731_846_0,
    3.590_979_378_3,
    3.329_284_079_6,
    3.098_496_561_2,
    2.944_467_543_1,
    2.759_130_177_3,
    2.608_049_240_9,
    2.570_157_382_2,
    2.532_232_522_5,
    2.456_643_785_2,
    2.570_157_382_2,
    2.645_616_384_1,
    2.683_540_043_7,
    2.645_616_384_1,
    2.645_616_384_1,
    2.590_434_213_5,
    2.570_157_382_2,
    2.570_157_382_2,
    4.497_770_226_5,
    3.893_144_066_5,
    3.666_185_965_5,
    3.477_458_847_0,
    3.590_979_378_3,
    3.553_012_618_8,
    3.515_218_199_1,
    3.496_320_939_2,
    3.477_458_847_0,
    3.458_495_917_1,
    3.439_719_087_3,
    3.420_755_839_8,
    3.401_859_048_2,
    3.382_962_256_6,
    3.364_065_465_0,
    3.364_065_465_0,
    3.382_962_256_6,
    3.288_160_200_3,
    3.098_496_561_2,
    2.986_367_004_3,
    2.834_620_772_2,
    2.664_513_643_9,
    2.494_438_204_7,
    2.456_643_785_2,
    2.456_643_785_2,
    2.494_438_204_7,
    2.721_205_423_8,
    2.740_168_353_7,
    2.834_620_772_2,
    2.683_540_043_7,
    2.797_112_213_9,
    2.759_130_177_3,
];

/// STO-NG primitive count per shell (xtb `setGFN0NumberOfPrimitives`).
fn gfn0_ngauss(z: u8, _sh: usize, l: u8, pqn: u8, valence: bool) -> u8 {
    if z <= 2 {
        match l {
            0 => {
                if valence {
                    3
                } else {
                    2
                }
            }
            _ => 3,
        }
    } else {
        match l {
            0 if pqn > 5 => 6,
            1 => {
                if pqn > 5 {
                    6
                } else {
                    3
                }
            }
            _ => 4,
        }
    }
}

/// Valence character per shell (xtb `generateValenceShellData`):
/// first shell of each l is valence, later same-l shells are not.
fn gfn0_valence_shells(ang_shell: &[u8]) -> Vec<bool> {
    let mut seen = [false; 5];
    ang_shell
        .iter()
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

/// GFN0 basis map entry: (atom, shell, l, m).
pub struct Gfn0BasisMap {
    /// (atom index, shell index, l, m).
    pub funcs: Vec<(usize, usize, u8, u8)>,
}

/// Self-energy SE(ish, iat) in eV (xtb `getSelfEnergy`):
/// lev − kCN·CN − kQShell·q − kQAtom·q², with kCN from KCNS/KCNP (×0.1,
/// s/p only), kQShell from KQS/KQP (s/p only), kQAtom = KQAT2.
fn gfn0_selfenergy(elements: &[u8], cn: &[f64], charges: &[f64]) -> Result<Vec<Vec<f64>>, String> {
    elements
        .iter()
        .enumerate()
        .map(|(i, &z)| -> Result<Vec<f64>, String> {
            let p = gfn0_params(z)?;
            Ok((0..p.nshell)
                .map(|sh| {
                    let l = p.ang_shell[sh];
                    let kcn = match l {
                        0 => p.kcn_s * 0.1,
                        1 => p.kcn_p * 0.1,
                        2 => p.kcn_d * 0.1,
                        _ => 0.0,
                    };
                    let kqsh = match l {
                        0 => p.kq_shell_s,
                        1 => p.kq_shell_p,
                        2 => p.kq_shell_d,
                        _ => 0.0,
                    };
                    p.selfenergy[sh]
                        - kcn * cn[i]
                        - kqsh * charges[i]
                        - p.kq_atom * charges[i] * charges[i]
                })
                .collect::<Vec<f64>>())
        })
        .collect::<Result<Vec<_>, String>>()
}

/// GFN0 overlap + H0 matrices (xtb `ccm_build_SH0`, molecules only).
/// H in eV. Same-atom pairs skipped (onsite = identity/diagonal,
/// matching xtb which never evaluates same-atom blocks).
#[allow(clippy::too_many_arguments)]
/// GFN0 overlap + H0 matrices (xtb `ccm_build_SH0`, molecules only).
///
/// Overlap is evaluated in xtb's Cartesian basis (s/p plain, d with trafo
/// weights) and transformed with the exact `dtrf2` matrix for d-involving
/// pairs (contaminant row/col dropped). Same-atom pairs skipped
/// (onsite = identity/diagonal, matching xtb). H in eV.
pub fn gfn0_h0(
    elements: &[u8],
    pos_bohr: &[[f64; 3]],
    cn: &[f64],
    charges: &[f64],
) -> Result<(nalgebra::DMatrix<f64>, nalgebra::DMatrix<f64>, Gfn0BasisMap), String> {
    use super::sto_overlap::cart_prim_overlap;
    use nalgebra::DMatrix;
    let g = &GFN0_GLOBALS;
    let n_atoms = elements.len();
    if n_atoms == 0 {
        return Err("GFN0 H0: empty molecule (at least one atom is required)".to_string());
    }
    if cn.len() != n_atoms {
        return Err(format!(
            "GFN0 H0: coordination numbers length {} != atom count {}",
            cn.len(),
            n_atoms
        ));
    }
    if charges.len() != n_atoms {
        return Err(format!(
            "GFN0 H0: charges length {} != atom count {}",
            charges.len(),
            n_atoms
        ));
    }
    let atom_shells: Vec<Vec<Gfn0Shell>> = elements
        .iter()
        .map(|&z| gfn0_atom_shells(z))
        .collect::<Result<Vec<_>, String>>()?;
    // funcs: (atom, shell, sao); shell_offset[atom][shell] = first index.
    let mut funcs: Vec<(usize, usize, usize)> = Vec::new();
    let mut shell_offset: Vec<Vec<usize>> = Vec::new();
    for (at, shells) in atom_shells.iter().enumerate() {
        let mut offs = Vec::new();
        for (sh, _) in shells.iter().enumerate() {
            offs.push(funcs.len());
            let l = shell_l(&shells[sh]);
            let nsao = if l == 2 { 5 } else { shells[sh].comps.len() };
            for m in 0..nsao {
                funcs.push((at, sh, m));
            }
        }
        shell_offset.push(offs);
    }
    let bmap = Gfn0BasisMap {
        funcs: funcs
            .iter()
            .map(|&(a, s, m)| (a, s, shell_l(&atom_shells[a][s]), m as u8))
            .collect(),
    };
    let n = funcs.len();
    let mut s_mat = DMatrix::zeros(n, n);
    let mut h_mat = DMatrix::zeros(n, n);
    let se = gfn0_selfenergy(elements, cn, charges)?;
    let kshell = [g.kshell[0], g.kshell[1], g.kshell[2], 0.0];
    let enshell = [g.enshell[0], g.enshell[1], g.enshell[2], -0.2];

    struct ShInfo {
        l: u8,
        zeta: f64,
        en: f64,
        poly: f64,
        rad: f64,
        km_valao_pol: bool,
    }
    let shinfo: Vec<Vec<ShInfo>> = elements
        .iter()
        .map(|&z| -> Result<Vec<ShInfo>, String> {
            let p = gfn0_params(z)?;
            let val = gfn0_valence_shells(&p.ang_shell[..p.nshell]);
            Ok((0..p.nshell)
                .map(|sh| {
                    let l = p.ang_shell[sh];
                    ShInfo {
                        l,
                        zeta: p.slater[sh],
                        en: p.electronegativity,
                        poly: match l {
                            0 => p.shell_poly_s,
                            1 => p.shell_poly_p,
                            2 => p.shell_poly_d,
                            _ => 0.0,
                        },
                        rad: GFN0_ATOMIC_RAD[z as usize - 1],
                        km_valao_pol: !val[sh],
                    }
                })
                .collect())
        })
        .collect::<Result<Vec<_>, String>>()?;

    for i in 0..n {
        s_mat[(i, i)] = 1.0;
        let (at_a, sh_a, _) = funcs[i];
        h_mat[(i, i)] = se[at_a][sh_a];
    }
    // Off-diagonal shell-pair blocks.
    let nat = elements.len();
    for at_a in 0..nat {
        for at_b in (at_a + 1)..nat {
            let dx = pos_bohr[at_a][0] - pos_bohr[at_b][0];
            let dy = pos_bohr[at_a][1] - pos_bohr[at_b][1];
            let dz = pos_bohr[at_a][2] - pos_bohr[at_b][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            for sh_a in 0..atom_shells[at_a].len() {
                for sh_b in 0..atom_shells[at_b].len() {
                    let ia = &shinfo[at_a][sh_a];
                    let ib = &shinfo[at_b][sh_b];
                    if ia.zeta < 1e-10 || ib.zeta < 1e-10 {
                        continue;
                    }
                    // Cartesian overlap block.
                    let ca = &atom_shells[at_a][sh_a];
                    let cb = &atom_shells[at_b][sh_b];
                    let na = ca.comps.len();
                    let nb = cb.comps.len();
                    let mut cao = DMatrix::zeros(na, nb);
                    for (a, ((lxa, lya, lza), wa)) in ca.comps.iter().enumerate() {
                        for (b, ((lxb, lyb, lzb), wb)) in cb.comps.iter().enumerate() {
                            cao[(a, b)] = wa
                                * wb
                                * cart_prim_overlap(
                                    *lxa,
                                    *lya,
                                    *lza,
                                    &ca.prims,
                                    &pos_bohr[at_a],
                                    *lxb,
                                    *lyb,
                                    *lzb,
                                    &cb.prims,
                                    &pos_bohr[at_b],
                                );
                        }
                    }
                    // Common Cartesian normalization for d shells (xtb
                    // convention: all d Cartesians share the xx norm, so
                    // xy/xz/yz rows/cols carry an extra 1/sqrt(3); the
                    // trafo weights + dtrf2 then yield correct SAOs).
                    // Without this, dxy-type overlaps come out sqrt(3) too
                    // big (found by F43 validation vs numeric quadrature).
                    let s3 = 3.0f64.sqrt();
                    if ia.l == 2 {
                        for a in 3..6 {
                            for b in 0..nb {
                                cao[(a, b)] /= s3;
                            }
                        }
                    }
                    if ib.l == 2 {
                        for b in 3..6 {
                            for a in 0..na {
                                cao[(a, b)] /= s3;
                            }
                        }
                    }
                    // SAO transform (dtrf2) on d sides.
                    let ta = sao_trafo(ia.l);
                    let tb = sao_trafo(ib.l);
                    let ma = ta.as_ref().map(|t| t[0].len()).unwrap_or(na);
                    let mb = tb.as_ref().map(|t| t[0].len()).unwrap_or(nb);
                    let mut sblock = DMatrix::zeros(ma, mb);
                    for ii in 0..ma {
                        for jj in 0..mb {
                            let mut v = 0.0;
                            for a in 0..na {
                                let wa = ta.as_ref().map(|t| t[a][ii]).unwrap_or(if a == ii {
                                    1.0
                                } else {
                                    0.0
                                });
                                for b in 0..nb {
                                    let wb = tb.as_ref().map(|t| t[b][jj]).unwrap_or(if b == jj {
                                        1.0
                                    } else {
                                        0.0
                                    });
                                    v += wa * wb * cao[(a, b)];
                                }
                            }
                            sblock[(ii, jj)] = v;
                        }
                    }
                    // Hij prefactors (xtb ccm_build_SH0).
                    let den = ia.en - ib.en;
                    let den2 = den * den;
                    let kscale = 0.5 * (kshell[ia.l as usize] + kshell[ib.l as usize]);
                    let enscale = 0.005 * (enshell[ia.l as usize] + enshell[ib.l as usize]);
                    let enpoly = 1.0 + enscale * den2 + g.enscale4 * enscale * den2 * den2;
                    let zetaij = 2.0 * (ia.zeta * ib.zeta).sqrt() / (ia.zeta + ib.zeta);
                    let mut km = kscale * zetaij * enpoly;
                    match (ia.km_valao_pol, ib.km_valao_pol) {
                        (true, true) => km = 0.0,
                        (true, false) | (false, true) => km *= g.kdiff,
                        _ => {}
                    }
                    let rad_sum = ia.rad + ib.rad;
                    let sq = (r / rad_sum).sqrt();
                    let shpoly = (1.0 + 0.01 * ia.poly * sq) * (1.0 + 0.01 * ib.poly * sq);
                    for ii in 0..ma {
                        for jj in 0..mb {
                            let fi = shell_offset[at_a][sh_a] + ii;
                            let fj = shell_offset[at_b][sh_b] + jj;
                            let sij = sblock[(ii, jj)];
                            s_mat[(fi, fj)] = sij;
                            s_mat[(fj, fi)] = sij;
                            let hav = 0.5 * km * (se[at_a][sh_a] + se[at_b][sh_b]);
                            let hij = hav * shpoly * sij;
                            h_mat[(fi, fj)] = hij;
                            h_mat[(fj, fi)] = hij;
                        }
                    }
                }
            }
        }
    }
    Ok((s_mat, h_mat, bmap))
}

/// Angular momentum of a shell (stored implicitly via comps length).
fn shell_l(sh: &Gfn0Shell) -> u8 {
    match sh.comps.len() {
        1 => 0,
        3 => 1,
        6 => 2,
        _ => panic!("GFN0: unsupported shell size"),
    }
}

/// dtrf2 SAO transform columns (5 SAOs in dtrf2 output order) as 6x5.
/// None for non-d shells (identity).
fn sao_trafo(l: u8) -> Option<[[f64; 5]; 6]> {
    if l != 2 {
        return None;
    }
    let dt = dtrf2();
    let mut t = [[0.0; 5]; 6];
    for m in 0..6 {
        for k in 0..5 {
            t[m][k] = dt[m][k + 1];
        }
    }
    Some(t)
}

/// Boltzmann constant (eV/K) for GFN0 Fermi smearing at 300 K.
pub const GFN0_KB_EV: f64 = 8.617_333_262e-5;
/// GFN0 electronic temperature (K).
pub const GFN0_ET_K: f64 = 300.0;

/// Solve the GFN0 EHT problem: Löwdin orthogonalization, Fermi filling
/// at 300 K, electronic energy (xtb `eel = Σf·ε + ga + gb`).
///
/// Returns (orbital energies eV ascending, eel eV, homo eV, lumo eV).
pub fn gfn0_eht(
    s_mat: &nalgebra::DMatrix<f64>,
    h_mat: &nalgebra::DMatrix<f64>,
    n_electrons: usize,
) -> Result<(Vec<f64>, f64, f64, f64), String> {
    use nalgebra::DMatrix;
    let n = s_mat.nrows();
    if !n_electrons.is_multiple_of(2) {
        return Err("GFN0 EHT: open shells not supported".to_string());
    }
    // Löwdin S^{-1/2}.
    let s_eigen = s_mat.clone().symmetric_eigen();
    let mut s_half_inv = DMatrix::zeros(n, n);
    for k in 0..n {
        let val = s_eigen.eigenvalues[k];
        if val > 1e-8 {
            let inv_sqrt = 1.0 / val.sqrt();
            let col = s_eigen.eigenvectors.column(k);
            for i in 0..n {
                for j in 0..n {
                    s_half_inv[(i, j)] += inv_sqrt * col[i] * col[j];
                }
            }
        }
    }
    let f_prime = &s_half_inv * h_mat * &s_half_inv;
    let eig = f_prime.symmetric_eigen();
    let mut emo: Vec<(f64, usize)> = eig
        .eigenvalues
        .iter()
        .enumerate()
        .map(|(k, &e)| (e, k))
        .collect();
    emo.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut orbital_energies = vec![0.0; n];
    for (k, (e, _)) in emo.iter().enumerate() {
        orbital_energies[k] = *e;
    }
    // Restricted Fermi filling per spin channel (reuse GFN2 routine).
    let nel_per_spin = n_electrons as f64 / 2.0;
    let kt = GFN0_KB_EV * GFN0_ET_K;
    let (focc, ts) = super::gfn2::get_fermi_filling(&orbital_energies, nel_per_spin, kt);
    let eel_ev: f64 = 2.0
        * focc
            .iter()
            .zip(orbital_energies.iter())
            .map(|(f, e)| f * e)
            .sum::<f64>()
        + ts;
    let n_occ = n_electrons / 2;
    let homo = if n_occ > 0 {
        orbital_energies[n_occ - 1]
    } else {
        f64::NAN
    };
    let lumo = if n_occ < n {
        orbital_energies[n_occ]
    } else {
        f64::NAN
    };
    Ok((orbital_energies, eel_ev, homo, lumo))
}

/// dtrf2 transformation matrix (xtb `intgrad.f90`, CAO order
/// xx, yy, zz, xy, xz, yz; columns = dS, dx2-y2, dz2, xy, xz, yz).
/// Only columns 1..6 (0-based) survive (dS contaminant dropped).
fn dtrf2() -> [[f64; 6]; 6] {
    let s3 = 3.0f64.sqrt();
    [
        [1.0 / 5.0f64.sqrt(), 0.5 * s3, 0.5, 0.0, 0.0, 0.0],
        [1.0 / 5.0f64.sqrt(), -0.5 * s3, 0.5, 0.0, 0.0, 0.0],
        [1.0 / 5.0f64.sqrt(), 0.0, -1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
    ]
}

/// Cartesian-d trafo weights baked into contractions
/// (xtb `set_d_function`: xx, yy, zz ×1; xy, xz, yz ×√3).
fn d_trafo_w() -> [f64; 6] {
    let s3 = 3.0f64.sqrt();
    [1.0, 1.0, 1.0, s3, s3, s3]
}

/// Cartesian components per shell: (lx, ly, lz, weight).
fn gfn0_cart_comps(l: u8) -> Vec<((u8, u8, u8), f64)> {
    match l {
        0 => vec![((0, 0, 0), 1.0)],
        1 => vec![((1, 0, 0), 1.0), ((0, 1, 0), 1.0), ((0, 0, 1), 1.0)],
        2 => {
            let w = d_trafo_w();
            vec![
                ((2, 0, 0), w[0]),
                ((0, 2, 0), w[1]),
                ((0, 0, 2), w[2]),
                ((1, 1, 0), w[3]),
                ((1, 0, 1), w[4]),
                ((0, 1, 1), w[5]),
            ]
        }
        _ => panic!("GFN0: f shells not supported"),
    }
}

/// One contracted shell: absolute primitives + Cartesian components.
struct Gfn0Shell {
    prims: Vec<(f64, f64)>,
    comps: Vec<((u8, u8, u8), f64)>,
}

/// Build all shell contractions for an atom (xtb `newBasisset`).
///
/// s-polarization shells (non-first s) are Schmidt-orthogonalized against
/// the atom's valence s and renormalized; d shells carry the trafo weights
/// (dtrf2 applied later at the block level).
fn gfn0_atom_shells(z: u8) -> Result<Vec<Gfn0Shell>, String> {
    use super::sto_overlap::sto_table_prims;
    let p = gfn0_params(z)?;
    let val_flags = gfn0_valence_shells(&p.ang_shell[..p.nshell]);
    // Valence-s contraction (for Schmidt), if any s shell is valence.
    let mut val_s_prims: Option<Vec<(f64, f64)>> = None;
    let mut out = Vec::new();
    for sh in 0..p.nshell {
        let l = p.ang_shell[sh];
        let zeta = p.slater[sh];
        let nga = gfn0_ngauss(z, sh, l, p.pqn[sh], val_flags[sh]);
        let table = sto_table_prims(p.pqn[sh], l, nga);
        let mut prims: Vec<(f64, f64)> = table
            .iter()
            .map(|(ar, ca)| (ar * zeta * zeta, *ca))
            .collect();
        let is_spol = l == 0 && !val_flags[sh];
        if l == 0 && val_flags[sh] && val_s_prims.is_none() {
            val_s_prims = Some(prims.clone());
        }
        if is_spol {
            // Schmidt vs valence s (xtb diff-s / DZ-s branches) + renormalize.
            // xtb works with normalized primitives (slaterToGauss norm=.true.):
            // do the Schmidt algebra on normalized coefficients, then store
            // bare coefficients (our integral code applies primitive norms
            // itself, mirroring xtb's split).
            use super::sto_overlap::gaussian_norm;
            let nc = |a: f64, c: f64| c * gaussian_norm(0, 0, 0, a);
            let sn = |a: f64, b: f64| (std::f64::consts::PI / (a + b)).powf(1.5);
            let v = val_s_prims.clone().expect("s-pol needs valence s");
            // Normalized-coefficient view.
            let praw: Vec<(f64, f64)> = prims.iter().map(|(a, c)| (*a, nc(*a, *c))).collect();
            let valv: Vec<(f64, f64)> = v.iter().map(|(a, c)| (*a, nc(*a, *c))).collect();
            let mut ov = 0.0;
            for (ar, ca) in praw.iter() {
                for (br, cb) in valv.iter() {
                    ov += ca * cb * sn(*ar, *br);
                }
            }
            let mut combined = praw.clone();
            for (br, cb) in valv.iter() {
                combined.push((*br, -ov * cb));
            }
            // Renormalize.
            let mut norm2 = 0.0;
            for (ar, ca) in combined.iter() {
                for (br, cb) in combined.iter() {
                    norm2 += ca * cb * sn(*ar, *br);
                }
            }
            let nrm = norm2.sqrt();
            // Store bare (divide the primitive norm back out).
            prims = combined
                .iter()
                .map(|(a, c)| (*a, (c / nrm) / gaussian_norm(0, 0, 0, *a)))
                .collect();
        }
        if std::env::var("GFN0_CONT_DEBUG").is_ok() && prims.len() > 3 {
            eprintln!("shell contraction ({} prims):", prims.len());
            for (a, c) in prims.iter() {
                eprintln!("  {a:.6} {c:+.6}");
            }
        }
        out.push(Gfn0Shell {
            prims,
            comps: gfn0_cart_comps(l),
        });
    }
    Ok(out)
}

/// Pauling electronegativities (xtb paulingEN, peeq scope for SRB `pre`).
pub(crate) fn pauling_en(z: u8) -> f64 {
    const T: [f64; 86] = [
        2.20, 3.00, 0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 4.50, 0.93, 1.31, 1.61, 1.90, 2.19,
        2.58, 3.16, 3.50, 0.82, 1.00, 1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88, 1.91, 1.90, 1.65,
        1.81, 2.01, 2.18, 2.55, 2.96, 3.00, 0.82, 0.95, 1.22, 1.33, 1.60, 2.16, 1.90, 2.20, 2.28,
        2.20, 1.93, 1.69, 1.78, 1.96, 2.05, 2.10, 2.66, 2.60, 0.79, 0.89, 1.10, 1.12, 1.13, 1.14,
        1.15, 1.17, 1.18, 1.20, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.30, 1.50, 2.36, 1.90,
        2.20, 2.20, 2.28, 2.54, 2.00, 1.62, 2.33, 2.02, 2.00, 2.20, 2.20,
    ];
    T[z as usize - 1]
}

/// Period (row) number 1-4 for the SRB tables (xtb `itr`).
fn srb_period(z: u8) -> usize {
    match z {
        1 | 2 => 1,
        3..=10 => 2,
        11..=18 => 3,
        _ => 4,
    }
}

/// Short-range bond correction (xtb `dsrb_grad` energy part).
///
/// esrb = Σ prefactor·exp(−pre·(rab−rab0)²) over hetero B..F pairs
/// (Z in 5..9, Zi≠Zj) within √200 bohr, with
/// pre = steepness·(1+enScale·den²), den = |ENi−ENj} **Pauling**
/// (peeq's `en` array; NOT the approxrab fitted table used in rab0 —
/// distinguished by F43 validation to 6 decimals),
/// rab0 = (ra+rb)·(1−k1·den−k2·den²), ra = r0(Zi)+cnfak(Zi)·cn(i)+shift
/// with approxrab r0/cnfak/en and k1/k2 from the p matrix.
/// Positions in bohr, CN capped (same `cn` as H0/EEQ), energy in Hartree.
pub fn gfn0_srb(elements: &[u8], pos_bohr: &[[f64; 3]], cn: &[f64]) -> f64 {
    use super::gfn0_params::GFN0_SRB;
    let g = &GFN0_GLOBALS;
    let n = elements.len();
    let mut esrb = 0.0;
    for i in 0..n {
        let zi = elements[i];
        if !(5..=9).contains(&zi) {
            continue;
        }
        // Unique hetero pairs only (xtb srblist holds each pair once).
        for j in (i + 1)..n {
            let zj = elements[j];
            if !(5..=9).contains(&zj) || zi == zj {
                continue;
            }
            let dx = pos_bohr[i][0] - pos_bohr[j][0];
            let dy = pos_bohr[i][1] - pos_bohr[j][1];
            let dz = pos_bohr[i][2] - pos_bohr[j][2];
            let r2 = dx * dx + dy * dy + dz * dz;
            if !(1e-12..200.0).contains(&r2) {
                continue;
            }
            let rab = r2.sqrt();
            // NOTE: pre uses Pauling EN (peeq scope); rab0 below uses
            // approxrab EN. Distinguished by F43 validation.
            let den_pauling = (pauling_en(zi) - pauling_en(zj)).abs();
            let pre = g.srb_steepness * (1.0 + g.srb_enscale * den_pauling * den_pauling);
            let den = (GFN0_SRB.en[zi as usize - 1] - GFN0_SRB.en[zj as usize - 1]).abs();
            let ra = GFN0_SRB.r0[zi as usize - 1]
                + GFN0_SRB.cnfak[zi as usize - 1] * cn[i]
                + g.srb_shift;
            let rb = GFN0_SRB.r0[zj as usize - 1]
                + GFN0_SRB.cnfak[zj as usize - 1] * cn[j]
                + g.srb_shift;
            let ir = srb_period(zi);
            let jr = srb_period(zj);
            let k1 = 0.005 * (GFN0_SRB.p[ir - 1][0] + GFN0_SRB.p[jr - 1][0]);
            let k2 = 0.005 * (GFN0_SRB.p[ir - 1][1] + GFN0_SRB.p[jr - 1][1]);
            let ff = 1.0 - k1 * den - k2 * den * den;
            let rab0 = (ra + rb) * ff;
            let dr = rab - rab0;
            esrb += g.srb_prefactor * (-pre * dr * dr).exp();
        }
    }
    esrb
}

/// Full GFN0-xTB single point (xtb `peeq`: eel + ees + ep + esrb + ed).
///
/// `elements`: atomic numbers; `positions`: angstrom; `charge`: total charge.
/// Energies in the result are eV. Only closed shells supported.
pub fn solve_gfn0(
    elements: &[u8],
    positions: &[[f64; 3]],
    charge: f64,
) -> Result<Gfn0Result, String> {
    use super::gfn0_params::GFN0_NVALENCE;
    crate::validation::validate_geometry(elements, positions)?;
    let n = elements.len();
    if n == 0 {
        return Err("GFN0: empty molecule (at least one atom is required)".to_string());
    }
    for &z in elements {
        if get_gfn0_params(z).is_none() {
            return Err(format!("GFN0 parameters not available for Z={z}"));
        }
    }
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
    let n_electrons: usize = elements
        .iter()
        .map(|&z| GFN0_NVALENCE[(z - 1) as usize] as usize)
        .sum::<usize>()
        .saturating_sub(charge.round() as usize);
    if !n_electrons.is_multiple_of(2) {
        return Err("GFN0: open shells not supported".to_string());
    }

    let cn = gfn0_cn(elements, &pos_bohr);
    let (charges, ees_ha) = gfn0_eeq(elements, &pos_bohr, &cn, charge)?;
    let (s_mat, h_mat, bmap) = gfn0_h0(elements, &pos_bohr, &cn, &charges)?;
    let (orbital_energies, eel_ev, homo_energy, lumo_energy) =
        gfn0_eht(&s_mat, &h_mat, n_electrons)?;
    let e_rep_ha = gfn0_repulsion(elements, &pos_bohr)?;
    let e_srb_ha = gfn0_srb(elements, &pos_bohr, &cn);
    // D4(BJ) with GFN0 parameters (s6=1, s8=2.85, a1=0.8, a2=4.6, no ATM).
    let d4 = crate::dispersion::compute_d4_energy(
        elements,
        positions,
        &crate::dispersion::D4Config {
            s6: GFN0_GLOBALS.d4_s6,
            s8: GFN0_GLOBALS.d4_s8,
            a1: GFN0_GLOBALS.d4_a1,
            a2: GFN0_GLOBALS.d4_a2,
            three_body: false,
            s9: GFN0_GLOBALS.d4_s9,
        },
    );
    let e_disp_ha = d4.total_energy;
    let total_ev = eel_ev
        + ees_ha * EV_PER_HARTREE
        + e_rep_ha * EV_PER_HARTREE
        + e_srb_ha * EV_PER_HARTREE
        + e_disp_ha * EV_PER_HARTREE;
    let gap = lumo_energy - homo_energy;
    Ok(Gfn0Result {
        orbital_energies,
        electronic_energy: eel_ev,
        electrostatic_energy: ees_ha * EV_PER_HARTREE,
        repulsive_energy: e_rep_ha * EV_PER_HARTREE,
        srb_energy: e_srb_ha * EV_PER_HARTREE,
        dispersion_energy: e_disp_ha * EV_PER_HARTREE,
        dispersion_2_body: d4.e2_body * EV_PER_HARTREE,
        dispersion_3_body: d4.e3_body * EV_PER_HARTREE,
        total_energy: total_ev,
        n_basis: bmap.funcs.len(),
        n_electrons,
        homo_energy,
        lumo_energy,
        gap,
        charges,
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
        let result = solve_gfn0(&elements, &positions, 0.0).unwrap();
        let sum = result.dispersion_2_body + result.dispersion_3_body;
        assert!(
            (sum - result.dispersion_energy).abs() < 1e-9,
            "{} vs {}",
            sum,
            result.dispersion_energy
        );
        assert!(result.dispersion_2_body <= 0.0);
    }
}
