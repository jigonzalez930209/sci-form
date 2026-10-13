//! Faithful PM3 (NDDO) solver — term-by-term port of MOPAC 23.2.5.
//!
//! Sources ported (openmopac/mopac, Apache-2.0):
//! - `src/integrals/diat.F90` + `set.F90` (analytic s/p Slater overlaps)
//! - `src/integrals/mndod.F90`: `reppd` (22 local two-electron integrals),
//!   `rotmat` (direction cosines), `rotatd`/`tx` (molecular-frame rotation),
//!   `spcore` (electron-nuclear attraction), `elenuc`
//! - `src/integrals/hcore.F90` + `h1elec.F90` (one-electron matrix)
//! - `src/integrals/ccrep.F90` (core-core repulsion with PM3 Gaussians)
//! - `src/models/calpar.F90` (derived α/ρ parameters)
//! - `src/SCF/fock2.F90`, `jab`/`kab` (NDDO Fock build)
//! - `src/SCF/iter.F90` (heat of formation)
//!
//! Conventions: distances in Å at the API boundary, bohr inside the integral
//! kernels; energies in eV; `kcal_per_ev` for the heat of formation.

use super::params::{get_pm3_params, Pm3Params};
use nalgebra::DMatrix;

/// eV per Hartree (MOPAC `ev`, 2018 CODATA).
pub(crate) const EV: f64 = 27.211_386_245_988;
/// Bohr radius in Å (MOPAC `a0`).
pub(crate) const A0: f64 = 0.529_177_210_903;
/// kcal/mol per eV (MOPAC `fpc_9`).
pub(crate) const KCAL_PER_EV: f64 = 23.060_547_830_619_029;

/// Derived per-element quantities (MOPAC `calpar` + `inid`).
#[derive(Debug, Clone, Copy)]
pub struct Pm3Derived {
    pub am: f64,
    pub ad: f64,
    pub aq: f64,
    pub dd: f64,
    pub qq: f64,
    pub po1: f64,
    pub po2: f64,
    pub po3: f64,
    pub po7: f64,
    pub po9: f64,
}

impl Pm3Derived {
    /// Port of the `calpar`/`inid` parameter derivation.
    pub fn new(p: &Pm3Params) -> Self {
        let zs = p.zeta_s;
        let zp = p.zeta_p.max(0.3);
        let hpp = (0.5 * (p.gpp - p.gp2)).max(0.1);
        let hsp = p.hsp.max(1e-7);
        let qn = p.nspqn as f64;

        let dd = (2.0 * qn + 1.0) * (4.0 * zs * zp).powf(qn + 0.5)
            / (zs + zp).powf(2.0 * qn + 2.0)
            / 3.0_f64.sqrt();
        let qq = ((4.0 * qn * qn + 6.0 * qn + 2.0) / 20.0).sqrt() / zp;

        // Newton iteration for the "d" additive term (MOPAC calpar).
        let gdd1 = (hsp / (EV * dd * dd)).powf(1.0 / 3.0);
        let (mut d1, mut d2) = (gdd1, gdd1 + 0.04);
        for _ in 0..5 {
            let df = d2 - d1;
            let hsp1 = 0.5 * d1 - 0.5 / (4.0 * dd * dd + 1.0 / (d1 * d1)).sqrt();
            let hsp2 = 0.5 * d2 - 0.5 / (4.0 * dd * dd + 1.0 / (d2 * d2)).sqrt();
            if (hsp2 - hsp1).abs() < 1e-25 {
                break;
            }
            let d3 = d1 + df * (hsp / EV - hsp1) / (hsp2 - hsp1);
            d1 = d2;
            d2 = d3;
        }
        let ad = d2;

        // Newton iteration for the "q" additive term.
        let gqq = (4.0 * hpp / (EV * 48.0 * qq.powi(4))).powf(0.2);
        let (mut q1, mut q2) = (gqq, gqq + 0.04);
        for _ in 0..5 {
            let qf = q2 - q1;
            let hpp1 = 0.25 * q1 - 0.5 / (4.0 * qq * qq + 1.0 / (q1 * q1)).sqrt()
                + 0.25 / (8.0 * qq * qq + 1.0 / (q1 * q1)).sqrt();
            let hpp2 = 0.25 * q2 - 0.5 / (4.0 * qq * qq + 1.0 / (q2 * q2)).sqrt()
                + 0.25 / (8.0 * qq * qq + 1.0 / (q2 * q2)).sqrt();
            if (hpp2 - hpp1).abs() < 1e-25 {
                break;
            }
            let q3 = q1 + qf * (hpp / EV - hpp1) / (hpp2 - hpp1);
            q1 = q2;
            q2 = q3;
        }
        let aq = q2;

        let am = p.gss / EV;
        let po1 = 0.5 / am;
        let mut po2 = 0.0;
        let mut po3 = 0.0;
        let mut dd = dd;
        let mut qq = qq;
        let mut ad = ad;
        let mut aq = aq;
        if p.z == 1 {
            // Hydrogen: MOPAC forces ad = aq = am and zeroes po2/po3.
            ad = am;
            aq = am;
            dd = 0.0;
            qq = 0.0;
        } else {
            po2 = 0.5 / ad;
            po3 = 0.5 / aq;
        }
        Pm3Derived {
            am,
            ad,
            aq,
            dd,
            qq,
            po1,
            po2,
            po3,
            po7: po1,
            po9: po1,
        }
    }
}

/// Heavy atom = has p functions (MOPAC `natorb >= 3`).
fn is_heavy(z: u8) -> bool {
    z != 1
}

// ─── Overlaps (MOPAC diat2 + set + aintgs/bintgs) ───────────────────────────

struct AbIntegrals {
    a: [f64; 8],
    b: [f64; 8],
    sa: f64,
    sb: f64,
    isp: usize,
    ips: usize,
}

fn aintgs(x: f64, k: usize, a: &mut [f64; 8]) {
    let c = (-x).exp();
    a[0] = c / x;
    for i in 1..=k {
        a[i] = (a[i - 1] * i as f64 + c) / x;
    }
}

fn bintgs(x: f64, k: usize, b: &mut [f64; 8]) {
    const FACT: [f64; 16] = [
        1.0,
        1.0,
        2.0,
        6.0,
        24.0,
        120.0,
        720.0,
        5040.0,
        40320.0,
        362880.0,
        3628800.0,
        39916800.0,
        479001600.0,
        6227020800.0,
        87178291200.0,
        1307674368000.0,
    ];
    let absx = x.abs();
    let direct =
        absx > 3.0 || (absx > 2.0 && k <= 10) || (absx > 1.0 && k <= 7) || (absx > 0.5 && k <= 5);
    if direct {
        let expx = x.exp();
        let expmx = 1.0 / expx;
        b[0] = (expx - expmx) / x;
        for i in 1..=k {
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            b[i] = (i as f64 * b[i - 1] + sign * expx - expmx) / x;
        }
        return;
    }
    if absx <= 1e-6 {
        for i in 0..=k {
            b[i] = if (i + 1) % 2 == 1 { 2.0 } else { 0.0 } / (i + 1) as f64;
        }
        return;
    }
    let last = if absx > 2.0 {
        15
    } else if absx > 1.0 {
        12
    } else if absx > 0.5 {
        7
    } else {
        6
    };
    for i in 0..=k {
        let mut y = 0.0;
        for m in 0..=last {
            let xf = if m == 0 { 1.0 } else { FACT[m.min(15)] };
            let two_mod = if (m + i + 1) % 2 == 1 { 2.0 } else { 0.0 };
            y += (-x).powi(m as i32) * two_mod / (xf * (m + i + 1) as f64);
        }
        b[i] = y;
    }
}

fn set_ab(esa: f64, esb: f64, na: u8, nb: u8, rab: f64, ii: usize) -> AbIntegrals {
    let (isp, ips, sa, sb) = if na <= nb {
        (1usize, 2usize, esa, esb)
    } else {
        (2, 1, esb, esa)
    };
    let j = if ii > 3 { ii + 1 } else { ii + 2 };
    let jcall = j - 1;
    let alpha = 0.5 * rab * (sa + sb);
    let beta = 0.5 * rab * (sb - sa);
    let mut a = [0.0; 8];
    let mut b = [0.0; 8];
    aintgs(alpha, jcall, &mut a);
    bintgs(beta, jcall, &mut b);
    AbIntegrals {
        a,
        b,
        sa,
        sb,
        isp,
        ips,
    }
}

/// Full molecular-frame overlap matrix for the sp basis (MOPAC `diat`/`diat2`).
pub(crate) fn overlap_matrix(elements: &[u8], positions: &[[f64; 3]]) -> DMatrix<f64> {
    let n_atoms = elements.len();
    let mut basis: Vec<(usize, usize)> = Vec::new(); // (atom, orbital 0=s,1=px,2=py,3=pz)
    for (i, &z) in elements.iter().enumerate() {
        basis.push((i, 0));
        if is_heavy(z) {
            for k in 1..4 {
                basis.push((i, k));
            }
        }
    }
    let n = basis.len();
    let mut s = DMatrix::zeros(n, n);
    for a in 0..n {
        s[(a, a)] = 1.0;
        for b in (a + 1)..n {
            let (ia, oa) = basis[a];
            let (ib, ob) = basis[b];
            if ia == ib {
                continue;
            }
            let pa = get_pm3_params(elements[ia]).unwrap();
            let pb = get_pm3_params(elements[ib]).unwrap();
            let r = ((positions[ia][0] - positions[ib][0]).powi(2)
                + (positions[ia][1] - positions[ib][1]).powi(2)
                + (positions[ia][2] - positions[ib][2]).powi(2))
            .sqrt();
            let val = overlap_sp(
                pa,
                elements[ia],
                pb,
                elements[ib],
                r,
                oa,
                ob,
                &positions[ia],
                &positions[ib],
            );
            s[(a, b)] = val;
            s[(b, a)] = val;
        }
    }
    s
}

/// Analytic s/p Slater overlaps in the diatomic frame (MOPAC `diat2`).
/// Fills `s[shell_a][shell_b][σ/σπ/ππ]` (0-based shells).
#[allow(clippy::too_many_arguments)]
fn diat2_overlaps(
    na: u8,
    esa: f64,
    epa: f64,
    nb: u8,
    esb: f64,
    epb: f64,
    r_ang: f64,
    ii: usize,
    s: &mut [[[f64; 3]; 3]; 3],
) {
    let rab = r_ang / A0;
    let mut ab = set_ab(esa, esb, na, nb, rab, ii);
    let mut isp = ab.isp;
    let mut ips = ab.ips;
    let s11;
    let mut s22a = 0.0;
    let mut s22b = 0.0;
    match ii {
        1 => {
            let (sa, sb, a, b) = (ab.sa, ab.sb, ab.a, ab.b);
            s11 = 0.25 * ((sa * sb * rab * rab).powi(3)).sqrt() * (a[2] * b[0] - b[2] * a[0]);
        }
        2 => {
            let rab4 = rab.powi(4) * 0.125;
            let w = (ab.sa.powi(3) * ab.sb.powi(5)).sqrt() * rab4;
            s11 = w
                * (1.0 / 3.0_f64.sqrt())
                * (ab.a[3] * ab.b[0] - ab.b[3] * ab.a[0] + ab.a[2] * ab.b[1] - ab.b[2] * ab.a[1]);
            if na > 1 {
                ab = set_ab(epa, esb, na, nb, rab, ii);
            }
            if nb > 1 {
                ab = set_ab(esa, epb, na, nb, rab, ii);
            }
            isp = ab.isp;
            ips = ab.ips;
            let w = (ab.sa.powi(3) * ab.sb.powi(5)).sqrt() * rab4;
            s[isp - 1][ips - 1][0] =
                w * (ab.a[2] * ab.b[0] - ab.b[2] * ab.a[0] + ab.a[3] * ab.b[1] - ab.b[3] * ab.a[1]);
        }
        3 => {
            let rab4 = rab.powi(5) * 0.0625;
            let w = (ab.sa.powi(3) * ab.sb.powi(7) / 22.5).sqrt() * rab4;
            s11 = w
                * (ab.a[4] * ab.b[0] - ab.b[4] * ab.a[0]
                    + (ab.a[3] * ab.b[1] - ab.b[3] * ab.a[1]) * 2.0);
            if na > 1 {
                ab = set_ab(epa, esb, na, nb, rab, ii);
            }
            if nb > 1 {
                ab = set_ab(esa, epb, na, nb, rab, ii);
            }
            isp = ab.isp;
            ips = ab.ips;
            let w = (ab.sa.powi(3) * ab.sb.powi(7) / 7.5).sqrt() * rab4;
            s[isp - 1][ips - 1][0] = w
                * (ab.a[3] * (ab.b[0] + ab.b[2]) - ab.b[3] * (ab.a[0] + ab.a[2])
                    + ab.b[1] * (ab.a[2] + ab.a[4])
                    - ab.a[1] * (ab.b[2] + ab.b[4]));
        }
        4 => {
            let rab4 = rab.powi(5) * 0.0625;
            let w = ((ab.sa * ab.sb).powi(5)).sqrt() * rab4;
            s11 = w * (ab.a[4] * ab.b[0] + ab.b[4] * ab.a[0] - 2.0 * ab.a[2] * ab.b[2]) / 3.0;
            // s-p
            ab = set_ab(esa, epb, na, nb, rab, ii);
            if na > nb {
                ab = set_ab(epa, esb, na, nb, rab, ii);
            }
            isp = ab.isp;
            ips = ab.ips;
            let w = ((ab.sa * ab.sb).powi(5)).sqrt() * rab4;
            let rt3 = 1.0 / 3.0_f64.sqrt();
            let d = ab.a[3] * (ab.b[0] - ab.b[2]) - ab.a[1] * (ab.b[2] - ab.b[4]);
            let e = ab.b[3] * (ab.a[0] - ab.a[2]) - ab.b[1] * (ab.a[2] - ab.a[4]);
            s[isp - 1][ips - 1][0] = w * rt3 * (d + e);
            // p-s
            ab = set_ab(epa, esb, na, nb, rab, ii);
            if na > nb {
                ab = set_ab(esa, epb, na, nb, rab, ii);
            }
            isp = ab.isp;
            ips = ab.ips;
            let w = ((ab.sa * ab.sb).powi(5)).sqrt() * rab4;
            let d = ab.a[3] * (ab.b[0] - ab.b[2]) - ab.a[1] * (ab.b[2] - ab.b[4]);
            let e = ab.b[3] * (ab.a[0] - ab.a[2]) - ab.b[1] * (ab.a[2] - ab.a[4]);
            s[ips - 1][isp - 1][0] = w * rt3 * (d - e);
            // p-p
            ab = set_ab(epa, epb, na, nb, rab, ii);
            let w = ((ab.sa * ab.sb).powi(5)).sqrt() * rab4;
            s22a = -w * (ab.b[2] * (ab.a[4] + ab.a[0]) - ab.a[2] * (ab.b[4] + ab.b[0]));
            s22b = 0.5
                * w
                * (ab.a[4] * (ab.b[0] - ab.b[2])
                    - ab.b[4] * (ab.a[0] - ab.a[2])
                    - ab.a[2] * ab.b[0]
                    + ab.b[2] * ab.a[0]);
        }
        5 => {
            let rab6 = rab.powi(6) * 0.03125 / 7.5_f64.sqrt();
            let w = (ab.sa.powi(5) * ab.sb.powi(7)).sqrt() * rab6;
            let rt3 = 1.0 / 3.0_f64.sqrt();
            s11 = w
                * (ab.a[5] * ab.b[0] + ab.a[4] * ab.b[1]
                    - 2.0 * (ab.a[3] * ab.b[2] + ab.a[2] * ab.b[3])
                    + ab.a[1] * ab.b[4]
                    + ab.a[0] * ab.b[5])
                / 3.0;
            ab = set_ab(esa, epb, na, nb, rab, ii);
            if na > nb {
                ab = set_ab(epa, esb, na, nb, rab, ii);
            }
            isp = ab.isp;
            ips = ab.ips;
            let w = (ab.sa.powi(5) * ab.sb.powi(7)).sqrt() * rab6;
            s[isp - 1][ips - 1][0] = w
                * rt3
                * (ab.a[5] * ab.b[1] + ab.a[4] * ab.b[0]
                    - 2.0 * (ab.a[3] * ab.b[3] + ab.a[2] * ab.b[2])
                    + ab.a[1] * ab.b[5]
                    + ab.a[0] * ab.b[4]);
            ab = set_ab(epa, esb, na, nb, rab, ii);
            if na > nb {
                ab = set_ab(esa, epb, na, nb, rab, ii);
            }
            isp = ab.isp;
            ips = ab.ips;
            let w = (ab.sa.powi(5) * ab.sb.powi(7)).sqrt() * rab6;
            s[ips - 1][isp - 1][0] = -w
                * rt3
                * (ab.a[4] * (2.0 * ab.b[2] - ab.b[0])
                    - ab.b[4] * (2.0 * ab.a[2] - ab.a[0])
                    - ab.a[1] * (ab.b[5] - 2.0 * ab.b[3])
                    + ab.b[1] * (ab.a[5] - 2.0 * ab.a[3]));
            ab = set_ab(epa, epb, na, nb, rab, ii);
            let w = (ab.sa.powi(5) * ab.sb.powi(7)).sqrt() * rab6;
            s22a = -w
                * (ab.b[3] * (ab.a[0] + ab.a[4]) - ab.a[3] * (ab.b[0] + ab.b[4])
                    + ab.b[2] * (ab.a[1] + ab.a[5])
                    - ab.a[2] * (ab.b[1] + ab.b[5]));
            s22b = 0.5
                * w
                * (ab.a[5] * (ab.b[0] - ab.b[2]) - ab.b[5] * (ab.a[0] - ab.a[2])
                    + ab.a[4] * (ab.b[1] - ab.b[3])
                    - ab.b[4] * (ab.a[1] - ab.a[3])
                    - ab.a[3] * ab.b[0]
                    + ab.b[3] * ab.a[0]
                    - ab.a[2] * ab.b[1]
                    + ab.b[2] * ab.a[1]);
        }
        _ => {
            // case 6: third row - third row
            let rab4 = rab.powi(7) / 480.0;
            let w = ((ab.sa * ab.sb).powi(7)).sqrt() * rab4;
            let rt3 = 1.0 / 3.0_f64.sqrt();
            s11 = w
                * (ab.a[6] * ab.b[0]
                    - 3.0 * (ab.a[4] * ab.b[2] - ab.a[2] * ab.b[4])
                    - ab.a[0] * ab.b[6])
                / 3.0;
            ab = set_ab(esa, epb, na, nb, rab, ii);
            if na > nb {
                ab = set_ab(epa, esb, na, nb, rab, ii);
            }
            isp = ab.isp;
            ips = ab.ips;
            let w = ((ab.sa * ab.sb).powi(7)).sqrt() * rab4;
            let d = ab.a[5] * (ab.b[0] - ab.b[2]) - 2.0 * ab.a[3] * (ab.b[2] - ab.b[4])
                + ab.a[1] * (ab.b[4] - ab.b[6]);
            let e = ab.b[5] * (ab.a[0] - ab.a[2]) - 2.0 * ab.b[3] * (ab.a[2] - ab.a[4])
                + ab.b[1] * (ab.a[4] - ab.a[6]);
            s[isp - 1][ips - 1][0] = w * rt3 * (d - e);
            ab = set_ab(epa, esb, na, nb, rab, ii);
            if na > nb {
                ab = set_ab(esa, epb, na, nb, rab, ii);
            }
            isp = ab.isp;
            ips = ab.ips;
            let w = ((ab.sa * ab.sb).powi(7)).sqrt() * rab4;
            let d = ab.a[5] * (ab.b[0] - ab.b[2]) - 2.0 * ab.a[3] * (ab.b[2] - ab.b[4])
                + ab.a[1] * (ab.b[4] - ab.b[6]);
            let e = ab.b[5] * (ab.a[0] - ab.a[2]) - 2.0 * ab.b[3] * (ab.a[2] - ab.a[4])
                + ab.b[1] * (ab.a[4] - ab.a[6]);
            s[ips - 1][isp - 1][0] = -w * rt3 * ((-d) - e);
            ab = set_ab(epa, epb, na, nb, rab, ii);
            let w = ((ab.sa * ab.sb).powi(7)).sqrt() * rab4;
            let d = ab.a[2] * (ab.b[6] + ab.b[2] + ab.b[2])
                - ab.a[4] * (ab.b[0] + ab.b[4] + ab.b[4])
                - ab.b[4] * ab.a[0]
                + ab.a[6] * ab.b[2];
            s22a = -w * d;
            let d = ab.a[6] * (ab.b[0] - ab.b[2]) + ab.b[6] * (ab.a[0] - ab.a[2]);
            let e = ab.a[4] * (ab.b[4] - ab.b[2] - ab.b[0])
                + ab.b[4] * (ab.a[4] - ab.a[2] - ab.a[0])
                + 2.0 * ab.a[2] * ab.b[2];
            s22b = 0.5 * w * (d + e);
        }
    }
    s[0][0][0] = s11;
    s[1][1][0] = s22a;
    s[1][1][1] = s22b;
}

/// Overlap between one orbital of atom A and one of atom B (sp basis),
/// following MOPAC `diat` + `diat2` + `coe`.
#[allow(clippy::too_many_arguments)]
fn overlap_sp(
    pa: &Pm3Params,
    za: u8,
    pb: &Pm3Params,
    zb: u8,
    r_ang: f64,
    oa: usize,
    ob: usize,
    pos_a: &[f64; 3],
    pos_b: &[f64; 3],
) -> f64 {
    // Direction cosines (MOPAC `coe`): c(i, m, n) with 0-based i.
    let x2 = pos_b[0] - pos_a[0];
    let y2 = pos_b[1] - pos_a[1];
    let z2 = pos_b[2] - pos_a[2];
    let xy = (x2 * x2 + y2 * y2).sqrt();
    let r = (xy * xy + z2 * z2).sqrt();
    let (ca, cb, sa, sb) = if xy >= 1e-10 {
        (x2 / xy, z2 / r, y2 / xy, xy / r)
    } else if z2 <= 0.0 {
        if z2 != 0.0 {
            (-1.0, -1.0, 0.0, 0.0)
        } else {
            (0.0, 0.0, 0.0, 0.0)
        }
    } else {
        (1.0, 1.0, 0.0, 0.0)
    };
    let mut c = [[[0.0f64; 5]; 5]; 3];
    c[0][2][2] = 1.0;
    if is_heavy(za) || is_heavy(zb) {
        c[1][1][1] = ca;
        c[1][1][2] = sa * sb;
        c[1][1][3] = sa * cb;
        c[1][2][1] = 0.0;
        c[1][2][2] = cb;
        c[1][2][3] = -sb;
        c[1][3][1] = -sa;
        c[1][3][2] = ca * sb;
        c[1][3][3] = ca * cb;
    }
    let ii = bond_type(za, zb);
    let mut s = [[[0.0f64; 3]; 3]; 3];
    diat2_overlaps(
        za, pa.zeta_s, pa.zeta_p, zb, pb.zeta_s, pb.zeta_p, r_ang, ii, &mut s,
    );
    let ia = if is_heavy(za) { 2 } else { 1 };
    let ib = if is_heavy(zb) { 2 } else { 1 };
    let shell_m = |o: usize| -> (usize, usize) {
        match o {
            0 => (1, 3),
            1 => (2, 4),
            2 => (2, 2),
            _ => (2, 3),
        }
    };
    let (i, k) = shell_m(oa);
    let (j, l) = shell_m(ob);
    if i > ia || j > ib {
        return 0.0;
    }
    let aa = if j == 2 { -1.0 } else { 1.0 };
    let bb = if j == 3 { -1.0 } else { 1.0 };
    let s1 = s[i - 1][j - 1][0];
    let s2 = s[i - 1][j - 1][1];
    let s3 = s[i - 1][j - 1][2];
    s1 * (c[i - 1][k - 1][2] * c[j - 1][l - 1][2]) * aa
        + s2 * (c[i - 1][k - 1][3] * c[j - 1][l - 1][3] + c[i - 1][k - 1][1] * c[j - 1][l - 1][1])
            * bb
        + s3 * (c[i - 1][k - 1][4] * c[j - 1][l - 1][4] + c[i - 1][k - 1][0] * c[j - 1][l - 1][0])
}

/// MOPAC bond-type code (`inmb`/`iii`): element-based, not orbital-based.
fn bond_type(za: u8, zb: u8) -> usize {
    const INMB: [usize; 17] = [1, 0, 2, 2, 3, 4, 5, 6, 7, 0, 8, 8, 8, 9, 10, 11, 12];
    const III: [usize; 78] = [
        1, 2, 4, 2, 4, 4, 2, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 3, 5,
        5, 5, 5, 5, 5, 6, 3, 5, 5, 5, 5, 5, 5, 6, 6, 3, 5, 5, 5, 5, 5, 5, 6, 6, 6, 3, 5, 5, 5, 5,
        5, 5, 6, 6, 6, 6, 3, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6,
    ];
    let inmb = |z: u8| -> usize {
        INMB.get((z as usize).saturating_sub(1))
            .copied()
            .unwrap_or(12)
    };
    let (ia, ib) = (inmb(za), inmb(zb));
    let (jmax, jmin) = (ia.max(ib), ia.min(ib));
    let nbond = jmax * (jmax - 1) / 2 + jmin;
    III.get(nbond.saturating_sub(1)).copied().unwrap_or(4)
}

// ─── Two-center two-electron integrals (MOPAC reppd, sp subset) ─────────────

/// 22 local-frame two-electron integrals + the `gab` used by the core-core term.
fn reppd(
    pa: &Pm3Params,
    pb: &Pm3Params,
    da: &Pm3Derived,
    db: &Pm3Derived,
    r_ang: f64,
) -> ([f64; 22], f64) {
    let r = r_ang / A0;
    let si = is_heavy(pa.z);
    let sj = is_heavy(pb.z);
    let (ev1, ev2, ev3, ev4) = (EV / 2.0, EV / 4.0, EV / 8.0, EV / 16.0);
    let mut ri = [0.0f64; 22];

    // <ss|ss> and the core-core `gab` (MOPAC keeps them consistent).
    let aee_gab = da.po9 + db.po9;
    let gab = EV / (r * r + aee_gab * aee_gab).sqrt();
    let aee = 0.5 / da.am + 0.5 / db.am;
    let aee = aee * aee;

    if !si && !sj {
        // H–H
        ri[0] = EV / (r * r + aee).sqrt();
    } else if si && !sj {
        // heavy–H
        let da_ = da.dd;
        let qa = da.qq * 2.0;
        let ade = (0.5 / da.ad + 0.5 / db.am).powi(2);
        let aqe = (0.5 / da.aq + 0.5 / db.am).powi(2);
        let rsq = r * r;
        let arg = [
            rsq + aee,
            (r + da_).powi(2) + ade,
            (r - da_).powi(2) + ade,
            (r + qa).powi(2) + aqe,
            (r - qa).powi(2) + aqe,
            rsq + aqe,
            rsq + aqe + qa * qa,
        ];
        let sqr: Vec<f64> = arg.iter().map(|x| x.sqrt()).collect();
        let ee = EV / sqr[0];
        ri[0] = ee;
        ri[1] = ev1 / sqr[1] - ev1 / sqr[2];
        ri[2] = ee + ev2 / sqr[3] + ev2 / sqr[4] - ev1 / sqr[5];
        ri[3] = ee + ev1 / sqr[6] - ev1 / sqr[5];
    } else if !si && sj {
        // H–heavy
        let db_ = db.dd;
        let qb = db.qq * 2.0;
        let aed = (0.5 / da.am + 0.5 / db.ad).powi(2);
        let aeq = (0.5 / da.am + 0.5 / db.aq).powi(2);
        let rsq = r * r;
        let arg = [
            rsq + aee,
            (r - db_).powi(2) + aed,
            (r + db_).powi(2) + aed,
            (r - qb).powi(2) + aeq,
            (r + qb).powi(2) + aeq,
            rsq + aeq,
            rsq + aeq + qb * qb,
        ];
        let sqr: Vec<f64> = arg.iter().map(|x| x.sqrt()).collect();
        let ee = EV / sqr[0];
        ri[0] = ee;
        ri[4] = ev1 / sqr[1] - ev1 / sqr[2];
        ri[10] = ee + ev2 / sqr[3] + ev2 / sqr[4] - ev1 / sqr[5];
        ri[11] = ee + ev1 / sqr[6] - ev1 / sqr[5];
    } else {
        // heavy–heavy
        let da_ = da.dd;
        let db_ = db.dd;
        let qa = da.qq * 2.0;
        let qb = db.qq * 2.0;
        let ade = (0.5 / da.ad + 0.5 / db.am).powi(2);
        let aqe = (0.5 / da.aq + 0.5 / db.am).powi(2);
        let aed = (0.5 / da.am + 0.5 / db.ad).powi(2);
        let aeq = (0.5 / da.am + 0.5 / db.aq).powi(2);
        let axx = (0.5 / da.ad + 0.5 / db.ad).powi(2);
        let adq = (0.5 / da.ad + 0.5 / db.aq).powi(2);
        let aqd = (0.5 / da.aq + 0.5 / db.ad).powi(2);
        let aqq = (0.5 / da.aq + 0.5 / db.aq).powi(2);
        let rsq = r * r;
        let mut arg = [0.0f64; 72];
        arg[0] = rsq + aee;
        arg[1] = (r + da_).powi(2) + ade;
        arg[2] = (r - da_).powi(2) + ade;
        arg[3] = (r - qa).powi(2) + aqe;
        arg[4] = (r + qa).powi(2) + aqe;
        arg[5] = rsq + aqe;
        arg[6] = rsq + aqe + qa * qa;
        arg[7] = (r - db_).powi(2) + aed;
        arg[8] = (r + db_).powi(2) + aed;
        arg[9] = (r - qb).powi(2) + aeq;
        arg[10] = (r + qb).powi(2) + aeq;
        arg[11] = rsq + aeq;
        arg[12] = rsq + aeq + qb * qb;
        arg[13] = rsq + axx + (da_ - db_).powi(2);
        arg[14] = rsq + axx + (da_ + db_).powi(2);
        arg[15] = (r + da_ - db_).powi(2) + axx;
        arg[16] = (r - da_ + db_).powi(2) + axx;
        arg[17] = (r - da_ - db_).powi(2) + axx;
        arg[18] = (r + da_ + db_).powi(2) + axx;
        arg[19] = (r + da_).powi(2) + adq;
        arg[20] = arg[19] + qb * qb;
        arg[21] = (r - da_).powi(2) + adq;
        arg[22] = arg[21] + qb * qb;
        arg[23] = (r - db_).powi(2) + aqd;
        arg[24] = arg[23] + qa * qa;
        arg[25] = (r + db_).powi(2) + aqd;
        arg[26] = arg[25] + qa * qa;
        arg[27] = (r + da_ - qb).powi(2) + adq;
        arg[28] = (r - da_ - qb).powi(2) + adq;
        arg[29] = (r + da_ + qb).powi(2) + adq;
        arg[30] = (r - da_ + qb).powi(2) + adq;
        arg[31] = (r + qa - db_).powi(2) + aqd;
        arg[32] = (r + qa + db_).powi(2) + aqd;
        arg[33] = (r - qa - db_).powi(2) + aqd;
        arg[34] = (r - qa + db_).powi(2) + aqd;
        arg[35] = rsq + aqq;
        arg[36] = arg[35] + (qa - qb).powi(2);
        arg[37] = arg[35] + (qa + qb).powi(2);
        arg[38] = arg[35] + qa * qa;
        arg[39] = arg[35] + qb * qb;
        arg[40] = arg[38] + qb * qb;
        arg[41] = (r - qb).powi(2) + aqq;
        arg[42] = arg[41] + qa * qa;
        arg[43] = (r + qb).powi(2) + aqq;
        arg[44] = arg[43] + qa * qa;
        arg[45] = (r + qa).powi(2) + aqq;
        arg[46] = arg[45] + qb * qb;
        arg[47] = (r - qa).powi(2) + aqq;
        arg[48] = arg[47] + qb * qb;
        arg[49] = (r + qa - qb).powi(2) + aqq;
        arg[50] = (r + qa + qb).powi(2) + aqq;
        arg[51] = (r - qa - qb).powi(2) + aqq;
        arg[52] = (r - qa + qb).powi(2) + aqq;
        arg[53] = (da_ - qb).powi(2) + (r - qb).powi(2) + adq;
        arg[54] = (da_ - qb).powi(2) + (r + qb).powi(2) + adq;
        arg[55] = (da_ + qb).powi(2) + (r - qb).powi(2) + adq;
        arg[56] = (da_ + qb).powi(2) + (r + qb).powi(2) + adq;
        arg[57] = (r + qa).powi(2) + (qa - db_).powi(2) + aqd;
        arg[58] = (r - qa).powi(2) + (qa - db_).powi(2) + aqd;
        arg[59] = (r + qa).powi(2) + (qa + db_).powi(2) + aqd;
        arg[60] = (r - qa).powi(2) + (qa + db_).powi(2) + aqd;
        arg[61] = arg[35] + 2.0 * (qa - qb).powi(2);
        arg[62] = arg[35] + 2.0 * (qa + qb).powi(2);
        arg[63] = arg[35] + 2.0 * (qa * qa + qb * qb);
        arg[64] = (r + qa - qb).powi(2) + (qa - qb).powi(2) + aqq;
        arg[65] = (r + qa - qb).powi(2) + (qa + qb).powi(2) + aqq;
        arg[66] = (r + qa + qb).powi(2) + (qa - qb).powi(2) + aqq;
        arg[67] = (r + qa + qb).powi(2) + (qa + qb).powi(2) + aqq;
        arg[68] = (r - qa - qb).powi(2) + (qa - qb).powi(2) + aqq;
        arg[69] = (r - qa - qb).powi(2) + (qa + qb).powi(2) + aqq;
        arg[70] = (r - qa + qb).powi(2) + (qa - qb).powi(2) + aqq;
        arg[71] = (r - qa + qb).powi(2) + (qa + qb).powi(2) + aqq;
        let sqr: Vec<f64> = arg.iter().map(|x| x.sqrt()).collect();
        let ee = EV / sqr[0];
        let dze = -ev1 / sqr[1] + ev1 / sqr[2];
        let qzze = ev2 / sqr[3] + ev2 / sqr[4] - ev1 / sqr[5];
        let qxxe = ev1 / sqr[6] - ev1 / sqr[5];
        let edz = -ev1 / sqr[7] + ev1 / sqr[8];
        let eqzz = ev2 / sqr[9] + ev2 / sqr[10] - ev1 / sqr[11];
        let eqxx = ev1 / sqr[12] - ev1 / sqr[11];
        let dxdx = ev1 / sqr[13] - ev1 / sqr[14];
        let dzdz = ev2 / sqr[15] + ev2 / sqr[16] - ev2 / sqr[17] - ev2 / sqr[18];
        let dzqxx = ev2 / sqr[19] - ev2 / sqr[20] - ev2 / sqr[21] + ev2 / sqr[22];
        let qxxdz = ev2 / sqr[23] - ev2 / sqr[24] - ev2 / sqr[25] + ev2 / sqr[26];
        let dzqzz = -ev3 / sqr[27] + ev3 / sqr[28] - ev3 / sqr[29] + ev3 / sqr[30] - ev2 / sqr[21]
            + ev2 / sqr[19];
        let qzzdz = -ev3 / sqr[31] + ev3 / sqr[32] - ev3 / sqr[33] + ev3 / sqr[34] + ev2 / sqr[23]
            - ev2 / sqr[25];
        let qxxqxx = ev3 / sqr[36] + ev3 / sqr[37] - ev2 / sqr[38] - ev2 / sqr[39] + ev2 / sqr[35];
        let qxxqyy = ev2 / sqr[40] - ev2 / sqr[38] - ev2 / sqr[39] + ev2 / sqr[35];
        let qxxqzz = ev3 / sqr[42] + ev3 / sqr[44] - ev3 / sqr[41] - ev3 / sqr[43] - ev2 / sqr[38]
            + ev2 / sqr[35];
        let qzzqxx = ev3 / sqr[46] + ev3 / sqr[48] - ev3 / sqr[45] - ev3 / sqr[47] - ev2 / sqr[39]
            + ev2 / sqr[35];
        let qzzqzz = ev4 / sqr[49] + ev4 / sqr[50] + ev4 / sqr[51] + ev4 / sqr[52]
            - ev3 / sqr[47]
            - ev3 / sqr[45]
            - ev3 / sqr[41]
            - ev3 / sqr[43]
            + ev2 / sqr[35];
        let dxqxz = -ev2 / sqr[53] + ev2 / sqr[54] + ev2 / sqr[55] - ev2 / sqr[56];
        let qxzdx = -ev2 / sqr[57] + ev2 / sqr[58] + ev2 / sqr[59] - ev2 / sqr[60];
        let qxzqxz = ev3 / sqr[64] - ev3 / sqr[66] - ev3 / sqr[68] + ev3 / sqr[70] - ev3 / sqr[65]
            + ev3 / sqr[67]
            + ev3 / sqr[69]
            - ev3 / sqr[71];
        ri[0] = ee;
        ri[1] = -dze;
        ri[2] = ee + qzze;
        ri[3] = ee + qxxe;
        ri[4] = -edz;
        ri[5] = dzdz;
        ri[6] = dxdx;
        ri[7] = -edz - qzzdz;
        ri[8] = -edz - qxxdz;
        ri[9] = -qxzdx;
        ri[10] = ee + eqzz;
        ri[11] = ee + eqxx;
        ri[12] = -dze - dzqzz;
        ri[13] = -dze - dzqxx;
        ri[14] = -dxqxz;
        ri[15] = ee + eqzz + qzze + qzzqzz;
        ri[16] = ee + eqzz + qxxe + qxxqzz;
        ri[17] = ee + eqxx + qzze + qzzqxx;
        ri[18] = ee + eqxx + qxxe + qxxqxx;
        ri[19] = qxzqxz;
        ri[20] = ee + eqxx + qxxe + qxxqyy;
        ri[21] = 0.5 * (qxxqxx - qxxqyy);
    }
    // MOPAC `nri` sign pattern.
    const NRI: [f64; 22] = [
        1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    ];
    for (v, s) in ri.iter_mut().zip(NRI.iter()) {
        *v *= s;
    }
    (ri, gab)
}

/// Direction cosines (MOPAC `rotmat`, sp block): `p[3][3]` and `pp[6][3][3]`.
fn rotmat(pos_a: &[f64; 3], pos_b: &[f64; 3]) -> ([f64; 3], [[f64; 3]; 3], [[[f64; 3]; 3]; 6]) {
    let (x11, x22, x33) = (
        pos_b[0] - pos_a[0],
        pos_b[1] - pos_a[1],
        pos_b[2] - pos_a[2],
    );
    let b = x11 * x11 + x22 * x22;
    let r = (b + x33 * x33).sqrt();
    let sqb = b.sqrt();
    let sb0 = sqb / r;
    let (ca, sa, cb, sb);
    if sb0 > 1e-7 {
        ca = x11 / sqb;
        sa = x22 / sqb;
        cb = x33 / r;
        sb = sb0;
    } else {
        sa = 0.0;
        sb = 0.0;
        if x33 < 0.0 {
            ca = -1.0;
            cb = -1.0;
        } else if x33 > 0.0 {
            ca = 1.0;
            cb = 1.0;
        } else {
            ca = 0.0;
            cb = 0.0;
        }
    }
    // p(k, m): k = local axis (1=σ, 2, 3), m = molecular p index.
    let p = [
        [ca * sb, sa * sb, cb],
        [ca * cb, sa * cb, -sb],
        [-sa, ca, 0.0],
    ];
    let mut pp = [[[0.0f64; 3]; 3]; 6];
    for k in 0..3 {
        for c in 0..3 {
            pp[0][k][k] += 0.0; // placeholder, filled below
            let _ = c;
        }
        pp[0][k][k] = p[k][0] * p[k][0];
        pp[1][k][k] = p[k][1] * p[k][1];
        pp[2][k][k] = p[k][2] * p[k][2];
        pp[3][k][k] = p[k][0] * p[k][1];
        pp[4][k][k] = p[k][0] * p[k][2];
        pp[5][k][k] = p[k][1] * p[k][2];
        for l in 0..k {
            pp[0][k][l] = 2.0 * p[k][0] * p[l][0];
            pp[1][k][l] = 2.0 * p[k][1] * p[l][1];
            pp[2][k][l] = 2.0 * p[k][2] * p[l][2];
            pp[3][k][l] = p[k][0] * p[l][1] + p[k][1] * p[l][0];
            pp[4][k][l] = p[k][0] * p[l][2] + p[k][2] * p[l][0];
            pp[5][k][l] = p[k][1] * p[l][2] + p[k][2] * p[l][1];
        }
    }
    (p[0], p, pp)
}

/// indexd convention (MOPAC `fordd`): (i, j) with i >= j, 1-based orbitals.
fn indexd(i: usize, j: usize) -> usize {
    i + 9 * (j - 1) - (j * (j - 1)) / 2
}

/// indx convention (MOPAC `fordd`): standard pair index.
fn indx(i: usize, j: usize) -> usize {
    (i * (i - 1)) / 2 + j
}

/// `met` classification (MOPAC `mndod_C` / `tx`).
const MET: [usize; 10] = [1, 2, 3, 2, 3, 3, 2, 3, 3, 3];

/// `ind2` mapping for the sp block (MOPAC `fordd`).
fn ind2(ij: usize, kl: usize) -> usize {
    const T: [(usize, usize, usize); 34] = [
        (1, 1, 1),
        (1, 2, 2),
        (1, 10, 3),
        (1, 18, 4),
        (1, 25, 5),
        (2, 1, 6),
        (2, 2, 7),
        (2, 10, 8),
        (2, 18, 9),
        (2, 25, 10),
        (10, 1, 11),
        (10, 2, 12),
        (10, 10, 13),
        (10, 18, 14),
        (10, 25, 15),
        (3, 3, 16),
        (3, 11, 17),
        (11, 3, 18),
        (11, 11, 19),
        (18, 1, 20),
        (18, 2, 21),
        (18, 10, 22),
        (18, 18, 23),
        (18, 25, 24),
        (4, 4, 25),
        (4, 12, 26),
        (12, 4, 27),
        (12, 12, 28),
        (19, 19, 29),
        (25, 1, 30),
        (25, 2, 31),
        (25, 10, 32),
        (25, 18, 33),
        (25, 25, 34),
    ];
    T.iter()
        .find(|(a, b, _)| *a == ij && *b == kl)
        .map(|(_, _, v)| *v)
        .unwrap_or(0)
}

/// Two-center two-electron block in the molecular frame, indexed by the `indx`
/// pair convention (1=s s, 2=px s, 3=px px, 4=py s, 5=py px, 6=py py, 7=pz s,
/// 8=pz px, 9=pz py, 10=pz pz) for both atoms.
fn two_center_block(
    pa: &Pm3Params,
    pb: &Pm3Params,
    da: &Pm3Derived,
    db: &Pm3Derived,
    pos_a: &[f64; 3],
    pos_b: &[f64; 3],
) -> ([[f64; 10]; 10], [[f64; 4]; 4], [[f64; 4]; 4], f64) {
    let r_ang = ((pos_a[0] - pos_b[0]).powi(2)
        + (pos_a[1] - pos_b[1]).powi(2)
        + (pos_a[2] - pos_b[2]).powi(2))
    .sqrt();
    let (_, p, pp) = rotmat(pos_a, pos_b);
    let (ri, gab) = reppd(pa, pb, da, db, r_ang);
    // MOPAC `reppd2`: rep(1..34) = ri(ipos) for the sp block.
    const IPOS: [usize; 34] = [
        1, 5, 11, 12, 12, 2, 6, 13, 14, 14, 3, 8, 16, 18, 18, 7, 15, 10, 20, 4, 9, 17, 19, 21, 7,
        15, 10, 20, 22, 4, 9, 17, 21, 19,
    ];
    let rep: Vec<f64> = IPOS.iter().map(|&k| ri[k - 1]).collect();
    let cored = spcore(pa, pb, da, db, r_ang);

    // Step 1 (MOPAC `tx`): local integrals -> v[ij][pair_on_B].
    let mut v = [[0.0f64; 10]; 26]; // first index = indexd value (max 25)
    for i1 in 1..=4usize {
        for j1 in 1..=i1 {
            let ij = indexd(i1, j1);
            for k1 in 1..=4usize {
                for l1 in 1..=k1 {
                    let kl = indexd(k1, l1);
                    let nd = ind2(ij, kl);
                    if nd == 0 {
                        continue;
                    }
                    let wrepp = rep[nd - 1]; // rep = reppd2's 34-entry mapping
                    let ll = indx(k1, l1);
                    match MET[ll - 1] {
                        1 => v[ij][0] = wrepp,
                        2 => {
                            let k = k1 - 2;
                            v[ij][1] += p[k][0] * wrepp;
                            v[ij][3] += p[k][1] * wrepp;
                            v[ij][6] += p[k][2] * wrepp;
                        }
                        3 => {
                            let k = k1 - 2;
                            let l = l1 - 2;
                            v[ij][2] += pp[0][k][l] * wrepp;
                            v[ij][5] += pp[1][k][l] * wrepp;
                            v[ij][9] += pp[2][k][l] * wrepp;
                            v[ij][4] += pp[3][k][l] * wrepp;
                            v[ij][7] += pp[4][k][l] * wrepp;
                            v[ij][8] += pp[5][k][l] * wrepp;
                        }
                        _ => {}
                    }
                }
            }
        }
    }
    // Step 2 (MOPAC `rotatd`): v -> ww[pair_on_A][pair_on_B].
    let mut ww = [[0.0f64; 10]; 10];
    for i1 in 1..=4usize {
        for j1 in 1..=i1 {
            let ij = indexd(i1, j1);
            let jj = indx(i1, j1);
            let mm = MET[jj - 1];
            for kl in 1..=10usize {
                if v[ij][kl - 1] == 0.0 {
                    continue;
                }
                let wrepp = v[ij][kl - 1];
                match mm {
                    1 => ww[0][kl - 1] = wrepp,
                    2 => {
                        for i in 1..=3usize {
                            ww[indx(i + 1, 1) - 1][kl - 1] += p[i1 - 2][i - 1] * wrepp;
                        }
                    }
                    3 => {
                        for i in 1..=3usize {
                            let cc = pp[i - 1][i1 - 2][j1 - 2];
                            ww[indx(i + 1, i + 1) - 1][kl - 1] += cc * wrepp;
                            for j in 1..i {
                                let cc = pp[1 + i + j - 1][i1 - 2][j1 - 2];
                                ww[indx(i + 1, j + 1) - 1][kl - 1] += cc * wrepp;
                            }
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    // Electron-nuclear attraction terms (MOPAC `elenuc` placement):
    // lower-triangular 4x4 blocks over (s, px, py, pz).
    let mut e1b = [[0.0f64; 4]; 4];
    let mut e2a = [[0.0f64; 4]; 4];
    e1b[0][0] = cored[0][0];
    e2a[0][0] = cored[0][1];
    for k in 1..=3usize {
        e1b[k][0] = p[0][k - 1] * cored[1][0];
        e2a[k][0] = p[0][k - 1] * cored[1][1];
    }
    for k in 1..=3usize {
        for l in 1..=k {
            let ipp = indpp(k, l);
            e1b[k][l] = cored[2][0] * pp[ipp - 1][0][0]
                + cored[3][0] * (pp[ipp - 1][1][1] + pp[ipp - 1][2][2]);
            e2a[k][l] = cored[2][1] * pp[ipp - 1][0][0]
                + cored[3][1] * (pp[ipp - 1][1][1] + pp[ipp - 1][2][2]);
        }
    }
    let enuc = ccrep(pa, pb, r_ang, gab);
    (ww, e1b, e2a, enuc)
}

/// `indpp` (MOPAC `rotmat` comment).
fn indpp(k: usize, l: usize) -> usize {
    const INDPP: [[usize; 3]; 3] = [[1, 4, 5], [4, 2, 6], [5, 6, 3]];
    INDPP[k - 1][l - 1]
}

/// Electron-nuclear attraction integrals (MOPAC `spcore`), per atom side.
fn spcore(
    pa: &Pm3Params,
    pb: &Pm3Params,
    da: &Pm3Derived,
    db: &Pm3Derived,
    r_ang: f64,
) -> [[f64; 2]; 4] {
    let r = r_ang / A0;
    let mut core = [[0.0f64; 2]; 4];
    let r2 = r * r;
    let aci = da.po9;
    let acj = db.po9;
    let ssi = (aci + db.po1).powi(2);
    let ssj = (acj + da.po1).powi(2);
    core[0][0] = -pb.core_charge * EV / (r2 + ssj).sqrt();
    core[0][1] = -pa.core_charge * EV / (r2 + ssi).sqrt();
    if pa.z >= 3 || pb.z >= 3 {
        if pa.z >= 3 {
            let ppj = (acj + da.po7).powi(2);
            let da_ = da.dd;
            let qa = da.qq;
            let twoqa = 2.0 * qa;
            let adj = (da.po2 + acj).powi(2);
            let aqj = (da.po3 + acj).powi(2);
            let xj = [
                r2 + ppj,
                r2 + aqj,
                (r + da_).powi(2) + adj,
                (r - da_).powi(2) + adj,
                (r - twoqa).powi(2) + aqj,
                (r + twoqa).powi(2) + aqj,
                r2 + twoqa * twoqa + aqj,
            ];
            let pxy = [1.0, -0.5, -0.5, 0.5, 0.25, 0.25, 0.5];
            let xj: Vec<f64> = xj
                .iter()
                .zip(pxy.iter())
                .map(|(x, w)| w / x.sqrt())
                .collect();
            let aj2 = (xj[2] + xj[3]) * EV;
            let aj3 = (xj[0] + xj[1] + xj[4] + xj[5]) * EV;
            let aj4 = (xj[0] + xj[1] + xj[6]) * EV;
            core[1][0] = -pb.core_charge * aj2;
            core[2][0] = -pb.core_charge * aj3;
            core[3][0] = -pb.core_charge * aj4;
        }
        if pb.z >= 3 {
            let ppi = (aci + db.po7).powi(2);
            let db_ = db.dd;
            let qb = db.qq;
            let twoqb = 2.0 * qb;
            let adi = (db.po2 + aci).powi(2);
            let aqi = (db.po3 + aci).powi(2);
            let xi = [
                r2 + ppi,
                r2 + aqi,
                (r + db_).powi(2) + adi,
                (r - db_).powi(2) + adi,
                (r - twoqb).powi(2) + aqi,
                (r + twoqb).powi(2) + aqi,
                r2 + twoqb * twoqb + aqi,
            ];
            let pxy = [1.0, -0.5, -0.5, 0.5, 0.25, 0.25, 0.5];
            let xi: Vec<f64> = xi
                .iter()
                .zip(pxy.iter())
                .map(|(x, w)| w / x.sqrt())
                .collect();
            let ai2 = -(xi[2] + xi[3]) * EV;
            let ai3 = (xi[0] + xi[1] + xi[4] + xi[5]) * EV;
            let ai4 = (xi[0] + xi[1] + xi[6]) * EV;
            core[1][1] = -pa.core_charge * ai2;
            core[2][1] = -pa.core_charge * ai3;
            core[3][1] = -pa.core_charge * ai4;
        }
    }
    core
}

/// Core-core repulsion with PM3 Gaussians (MOPAC `ccrep`, PM3 branch).
fn ccrep(pa: &Pm3Params, pb: &Pm3Params, r_ang: f64, gab: f64) -> f64 {
    let r = r_ang;
    let enuc = pa.core_charge * pb.core_charge * gab;
    let mut enuclr;
    let n_gauss;
    // PM3 has no bond parameters (`xfac`) for these pairs: exponential branch.
    let eni = (-pa.alpha * r).exp();
    let enj = (-pb.alpha * r).exp();
    let mut scale0 = eni + enj;
    // MOPAC's N-H/O-H correction (the `nt == 8 or 9` block in ccrep).
    let nt = pa.z as u32 + pb.z as u32;
    if nt == 8 || nt == 9 {
        if pa.z == 7 || pa.z == 8 {
            scale0 += (r - 1.0) * eni;
        }
        if pb.z == 7 || pb.z == 8 {
            scale0 += (r - 1.0) * enj;
        }
    }
    enuclr = scale0 * enuc.abs() + enuc;
    n_gauss = 4;
    let mut scale = 0.0;
    for ig in 0..n_gauss {
        if let Some(&(g1, g2, g3)) = pa.gaussians.get(ig) {
            if g1.abs() > 0.0 {
                let ax = g2 * (r - g3).powi(2);
                if ax <= 25.0 {
                    scale += pa.core_charge * pb.core_charge / r * g1 * (-ax).exp();
                }
            }
        }
        if let Some(&(g1, g2, g3)) = pb.gaussians.get(ig) {
            if g1.abs() > 0.0 {
                let ax = g2 * (r - g3).powi(2);
                if ax <= 25.0 {
                    scale += pa.core_charge * pb.core_charge / r * g1 * (-ax).exp();
                }
            }
        }
    }
    enuclr += scale;
    enuclr
}

// ─── One-electron matrix, Fock build, SCF and heat of formation ─────────────

/// Basis layout: (atom index, orbital type 0=s,1=px,2=py,3=pz).
fn basis_layout(elements: &[u8]) -> Vec<(usize, usize)> {
    let mut basis = Vec::new();
    for (i, &z) in elements.iter().enumerate() {
        basis.push((i, 0));
        if is_heavy(z) {
            for k in 1..4 {
                basis.push((i, k));
            }
        }
    }
    basis
}

/// Pair index (`indx` convention) from two orbital types.
fn pair_index(o1: usize, o2: usize) -> usize {
    let (i, j) = if o1 >= o2 {
        (o1 + 1, o2 + 1)
    } else {
        (o2 + 1, o1 + 1)
    };
    indx(i, j) - 1
}

/// Eigenproblem as MOPAC solves it: the AO-basis Fock is diagonalized
/// directly (`eigenvectors_LAPACK`) and the density uses `P = 2 c cᵀ`.
fn solve_generalized(fock: &DMatrix<f64>, _overlap: &DMatrix<f64>) -> (Vec<f64>, DMatrix<f64>) {
    let n = fock.nrows();
    let eig = fock.clone().symmetric_eigen();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| eig.eigenvalues[a].partial_cmp(&eig.eigenvalues[b]).unwrap());
    let mut energies = vec![0.0; n];
    let mut coefficients = DMatrix::zeros(n, n);
    for (col, &idx) in order.iter().enumerate() {
        energies[col] = eig.eigenvalues[idx];
        for row in 0..n {
            coefficients[(row, col)] = eig.eigenvectors[(row, idx)];
        }
    }
    (energies, coefficients)
}

/// Result of the faithful solver.
pub struct FaithfulResult {
    pub orbital_energies: Vec<f64>,
    pub electronic_energy: f64,
    pub nuclear_repulsion: f64,
    pub total_energy: f64,
    pub heat_of_formation: f64,
    pub homo_energy: f64,
    pub lumo_energy: f64,
    pub gap: f64,
    pub n_basis: usize,
    pub n_electrons: usize,
    pub mulliken_charges: Vec<f64>,
    pub scf_iterations: usize,
    pub converged: bool,
    /// Converged Fock matrix (validation/debug).
    pub fock: DMatrix<f64>,
    /// Converged density matrix (validation/debug).
    pub density: DMatrix<f64>,
}

/// Two-center integral block for one atom pair.
struct PairBlock {
    a: usize,
    b: usize,
    ww: [[f64; 10]; 10],
}

/// Build the NDDO G matrix for a density (MOPAC `fock1` + `fock2`).
fn build_g_matrix(
    params: &[&Pm3Params],
    basis: &[(usize, usize)],
    blocks: &[PairBlock],
    density: &DMatrix<f64>,
) -> DMatrix<f64> {
    let n_basis = basis.len();

    let mut g = DMatrix::zeros(n_basis, n_basis);
    // One-center contributions (MOPAC `fock1`, exact):
    //   G(i,j) += Σ_{k,l} [ P(k,l)<i,j|k,l> - ½ P(k,l)<i,k|j,l> ].
    for (atom_idx, &p) in params.iter().enumerate() {
        let orbs: Vec<usize> = (0..n_basis).filter(|&i| basis[i].0 == atom_idx).collect();
        let mut one = [[0.0f64; 10]; 10];
        let pair = |o1: usize, o2: usize| pair_index(o1, o2);
        one[0][0] = p.gss;
        for k in 1..=3usize {
            one[0][pair(k, k)] = p.gsp;
            one[pair(k, k)][0] = p.gsp;
            one[pair(k, k)][pair(k, k)] = p.gpp;
        }
        for k in 1..=3usize {
            // <s p_k|s p_k> = hsp.
            one[pair(0, k)][pair(0, k)] = p.hsp;
            for l in (k + 1)..=3usize {
                one[pair(k, k)][pair(l, l)] = p.gp2;
                one[pair(l, l)][pair(k, k)] = p.gp2;
                // <p_k p_l|p_k p_l> = hpp = ½(gpp − gp2).
                one[pair(l, k)][pair(l, k)] = 0.5 * (p.gpp - p.gp2);
            }
        }
        for &i in &orbs {
            for &j in &orbs {
                let mut sum = 0.0;
                for &k in &orbs {
                    for &l in &orbs {
                        let pij_kl =
                            one[pair(basis[i].1, basis[j].1)][pair(basis[k].1, basis[l].1)];
                        let pik_jl =
                            one[pair(basis[i].1, basis[k].1)][pair(basis[j].1, basis[l].1)];
                        sum += density[(k, l)] * pij_kl - 0.5 * density[(k, l)] * pik_jl;
                    }
                }
                g[(i, j)] += sum;
            }
        }
    }
    // Two-center contributions.
    for blk in blocks {
        let a_orbs: Vec<usize> = (0..n_basis).filter(|&i| basis[i].0 == blk.a).collect();
        let b_orbs: Vec<usize> = (0..n_basis).filter(|&i| basis[i].0 == blk.b).collect();
        // Coulomb.
        for &mu in &a_orbs {
            for &nu in &a_orbs {
                let pa = pair_index(basis[mu].1, basis[nu].1);
                let mut sum = 0.0;
                for &lam in &b_orbs {
                    for &sig in &b_orbs {
                        let pb = pair_index(basis[lam].1, basis[sig].1);
                        sum += density[(lam, sig)] * blk.ww[pa][pb];
                    }
                }
                g[(mu, nu)] += sum;
            }
        }
        for &lam in &b_orbs {
            for &sig in &b_orbs {
                let pb = pair_index(basis[lam].1, basis[sig].1);
                let mut sum = 0.0;
                for &mu in &a_orbs {
                    for &nu in &a_orbs {
                        let pa = pair_index(basis[mu].1, basis[nu].1);
                        sum += density[(mu, nu)] * blk.ww[pa][pb];
                    }
                }
                g[(lam, sig)] += sum;
            }
        }
        // Exchange: G_μν (μ∈A, ν∈B) -= ½ Σ P_λσ (μλ|νσ).
        for &mu in &a_orbs {
            for &nu in &b_orbs {
                let mut sum = 0.0;
                for &lam in &a_orbs {
                    for &sig in &b_orbs {
                        let pa = pair_index(basis[mu].1, basis[lam].1);
                        let pb = pair_index(basis[nu].1, basis[sig].1);
                        sum += density[(lam, sig)] * blk.ww[pa][pb];
                    }
                }
                let v = -0.5 * sum;
                g[(mu, nu)] += v;
                g[(nu, mu)] += v;
            }
        }
    }
    g
}

/// Run the faithful PM3 SCF (sp elements, Z <= 17).
pub fn solve(elements: &[u8], positions: &[[f64; 3]]) -> Result<FaithfulResult, String> {
    let n_atoms = elements.len();
    let params: Vec<&Pm3Params> = elements
        .iter()
        .map(|&z| {
            get_pm3_params(z).ok_or_else(|| format!("PM3 parameters not available for Z={z}"))
        })
        .collect::<Result<_, _>>()?;
    for &z in elements {
        if z > 17 {
            return Err(format!(
                "faithful PM3 currently covers Z <= 17 (got Z={z}); the generic Slater-overlap path is pending"
            ));
        }
    }
    let derived: Vec<Pm3Derived> = params.iter().map(|p| Pm3Derived::new(p)).collect();

    let basis = basis_layout(elements);
    let n_basis = basis.len();
    let n_electrons: usize = params.iter().map(|p| p.core_charge as usize).sum();
    if n_electrons % 2 != 0 {
        return Err(
            "faithful PM3 currently supports closed-shell (even electron) systems".to_string(),
        );
    }
    let n_occ = n_electrons / 2;

    let s_mat = overlap_matrix(elements, positions);

    // One-electron matrix + pair integral blocks.
    let mut h_core = DMatrix::zeros(n_basis, n_basis);
    for (idx, (atom, orb)) in basis.iter().enumerate() {
        let p = params[*atom];
        h_core[(idx, idx)] = if *orb == 0 { p.uss } else { p.upp };
    }
    let mut blocks: Vec<PairBlock> = Vec::new();
    let mut e_nuc = 0.0;
    for a in 0..n_atoms {
        for b in (a + 1)..n_atoms {
            let (ww, e1b, e2a, enuc) = two_center_block(
                params[a],
                params[b],
                &derived[a],
                &derived[b],
                &positions[a],
                &positions[b],
            );
            e_nuc += enuc;
            blocks.push(PairBlock { a, b, ww });
            // Attraction contributions to H.
            for (i1, &(atom_i, orb_i)) in basis.iter().enumerate() {
                if atom_i == a {
                    for (j1, &(atom_j, orb_j)) in basis.iter().enumerate() {
                        if atom_j == a {
                            let (r, c) = if orb_i >= orb_j {
                                (orb_i, orb_j)
                            } else {
                                (orb_j, orb_i)
                            };
                            h_core[(i1, j1)] += e1b[r][c];
                        }
                    }
                } else if atom_i == b {
                    for (j1, &(atom_j, orb_j)) in basis.iter().enumerate() {
                        if atom_j == b {
                            let (r, c) = if orb_i >= orb_j {
                                (orb_i, orb_j)
                            } else {
                                (orb_j, orb_i)
                            };
                            h_core[(i1, j1)] += e2a[r][c];
                        }
                    }
                }
            }
        }
    }
    // Resonance (off-diagonal) terms: h_μν = ½(β_μ + β_ν) S_μν.
    for i in 0..n_basis {
        for j in (i + 1)..n_basis {
            let (ai, oi) = basis[i];
            let (aj, oj) = basis[j];
            if ai == aj {
                continue;
            }
            let bi = if oi == 0 {
                params[ai].beta_s
            } else {
                params[ai].beta_p
            };
            let bj = if oj == 0 {
                params[aj].beta_s
            } else {
                params[aj].beta_p
            };
            let h = 0.5 * (bi + bj) * s_mat[(i, j)];
            h_core[(i, j)] = h;
            h_core[(j, i)] = h;
        }
    }

    // SCF with adaptive damping on the density (MOPAC converges with
    // Camp-King/Pulay; a damped fixed point is enough for the validated set).
    let mut density = DMatrix::zeros(n_basis, n_basis);
    let mut fock = h_core.clone();
    let mut converged = false;
    let mut iterations = 0;
    let mut prev_energy = f64::MAX;
    let mut e_elec = 0.0;
    let mut alpha: f64 = 0.4;
    let build_g = |density: &DMatrix<f64>| build_g_matrix(&params, &basis, &blocks, density);
    let density_from = |coeff: &DMatrix<f64>| -> DMatrix<f64> {
        let mut d = DMatrix::zeros(n_basis, n_basis);
        for i in 0..n_basis {
            for j in 0..n_basis {
                let mut sum = 0.0;
                for k in 0..n_occ {
                    sum += coeff[(i, k)] * coeff[(j, k)];
                }
                d[(i, j)] = 2.0 * sum;
            }
        }
        d
    };
    let trace_p_hf = |p: &DMatrix<f64>, f: &DMatrix<f64>| -> f64 {
        let mut e = 0.0;
        for i in 0..n_basis {
            for j in 0..n_basis {
                e += 0.5 * p[(i, j)] * (h_core[(i, j)] + f[(i, j)]);
            }
        }
        e
    };

    for iter in 0..200 {
        iterations = iter + 1;
        let (_, coeff) = solve_generalized(&fock, &s_mat);
        let p_new = density_from(&coeff);
        let g_new = build_g(&p_new);
        let f_new = &h_core + &g_new;
        e_elec = trace_p_hf(&p_new, &f_new);
        if std::env::var("PM3_TRACE").is_ok() {
            eprintln!(
                "iter {iter}: E_elec={e_elec:.6} P_diag={:?}",
                (0..n_basis).map(|i| p_new[(i, i)]).collect::<Vec<_>>()
            );
        }
        if (e_elec - prev_energy).abs() < 1e-9 && iter > 1 {
            converged = true;
            density = p_new;
            fock = f_new;
            break;
        }
        // Adaptive damping: shrink on energy increase, grow otherwise.
        if iter > 0 && e_elec > prev_energy {
            alpha = (alpha * 0.5).max(0.02);
        } else {
            alpha = (alpha * 1.05).min(0.9);
        }
        prev_energy = e_elec;
        density = &density * (1.0 - alpha) + &p_new * alpha;
        fock = &h_core + &build_g(&density);
    }

    let (orbital_energies, coeff) = solve_generalized(&fock, &s_mat);
    // Final density from the converged coefficients.
    let mut density = DMatrix::zeros(n_basis, n_basis);
    for i in 0..n_basis {
        for j in 0..n_basis {
            let mut sum = 0.0;
            for k in 0..n_occ {
                sum += coeff[(i, k)] * coeff[(j, k)];
            }
            density[(i, j)] = 2.0 * sum;
        }
    }
    e_elec = 0.0;
    for i in 0..n_basis {
        for j in 0..n_basis {
            e_elec += 0.5 * density[(i, j)] * (h_core[(i, j)] + fock[(i, j)]);
        }
    }

    let total_energy = e_elec + e_nuc;
    // MOPAC: ΔHf = (E_elec + E_nuc - Σ eisol)·23.0605 + Σ eheat.
    let atheat: f64 = params.iter().map(|p| p.heat_of_atomization).sum();
    let eisol_sum: f64 = params.iter().map(|p| p.eisol).sum();
    let heat_of_formation = (total_energy - eisol_sum) * KCAL_PER_EV + atheat;

    // Mulliken charges.
    let mut mulliken = vec![0.0f64; n_atoms];
    for i in 0..n_basis {
        let (atom, _) = basis[i];
        mulliken[atom] += density[(i, i)];
    }
    for (atom, &p) in params.iter().enumerate() {
        mulliken[atom] = p.core_charge - mulliken[atom];
    }

    let homo_energy = orbital_energies
        .get(n_occ.saturating_sub(1))
        .copied()
        .unwrap_or(0.0);
    let lumo_energy = orbital_energies.get(n_occ).copied().unwrap_or(0.0);
    Ok(FaithfulResult {
        orbital_energies,
        electronic_energy: e_elec,
        nuclear_repulsion: e_nuc,
        total_energy,
        heat_of_formation,
        homo_energy,
        lumo_energy,
        gap: lumo_energy - homo_energy,
        n_basis,
        n_electrons,
        mulliken_charges: mulliken,
        scf_iterations: iterations,
        converged,
        fock: fock.clone(),
        density: density.clone(),
    })
}

/// Debug helper for validation: overlap and one-electron matrices.
#[doc(hidden)]
pub fn debug_matrices(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<(DMatrix<f64>, DMatrix<f64>, f64, Vec<f64>), String> {
    let params: Vec<&Pm3Params> = elements
        .iter()
        .map(|&z| get_pm3_params(z).ok_or_else(|| format!("no params Z={z}")))
        .collect::<Result<_, _>>()?;
    let derived: Vec<Pm3Derived> = params.iter().map(|p| Pm3Derived::new(p)).collect();
    let basis = basis_layout(elements);
    let n_basis = basis.len();
    let s_mat = overlap_matrix(elements, positions);
    let mut h = DMatrix::zeros(n_basis, n_basis);
    for (idx, (atom, orb)) in basis.iter().enumerate() {
        h[(idx, idx)] = if *orb == 0 {
            params[*atom].uss
        } else {
            params[*atom].upp
        };
    }
    let mut e_nuc = 0.0;
    let mut attrs = vec![0.0f64; n_basis];
    for a in 0..elements.len() {
        for b in (a + 1)..elements.len() {
            let (_, e1b, e2a, enuc) = two_center_block(
                params[a],
                params[b],
                &derived[a],
                &derived[b],
                &positions[a],
                &positions[b],
            );
            e_nuc += enuc;
            for (i, &(ai, oi)) in basis.iter().enumerate() {
                if ai == a {
                    for (j, &(aj, oj)) in basis.iter().enumerate() {
                        if aj == a {
                            let (r, c) = if oi >= oj { (oi, oj) } else { (oj, oi) };
                            h[(i, j)] += e1b[r][c];
                            if i == j {
                                attrs[i] += e1b[r][c];
                            }
                        }
                    }
                } else if ai == b {
                    for (j, &(aj, oj)) in basis.iter().enumerate() {
                        if aj == b {
                            let (r, c) = if oi >= oj { (oi, oj) } else { (oj, oi) };
                            h[(i, j)] += e2a[r][c];
                            if i == j {
                                attrs[i] += e2a[r][c];
                            }
                        }
                    }
                }
            }
        }
    }
    for i in 0..n_basis {
        for j in (i + 1)..n_basis {
            let (ai, oi) = basis[i];
            let (aj, oj) = basis[j];
            if ai == aj {
                continue;
            }
            let bi = if oi == 0 {
                params[ai].beta_s
            } else {
                params[ai].beta_p
            };
            let bj = if oj == 0 {
                params[aj].beta_s
            } else {
                params[aj].beta_p
            };
            let v = 0.5 * (bi + bj) * s_mat[(i, j)];
            h[(i, j)] = v;
            h[(j, i)] = v;
        }
    }
    Ok((s_mat, h, e_nuc, attrs))
}

/// Debug: derived parameters (validation harness).
#[doc(hidden)]
pub fn debug_derived(p: &Pm3Params) -> Pm3Derived {
    Pm3Derived::new(p)
}

/// Debug: two-center block for one pair (validation harness).
#[doc(hidden)]
pub fn debug_block(
    za: u8,
    zb: u8,
    pos_a: &[f64; 3],
    pos_b: &[f64; 3],
) -> Result<([[f64; 10]; 10], [[f64; 4]; 4], [[f64; 4]; 4], f64), String> {
    let pa = get_pm3_params(za).ok_or("no params")?;
    let pb = get_pm3_params(zb).ok_or("no params")?;
    let da = Pm3Derived::new(pa);
    let db = Pm3Derived::new(pb);
    Ok(two_center_block(pa, pb, &da, &db, pos_a, pos_b))
}

/// Debug: G matrix for a given density (validation harness).
#[doc(hidden)]
pub fn debug_g(
    elements: &[u8],
    positions: &[[f64; 3]],
    density: &DMatrix<f64>,
) -> Result<DMatrix<f64>, String> {
    let params: Vec<&Pm3Params> = elements
        .iter()
        .map(|&z| get_pm3_params(z).ok_or_else(|| format!("no params Z={z}")))
        .collect::<Result<_, _>>()?;
    let derived: Vec<Pm3Derived> = params.iter().map(|p| Pm3Derived::new(p)).collect();
    let basis = basis_layout(elements);
    let mut blocks: Vec<PairBlock> = Vec::new();
    for a in 0..elements.len() {
        for b in (a + 1)..elements.len() {
            let (ww, _, _, _) = two_center_block(
                params[a],
                params[b],
                &derived[a],
                &derived[b],
                &positions[a],
                &positions[b],
            );
            blocks.push(PairBlock { a, b, ww });
        }
    }
    Ok(build_g_matrix(&params, &basis, &blocks, density))
}

#[cfg(test)]
mod mopac_validation_tests {
    use super::*;

    /// Core-core repulsion per X-H pair at 1.1 Å, vs MOPAC 23.2.5 (`PM3 1SCF DISP`).
    #[test]
    fn core_core_matches_mopac_for_x_h_pairs() {
        let cases: [(u8, f64); 8] = [
            (6, 40.00082),
            (7, 50.60331),
            (8, 63.89158),
            (9, 66.73093),
            (14, 29.38527),
            (15, 44.27043),
            (16, 55.20506),
            (17, 75.71863),
        ];
        for (z, mopac) in cases {
            let (_, _, _, enuc) = debug_block(z, 1, &[0.0, 0.0, 0.0], &[1.1, 0.0, 0.0]).unwrap();
            assert!(
                (enuc - mopac).abs() < 1e-4,
                "Z={z}: {enuc} vs MOPAC {mopac}"
            );
        }
    }

    /// Water one-electron matrix vs MOPAC 23.2.5 (`PM3 1SCF HCORE`).
    #[test]
    fn water_h_matrix_matches_mopac() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let (_, h, enuc, _) = debug_matrices(&elements, &positions).unwrap();
        let mopac_h = [
            [-108.417811, 0.0, -1.395684, 0.0, -7.341166, -7.341166],
            [0.0, -93.358540, 0.0, 0.0, -3.832985, 3.832985],
            [-1.395684, 0.0, -92.958650, 0.0, -2.967146, -2.967146],
            [0.0, 0.0, 0.0, -92.360702, 0.0, 0.0],
            [-7.341166, -3.832985, -2.967146, 0.0, -85.348085, -2.232226],
            [-7.341166, 3.832985, -2.967146, 0.0, -2.232226, -85.348085],
        ];
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (h[(i, j)] - mopac_h[i][j]).abs() < 1e-4,
                    "H[{i},{j}] = {} vs MOPAC {}",
                    h[(i, j)],
                    mopac_h[i][j]
                );
            }
        }
        // MOPAC: CORE-CORE REPULSION = 148.85535 EV.
        assert!((enuc - 148.85535).abs() < 1e-4, "E_nuc = {enuc}");
    }

    /// Water SCF state vs MOPAC: eigenvalues and Mulliken charges are exact.
    #[test]
    fn water_scf_matches_mopac() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let r = solve(&elements, &positions).unwrap();
        assert!(r.converged);
        let mopac_eigs = [-36.73327, -17.40836, -14.67440, -12.33356, 4.08042, 5.16984];
        for (e, m) in r.orbital_energies.iter().zip(mopac_eigs.iter()) {
            assert!((e - m).abs() < 1e-4, "eig {e} vs {m}");
        }
        let mopac_q = [-0.350797, 0.175399, 0.175399];
        for (q, m) in r.mulliken_charges.iter().zip(mopac_q.iter()) {
            assert!((q - m).abs() < 1e-4, "charge {q} vs {m}");
        }
        // MOPAC: FINAL HEAT OF FORMATION = -53.24348 (same geometry).
        assert!(
            (r.heat_of_formation - (-53.24348)).abs() < 1e-3,
            "HoF = {}",
            r.heat_of_formation
        );
    }

    /// H2 is exact end to end (HoF, eigenvalues, Fock).
    #[test]
    fn h2_matches_mopac_exactly() {
        let r = solve(&[1, 1], &[[0.0, 0.0, 0.0], [0.0, 0.0, 0.740]]).unwrap();
        assert!((r.heat_of_formation - (-12.71113)).abs() < 1e-4);
        assert!((r.orbital_energies[0] - (-15.87338)).abs() < 1e-4);
        assert!((r.fock[(0, 0)] - (-5.676214)).abs() < 1e-4);
    }
}
