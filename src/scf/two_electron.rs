//! Two-electron repulsion integrals (μν|λσ).
//!
//! (μν|λσ) = ∫∫ χ_μ(r₁)χ_ν(r₁) (1/r₁₂) χ_λ(r₂)χ_σ(r₂) d³r₁ d³r₂
//!
//! These are the most expensive integrals in quantum chemistry, scaling
//! as O(N⁴) for N basis functions.
//!
//! # Algorithm
//!
//! Uses the Obara-Saika scheme for electron repulsion integrals (ERIs).

use std::collections::HashMap;
use std::f64::consts::PI;

use super::basis::{BasisFunction, BasisSet};
use super::gaussian_integrals::{boys_function, distance_squared, gaussian_product_center};

/// Store two-electron integrals in a compact format.
///
/// For N basis functions, there are N⁴ integrals but only
/// N(N+1)/2 × (N(N+1)/2 + 1)/2 unique ones due to symmetry:
///   (μν|λσ) = (νμ|λσ) = (μν|σλ) = (λσ|μν) = ...
#[derive(Debug, Clone)]
pub struct TwoElectronIntegrals {
    /// Flat storage of integrals indexed by compound index.
    data: Vec<f64>,
    /// Number of basis functions.
    n_basis: usize,
}

impl TwoElectronIntegrals {
    /// Compute all unique two-electron integrals for the basis set.
    pub fn compute(basis: &BasisSet) -> Self {
        let n = basis.n_basis;
        let n2 = n * n;
        let mut data = vec![0.0f64; n2 * n2];

        for i in 0..n {
            for j in 0..=i {
                let ij = i * n + j;
                for k in 0..n {
                    for l in 0..=k {
                        let kl = k * n + l;
                        if ij < kl {
                            continue;
                        }

                        let eri = contracted_eri(
                            &basis.functions[i],
                            &basis.functions[j],
                            &basis.functions[k],
                            &basis.functions[l],
                        );

                        // Store all 8-fold symmetry permutations
                        data[i * n * n2 + j * n2 + k * n + l] = eri;
                        data[j * n * n2 + i * n2 + k * n + l] = eri;
                        data[i * n * n2 + j * n2 + l * n + k] = eri;
                        data[j * n * n2 + i * n2 + l * n + k] = eri;
                        data[k * n * n2 + l * n2 + i * n + j] = eri;
                        data[l * n * n2 + k * n2 + i * n + j] = eri;
                        data[k * n * n2 + l * n2 + j * n + i] = eri;
                        data[l * n * n2 + k * n2 + j * n + i] = eri;
                    }
                }
            }
        }

        Self { data, n_basis: n }
    }

    /// Get integral (μν|λσ).
    #[inline]
    pub fn get(&self, mu: usize, nu: usize, lam: usize, sig: usize) -> f64 {
        let n = self.n_basis;
        self.data[mu * n * n * n + nu * n * n + lam * n + sig]
    }

    /// Number of basis functions.
    pub fn n_basis(&self) -> usize {
        self.n_basis
    }

    /// Construct from pre-computed raw data (e.g. from GPU dispatch).
    pub fn from_raw(data: Vec<f64>, n_basis: usize) -> Self {
        debug_assert_eq!(data.len(), n_basis * n_basis * n_basis * n_basis);
        Self { data, n_basis }
    }

    /// Compute two-electron integrals using rayon parallelism.
    ///
    /// Parallelizes the outer `i` loop. Each thread writes to a
    /// disjoint sub-region via a Mutex-protected accumulation step.
    #[cfg(feature = "parallel")]
    pub fn compute_parallel(basis: &BasisSet) -> Self {
        use rayon::prelude::*;
        use std::sync::Mutex;

        let n = basis.n_basis;
        let n2 = n * n;
        let data = Mutex::new(vec![0.0f64; n2 * n2]);

        (0..n).into_par_iter().for_each(|i| {
            let mut local: Vec<(usize, f64)> = Vec::new();
            for j in 0..=i {
                let ij = i * n + j;
                for k in 0..n {
                    for l in 0..=k {
                        let kl = k * n + l;
                        if ij < kl {
                            continue;
                        }
                        let eri = contracted_eri(
                            &basis.functions[i],
                            &basis.functions[j],
                            &basis.functions[k],
                            &basis.functions[l],
                        );
                        local.push((i * n * n2 + j * n2 + k * n + l, eri));
                        local.push((j * n * n2 + i * n2 + k * n + l, eri));
                        local.push((i * n * n2 + j * n2 + l * n + k, eri));
                        local.push((j * n * n2 + i * n2 + l * n + k, eri));
                        local.push((k * n * n2 + l * n2 + i * n + j, eri));
                        local.push((l * n * n2 + k * n2 + i * n + j, eri));
                        local.push((k * n * n2 + l * n2 + j * n + i, eri));
                        local.push((l * n * n2 + k * n2 + j * n + i, eri));
                    }
                }
            }
            let mut d = data.lock().unwrap();
            for (idx, val) in local {
                d[idx] = val;
            }
        });

        Self {
            data: data.into_inner().unwrap(),
            n_basis: n,
        }
    }
}

/// Contracted ERI between four basis functions.
fn contracted_eri(
    bf_a: &BasisFunction,
    bf_b: &BasisFunction,
    bf_c: &BasisFunction,
    bf_d: &BasisFunction,
) -> f64 {
    let mut eri = 0.0;

    for pa in &bf_a.primitives {
        let na = BasisFunction::normalization(
            pa.alpha,
            bf_a.angular[0],
            bf_a.angular[1],
            bf_a.angular[2],
        );
        for pb in &bf_b.primitives {
            let nb = BasisFunction::normalization(
                pb.alpha,
                bf_b.angular[0],
                bf_b.angular[1],
                bf_b.angular[2],
            );
            for pc in &bf_c.primitives {
                let nc = BasisFunction::normalization(
                    pc.alpha,
                    bf_c.angular[0],
                    bf_c.angular[1],
                    bf_c.angular[2],
                );
                for pd in &bf_d.primitives {
                    let nd = BasisFunction::normalization(
                        pd.alpha,
                        bf_d.angular[0],
                        bf_d.angular[1],
                        bf_d.angular[2],
                    );

                    let prim_eri = eri_primitive(
                        pa.alpha,
                        &bf_a.center,
                        bf_a.angular,
                        pb.alpha,
                        &bf_b.center,
                        bf_b.angular,
                        pc.alpha,
                        &bf_c.center,
                        bf_c.angular,
                        pd.alpha,
                        &bf_d.center,
                        bf_d.angular,
                    );

                    eri += na
                        * pa.coefficient
                        * nb
                        * pb.coefficient
                        * nc
                        * pc.coefficient
                        * nd
                        * pd.coefficient
                        * prim_eri;
                }
            }
        }
    }

    eri
}

/// ERI between four primitive Gaussians.
fn eri_primitive(
    alpha: f64,
    center_a: &[f64; 3],
    la: [u32; 3],
    beta: f64,
    center_b: &[f64; 3],
    lb: [u32; 3],
    gamma: f64,
    center_c: &[f64; 3],
    lc: [u32; 3],
    delta: f64,
    center_d: &[f64; 3],
    ld: [u32; 3],
) -> f64 {
    let p = alpha + beta;
    let q = gamma + delta;
    let alpha_pq = p * q / (p + q);

    let mu_ab = alpha * beta / p;
    let mu_cd = gamma * delta / q;

    let ab2 = distance_squared(center_a, center_b);
    let cd2 = distance_squared(center_c, center_d);

    let px = gaussian_product_center(alpha, center_a[0], beta, center_b[0]);
    let py = gaussian_product_center(alpha, center_a[1], beta, center_b[1]);
    let pz = gaussian_product_center(alpha, center_a[2], beta, center_b[2]);

    let qx = gaussian_product_center(gamma, center_c[0], delta, center_d[0]);
    let qy = gaussian_product_center(gamma, center_c[1], delta, center_d[1]);
    let qz = gaussian_product_center(gamma, center_c[2], delta, center_d[2]);

    let pq2 = (px - qx).powi(2) + (py - qy).powi(2) + (pz - qz).powi(2);

    let prefactor = 2.0 * PI.powf(2.5) / (p * q * (p + q).sqrt());
    let k_ab = (-mu_ab * ab2).exp();
    let k_cd = (-mu_cd * cd2).exp();
    let pref = prefactor * k_ab * k_cd;
    let t = alpha_pq * pq2;

    let p_vec = [px, py, pz];
    let q_vec = [qx, qy, qz];
    let w = [
        gaussian_product_center(p, px, q, qx),
        gaussian_product_center(p, py, q, qy),
        gaussian_product_center(p, pz, q, qz),
    ];

    let ctx = OsCtx {
        a: *center_a,
        b: *center_b,
        c: *center_c,
        d: *center_d,
        p: p_vec,
        q: q_vec,
        w,
        zeta: p,
        eta: q,
        rho: alpha_pq,
        pref,
        t,
    };
    let mut cache = HashMap::new();
    os_cartesian(&ctx, la, lb, lc, ld, 0, &mut cache)
}

/// Shared geometry for the Obara–Saika Cartesian ERI recurrence.
struct OsCtx {
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
    p: [f64; 3],
    q: [f64; 3],
    w: [f64; 3],
    zeta: f64,
    eta: f64,
    rho: f64,
    pref: f64,
    t: f64,
}

#[inline]
fn dec(l: [u32; 3], axis: usize) -> [u32; 3] {
    let mut out = l;
    out[axis] -= 1;
    out
}

#[inline]
fn pack_os_key(la: [u32; 3], lb: [u32; 3], lc: [u32; 3], ld: [u32; 3], m: u32) -> u64 {
    let pack = |l: [u32; 3], shift: u32| {
        (l[0] as u64) << shift | (l[1] as u64) << (shift + 4) | (l[2] as u64) << (shift + 8)
    };
    pack(la, 0) | pack(lb, 12) | pack(lc, 24) | pack(ld, 36) | (m as u64) << 48
}

/// Cartesian Obara–Saika ERI. Valid for s/p/d Cartesian Gaussians.
fn os_cartesian(
    ctx: &OsCtx,
    la: [u32; 3],
    lb: [u32; 3],
    lc: [u32; 3],
    ld: [u32; 3],
    m: u32,
    cache: &mut HashMap<u64, f64>,
) -> f64 {
    let key = pack_os_key(la, lb, lc, ld, m);
    if let Some(&hit) = cache.get(&key) {
        return hit;
    }
    let val = os_cartesian_uncached(ctx, la, lb, lc, ld, m, cache);
    cache.insert(key, val);
    val
}

fn os_cartesian_uncached(
    ctx: &OsCtx,
    la: [u32; 3],
    lb: [u32; 3],
    lc: [u32; 3],
    ld: [u32; 3],
    m: u32,
    cache: &mut HashMap<u64, f64>,
) -> f64 {
    let lsum = la[0]
        + la[1]
        + la[2]
        + lb[0]
        + lb[1]
        + lb[2]
        + lc[0]
        + lc[1]
        + lc[2]
        + ld[0]
        + ld[1]
        + ld[2];
    if lsum == 0 {
        return ctx.pref * boys_function(m, ctx.t);
    }

    let zeta_eta = ctx.zeta + ctx.eta;
    if let Some(i) = (0..3).find(|&i| la[i] > 0) {
        let lam = dec(la, i);
        let mut val = (ctx.p[i] - ctx.a[i]) * os_cartesian(ctx, lam, lb, lc, ld, m, cache)
            + (ctx.w[i] - ctx.p[i]) * os_cartesian(ctx, lam, lb, lc, ld, m + 1, cache);
        if lam[i] > 0 {
            val += lam[i] as f64 / (2.0 * ctx.zeta)
                * (os_cartesian(ctx, dec(lam, i), lb, lc, ld, m, cache)
                    - ctx.rho / ctx.zeta
                        * os_cartesian(ctx, dec(lam, i), lb, lc, ld, m + 1, cache));
        }
        if lb[i] > 0 {
            val += lb[i] as f64 / (2.0 * ctx.zeta)
                * (os_cartesian(ctx, lam, dec(lb, i), lc, ld, m, cache)
                    - ctx.rho / ctx.zeta
                        * os_cartesian(ctx, lam, dec(lb, i), lc, ld, m + 1, cache));
        }
        if lc[i] > 0 {
            val += lc[i] as f64 / (2.0 * zeta_eta)
                * os_cartesian(ctx, lam, lb, dec(lc, i), ld, m + 1, cache);
        }
        if ld[i] > 0 {
            val += ld[i] as f64 / (2.0 * zeta_eta)
                * os_cartesian(ctx, lam, lb, lc, dec(ld, i), m + 1, cache);
        }
        return val;
    }
    if let Some(i) = (0..3).find(|&i| lb[i] > 0) {
        let lbm = dec(lb, i);
        let mut val = (ctx.p[i] - ctx.b[i]) * os_cartesian(ctx, la, lbm, lc, ld, m, cache)
            + (ctx.w[i] - ctx.p[i]) * os_cartesian(ctx, la, lbm, lc, ld, m + 1, cache);
        if la[i] > 0 {
            val += la[i] as f64 / (2.0 * ctx.zeta)
                * (os_cartesian(ctx, dec(la, i), lbm, lc, ld, m, cache)
                    - ctx.rho / ctx.zeta
                        * os_cartesian(ctx, dec(la, i), lbm, lc, ld, m + 1, cache));
        }
        if lbm[i] > 0 {
            val += lbm[i] as f64 / (2.0 * ctx.zeta)
                * (os_cartesian(ctx, la, dec(lbm, i), lc, ld, m, cache)
                    - ctx.rho / ctx.zeta
                        * os_cartesian(ctx, la, dec(lbm, i), lc, ld, m + 1, cache));
        }
        if lc[i] > 0 {
            val += lc[i] as f64 / (2.0 * zeta_eta)
                * os_cartesian(ctx, la, lbm, dec(lc, i), ld, m + 1, cache);
        }
        if ld[i] > 0 {
            val += ld[i] as f64 / (2.0 * zeta_eta)
                * os_cartesian(ctx, la, lbm, lc, dec(ld, i), m + 1, cache);
        }
        return val;
    }
    if let Some(i) = (0..3).find(|&i| lc[i] > 0) {
        let lcm = dec(lc, i);
        let mut val = (ctx.q[i] - ctx.c[i]) * os_cartesian(ctx, la, lb, lcm, ld, m, cache)
            + (ctx.w[i] - ctx.q[i]) * os_cartesian(ctx, la, lb, lcm, ld, m + 1, cache);
        if lcm[i] > 0 {
            val += lcm[i] as f64 / (2.0 * ctx.eta)
                * (os_cartesian(ctx, la, lb, dec(lcm, i), ld, m, cache)
                    - ctx.rho / ctx.eta * os_cartesian(ctx, la, lb, dec(lcm, i), ld, m + 1, cache));
        }
        if ld[i] > 0 {
            val += ld[i] as f64 / (2.0 * ctx.eta)
                * (os_cartesian(ctx, la, lb, lcm, dec(ld, i), m, cache)
                    - ctx.rho / ctx.eta * os_cartesian(ctx, la, lb, lcm, dec(ld, i), m + 1, cache));
        }
        if la[i] > 0 {
            val += la[i] as f64 / (2.0 * zeta_eta)
                * os_cartesian(ctx, dec(la, i), lb, lcm, ld, m + 1, cache);
        }
        if lb[i] > 0 {
            val += lb[i] as f64 / (2.0 * zeta_eta)
                * os_cartesian(ctx, la, dec(lb, i), lcm, ld, m + 1, cache);
        }
        return val;
    }
    let i = (0..3)
        .find(|&i| ld[i] > 0)
        .expect("ld has angular momentum");
    let ldm = dec(ld, i);
    let mut val = (ctx.q[i] - ctx.d[i]) * os_cartesian(ctx, la, lb, lc, ldm, m, cache)
        + (ctx.w[i] - ctx.q[i]) * os_cartesian(ctx, la, lb, lc, ldm, m + 1, cache);
    if lc[i] > 0 {
        val += lc[i] as f64 / (2.0 * ctx.eta)
            * (os_cartesian(ctx, la, lb, dec(lc, i), ldm, m, cache)
                - ctx.rho / ctx.eta * os_cartesian(ctx, la, lb, dec(lc, i), ldm, m + 1, cache));
    }
    if ldm[i] > 0 {
        val += ldm[i] as f64 / (2.0 * ctx.eta)
            * (os_cartesian(ctx, la, lb, lc, dec(ldm, i), m, cache)
                - ctx.rho / ctx.eta * os_cartesian(ctx, la, lb, lc, dec(ldm, i), m + 1, cache));
    }
    if la[i] > 0 {
        val += la[i] as f64 / (2.0 * zeta_eta)
            * os_cartesian(ctx, dec(la, i), lb, lc, ldm, m + 1, cache);
    }
    if lb[i] > 0 {
        val += lb[i] as f64 / (2.0 * zeta_eta)
            * os_cartesian(ctx, la, dec(lb, i), lc, ldm, m + 1, cache);
    }
    val
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eri_h2_computed() {
        let basis = BasisSet::sto3g(&[1, 1], &[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]);
        let eris = TwoElectronIntegrals::compute(&basis);

        assert!(eris.get(0, 0, 0, 0) > 0.0);
        assert!((eris.get(0, 1, 0, 0) - eris.get(1, 0, 0, 0)).abs() < 1e-14);
    }

    #[test]
    fn test_eri_symmetry() {
        let basis = BasisSet::sto3g(&[1], &[[0.0, 0.0, 0.0]]);
        let eris = TwoElectronIntegrals::compute(&basis);
        assert!(eris.get(0, 0, 0, 0) > 0.0);
    }

    #[test]
    fn same_center_pxpx_ss_matches_os_closed_form() {
        let origin = [0.0, 0.0, 0.0];
        let alpha = 0.5;
        let beta = 0.6;
        let gamma = 0.7;
        let delta = 0.8;
        let s = [0u32, 0, 0];
        let px = [1u32, 0, 0];
        let got = eri_primitive(
            alpha, &origin, px, beta, &origin, px, gamma, &origin, s, delta, &origin, s,
        );

        let zeta = alpha + beta;
        let eta = gamma + delta;
        let rho = zeta * eta / (zeta + eta);
        let pref = 2.0 * PI.powf(2.5) / (zeta * eta * (zeta + eta).sqrt());
        let f0 = boys_function(0, 0.0);
        let f1 = boys_function(1, 0.0);
        let expected = pref / (2.0 * zeta) * (f0 - rho / zeta * f1);
        assert!(
            (got - expected).abs() < 1e-12,
            "same-center (px px|ss): got {got}, expected {expected}"
        );
        let ssss = eri_primitive(
            alpha, &origin, s, beta, &origin, s, gamma, &origin, s, delta, &origin, s,
        );
        assert!(
            (got - ssss).abs() > 1e-8,
            "(px px|ss) must differ from (ss|ss)"
        );
    }

    #[test]
    fn primitive_ps_ss_and_pp_pp_are_finite_and_asymmetric() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.2, 0.0, 0.0];
        let s = [0u32, 0, 0];
        let px = [1u32, 0, 0];
        let py = [0u32, 1, 0];
        let psss = eri_primitive(0.4, &a, px, 0.5, &b, s, 0.6, &a, s, 0.7, &b, s);
        let ssss = eri_primitive(0.4, &a, s, 0.5, &b, s, 0.6, &a, s, 0.7, &b, s);
        let pppp = eri_primitive(0.4, &a, px, 0.5, &b, py, 0.6, &a, px, 0.7, &b, py);
        assert!(psss.is_finite() && pppp.is_finite());
        assert!((psss - ssss).abs() > 1e-10);
        let ppss = eri_primitive(0.4, &a, px, 0.5, &b, px, 0.6, &a, s, 0.7, &b, s);
        let psps = eri_primitive(0.4, &a, px, 0.5, &b, s, 0.6, &a, px, 0.7, &b, s);
        assert!(ppss.is_finite() && psps.is_finite());
        assert!((ppss - ssss).abs() > 1e-10);
    }

    #[test]
    fn same_center_dxx_s_ss_matches_os_closed_form() {
        let origin = [0.0, 0.0, 0.0];
        let alpha = 0.5;
        let beta = 0.6;
        let gamma = 0.7;
        let delta = 0.8;
        let s = [0u32, 0, 0];
        let dxx = [2u32, 0, 0];
        let dxy = [1u32, 1, 0];
        let got = eri_primitive(
            alpha, &origin, dxx, beta, &origin, s, gamma, &origin, s, delta, &origin, s,
        );
        let zeta = alpha + beta;
        let eta = gamma + delta;
        let rho = zeta * eta / (zeta + eta);
        let pref = 2.0 * PI.powf(2.5) / (zeta * eta * (zeta + eta).sqrt());
        let expected =
            pref / (2.0 * zeta) * (boys_function(0, 0.0) - rho / zeta * boys_function(1, 0.0));
        assert!(
            (got - expected).abs() < 1e-12,
            "same-center (dxx s|ss): got {got}, expected {expected}"
        );
        let odd = eri_primitive(
            alpha, &origin, dxy, beta, &origin, s, gamma, &origin, s, delta, &origin, s,
        );
        assert!(
            odd.abs() < 1e-12,
            "(dxy s|ss) same-center must vanish, got {odd}"
        );
        let dddd = eri_primitive(
            alpha, &origin, dxx, beta, &origin, dxx, gamma, &origin, dxx, delta, &origin, dxx,
        );
        assert!(dddd.is_finite() && dddd > 0.0);
        let dxx_dxx_ss = eri_primitive(
            alpha, &origin, dxx, beta, &origin, dxx, gamma, &origin, s, delta, &origin, s,
        );
        let ss_dxx_dxx = eri_primitive(
            gamma, &origin, s, delta, &origin, s, alpha, &origin, dxx, beta, &origin, dxx,
        );
        assert!((dxx_dxx_ss - ss_dxx_dxx).abs() < 1e-12);
        assert!((dxx_dxx_ss - got).abs() > 1e-8);
    }

    #[test]
    fn cartesian_d_shell_contracted_eris() {
        use super::super::basis::{ContractedShell, GaussianPrimitive};
        let origin = [0.0, 0.0, 0.0];
        let prim = |alpha, coefficient| GaussianPrimitive { alpha, coefficient };
        let s = ContractedShell {
            atom_index: 0,
            center: origin,
            l: 0,
            primitives: vec![prim(0.5, 1.0)],
        };
        let d = ContractedShell {
            atom_index: 0,
            center: origin,
            l: 2,
            primitives: vec![prim(0.4, 1.0)],
        };
        let basis = BasisSet::from_shells(vec![s, d]);
        assert_eq!(basis.n_basis, 7);
        let eris = TwoElectronIntegrals::compute(&basis);
        let ssss = eris.get(0, 0, 0, 0);
        let dxx = 1; // xx after the s function
        let dxx_dxx_ss = eris.get(dxx, dxx, 0, 0);
        assert!(ssss > 0.0 && dxx_dxx_ss.is_finite());
        assert!((dxx_dxx_ss - ssss).abs() > 1e-6);
        assert!((eris.get(dxx, dxx, 0, 0) - eris.get(0, 0, dxx, dxx)).abs() < 1e-12);
        let dxy = 2;
        assert!(eris.get(dxy, 0, 0, 0).abs() < 1e-10);
    }

    #[test]
    fn water_sto3g_eri_symmetry_and_p_shells() {
        // O at origin; H in Bohr (0.757 Å, 0.586 Å).
        let ang = 1.889_726_124_6;
        let basis = BasisSet::sto3g(
            &[8, 1, 1],
            &[
                [0.0, 0.0, 0.0],
                [0.757 * ang, 0.586 * ang, 0.0],
                [-0.757 * ang, 0.586 * ang, 0.0],
            ],
        );
        assert_eq!(basis.n_basis, 7);
        let eris = TwoElectronIntegrals::compute(&basis);
        let n = eris.n_basis();
        let mut rms = 0.0;
        let mut count = 0u64;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    for l in 0..n {
                        let v1 = eris.get(i, j, k, l);
                        assert!(v1.is_finite());
                        assert!((v1 - eris.get(j, i, k, l)).abs() < 1e-12);
                        assert!((v1 - eris.get(i, j, l, k)).abs() < 1e-12);
                        assert!((v1 - eris.get(k, l, i, j)).abs() < 1e-12);
                        rms += v1 * v1;
                        count += 1;
                    }
                }
            }
        }
        rms = (rms / count as f64).sqrt();
        // STO-3G water ERIs are O(1) Ha; a vanishing RMS would mean the p stub collapsed.
        assert!(rms > 0.1 && rms < 20.0, "H2O ERI RMS = {rms}");
        // Oxygen 2px is basis index 2 in STO-3G (1s, 2s, 2px, 2py, 2pz, H1s, H1s).
        let pxpx_ss = eris.get(2, 2, 0, 0);
        let ssss = eris.get(0, 0, 0, 0);
        assert!(
            (pxpx_ss - ssss).abs() > 1e-4,
            "(2px 2px|1s 1s) must differ from (1s 1s|1s 1s)"
        );
    }

    #[test]
    fn test_eri_sequential_vs_sequential_consistency() {
        let basis = BasisSet::sto3g(&[1, 1], &[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]);
        let eris = TwoElectronIntegrals::compute(&basis);

        // Test 8-fold symmetry
        let n = eris.n_basis();
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    for l in 0..n {
                        let v1 = eris.get(i, j, k, l);
                        let v2 = eris.get(j, i, k, l);
                        let v3 = eris.get(i, j, l, k);
                        let v4 = eris.get(k, l, i, j);
                        assert!((v1 - v2).abs() < 1e-14, "Symmetry (ij) failed");
                        assert!((v1 - v3).abs() < 1e-14, "Symmetry (kl) failed");
                        assert!((v1 - v4).abs() < 1e-14, "Symmetry (ij|kl)↔(kl|ij) failed");
                    }
                }
            }
        }
    }
}
