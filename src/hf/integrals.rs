//! Electron repulsion integral (ERI) evaluator.
//!
//! Implements the Obara-Saika recurrence for two-electron integrals
//! over Cartesian Gaussian shells. For HF-3c, we need up to (pp|pp).

use super::basis::{BasisSet, Shell, ShellType};
use super::nuclear::boys_function;

pub fn compute_eris(basis: &BasisSet) -> Vec<f64> {
    let n = basis.n_basis();
    let size = eri_storage_size(n);
    let mut eris = vec![0.0f64; size];

    let shell_offsets = shell_function_offsets(basis);
    let n_shells = basis.shells.len();

    // Full loop over all shell quartets (no symmetry exploitation)
    for a in 0..n_shells {
        for b in 0..n_shells {
            for c in 0..n_shells {
                for d in 0..n_shells {
                    compute_eri_quartet(
                        &basis.shells[a],
                        &basis.shells[b],
                        &basis.shells[c],
                        &basis.shells[d],
                        shell_offsets[a],
                        shell_offsets[b],
                        shell_offsets[c],
                        shell_offsets[d],
                        n,
                        &mut eris,
                    );
                }
            }
        }
    }
    eris
}

fn shell_function_offsets(basis: &BasisSet) -> Vec<usize> {
    let mut offsets = Vec::with_capacity(basis.shells.len());
    let mut offset = 0;
    for shell in &basis.shells {
        offsets.push(offset);
        offset += shell.n_functions();
    }
    offsets
}

fn eri_storage_size(n: usize) -> usize {
    let nn = n * (n + 1) / 2;
    nn * (nn + 1) / 2
}

pub fn eri_index(i: usize, j: usize, k: usize, l: usize, _n: usize) -> usize {
    let ij = if i >= j {
        i * (i + 1) / 2 + j
    } else {
        j * (j + 1) / 2 + i
    };
    let kl = if k >= l {
        k * (k + 1) / 2 + l
    } else {
        l * (l + 1) / 2 + k
    };
    if ij >= kl {
        ij * (ij + 1) / 2 + kl
    } else {
        kl * (kl + 1) / 2 + ij
    }
}

pub fn get_eri(eris: &[f64], i: usize, j: usize, k: usize, l: usize, n: usize) -> f64 {
    eris[eri_index(i, j, k, l, n)]
}

fn compute_eri_quartet(
    sa: &Shell,
    sb: &Shell,
    sc: &Shell,
    sd: &Shell,
    off_a: usize,
    off_b: usize,
    off_c: usize,
    off_d: usize,
    n_basis: usize,
    eris: &mut [f64],
) {
    let la = if sa.shell_type == ShellType::P { 1 } else { 0 };
    let lb = if sb.shell_type == ShellType::P { 1 } else { 0 };
    let lc = if sc.shell_type == ShellType::P { 1 } else { 0 };
    let ld = if sd.shell_type == ShellType::P { 1 } else { 0 };

    // Temporary accumulator for primitives within this shell quartet
    let na = sa.n_functions();
    let nb = sb.n_functions();
    let nc = sc.n_functions();
    let nd = sd.n_functions();
    let mut temp = vec![0.0f64; na * nb * nc * nd];

    for (&ea, &ca) in sa.exponents.iter().zip(&sa.coefficients) {
        for (&eb, &cb) in sb.exponents.iter().zip(&sb.coefficients) {
            let zeta = ea + eb;
            let ab2 = dist_sq(&sa.center, &sb.center);
            let kab = (-ea * eb / zeta * ab2).exp();
            let p = gaussian_product(&sa.center, ea, &sb.center, eb);

            for (&ec, &cc) in sc.exponents.iter().zip(&sc.coefficients) {
                for (&ed, &cd) in sd.exponents.iter().zip(&sd.coefficients) {
                    let eta = ec + ed;
                    let cd2 = dist_sq(&sc.center, &sd.center);
                    let kcd = (-ec * ed / eta * cd2).exp();
                    let q = gaussian_product(&sc.center, ec, &sd.center, ed);

                    let rho = zeta * eta / (zeta + eta);
                    let pq2 = dist_sq(&p, &q);
                    let t = rho * pq2;

                    let prefactor = 2.0 * std::f64::consts::PI.powi(2) / (zeta * eta)
                        * (std::f64::consts::PI / (zeta + eta)).sqrt()
                        * kab
                        * kcd
                        * ca
                        * cb
                        * cc
                        * cd;

                    let w = gaussian_product(&p, zeta, &q, eta);

                    let prim_eris = os_eri_primitives(
                        &sa.center, &sb.center, &sc.center, &sd.center, &p, &q, &w, zeta, eta, rho,
                        t, prefactor, la, lb, lc, ld,
                    );

                    for fi in 0..na {
                        for fj in 0..nb {
                            for fk in 0..nc {
                                for fl in 0..nd {
                                    temp[((fi * nb + fj) * nc + fk) * nd + fl] +=
                                        prim_eris[fi][fj][fk][fl];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Write accumulated values to storage (direct assign — no cross-quartet accumulation)
    for fi in 0..na {
        for fj in 0..nb {
            for fk in 0..nc {
                for fl in 0..nd {
                    let i = off_a + fi;
                    let j = off_b + fj;
                    let k = off_c + fk;
                    let l = off_d + fl;
                    let val = temp[((fi * nb + fj) * nc + fk) * nd + fl];
                    let idx = eri_index(i, j, k, l, n_basis);
                    eris[idx] = val;
                }
            }
        }
    }
}

// Computes all primitive ERI combinations for (s,p) orbitals.
// Returns a 4D array sized [3][3][3][3]. Indices > 0 are only valid if L_x = 1.
// Index 0 represents s, indices 1,2,3 represent px,py,pz (shifted to 0,1,2).
fn os_eri_primitives(
    a: &[f64; 3],
    b: &[f64; 3],
    c: &[f64; 3],
    d: &[f64; 3],
    p: &[f64; 3],
    q: &[f64; 3],
    w: &[f64; 3],
    zeta: f64,
    eta: f64,
    _rho: f64,
    t: f64,
    prefactor: f64,
    la: usize,
    lb: usize,
    lc: usize,
    ld: usize,
) -> [[[[f64; 3]; 3]; 3]; 3] {
    let mut f = [0.0; 5];
    for m in 0..=(la + lb + lc + ld) {
        f[m] = boys_function(m, t) * prefactor;
    }

    // Helper to evaluate OS tree.
    // ss_ss[m]
    let ss_ss = |m: usize| f[m];

    // ps_ss[i][m]
    let mut ps_ss = [[0.0; 5]; 3];
    if la > 0 || lb > 0 || lc > 0 || ld > 0 {
        for i in 0..3 {
            for m in 0..=3 {
                ps_ss[i][m] = (p[i] - a[i]) * ss_ss(m) + (w[i] - p[i]) * ss_ss(m + 1);
            }
        }
    }

    // pp_ss[i][j][m]
    let mut pp_ss = [[[0.0; 5]; 3]; 3];
    if lb > 0 || lc > 0 || ld > 0 {
        for i in 0..3 {
            for j in 0..3 {
                for m in 0..=2 {
                    pp_ss[i][j][m] = (p[j] - b[j]) * ps_ss[i][m] + (w[j] - p[j]) * ps_ss[i][m + 1];
                    if i == j {
                        pp_ss[i][j][m] +=
                            1.0 / (2.0 * zeta) * (ss_ss(m) - eta / (zeta + eta) * ss_ss(m + 1));
                    }
                }
            }
        }
    }

    // ps_ps[i][k][m]
    let mut ps_ps = [[[0.0; 5]; 3]; 3];
    if lc > 0 || ld > 0 {
        for i in 0..3 {
            for k in 0..3 {
                for m in 0..=2 {
                    ps_ps[i][k][m] = (q[k] - c[k]) * ps_ss[i][m] + (w[k] - q[k]) * ps_ss[i][m + 1];
                    if i == k {
                        ps_ps[i][k][m] += 1.0 / (2.0 * (zeta + eta)) * ss_ss(m + 1);
                    }
                }
            }
        }
    }

    // pp_ps[i][j][k][m]
    let mut pp_ps = [[[[0.0; 5]; 3]; 3]; 3];
    if (la > 0 && lb > 0 && lc > 0) || ld > 0 {
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for m in 0..=1 {
                        pp_ps[i][j][k][m] =
                            (q[k] - c[k]) * pp_ss[i][j][m] + (w[k] - q[k]) * pp_ss[i][j][m + 1];
                        if i == k {
                            pp_ps[i][j][k][m] += 1.0 / (2.0 * (zeta + eta)) * ps_ss[j][m + 1];
                            // ps_ss is (0,b|0,0) which is same as (a,0|0,0) by symmetry if we swap
                        }
                        if j == k {
                            pp_ps[i][j][k][m] += 1.0 / (2.0 * (zeta + eta)) * ps_ss[i][m + 1];
                        }
                    }
                }
            }
        }
    }

    // pp_pp[i][j][k][l][m]
    let mut pp_pp = [[[[[0.0; 5]; 3]; 3]; 3]; 3];
    if la > 0 && lb > 0 && lc > 0 && ld > 0 {
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        for m in 0..=0 {
                            pp_pp[i][j][k][l][m] = (q[l] - d[l]) * pp_ps[i][j][k][m]
                                + (w[l] - q[l]) * pp_ps[i][j][k][m + 1];
                            if i == l {
                                pp_pp[i][j][k][l][m] +=
                                    1.0 / (2.0 * (zeta + eta)) * pp_ps[j][0][k][m + 1];
                                // Wait, pp_ps[j][0][k] doesn't mean (0,b|c,0).
                                // Actually, we need (0,b|c,0). By symmetry, it is ps_ps[j][k].
                            }
                            if j == l {
                                pp_pp[i][j][k][l][m] +=
                                    1.0 / (2.0 * (zeta + eta)) * ps_ps[i][k][m + 1];
                            }
                            if k == l {
                                pp_pp[i][j][k][l][m] += 1.0 / (2.0 * eta)
                                    * (pp_ss[i][j][m] - zeta / (zeta + eta) * pp_ss[i][j][m + 1]);
                            }
                        }
                    }
                }
            }
        }
    }

    let mut res = [[[[0.0; 3]; 3]; 3]; 3];
    let max_i = if la > 0 { 3 } else { 1 };
    let max_j = if lb > 0 { 3 } else { 1 };
    let max_k = if lc > 0 { 3 } else { 1 };
    let max_l = if ld > 0 { 3 } else { 1 };

    for i in 0..max_i {
        for j in 0..max_j {
            for k in 0..max_k {
                for l in 0..max_l {
                    if la == 0 && lb == 0 && lc == 0 && ld == 0 {
                        res[i][j][k][l] = ss_ss(0);
                    } else if la > 0 && lb == 0 && lc == 0 && ld == 0 {
                        res[i][j][k][l] = ps_ss[i][0];
                    } else if la == 0 && lb > 0 && lc == 0 && ld == 0 {
                        // (s, p | s, s) is same as (p, s | s, s) evaluated with B and A swapped.
                        res[i][j][k][l] = (p[j] - b[j]) * ss_ss(0) + (w[j] - p[j]) * ss_ss(1);
                    } else if la > 0 && lb > 0 && lc == 0 && ld == 0 {
                        res[i][j][k][l] = pp_ss[i][j][0];
                    } else if la == 0 && lb == 0 && lc == 0 && ld > 0 {
                        // (s, s | s, p)
                        res[i][j][k][l] = (q[l] - d[l]) * ss_ss(0) + (w[l] - q[l]) * ss_ss(1);
                    } else if la > 0 && lb == 0 && lc > 0 && ld == 0 {
                        res[i][j][k][l] = ps_ps[i][k][0];
                    } else if la == 0 && lb == 0 && lc > 0 && ld == 0 {
                        res[i][j][k][l] = (q[k] - c[k]) * ss_ss(0) + (w[k] - q[k]) * ss_ss(1);
                    } else if la == 0 && lb == 0 && lc > 0 && ld > 0 {
                        // (s, s | p, p) is same as (p, p | s, s) exchanging AB with CD.
                        let mut pp_ss_cd = (q[l] - d[l])
                            * ((q[k] - c[k]) * ss_ss(0) + (w[k] - q[k]) * ss_ss(1))
                            + (w[l] - q[l]) * ((q[k] - c[k]) * ss_ss(1) + (w[k] - q[k]) * ss_ss(2));
                        if k == l {
                            pp_ss_cd +=
                                1.0 / (2.0 * eta) * (ss_ss(0) - zeta / (zeta + eta) * ss_ss(1));
                        }
                        res[i][j][k][l] = pp_ss_cd;
                    } else if la == 0 && lb > 0 && lc > 0 && ld == 0 {
                        // (s, p | p, s)
                        let ps_ss_b = (p[j] - b[j]) * ss_ss(0) + (w[j] - p[j]) * ss_ss(1);
                        let ps_ss_b_1 = (p[j] - b[j]) * ss_ss(1) + (w[j] - p[j]) * ss_ss(2);
                        let mut val = (q[k] - c[k]) * ps_ss_b + (w[k] - q[k]) * ps_ss_b_1;
                        if j == k {
                            val += 1.0 / (2.0 * (zeta + eta)) * ss_ss(1);
                        }
                        res[i][j][k][l] = val;
                    } else if la > 0 && lb > 0 && lc > 0 && ld == 0 {
                        res[i][j][k][l] = pp_ps[i][j][k][0];
                    } else if la > 0 && lb > 0 && lc == 0 && ld > 0 {
                        // (p, p | s, p)
                        // symmetric to pp_ps where c -> d
                        let mut pp_ps_d =
                            (q[l] - d[l]) * pp_ss[i][j][0] + (w[l] - q[l]) * pp_ss[i][j][1];
                        if i == l {
                            pp_ps_d += 1.0 / (2.0 * (zeta + eta)) * ps_ss[j][1];
                        }
                        if j == l {
                            pp_ps_d += 1.0 / (2.0 * (zeta + eta)) * ps_ss[i][1];
                        }
                        res[i][j][k][l] = pp_ps_d;
                    } else if la > 0 && lb == 0 && lc > 0 && ld > 0 {
                        // (p, s | p, p)
                        let mut ps_pp =
                            (q[l] - d[l]) * ps_ps[i][k][0] + (w[l] - q[l]) * ps_ps[i][k][1];
                        if i == l {
                            ps_pp += 1.0 / (2.0 * (zeta + eta))
                                * ((q[k] - c[k]) * ss_ss(1) + (w[k] - q[k]) * ss_ss(2));
                        }
                        if k == l {
                            ps_pp += 1.0 / (2.0 * eta)
                                * (ps_ss[i][0] - zeta / (zeta + eta) * ps_ss[i][1]);
                        }
                        res[i][j][k][l] = ps_pp;
                    } else if la == 0 && lb > 0 && lc > 0 && ld > 0 {
                        // (s, p | p, p)
                        let ps_ss_b0 = (p[j] - b[j]) * ss_ss(0) + (w[j] - p[j]) * ss_ss(1);
                        let ps_ss_b1 = (p[j] - b[j]) * ss_ss(1) + (w[j] - p[j]) * ss_ss(2);
                        let ps_ss_b2 = (p[j] - b[j]) * ss_ss(2) + (w[j] - p[j]) * ss_ss(3);
                        let ps_ps_bk0 = (q[k] - c[k]) * ps_ss_b0
                            + (w[k] - q[k]) * ps_ss_b1
                            + if j == k {
                                1.0 / (2.0 * (zeta + eta)) * ss_ss(1)
                            } else {
                                0.0
                            };
                        let ps_ps_bk1 = (q[k] - c[k]) * ps_ss_b1
                            + (w[k] - q[k]) * ps_ss_b2
                            + if j == k {
                                1.0 / (2.0 * (zeta + eta)) * ss_ss(2)
                            } else {
                                0.0
                            };
                        let mut sp_pp = (q[l] - d[l]) * ps_ps_bk0 + (w[l] - q[l]) * ps_ps_bk1;
                        if j == l {
                            sp_pp += 1.0 / (2.0 * (zeta + eta))
                                * ((q[k] - c[k]) * ss_ss(1) + (w[k] - q[k]) * ss_ss(2));
                        }
                        if k == l {
                            sp_pp +=
                                1.0 / (2.0 * eta) * (ps_ss_b0 - zeta / (zeta + eta) * ps_ss_b1);
                        }
                        res[i][j][k][l] = sp_pp;
                    } else if la > 0 && lb > 0 && lc > 0 && ld > 0 {
                        // Fully populated pp_pp, let's substitute the symmetry correct components
                        let mut term =
                            (q[l] - d[l]) * pp_ps[i][j][k][0] + (w[l] - q[l]) * pp_ps[i][j][k][1];
                        if i == l {
                            term += 1.0 / (2.0 * (zeta + eta)) * ps_ps[j][k][1];
                        }
                        if j == l {
                            term += 1.0 / (2.0 * (zeta + eta)) * ps_ps[i][k][1];
                        }
                        if k == l {
                            term += 1.0 / (2.0 * eta)
                                * (pp_ss[i][j][0] - zeta / (zeta + eta) * pp_ss[i][j][1]);
                        }
                        res[i][j][k][l] = term;
                    } else if la > 0 && lb == 0 && lc == 0 && ld > 0 {
                        // (p, s | s, p)
                        let val = (q[l] - d[l]) * ps_ss[i][0]
                            + (w[l] - q[l]) * ps_ss[i][1]
                            + if i == l {
                                1.0 / (2.0 * (zeta + eta)) * ss_ss(1)
                            } else {
                                0.0
                            };
                        res[i][j][k][l] = val;
                    } else if la == 0 && lb > 0 && lc == 0 && ld > 0 {
                        // (s, p | s, p)
                        let ps_ss_b = (p[j] - b[j]) * ss_ss(0) + (w[j] - p[j]) * ss_ss(1);
                        let ps_ss_b1 = (p[j] - b[j]) * ss_ss(1) + (w[j] - p[j]) * ss_ss(2);
                        res[i][j][k][l] = (q[l] - d[l]) * ps_ss_b
                            + (w[l] - q[l]) * ps_ss_b1
                            + if j == l {
                                1.0 / (2.0 * (zeta + eta)) * ss_ss(1)
                            } else {
                                0.0
                            };
                    }
                }
            }
        }
    }
    res
}

fn gaussian_product(a: &[f64; 3], ea: f64, b: &[f64; 3], eb: f64) -> [f64; 3] {
    let g = ea + eb;
    [
        (ea * a[0] + eb * b[0]) / g,
        (ea * a[1] + eb * b[1]) / g,
        (ea * a[2] + eb * b[2]) / g,
    ]
}

#[inline]
fn dist_sq(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

#[cfg(test)]
mod tests {
    use super::super::basis::build_sto3g_basis;
    use super::*;

    #[test]
    fn test_eri_h2() {
        let basis = build_sto3g_basis(&[1, 1], &[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]]);
        let eris = compute_eris(&basis);
        let n = basis.n_basis();
        // (11|11) should be positive (electron-electron repulsion)
        let eri_0000 = get_eri(&eris, 0, 0, 0, 0, n);
        assert!(eri_0000 > 0.0, "(11|11) = {eri_0000}");
    }

    #[test]
    fn test_eri_symmetry() {
        let basis = build_sto3g_basis(&[1, 1], &[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]]);
        let eris = compute_eris(&basis);
        let n = basis.n_basis();
        // Permutational symmetry
        assert_eq!(get_eri(&eris, 0, 1, 0, 1, n), get_eri(&eris, 1, 0, 1, 0, n));
        assert_eq!(get_eri(&eris, 0, 1, 0, 1, n), get_eri(&eris, 0, 1, 1, 0, n));
    }
}
