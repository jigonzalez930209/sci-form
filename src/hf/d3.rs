//! Grimme DFT-D3 dispersion correction with Becke-Johnson damping.
//!
//! Adds long-range van der Waals interactions missing from HF:
//! $$E_{disp} = -\sum_{n=6,8} s_n \sum_{A>B} \frac{C_n^{AB}}{R_{AB}^n + f_{damp}^n}$$
//!
//! Reference: Grimme, S. et al. J. Chem. Phys. 132 (2010): 154104.

use serde::{Deserialize, Serialize};

/// D3 correction result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct D3Result {
    /// Total D3 dispersion energy (Hartree).
    pub energy: f64,
    /// C6 contribution.
    pub e6: f64,
    /// C8 contribution.
    pub e8: f64,
}

/// D3 parameters for HF-3c (BJ damping).
const S6: f64 = 1.0;
const S8: f64 = 0.8777;
const A1: f64 = 0.4171;
const A2: f64 = 2.9149;

/// D3 reference data (generated tables, simple-dftd3 reference).
#[path = "d3_data.rs"]
mod d3_data;
use d3_data::{D3_C6, D3_COVRAD_BOHR, D3_NREFS, D3_R4R2, D3_REFCN, D3_VDWRAD_BOHR};

/// D3 model dimensions (reference states per element at most).
const D3_MAX_REF: usize = 7;
/// Elements covered by the D3 tables.
const D3_MAX_ELEM: usize = 103;
/// Gaussian interpolation width for CN weighting (simple-dftd3 `wf`).
const D3_WF: f64 = 4.0;
/// Steepness of the logistic CN counting function (simple-dftd3 `kcn`).
const D3_KCN: f64 = 16.0;

/// Triangular pair index (0-based) for 1-based atomic numbers with `za <= zb`.
fn d3_pair_index(za: u8, zb: u8) -> Option<usize> {
    let (a, b) = (za as usize, zb as usize);
    if a == 0 || b == 0 || a > D3_MAX_ELEM || b > D3_MAX_ELEM {
        return None;
    }
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    Some(lo + hi * (hi - 1) / 2 - 1)
}

/// Reference C6 for reference states (`iref`, `jref`) of elements
/// (`za`, `zb`), mirroring simple-dftd3's canonical index ordering.
fn d3_c6_ref(iref: usize, jref: usize, za: u8, zb: u8) -> f64 {
    let ic = match d3_pair_index(za.min(zb), za.max(zb)) {
        Some(ic) => ic,
        None => return 0.0,
    };
    // Fortran: ati > atj → ref(iref,jref,ic) else ref(jref,iref,ic).
    let (r1, r2) = if za > zb { (iref, jref) } else { (jref, iref) };
    if r1 >= D3_MAX_REF || r2 >= D3_MAX_REF {
        return 0.0;
    }
    D3_C6[r1 + D3_MAX_REF * (r2 + D3_MAX_REF * ic)]
}

/// Fractional coordination numbers (logistic counting, D3 covalent radii).
fn d3_cn(elements: &[u8], positions_bohr: &[[f64; 3]]) -> Vec<f64> {
    let n = elements.len();
    let mut cn = vec![0.0; n];
    for i in 0..n {
        let zi = elements[i] as usize;
        if zi == 0 || zi > D3_COVRAD_BOHR.len() {
            continue;
        }
        for j in 0..n {
            if i == j {
                continue;
            }
            let zj = elements[j] as usize;
            if zj == 0 || zj > D3_COVRAD_BOHR.len() {
                continue;
            }
            let dx = positions_bohr[i][0] - positions_bohr[j][0];
            let dy = positions_bohr[i][1] - positions_bohr[j][1];
            let dz = positions_bohr[i][2] - positions_bohr[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < 1e-10 {
                continue;
            }
            let rc = D3_COVRAD_BOHR[zi - 1] + D3_COVRAD_BOHR[zj - 1];
            cn[i] += 1.0 / (1.0 + (-D3_KCN * (rc / r - 1.0)).exp());
        }
    }
    cn
}

/// Normalized Gaussian reference weights for one atom.
fn d3_weights(cn_atom: f64, z: u8) -> Vec<f64> {
    let zi = z as usize;
    if zi == 0 || zi > D3_NREFS.len() {
        return vec![];
    }
    let nref = D3_NREFS[zi - 1] as usize;
    let mut gw = vec![0.0; nref];
    let mut norm = 0.0;
    for iref in 0..nref {
        let cnref = D3_REFCN[iref + D3_MAX_REF * (zi - 1)];
        let g = (-D3_WF * (cn_atom - cnref) * (cn_atom - cnref)).exp();
        gw[iref] = g;
        norm += g;
    }
    if norm.is_nan() || norm <= 1e-300 {
        // Exceptional fallback: weight 1 on a maximal-CN reference.
        let mut max_cn = f64::NEG_INFINITY;
        for iref in 0..nref {
            max_cn = max_cn.max(D3_REFCN[iref + D3_MAX_REF * (zi - 1)]);
        }
        for iref in 0..nref {
            gw[iref] = if (D3_REFCN[iref + D3_MAX_REF * (zi - 1)] - max_cn).abs() < 1e-12 {
                1.0
            } else {
                0.0
            };
        }
        return gw;
    }
    for g in gw.iter_mut() {
        *g /= norm;
    }
    gw
}

/// Interpolated pair C6 from per-atom reference weights.
fn interpolated_c6(wa: &[f64], wb: &[f64], za: u8, zb: u8) -> f64 {
    let mut c6 = 0.0;
    for (iref, &gwa) in wa.iter().enumerate() {
        if gwa == 0.0 {
            continue;
        }
        for (jref, &gwb) in wb.iter().enumerate() {
            if gwb == 0.0 {
                continue;
            }
            c6 += gwa * gwb * d3_c6_ref(iref, jref, za, zb);
        }
    }
    c6
}

/// Pairwise vdW radius from the D3 table (triangular storage like C6).
fn d3_vdwrad(za: u8, zb: u8) -> Option<f64> {
    let (a, b) = (za as usize, zb as usize);
    if a == 0 || b == 0 {
        return None;
    }
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    let idx = lo + hi * (hi - 1) / 2 - 1;
    D3_VDWRAD_BOHR.get(idx).copied()
}

/// Axilrod-Teller-Muto three-body dispersion energy (Hartree).
///
/// Uses the interpolated pair C6 coefficients, `rs9 = 4/3`, `alp = 16` and
/// the given `s9` scaling (simple-dftd3 convention). HF-3c does not include
/// this term; it is needed to match s-dftd3 totals (F46).
pub fn compute_d3_atm_energy(elements: &[u8], positions_bohr: &[[f64; 3]], s9: f64) -> f64 {
    const RS9: f64 = 4.0 / 3.0;
    const ALP: f64 = 16.0;
    let n = elements.len();
    if n < 3 || s9.abs() < 1e-15 {
        return 0.0;
    }
    let cn = d3_cn(elements, positions_bohr);
    let weights: Vec<Vec<f64>> = elements
        .iter()
        .zip(cn.iter())
        .map(|(&z, &c)| d3_weights(c, z))
        .collect();
    // Interpolated pair C6 matrix (symmetric).
    let mut c6mat = vec![vec![0.0; n]; n];
    for a in 0..n {
        for b in (a + 1)..n {
            let c = interpolated_c6(&weights[a], &weights[b], elements[a], elements[b]);
            c6mat[a][b] = c;
            c6mat[b][a] = c;
        }
    }

    let alp3 = ALP / 3.0;
    let mut energy = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            let r2ij = dist2(positions_bohr, i, j);
            if r2ij < 1e-15 {
                continue;
            }
            for k in (j + 1)..n {
                let r2ik = dist2(positions_bohr, i, k);
                let r2jk = dist2(positions_bohr, j, k);
                if r2ik < 1e-15 || r2jk < 1e-15 {
                    continue;
                }
                let c6ij = c6mat[i][j];
                let c6ik = c6mat[i][k];
                let c6jk = c6mat[j][k];
                if c6ij.abs() < 1e-300 || c6ik.abs() < 1e-300 || c6jk.abs() < 1e-300 {
                    continue;
                }
                let Some(rv_ij) = d3_vdwrad(elements[i], elements[j]) else {
                    continue;
                };
                let Some(rv_ik) = d3_vdwrad(elements[i], elements[k]) else {
                    continue;
                };
                let Some(rv_jk) = d3_vdwrad(elements[j], elements[k]) else {
                    continue;
                };
                let c9 = -s9 * (c6ij * c6ik * c6jk).abs().sqrt();
                let r0 = RS9 * rv_ij * RS9 * rv_ik * RS9 * rv_jk;
                let r2 = r2ij * r2ik * r2jk;
                let r1 = r2.sqrt();
                let r3 = r2 * r1;
                let r5 = r3 * r2;
                let fdmp = 1.0 / (1.0 + 6.0 * (r0 / r1).powf(alp3));
                let ang =
                    0.375 * (r2ij + r2jk - r2ik) * (r2ij - r2jk + r2ik) * (-r2ij + r2jk + r2ik)
                        / r5
                        + 1.0 / r3;
                energy -= ang * fdmp * c9;
            }
        }
    }
    energy
}

/// Coordination numbers from the exponential counting function
/// (`kcn = 16`, D3 covalent radii, simple-dftd3 convention).
///
/// This is the same counting function and radius set the GFN1-xTB Hamiltonian
/// uses through tblite's `cn_count%exp` (validation F44).
pub fn coordination_numbers_exp(elements: &[u8], positions_bohr: &[[f64; 3]]) -> Vec<f64> {
    d3_cn(elements, positions_bohr)
}

fn dist2(positions_bohr: &[[f64; 3]], a: usize, b: usize) -> f64 {
    let dx = positions_bohr[a][0] - positions_bohr[b][0];
    let dy = positions_bohr[a][1] - positions_bohr[b][1];
    let dz = positions_bohr[a][2] - positions_bohr[b][2];
    dx * dx + dy * dy + dz * dz
}

/// Compute D3-BJ dispersion correction with CN-interpolated C6 coefficients.
///
/// Same signature and HF-3c damping parameters as before; the C6 coefficients
/// are now interpolated over coordination-number-dependent reference systems
/// (Grimme D3), validated against s-dftd3 (F46).
/// D3-BJ damping parameters (functional-dependent).
#[derive(Debug, Clone, Copy)]
pub struct D3Params {
    /// s6 scaling (usually 1.0).
    pub s6: f64,
    /// s8 scaling.
    pub s8: f64,
    /// BJ a1.
    pub a1: f64,
    /// BJ a2 (Bohr).
    pub a2: f64,
}

impl D3Params {
    /// HF-3c D3(BJ) parameters (validated F46 vs s-dftd3).
    pub fn hf3c() -> Self {
        Self {
            s6: S6,
            s8: S8,
            a1: A1,
            a2: A2,
        }
    }
    /// GFN1-xTB D3(BJ) parameters (tblite; validation F44).
    pub fn gfn1() -> Self {
        Self {
            s6: 1.0,
            s8: 2.4,
            a1: 0.63,
            a2: 5.0,
        }
    }
}

pub fn compute_d3_energy(elements: &[u8], positions_bohr: &[[f64; 3]]) -> D3Result {
    compute_d3_energy_with_params(elements, positions_bohr, D3Params::hf3c())
}

/// Compute D3-BJ dispersion correction with explicit damping parameters.
pub fn compute_d3_energy_with_params(
    elements: &[u8],
    positions_bohr: &[[f64; 3]],
    params: D3Params,
) -> D3Result {
    let n = elements.len();
    let cn = d3_cn(elements, positions_bohr);
    // Per-atom reference weights (empty when the element is not tabulated).
    let weights: Vec<Vec<f64>> = elements
        .iter()
        .zip(cn.iter())
        .map(|(&z, &c)| d3_weights(c, z))
        .collect();

    let mut e6 = 0.0;
    let mut e8 = 0.0;

    for a in 0..n {
        let za = elements[a];
        if weights[a].is_empty() {
            continue;
        }
        for b in (a + 1)..n {
            let zb = elements[b];
            if weights[b].is_empty() {
                continue;
            }
            let r = distance(positions_bohr, a, b);
            if r < 1e-10 {
                continue;
            }

            // Interpolated C6 over the reference grid.
            let c6 = interpolated_c6(&weights[a], &weights[b], za, zb);
            if c6 <= 0.0 {
                continue;
            }

            // C8/C6 ratio and BJ radius from the D3 r4r2 table. Note the
            // effective quadrupole is sqrt(0.5 * R * sqrt(Z)), mirroring
            // simple-dftd3's `sqrt_z_r4_over_r2` (validation F46).
            let (ra, rb) = (za as usize, zb as usize);
            if ra == 0 || rb == 0 || ra > D3_R4R2.len() || rb > D3_R4R2.len() {
                continue;
            }
            let qa = (0.5 * D3_R4R2[ra - 1] * (ra as f64).sqrt()).sqrt();
            let qb = (0.5 * D3_R4R2[rb - 1] * (rb as f64).sqrt()).sqrt();
            let rrij: f64 = 3.0 * qa * qb;
            let r0 = params.a1 * rrij.sqrt() + params.a2;

            let t6 = 1.0 / (r.powi(6) + r0.powi(6));
            let t8 = 1.0 / (r.powi(8) + r0.powi(8));
            e6 -= params.s6 * c6 * t6;
            e8 -= params.s8 * c6 * rrij * t8;
        }
    }

    D3Result {
        energy: e6 + e8,
        e6,
        e8,
    }
}

fn distance(positions: &[[f64; 3]], a: usize, b: usize) -> f64 {
    let dx = positions[a][0] - positions[b][0];
    let dy = positions[a][1] - positions[b][1];
    let dz = positions[a][2] - positions[b][2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}
#[cfg(test)]
mod tests {
    use super::super::basis::ANG_TO_BOHR;
    use super::*;

    #[test]
    fn test_d3_h2() {
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74 * ANG_TO_BOHR]];
        let result = compute_d3_energy(&[1, 1], &positions);
        // D3 energy should be small and negative
        assert!(result.energy < 0.0, "D3 should be attractive");
        assert!(result.energy.abs() < 0.01, "D3 for H2 should be tiny");
    }

    #[test]
    fn test_d3_ethane_larger() {
        // Two carbons + H atoms → more D3
        let r = 1.54 * ANG_TO_BOHR;
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, r]];
        let result = compute_d3_energy(&[6, 6], &positions);
        assert!(result.energy < 0.0);
    }
}
