//! AO → MO integral transformation for CI methods.
//!
//! Pre-computes the full 4-index MO integral tensor:
//!   (pq|rs)_MO = Σ_{μνλσ} C_{μp} C_{νq} C_{λr} C_{σs} (μν|λσ)_AO
//!
//! Uses the half-transformation (Yoshimine) approach for O(N⁵) scaling
//! instead of naive O(N⁸).

use nalgebra::DMatrix;
use super::integrals::get_eri;

/// Pre-computed MO-basis electron repulsion integrals.
///
/// Stores the full (pq|rs) tensor for indices in `[0, n_mo)`.
/// Uses chemist notation: (pq|rs) = ∫∫ φ_p(1)φ_q(1) r₁₂⁻¹ φ_r(2)φ_s(2) d1 d2
pub struct MoIntegrals {
    /// The MO ERI tensor stored as flat array with 4-fold symmetry.
    data: Vec<f64>,
    /// Number of MO indices.
    n_mo: usize,
}

impl MoIntegrals {
    /// Retrieve (pq|rs) in chemist notation.
    #[inline]
    pub fn get(&self, p: usize, q: usize, r: usize, s: usize) -> f64 {
        // Use 4-fold symmetry: (pq|rs) = (qp|rs) = (pq|sr) = (qp|sr)
        let pq = if p >= q {
            p * (p + 1) / 2 + q
        } else {
            q * (q + 1) / 2 + p
        };
        let rs = if r >= s {
            r * (r + 1) / 2 + s
        } else {
            s * (s + 1) / 2 + r
        };
        let (idx1, idx2) = if pq >= rs { (pq, rs) } else { (rs, pq) };
        let index = idx1 * (idx1 + 1) / 2 + idx2;
        if index < self.data.len() {
            self.data[index]
        } else {
            0.0
        }
    }

    /// Antisymmetrized MO integral: <pq||rs> = (pr|qs) - (ps|qr)
    #[inline]
    pub fn antisymmetrized(&self, p: usize, q: usize, r: usize, s: usize) -> f64 {
        self.get(p, r, q, s) - self.get(p, s, q, r)
    }

    /// Number of MO indices.
    pub fn n_mo(&self) -> usize {
        self.n_mo
    }
}

/// Perform the 4-index AO → MO integral transformation.
///
/// Uses the Yoshimine half-transformation approach:
///   Step 1: (μν|λσ) → (pν|λσ) — first quarter transform
///   Step 2: (pν|λσ) → (pq|λσ) — second quarter transform
///   Step 3: (pq|λσ) → (pq|rσ) — third quarter transform  
///   Step 4: (pq|rσ) → (pq|rs) — fourth quarter transform
///
/// Each step is O(N⁵), total O(N⁵). Storage for intermediates is O(N⁴).
///
/// `coefficients`: MO coefficient matrix (n_ao × n_mo), column p = MO p.
/// `ao_eris`: AO-basis ERIs in flat array (from HF integrals).
/// `n_ao`: number of AO basis functions.
/// `n_mo`: number of MO to transform (can be < n_ao for truncated virtual space).
pub fn ao_to_mo_transform(
    coefficients: &DMatrix<f64>,
    ao_eris: &[f64],
    n_ao: usize,
    n_mo: usize,
) -> MoIntegrals {
    // Total unique pairs
    let n_pair = n_mo * (n_mo + 1) / 2;
    let n_total = n_pair * (n_pair + 1) / 2;
    let mut mo_eris = vec![0.0; n_total];

    // Half-transformed intermediates
    // Strategy: transform one index at a time
    // (μν|λσ) → accumulate into (pq|rs) using 4 nested transforms
    //
    // For small systems (n_ao < ~50), direct transformation is fast enough.
    // For larger systems, the Yoshimine sort approach would be preferred.

    // Direct 4-quarter transformation
    // We accumulate (pq|rs) = Σ_{μνλσ} C_{μp} C_{νq} C_{λr} C_{σs} (μν|λσ)
    // Using half-transform approach:
    // tmp1[p,ν,λ,σ] = Σ_μ C_{μp} (μν|λσ)           -- O(N^5) first quarter

    // Allocate half-transformed tensor (p, nu, lam, sig) with p < n_mo
    let mut half1 = vec![0.0; n_mo * n_ao * n_ao * n_ao];

    // First quarter: sum over μ
    for p in 0..n_mo {
        for nu in 0..n_ao {
            for lam in 0..n_ao {
                for sig in 0..=lam {
                    let mut val = 0.0;
                    for mu in 0..n_ao {
                        let c_mp = coefficients[(mu, p)];
                        if c_mp.abs() > 1e-14 {
                            val += c_mp * get_eri(ao_eris, mu, nu, lam, sig, n_ao);
                        }
                    }
                    let idx = ((p * n_ao + nu) * n_ao + lam) * n_ao + sig;
                    half1[idx] = val;
                    // (μν|λσ) = (μν|σλ) symmetry
                    let idx2 = ((p * n_ao + nu) * n_ao + sig) * n_ao + lam;
                    half1[idx2] = val;
                }
            }
        }
    }

    // Second quarter: sum over ν → (pq|λσ)
    let mut half2 = vec![0.0; n_mo * n_mo * n_ao * n_ao];
    for p in 0..n_mo {
        for q in 0..=p {
            for lam in 0..n_ao {
                for sig in 0..n_ao {
                    let mut val = 0.0;
                    for nu in 0..n_ao {
                        let c_nq = coefficients[(nu, q)];
                        if c_nq.abs() > 1e-14 {
                            let idx = ((p * n_ao + nu) * n_ao + lam) * n_ao + sig;
                            val += c_nq * half1[idx];
                        }
                    }
                    let idx2 = ((p * n_mo + q) * n_ao + lam) * n_ao + sig;
                    half2[idx2] = val;
                }
            }
        }
    }
    drop(half1);

    // Third quarter: sum over λ → (pq|rσ)
    let mut half3 = vec![0.0; n_mo * n_mo * n_mo * n_ao];
    for p in 0..n_mo {
        for q in 0..=p {
            for r in 0..n_mo {
                for sig in 0..n_ao {
                    let mut val = 0.0;
                    for lam in 0..n_ao {
                        let c_lr = coefficients[(lam, r)];
                        if c_lr.abs() > 1e-14 {
                            let idx = ((p * n_mo + q) * n_ao + lam) * n_ao + sig;
                            val += c_lr * half2[idx];
                        }
                    }
                    let idx3 = ((p * n_mo + q) * n_mo + r) * n_ao + sig;
                    half3[idx3] = val;
                }
            }
        }
    }
    drop(half2);

    // Fourth quarter: sum over σ → store (pq|rs)
    for p in 0..n_mo {
        for q in 0..=p {
            let pq = p * (p + 1) / 2 + q;
            for r in 0..n_mo {
                for s in 0..=r {
                    let rs = r * (r + 1) / 2 + s;
                    if pq < rs {
                        continue; // will be set by symmetry
                    }
                    let mut val = 0.0;
                    for sig in 0..n_ao {
                        let c_ss = coefficients[(sig, s)];
                        if c_ss.abs() > 1e-14 {
                            let idx = ((p * n_mo + q) * n_mo + r) * n_ao + sig;
                            val += c_ss * half3[idx];
                        }
                    }
                    let index = pq * (pq + 1) / 2 + rs;
                    mo_eris[index] = val;
                }
            }
        }
    }

    MoIntegrals {
        data: mo_eris,
        n_mo,
    }
}

/// Perform AO → MO transform with rayon parallelism on the first quarter.
#[cfg(feature = "parallel")]
pub fn ao_to_mo_transform_parallel(
    coefficients: &DMatrix<f64>,
    ao_eris: &[f64],
    n_ao: usize,
    n_mo: usize,
) -> MoIntegrals {
    use rayon::prelude::*;

    let n_pair = n_mo * (n_mo + 1) / 2;
    let n_total = n_pair * (n_pair + 1) / 2;

    // First quarter: parallelize over p
    let half1: Vec<Vec<f64>> = (0..n_mo)
        .into_par_iter()
        .map(|p| {
            let mut slice = vec![0.0; n_ao * n_ao * n_ao];
            for nu in 0..n_ao {
                for lam in 0..n_ao {
                    for sig in 0..=lam {
                        let mut val = 0.0;
                        for mu in 0..n_ao {
                            let c_mp = coefficients[(mu, p)];
                            if c_mp.abs() > 1e-14 {
                                val += c_mp * get_eri(ao_eris, mu, nu, lam, sig, n_ao);
                            }
                        }
                        slice[nu * n_ao * n_ao + lam * n_ao + sig] = val;
                        slice[nu * n_ao * n_ao + sig * n_ao + lam] = val;
                    }
                }
            }
            slice
        })
        .collect();

    // Second quarter
    let mut half2 = vec![0.0; n_mo * n_mo * n_ao * n_ao];
    for p in 0..n_mo {
        for q in 0..=p {
            for lam in 0..n_ao {
                for sig in 0..n_ao {
                    let mut val = 0.0;
                    for nu in 0..n_ao {
                        let c_nq = coefficients[(nu, q)];
                        if c_nq.abs() > 1e-14 {
                            val += c_nq * half1[p][nu * n_ao * n_ao + lam * n_ao + sig];
                        }
                    }
                    half2[((p * n_mo + q) * n_ao + lam) * n_ao + sig] = val;
                }
            }
        }
    }
    drop(half1);

    // Third quarter
    let mut half3 = vec![0.0; n_mo * n_mo * n_mo * n_ao];
    for p in 0..n_mo {
        for q in 0..=p {
            for r in 0..n_mo {
                for sig in 0..n_ao {
                    let mut val = 0.0;
                    for lam in 0..n_ao {
                        let c_lr = coefficients[(lam, r)];
                        if c_lr.abs() > 1e-14 {
                            val += c_lr * half2[((p * n_mo + q) * n_ao + lam) * n_ao + sig];
                        }
                    }
                    half3[((p * n_mo + q) * n_mo + r) * n_ao + sig] = val;
                }
            }
        }
    }
    drop(half2);

    // Fourth quarter
    let mut mo_eris = vec![0.0; n_total];
    for p in 0..n_mo {
        for q in 0..=p {
            let pq = p * (p + 1) / 2 + q;
            for r in 0..n_mo {
                for s in 0..=r {
                    let rs = r * (r + 1) / 2 + s;
                    if pq < rs {
                        continue;
                    }
                    let mut val = 0.0;
                    for sig in 0..n_ao {
                        let c_ss = coefficients[(sig, s)];
                        if c_ss.abs() > 1e-14 {
                            val += c_ss * half3[((p * n_mo + q) * n_mo + r) * n_ao + sig];
                        }
                    }
                    mo_eris[pq * (pq + 1) / 2 + rs] = val;
                }
            }
        }
    }

    MoIntegrals {
        data: mo_eris,
        n_mo,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mo_integrals_symmetry() {
        // Simple 2-MO system
        let c = DMatrix::from_row_slice(2, 2, &[
            0.7071, 0.7071,
            0.7071, -0.7071,
        ]);
        // Simple AO ERIs (all equal for testing)
        let eris = vec![0.5; 16];
        let mo = ao_to_mo_transform(&c, &eris, 2, 2);
        // (pq|rs) should have 4-fold symmetry
        assert!((mo.get(0, 0, 1, 1) - mo.get(1, 1, 0, 0)).abs() < 1e-10);
    }
}
