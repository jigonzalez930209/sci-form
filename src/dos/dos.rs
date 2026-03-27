//! Density of States (DOS) and Projected DOS (PDOS).
//!
//! Computes total DOS by Gaussian-smearing EHT orbital energies
//! and atom-projected DOS by weighting with Mulliken orbital populations.

use crate::eht::basis::{build_basis, AtomicOrbital};
use crate::eht::overlap::build_overlap_matrix;
use serde::{Deserialize, Serialize};

#[allow(dead_code)]
fn gaussian_value(energy: f64, center: f64, norm: f64, inv_2s2: f64) -> f64 {
    norm * (-(energy - center).powi(2) * inv_2s2).exp()
}

/// Result of a DOS/PDOS calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DosResult {
    /// Energy grid values (eV).
    pub energies: Vec<f64>,
    /// Total DOS values (states/eV).
    pub total_dos: Vec<f64>,
    /// Per-atom PDOS: pdos\[atom_idx\]\[grid_idx\].
    pub pdos: Vec<Vec<f64>>,
    /// Smearing width used (eV).
    pub sigma: f64,
}

/// Compute total density of states from EHT orbital energies.
///
/// `orbital_energies`: eigenvalues from EHT (eV).
/// `sigma`: Gaussian smearing width (eV).
/// `e_min`, `e_max`: energy window.
/// `n_points`: grid resolution.
pub fn compute_dos(
    orbital_energies: &[f64],
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> DosResult {
    let step = (e_max - e_min) / (n_points - 1).max(1) as f64;
    let energies: Vec<f64> = (0..n_points).map(|i| e_min + i as f64 * step).collect();

    let norm = 1.0 / (sigma * (2.0 * std::f64::consts::PI).sqrt());
    let inv_2s2 = 1.0 / (2.0 * sigma * sigma);

    let total_dos: Vec<f64> = energies
        .iter()
        .map(|&e| {
            orbital_energies
                .iter()
                .map(|&ei| norm * (-(e - ei).powi(2) * inv_2s2).exp())
                .sum()
        })
        .collect();

    DosResult {
        energies,
        total_dos,
        pdos: Vec::new(),
        sigma,
    }
}

/// Compute total DOS using rayon over the energy grid.
#[cfg(feature = "parallel")]
pub fn compute_dos_parallel(
    orbital_energies: &[f64],
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> DosResult {
    use rayon::prelude::*;

    let step = (e_max - e_min) / (n_points - 1).max(1) as f64;
    let energies: Vec<f64> = (0..n_points).map(|i| e_min + i as f64 * step).collect();

    let norm = 1.0 / (sigma * (2.0 * std::f64::consts::PI).sqrt());
    let inv_2s2 = 1.0 / (2.0 * sigma * sigma);

    let total_dos: Vec<f64> = energies
        .par_iter()
        .map(|&energy| {
            orbital_energies
                .iter()
                .map(|&center| gaussian_value(energy, center, norm, inv_2s2))
                .sum()
        })
        .collect();

    DosResult {
        energies,
        total_dos,
        pdos: Vec::new(),
        sigma,
    }
}

/// Compute atom-projected DOS from EHT results.
///
/// `elements`: atomic numbers per atom.
/// `positions`: flat \[x0,y0,z0, x1,y1,z1,...\] in Å.
/// `orbital_energies`: eigenvalues from EHT (eV).
/// `coefficients`: coefficients\[orbital\]\[basis\] from EHT.
/// `n_electrons`: number of electrons.
/// `sigma`: Gaussian smearing width (eV).
/// `e_min`, `e_max`: energy window.
/// `n_points`: grid resolution.
#[allow(clippy::too_many_arguments)]
pub fn compute_pdos(
    elements: &[u8],
    positions: &[f64],
    orbital_energies: &[f64],
    coefficients: &[Vec<f64>],
    n_electrons: usize,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> DosResult {
    let n_atoms = elements.len();
    let pos_arr: Vec<[f64; 3]> = positions
        .chunks_exact(3)
        .map(|c| [c[0], c[1], c[2]])
        .collect();
    let basis: Vec<AtomicOrbital> = build_basis(elements, &pos_arr);
    let overlap = build_overlap_matrix(&basis);
    let n_basis = basis.len();

    // Build density-like Mulliken weight per (orbital, atom):
    // w_{k,A} = Σ_{μ∈A} Σ_ν c_{k,μ} S_{μν} c_{k,ν}
    let n_orb = orbital_energies.len().min(coefficients.len());
    let mut orbital_atom_weight = vec![vec![0.0f64; n_atoms]; n_orb];

    for k in 0..n_orb {
        for mu in 0..n_basis {
            if coefficients.len() <= mu || coefficients[mu].len() <= k {
                continue;
            }
            let atom_mu = basis[mu].atom_index;
            let mut w = 0.0;
            for nu in 0..n_basis {
                if coefficients.len() <= nu || coefficients[nu].len() <= k {
                    continue;
                }
                // coefficients[mu][k] = AO mu, MO k
                w += coefficients[mu][k] * overlap[(mu, nu)] * coefficients[nu][k];
            }
            orbital_atom_weight[k][atom_mu] += w;
        }
        // Normalize weights so they sum to 1 per orbital
        let total_w: f64 = orbital_atom_weight[k].iter().sum();
        if total_w.abs() > 1e-12 {
            for a in 0..n_atoms {
                orbital_atom_weight[k][a] /= total_w;
            }
        }
    }

    let step = (e_max - e_min) / (n_points - 1).max(1) as f64;
    let energies: Vec<f64> = (0..n_points).map(|i| e_min + i as f64 * step).collect();

    let norm = 1.0 / (sigma * (2.0 * std::f64::consts::PI).sqrt());
    let inv_2s2 = 1.0 / (2.0 * sigma * sigma);

    // Total DOS
    let total_dos: Vec<f64> = energies
        .iter()
        .map(|&e| {
            (0..n_orb)
                .map(|k| norm * (-(e - orbital_energies[k]).powi(2) * inv_2s2).exp())
                .sum()
        })
        .collect();

    // Per-atom PDOS
    let mut pdos = vec![vec![0.0f64; n_points]; n_atoms];
    for a in 0..n_atoms {
        for (gi, &e) in energies.iter().enumerate() {
            let mut val = 0.0;
            for k in 0..n_orb {
                let gauss = norm * (-(e - orbital_energies[k]).powi(2) * inv_2s2).exp();
                val += orbital_atom_weight[k][a] * gauss;
            }
            pdos[a][gi] = val;
        }
    }

    let _ = n_electrons; // used contextually; weight already normalized via Mulliken

    DosResult {
        energies,
        total_dos,
        pdos,
        sigma,
    }
}

/// Compute atom-projected DOS using rayon for orbital weights and atom-grid accumulation.
#[cfg(feature = "parallel")]
#[allow(clippy::too_many_arguments)]
pub fn compute_pdos_parallel(
    elements: &[u8],
    positions: &[f64],
    orbital_energies: &[f64],
    coefficients: &[Vec<f64>],
    n_electrons: usize,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> DosResult {
    use rayon::prelude::*;

    let n_atoms = elements.len();
    let pos_arr: Vec<[f64; 3]> = positions
        .chunks_exact(3)
        .map(|c| [c[0], c[1], c[2]])
        .collect();
    let basis: Vec<AtomicOrbital> = build_basis(elements, &pos_arr);
    let overlap = build_overlap_matrix(&basis);
    let n_basis = basis.len();
    let n_orb = orbital_energies.len().min(coefficients.len());

    let orbital_atom_weight: Vec<Vec<f64>> = (0..n_orb)
        .into_par_iter()
        .map(|k| {
            let mut weights = vec![0.0f64; n_atoms];
            for mu in 0..n_basis {
                if coefficients.len() <= mu || coefficients[mu].len() <= k {
                    continue;
                }
                let atom_mu = basis[mu].atom_index;
                let mut weight = 0.0;
                for nu in 0..n_basis {
                    if coefficients.len() <= nu || coefficients[nu].len() <= k {
                        continue;
                    }
                    weight += coefficients[mu][k] * overlap[(mu, nu)] * coefficients[nu][k];
                }
                weights[atom_mu] += weight;
            }

            let total_weight: f64 = weights.iter().sum();
            if total_weight.abs() > 1e-12 {
                for weight in &mut weights {
                    *weight /= total_weight;
                }
            }
            weights
        })
        .collect();

    let step = (e_max - e_min) / (n_points - 1).max(1) as f64;
    let energies: Vec<f64> = (0..n_points).map(|i| e_min + i as f64 * step).collect();

    let norm = 1.0 / (sigma * (2.0 * std::f64::consts::PI).sqrt());
    let inv_2s2 = 1.0 / (2.0 * sigma * sigma);

    let total_dos: Vec<f64> = energies
        .par_iter()
        .map(|&energy| {
            (0..n_orb)
                .map(|k| gaussian_value(energy, orbital_energies[k], norm, inv_2s2))
                .sum()
        })
        .collect();

    let pdos: Vec<Vec<f64>> = (0..n_atoms)
        .into_par_iter()
        .map(|atom_index| {
            energies
                .iter()
                .map(|&energy| {
                    (0..n_orb)
                        .map(|k| {
                            orbital_atom_weight[k][atom_index]
                                * gaussian_value(energy, orbital_energies[k], norm, inv_2s2)
                        })
                        .sum()
                })
                .collect()
        })
        .collect();

    let _ = n_electrons;

    DosResult {
        energies,
        total_dos,
        pdos,
        sigma,
    }
}

/// Compute mean-squared error between two DOS curves.
///
/// Both curves must have the same length.  Useful for comparing against
/// reference DOS (e.g. Multiwfn output).
pub fn dos_mse(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "DOS curves must have same length");
    let n = a.len() as f64;
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).powi(2))
        .sum::<f64>()
        / n
}

/// Serialize DOS/PDOS result to JSON for web visualization.
///
/// Format:
/// ```json
/// {
///   "energies": [...],
///   "total_dos": [...],
///   "sigma": 0.3,
///   "pdos": { "0": [...], "1": [...], ... }
/// }
/// ```
pub fn export_dos_json(result: &DosResult) -> String {
    let mut json = String::from("{");
    json.push_str("\"energies\":[");
    for (i, e) in result.energies.iter().enumerate() {
        if i > 0 {
            json.push(',');
        }
        json.push_str(&format!("{:.6}", e));
    }
    json.push_str("],\"total_dos\":[");
    for (i, d) in result.total_dos.iter().enumerate() {
        if i > 0 {
            json.push(',');
        }
        json.push_str(&format!("{:.6}", d));
    }
    json.push_str(&format!("],\"sigma\":{:.6}", result.sigma));
    if !result.pdos.is_empty() {
        json.push_str(",\"pdos\":{");
        for (a, pdos_a) in result.pdos.iter().enumerate() {
            if a > 0 {
                json.push(',');
            }
            json.push_str(&format!("\"{}\":[", a));
            for (i, v) in pdos_a.iter().enumerate() {
                if i > 0 {
                    json.push(',');
                }
                json.push_str(&format!("{:.6}", v));
            }
            json.push(']');
        }
        json.push('}');
    }
    json.push('}');
    json
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dos_single_level() {
        // Single orbital at 0 eV → DOS should peak at 0.
        let res = compute_dos(&[0.0], 0.1, -1.0, 1.0, 201);
        assert_eq!(res.energies.len(), 201);
        assert_eq!(res.total_dos.len(), 201);
        // Peak should be at mid-point (index 100)
        let peak_idx = res
            .total_dos
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap()
            .0;
        assert_eq!(peak_idx, 100);
    }

    #[test]
    fn test_dos_integral_approx_one() {
        // Integral of DOS for one level ≈ 1 (normalized Gaussian).
        let res = compute_dos(&[0.0], 0.2, -3.0, 3.0, 1001);
        let de = (3.0 - (-3.0)) / 1000.0;
        let integral: f64 = res.total_dos.iter().sum::<f64>() * de;
        assert!((integral - 1.0).abs() < 0.01, "integral = {integral}");
    }

    #[test]
    fn test_dos_two_peaks() {
        let res = compute_dos(&[-5.0, 5.0], 0.3, -10.0, 10.0, 501);
        // Should have two peaks, one near index ~125, one near ~375
        let mid = res.total_dos[250];
        let left_peak = res.total_dos[125];
        let right_peak = res.total_dos[375];
        assert!(left_peak > mid * 5.0);
        assert!(right_peak > mid * 5.0);
    }

    #[test]
    fn test_pdos_h2() {
        // H₂: two atoms should have symmetric PDOS.
        let elements = vec![1u8, 1];
        let pos_arr = vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let positions: Vec<f64> = pos_arr.iter().flat_map(|p| p.iter().copied()).collect();
        let eht = crate::eht::solve_eht(&elements, &pos_arr, None).unwrap();
        let res = compute_pdos(
            &elements,
            &positions,
            &eht.energies,
            &eht.coefficients,
            eht.n_electrons,
            0.2,
            -20.0,
            5.0,
            201,
        );
        assert_eq!(res.pdos.len(), 2);
        // Both H atoms should have nearly equal PDOS at peaks.
        // The peak region is where DOS > 10% of max.
        let peak_val = res.pdos[0].iter().cloned().fold(0.0f64, f64::max);
        let threshold = peak_val * 0.1;
        for i in 0..201 {
            if res.pdos[0][i].abs() > threshold || res.pdos[1][i].abs() > threshold {
                let diff = (res.pdos[0][i] - res.pdos[1][i]).abs();
                let avg = (res.pdos[0][i].abs() + res.pdos[1][i].abs()) / 2.0;
                assert!(
                    diff < avg * 0.05 + 1e-6,
                    "PDOS mismatch at grid point {i}: {} vs {} (peak={})",
                    res.pdos[0][i],
                    res.pdos[1][i],
                    peak_val
                );
            }
        }
    }

    #[test]
    fn test_pdos_sums_to_total() {
        // Sum of PDOS over all atoms ≈ total DOS.
        let elements = vec![8u8, 1, 1];
        let pos_arr = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
        let positions: Vec<f64> = pos_arr.iter().flat_map(|p| p.iter().copied()).collect();
        let eht = crate::eht::solve_eht(&elements, &pos_arr, None).unwrap();
        let res = compute_pdos(
            &elements,
            &positions,
            &eht.energies,
            &eht.coefficients,
            eht.n_electrons,
            0.3,
            -30.0,
            5.0,
            201,
        );
        for i in 0..201 {
            let pdos_sum: f64 = res.pdos.iter().map(|p| p[i]).sum();
            let diff = (pdos_sum - res.total_dos[i]).abs();
            assert!(
                diff < res.total_dos[i].abs() * 0.05 + 1e-10,
                "PDOS sum {pdos_sum} vs total {} at grid {i}",
                res.total_dos[i]
            );
        }
    }

    #[test]
    fn test_dos_mse_identical() {
        let a = vec![1.0, 2.0, 3.0, 4.0];
        assert!((dos_mse(&a, &a)) < 1e-15);
    }

    #[test]
    fn test_dos_mse_known() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![1.1, 1.9, 3.2];
        // MSE = (0.01 + 0.01 + 0.04) / 3 = 0.02
        assert!((dos_mse(&a, &b) - 0.02).abs() < 1e-10);
    }

    #[test]
    fn test_export_dos_json_roundtrip() {
        let res = compute_dos(&[0.0, -5.0], 0.3, -10.0, 5.0, 51);
        let json = export_dos_json(&res);

        // Should be valid JSON
        let parsed: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");
        assert!(parsed["energies"].is_array());
        assert!(parsed["total_dos"].is_array());
        assert_eq!(parsed["energies"].as_array().unwrap().len(), 51);
        assert_eq!(parsed["total_dos"].as_array().unwrap().len(), 51);
    }

    #[test]
    fn test_export_pdos_json() {
        let elements = vec![1u8, 1];
        let pos_arr = vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let positions: Vec<f64> = pos_arr.iter().flat_map(|p| p.iter().copied()).collect();
        let eht = crate::eht::solve_eht(&elements, &pos_arr, None).unwrap();
        let res = compute_pdos(
            &elements,
            &positions,
            &eht.energies,
            &eht.coefficients,
            eht.n_electrons,
            0.2,
            -20.0,
            5.0,
            51,
        );
        let json = export_dos_json(&res);
        let parsed: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");
        assert!(parsed["pdos"].is_object());
        assert!(parsed["pdos"]["0"].is_array());
        assert!(parsed["pdos"]["1"].is_array());
    }
}
