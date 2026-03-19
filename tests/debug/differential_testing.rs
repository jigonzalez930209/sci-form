#![allow(
    unused_imports,
    unused_variables,
    dead_code,
    clippy::unnecessary_cast,
    clippy::needless_range_loop,
    clippy::manual_repeat_n,
    clippy::manual_str_repeat,
    clippy::manual_is_multiple_of,
    clippy::redundant_field_names,
    clippy::useless_vec,
    clippy::single_range_in_vec_init
)]
use rand::SeedableRng;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs;

#[derive(Deserialize, Debug)]
struct OracleAtom {
    element: u8,
    x: f32,
    y: f32,
    z: f32,
    formal_charge: i8,
    hybridization: String,
}

#[derive(Deserialize, Debug)]
struct OracleBond {
    start: usize,
    end: usize,
    order: String,
}

#[derive(Deserialize, Debug)]
struct OracleMolecule {
    smiles: String,
    atoms: Vec<OracleAtom>,
    bonds: Vec<OracleBond>,
}

#[derive(Deserialize, Debug)]
struct CsdTorsion {
    atoms: Vec<usize>,
    v: Vec<f64>,
    signs: Vec<i32>,
}

#[test]
fn test_parse_reference_data() {
    let ref_file = std::env::var("REF_FILE")
        .unwrap_or_else(|_| "tests/fixtures/reference_coords.json".to_string());
    let data = fs::read_to_string(&ref_file).expect("Should be able to read reference JSON");
    let mut molecules: Vec<OracleMolecule> =
        serde_json::from_str(&data).expect("JSON was not well-formatted");

    assert!(!molecules.is_empty());

    // Shuffle and pick 100 random molecules (deterministic seed for reproducibility)
    use rand::seq::SliceRandom;
    let mut rng = rand::rngs::StdRng::seed_from_u64(123);
    molecules.shuffle(&mut rng);
    molecules.truncate(100);

    // Load CSD torsion parameters extracted from RDKit
    let csd_torsions: HashMap<String, Vec<CsdTorsion>> = {
        let path = "tests/fixtures/torsion_params.json";
        match fs::read_to_string(path) {
            Ok(data) => serde_json::from_str(&data).unwrap_or_default(),
            Err(_) => HashMap::new(),
        }
    };
    println!(
        "Loaded CSD torsion params for {} molecules",
        csd_torsions.len()
    );

    let mut total_rmsd = 0.0;
    let mut count = 0;
    let mut oracle_total = 0.0f32;
    let mut oracle_max = 0.0f32;
    let mut oracle_above = 0u32;
    let mut select_max = 0.0f32;
    let mut select_above = 0u32;
    // Track multiple selection strategies
    let mut strat_above = [0u32; 6]; // above 0.5 count for each strategy
    let mut strat_total = [0.0f32; 6]; // total RMSD
    let strat_names = [
        "BV",
        "FF_E",
        "CSD_full",
        "BV+CSD_full",
        "BV+FF_E",
        "CSD_full2",
    ];

    for mol in molecules {
        assert!(!mol.atoms.is_empty());
        let n = mol.atoms.len();

        // Convert OracleMolecule to our format
        let mut our_mol = sci_form::graph::Molecule::new(&mol.smiles);

        let mut node_indices = Vec::new();

        for atom in &mol.atoms {
            let hybridization = match atom.hybridization.as_str() {
                "SP" => sci_form::graph::Hybridization::SP,
                "SP2" => sci_form::graph::Hybridization::SP2,
                "SP3" => sci_form::graph::Hybridization::SP3,
                "SP3D" => sci_form::graph::Hybridization::SP3D,
                "SP3D2" => sci_form::graph::Hybridization::SP3D2,
                _ => sci_form::graph::Hybridization::Unknown,
            };

            let new_atom = sci_form::graph::Atom {
                element: atom.element,
                position: nalgebra::Vector3::zeros(), // Ignore RDKit coords when building initial graph
                charge: 0.0,                          // Not parsed yet
                formal_charge: atom.formal_charge,
                hybridization: hybridization,
                chiral_tag: sci_form::graph::ChiralType::Unspecified,
                explicit_h: if atom.element == 1 || atom.element == 0 {
                    1
                } else {
                    0
                },
            };
            let n_idx = our_mol.add_atom(new_atom);
            node_indices.push(n_idx);
        }

        for bond in &mol.bonds {
            let order = match bond.order.as_str() {
                "DOUBLE" => sci_form::graph::BondOrder::Double,
                "TRIPLE" => sci_form::graph::BondOrder::Triple,
                "AROMATIC" => sci_form::graph::BondOrder::Aromatic,
                _ => sci_form::graph::BondOrder::Single,
            };

            our_mol.add_bond(
                node_indices[bond.start],
                node_indices[bond.end],
                sci_form::graph::Bond {
                    order,
                    stereo: sci_form::graph::BondStereo::None,
                },
            );
        }

        // 1. Calculate Bounds
        let bounds = sci_form::distgeom::calculate_bounds_matrix(&our_mol);
        let smoothed = sci_form::distgeom::smooth_bounds_matrix(bounds.clone());

        // 2. Try seeds with per-seed 3D FF
        let mut all_refined: Vec<nalgebra::DMatrix<f32>> = Vec::new();
        let mut all_bv: Vec<f32> = Vec::new();
        let mut all_3dff_e: Vec<f32> = Vec::new();
        let mut all_csd_full_e: Vec<f32> = Vec::new();
        let mut best_oracle_rmsd = f32::MAX;

        // Build CSD torsion contribs for this molecule (if available)
        let csd_contribs: Vec<sci_form::forcefield::M6TorsionContrib> =
            if let Some(torsions) = csd_torsions.get(&mol.smiles) {
                torsions
                    .iter()
                    .map(|t| {
                        let mut signs = [0.0f64; 6];
                        let mut v = [0.0f64; 6];
                        for k in 0..6 {
                            signs[k] = t.signs[k] as f64;
                            v[k] = t.v[k] as f64;
                        }
                        sci_form::forcefield::M6TorsionContrib {
                            i: t.atoms[0],
                            j: t.atoms[1],
                            k: t.atoms[2],
                            l: t.atoms[3],
                            signs,
                            v,
                        }
                    })
                    .collect()
            } else {
                Vec::new()
            };
        let has_csd = !csd_contribs.is_empty();

        for seed in 0..500u64 {
            let mut rng = rand::rngs::StdRng::seed_from_u64(seed * 7 + 13);
            let dists = sci_form::distgeom::pick_random_distances(&mut rng, &smoothed);
            let metric = sci_form::distgeom::compute_metric_matrix(&dists);
            let mut coords4d = sci_form::distgeom::generate_nd_coordinates(&mut rng, &metric, 4);

            sci_form::forcefield::minimize_bounds_lbfgs_ex(
                &mut coords4d,
                &smoothed,
                &[],
                400,
                1e-4,
                5.0,
                0.1,
            );
            sci_form::forcefield::minimize_bounds_lbfgs_ex(
                &mut coords4d,
                &smoothed,
                &[],
                200,
                1e-4,
                5.0,
                1.0,
            );

            let coords3d = coords4d.columns(0, 3).into_owned();
            // 3D FF minimization: use CSD torsion preferences if available
            let refined = if has_csd {
                let mut ff = sci_form::forcefield::build_etkdg_3d_ff(
                    &our_mol,
                    &coords3d.map(|v| v as f64),
                    &smoothed,
                );
                ff.torsion_contribs = csd_contribs.clone();
                ff.use_m6_torsions = false;
                sci_form::forcefield::minimize_etkdg_3d_with_ff(&our_mol, &coords3d, &ff, 200, 1e-4)
            } else {
                sci_form::forcefield::minimize_etkdg_3d(&our_mol, &coords3d, &smoothed, 200, 1e-4)
            };

            let e = sci_form::forcefield::bounds_violation_energy(&refined, &smoothed);
            let ff_e = {
                let ff = sci_form::forcefield::build_etkdg_3d_ff(
                    &our_mol,
                    &refined.map(|v| v as f64),
                    &smoothed,
                );
                sci_form::forcefield::etkdg_3d_energy(&refined, &our_mol, &ff)
            };

            // Full CSD FF energy (distances + OOP + torsions)
            let csd_full_e = if has_csd {
                let mut ff = sci_form::forcefield::build_etkdg_3d_ff(
                    &our_mol,
                    &refined.map(|v| v as f64),
                    &smoothed,
                );
                ff.torsion_contribs = csd_contribs.clone();
                ff.use_m6_torsions = false;
                sci_form::forcefield::etkdg_3d_energy(&refined, &our_mol, &ff)
            } else {
                ff_e
            };

            all_refined.push(refined.clone());
            all_bv.push(e);
            all_3dff_e.push(ff_e);
            all_csd_full_e.push(csd_full_e);

            // Oracle
            let mut sq = 0.0;
            let mut np = 0;
            for i in 0..n {
                for j in (i + 1)..n {
                    let dr = ((mol.atoms[i].x - mol.atoms[j].x).powi(2)
                        + (mol.atoms[i].y - mol.atoms[j].y).powi(2)
                        + (mol.atoms[i].z - mol.atoms[j].z).powi(2))
                    .sqrt();
                    let du = ((refined[(i, 0)] - refined[(j, 0)]).powi(2)
                        + (refined[(i, 1)] - refined[(j, 1)]).powi(2)
                        + (refined[(i, 2)] - refined[(j, 2)]).powi(2))
                    .sqrt();
                    sq += (dr - du).powi(2);
                    np += 1;
                }
            }
            let r = if np > 0 { (sq / np as f32).sqrt() } else { 0.0 };
            if r < best_oracle_rmsd {
                best_oracle_rmsd = r;
            }
        }

        // === Multiple selection strategies ===
        let num_confs = all_refined.len();
        let npairs = n * (n - 1) / 2;

        // Precompute pairwise distance matrices for all conformers
        let mut all_pair_dists: Vec<Vec<f32>> = Vec::with_capacity(num_confs);
        for ci in 0..num_confs {
            let mut pd = vec![0.0f32; npairs];
            let mut idx = 0;
            for ai in 0..n {
                for aj in (ai + 1)..n {
                    pd[idx] = ((all_refined[ci][(ai, 0)] - all_refined[ci][(aj, 0)]).powi(2)
                        + (all_refined[ci][(ai, 1)] - all_refined[ci][(aj, 1)]).powi(2)
                        + (all_refined[ci][(ai, 2)] - all_refined[ci][(aj, 2)]).powi(2))
                    .sqrt();
                    idx += 1;
                }
            }
            all_pair_dists.push(pd);
        }

        // BV threshold
        let bv_threshold = {
            let mut sorted_bv = all_bv.clone();
            sorted_bv.sort_by(|a, b| a.partial_cmp(b).unwrap());
            sorted_bv[num_confs / 4] * 3.0
        };

        // Helper: compute RMSD for a given conformer against reference
        let compute_rmsd = |ci: usize| -> f32 {
            let mut sq = 0.0f32;
            let mut np = 0;
            for i in 0..n {
                for j in (i + 1)..n {
                    let dr = ((mol.atoms[i].x - mol.atoms[j].x).powi(2)
                        + (mol.atoms[i].y - mol.atoms[j].y).powi(2)
                        + (mol.atoms[i].z - mol.atoms[j].z).powi(2))
                    .sqrt();
                    let du = all_pair_dists[ci][np];
                    sq += (dr - du).powi(2);
                    np += 1;
                }
            }
            if np > 0 {
                (sq / np as f32).sqrt()
            } else {
                0.0
            }
        };

        // Strategy 0: lowest BV
        let strat0_idx = all_bv
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap()
            .0;

        // Strategy 1: lowest FF energy (no torsions)
        let strat1_idx = all_3dff_e
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap()
            .0;

        // Strategy 2: CSD full FF energy (includes CSD torsions + distances + OOP)
        let strat2_idx = all_csd_full_e
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap()
            .0;

        // Strategy 3: BV-filtered + lowest CSD full energy
        let strat3_idx = {
            let mut best_i = 0;
            let mut best_e = f32::MAX;
            for ci in 0..num_confs {
                if all_bv[ci] > bv_threshold {
                    continue;
                }
                if all_csd_full_e[ci] < best_e {
                    best_e = all_csd_full_e[ci];
                    best_i = ci;
                }
            }
            best_i
        };

        // Strategy 4: FF_E filtered by BV threshold
        let strat4_idx = {
            let mut best_ci = 0;
            let mut best_e = f32::MAX;
            for ci in 0..num_confs {
                if all_bv[ci] > bv_threshold {
                    continue;
                }
                if all_3dff_e[ci] < best_e {
                    best_e = all_3dff_e[ci];
                    best_ci = ci;
                }
            }
            best_ci
        };

        // Strategy 5: lowest CSD full energy (same as strat2, placeholder)
        let strat5_idx = strat2_idx;

        let strat_indices = [
            strat0_idx, strat1_idx, strat2_idx, strat3_idx, strat4_idx, strat5_idx,
        ];

        // Evaluate all strategies
        for (si, &sel_idx) in strat_indices.iter().enumerate() {
            let rmsd = compute_rmsd(sel_idx);
            strat_total[si] += rmsd;
            if rmsd > 0.5 {
                strat_above[si] += 1;
            }
        }

        // Use strategy 0 (BV) for detailed per-molecule output
        let best_coords = all_refined[strat0_idx].clone();
        let minimized_coords = best_coords.clone();
        assert_eq!(minimized_coords.nrows(), n);

        // Compare pairwise distances decomposed by topological distance
        let mut diff_min_sq_sum = 0.0f32;
        let mut pairs = 0;
        let mut short_sq_sum = 0.0f32;
        let mut short_pairs = 0;
        let mut long_sq_sum = 0.0f32;
        let mut long_pairs = 0;
        for i in 0..n {
            for j in (i + 1)..n {
                let dx_r = mol.atoms[i].x - mol.atoms[j].x;
                let dy_r = mol.atoms[i].y - mol.atoms[j].y;
                let dz_r = mol.atoms[i].z - mol.atoms[j].z;
                let rdkit_dist = (dx_r * dx_r + dy_r * dy_r + dz_r * dz_r).sqrt();

                let dx_min = minimized_coords[(i, 0)] - minimized_coords[(j, 0)];
                let dy_min = minimized_coords[(i, 1)] - minimized_coords[(j, 1)];
                let dz_min = minimized_coords[(i, 2)] - minimized_coords[(j, 2)];
                let min_dist = (dx_min * dx_min + dy_min * dy_min + dz_min * dz_min).sqrt();

                let diff_min = rdkit_dist - min_dist;
                diff_min_sq_sum += diff_min * diff_min;
                pairs += 1;
                let upper = smoothed[(i, j)];
                let lower = smoothed[(j, i)];
                let range = upper - lower;
                if range < 0.3 {
                    // Tightly constrained: bonded or 1-3
                    short_sq_sum += diff_min * diff_min;
                    short_pairs += 1;
                } else {
                    long_sq_sum += diff_min * diff_min;
                    long_pairs += 1;
                }
            }
        }

        let dist_min_rmsd = if pairs > 0 {
            (diff_min_sq_sum / pairs as f32).sqrt()
        } else {
            0.0
        };
        let short_rmsd = if short_pairs > 0 {
            (short_sq_sum / short_pairs as f32).sqrt()
        } else {
            0.0
        };
        let long_rmsd = if long_pairs > 0 {
            (long_sq_sum / long_pairs as f32).sqrt()
        } else {
            0.0
        };

        total_rmsd += dist_min_rmsd;
        oracle_total += best_oracle_rmsd;
        if best_oracle_rmsd > oracle_max {
            oracle_max = best_oracle_rmsd;
        }
        if best_oracle_rmsd > 0.5 {
            oracle_above += 1;
        }
        if dist_min_rmsd > select_max {
            select_max = dist_min_rmsd;
        }
        if dist_min_rmsd > 0.5 {
            select_above += 1;
        }
        // Print every molecule's RMSD with decomposition
        {
            println!(
                "Mol {:3} ({:40}) -> RMSD: {:.3} Å  oracle: {:.3}  short: {:.3} ({}) long: {:.3} ({})",
                count, &mol.smiles[..mol.smiles.len().min(40)], dist_min_rmsd,
                best_oracle_rmsd,
                short_rmsd, short_pairs, long_rmsd, long_pairs
            );
            // No per-pair diagnostic to keep output manageable
            // But show diagnostic for oracle failures (oracle > 0.48)
            if best_oracle_rmsd > 0.48 {
                println!(
                    "  ORACLE FAIL: best_oracle={:.3}, analyzing bounds vs ref:",
                    best_oracle_rmsd
                );
                let mut viol_count = 0;
                for ii in 0..n {
                    for jj in (ii + 1)..n {
                        let lb = smoothed[(jj, ii)];
                        let ub = smoothed[(ii, jj)];
                        let dx_r = mol.atoms[ii].x - mol.atoms[jj].x;
                        let dy_r = mol.atoms[ii].y - mol.atoms[jj].y;
                        let dz_r = mol.atoms[ii].z - mol.atoms[jj].z;
                        let rdkit_d = (dx_r * dx_r + dy_r * dy_r + dz_r * dz_r).sqrt() as f64;
                        let viol = if rdkit_d < lb {
                            lb - rdkit_d
                        } else if rdkit_d > ub {
                            rdkit_d - ub
                        } else {
                            0.0
                        };
                        if viol > 0.03 {
                            viol_count += 1;
                            let range = ub - lb;
                            println!("    pair ({},{}) e({},{}) ref={:.3} bounds=[{:.3},{:.3}] range={:.3} viol={:.3}",
                                ii, jj, mol.atoms[ii].element, mol.atoms[jj].element,
                                rdkit_d, lb, ub, range, viol);
                        }
                    }
                }
                if viol_count == 0 {
                    println!("    No ref violations > 0.03 — oracle failure is from torsion variance, not bounds");
                }
            }
        }
        count += 1;

        if count % 1000 == 0 {
            println!(
                "Processed {} molecules... Current Avg Error: {:.3} Å",
                count,
                total_rmsd / count as f32
            );
        }
    }

    let final_avg_rmsd = total_rmsd / count as f32;
    let oracle_avg = oracle_total / count as f32;
    println!("=== TEST COMPLETE ===");
    println!("Successfully processed {} molecules.", count);
    println!(
        "Selection(BV): Avg {:.3} Max {:.3} Above0.5 {}",
        final_avg_rmsd, select_max, select_above
    );
    println!(
        "Oracle:        Avg {:.3} Max {:.3} Above0.5 {}",
        oracle_avg, oracle_max, oracle_above
    );
    println!("--- Strategy comparison (Above 0.5 / Avg RMSD) ---");
    for si in 0..strat_names.len() {
        println!(
            "  {:20} : Above0.5 {:3}  Avg {:.3}",
            strat_names[si],
            strat_above[si],
            strat_total[si] / count as f32
        );
    }
}

// Function to calculate RMSD between two sets of coordinates.
// Future implementation when we generate 3D coordinates.
pub fn calculate_rmsd(coords1: &[(f32, f32, f32)], coords2: &[(f32, f32, f32)]) -> f32 {
    assert_eq!(coords1.len(), coords2.len());
    let mut sum_sq_diff = 0.0;

    for (p1, p2) in coords1.iter().zip(coords2.iter()) {
        let dx = p1.0 - p2.0;
        let dy = p1.1 - p2.1;
        let dz = p1.2 - p2.2;
        sum_sq_diff += dx * dx + dy * dy + dz * dz;
    }

    (sum_sq_diff / coords1.len() as f32).sqrt()
}

#[test]
fn test_rmsd_calculation_dummy() {
    let raw1 = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)];
    let raw2 = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.1)];
    let rmsd = calculate_rmsd(&raw1, &raw2);
    assert!(rmsd < 0.1);
}
