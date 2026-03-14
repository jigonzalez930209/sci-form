#![allow(unused_imports, unused_variables, dead_code, clippy::unnecessary_cast, clippy::needless_range_loop, clippy::manual_repeat_n, clippy::manual_str_repeat, clippy::manual_is_multiple_of, clippy::redundant_field_names, clippy::useless_vec, clippy::single_range_in_vec_init)]
//! Ensemble RMSD comparison test.
//!
//! Compares our conformers against an ensemble of RDKit conformers (21 seeds).
//! For each molecule, computes the MINIMUM pairwise-distance RMSD across all
//! RDKit seeds. This is a fair comparison because flexible molecules can have
//! multiple valid conformations with very different geometries.
//!
//! Run: GDB20_LIMIT=1000 cargo test --release --test test_ensemble_rmsd -- --nocapture

use rayon::prelude::*;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs;
use std::sync::atomic::{AtomicU32, Ordering};
use std::time::Instant;

#[derive(Deserialize)]
struct EnsembleAtom {
    element: u8,
    hybridization: String,
    formal_charge: i8,
}

#[derive(Deserialize)]
struct EnsembleBond {
    start: usize,
    end: usize,
    order: String,
}

#[derive(Deserialize)]
struct EnsembleTorsion {
    atoms: Vec<usize>,
    v: Vec<f64>,
    signs: Vec<i32>,
}

#[derive(Deserialize)]
#[allow(dead_code)]
struct EnsembleMolecule {
    smiles: String,
    atoms: Vec<EnsembleAtom>,
    bonds: Vec<EnsembleBond>,
    torsions: Vec<EnsembleTorsion>,
    conformers: HashMap<String, Vec<Vec<f64>>>,
}

/// Original reference molecule (has correct torsion data)
#[derive(Deserialize)]
#[allow(dead_code)]
struct RefMolecule {
    smiles: String,
    atoms: Vec<EnsembleAtom>,
    bonds: Vec<EnsembleBond>,
    torsions: Vec<EnsembleTorsion>,
}

#[allow(dead_code)]
struct EnsembleResult {
    smiles: String,
    n_atoms: usize,
    min_rmsd: Option<f32>,       // min RMSD across all RDKit seeds (all atoms)
    min_rmsd_heavy: Option<f32>, // min RMSD across all RDKit seeds (heavy atoms only)
    rmsd_seed42: Option<f32>,    // RMSD vs seed=42 specifically
    best_seed: Option<String>,   // which RDKit seed gives min RMSD
    n_rdkit_seeds: usize,
    time_ms: f64,
    error: Option<String>,
}

fn pairwise_rmsd(coords: &nalgebra::DMatrix<f32>, ref_coords: &[Vec<f64>]) -> f32 {
    let n = ref_coords.len();
    let mut sq_sum = 0.0f64;
    let mut npairs = 0u64;
    for a in 0..n {
        for b in (a + 1)..n {
            let dr = ((ref_coords[a][0] - ref_coords[b][0]).powi(2)
                + (ref_coords[a][1] - ref_coords[b][1]).powi(2)
                + (ref_coords[a][2] - ref_coords[b][2]).powi(2))
            .sqrt();
            let du = ((coords[(a, 0)] as f64 - coords[(b, 0)] as f64).powi(2)
                + (coords[(a, 1)] as f64 - coords[(b, 1)] as f64).powi(2)
                + (coords[(a, 2)] as f64 - coords[(b, 2)] as f64).powi(2))
            .sqrt();
            sq_sum += (dr - du).powi(2);
            npairs += 1;
        }
    }
    if npairs > 0 {
        (sq_sum / npairs as f64).sqrt() as f32
    } else {
        0.0
    }
}

/// Heavy-atom-only pairwise distance RMSD (skip hydrogens)
fn pairwise_rmsd_heavy(
    coords: &nalgebra::DMatrix<f32>,
    ref_coords: &[Vec<f64>],
    atoms: &[EnsembleAtom],
) -> f32 {
    let heavy_indices: Vec<usize> = atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| a.element != 1)
        .map(|(i, _)| i)
        .collect();
    let n = heavy_indices.len();
    let mut sq_sum = 0.0f64;
    let mut npairs = 0u64;
    for ai in 0..n {
        for bi in (ai + 1)..n {
            let a = heavy_indices[ai];
            let b = heavy_indices[bi];
            let dr = ((ref_coords[a][0] - ref_coords[b][0]).powi(2)
                + (ref_coords[a][1] - ref_coords[b][1]).powi(2)
                + (ref_coords[a][2] - ref_coords[b][2]).powi(2))
            .sqrt();
            let du = ((coords[(a, 0)] as f64 - coords[(b, 0)] as f64).powi(2)
                + (coords[(a, 1)] as f64 - coords[(b, 1)] as f64).powi(2)
                + (coords[(a, 2)] as f64 - coords[(b, 2)] as f64).powi(2))
            .sqrt();
            sq_sum += (dr - du).powi(2);
            npairs += 1;
        }
    }
    if npairs > 0 {
        (sq_sum / npairs as f64).sqrt() as f32
    } else {
        0.0
    }
}

fn build_mol(m: &EnsembleMolecule) -> sci_form::graph::Molecule {
    let mut mol = sci_form::graph::Molecule::new(&m.smiles);
    let mut nodes = Vec::with_capacity(m.atoms.len());
    for a in &m.atoms {
        let hyb = match a.hybridization.as_str() {
            "SP" => sci_form::graph::Hybridization::SP,
            "SP2" => sci_form::graph::Hybridization::SP2,
            "SP3" => sci_form::graph::Hybridization::SP3,
            "SP3D" => sci_form::graph::Hybridization::SP3D,
            "SP3D2" => sci_form::graph::Hybridization::SP3D2,
            _ => sci_form::graph::Hybridization::Unknown,
        };
        nodes.push(mol.add_atom(sci_form::graph::Atom {
            element: a.element,
            position: nalgebra::Vector3::zeros(),
            charge: 0.0,
            formal_charge: a.formal_charge,
            hybridization: hyb,
            chiral_tag: sci_form::graph::ChiralType::Unspecified,
            explicit_h: if a.element == 1 || a.element == 0 {
                1
            } else {
                0
            },
        }));
    }
    for b in &m.bonds {
        let order = match b.order.as_str() {
            "DOUBLE" => sci_form::graph::BondOrder::Double,
            "TRIPLE" => sci_form::graph::BondOrder::Triple,
            "AROMATIC" => sci_form::graph::BondOrder::Aromatic,
            _ => sci_form::graph::BondOrder::Single,
        };
        mol.add_bond(
            nodes[b.start],
            nodes[b.end],
            sci_form::graph::Bond {
                order,
                stereo: sci_form::graph::BondStereo::None,
            },
        );
    }
    mol
}

fn build_torsions(t: &[EnsembleTorsion]) -> Vec<sci_form::forcefield::etkdg_3d::M6TorsionContrib> {
    t.iter()
        .filter_map(|t| {
            if t.atoms.len() < 4 || t.v.len() < 6 || t.signs.len() < 6 {
                return None;
            }
            let mut signs = [0.0f64; 6];
            let mut v = [0.0f64; 6];
            for k in 0..6 {
                signs[k] = t.signs[k] as f64;
                v[k] = t.v[k];
            }
            Some(sci_form::forcefield::etkdg_3d::M6TorsionContrib {
                i: t.atoms[0],
                j: t.atoms[1],
                k: t.atoms[2],
                l: t.atoms[3],
                signs,
                v,
            })
        })
        .collect()
}

fn process(
    m: &EnsembleMolecule,
    torsion_data: &[EnsembleTorsion],
    our_seeds: &[u64],
    use_own_torsions: bool,
) -> EnsembleResult {
    let start = Instant::now();
    let mol = build_mol(m);
    let torsions = if use_own_torsions {
        sci_form::smarts::match_experimental_torsions(&mol)
    } else {
        build_torsions(torsion_data)
    };

    let mut global_min_rmsd = f32::MAX;
    let mut global_min_rmsd_heavy = f32::MAX;
    let mut global_best_seed = String::new();
    let mut rmsd_42_42 = None;
    let mut any_ok = false;

    for &our_seed in our_seeds {
        let result =
            sci_form::conformer::generate_3d_conformer_with_torsions(&mol, our_seed, &torsions);
        if let Ok(coords) = result {
            any_ok = true;
            for (rdkit_seed, ref_coords) in &m.conformers {
                if ref_coords.len() != coords.nrows() {
                    continue;
                }
                let r = pairwise_rmsd(&coords, ref_coords);
                let r_heavy = pairwise_rmsd_heavy(&coords, ref_coords, &m.atoms);
                if our_seed == 42 && rdkit_seed == "42" {
                    rmsd_42_42 = Some(r);
                }
                if r < global_min_rmsd {
                    global_min_rmsd = r;
                    global_best_seed = format!("ours={}/rdkit={}", our_seed, rdkit_seed);
                }
                if r_heavy < global_min_rmsd_heavy {
                    global_min_rmsd_heavy = r_heavy;
                }
            }
        }
    }
    let time_ms = start.elapsed().as_secs_f64() * 1000.0;

    if any_ok {
        EnsembleResult {
            smiles: m.smiles.clone(),
            n_atoms: m.atoms.len(),
            min_rmsd: if global_min_rmsd < f32::MAX {
                Some(global_min_rmsd)
            } else {
                None
            },
            min_rmsd_heavy: if global_min_rmsd_heavy < f32::MAX {
                Some(global_min_rmsd_heavy)
            } else {
                None
            },
            rmsd_seed42: rmsd_42_42,
            best_seed: Some(global_best_seed),
            n_rdkit_seeds: m.conformers.len(),
            time_ms,
            error: None,
        }
    } else {
        EnsembleResult {
            smiles: m.smiles.clone(),
            n_atoms: m.atoms.len(),
            min_rmsd: None,
            min_rmsd_heavy: None,
            rmsd_seed42: None,
            best_seed: None,
            n_rdkit_seeds: m.conformers.len(),
            time_ms,
            error: Some("all seeds failed".to_string()),
        }
    }
}

#[test]
fn test_ensemble_rmsd() {
    let ensemble_fixture = "tests/fixtures/gdb20_ensemble.json";
    let reference_fixture = "tests/fixtures/gdb20_reference.json";
    for fixture in [ensemble_fixture, reference_fixture] {
        if !std::path::Path::new(fixture).exists() {
            eprintln!("SKIP {fixture}: run scripts/generate_ensemble_reference.py to generate it");
            return;
        }
    }
    let data = fs::read_to_string(ensemble_fixture).expect("Failed to read gdb20_ensemble.json");
    let mut mols: Vec<EnsembleMolecule> =
        serde_json::from_str(&data).expect("Invalid gdb20_ensemble.json");

    // Sort by heaviest molecules first (most atoms) for a representative sample
    mols.sort_by(|a, b| b.atoms.len().cmp(&a.atoms.len()));

    // Load original reference for correct torsion data
    let ref_data =
        fs::read_to_string(reference_fixture).expect("Failed to read gdb20_reference.json");
    let ref_mols: Vec<RefMolecule> =
        serde_json::from_str(&ref_data).expect("Invalid gdb20_reference.json");

    // Build SMILES → torsion index
    let ref_torsion_map: HashMap<&str, &Vec<EnsembleTorsion>> = ref_mols
        .iter()
        .map(|m| (m.smiles.as_str(), &m.torsions))
        .collect();

    // Optional limit — default 300 heaviest; use GDB20_LIMIT=0 for all
    let limit: usize = std::env::var("GDB20_LIMIT")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(300);
    let limit = if limit == 0 { mols.len() } else { limit };
    let mols = &mols[..limit.min(mols.len())];

    // Our seeds: OUR_SEEDS=1 means just seed 42, OUR_SEEDS=21 means seeds 42,0-19
    let num_our_seeds: usize = std::env::var("OUR_SEEDS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1);
    let our_seeds: Vec<u64> = if num_our_seeds <= 1 {
        vec![42]
    } else {
        let mut s = vec![42u64];
        for i in 0..((num_our_seeds - 1) as u64) {
            s.push(i);
        }
        s
    };

    // OWN_TORSIONS=1 uses sci-form's own SMARTS matcher instead of RDKit oracle torsions
    let use_own_torsions = std::env::var("OWN_TORSIONS")
        .ok()
        .map(|s| s == "1" || s == "true")
        .unwrap_or(false);

    let total = mols.len();
    let ncpus = rayon::current_num_threads();
    eprintln!(
        "\n=== ENSEMBLE RMSD TEST ({} molecules, {} RDKit seeds, {} our seeds, torsions={}, {} threads) ===\n",
        total, mols[0].conformers.len(), our_seeds.len(),
        if use_own_torsions { "own" } else { "oracle" }, ncpus
    );

    let done = AtomicU32::new(0);
    let start = Instant::now();

    let results: Vec<EnsembleResult> = mols
        .par_iter()
        .map(|m| {
            let torsion_data = ref_torsion_map
                .get(m.smiles.as_str())
                .map(|t| t.as_slice())
                .unwrap_or(&[]);
            let r = process(m, torsion_data, &our_seeds, use_own_torsions);
            let d = done.fetch_add(1, Ordering::Relaxed) + 1;
            if d.is_multiple_of(500) {
                eprintln!("  Progress: {}/{}", d, total);
            }
            r
        })
        .collect();

    let elapsed = start.elapsed();

    // === Statistics ===
    let embed_ok: Vec<&EnsembleResult> = results.iter().filter(|r| r.min_rmsd.is_some()).collect();
    let embed_fail = results.len() - embed_ok.len();

    // Single-seed comparison (vs seed=42 only, current method)
    let rmsds_42: Vec<f32> = embed_ok.iter().filter_map(|r| r.rmsd_seed42).collect();
    let above_05_seed42 = rmsds_42.iter().filter(|&&r| r > 0.5).count();

    // Ensemble comparison (min RMSD across all seeds)
    let min_rmsds: Vec<f32> = embed_ok.iter().filter_map(|r| r.min_rmsd).collect();
    let above_05_ensemble = min_rmsds.iter().filter(|&&r| r > 0.5).count();
    let above_03_ensemble = min_rmsds.iter().filter(|&&r| r > 0.3).count();
    let above_01_ensemble = min_rmsds.iter().filter(|&&r| r > 0.1).count();

    let avg_min = min_rmsds.iter().map(|&r| r as f64).sum::<f64>() / min_rmsds.len().max(1) as f64;
    let avg_42 = rmsds_42.iter().map(|&r| r as f64).sum::<f64>() / rmsds_42.len().max(1) as f64;

    println!(
        "\n=== ENSEMBLE RMSD RESULTS ({} molecules, {:.1}s) ===",
        total,
        elapsed.as_secs_f64()
    );
    println!("Embed OK: {}, Embed FAIL: {}", embed_ok.len(), embed_fail);
    println!();
    println!("--- Single-seed comparison (vs seed=42) ---");
    println!("  Avg RMSD: {:.4} Å", avg_42);
    println!(
        "  Above 0.5 Å: {} ({:.2}%)",
        above_05_seed42,
        100.0 * above_05_seed42 as f64 / rmsds_42.len().max(1) as f64
    );
    println!();
    println!(
        "--- Ensemble comparison (min RMSD over {} seeds) ---",
        mols[0].conformers.len()
    );
    println!("  Avg min-RMSD: {:.4} Å", avg_min);
    println!(
        "  Above 0.1 Å: {} ({:.2}%)",
        above_01_ensemble,
        100.0 * above_01_ensemble as f64 / min_rmsds.len().max(1) as f64
    );
    println!(
        "  Above 0.3 Å: {} ({:.2}%)",
        above_03_ensemble,
        100.0 * above_03_ensemble as f64 / min_rmsds.len().max(1) as f64
    );
    println!(
        "  Above 0.5 Å: {} ({:.2}%)",
        above_05_ensemble,
        100.0 * above_05_ensemble as f64 / min_rmsds.len().max(1) as f64
    );

    // Heavy-atom ensemble comparison
    let min_rmsds_heavy: Vec<f32> = embed_ok.iter().filter_map(|r| r.min_rmsd_heavy).collect();
    let above_05_heavy = min_rmsds_heavy.iter().filter(|&&r| r > 0.5).count();
    let above_03_heavy = min_rmsds_heavy.iter().filter(|&&r| r > 0.3).count();
    let above_01_heavy = min_rmsds_heavy.iter().filter(|&&r| r > 0.1).count();
    let avg_heavy = min_rmsds_heavy.iter().map(|&r| r as f64).sum::<f64>()
        / min_rmsds_heavy.len().max(1) as f64;

    println!();
    println!("--- Heavy-atom ensemble comparison (min RMSD, heavy atoms only) ---");
    println!("  Avg min-RMSD: {:.4} Å", avg_heavy);
    println!(
        "  Above 0.1 Å: {} ({:.2}%)",
        above_01_heavy,
        100.0 * above_01_heavy as f64 / min_rmsds_heavy.len().max(1) as f64
    );
    println!(
        "  Above 0.3 Å: {} ({:.2}%)",
        above_03_heavy,
        100.0 * above_03_heavy as f64 / min_rmsds_heavy.len().max(1) as f64
    );
    println!(
        "  Above 0.5 Å: {} ({:.2}%)",
        above_05_heavy,
        100.0 * above_05_heavy as f64 / min_rmsds_heavy.len().max(1) as f64
    );

    // Distribution (all atoms)
    println!("\n--- min-RMSD Distribution (all atoms) ---");
    let buckets_def = [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, f32::MAX];
    for w in buckets_def.windows(2) {
        let lo = w[0];
        let hi = w[1];
        let count = min_rmsds.iter().filter(|&&r| r >= lo && r < hi).count();
        let pct = 100.0 * count as f64 / min_rmsds.len().max(1) as f64;
        let bar: String = std::iter::repeat_n('#', (pct * 0.5) as usize).collect();
        let label = if hi == f32::MAX {
            format!("{:.2}+", lo)
        } else {
            format!("{:.2}-{:.2}", lo, hi)
        };
        println!("  {:12}: {:5} ({:5.2}%) {}", label, count, pct, bar);
    }

    // Distribution (heavy atoms only)
    println!("\n--- min-RMSD Distribution (heavy atoms) ---");
    for w in buckets_def.windows(2) {
        let lo = w[0];
        let hi = w[1];
        let count = min_rmsds_heavy
            .iter()
            .filter(|&&r| r >= lo && r < hi)
            .count();
        let pct = 100.0 * count as f64 / min_rmsds_heavy.len().max(1) as f64;
        let bar: String = std::iter::repeat_n('#', (pct * 0.5) as usize).collect();
        let label = if hi == f32::MAX {
            format!("{:.2}+", lo)
        } else {
            format!("{:.2}-{:.2}", lo, hi)
        };
        println!("  {:12}: {:5} ({:5.2}%) {}", label, count, pct, bar);
    }

    // Worst molecules (highest min-RMSD)
    let mut worst: Vec<(&EnsembleResult, f32)> = embed_ok
        .iter()
        .filter_map(|r| r.min_rmsd.map(|mr| (*r, mr)))
        .collect();
    worst.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    if !worst.is_empty() {
        println!("\n--- Worst molecules (highest min-RMSD all-atom, first 20) ---");
        for (r, mr) in worst.iter().take(20) {
            let s42 = r
                .rmsd_seed42
                .map_or("N/A".to_string(), |v| format!("{:.3}", v));
            let heavy = r
                .min_rmsd_heavy
                .map_or("N/A".to_string(), |v| format!("{:.3}", v));
            println!(
                "  all={:.3} heavy={:6}  seed42={:6}  best={}  {}",
                mr,
                heavy,
                s42,
                r.best_seed.as_deref().unwrap_or("?"),
                &r.smiles[..r.smiles.len().min(50)]
            );
        }
    }

    // How much does ensemble help?
    let improved: Vec<(f32, f32)> = embed_ok
        .iter()
        .filter_map(|r| match (r.rmsd_seed42, r.min_rmsd) {
            (Some(s42), Some(min_r)) if s42 > 0.5 => Some((s42, min_r)),
            _ => None,
        })
        .collect();

    if !improved.is_empty() {
        let rescued = improved.iter().filter(|(_, min_r)| *min_r <= 0.5).count();
        println!("\n--- Ensemble rescue rate (seed42 > 0.5 → min-RMSD ≤ 0.5) ---");
        println!(
            "  {}/{} molecules rescued ({:.1}%)",
            rescued,
            improved.len(),
            100.0 * rescued as f64 / improved.len() as f64
        );
        let still_failing: Vec<&(f32, f32)> =
            improved.iter().filter(|(_, min_r)| *min_r > 0.5).collect();
        if !still_failing.is_empty() {
            println!("  Still failing: {} molecules", still_failing.len());
        }
    }

    println!();
}
