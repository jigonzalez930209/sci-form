#![allow(unused_imports, unused_variables, dead_code, clippy::unnecessary_cast, clippy::needless_range_loop, clippy::manual_repeat_n, clippy::manual_str_repeat, clippy::manual_is_multiple_of, clippy::redundant_field_names, clippy::useless_vec, clippy::single_range_in_vec_init)]
//! GDB-20 50k parallel RMSD test.
//!
//! Loads RDKit reference coordinates for ~50k GDB-20 molecules,
//! generates conformers in parallel using all CPU cores (1 molecule per thread),
//! and reports pairwise-distance RMSD vs RDKit reference.
//!
//! Run:  cargo test --release --test test_gdb20_rmsd test_gdb20_parallel -- --nocapture
//! Limit: GDB20_LIMIT=1000 cargo test ...

use rayon::prelude::*;
use serde::Deserialize;
use std::fs;
use std::sync::atomic::{AtomicU32, Ordering};
use std::time::Instant;

#[derive(Deserialize)]
struct RefAtom {
    element: u8,
    x: f32,
    y: f32,
    z: f32,
    formal_charge: i8,
    hybridization: String,
}

#[derive(Deserialize)]
struct RefBond {
    start: usize,
    end: usize,
    order: String,
}

#[derive(Deserialize)]
struct RefTorsion {
    atoms: Vec<usize>,
    v: Vec<f64>,
    signs: Vec<i32>,
}

#[derive(Deserialize)]
#[allow(dead_code)]
struct RefMolecule {
    smiles: String,
    atoms: Vec<RefAtom>,
    bonds: Vec<RefBond>,
    torsions: Vec<RefTorsion>,
}

/// Per-molecule result from parallel processing
#[allow(dead_code)]
struct MolResult {
    smiles: String,
    n_atoms: usize,
    rmsd: Option<f32>, // None = embed failure
    time_ms: f64,
    error: Option<String>,
}

/// Compute pairwise-distance RMSD between our coords and reference atoms.
fn pairwise_rmsd(coords: &nalgebra::DMatrix<f32>, ref_atoms: &[RefAtom]) -> f32 {
    let n = ref_atoms.len();
    let mut sq_sum = 0.0f64;
    let mut npairs = 0u64;
    for a in 0..n {
        for b in (a + 1)..n {
            let dr = ((ref_atoms[a].x - ref_atoms[b].x).powi(2)
                + (ref_atoms[a].y - ref_atoms[b].y).powi(2)
                + (ref_atoms[a].z - ref_atoms[b].z).powi(2))
            .sqrt() as f64;
            let du = ((coords[(a, 0)] - coords[(b, 0)]).powi(2)
                + (coords[(a, 1)] - coords[(b, 1)]).powi(2)
                + (coords[(a, 2)] - coords[(b, 2)]).powi(2))
            .sqrt() as f64;
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

/// Build a Molecule from reference data (using oracle atoms/bonds directly).
fn build_mol_from_ref(ref_mol: &RefMolecule) -> sci_form::graph::Molecule {
    let mut mol = sci_form::graph::Molecule::new(&ref_mol.smiles);
    let mut node_indices = Vec::with_capacity(ref_mol.atoms.len());

    for atom in &ref_mol.atoms {
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
            position: nalgebra::Vector3::zeros(),
            charge: 0.0,
            formal_charge: atom.formal_charge,
            hybridization,
            chiral_tag: sci_form::graph::ChiralType::Unspecified,
            explicit_h: if atom.element == 1 || atom.element == 0 {
                1
            } else {
                0
            },
        };
        node_indices.push(mol.add_atom(new_atom));
    }

    for bond in &ref_mol.bonds {
        let order = match bond.order.as_str() {
            "DOUBLE" => sci_form::graph::BondOrder::Double,
            "TRIPLE" => sci_form::graph::BondOrder::Triple,
            "AROMATIC" => sci_form::graph::BondOrder::Aromatic,
            _ => sci_form::graph::BondOrder::Single,
        };
        mol.add_bond(
            node_indices[bond.start],
            node_indices[bond.end],
            sci_form::graph::Bond {
                order,
                stereo: sci_form::graph::BondStereo::None,
            },
        );
    }
    mol
}

/// Build CSD torsion contribs from reference data.
fn build_csd_torsions(
    ref_torsions: &[RefTorsion],
) -> Vec<sci_form::forcefield::etkdg_3d::M6TorsionContrib> {
    ref_torsions
        .iter()
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

/// Process a single molecule: parse → embed → compute RMSD.
fn process_molecule(ref_mol: &RefMolecule) -> MolResult {
    let start = Instant::now();

    let mol = build_mol_from_ref(ref_mol);
    let csd_torsions = build_csd_torsions(&ref_mol.torsions);

    let num_seeds: usize = std::env::var("BEST_OF_K")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1);

    let result = if num_seeds > 1 {
        sci_form::conformer::generate_3d_conformer_best_of_k(&mol, 42, &csd_torsions, num_seeds)
    } else {
        sci_form::conformer::generate_3d_conformer_with_torsions(&mol, 42, &csd_torsions)
    };

    let time_ms = start.elapsed().as_secs_f64() * 1000.0;

    match result {
        Ok(coords) => {
            let rmsd = pairwise_rmsd(&coords, &ref_mol.atoms);
            MolResult {
                smiles: ref_mol.smiles.clone(),
                n_atoms: ref_mol.atoms.len(),
                rmsd: Some(rmsd),
                time_ms,
                error: None,
            }
        }
        Err(e) => MolResult {
            smiles: ref_mol.smiles.clone(),
            n_atoms: ref_mol.atoms.len(),
            rmsd: None,
            time_ms,
            error: Some(e),
        },
    }
}

#[test]
fn test_gdb20_parallel() {
    // Load reference data (skip gracefully when the large fixture is absent)
    let fixture = "tests/fixtures/gdb20_reference.json";
    if !std::path::Path::new(fixture).exists() {
        eprintln!("SKIP {fixture}: run scripts/generate_gdb20_reference.py to generate it");
        return;
    }
    let ref_data = fs::read_to_string(fixture).expect("Failed to read gdb20_reference.json");
    let mut ref_mols: Vec<RefMolecule> =
        serde_json::from_str(&ref_data).expect("Invalid gdb20_reference.json");

    // Sort by heaviest molecules first (most atoms) for a representative sample
    ref_mols.sort_by(|a, b| b.atoms.len().cmp(&a.atoms.len()));

    // Optional limit — default 300 heaviest; use GDB20_LIMIT=0 for full run
    let limit: usize = std::env::var("GDB20_LIMIT")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(300);
    let limit = if limit == 0 { ref_mols.len() } else { limit };
    let ref_mols = &ref_mols[..limit.min(ref_mols.len())];

    let total = ref_mols.len();
    let ncpus = rayon::current_num_threads();
    eprintln!("\n=== GDB-20 PARALLEL RMSD TEST ===");
    eprintln!("Molecules: {}, Threads: {}\n", total, ncpus);

    // Atomic progress counter
    let done = AtomicU32::new(0);
    let start = Instant::now();

    // Process all molecules in parallel
    let results: Vec<MolResult> = ref_mols
        .par_iter()
        .map(|ref_mol| {
            let result = process_molecule(ref_mol);
            let d = done.fetch_add(1, Ordering::Relaxed) + 1;
            if d.is_multiple_of(2000) {
                let elapsed = start.elapsed().as_secs_f64();
                let rate = d as f64 / elapsed;
                let eta = (total as f64 - d as f64) / rate;
                eprintln!(
                    "  {}/{} done ({:.0} mol/s, ETA {:.0}s)",
                    d, total, rate, eta
                );
            }
            result
        })
        .collect();

    let elapsed = start.elapsed();

    // Collect statistics
    let mut embed_ok = 0u32;
    let mut embed_fail = 0u32;
    let mut above_05 = 0u32;
    let mut above_1 = 0u32;
    let mut total_rmsd = 0.0f64;
    let mut max_rmsd = 0.0f32;
    let mut max_rmsd_smi = String::new();
    let mut rmsd_buckets = [0u32; 10];
    let mut failures: Vec<&MolResult> = Vec::new();
    let mut high_rmsd: Vec<&MolResult> = Vec::new();

    for r in &results {
        match r.rmsd {
            Some(rmsd) => {
                embed_ok += 1;
                total_rmsd += rmsd as f64;
                if rmsd > max_rmsd {
                    max_rmsd = rmsd;
                    max_rmsd_smi = r.smiles.clone();
                }
                if rmsd > 0.5 {
                    above_05 += 1;
                    high_rmsd.push(r);
                }
                if rmsd > 1.0 {
                    above_1 += 1;
                }
                let bucket = (rmsd * 10.0).min(9.0) as usize;
                rmsd_buckets[bucket] += 1;
            }
            None => {
                embed_fail += 1;
                failures.push(r);
            }
        }
    }

    let avg_rmsd = if embed_ok > 0 {
        total_rmsd / embed_ok as f64
    } else {
        0.0
    };

    // Sort high RMSD by descending RMSD
    let mut high_rmsd_sorted: Vec<(&str, f32)> = high_rmsd
        .iter()
        .map(|r| (r.smiles.as_str(), r.rmsd.unwrap()))
        .collect();
    high_rmsd_sorted.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    // Print first 10 failing molecules for debugging
    println!("\n=== FIRST 10 FAILING MOLECULES (RMSD > 0.5 Å) ===");
    for (idx, (smi, rmsd)) in high_rmsd_sorted.iter().take(10).enumerate() {
        println!("  #{}: RMSD={:.4} SMILES={}", idx + 1, rmsd, smi);
    }

    // Print results
    println!("\n=== GDB-20 50K RESULTS ===");
    println!("Total molecules: {}", total);
    println!("Embed OK: {}, Embed FAIL: {}", embed_ok, embed_fail);
    println!("Avg RMSD: {:.4} Å", avg_rmsd);
    println!("Max RMSD: {:.3} Å ({})", max_rmsd, max_rmsd_smi);
    println!(
        "Above 0.5 Å: {} ({:.2}%)",
        above_05,
        above_05 as f64 / embed_ok.max(1) as f64 * 100.0
    );
    println!(
        "Above 1.0 Å: {} ({:.2}%)",
        above_1,
        above_1 as f64 / embed_ok.max(1) as f64 * 100.0
    );

    println!("\nRMSD Distribution:");
    let labels = [
        "0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8",
        "0.8-0.9", "0.9+",
    ];
    for (bucket, label) in rmsd_buckets.iter().zip(labels.iter()) {
        let pct = *bucket as f64 / embed_ok.max(1) as f64 * 100.0;
        let bar: String = std::iter::repeat_n('#', (pct * 0.5) as usize).collect();
        println!("  {:8}: {:6} ({:5.2}%) {}", label, bucket, pct, bar);
    }

    println!(
        "\nTime: {:.1}s ({:.1} ms/mol, {} threads)",
        elapsed.as_secs_f64(),
        elapsed.as_secs_f64() * 1000.0 / total as f64,
        ncpus
    );

    if !failures.is_empty() {
        println!("\n--- Embed Failures (first 30) ---");
        for r in failures.iter().take(30) {
            println!(
                "  {} (n={}) → {}",
                r.smiles,
                r.n_atoms,
                r.error.as_deref().unwrap_or("?")
            );
        }
    }

    if !high_rmsd_sorted.is_empty() {
        println!("\n--- Highest RMSD (first 50) ---");
        for (smi, rmsd) in high_rmsd_sorted.iter().take(50) {
            println!("  {:.3} Å  {}", rmsd, smi);
        }
    }

    // === Structural feature analysis ===
    println!("\n=== STRUCTURAL FEATURE ANALYSIS ===");

    // Track feature correlations with failure
    struct FeatureStats {
        total: u32,
        above_05: u32,
        sum_rmsd: f64,
    }
    impl FeatureStats {
        fn new() -> Self {
            Self {
                total: 0,
                above_05: 0,
                sum_rmsd: 0.0,
            }
        }
        fn add(&mut self, rmsd: f32) {
            self.total += 1;
            self.sum_rmsd += rmsd as f64;
            if rmsd > 0.5 {
                self.above_05 += 1;
            }
        }
        fn fail_pct(&self) -> f64 {
            if self.total == 0 {
                0.0
            } else {
                self.above_05 as f64 / self.total as f64 * 100.0
            }
        }
        fn avg_rmsd(&self) -> f64 {
            if self.total == 0 {
                0.0
            } else {
                self.sum_rmsd / self.total as f64
            }
        }
    }

    let mut has_triple = FeatureStats::new();
    let mut no_triple = FeatureStats::new();
    let mut has_double = FeatureStats::new();
    let mut no_double = FeatureStats::new();
    let mut has_aromatic = FeatureStats::new();
    let mut no_aromatic = FeatureStats::new();
    // By heavy atom count ranges
    let mut heavy_small = FeatureStats::new(); // <=10
    let mut heavy_medium = FeatureStats::new(); // 11-15
    let mut heavy_large = FeatureStats::new(); // 16-20
                                               // By CSD torsion count
    let mut csd_0 = FeatureStats::new();
    let mut csd_1_3 = FeatureStats::new();
    let mut csd_4plus = FeatureStats::new();
    // By total atom count
    let mut atoms_lt20 = FeatureStats::new();
    let mut atoms_20_35 = FeatureStats::new();
    let mut atoms_36plus = FeatureStats::new();
    // SMILES patterns
    let mut has_ring_smi = FeatureStats::new(); // contains digits (ring closure)
    let mut no_ring_smi = FeatureStats::new();
    // Combined: triple + ring
    let mut triple_with_ring = FeatureStats::new();
    let mut triple_no_ring = FeatureStats::new();
    let mut no_triple_with_ring = FeatureStats::new();
    let mut no_triple_no_ring = FeatureStats::new();

    for (i, r) in results.iter().enumerate() {
        let rmsd = match r.rmsd {
            Some(v) => v,
            None => continue,
        };
        let ref_mol = &ref_mols[i];

        // Check bond types
        let has_triple_bond = ref_mol.bonds.iter().any(|b| b.order == "TRIPLE");
        let has_double_bond = ref_mol.bonds.iter().any(|b| b.order == "DOUBLE");
        let has_aromatic_bond = ref_mol.bonds.iter().any(|b| b.order == "AROMATIC");
        let has_ring_closure = ref_mol.smiles.chars().any(|c| c.is_ascii_digit());

        if has_triple_bond {
            has_triple.add(rmsd);
        } else {
            no_triple.add(rmsd);
        }
        if has_double_bond {
            has_double.add(rmsd);
        } else {
            no_double.add(rmsd);
        }
        if has_aromatic_bond {
            has_aromatic.add(rmsd);
        } else {
            no_aromatic.add(rmsd);
        }
        if has_ring_closure {
            has_ring_smi.add(rmsd);
        } else {
            no_ring_smi.add(rmsd);
        }

        // Combined
        match (has_triple_bond, has_ring_closure) {
            (true, true) => triple_with_ring.add(rmsd),
            (true, false) => triple_no_ring.add(rmsd),
            (false, true) => no_triple_with_ring.add(rmsd),
            (false, false) => no_triple_no_ring.add(rmsd),
        }

        // Heavy atom count
        let n_heavy = ref_mol.atoms.iter().filter(|a| a.element != 1).count();
        if n_heavy <= 10 {
            heavy_small.add(rmsd);
        } else if n_heavy <= 15 {
            heavy_medium.add(rmsd);
        } else {
            heavy_large.add(rmsd);
        }

        // Total atom count
        let n_total = ref_mol.atoms.len();
        if n_total < 20 {
            atoms_lt20.add(rmsd);
        } else if n_total <= 35 {
            atoms_20_35.add(rmsd);
        } else {
            atoms_36plus.add(rmsd);
        }

        // CSD torsion count
        let n_csd = ref_mol.torsions.len();
        if n_csd == 0 {
            csd_0.add(rmsd);
        } else if n_csd <= 3 {
            csd_1_3.add(rmsd);
        } else {
            csd_4plus.add(rmsd);
        }
    }

    println!("\nBond type correlation:");
    println!(
        "  Has triple bond:  n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        has_triple.total,
        has_triple.fail_pct(),
        has_triple.avg_rmsd()
    );
    println!(
        "  No triple bond:   n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        no_triple.total,
        no_triple.fail_pct(),
        no_triple.avg_rmsd()
    );
    println!(
        "  Has double bond:  n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        has_double.total,
        has_double.fail_pct(),
        has_double.avg_rmsd()
    );
    println!(
        "  No double bond:   n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        no_double.total,
        no_double.fail_pct(),
        no_double.avg_rmsd()
    );
    println!(
        "  Has aromatic:     n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        has_aromatic.total,
        has_aromatic.fail_pct(),
        has_aromatic.avg_rmsd()
    );
    println!(
        "  No aromatic:      n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        no_aromatic.total,
        no_aromatic.fail_pct(),
        no_aromatic.avg_rmsd()
    );

    println!("\nRing presence (SMILES digits):");
    println!(
        "  Has rings:        n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        has_ring_smi.total,
        has_ring_smi.fail_pct(),
        has_ring_smi.avg_rmsd()
    );
    println!(
        "  No rings:         n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        no_ring_smi.total,
        no_ring_smi.fail_pct(),
        no_ring_smi.avg_rmsd()
    );

    println!("\nTriple × Ring cross-tab:");
    println!(
        "  Triple + Ring:    n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        triple_with_ring.total,
        triple_with_ring.fail_pct(),
        triple_with_ring.avg_rmsd()
    );
    println!(
        "  Triple, no Ring:  n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        triple_no_ring.total,
        triple_no_ring.fail_pct(),
        triple_no_ring.avg_rmsd()
    );
    println!(
        "  No Triple + Ring: n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        no_triple_with_ring.total,
        no_triple_with_ring.fail_pct(),
        no_triple_with_ring.avg_rmsd()
    );
    println!(
        "  No Triple, no Ring: n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        no_triple_no_ring.total,
        no_triple_no_ring.fail_pct(),
        no_triple_no_ring.avg_rmsd()
    );

    println!("\nHeavy atom count:");
    println!(
        "  <=10 atoms:       n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        heavy_small.total,
        heavy_small.fail_pct(),
        heavy_small.avg_rmsd()
    );
    println!(
        "  11-15 atoms:      n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        heavy_medium.total,
        heavy_medium.fail_pct(),
        heavy_medium.avg_rmsd()
    );
    println!(
        "  16-20 atoms:      n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        heavy_large.total,
        heavy_large.fail_pct(),
        heavy_large.avg_rmsd()
    );

    println!("\nTotal atom count:");
    println!(
        "  <20 atoms:        n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        atoms_lt20.total,
        atoms_lt20.fail_pct(),
        atoms_lt20.avg_rmsd()
    );
    println!(
        "  20-35 atoms:      n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        atoms_20_35.total,
        atoms_20_35.fail_pct(),
        atoms_20_35.avg_rmsd()
    );
    println!(
        "  36+ atoms:        n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        atoms_36plus.total,
        atoms_36plus.fail_pct(),
        atoms_36plus.avg_rmsd()
    );

    println!("\nCSD torsion count:");
    println!(
        "  0 torsions:       n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        csd_0.total,
        csd_0.fail_pct(),
        csd_0.avg_rmsd()
    );
    println!(
        "  1-3 torsions:     n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        csd_1_3.total,
        csd_1_3.fail_pct(),
        csd_1_3.avg_rmsd()
    );
    println!(
        "  4+ torsions:      n={:5}, fail={:5.2}%, avg_rmsd={:.4}",
        csd_4plus.total,
        csd_4plus.fail_pct(),
        csd_4plus.avg_rmsd()
    );

    // Quality gate
    let fail_pct = above_05 as f64 / embed_ok.max(1) as f64 * 100.0;
    println!("\nQuality gate: {:.2}% above 0.5 Å (target: 0%)", fail_pct);
}
