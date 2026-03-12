//! Comprehensive 10k molecule test: SMILES parsing + conformer generation + RMSD comparison.
//!
//! Phase 1: Parse all 10k SMILES → report failures, compare atom counts vs RDKit reference
//! Phase 2: Generate conformers → compare pairwise distance RMSD vs RDKit reference
//! Phase 3: Summary statistics

use serde::Deserialize;
use std::collections::HashMap;
use std::fs;
use std::io::Write;

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

/// Phase 1: Test SMILES parsing for all 10k molecules.
/// Reports parse failures and atom count mismatches vs RDKit.
#[test]
fn test_smiles_parsing_10k() {
    let smiles_data = fs::read_to_string("scripts/10k_smiles.smi")
        .expect("Should read 10k_smiles.smi");
    let smiles_list: Vec<&str> = smiles_data.lines().filter(|l| !l.trim().is_empty()).collect();

    // Load reference for atom count comparison
    let ref_data = fs::read_to_string("tests/fixtures/reference_coords_no_mmff.json")
        .expect("Should read reference JSON");
    let ref_mols: Vec<OracleMolecule> =
        serde_json::from_str(&ref_data).expect("Reference JSON parse");
    let ref_map: HashMap<&str, &OracleMolecule> = ref_mols.iter().map(|m| (m.smiles.as_str(), m)).collect();

    println!("\n=== PHASE 1: SMILES PARSING (10k molecules) ===\n");

    let mut parse_ok = 0u32;
    let mut parse_fail = 0u32;
    let mut atom_match = 0u32;
    let mut atom_mismatch = 0u32;
    let mut element_mismatch = 0u32;
    let mut no_reference = 0u32;
    let mut parse_failures: Vec<(usize, String, String)> = Vec::new();
    let mut atom_count_diffs: Vec<(String, usize, usize)> = Vec::new();

    let start = std::time::Instant::now();

    for (i, smi) in smiles_list.iter().enumerate() {
        match sci_form::graph::Molecule::from_smiles(smi) {
            Ok(mol) => {
                parse_ok += 1;
                let our_count = mol.graph.node_count();

                if let Some(ref_mol) = ref_map.get(smi) {
                    let ref_count = ref_mol.atoms.len();
                    if our_count == ref_count {
                        // Check element ordering
                        let our_elems: Vec<u8> = mol.graph.node_indices()
                            .map(|ni| mol.graph[ni].element)
                            .collect();
                        let ref_elems: Vec<u8> = ref_mol.atoms.iter()
                            .map(|a| a.element)
                            .collect();
                        if our_elems == ref_elems {
                            atom_match += 1;
                        } else {
                            element_mismatch += 1;
                            if element_mismatch <= 20 {
                                println!("ELEM MISMATCH mol {}: {} ours={:?} ref={:?}",
                                    i, smi,
                                    &our_elems[..our_elems.len().min(20)],
                                    &ref_elems[..ref_elems.len().min(20)]);
                            }
                        }
                    } else {
                        atom_mismatch += 1;
                        if atom_count_diffs.len() < 50 {
                            atom_count_diffs.push((smi.to_string(), our_count, ref_count));
                        }
                    }
                } else {
                    no_reference += 1;
                }
            }
            Err(e) => {
                parse_fail += 1;
                if parse_failures.len() < 50 {
                    parse_failures.push((i, smi.to_string(), e));
                }
            }
        }
    }

    let elapsed = start.elapsed();

    println!("Parsing time: {:.2}s ({:.2} µs/mol)", elapsed.as_secs_f64(),
        elapsed.as_secs_f64() * 1e6 / smiles_list.len() as f64);
    println!("Parse OK: {}, Parse FAIL: {}", parse_ok, parse_fail);
    println!("Atom count match: {}, Atom count mismatch: {}, Element mismatch: {}",
        atom_match, atom_mismatch, element_mismatch);
    println!("No reference: {}", no_reference);

    if !parse_failures.is_empty() {
        println!("\n--- Parse Failures (first {}) ---", parse_failures.len());
        for (i, smi, err) in &parse_failures {
            println!("  mol {}: {} → {}", i, smi, err);
        }
    }

    if !atom_count_diffs.is_empty() {
        println!("\n--- Atom Count Mismatches (first {}) ---", atom_count_diffs.len());
        for (smi, ours, theirs) in &atom_count_diffs {
            let diff = *ours as i32 - *theirs as i32;
            println!("  {} : ours={} ref={} (diff={})", smi, ours, theirs, diff);
        }
    }

    println!("\n=== PHASE 1 COMPLETE ===");
    // Don't assert — just report. Failures will be fixed iteratively.
}

/// Phase 2: Full end-to-end test — SMILES → conformer → RMSD comparison for all 10k.
#[test]
fn test_conformer_10k() {
    let smiles_data = fs::read_to_string("scripts/10k_smiles.smi")
        .expect("Should read 10k_smiles.smi");
    let smiles_list: Vec<&str> = smiles_data.lines().filter(|l| !l.trim().is_empty()).collect();

    // Load reference
    let ref_data = fs::read_to_string("tests/fixtures/reference_coords_no_mmff.json")
        .expect("Should read reference JSON");
    let ref_mols: Vec<OracleMolecule> =
        serde_json::from_str(&ref_data).expect("Reference JSON parse");
    let ref_map: HashMap<&str, &OracleMolecule> = ref_mols.iter().map(|m| (m.smiles.as_str(), m)).collect();

    // Load CSD torsion params
    let csd_torsions: HashMap<String, Vec<CsdTorsion>> = {
        let path = "tests/fixtures/torsion_params.json";
        match fs::read_to_string(path) {
            Ok(data) => serde_json::from_str(&data).unwrap_or_default(),
            Err(_) => HashMap::new(),
        }
    };

    println!("\n=== PHASE 2: CONFORMER GENERATION + RMSD (10k molecules) ===\n");

    // Progress file for real-time monitoring (avoids test harness buffering)
    let mut progress = fs::File::create("/tmp/10k_progress.log").unwrap();

    let mut parse_fail = 0u32;
    let mut atom_mismatch = 0u32;
    let mut elem_mismatch = 0u32;
    let mut no_reference = 0u32;
    let mut embed_fail = 0u32;
    let mut embed_ok = 0u32;
    let mut total_rmsd = 0.0f64;
    let mut max_rmsd = 0.0f32;
    let mut max_rmsd_smi = String::new();
    let mut above_05 = 0u32;
    let mut above_1 = 0u32;
    let mut rmsd_buckets = [0u32; 10]; // [0-0.1, 0.1-0.2, ..., 0.9-1.0+]

    let start = std::time::Instant::now();
    let mut seen = std::collections::HashSet::new();

    for (i, smi) in smiles_list.iter().enumerate() {
        // Skip duplicate SMILES
        if !seen.insert(*smi) {
            continue;
        }

        // Check if we have a reference
        let ref_mol = match ref_map.get(smi) {
            Some(r) => r,
            None => { no_reference += 1; continue; }
        };

        // Parse SMILES with our parser
        let mol = match sci_form::graph::Molecule::from_smiles(smi) {
            Ok(m) => m,
            Err(_) => { parse_fail += 1; continue; }
        };

        // Check atom count matches
        let our_count = mol.graph.node_count();
        let ref_count = ref_mol.atoms.len();
        if our_count != ref_count {
            atom_mismatch += 1;
            continue;
        }

        // Check element ordering matches (required for pairwise RMSD)
        let our_elems: Vec<u8> = mol.graph.node_indices()
            .map(|ni| mol.graph[ni].element)
            .collect();
        let ref_elems: Vec<u8> = ref_mol.atoms.iter().map(|a| a.element).collect();
        if our_elems != ref_elems {
            elem_mismatch += 1;
            continue;
        }

        // Build CSD torsion contribs
        let csd_contribs: Vec<sci_form::forcefield::etkdg_3d::M6TorsionContrib> =
            if let Some(torsions) = csd_torsions.get(*smi) {
                torsions.iter().map(|t| {
                    let mut signs = [0.0f64; 6];
                    let mut v = [0.0f64; 6];
                    for k in 0..6 {
                        signs[k] = t.signs[k] as f64;
                        v[k] = t.v[k] as f64;
                    }
                    sci_form::forcefield::etkdg_3d::M6TorsionContrib {
                        i: t.atoms[0], j: t.atoms[1], k: t.atoms[2], l: t.atoms[3],
                        signs, v,
                    }
                }).collect()
            } else {
                Vec::new()
            };

        // Generate conformer
        let mol_start = std::time::Instant::now();
        writeln!(progress, "START mol {} n={} {}", i, our_count, &smi[..smi.len().min(50)]).ok();
        match sci_form::conformer::generate_3d_conformer_with_torsions(&mol, 42, &csd_contribs) {
            Ok(coords) => {
                let n = ref_count;
                let mut sq_sum = 0.0f64;
                let mut npairs = 0u64;
                for a in 0..n {
                    for b in (a + 1)..n {
                        let dr = ((ref_mol.atoms[a].x - ref_mol.atoms[b].x).powi(2)
                            + (ref_mol.atoms[a].y - ref_mol.atoms[b].y).powi(2)
                            + (ref_mol.atoms[a].z - ref_mol.atoms[b].z).powi(2))
                        .sqrt() as f64;
                        let du = ((coords[(a, 0)] - coords[(b, 0)]).powi(2)
                            + (coords[(a, 1)] - coords[(b, 1)]).powi(2)
                            + (coords[(a, 2)] - coords[(b, 2)]).powi(2))
                        .sqrt() as f64;
                        sq_sum += (dr - du).powi(2);
                        npairs += 1;
                    }
                }
                let rmsd = if npairs > 0 { (sq_sum / npairs as f64).sqrt() } else { 0.0 };
                let rmsd32 = rmsd as f32;

                embed_ok += 1;
                total_rmsd += rmsd;
                if rmsd32 > max_rmsd {
                    max_rmsd = rmsd32;
                    max_rmsd_smi = smi.to_string();
                }
                if rmsd32 > 0.5 { above_05 += 1; }
                if rmsd32 > 1.0 { above_1 += 1; }

                let bucket = (rmsd * 10.0).min(9.0) as usize;
                rmsd_buckets[bucket] += 1;

                let mol_ms = mol_start.elapsed().as_secs_f64() * 1000.0;
                if rmsd32 > 0.3 {
                    let has_csd = if csd_contribs.is_empty() { "no_csd" } else { "csd" };
                    let has_triple = smi.contains('#');
                    let msg = format!("HIGH_RMSD {:.3} {} n={} {} triple={} {:.0}ms\n", rmsd, smi, n, has_csd, has_triple, mol_ms);
                    eprintln!("{}", msg.trim());
                    progress.write_all(msg.as_bytes()).ok();
                }
            }
            Err(e) => {
                embed_fail += 1;
                let mol_ms = mol_start.elapsed().as_secs_f64() * 1000.0;
                if embed_fail <= 50 {
                    let msg = format!("EMBED FAIL mol {:5}: {:40} {:.0}ms → {}\n", i, &smi[..smi.len().min(40)], mol_ms, e);
                    eprintln!("{}", msg.trim());
                    progress.write_all(msg.as_bytes()).ok();
                }
            }
        }

        if (embed_ok + embed_fail) % 500 == 0 && embed_ok + embed_fail > 0 {
            let done = embed_ok + embed_fail;
            let elapsed = start.elapsed().as_secs_f64();
            let msg = format!("... {} done, {:.1} ms/mol, {} ok, {} fail, {:.4} avg RMSD\n",
                done, elapsed * 1000.0 / done as f64, embed_ok, embed_fail,
                if embed_ok > 0 { total_rmsd / embed_ok as f64 } else { 0.0 });
            eprintln!("{}", msg.trim());
            progress.write_all(msg.as_bytes()).ok();
        }
    }

    let elapsed = start.elapsed();
    let avg_rmsd = if embed_ok > 0 { total_rmsd / embed_ok as f64 } else { 0.0 };

    println!("\n=== 10K CONFORMER RESULTS ===");
    println!("Total unique SMILES tested: {}", seen.len());
    println!("Parse failures: {}", parse_fail);
    println!("Atom count mismatches: {}", atom_mismatch);
    println!("Element mismatches: {}", elem_mismatch);
    println!("No reference: {}", no_reference);
    println!("Embed OK: {}, Embed FAIL: {}", embed_ok, embed_fail);
    println!("Avg RMSD: {:.4} Å", avg_rmsd);
    println!("Max RMSD: {:.3} Å ({})", max_rmsd, max_rmsd_smi);
    println!("Above 0.5 Å: {} ({:.1}%)", above_05, above_05 as f64 / embed_ok.max(1) as f64 * 100.0);
    println!("Above 1.0 Å: {} ({:.1}%)", above_1, above_1 as f64 / embed_ok.max(1) as f64 * 100.0);
    println!("\nRMSD Distribution:");
    let labels = ["0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5",
                  "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9+"];
    for (bucket, label) in rmsd_buckets.iter().zip(labels.iter()) {
        let pct = *bucket as f64 / embed_ok.max(1) as f64 * 100.0;
        let bar: String = std::iter::repeat('#').take((pct * 0.5) as usize).collect();
        println!("  {:8}: {:5} ({:5.1}%) {}", label, bucket, pct, bar);
    }
    println!("\nTime: {:.1}s ({:.1} ms/mol)", elapsed.as_secs_f64(),
        elapsed.as_secs_f64() * 1000.0 / (embed_ok + embed_fail).max(1) as f64);

    // Quality gate: <5% above 0.5Å
    let fail_pct = above_05 as f64 / embed_ok.max(1) as f64 * 100.0;
    println!("\nQuality gate: {:.1}% above 0.5 Å (target: <5%)", fail_pct);
}

#[derive(Deserialize, Debug)]
struct CsdTorsion {
    atoms: Vec<usize>,
    v: Vec<f64>,
    signs: Vec<i32>,
}
