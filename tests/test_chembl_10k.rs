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
//! ChemBL 10K Benchmark — intrinsic geometry quality test for large molecules.
//!
//! Tests the 10,000 largest molecules from ChemBL 36 (110-349 heavy atoms).
//! Validates:
//! - Embed success rate
//! - Bond lengths (within expected covalent radii ranges)
//! - Bond angles (no degenerate angles < 30°)
//! - No steric clashes (non-bonded atoms closer than 0.5 Å)
//! - Performance (ms/molecule)
//!
//! Run:  CHEMBL_LIMIT=100 cargo test --release --test test_chembl_10k -- --nocapture
//! Full: cargo test --release --test test_chembl_10k -- --nocapture

use rayon::prelude::*;
use std::fs;
use std::sync::atomic::{AtomicU32, Ordering};
use std::time::Instant;

struct ChemblResult {
    smiles: String,
    chembl_id: String,
    heavy_atoms: usize,
    total_atoms: usize,
    time_ms: f64,
    embed_ok: bool,
    bad_bonds: usize,
    bad_angles: usize,
    steric_clashes: usize,
    error: Option<String>,
}

fn validate_chembl_molecule(smiles: &str, chembl_id: &str, heavy_atoms: usize) -> ChemblResult {
    let start = Instant::now();

    let mol = match sci_form::graph::Molecule::from_smiles(smiles) {
        Ok(m) => m,
        Err(e) => {
            return ChemblResult {
                smiles: smiles.to_string(),
                chembl_id: chembl_id.to_string(),
                heavy_atoms,
                total_atoms: 0,
                time_ms: start.elapsed().as_secs_f64() * 1000.0,
                embed_ok: false,
                bad_bonds: 0,
                bad_angles: 0,
                steric_clashes: 0,
                error: Some(format!("Parse error: {}", e)),
            };
        }
    };

    let total_atoms = mol.graph.node_count();

    let result = sci_form::conformer::generate_3d_conformer(&mol, 42);
    let time_ms = start.elapsed().as_secs_f64() * 1000.0;

    match result {
        Ok(coords) => {
            let n = coords.nrows();
            let mut bad_bonds = 0usize;
            let mut bad_angles = 0usize;
            let mut steric_clashes = 0usize;

            // Bond length validation
            use petgraph::visit::EdgeRef;
            for edge in mol.graph.edge_references() {
                let (a, b) = (edge.source().index(), edge.target().index());
                if a >= n || b >= n {
                    continue;
                }
                let dx = coords[(a, 0)] - coords[(b, 0)];
                let dy = coords[(a, 1)] - coords[(b, 1)];
                let dz = coords[(a, 2)] - coords[(b, 2)];
                let d = (dx * dx + dy * dy + dz * dz).sqrt();
                if !(0.7..=2.5).contains(&d) {
                    bad_bonds += 1;
                }
            }

            // Angle validation (no degenerate angles < 30°)
            for idx in 0..n {
                let node = petgraph::graph::NodeIndex::new(idx);
                let neighbors: Vec<usize> = mol.graph.neighbors(node).map(|n| n.index()).collect();
                if neighbors.len() < 2 {
                    continue;
                }
                for i in 0..neighbors.len() {
                    for j in (i + 1)..neighbors.len() {
                        let a = neighbors[i];
                        let b = neighbors[j];
                        if a >= n || b >= n {
                            continue;
                        }
                        let vax = coords[(a, 0)] - coords[(idx, 0)];
                        let vay = coords[(a, 1)] - coords[(idx, 1)];
                        let vaz = coords[(a, 2)] - coords[(idx, 2)];
                        let vbx = coords[(b, 0)] - coords[(idx, 0)];
                        let vby = coords[(b, 1)] - coords[(idx, 1)];
                        let vbz = coords[(b, 2)] - coords[(idx, 2)];
                        let dot = vax * vbx + vay * vby + vaz * vbz;
                        let la = (vax * vax + vay * vay + vaz * vaz).sqrt();
                        let lb = (vbx * vbx + vby * vby + vbz * vbz).sqrt();
                        if la < 1e-6 || lb < 1e-6 {
                            continue;
                        }
                        let cos_ang = (dot / (la * lb)).clamp(-1.0, 1.0);
                        let ang_deg = cos_ang.acos().to_degrees();
                        if ang_deg < 30.0 {
                            bad_angles += 1;
                        }
                    }
                }
            }

            // Steric clash validation (non-bonded atoms < 0.5 Å)
            let mut bonded = std::collections::HashSet::new();
            for edge in mol.graph.edge_references() {
                let (a, b) = (edge.source().index(), edge.target().index());
                bonded.insert((a, b));
                bonded.insert((b, a));
            }
            for a in 0..n {
                for b in (a + 1)..n {
                    if bonded.contains(&(a, b)) {
                        continue;
                    }
                    let dx = coords[(a, 0)] - coords[(b, 0)];
                    let dy = coords[(a, 1)] - coords[(b, 1)];
                    let dz = coords[(a, 2)] - coords[(b, 2)];
                    let d = (dx * dx + dy * dy + dz * dz).sqrt();
                    if d < 0.5 {
                        steric_clashes += 1;
                    }
                }
            }

            ChemblResult {
                smiles: smiles.to_string(),
                chembl_id: chembl_id.to_string(),
                heavy_atoms,
                total_atoms,
                time_ms,
                embed_ok: true,
                bad_bonds,
                bad_angles,
                steric_clashes,
                error: None,
            }
        }
        Err(e) => ChemblResult {
            smiles: smiles.to_string(),
            chembl_id: chembl_id.to_string(),
            heavy_atoms,
            total_atoms,
            time_ms,
            embed_ok: false,
            bad_bonds: 0,
            bad_angles: 0,
            steric_clashes: 0,
            error: Some(e),
        },
    }
}

#[test]
fn test_chembl_10k() {
    let chembl_file = std::env::var("CHEMBL_FILE")
        .unwrap_or_else(|_| "data/chembl_10k_practical.smi".to_string());
    let data = fs::read_to_string(&chembl_file).unwrap_or_else(|_| {
        panic!(
            "File not found: {}. Run scripts/extract_chembl_practical.py first",
            chembl_file
        )
    });

    let molecules: Vec<(&str, &str, usize)> = data
        .lines()
        .filter_map(|line| {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 {
                let heavy: usize = parts[2].trim().parse().ok()?;
                Some((parts[0], parts[1], heavy))
            } else {
                None
            }
        })
        .collect();

    let limit: usize = std::env::var("CHEMBL_LIMIT")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(molecules.len());
    let molecules = &molecules[..limit.min(molecules.len())];

    let total = molecules.len();
    let ncpus = rayon::current_num_threads();
    eprintln!("\n=== CHEMBL 10K BENCHMARK ===");
    eprintln!(
        "Molecules: {}, Threads: {}, Heavy atoms: {}-{}",
        total,
        ncpus,
        molecules.last().map(|m| m.2).unwrap_or(0),
        molecules.first().map(|m| m.2).unwrap_or(0)
    );
    eprintln!();

    let done = AtomicU32::new(0);
    let start = Instant::now();

    let results: Vec<ChemblResult> = molecules
        .par_iter()
        .map(|(smiles, id, heavy)| {
            let result = validate_chembl_molecule(smiles, id, *heavy);
            let d = done.fetch_add(1, Ordering::Relaxed) + 1;
            if d.is_multiple_of(100) || d == total as u32 {
                let elapsed = start.elapsed().as_secs_f64();
                let rate = d as f64 / elapsed;
                eprintln!("  Progress: {}/{} ({:.1} mol/s)", d, total, rate);
            }
            result
        })
        .collect();

    let elapsed = start.elapsed().as_secs_f64();

    // === Summary Statistics ===
    let parse_fail = results
        .iter()
        .filter(|r| r.error.as_ref().is_some_and(|e| e.starts_with("Parse")))
        .count();
    let embed_ok = results.iter().filter(|r| r.embed_ok).count();
    let embed_fail = results
        .iter()
        .filter(|r| !r.embed_ok && !r.error.as_ref().is_some_and(|e| e.starts_with("Parse")))
        .count();
    let perfect = results
        .iter()
        .filter(|r| r.embed_ok && r.bad_bonds == 0 && r.bad_angles == 0 && r.steric_clashes == 0)
        .count();
    let with_bad_bonds = results.iter().filter(|r| r.bad_bonds > 0).count();
    let with_bad_angles = results.iter().filter(|r| r.bad_angles > 0).count();
    let with_clashes = results.iter().filter(|r| r.steric_clashes > 0).count();
    let total_bad_bonds: usize = results.iter().map(|r| r.bad_bonds).sum();
    let total_bad_angles: usize = results.iter().map(|r| r.bad_angles).sum();
    let total_clashes: usize = results.iter().map(|r| r.steric_clashes).sum();

    let embedded_results: Vec<&ChemblResult> = results.iter().filter(|r| r.embed_ok).collect();
    let avg_time: f64 = if embedded_results.is_empty() {
        0.0
    } else {
        embedded_results.iter().map(|r| r.time_ms).sum::<f64>() / embedded_results.len() as f64
    };
    let p95_time = {
        let mut times: Vec<f64> = embedded_results.iter().map(|r| r.time_ms).collect();
        times.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if times.is_empty() {
            0.0
        } else {
            times[(times.len() as f64 * 0.95) as usize]
        }
    };
    let max_time = embedded_results
        .iter()
        .map(|r| r.time_ms)
        .fold(0.0f64, f64::max);

    eprintln!(
        "\n=== CHEMBL BENCHMARK RESULTS ({} molecules, {:.1}s) ===",
        total, elapsed
    );
    eprintln!(
        "Parse failures:    {}/{} ({:.2}%)",
        parse_fail,
        total,
        100.0 * parse_fail as f64 / total as f64
    );
    eprintln!(
        "Embed success:     {}/{} ({:.2}%)",
        embed_ok,
        total - parse_fail,
        100.0 * embed_ok as f64 / (total - parse_fail).max(1) as f64
    );
    eprintln!(
        "Embed failures:    {}/{} ({:.2}%)",
        embed_fail,
        total - parse_fail,
        100.0 * embed_fail as f64 / (total - parse_fail).max(1) as f64
    );
    eprintln!(
        "Perfect geometry:  {}/{} ({:.2}%)",
        perfect,
        embed_ok,
        100.0 * perfect as f64 / embed_ok.max(1) as f64
    );
    eprintln!(
        "Bad bonds:         {} issues in {} molecules",
        total_bad_bonds, with_bad_bonds
    );
    eprintln!(
        "Bad angles:        {} issues in {} molecules",
        total_bad_angles, with_bad_angles
    );
    eprintln!(
        "Steric clashes:    {} issues in {} molecules",
        total_clashes, with_clashes
    );
    eprintln!();
    eprintln!("Performance (embedded only):");
    eprintln!(
        "  Avg: {:.1}ms, P95: {:.1}ms, Max: {:.1}ms",
        avg_time, p95_time, max_time
    );
    eprintln!(
        "  Rate: {:.1} mol/s (wall, {} threads)",
        total as f64 / elapsed,
        ncpus
    );

    // === By size bucket ===
    eprintln!("\n=== BY HEAVY ATOM COUNT ===");
    let buckets = [
        (30, 39),
        (40, 49),
        (50, 59),
        (60, 69),
        (70, 79),
        (80, 89),
        (90, 100),
    ];
    for (lo, hi) in buckets {
        let in_bucket: Vec<&ChemblResult> = results
            .iter()
            .filter(|r| r.heavy_atoms >= lo && r.heavy_atoms <= hi)
            .collect();
        if in_bucket.is_empty() {
            continue;
        }
        let n = in_bucket.len();
        let ok = in_bucket.iter().filter(|r| r.embed_ok).count();
        let perf = in_bucket
            .iter()
            .filter(|r| {
                r.embed_ok && r.bad_bonds == 0 && r.bad_angles == 0 && r.steric_clashes == 0
            })
            .count();
        let avg_t: f64 = if ok > 0 {
            in_bucket
                .iter()
                .filter(|r| r.embed_ok)
                .map(|r| r.time_ms)
                .sum::<f64>()
                / ok as f64
        } else {
            0.0
        };
        eprintln!(
            "  {}-{}: n={:>5}, embed={:.1}%, perfect={:.1}%, avg={:.0}ms",
            lo,
            hi,
            n,
            100.0 * ok as f64 / n as f64,
            100.0 * perf as f64 / ok.max(1) as f64,
            avg_t
        );
    }

    // === Top 10 worst failures ===
    let mut failures: Vec<&ChemblResult> = results.iter().filter(|r| !r.embed_ok).collect();
    failures.sort_by_key(|r| r.heavy_atoms);
    if !failures.is_empty() {
        eprintln!("\n=== SAMPLE FAILURES (smallest first) ===");
        for f in failures.iter().take(10) {
            eprintln!(
                "  {} ({}): heavy={}, total={}, err={}",
                f.chembl_id,
                &f.smiles[..f.smiles.len().min(60)],
                f.heavy_atoms,
                f.total_atoms,
                f.error.as_deref().unwrap_or("?")
            );
        }
    }

    // === Top 10 worst geometry ===
    let mut bad_geom: Vec<&ChemblResult> = results
        .iter()
        .filter(|r| r.embed_ok && (r.bad_bonds > 0 || r.bad_angles > 0 || r.steric_clashes > 0))
        .collect();
    bad_geom.sort_by(|a, b| {
        (b.bad_bonds + b.bad_angles + b.steric_clashes)
            .cmp(&(a.bad_bonds + a.bad_angles + a.steric_clashes))
    });
    if !bad_geom.is_empty() {
        eprintln!("\n=== WORST GEOMETRY (most issues) ===");
        for g in bad_geom.iter().take(10) {
            eprintln!(
                "  {} (heavy={}): bonds={}, angles={}, clashes={}",
                g.chembl_id, g.heavy_atoms, g.bad_bonds, g.bad_angles, g.steric_clashes
            );
        }
    }

    eprintln!();
}
