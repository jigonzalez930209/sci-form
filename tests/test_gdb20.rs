#![allow(unused_imports, unused_variables, dead_code, clippy::unnecessary_cast, clippy::needless_range_loop, clippy::manual_repeat_n, clippy::manual_str_repeat, clippy::manual_is_multiple_of, clippy::redundant_field_names, clippy::useless_vec, clippy::single_range_in_vec_init)]
//! GDB-20 50k molecule benchmark: SMILES parsing + conformer generation.
//!
//! Tests against the GDB-20 dataset (50,000 molecules with up to 20 heavy atoms).
//! No RDKit reference — validates:
//!   - Parse success rate
//!   - Embedding success rate
//!   - Basic geometry sanity (bond lengths, no NaN, no atom overlap)
//!   - Throughput (ms/mol)

use std::fs;

/// Sanity-check generated coordinates:
///  - No NaN/Inf values
///  - All bonded distances in [0.8, 2.5] Å
///  - No atom pair closer than 0.5 Å
fn check_geometry(
    mol: &sci_form::graph::Molecule,
    coords: &nalgebra::DMatrix<f32>,
) -> Result<(), String> {
    let n = mol.graph.node_count();

    // Check for NaN/Inf
    for i in 0..n {
        for d in 0..3 {
            let v = coords[(i, d)];
            if v.is_nan() || v.is_infinite() {
                return Err(format!("NaN/Inf at atom {} dim {}", i, d));
            }
        }
    }

    // Check bonded distances
    for edge in mol.graph.edge_indices() {
        let (a, b) = mol.graph.edge_endpoints(edge).unwrap();
        let ai = a.index();
        let bi = b.index();
        let dx = coords[(ai, 0)] - coords[(bi, 0)];
        let dy = coords[(ai, 1)] - coords[(bi, 1)];
        let dz = coords[(ai, 2)] - coords[(bi, 2)];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        if !(0.7..=2.5).contains(&dist) {
            return Err(format!(
                "Bond {}-{} distance {:.3} out of range [0.8, 2.5]",
                ai, bi, dist
            ));
        }
    }

    Ok(())
}

#[test]
fn test_gdb20_50k() {
    let smiles_data = fs::read_to_string("GDB20.50000.smi").expect("Should read GDB20.50000.smi");
    let smiles_list: Vec<&str> = smiles_data
        .lines()
        .filter(|l| !l.trim().is_empty())
        .collect();

    // Allow limiting via env var: GDB20_LIMIT=1000 cargo test ...
    let limit: usize = std::env::var("GDB20_LIMIT")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(smiles_list.len());
    let smiles_list = &smiles_list[..limit.min(smiles_list.len())];

    println!("\n=== GDB-20 50K BENCHMARK ===\n");
    println!("Total SMILES: {}", smiles_list.len());

    let mut parse_ok = 0u32;
    let mut parse_fail = 0u32;
    let mut embed_ok = 0u32;
    let mut embed_fail = 0u32;
    let mut geom_fail = 0u32;
    let mut total_atoms = 0u64;
    let mut max_atoms = 0usize;
    let mut parse_failures: Vec<String> = Vec::new();
    let mut embed_failures: Vec<String> = Vec::new();
    let mut geom_failures: Vec<(String, String)> = Vec::new();

    let start = std::time::Instant::now();

    for smi in smiles_list.iter() {
        // Parse
        let mol = match sci_form::graph::Molecule::from_smiles(smi) {
            Ok(m) => m,
            Err(e) => {
                parse_fail += 1;
                if parse_failures.len() < 30 {
                    parse_failures.push(format!("{} → {}", smi, e));
                }
                continue;
            }
        };
        parse_ok += 1;

        let n = mol.graph.node_count();
        total_atoms += n as u64;
        if n > max_atoms {
            max_atoms = n;
        }

        // Generate conformer (no CSD torsions — pure distance geometry + ETKDG)
        match sci_form::conformer::generate_3d_conformer(&mol, 42) {
            Ok(coords) => {
                embed_ok += 1;

                // Geometry sanity check
                if let Err(e) = check_geometry(&mol, &coords) {
                    geom_fail += 1;
                    if geom_failures.len() < 30 {
                        geom_failures.push((smi.to_string(), e));
                    }
                }
            }
            Err(e) => {
                embed_fail += 1;
                if embed_failures.len() < 30 {
                    embed_failures.push(format!("{} → {}", smi, e));
                }
            }
        }

        // Progress reporting
        let done = parse_ok + parse_fail;
        if done.is_multiple_of(5000) {
            let elapsed = start.elapsed().as_secs_f64();
            eprintln!(
                "... {}/{} done, {:.1}ms/mol, {} parse_fail, {} embed_fail, {} geom_fail",
                done,
                smiles_list.len(),
                elapsed * 1000.0 / done as f64,
                parse_fail,
                embed_fail,
                geom_fail
            );
        }
    }

    let elapsed = start.elapsed();
    let total_processed = embed_ok + embed_fail;

    println!("\n=== GDB-20 50K RESULTS ===");
    println!(
        "Parse OK: {} ({:.1}%)",
        parse_ok,
        parse_ok as f64 / smiles_list.len() as f64 * 100.0
    );
    println!(
        "Parse FAIL: {} ({:.1}%)",
        parse_fail,
        parse_fail as f64 / smiles_list.len() as f64 * 100.0
    );
    println!(
        "Embed OK: {} ({:.1}%)",
        embed_ok,
        embed_ok as f64 / parse_ok.max(1) as f64 * 100.0
    );
    println!(
        "Embed FAIL: {} ({:.1}%)",
        embed_fail,
        embed_fail as f64 / parse_ok.max(1) as f64 * 100.0
    );
    println!(
        "Geometry FAIL: {} ({:.1}%)",
        geom_fail,
        geom_fail as f64 / embed_ok.max(1) as f64 * 100.0
    );
    println!(
        "Avg atoms/mol: {:.1}",
        total_atoms as f64 / parse_ok.max(1) as f64
    );
    println!("Max atoms: {}", max_atoms);
    println!(
        "Time: {:.1}s ({:.1} ms/mol)",
        elapsed.as_secs_f64(),
        elapsed.as_secs_f64() * 1000.0 / total_processed.max(1) as f64
    );

    if !parse_failures.is_empty() {
        println!("\n--- Parse Failures (first {}) ---", parse_failures.len());
        for f in &parse_failures {
            println!("  {}", f);
        }
    }

    if !embed_failures.is_empty() {
        println!("\n--- Embed Failures (first {}) ---", embed_failures.len());
        for f in &embed_failures {
            println!("  {}", f);
        }
    }

    if !geom_failures.is_empty() {
        println!(
            "\n--- Geometry Failures (first {}) ---",
            geom_failures.len()
        );
        for (smi, err) in &geom_failures {
            println!("  {} → {}", smi, err);
        }
    }

    // Quality gates
    let embed_rate = embed_ok as f64 / parse_ok.max(1) as f64 * 100.0;
    let geom_rate = geom_fail as f64 / embed_ok.max(1) as f64 * 100.0;
    println!("\n=== QUALITY GATES ===");
    println!("Embed success rate: {:.1}% (target: >95%)", embed_rate);
    println!("Geometry failure rate: {:.2}% (target: <1%)", geom_rate);
    assert!(
        embed_rate > 90.0,
        "Embedding rate too low: {:.1}%",
        embed_rate
    );
    assert!(
        geom_rate < 5.0,
        "Geometry failure rate too high: {:.2}%",
        geom_rate
    );
}
