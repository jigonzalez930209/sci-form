//! Embed, batch, parse, info, charges, uff, rmsd command handlers.

use crate::format::{format_sdf, format_xyz};
use std::io::{self, BufRead, BufReader, Write};
use std::time::Instant;

/// Molecules per streaming batch (bounded memory for large inputs).
const BATCH_CHUNK: usize = 512;

pub fn cmd_embed(smiles: &str, seed: u64, format: &str) {
    let result = sci_form::embed(smiles, seed);
    match format {
        "xyz" => print!("{}", format_xyz(&result)),
        "sdf" => print!("{}", format_sdf(&result)),
        _ => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
    }
    if result.error.is_some() {
        std::process::exit(1);
    }
}

/// Batch-process SMILES in fixed-size chunks, streaming results as they are
/// produced instead of materializing the whole input and output first.
pub fn cmd_batch(
    input: Option<String>,
    output: Option<String>,
    seed: u64,
    threads: usize,
    format: &str,
) {
    let source_label = input.as_deref().unwrap_or("stdin").to_string();
    let reader: Box<dyn BufRead> = match input {
        Some(path) => {
            let file = std::fs::File::open(&path).unwrap_or_else(|e| {
                eprintln!("Error reading {}: {}", path, e);
                std::process::exit(1);
            });
            Box::new(BufReader::new(file))
        }
        None => Box::new(io::stdin().lock()),
    };

    #[cfg(not(feature = "parallel"))]
    if threads > 1 {
        eprintln!(
            "Warning: --threads {} requested but binary was compiled without 'parallel' feature. \
             Running single-threaded.",
            threads
        );
    }

    eprintln!(
        "Processing molecules from {} ({} threads)...",
        source_label,
        if threads == 0 {
            "auto".to_string()
        } else {
            threads.to_string()
        }
    );

    let mut writer: Box<dyn Write> = match output {
        Some(path) => Box::new(std::fs::File::create(&path).unwrap_or_else(|e| {
            eprintln!("Error creating {}: {}", path, e);
            std::process::exit(1);
        })),
        None => Box::new(io::stdout().lock()),
    };

    let config = sci_form::ConformerConfig {
        seed,
        num_threads: threads,
    };
    let start = Instant::now();
    let json = !matches!(format, "xyz" | "sdf");
    let mut first_json = true;
    if json {
        writer.write_all(b"[").unwrap();
    }

    let mut total = 0usize;
    let mut ok = 0usize;
    let mut fail = 0usize;
    let mut chunk: Vec<String> = Vec::with_capacity(BATCH_CHUNK);

    let flush_chunk = |chunk: &mut Vec<String>,
                       writer: &mut dyn Write,
                       total: &mut usize,
                       ok: &mut usize,
                       fail: &mut usize,
                       first_json: &mut bool| {
        if chunk.is_empty() {
            return;
        }
        let smiles_refs: Vec<&str> = chunk.iter().map(|s| s.as_str()).collect();
        let results = sci_form::embed_batch(&smiles_refs, &config);
        for r in &results {
            if r.error.is_none() {
                *ok += 1;
            } else {
                *fail += 1;
            }
        }
        *total += results.len();
        match format {
            "xyz" => {
                for r in &results {
                    write!(writer, "{}", format_xyz(r)).unwrap();
                }
            }
            "sdf" => {
                for r in &results {
                    write!(writer, "{}", format_sdf(r)).unwrap();
                }
            }
            _ => {
                for r in &results {
                    if !*first_json {
                        writer.write_all(b",").unwrap();
                    }
                    *first_json = false;
                    write!(writer, "{}", serde_json::to_string(r).unwrap()).unwrap();
                }
            }
        }
        writer.flush().unwrap();
        chunk.clear();
    };

    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                eprintln!("Error reading input: {}", e);
                break;
            }
        };
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        chunk.push(trimmed.to_string());
        if chunk.len() >= BATCH_CHUNK {
            flush_chunk(
                &mut chunk,
                &mut *writer,
                &mut total,
                &mut ok,
                &mut fail,
                &mut first_json,
            );
        }
    }
    flush_chunk(
        &mut chunk,
        &mut *writer,
        &mut total,
        &mut ok,
        &mut fail,
        &mut first_json,
    );
    if json {
        writer.write_all(b"]\n").unwrap();
        writer.flush().unwrap();
    }

    let elapsed = start.elapsed();
    let avg_ms = if total > 0 {
        elapsed.as_secs_f64() * 1000.0 / total as f64
    } else {
        0.0
    };
    eprintln!(
        "Done: {}/{} OK, {} failed, {:.1}s total ({:.1} ms/mol)",
        ok,
        total,
        fail,
        elapsed.as_secs_f64(),
        avg_ms
    );
}

pub fn cmd_parse(smiles: &str) {
    match sci_form::parse(smiles) {
        Ok(mol) => {
            let n = mol.graph.node_count();
            let nb = mol.graph.edge_count();
            println!("Atoms: {}, Bonds: {}", n, nb);
            for i in 0..n {
                let idx = sci_form::graph::NodeIndex::new(i);
                let atom = &mol.graph[idx];
                println!(
                    "  {:3}: elem={:2} hyb={:?} charge={} h={}",
                    i, atom.element, atom.hybridization, atom.formal_charge, atom.explicit_h
                );
            }
        }
        Err(e) => {
            eprintln!("Parse error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_info() {
    println!("{}", sci_form::version());
    println!("Features: ETKDG distance geometry, CSD torsion patterns, SMARTS matching");
    println!("  EHT, PM3, xTB, HF-3c, ANI, Gasteiger charges, SASA, IR, NMR, UV-Vis");
    println!("Formats: JSON, XYZ, SDF");
    println!("Bindings: CLI, Python (PyO3), TypeScript (WASM)");
}

pub fn cmd_charges(smiles: &str) {
    match sci_form::compute_charges(smiles) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("Charges error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_uff(smiles: &str, coords: &str) {
    let flat: Vec<f64> = serde_json::from_str(coords).unwrap_or_else(|e| {
        eprintln!("Bad coords JSON: {}", e);
        std::process::exit(1);
    });
    match sci_form::compute_uff_energy(smiles, &flat) {
        Ok(energy) => println!("{{\"energy\":{:.6},\"unit\":\"kcal/mol\"}}", energy),
        Err(e) => {
            eprintln!("UFF error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_rmsd(coords: &str, reference: &str) {
    let c: Vec<f64> = serde_json::from_str(coords).unwrap_or_else(|e| {
        eprintln!("Bad coords JSON: {}", e);
        std::process::exit(1);
    });
    let r: Vec<f64> = serde_json::from_str(reference).unwrap_or_else(|e| {
        eprintln!("Bad reference JSON: {}", e);
        std::process::exit(1);
    });
    let result = sci_form::alignment::align_coordinates(&c, &r);
    println!(
        "{{\"rmsd\":{:.6},\"rotation\":{:?}}}",
        result.rmsd, result.rotation
    );
}
