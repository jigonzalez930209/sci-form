//! Embed, batch, parse, info, charges, uff, rmsd command handlers.

use crate::format::{format_sdf, format_xyz};
use std::io::{self, BufRead, Write};
use std::time::Instant;

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

pub fn cmd_batch(
    input: Option<String>,
    output: Option<String>,
    seed: u64,
    threads: usize,
    format: &str,
) {
    let lines: Vec<String> = match input {
        Some(path) => {
            let file = std::fs::File::open(&path).unwrap_or_else(|e| {
                eprintln!("Error reading {}: {}", path, e);
                std::process::exit(1);
            });
            io::BufReader::new(file)
                .lines()
                .filter_map(|l| l.ok())
                .filter(|l| !l.trim().is_empty())
                .map(|l| l.trim().to_string())
                .collect()
        }
        None => {
            let stdin = io::stdin();
            stdin
                .lock()
                .lines()
                .filter_map(|l| l.ok())
                .filter(|l| !l.trim().is_empty())
                .map(|l| l.trim().to_string())
                .collect()
        }
    };
    let total = lines.len();
    let start = Instant::now();

    #[cfg(not(feature = "parallel"))]
    if threads > 1 {
        eprintln!(
            "Warning: --threads {} requested but binary was compiled without 'parallel' feature. \
             Running single-threaded.",
            threads
        );
    }

    eprintln!(
        "Processing {} molecules ({} threads)...",
        total,
        if threads == 0 {
            "auto".to_string()
        } else {
            threads.to_string()
        }
    );

    let smiles_refs: Vec<&str> = lines.iter().map(|s| s.as_str()).collect();
    let config = sci_form::ConformerConfig {
        seed,
        num_threads: threads,
    };
    let results = sci_form::embed_batch(&smiles_refs, &config);
    let elapsed = start.elapsed();

    let mut writer: Box<dyn Write> = match output {
        Some(path) => Box::new(std::fs::File::create(&path).unwrap_or_else(|e| {
            eprintln!("Error creating {}: {}", path, e);
            std::process::exit(1);
        })),
        None => Box::new(io::stdout()),
    };
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
            writeln!(writer, "{}", serde_json::to_string(&results).unwrap()).unwrap();
        }
    }
    let ok = results.iter().filter(|r| r.error.is_none()).count();
    let fail = results.iter().filter(|r| r.error.is_some()).count();
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
