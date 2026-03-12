use clap::{Parser, Subcommand};
use std::io::{self, BufRead, Write};
use std::time::Instant;

#[derive(Parser)]
#[command(name = "sci-form", version, about = "3D molecular conformer generation")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Generate 3D conformer for a single SMILES
    Embed {
        /// SMILES string
        smiles: String,
        /// RNG seed
        #[arg(short, long, default_value_t = 42)]
        seed: u64,
        /// Output format: json, xyz, sdf
        #[arg(short, long, default_value = "json")]
        format: String,
    },
    /// Batch-process SMILES from stdin or file (one per line)
    Batch {
        /// Input file (omit for stdin)
        #[arg(short, long)]
        input: Option<String>,
        /// Output file (omit for stdout)
        #[arg(short, long)]
        output: Option<String>,
        /// RNG seed
        #[arg(short, long, default_value_t = 42)]
        seed: u64,
        /// Number of threads (0 = auto)
        #[arg(short = 't', long, default_value_t = 0)]
        threads: usize,
        /// Output format: json, xyz, sdf
        #[arg(short, long, default_value = "json")]
        format: String,
    },
    /// Parse SMILES and show molecular structure (no 3D)
    Parse {
        /// SMILES string
        smiles: String,
    },
    /// Show version info
    Info,
}

fn format_xyz(result: &sci_form::ConformerResult) -> String {
    if result.error.is_some() || result.coords.is_empty() {
        return format!(
            "0\nERROR: {}\n",
            result.error.as_deref().unwrap_or("unknown")
        );
    }
    let n = result.num_atoms;
    let mut out = format!("{}\n{}\n", n, result.smiles);
    for i in 0..n {
        let elem = element_symbol(result.elements[i]);
        let x = result.coords[i * 3];
        let y = result.coords[i * 3 + 1];
        let z = result.coords[i * 3 + 2];
        out.push_str(&format!("{:2} {:12.6} {:12.6} {:12.6}\n", elem, x, y, z));
    }
    out
}

fn format_sdf(result: &sci_form::ConformerResult) -> String {
    if result.error.is_some() || result.coords.is_empty() {
        return format!(
            "\n  sci-form\n\nERROR: {}\n$$$$\n",
            result.error.as_deref().unwrap_or("unknown")
        );
    }
    let n = result.num_atoms;
    let nb = result.bonds.len();
    let mut out = format!("{}\n  sci-form\n\n", result.smiles);
    // Counts line: V2000
    out.push_str(&format!("{:3}{:3}  0  0  0  0  0  0  0  0999 V2000\n", n, nb));
    // Atom block
    for i in 0..n {
        let x = result.coords[i * 3];
        let y = result.coords[i * 3 + 1];
        let z = result.coords[i * 3 + 2];
        let elem = element_symbol(result.elements[i]);
        out.push_str(&format!(
            "{:10.4}{:10.4}{:10.4} {:3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
            x, y, z, elem
        ));
    }
    // Bond block
    for (a, b, order) in &result.bonds {
        let bond_type = match order.as_str() {
            "SINGLE" => 1,
            "DOUBLE" => 2,
            "TRIPLE" => 3,
            "AROMATIC" => 4,
            _ => 1,
        };
        out.push_str(&format!("{:3}{:3}{:3}  0\n", a + 1, b + 1, bond_type));
    }
    out.push_str("M  END\n$$$$\n");
    out
}

fn element_symbol(atomic_num: u8) -> &'static str {
    match atomic_num {
        1 => "H",
        6 => "C",
        7 => "N",
        8 => "O",
        9 => "F",
        15 => "P",
        16 => "S",
        17 => "Cl",
        35 => "Br",
        53 => "I",
        5 => "B",
        14 => "Si",
        34 => "Se",
        _ => "X",
    }
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Embed {
            smiles,
            seed,
            format,
        } => {
            let result = sci_form::embed(&smiles, seed);
            match format.as_str() {
                "xyz" => print!("{}", format_xyz(&result)),
                "sdf" => print!("{}", format_sdf(&result)),
                _ => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
            }
            if result.error.is_some() {
                std::process::exit(1);
            }
        }

        Commands::Batch {
            input,
            output,
            seed,
            threads,
            format,
        } => {
            let lines: Vec<String> = match input {
                Some(path) => {
                    let content = std::fs::read_to_string(&path)
                        .unwrap_or_else(|e| {
                            eprintln!("Error reading {}: {}", path, e);
                            std::process::exit(1);
                        });
                    content
                        .lines()
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
            eprintln!("Processing {} molecules ({} threads)...", total, if threads == 0 { "auto".to_string() } else { threads.to_string() });

            let smiles_refs: Vec<&str> = lines.iter().map(|s| s.as_str()).collect();
            let config = sci_form::ConformerConfig {
                seed,
                num_threads: threads,
            };
            let results = sci_form::embed_batch(&smiles_refs, &config);

            let elapsed = start.elapsed();

            // Write output
            let mut writer: Box<dyn Write> = match output {
                Some(path) => Box::new(
                    std::fs::File::create(&path).unwrap_or_else(|e| {
                        eprintln!("Error creating {}: {}", path, e);
                        std::process::exit(1);
                    }),
                ),
                None => Box::new(io::stdout()),
            };

            match format.as_str() {
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

            // Stats to stderr
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

        Commands::Parse { smiles } => match sci_form::parse(&smiles) {
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
        },

        Commands::Info => {
            println!("{}", sci_form::version());
            println!("Features: ETKDG distance geometry, CSD torsion patterns, SMARTS matching");
            println!("Formats: JSON, XYZ, SDF");
            println!("Bindings: CLI, Python (PyO3), TypeScript (WASM)");
        }
    }
}
