use clap::{Parser, Subcommand};
use std::io::{self, BufRead, Write};
use std::time::Instant;

#[derive(Parser)]
#[command(
    name = "sci-form",
    version,
    about = "3D molecular conformer generation"
)]
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
    /// Run Extended Hückel Theory calculation
    Eht {
        /// JSON array of atomic numbers, e.g. "[8,1,1]"
        elements: String,
        /// JSON array of flat xyz coords, e.g. "[0,0,0, 0.96,0,0, -0.24,0.93,0]"
        coords: String,
        /// Wolfsberg-Helmholtz constant (0 = default 1.75)
        #[arg(short, long, default_value_t = 0.0)]
        k: f64,
    },
    /// Compute Gasteiger-Marsili partial charges from SMILES
    Charges {
        /// SMILES string
        smiles: String,
    },
    /// Compute solvent-accessible surface area
    Sasa {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords
        coords: String,
        /// Probe radius in Å (default 1.4)
        #[arg(short, long, default_value_t = 1.4)]
        probe_radius: f64,
    },
    /// Mulliken & Löwdin population analysis (requires EHT)
    Population {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords
        coords: String,
    },
    /// Compute molecular dipole moment (Debye)
    Dipole {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords
        coords: String,
    },
    /// Compute electrostatic potential on a 3D grid
    Esp {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords
        coords: String,
        /// Grid spacing in Å
        #[arg(short, long, default_value_t = 0.5)]
        spacing: f64,
        /// Padding around molecule in Å
        #[arg(short, long, default_value_t = 3.0)]
        padding: f64,
    },
    /// Compute density of states (DOS/PDOS)
    Dos {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords
        coords: String,
        /// Gaussian smearing width (eV)
        #[arg(short, long, default_value_t = 0.3)]
        sigma: f64,
        /// Energy window minimum (eV)
        #[arg(long, default_value_t = -30.0)]
        e_min: f64,
        /// Energy window maximum (eV)
        #[arg(long, default_value_t = 5.0)]
        e_max: f64,
        /// Number of grid points
        #[arg(short, long, default_value_t = 500)]
        n_points: usize,
    },
    /// Compute RMSD between two coordinate sets (Kabsch alignment)
    Rmsd {
        /// JSON array of flat xyz coords (mobile)
        coords: String,
        /// JSON array of flat xyz coords (reference)
        reference: String,
    },
    /// Compute UFF force field energy
    Uff {
        /// SMILES string
        smiles: String,
        /// JSON array of flat xyz coords
        coords: String,
    },
    /// Create a unit cell and show parameters
    Cell {
        /// Lattice parameter a (Å)
        #[arg(long)]
        a: f64,
        /// Lattice parameter b (Å)
        #[arg(long)]
        b: f64,
        /// Lattice parameter c (Å)
        #[arg(long)]
        c: f64,
        /// Angle α (degrees)
        #[arg(long, default_value_t = 90.0)]
        alpha: f64,
        /// Angle β (degrees)
        #[arg(long, default_value_t = 90.0)]
        beta: f64,
        /// Angle γ (degrees)
        #[arg(long, default_value_t = 90.0)]
        gamma: f64,
    },
    /// Assemble a framework from a topology
    Assemble {
        /// Topology name: pcu, dia, sql
        #[arg(short, long, default_value = "pcu")]
        topology: String,
        /// Lattice parameter a (Å)
        #[arg(long, default_value_t = 10.0)]
        a: f64,
        /// Metal center atomic number
        #[arg(long, default_value_t = 30)]
        metal: u8,
        /// Coordination geometry: linear, tetrahedral, octahedral, square_planar
        #[arg(long, default_value = "octahedral")]
        geometry: String,
        /// Supercell replication (e.g., 2 for 2x2x2)
        #[arg(long, default_value_t = 1)]
        supercell: usize,
    },
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
    out.push_str(&format!(
        "{:3}{:3}  0  0  0  0  0  0  0  0999 V2000\n",
        n, nb
    ));
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
                    let content = std::fs::read_to_string(&path).unwrap_or_else(|e| {
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

            // Write output
            let mut writer: Box<dyn Write> = match output {
                Some(path) => Box::new(std::fs::File::create(&path).unwrap_or_else(|e| {
                    eprintln!("Error creating {}: {}", path, e);
                    std::process::exit(1);
                })),
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
            println!("  EHT, Gasteiger charges, SASA");
            println!("Formats: JSON, XYZ, SDF");
            println!("Bindings: CLI, Python (PyO3), TypeScript (WASM)");
        }

        Commands::Eht {
            elements,
            coords,
            k,
        } => {
            let elems: Vec<u8> = serde_json::from_str(&elements).unwrap_or_else(|e| {
                eprintln!("Bad elements JSON: {}", e);
                std::process::exit(1);
            });
            let flat: Vec<f64> = serde_json::from_str(&coords).unwrap_or_else(|e| {
                eprintln!("Bad coords JSON: {}", e);
                std::process::exit(1);
            });
            if flat.len() != elems.len() * 3 {
                eprintln!(
                    "coords length {} != 3 * elements {}",
                    flat.len(),
                    elems.len()
                );
                std::process::exit(1);
            }
            let positions: Vec<[f64; 3]> =
                flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
            let k_opt = if k <= 0.0 { None } else { Some(k) };
            match sci_form::eht::solve_eht(&elems, &positions, k_opt) {
                Ok(result) => {
                    println!("{}", serde_json::to_string_pretty(&result).unwrap());
                }
                Err(e) => {
                    eprintln!("EHT error: {}", e);
                    std::process::exit(1);
                }
            }
        }

        Commands::Charges { smiles } => {
            match sci_form::compute_charges(&smiles) {
                Ok(result) => {
                    println!("{}", serde_json::to_string_pretty(&result).unwrap());
                }
                Err(e) => {
                    eprintln!("Charges error: {}", e);
                    std::process::exit(1);
                }
            }
        }

        Commands::Sasa {
            elements,
            coords,
            probe_radius,
        } => {
            let elems: Vec<u8> = serde_json::from_str(&elements).unwrap_or_else(|e| {
                eprintln!("Bad elements JSON: {}", e);
                std::process::exit(1);
            });
            let flat: Vec<f64> = serde_json::from_str(&coords).unwrap_or_else(|e| {
                eprintln!("Bad coords JSON: {}", e);
                std::process::exit(1);
            });
            match sci_form::compute_sasa(&elems, &flat, Some(probe_radius)) {
                Ok(result) => {
                    println!("{}", serde_json::to_string_pretty(&result).unwrap());
                }
                Err(e) => {
                    eprintln!("SASA error: {}", e);
                    std::process::exit(1);
                }
            }
        }

        Commands::Population { elements, coords } => {
            let elems: Vec<u8> = serde_json::from_str(&elements).unwrap_or_else(|e| {
                eprintln!("Bad elements JSON: {}", e);
                std::process::exit(1);
            });
            let flat: Vec<f64> = serde_json::from_str(&coords).unwrap_or_else(|e| {
                eprintln!("Bad coords JSON: {}", e);
                std::process::exit(1);
            });
            let positions: Vec<[f64; 3]> =
                flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
            match sci_form::compute_population(&elems, &positions) {
                Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
                Err(e) => {
                    eprintln!("Population error: {}", e);
                    std::process::exit(1);
                }
            }
        }

        Commands::Dipole { elements, coords } => {
            let elems: Vec<u8> = serde_json::from_str(&elements).unwrap_or_else(|e| {
                eprintln!("Bad elements JSON: {}", e);
                std::process::exit(1);
            });
            let flat: Vec<f64> = serde_json::from_str(&coords).unwrap_or_else(|e| {
                eprintln!("Bad coords JSON: {}", e);
                std::process::exit(1);
            });
            let positions: Vec<[f64; 3]> =
                flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
            match sci_form::compute_dipole(&elems, &positions) {
                Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
                Err(e) => {
                    eprintln!("Dipole error: {}", e);
                    std::process::exit(1);
                }
            }
        }

        Commands::Esp {
            elements,
            coords,
            spacing,
            padding,
        } => {
            let elems: Vec<u8> = serde_json::from_str(&elements).unwrap_or_else(|e| {
                eprintln!("Bad elements JSON: {}", e);
                std::process::exit(1);
            });
            let flat: Vec<f64> = serde_json::from_str(&coords).unwrap_or_else(|e| {
                eprintln!("Bad coords JSON: {}", e);
                std::process::exit(1);
            });
            let positions: Vec<[f64; 3]> =
                flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
            match sci_form::compute_esp(&elems, &positions, spacing, padding) {
                Ok(grid) => {
                    // Print dimensions + origin only (values too large for stdout)
                    println!(
                        "{{\"origin\":{:?},\"spacing\":{},\"dims\":{:?},\"n_values\":{}}}",
                        grid.origin,
                        grid.spacing,
                        grid.dims,
                        grid.values.len()
                    );
                }
                Err(e) => {
                    eprintln!("ESP error: {}", e);
                    std::process::exit(1);
                }
            }
        }

        Commands::Dos {
            elements,
            coords,
            sigma,
            e_min,
            e_max,
            n_points,
        } => {
            let elems: Vec<u8> = serde_json::from_str(&elements).unwrap_or_else(|e| {
                eprintln!("Bad elements JSON: {}", e);
                std::process::exit(1);
            });
            let flat: Vec<f64> = serde_json::from_str(&coords).unwrap_or_else(|e| {
                eprintln!("Bad coords JSON: {}", e);
                std::process::exit(1);
            });
            let positions: Vec<[f64; 3]> =
                flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
            match sci_form::compute_dos(&elems, &positions, sigma, e_min, e_max, n_points) {
                Ok(result) => {
                    // Print energy-DOS pairs as JSON array
                    let pairs: Vec<_> = result
                        .energies
                        .iter()
                        .zip(result.total_dos.iter())
                        .map(|(e, d)| format!("[{:.4},{:.6}]", e, d))
                        .collect();
                    println!("{{\"sigma\":{},\"data\":[{}]}}", result.sigma, pairs.join(","));
                }
                Err(e) => {
                    eprintln!("DOS error: {}", e);
                    std::process::exit(1);
                }
            }
        }

        Commands::Rmsd { coords, reference } => {
            let c: Vec<f64> = serde_json::from_str(&coords).unwrap_or_else(|e| {
                eprintln!("Bad coords JSON: {}", e);
                std::process::exit(1);
            });
            let r: Vec<f64> = serde_json::from_str(&reference).unwrap_or_else(|e| {
                eprintln!("Bad reference JSON: {}", e);
                std::process::exit(1);
            });
            let result = sci_form::alignment::align_coordinates(&c, &r);
            println!(
                "{{\"rmsd\":{:.6},\"rotation\":{:?}}}",
                result.rmsd, result.rotation
            );
        }

        Commands::Uff { smiles, coords } => {
            let flat: Vec<f64> = serde_json::from_str(&coords).unwrap_or_else(|e| {
                eprintln!("Bad coords JSON: {}", e);
                std::process::exit(1);
            });
            match sci_form::compute_uff_energy(&smiles, &flat) {
                Ok(energy) => println!("{{\"energy\":{:.6},\"unit\":\"kcal/mol\"}}", energy),
                Err(e) => {
                    eprintln!("UFF error: {}", e);
                    std::process::exit(1);
                }
            }
        }

        Commands::Cell { a, b, c, alpha, beta, gamma } => {
            let cell = sci_form::create_unit_cell(a, b, c, alpha, beta, gamma);
            let vol = cell.volume();
            let p = cell.parameters();
            println!(
                "{{\"a\":{:.4},\"b\":{:.4},\"c\":{:.4},\"alpha\":{:.2},\"beta\":{:.2},\"gamma\":{:.2},\"volume\":{:.4},\"lattice\":{:?}}}",
                p.a, p.b, p.c, p.alpha, p.beta, p.gamma, vol, cell.lattice
            );
        }

        Commands::Assemble { topology, a, metal, geometry, supercell } => {
            let geom = match geometry.as_str() {
                "linear" => sci_form::materials::CoordinationGeometry::Linear,
                "trigonal" => sci_form::materials::CoordinationGeometry::Trigonal,
                "tetrahedral" => sci_form::materials::CoordinationGeometry::Tetrahedral,
                "square_planar" => sci_form::materials::CoordinationGeometry::SquarePlanar,
                "octahedral" => sci_form::materials::CoordinationGeometry::Octahedral,
                _ => {
                    eprintln!("Unknown geometry: {}", geometry);
                    std::process::exit(1);
                }
            };
            let topo = match topology.as_str() {
                "pcu" => sci_form::materials::Topology::pcu(),
                "dia" => sci_form::materials::Topology::dia(),
                "sql" => sci_form::materials::Topology::sql(),
                _ => {
                    eprintln!("Unknown topology: {}", topology);
                    std::process::exit(1);
                }
            };
            let node = sci_form::materials::Sbu::metal_node(metal, 0.0, geom);
            let linker = sci_form::materials::Sbu::linear_linker(&[6, 6], 1.4, "carboxylate");
            let cell = sci_form::materials::UnitCell::cubic(a);
            let mut structure = sci_form::assemble_framework(&node, &linker, &topo, &cell);
            if supercell > 1 {
                structure = structure.make_supercell(supercell, supercell, supercell);
            }
            println!("{}", serde_json::to_string_pretty(&structure).unwrap());
        }
    }
}
