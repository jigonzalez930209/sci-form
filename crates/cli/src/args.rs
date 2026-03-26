//! CLI argument definitions.

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "sci-form",
    version,
    about = "3D molecular conformer generation"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
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
    Parse { smiles: String },
    /// Show version info
    Info,
    /// Run Extended Hückel Theory calculation
    Eht {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords
        coords: String,
        /// Wolfsberg-Helmholtz constant (0 = default 1.75)
        #[arg(short, long, default_value_t = 0.0)]
        k: f64,
    },
    /// Compute Gasteiger-Marsili partial charges from SMILES
    Charges { smiles: String },
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
    Population { elements: String, coords: String },
    /// Compute molecular dipole moment (Debye)
    Dipole { elements: String, coords: String },
    /// Compute electrostatic potential on a 3D grid
    Esp {
        elements: String,
        coords: String,
        #[arg(short, long, default_value_t = 0.5)]
        spacing: f64,
        #[arg(short, long, default_value_t = 3.0)]
        padding: f64,
    },
    /// Compute density of states (DOS/PDOS)
    Dos {
        elements: String,
        coords: String,
        #[arg(short, long, default_value_t = 0.3)]
        sigma: f64,
        #[arg(long, default_value_t = -30.0)]
        e_min: f64,
        #[arg(long, default_value_t = 5.0)]
        e_max: f64,
        #[arg(short, long, default_value_t = 500)]
        n_points: usize,
    },
    /// Compute RMSD between two coordinate sets (Kabsch alignment)
    Rmsd { coords: String, reference: String },
    /// Compute UFF force field energy
    Uff { smiles: String, coords: String },
    /// Create a unit cell and show parameters
    Cell {
        #[arg(long)]
        a: f64,
        #[arg(long)]
        b: f64,
        #[arg(long)]
        c: f64,
        #[arg(long, default_value_t = 90.0)]
        alpha: f64,
        #[arg(long, default_value_t = 90.0)]
        beta: f64,
        #[arg(long, default_value_t = 90.0)]
        gamma: f64,
    },
    /// Assemble a framework from a topology
    Assemble {
        #[arg(short, long, default_value = "pcu")]
        topology: String,
        #[arg(long, default_value_t = 10.0)]
        a: f64,
        #[arg(long, default_value_t = 30)]
        metal: u8,
        #[arg(long, default_value = "octahedral")]
        geometry: String,
        #[arg(long, default_value_t = 1)]
        supercell: usize,
    },
    /// Run ANI ML Potential calculation
    Ani { elements: String, coords: String },
    /// Run PM3 semi-empirical quantum calculation
    Pm3 { elements: String, coords: String },
    /// Run GFN0-xTB tight-binding quantum calculation
    Xtb { elements: String, coords: String },
    /// Run GFN1-xTB tight-binding quantum calculation
    Gfn1 { elements: String, coords: String },
    /// Run GFN2-xTB tight-binding quantum calculation
    Gfn2 { elements: String, coords: String },
    /// Run HF-3c quantum calculation
    Hf3c { elements: String, coords: String },
    /// Analyze stereochemistry (R/S stereocenters, E/Z double bonds)
    Stereo {
        /// SMILES string
        smiles: String,
        /// JSON array of flat xyz coords (optional, empty for topology-only)
        #[arg(short, long, default_value = "[]")]
        coords: String,
    },
    /// Compute non-polar solvation energy (SASA + ASP)
    Solvation {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords
        coords: String,
        /// JSON array of partial charges (for GB; omit for non-polar only)
        #[arg(short, long, default_value = "")]
        charges: String,
        /// Probe radius in Å
        #[arg(short, long, default_value_t = 1.4)]
        probe_radius: f64,
    },
    /// Compute SSSR (Smallest Set of Smallest Rings)
    Sssr { smiles: String },
    /// Compute ECFP fingerprint
    Ecfp {
        smiles: String,
        /// Fingerprint radius (2 = ECFP4, 3 = ECFP6)
        #[arg(short, long, default_value_t = 2)]
        radius: usize,
        /// Bit vector length
        #[arg(short, long, default_value_t = 2048)]
        n_bits: usize,
    },
    /// Compute Tanimoto similarity between two SMILES
    Tanimoto {
        smiles1: String,
        smiles2: String,
        #[arg(short, long, default_value_t = 2)]
        radius: usize,
        #[arg(short, long, default_value_t = 2048)]
        n_bits: usize,
    },

    // ─── Experimental ──────────────────────────────────────────────────
    /// [experimental] Compute EEQ geometry-dependent charges
    #[cfg(feature = "experimental-eeq")]
    Eeq {
        elements: String,
        coords: String,
        #[arg(short, long, default_value_t = 0.0)]
        total_charge: f64,
    },

    /// [experimental] Compute ALPB implicit solvation energy
    #[cfg(feature = "experimental-alpb")]
    Alpb {
        elements: String,
        coords: String,
        #[arg(short, long, default_value = "")]
        charges: String,
        #[arg(short, long, default_value_t = 78.5)]
        dielectric: f64,
    },

    /// [experimental] Compute DFT-D4 dispersion correction
    #[cfg(feature = "experimental-d4")]
    D4 {
        elements: String,
        coords: String,
        #[arg(long, default_value_t = false)]
        three_body: bool,
    },

    /// [experimental] Compute CPM charges at fixed electrochemical potential
    #[cfg(feature = "experimental-cpm")]
    Cpm {
        elements: String,
        coords: String,
        #[arg(short, long, default_value_t = -4.44)]
        mu_ev: f64,
        #[arg(short, long, default_value_t = 78.5)]
        dielectric: f64,
    },
}
