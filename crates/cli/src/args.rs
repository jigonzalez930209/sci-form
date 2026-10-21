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
    /// Compute a simplified NEB-like path between two geometries
    SimplifiedNebPath {
        /// SMILES string for topology and atom typing
        smiles: String,
        /// JSON array of flat start coords
        start_coords: String,
        /// JSON array of flat end coords
        end_coords: String,
        /// Number of images along the path
        #[arg(long, default_value_t = 8)]
        n_images: usize,
        /// Number of relaxation iterations
        #[arg(long, default_value_t = 50)]
        n_iter: usize,
        /// Harmonic spring constant between adjacent images
        #[arg(long, default_value_t = 0.1)]
        spring_k: f64,
        /// Cartesian update step size
        #[arg(long, default_value_t = 0.01)]
        step_size: f64,
        /// Energy backend: uff, mmff94, pm3, xtb, gfn1, gfn2, hf3c
        #[arg(long, default_value = "uff")]
        method: String,
    },
    /// Compute single-point energy with any NEB-capable backend
    NebEnergy {
        /// SMILES string for topology
        smiles: String,
        /// JSON array of flat coords
        coords: String,
        /// Energy backend: uff, mmff94, pm3, xtb, gfn1, gfn2, hf3c
        #[arg(long, default_value = "uff")]
        method: String,
    },
    /// Compute energy + gradient with any NEB-capable backend
    NebGradient {
        /// SMILES string for topology
        smiles: String,
        /// JSON array of flat coords
        coords: String,
        /// Energy backend: uff, mmff94, pm3, xtb, gfn1, gfn2, hf3c
        #[arg(long, default_value = "uff")]
        method: String,
    },
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
    /// Compute EEQ geometry-dependent charges
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

    /// [experimental-alpha] Compute an EDL profile
    #[cfg(feature = "alpha-edl")]
    EdlProfile {
        /// Surface potential (or surface charge) in V
        #[arg(short, long, default_value_t = 0.05)]
        potential: f64,
        /// EDL model: helmholtz, gouy-chapman, gouy-chapman-stern
        #[arg(short, long, default_value = "gouy-chapman-stern")]
        model: String,
        /// Bulk ionic strength in mol/L
        #[arg(short, long, default_value_t = 0.1)]
        ionic_strength: f64,
        /// Temperature in K
        #[arg(short, long, default_value_t = 298.15)]
        temperature: f64,
        /// Number of grid points
        #[arg(short, long, default_value_t = 128)]
        n_points: usize,
        /// Maximum profile extent in Å
        #[arg(short, long, default_value_t = 12.0)]
        extent: f64,
    },

    /// [experimental-alpha] Scan EDL capacitance over a potential range
    #[cfg(feature = "alpha-edl")]
    EdlCapacitanceScan {
        #[arg(long, default_value_t = -0.5)]
        v_min: f64,
        #[arg(long, default_value_t = 0.5)]
        v_max: f64,
        #[arg(short, long, default_value_t = 50)]
        n_points: usize,
        #[arg(short, long, default_value_t = 0.1)]
        ionic_strength: f64,
        /// EDL model: helmholtz, gouy-chapman, gouy-chapman-stern
        #[arg(short, long, default_value = "gouy-chapman-stern")]
        model: String,
    },

    /// [experimental-alpha] Compute a Helmholtz-only EDL profile
    #[cfg(feature = "alpha-edl")]
    EdlHelmholtz {
        #[arg(short, long, default_value_t = 0.05)]
        potential: f64,
        #[arg(short, long, default_value_t = 0.1)]
        ionic_strength: f64,
        #[arg(short, long, default_value_t = 298.15)]
        temperature: f64,
    },

    /// [experimental-alpha] Compute a Gouy-Chapman EDL profile
    #[cfg(feature = "alpha-edl")]
    EdlGouyChapman {
        #[arg(short, long, default_value_t = 0.05)]
        potential: f64,
        #[arg(short, long, default_value_t = 0.1)]
        ionic_strength: f64,
        #[arg(short, long, default_value_t = 298.15)]
        temperature: f64,
    },

    /// [experimental-alpha] Compute a Gouy-Chapman-Stern EDL profile
    #[cfg(feature = "alpha-edl")]
    EdlGcs {
        #[arg(short, long, default_value_t = 0.05)]
        potential: f64,
        #[arg(short, long, default_value_t = 0.1)]
        ionic_strength: f64,
        #[arg(short, long, default_value_t = 298.15)]
        temperature: f64,
        /// Stern-layer thickness in Å
        #[arg(short, long, default_value_t = 3.0)]
        stern_thickness: f64,
    },

    /// [experimental-alpha] Evaluate an HTST transition rate
    #[cfg(feature = "alpha-kinetics")]
    HtstRate {
        /// Activation free energy in eV
        #[arg(long)]
        activation_free_energy: f64,
        /// Reaction free energy in eV
        #[arg(long)]
        reaction_free_energy: f64,
        /// Temperature in K
        #[arg(short, long, default_value_t = 298.15)]
        temperature: f64,
    },

    /// [experimental-alpha] Generate a k-mesh
    #[cfg(feature = "alpha-periodic-linear")]
    Kmesh {
        /// Mesh divisions along the first reciprocal axis
        #[arg(long, default_value_t = 4)]
        nx: usize,
        /// Mesh divisions along the second reciprocal axis
        #[arg(long, default_value_t = 4)]
        ny: usize,
        /// Mesh divisions along the third reciprocal axis
        #[arg(long, default_value_t = 1)]
        nz: usize,
        /// Centering: monkhorst-pack or gamma
        #[arg(short, long, default_value = "monkhorst-pack")]
        centering: String,
    },

    /// [experimental-alpha] Generate EDL profile chart visualization
    #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-edl"))]
    EdlChart {
        #[arg(short, long, default_value_t = 0.1)]
        potential: f64,
        #[arg(short, long, default_value_t = 0.1)]
        ionic_strength: f64,
    },

    /// [experimental-alpha] Generate capacitance scan chart visualization
    #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-edl"))]
    CapacitanceChart {
        #[arg(long, default_value_t = -0.5)]
        v_min: f64,
        #[arg(long, default_value_t = 0.5)]
        v_max: f64,
        #[arg(short, long, default_value_t = 50)]
        n_points: usize,
        #[arg(short, long, default_value_t = 0.1)]
        ionic_strength: f64,
        #[arg(short, long, default_value = "gouy-chapman")]
        model: String,
    },

    /// [experimental-alpha] Generate Arrhenius chart
    #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-kinetics"))]
    ArrheniusChart {
        #[arg(short, long)]
        temperatures: String,
        #[arg(short, long)]
        rates: String,
    },

    /// [experimental-alpha] Generate thermal properties chart
    #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-kinetics"))]
    ThermalPropChart {
        /// Activation free energy in eV
        #[arg(long)]
        activation_free_energy: f64,
        /// Reaction free energy in eV
        #[arg(long)]
        reaction_free_energy: f64,
        /// JSON array of temperatures in K
        #[arg(short, long)]
        temperatures: String,
    },

    /// [experimental-alpha] Generate a phase portrait chart from a microkinetic network
    #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-kinetics"))]
    PhasePortraitChart {
        /// `MicrokineticNetworkConfig` JSON (file path or inline). Omit for the demo network.
        #[arg(long)]
        network: Option<String>,
        /// Number of integration frames
        #[arg(long, default_value_t = 50)]
        n_steps: usize,
    },

    /// [experimental-alpha] Generate band structure chart
    #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
    BandStructureChart {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords for the unit cell
        coords: String,
        /// JSON array of nine lattice-vector components
        lattice: String,
        /// k-points per high-symmetry segment
        #[arg(long, default_value_t = 20)]
        n_kpoints_per_segment: usize,
    },

    /// [experimental-alpha] Generate DOS chart
    #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
    DosChart {
        #[arg(short, long)]
        energies: String,
        #[arg(short, long)]
        dos_values: String,
    },

    /// [experimental-alpha] Generate k-point path chart
    #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
    KpointPathChart {
        /// Points per high-symmetry segment
        #[arg(long, default_value_t = 20)]
        n_points_per_segment: usize,
    },

    /// [experimental-alpha] Generate Fermi-surface slice chart
    #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
    FermiSurfaceChart {
        /// JSON array of real-space operator blocks
        blocks: String,
        /// Fermi energy in eV (estimated from the sampled band structure when omitted)
        #[arg(long)]
        fermi_energy_ev: Option<f64>,
        /// k-points along the first reciprocal axis
        #[arg(long, default_value_t = 24)]
        n_kx: usize,
        /// k-points along the second reciprocal axis
        #[arg(long, default_value_t = 24)]
        n_ky: usize,
    },

    /// [experimental-alpha] Generate trajectory chart from a GSM reaction path
    #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-gsm"))]
    TrajectoryChart {
        /// SMILES string for topology and atom typing
        smiles: String,
        /// JSON array of flat reactant coords
        reactant_coords: String,
        /// JSON array of flat product coords
        product_coords: String,
        /// GSM backend: uff, mmff94, pm3, xtb, gfn1, gfn2, hf3c
        #[arg(short, long, default_value = "uff")]
        method: String,
    },

    /// [experimental-alpha] Generate reaction coordinate chart from a GSM reaction path
    #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-gsm"))]
    ReactionCoordinateChart {
        /// SMILES string for topology and atom typing
        smiles: String,
        /// JSON array of flat reactant coords
        reactant_coords: String,
        /// JSON array of flat product coords
        product_coords: String,
        /// GSM backend: uff, mmff94, pm3, xtb, gfn1, gfn2, hf3c
        #[arg(short, long, default_value = "uff")]
        method: String,
    },

    // ─── PERIODIC LINEAR COMMANDS ──────────────
    /// [experimental-alpha] Compute the k-averaged periodic DOS with the KPM solver
    #[cfg(feature = "alpha-periodic-linear")]
    ComputePeriodicDos {
        /// JSON array of real-space Hamiltonian blocks
        hamiltonian_blocks: String,
        /// Mesh divisions along the first reciprocal axis
        #[arg(long, default_value_t = 2)]
        nx: usize,
        /// Mesh divisions along the second reciprocal axis
        #[arg(long, default_value_t = 2)]
        ny: usize,
        /// Mesh divisions along the third reciprocal axis
        #[arg(long, default_value_t = 2)]
        nz: usize,
        /// Centering: monkhorst-pack or gamma
        #[arg(short, long, default_value = "monkhorst-pack")]
        centering: String,
        /// Chebyshev expansion order
        #[arg(long, default_value_t = 100)]
        order: usize,
        /// Output minimum energy in eV
        #[arg(long, default_value_t = -20.0)]
        e_min: f64,
        /// Output maximum energy in eV
        #[arg(long, default_value_t = 10.0)]
        e_max: f64,
        /// Number of output grid points
        #[arg(long, default_value_t = 512)]
        n_points: usize,
        /// Fermi energy in eV
        #[arg(long)]
        fermi_energy_ev: Option<f64>,
    },

    /// [experimental-alpha] Solve periodic RandNLA
    #[cfg(feature = "alpha-periodic-linear")]
    SolvePeriodicRandnla {
        /// JSON array of real-space Hamiltonian blocks
        hamiltonian_blocks: String,
        /// JSON array of real-space overlap blocks
        overlap_blocks: String,
        /// Number of valence electrons
        #[arg(long, default_value_t = 2)]
        n_electrons: usize,
        /// Mesh divisions along the first reciprocal axis
        #[arg(long, default_value_t = 2)]
        nx: usize,
        /// Mesh divisions along the second reciprocal axis
        #[arg(long, default_value_t = 2)]
        ny: usize,
        /// Mesh divisions along the third reciprocal axis
        #[arg(long, default_value_t = 2)]
        nz: usize,
        /// Centering: monkhorst-pack or gamma
        #[arg(short, long, default_value = "monkhorst-pack")]
        centering: String,
        /// Nyström sketch rank
        #[arg(long)]
        sketch_size: Option<usize>,
    },

    /// [experimental-alpha] Compute Bloch phase
    #[cfg(feature = "alpha-periodic-linear")]
    BlochPhase {
        #[arg(long, default_value_t = 0.0)]
        k_x: f64,
        #[arg(long, default_value_t = 0.0)]
        k_y: f64,
        #[arg(long, default_value_t = 0.0)]
        k_z: f64,
        /// Lattice translation along the first axis
        #[arg(long, default_value_t = 1)]
        t_x: i32,
        /// Lattice translation along the second axis
        #[arg(long, default_value_t = 0)]
        t_y: i32,
        /// Lattice translation along the third axis
        #[arg(long, default_value_t = 0)]
        t_z: i32,
    },

    /// [experimental-alpha] Build Bloch Hamiltonian
    #[cfg(feature = "alpha-periodic-linear")]
    BuildBlochHamiltonian {
        /// JSON array of real-space operator blocks
        blocks: String,
        #[arg(long, default_value_t = 0.0)]
        k_x: f64,
        #[arg(long, default_value_t = 0.0)]
        k_y: f64,
        #[arg(long, default_value_t = 0.0)]
        k_z: f64,
    },

    /// [experimental-alpha] Assemble periodic operators for every k-point
    #[cfg(feature = "alpha-periodic-linear")]
    AssemblePeriodicOperators {
        /// JSON array of real-space Hamiltonian blocks
        hamiltonian_blocks: String,
        /// JSON array of real-space overlap blocks (optional)
        #[arg(default_value = "")]
        overlap_blocks: String,
        /// Mesh divisions along the first reciprocal axis
        #[arg(long, default_value_t = 2)]
        nx: usize,
        /// Mesh divisions along the second reciprocal axis
        #[arg(long, default_value_t = 2)]
        ny: usize,
        /// Mesh divisions along the third reciprocal axis
        #[arg(long, default_value_t = 2)]
        nz: usize,
        /// Centering: monkhorst-pack or gamma
        #[arg(short, long, default_value = "monkhorst-pack")]
        centering: String,
    },

    /// [experimental-alpha] Validate electron count
    #[cfg(feature = "alpha-periodic-linear")]
    ValidateElectronCount {
        #[arg(long)]
        n_electrons: usize,
        #[arg(long)]
        n_bands: usize,
        #[arg(long)]
        n_kpoints: usize,
    },

    /// [experimental-alpha] Show available GSM reaction-path backends for a molecule
    #[cfg(feature = "alpha-gsm")]
    GsmBackendPlan {
        /// SMILES string used to derive the element set
        smiles: String,
    },

    /// [experimental-alpha] Compare GSM backends on the same geometry
    #[cfg(feature = "alpha-gsm")]
    GsmCompareBackends {
        /// SMILES string for topology and atom typing
        smiles: String,
        /// JSON array of flat xyz coords
        coords: String,
        /// JSON array of backend names, or [] for all
        #[arg(long, default_value = "[]")]
        methods: String,
    },

    // ─── KINETICS COMMANDS ────────────────────
    /// [experimental-alpha] Extract kinetics diagnostics from a microkinetic trace
    #[cfg(feature = "alpha-kinetics")]
    ExtractKineticsDiagnostics {
        /// `MicrokineticNetworkConfig` as JSON (file path or inline). Omit for
        /// the built-in demo network.
        #[arg(long)]
        network: Option<String>,
    },

    /// [experimental-alpha] Solve microkinetic network
    #[cfg(feature = "alpha-kinetics")]
    SolveMicrokineticNetwork {
        /// `MicrokineticNetworkConfig` as JSON (file path or inline). Omit for
        /// the built-in demo network.
        #[arg(long)]
        network: Option<String>,
        /// Number of integration frames.
        #[arg(long, default_value_t = 5)]
        n_steps: usize,
    },

    /// [experimental-alpha] Solve microkinetic steady state
    #[cfg(feature = "alpha-kinetics")]
    SolveMicrokineticSteadyState {
        /// `MicrokineticNetworkConfig` as JSON (file path or inline). Omit for
        /// the built-in demo network.
        #[arg(long)]
        network: Option<String>,
        /// Number of integration frames.
        #[arg(long, default_value_t = 5)]
        n_steps: usize,
        /// Convergence tolerance on the population change.
        #[arg(long, default_value_t = 1e-10)]
        tol: f64,
    },

    /// [experimental-alpha] Analyze GSM+MBH+HTST step
    #[cfg(feature = "alpha-kinetics")]
    AnalyzeGsmMbhHtstStep {
        /// SMILES string for topology, ring perception, and backend capability checks
        smiles: String,
        /// JSON array of flat reactant coords
        reactant_coords: String,
        /// JSON array of flat product coords
        product_coords: String,
        /// GSM backend: uff, mmff94, pm3, xtb, gfn1, gfn2, hf3c, reaxff
        #[arg(long, default_value = "uff")]
        method: String,
        /// Stable identifier for the elementary step
        #[arg(long, default_value = "elementary-step")]
        step_id: String,
        /// Temperature in K
        #[arg(long, default_value_t = 298.15)]
        temperature_k: f64,
        /// Pressure in bar
        #[arg(long, default_value_t = 1.0)]
        pressure_bar: f64,
        /// Maximum number of GSM nodes
        #[arg(long, default_value_t = 11)]
        n_nodes: usize,
        /// MBH finite-difference step in Å
        #[arg(long, default_value_t = 0.005)]
        mbh_fd_step: f64,
    },

    // ──── RENDER BRIDGE: chart utilities ───────────────────────────────
    /// [render-bridge] Serialize chart to JSON
    #[command(args_conflicts_with_subcommands = true)]
    ChartToJson {
        /// JSON string input (a chart payload or the chart wrapper)
        #[arg(value_name = "JSON")]
        json_input: Option<String>,
    },

    /// [render-bridge] Deserialize chart from JSON
    #[command(args_conflicts_with_subcommands = true)]
    ChartFromJson {
        /// JSON string input
        #[arg(value_name = "JSON")]
        json_input: Option<String>,
    },

    /// [render-bridge] Validate chart schema
    #[command(args_conflicts_with_subcommands = true)]
    ValidateChartSchema {
        /// JSON string to validate
        #[arg(value_name = "JSON")]
        json_input: Option<String>,
    },

    /// [alpha-reaction-dynamics] Compute full 3D reaction dynamics path
    #[cfg(feature = "alpha-reaction-dynamics")]
    ReactionDynamics3d {
        /// JSON array of reactant SMILES, e.g. '["[Cl-]", "CBr"]'
        reactant_smiles: String,
        /// JSON array of product SMILES, e.g. '["ClC", "[Br-]"]'
        product_smiles: String,
        /// Energy backend: uff, xtb, gfn1, gfn2 (default), pm3, hf3c
        #[arg(long, default_value = "gfn2")]
        method: String,
        /// Number of NEB images
        #[arg(long, default_value_t = 30)]
        n_images: usize,
        /// Number of approach frames
        #[arg(long, default_value_t = 15)]
        n_approach: usize,
        /// Number of departure frames
        #[arg(long, default_value_t = 15)]
        n_departure: usize,
        /// Optional SMIRKS reaction string
        #[arg(long)]
        smirks: Option<String>,
        /// Random seed
        #[arg(long, default_value_t = 42)]
        seed: u64,
    },
}
