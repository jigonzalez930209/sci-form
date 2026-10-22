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
        /// Wolfsberg-Helmholtz constant (omit the flag for the default 1.75)
        #[arg(short, long)]
        k: Option<f64>,
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
        /// JSON array of flat xyz coords; omit the flag for a topology-only
        /// analysis
        #[arg(short, long)]
        coords: Option<String>,
    },
    /// Compute non-polar solvation energy (SASA + ASP)
    Solvation {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords
        coords: String,
        /// JSON array of partial charges (for GB); omit the flag for the
        /// non-polar contribution only
        #[arg(short, long)]
        charges: Option<String>,
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

    // ─── Cross-surface parity commands ─────────────────────────────────
    /// Generate a diverse conformer ensemble and cluster it by Butina RMSD
    EmbedDiverse {
        /// SMILES string
        smiles: String,
        /// Number of embedding attempts (different seeds)
        #[arg(long, default_value_t = 10)]
        n_conformers: usize,
        /// RMSD clustering cutoff in Å
        #[arg(long, default_value_t = 1.0)]
        rmsd_cutoff: f64,
        /// Base RNG seed (seeds = base_seed..base_seed+n_conformers)
        #[arg(long, default_value_t = 42)]
        base_seed: u64,
    },
    /// Print the canonical SMILES representation of a molecule
    CanonicalSmiles { smiles: String },
    /// Compute ECFP fingerprints for a batch of SMILES (JSON array)
    EcfpBatch {
        /// JSON array of SMILES strings
        smiles_list: String,
        /// Fingerprint radius (2 = ECFP4, 3 = ECFP6)
        #[arg(short, long, default_value_t = 2)]
        radius: usize,
        /// Bit vector length
        #[arg(short, long, default_value_t = 2048)]
        n_bits: usize,
    },
    /// Compute Gasteiger-Marsili charges with configurable parameters
    ChargesConfigured {
        /// SMILES string
        smiles: String,
        /// Number of electronegativity-equalization iterations
        #[arg(long, default_value_t = 6)]
        max_iter: usize,
        /// Initial damping factor (halved every iteration)
        #[arg(long, default_value_t = 0.5)]
        initial_damping: f64,
        /// Convergence threshold for max charge transfer
        #[arg(long, default_value_t = 1e-10)]
        convergence_threshold: f64,
    },
    /// Predict J-couplings averaged over a conformer ensemble (Boltzmann weighting)
    EnsembleJCouplings {
        /// SMILES string
        smiles: String,
        /// JSON array of conformers, each a flat xyz coordinate array
        conformer_coords: String,
        /// JSON array of conformer energies in kcal/mol
        energies_kcal: String,
        /// Temperature in K
        #[arg(long, default_value_t = 298.15)]
        temperature_k: f64,
    },
    /// Compute a broadened IR spectrum (runs the vibrational analysis internally)
    IrSpectrumBroadened {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords
        coords: String,
        /// Hessian method: eht, pm3, xtb, uff, ...
        #[arg(long, default_value = "eht")]
        method: String,
        /// Broadening width in cm⁻¹
        #[arg(long, default_value_t = 15.0)]
        gamma: f64,
        /// Minimum wavenumber in cm⁻¹
        #[arg(long, default_value_t = 400.0)]
        wn_min: f64,
        /// Maximum wavenumber in cm⁻¹
        #[arg(long, default_value_t = 4000.0)]
        wn_max: f64,
        /// Number of grid points
        #[arg(long, default_value_t = 1000)]
        n_points: usize,
        /// Broadening type: gaussian or lorentzian
        #[arg(long, default_value = "lorentzian")]
        broadening: String,
        /// Finite-difference step in Å; omit the flag for the engine default
        #[arg(long)]
        step_size: Option<f64>,
        /// SMILES string (required for method=uff)
        #[arg(long)]
        smiles: Option<String>,
    },
    /// Compute CISD excited states from an internal STO-3G RHF calculation
    Cisd {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords
        coords: String,
        /// Number of excited states to report
        #[arg(long, default_value_t = 5)]
        n_states: usize,
    },
    /// Compute DOS from the orbital energies of a chosen electronic-structure method
    DosMultimethod {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords
        coords: String,
        /// Method: eht, pm3, xtb, gfn1, gfn2 or hf3c
        #[arg(long, default_value = "eht")]
        method: String,
        #[arg(short, long, default_value_t = 0.3)]
        sigma: f64,
        #[arg(long, default_value_t = -30.0)]
        e_min: f64,
        #[arg(long, default_value_t = 5.0)]
        e_max: f64,
        #[arg(short, long, default_value_t = 500)]
        n_points: usize,
    },
    /// Compute EHT energy-surface gradients (central finite differences)
    EhtGradient {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords
        coords: String,
    },
    /// Optimize a geometry on the EHT energy surface (shared BFGS engine)
    OptimizeGeometryEht {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords
        coords: String,
        /// Maximum number of optimization steps
        #[arg(long, default_value_t = 200)]
        max_steps: usize,
        /// RMS gradient convergence threshold in eV/Å
        #[arg(long, default_value_t = 0.01)]
        grad_threshold: f64,
        /// Energy convergence threshold in eV
        #[arg(long, default_value_t = 1e-6)]
        energy_threshold: f64,
        /// Initial step size for steepest descent in Å
        #[arg(long, default_value_t = 0.05)]
        step_size: f64,
    },
    /// Compute ANI-TM atomic environment vectors (24-element extended set)
    AevsTm {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords
        coords: String,
    },
    /// Compute WHIM 3D molecular descriptors
    Whim {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords
        coords: String,
        /// Weighting: unit, mass, volume, electronegativity or polarizability
        #[arg(long, default_value = "mass")]
        weighting: String,
    },
    /// Compute radial distribution function (RDF) descriptors
    Rdf {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords
        coords: String,
        /// JSON array of partial charges; omit the flag to derive them with
        /// Gasteiger
        #[arg(long)]
        charges: Option<String>,
        /// Gaussian smoothing factor in Å⁻²
        #[arg(long, default_value_t = 100.0)]
        beta: f64,
        /// Maximum radius in Å
        #[arg(long, default_value_t = 12.0)]
        r_max: f64,
        /// Number of sampling points
        #[arg(long, default_value_t = 128)]
        n_points: usize,
    },
    /// Compute GETAWAY 3D molecular descriptors
    Getaway {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords
        coords: String,
        /// JSON array of [i, j] bond pairs; omit the flag to perceive bonds
        /// from the geometry
        #[arg(long)]
        bonds: Option<String>,
        /// Maximum topological distance for the autocorrelations
        #[arg(long, default_value_t = 8)]
        max_lag: usize,
    },
    /// Train a Random Forest regressor on a feature matrix
    TrainRandomForest {
        /// JSON array of feature vectors
        features: String,
        /// JSON array of target values
        targets: String,
        /// Number of trees
        #[arg(long, default_value_t = 100)]
        n_trees: usize,
        /// Maximum tree depth
        #[arg(long, default_value_t = 10)]
        max_depth: usize,
        /// Minimum samples required to split a node
        #[arg(long, default_value_t = 5)]
        min_samples_split: usize,
        /// Minimum samples required for a leaf
        #[arg(long, default_value_t = 2)]
        min_samples_leaf: usize,
        /// Bootstrap sample fraction (0.0-1.0)
        #[arg(long, default_value_t = 1.0)]
        sample_fraction: f64,
        /// RNG seed
        #[arg(long, default_value_t = 42)]
        seed: u64,
    },
    /// Train a Gradient Boosting regressor on a feature matrix
    TrainGradientBoosting {
        /// JSON array of feature vectors
        features: String,
        /// JSON array of target values
        targets: String,
        /// Number of boosting stages
        #[arg(long, default_value_t = 100)]
        n_estimators: usize,
        /// Learning rate
        #[arg(long, default_value_t = 0.1)]
        learning_rate: f64,
        /// Maximum tree depth
        #[arg(long, default_value_t = 5)]
        max_depth: usize,
        /// Minimum samples required to split a node
        #[arg(long, default_value_t = 5)]
        min_samples_split: usize,
        /// Minimum samples required for a leaf
        #[arg(long, default_value_t = 2)]
        min_samples_leaf: usize,
        /// Row subsample fraction (0.0-1.0)
        #[arg(long, default_value_t = 0.8)]
        subsample: f64,
        /// RNG seed
        #[arg(long, default_value_t = 42)]
        seed: u64,
    },
    /// Assign IR peaks to functional groups (group-frequency tables)
    AssignPeaks {
        /// JSON array of vibrational frequencies in cm⁻¹
        frequencies: String,
        /// JSON array of IR intensities
        intensities: String,
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of modes, each an array of [x,y,z] triples; omit the
        /// flag when the modes are unknown
        #[arg(long)]
        mode_displacements: Option<String>,
    },
    /// Look up a crystallographic space group by number or symbol
    SpaceGroup {
        /// ITC number (1-230)
        #[arg(long)]
        number: Option<u16>,
        /// Hermann-Mauguin symbol, e.g. "P 21/c"
        #[arg(long)]
        symbol: Option<String>,
    },
    /// List all 230 space groups
    AllSpaceGroups,
    /// Build a periodic molecular graph from fractional coordinates and a unit cell
    PeriodicMolecule {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of fractional coordinates
        frac_coords: String,
        /// Lattice length a in Å
        #[arg(long)]
        a: f64,
        /// Lattice length b in Å
        #[arg(long)]
        b: f64,
        /// Lattice length c in Å
        #[arg(long)]
        c: f64,
        /// Lattice angle alpha in degrees
        #[arg(long, default_value_t = 90.0)]
        alpha: f64,
        /// Lattice angle beta in degrees
        #[arg(long, default_value_t = 90.0)]
        beta: f64,
        /// Lattice angle gamma in degrees
        #[arg(long, default_value_t = 90.0)]
        gamma: f64,
        /// Extra bonding tolerance in Å; omit the flag for the default 0.3
        #[arg(long)]
        bond_tolerance: Option<f64>,
    },
    /// Detect haptic (η) metal-ring interactions and metallocene motifs
    Hapticity {
        /// JSON array of atomic numbers
        elements: String,
        /// JSON array of flat xyz coords
        coords: String,
        /// JSON array of [i, j] bond pairs; omit the flag to perceive bonds
        /// from the geometry
        #[arg(long)]
        bonds: Option<String>,
    },
    /// Optimize a geometry with the shared BFGS engine (UFF backend)
    OptimizeFramework {
        /// SMILES string (used to build the UFF force field)
        smiles: String,
        /// JSON array of flat xyz coords
        coords: String,
        /// Optimization method: bfgs or steepest_descent
        #[arg(long, default_value = "bfgs")]
        method: String,
        /// Maximum number of iterations
        #[arg(long, default_value_t = 200)]
        max_iter: usize,
        /// Force convergence threshold in eV/Å
        #[arg(long, default_value_t = 0.05)]
        force_tol: f64,
        /// Energy convergence threshold in eV
        #[arg(long, default_value_t = 1e-6)]
        energy_tol: f64,
        /// Maximum step size in Å
        #[arg(long, default_value_t = 0.2)]
        max_step: f64,
        /// JSON array of fixed atom indices; omit the flag to fix no atom
        #[arg(long)]
        fixed_atoms: Option<String>,
        /// Lattice length a in Å (enables periodic optimization)
        #[arg(long)]
        a: Option<f64>,
        /// Lattice length b in Å
        #[arg(long)]
        b: Option<f64>,
        /// Lattice length c in Å
        #[arg(long)]
        c: Option<f64>,
        /// Lattice angle alpha in degrees
        #[arg(long, default_value_t = 90.0)]
        alpha: f64,
        /// Lattice angle beta in degrees
        #[arg(long, default_value_t = 90.0)]
        beta: f64,
        /// Lattice angle gamma in degrees
        #[arg(long, default_value_t = 90.0)]
        gamma: f64,
    },
    /// Filter a conformer ensemble to keep only the diverse cluster centroids
    FilterDiverse {
        /// JSON array of conformers, each a flat xyz coordinate array
        conformers: String,
        /// RMSD clustering cutoff in Å
        #[arg(long, default_value_t = 1.0)]
        rmsd_cutoff: f64,
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
        #[arg(short, long)]
        charges: Option<String>,
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
        /// JSON array of real-space overlap blocks; omit the flag for a
        /// non-overlapping operator
        overlap_blocks: Option<String>,
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
        /// JSON array of backend names; omit the flag for all backends
        #[arg(long)]
        methods: Option<String>,
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
