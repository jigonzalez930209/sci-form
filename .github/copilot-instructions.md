# sci-form — Agent Instructions

**sci-form** is a Rust library (v0.13.0) for computational chemistry: 3D conformer generation, semi-empirical quantum chemistry (EHT, PM3, PM3(tm), GFN0/GFN1/GFN2-xTB), ab-initio (HF-3c, CISD, UHF/ROHF), neural network potentials (ANI-2x, ANI-TM), molecular properties, ML property prediction, 3D molecular descriptors (WHIM, RDF, GETAWAY), machine learning models (Random Forest, Gradient Boosting), stereochemistry (R/S, E/Z, helical, atropisomeric), solvation, ring perception, fingerprints, clustering, NMR/IR/UV-Vis spectroscopy, band structure, crystallographic materials (230 space groups, CIF import/export), and framework geometry optimization. Available as a Rust crate, Python package, WebAssembly/TypeScript module, and a CLI binary.

---

## Quick Reference — What the Library Does

| Module | Function |
|--------|----------|
| Conformer (ETKDG) | SMILES → 3D coordinates (ETKDGv2) |
| EHT | Semi-empirical Hückel MO theory, orbital grids, isosurfaces |
| EHT Band Structure | k-path band structure for periodic systems |
| EHT Gradients | Analytical gradients + geometry optimization |
| PM3 | NDDO PM3 semi-empirical SCF (heat of formation, HOMO/LUMO) |
| PM3(tm) | PM3 with transition metal parameters (Ti–Au) |
| GFN0-xTB | GFN0 tight-binding SCC (25 elements incl. transition metals) |
| GFN1-xTB | GFN1 tight-binding with shell-resolved charges + D3 dispersion |
| GFN2-xTB | GFN2 tight-binding with multipole electrostatics + D4 + XB |
| HF-3c | Minimal-basis Hartree-Fock with D3, gCP, SRB corrections |
| UHF/ROHF | Unrestricted & restricted open-shell Hartree-Fock SCF |
| CISD | Configuration Interaction Singles+Doubles (excited states) |
| AO→MO Transform | 4-index integral transform for post-HF methods |
| ANI-2x | Neural network potential for H, C, N, O, F, S, Cl |
| ANI-TM | ANI extended to 24 elements including transition metals |
| Charges | Gasteiger-Marsili partial charges |
| SASA | Solvent-accessible surface area (Shrake-Rupley) |
| Population | Mulliken & Löwdin population analysis (parallel, Z=1–86) |
| NPA/NBO | Natural Population Analysis & Natural Bond Orbital analysis |
| Dipole | Molecular dipole moment (Debye) |
| ESP | Electrostatic potential grid |
| DOS | Density of states with Gaussian smearing (multi-method) |
| Alignment | Kabsch SVD & quaternion RMSD alignment |
| Force Field | UFF and MMFF94 energy evaluation (+ GPU MMFF94) |
| ML Properties | LogP, MR, solubility, Lipinski Ro5, druglikeness |
| ML Descriptors | MW, Wiener index, Balaban J, FSP3, 17+ descriptors |
| 3D Descriptors | WHIM, RDF, GETAWAY molecular descriptors |
| ML Models | Decision Trees, Random Forest, Gradient Boosting, cross-validation |
| Materials | Unit cell, MOF framework assembly, 230 space groups, CIF import/export |
| Framework Opt | BFGS & steepest descent geometry optimization with PBC |
| Transport | Arrow columnar batch + Web Worker task splitting (parallelizable) |
| Stereochemistry | CIP priorities, R/S, E/Z, atropisomeric M/P, helical chirality |
| Solvation | Non-polar SASA solvation, Generalized Born (HCT) electrostatic |
| SSSR | Smallest Set of Smallest Rings (Horton's algorithm) |
| ECFP | Extended-Connectivity Fingerprints (Morgan), Tanimoto similarity |
| Clustering | Butina (Taylor-Butina) RMSD clustering, diversity filtering |
| NMR | Chemical shifts (HOSE codes), J-coupling (Karplus, 2J–5J), ensemble averaging |
| IR | Vibrational analysis, IR spectrum broadening, thermochemistry (RRHO) |
| Periodic Systems | Periodic molecular graphs with PBC, hapticity/metallocene detection |
| SMIRKS | Reaction transforms via atom-mapped reactant→product patterns (multi-component) |

---

## Rust

### Setup

```toml
[dependencies]
sci-form = "0.9"
# parallel batch (rayon):
sci-form = { version = "0.9", features = ["parallel"] }
# GPU acceleration:
sci-form = { version = "0.9", features = ["experimental-gpu"] }
```

### Core types

```rust
// ConformerResult — always returned, check .error
pub struct ConformerResult {
    pub smiles: String,
    pub num_atoms: usize,
    pub coords: Vec<f64>,                   // flat [x0,y0,z0, x1,y1,z1, ...]
    pub elements: Vec<u8>,                  // atomic numbers
    pub bonds: Vec<(usize, usize, String)>, // (a, b, "SINGLE"|"DOUBLE"|"TRIPLE"|"AROMATIC")
    pub error: Option<String>,
    pub time_ms: f64,
}

pub struct ConformerConfig { pub seed: u64, pub num_threads: usize }
```

### All public functions

```rust
// ── Conformer ──────────────────────────────────────────────
sci_form::embed(smiles: &str, seed: u64) -> ConformerResult
sci_form::embed_batch(smiles: &[&str], config: &ConformerConfig) -> Vec<ConformerResult>
sci_form::parse(smiles: &str) -> Result<graph::Molecule, String>

// ── Properties ─────────────────────────────────────────────
sci_form::compute_charges(smiles: &str) -> Result<ChargeResult, String>
sci_form::compute_charges_configured(smiles: &str, config: &GasteigerConfig) -> Result<ChargeResult, String>
sci_form::compute_sasa(elements: &[u8], coords_flat: &[f64], probe_radius: Option<f64>) -> Result<SasaResult, String>
sci_form::compute_population(elements: &[u8], positions: &[[f64;3]]) -> Result<PopulationResult, String>
sci_form::compute_dipole(elements: &[u8], positions: &[[f64;3]]) -> Result<DipoleResult, String>
sci_form::compute_esp(elements: &[u8], positions: &[[f64;3]], spacing: f64, padding: f64) -> Result<EspGrid, String>
sci_form::compute_dos(elements: &[u8], positions: &[[f64;3]], sigma: f64, e_min: f64, e_max: f64, n_points: usize) -> Result<DosResult, String>
sci_form::compute_rmsd(coords: &[f64], reference: &[f64]) -> f64

// ── Force Fields ───────────────────────────────────────────
sci_form::compute_uff_energy(smiles: &str, coords: &[f64]) -> Result<f64, String>
sci_form::compute_mmff94_energy(smiles: &str, coords: &[f64]) -> Result<f64, String>

// ── PM3 — NDDO semi-empirical SCF ─────────────────────────
sci_form::compute_pm3(elements: &[u8], positions: &[[f64;3]]) -> Result<pm3::Pm3Result, String>
// Pm3Result { orbital_energies, electronic_energy, nuclear_repulsion, total_energy,
//             heat_of_formation, n_basis, n_electrons, homo_energy, lumo_energy, gap,
//             mulliken_charges, scf_iterations, converged }
// PM3(tm) supports: H,C,N,O,F,P,S,Cl,Br,I + Al,Si,Ge + Ti,Cr,Mn,Fe,Co,Ni,Cu,Zn

// ── GFN-xTB family (GFN0, GFN1, GFN2) ────────────────────
sci_form::compute_xtb(elements: &[u8], positions: &[[f64;3]]) -> Result<xtb::XtbResult, String>
// GFN0: XtbResult { orbital_energies, electronic_energy, repulsive_energy, total_energy,
//        n_basis, n_electrons, homo_energy, lumo_energy, gap, mulliken_charges, scc_iterations, converged }

sci_form::xtb::gfn1::solve_gfn1(elements: &[u8], positions: &[[f64;3]]) -> Result<Gfn1Result, String>
// Gfn1Result { ..., dispersion_energy, shell_charges, ... }

sci_form::xtb::gfn2::solve_gfn2(elements: &[u8], positions: &[[f64;3]]) -> Result<Gfn2Result, String>
// Gfn2Result { ..., dispersion_energy, halogen_bond_energy, atomic_dipoles, atomic_quadrupoles, ... }

// ── ANI Neural Network Potentials ──────────────────────────
sci_form::compute_ani(elements: &[u8], positions: &[[f64;3]], config: &AniConfig, models: &HashMap<u8, FeedForwardNet>) -> Result<AniResult, String>
// AniResult { energy, forces: Vec<[f64;3]>, species, atomic_energies, aevs }

// ANI-TM — 24 elements including transition metals
sci_form::ani::ani_tm::compute_aevs_tm(elements: &[u8], positions: &[[f64;3]], params: &AevParams) -> Vec<Vec<f64>>
// Supported: H,C,N,O,F,Si,P,S,Cl,Ti,Cr,Mn,Fe,Co,Ni,Cu,Zn,Br,Ru,Pd,Ag,I,Pt,Au

// ── HF-3c — Minimal basis Hartree-Fock ────────────────────
sci_form::solve_hf3c(elements: &[u8], positions: &[[f64;3]], config: &HfConfig) -> Result<Hf3cResult, String>
// Hf3cResult { energy, hf_energy, nuclear_repulsion, d3_energy, gcp_energy, srb_energy,
//              orbital_energies, scf_iterations, converged, cis }

// ── UHF / ROHF — Open-shell Hartree-Fock ──────────────────
sci_form::compute_uhf(elements: &[u8], positions: &[[f64;3]], charge: i32, multiplicity: u32) -> Result<UhfResult, String>
sci_form::compute_rohf(elements: &[u8], positions: &[[f64;3]], charge: i32, multiplicity: u32) -> Result<UhfResult, String>
sci_form::compute_uhf_configured(elements: &[u8], positions: &[[f64;3]], charge: i32, multiplicity: u32, config: &UhfConfig) -> Result<UhfResult, String>
// UhfResult { total_energy, electronic_energy, nuclear_repulsion, alpha_orbital_energies,
//             beta_orbital_energies, n_alpha, n_beta, n_basis, scf_iterations, converged,
//             s2_expectation, spin_contamination, mulliken_charges }

// ── CISD — Excited states ─────────────────────────────────
sci_form::hf::cisd::compute_cisd(orbital_energies: &[f64], coefficients: &DMatrix<f64>,
    eris: &[f64], n_basis: usize, n_occupied: usize, n_states: usize) -> CisdResult
// CisdResult { excitations: Vec<CisdExcitation>, correlation_energy, n_csfs }
// CisdExcitation { energy, energy_ev, wavelength_nm, oscillator_strength, dominant_transition, character }

// ── AO→MO Integral Transform ─────────────────────────────
sci_form::hf::mo_transform::ao_to_mo_transform(eris_ao: &[f64], coefficients: &DMatrix<f64>,
    n_basis: usize) -> MoIntegrals
// MoIntegrals { get(p, q, r, s) -> f64, 4-fold symmetry: (pq|rs) = (qp|rs) = (pq|sr) = (rs|pq) }

// ── CIF Import/Export ─────────────────────────────────────
sci_form::parse_cif(cif_text: &str) -> Result<CifStructure, String>
sci_form::write_cif(structure: &CifStructure) -> String
// CifStructure { data_block, cell, cell_params, space_group_hm, space_group_number, atom_sites }

// ── EHT Band Structure (periodic systems) ─────────────────
sci_form::eht::band_structure::compute_band_structure(elements: &[u8], positions: &[[f64;3]],
    lattice: &[[f64;3]; 3], config: &BandStructureConfig, n_electrons: usize) -> Result<BandStructure, String>
// BandStructure { kpoints, bands, n_bands, n_kpoints, fermi_energy, direct_gap, indirect_gap }

// ── EHT Gradients & Geometry Optimization ─────────────────
sci_form::eht::gradients::compute_eht_gradient(elements: &[u8], positions: &[[f64;3]]) -> Result<EhtGradient, String>
// EhtGradient { gradients: Vec<[f64;3]>, rms_gradient, max_gradient, energy }

sci_form::eht::gradients::optimize_geometry_eht(elements: &[u8], initial_positions: &[[f64;3]],
    config: Option<EhtOptConfig>) -> Result<EhtOptResult, String>
// EhtOptResult { positions, energy, n_steps, converged, rms_gradient, energies }

// ── Multi-method DOS ──────────────────────────────────────
sci_form::dos::multi_method::compute_dos_multimethod(elements: &[u8], positions: &[[f64;3]],
    method: DosMethod, sigma: f64, e_min: f64, e_max: f64, n_points: usize) -> Result<MultiMethodDosResult, String>
// DosMethod: Eht, Pm3, Xtb, Gfn1, Gfn2, Hf3c
// MultiMethodDosResult { energies, total_dos, method, homo_energy, lumo_energy, gap, orbital_energies }

// ── NPA/NBO — Natural Population & Bond Orbital Analysis ──
sci_form::population::npa::compute_npa(elements: &[u8], overlap: &DMatrix<f64>,
    density: &DMatrix<f64>, basis_atom_map: &[usize]) -> Result<NpaResult, String>
// NpaResult { natural_charges, natural_config: Vec<NaturalConfig>, rydberg_population }

sci_form::population::npa::compute_nbo(elements: &[u8], overlap: &DMatrix<f64>,
    density: &DMatrix<f64>, basis_atom_map: &[usize], bonds: &[(usize, usize)]) -> Result<NboResult, String>
// NboResult { npa, bond_orbitals: Vec<NboBond>, lone_pairs: Vec<NboLonePair>, lewis_population_pct }

// ── Extended Basis Sets ───────────────────────────────────
sci_form::scf::extended_basis::build_basis_set(basis_type: BasisSetType,
    elements: &[u8], positions_bohr: &[[f64;3]]) -> BasisSet
// BasisSetType: Sto3g, Basis321g, Basis631g

// ── IR Spectroscopy ───────────────────────────────────────
sci_form::compute_ir_spectrum(elements: &[u8], positions: &[[f64;3]], method: ScientificMethod) -> Result<IrSpectrum, String>
sci_form::compute_vibrational_analysis(elements: &[u8], positions: &[[f64;3]], method: ScientificMethod) -> Result<VibrationalAnalysis, String>
sci_form::compute_ir_spectrum_broadened(analysis: &VibrationalAnalysis, gamma: f64, wn_min: f64, wn_max: f64, n_points: usize, broadening: &str) -> IrSpectrum

// IR Peak Assignment — functional group recognition
sci_form::ir::peak_assignment::assign_peaks(frequencies: &[f64], intensities: &[f64],
    elements: &[u8], mode_displacements: Option<&[Vec<[f64;3]>]>) -> AssignmentResult
// AssignmentResult { assignments: Vec<PeakAssignment>, functional_groups }

// ── NMR ───────────────────────────────────────────────────
sci_form::compute_ensemble_j_couplings(smiles: &str, conformer_coords: &[Vec<f64>],
    energies_kcal: &[f64], temperature_k: f64) -> Result<Vec<nmr::JCoupling>, String>

// ── ML Descriptors + Property Prediction ──────────────────
sci_form::compute_ml_descriptors(elements: &[u8], bonds: &[(usize, usize, String)],
    charges: &[f64], aromatic_atoms: &[bool]) -> ml::MolecularDescriptors
sci_form::predict_ml_properties(desc: &ml::MolecularDescriptors) -> ml::MlPropertyResult

// ── 3D Molecular Descriptors ──────────────────────────────
sci_form::ml::whim::compute_whim(elements: &[u8], positions: &[[f64;3]],
    weighting: WhimWeighting) -> WhimDescriptors
// WhimWeighting: Unit, Mass, Volume, Electronegativity, Polarizability
// WhimDescriptors { eigenvalues, lambda, nu, gamma, eta, total_spread, anisotropy }

sci_form::ml::rdf_descriptors::compute_rdf(elements: &[u8], positions: &[[f64;3]],
    charges: &[f64], beta: f64, r_max: f64, n_points: usize) -> RdfDescriptors
// RdfDescriptors { radii, rdf_unweighted, rdf_mass, rdf_electronegativity, rdf_charge }

sci_form::ml::getaway::compute_getaway(elements: &[u8], positions: &[[f64;3]],
    bonds: &[(usize, usize)], max_lag: usize) -> GetawayDescriptors
// GetawayDescriptors { leverages, hat_autocorrelation, r_autocorrelation, h_information, h_max, h_mean }

// ── ML Models (Random Forest, Gradient Boosting) ──────────
sci_form::ml::advanced_models::train_random_forest(features: &[Vec<f64>], targets: &[f64],
    config: &RandomForestConfig) -> RandomForest
sci_form::ml::advanced_models::train_gradient_boosting(features: &[Vec<f64>], targets: &[f64],
    config: &GradientBoostingConfig) -> GradientBoosting
// RandomForest { predict(&self, features) -> f64, predict_with_variance(..) -> (f64, f64) }
// GradientBoosting { predict(&self, features) -> f64 }

// ── Stereochemistry ───────────────────────────────────────
sci_form::analyze_stereo(smiles: &str, coords: &[f64]) -> Result<stereo::StereoAnalysis, String>
// StereoAnalysis { stereocenters, double_bonds, atropisomeric_axes, prochiral_centers,
//                  helical_chirality, n_stereocenters, n_double_bonds }
// HelicalChirality { backbone_atoms, configuration: "M"|"P", helix_type, torsion_sum }

// ── Solvation ─────────────────────────────────────────────
sci_form::compute_nonpolar_solvation(elements: &[u8], positions: &[[f64;3]], probe_radius: Option<f64>) -> NonPolarSolvation
sci_form::compute_gb_solvation(elements: &[u8], positions: &[[f64;3]], charges: &[f64],
    solvent_dielectric: Option<f64>, solute_dielectric: Option<f64>, probe_radius: Option<f64>) -> GbSolvation

// ── Ring Perception & Fingerprints ────────────────────────
sci_form::compute_sssr(smiles: &str) -> Result<rings::sssr::SssrResult, String>
sci_form::compute_ecfp(smiles: &str, radius: usize, n_bits: usize) -> Result<ECFPFingerprint, String>
sci_form::compute_tanimoto(fp1: &ECFPFingerprint, fp2: &ECFPFingerprint) -> f64

// ── Clustering ────────────────────────────────────────────
sci_form::butina_cluster(conformers: &[Vec<f64>], rmsd_cutoff: f64) -> clustering::ClusterResult
sci_form::compute_rmsd_matrix(conformers: &[Vec<f64>]) -> Vec<Vec<f64>>

// ── Materials & Crystallography ───────────────────────────
sci_form::create_unit_cell(a,b,c,alpha,beta,gamma: f64) -> materials::UnitCell
sci_form::assemble_framework(node: &Sbu, linker: &Sbu, topology: &Topology, cell: &UnitCell) -> CrystalStructure

// Space groups — all 230 ITC space groups
sci_form::materials::space_groups::space_group_by_number(number: u16) -> Option<SpaceGroup>
sci_form::materials::space_groups::space_group_by_symbol(symbol: &str) -> Option<SpaceGroup>
// SpaceGroup { number, symbol, crystal_system, lattice_type, point_group, operations }
// CrystalSystem: Triclinic, Monoclinic, Orthorhombic, Tetragonal, Trigonal, Hexagonal, Cubic
// SpaceGroup::generate_equivalent_positions(site: &[f64;3]) -> Vec<[f64;3]>

// Framework geometry optimization (BFGS / steepest descent with PBC)
sci_form::materials::geometry_opt::optimize_framework(elements: &[u8], initial_positions: &[[f64;3]],
    energy_force_fn: &EnergyForceFn, config: &FrameworkOptConfig) -> Result<FrameworkOptResult, String>
// OptMethod: SteepestDescent, Bfgs
// FrameworkOptResult { positions, energy, forces, max_force, n_iterations, converged, energy_history }

// ── Periodic Systems ──────────────────────────────────────
sci_form::periodic::build_periodic_molecule(elements: &[u8], frac_coords: &[[f64;3]],
    cell: &UnitCell, bond_tolerance: Option<f64>) -> PeriodicMolecule
sci_form::periodic::detect_hapticity(elements: &[u8], coords: &[[f64;3]],
    bonds: &[(usize, usize)]) -> HapticityAnalysis
// HapticityAnalysis { interactions: Vec<HapticInteraction>, is_metallocene, n_interactions }

// ── SMIRKS Reactions ──────────────────────────────────────
sci_form::smirks::parse_smirks(smirks: &str) -> Result<SmirksTransform, String>
sci_form::smirks::apply_smirks(smirks: &str, smiles: &str) -> Result<SmirksResult, String>

// ── Parallel Distance Bounds (requires feature "parallel") ─
sci_form::distgeom::parallel::compute_bounds_parallel(elements: &[u8],
    bonds: &[(usize, usize, String)]) -> ParallelBoundsResult

// ── GPU MMFF94 (requires feature "experimental-gpu") ──────
sci_form::gpu::mmff94_gpu::compute_mmff94_nonbonded_gpu(ctx: &GpuContext,
    coords: &[f64], charges: &[f64], vdw_radii: &[f64], vdw_epsilon: &[f64],
    exclusions_14: &[(usize, usize)]) -> Result<Mmff94GpuResult, String>
```

### Coordinate convention

`coords` is always a **flat `Vec<f64>`**: `[x0, y0, z0, x1, y1, z1, ...]` of length `num_atoms × 3`.

Convert to `&[[f64;3]]` when needed:
```rust
let pos: Vec<[f64;3]> = conf.coords.chunks(3).map(|c| [c[0],c[1],c[2]]).collect();
```

### Typical pipeline

```rust
let conf = sci_form::embed("c1ccccc1", 42);
assert!(conf.error.is_none());
let pos: Vec<[f64;3]> = conf.coords.chunks(3).map(|c| [c[0],c[1],c[2]]).collect();

// Population analysis
let pop = sci_form::compute_population(&conf.elements, &pos).unwrap();
println!("Gap: {:.3} eV", pop.homo_lumo_gap);

// Dipole moment
let dipole = sci_form::compute_dipole(&conf.elements, &pos).unwrap();
println!("|μ| = {:.3} D", dipole.magnitude);

// PM3 semi-empirical SCF (now supports transition metals via PM3(tm))
let pm3 = sci_form::compute_pm3(&conf.elements, &pos).unwrap();
println!("PM3 HOF: {:.2} kcal/mol, gap: {:.3} eV", pm3.heat_of_formation, pm3.gap);

// GFN-xTB family
let xtb = sci_form::compute_xtb(&conf.elements, &pos).unwrap(); // GFN0
let gfn1 = sci_form::xtb::gfn1::solve_gfn1(&conf.elements, &pos).unwrap();
let gfn2 = sci_form::xtb::gfn2::solve_gfn2(&conf.elements, &pos).unwrap();
println!("GFN2 gap: {:.3} eV, D4: {:.4} eV", gfn2.gap, gfn2.dispersion_energy);

// ML properties (no 3D needed — only SMILES topology)
let desc = sci_form::compute_ml_descriptors(&conf.elements, &conf.bonds, &[], &[]);
let props = sci_form::predict_ml_properties(&desc);
println!("LogP: {:.2}, Druglikeness: {:.3}", props.logp, props.druglikeness);

// 3D Descriptors (WHIM, RDF, GETAWAY)
use sci_form::ml::whim::{compute_whim, WhimWeighting};
let whim = compute_whim(&conf.elements, &pos, WhimWeighting::Mass);
println!("WHIM spread: {:.3}, anisotropy: {:.3}", whim.total_spread, whim.anisotropy);

// Stereochemistry (R/S, E/Z, helical, atropisomeric)
let stereo = sci_form::analyze_stereo("C(F)(Cl)(Br)I", &conf.coords).unwrap();
println!("Stereocenters: {}, helical: {:?}", stereo.n_stereocenters, stereo.helical_chirality);

// Ring perception + fingerprints
let fp1 = sci_form::compute_ecfp("c1ccccc1", 2, 2048).unwrap();
let fp2 = sci_form::compute_ecfp("Cc1ccccc1", 2, 2048).unwrap();
println!("Tanimoto: {:.3}", sci_form::compute_tanimoto(&fp1, &fp2));

// Solvation
let solv = sci_form::compute_nonpolar_solvation(&conf.elements, &pos, None);
println!("Non-polar solvation: {:.2} kcal/mol", solv.energy_kcal_mol);

// Multi-method DOS
use sci_form::dos::multi_method::{compute_dos_multimethod, DosMethod};
let dos = compute_dos_multimethod(&conf.elements, &pos, DosMethod::Pm3, 0.3, -30.0, 5.0, 500).unwrap();
println!("PM3 DOS gap: {:.3} eV", dos.gap);

// Space groups (230 ITC groups)
use sci_form::materials::space_groups::space_group_by_number;
let sg = space_group_by_number(225).unwrap(); // Fm-3m (NaCl-type)
println!("SG {}: {} ({})", sg.number, sg.symbol, sg.crystal_system);

// ML Models
use sci_form::ml::advanced_models::*;
let rf_config = RandomForestConfig { n_trees: 100, tree_config: TreeConfig::default(),
    sample_fraction: 0.8, seed: 42 };
let rf = train_random_forest(&features, &targets, &rf_config);
let prediction = rf.predict(&new_sample);

// Parallel batch (requires features = ["parallel"])
let config = sci_form::ConformerConfig { seed: 42, num_threads: 0 };
let results = sci_form::embed_batch(&["CCO", "c1ccccc1", "CC(=O)O"], &config);
```

---

## Python

### Setup

```bash
pip install sciforma      # PyPI package name
import sci_form           # import name
```

### All public functions

```python
from sci_form import (
    embed, embed_batch, parse,
    charges, sasa, eht_calculate, eht_orbital_mesh,
    population, dipole, dos, rmsd,
    uff_energy, mmff94_energy,
    pm3_calculate, xtb_calculate,
    ml_descriptors, ml_predict,
    unit_cell, assemble_framework,
    pack_conformers, split_worker_tasks, estimate_workers,
    stereo_analysis, nonpolar_solvation, gb_solvation,
    sssr, ecfp, tanimoto, butina_cluster, filter_diverse,
    hf3c_calculate, ani_calculate, compute_esp,
    stda_uvvis, vibrational_analysis, ir_spectrum,
    nmr_shift, j_coupling, nmr_spectrum, ensemble_j_couplings,
    fukui_descriptors, reactivity_ranking, empirical_pka,
    frontier_descriptors, bond_orders,
)

# ── Conformer ──────────────────────────────────────────────
embed(smiles: str, seed: int = 42) -> ConformerResult
embed_batch(smiles_list: list[str], seed: int = 42, num_threads: int = 0) -> list[ConformerResult]
parse(smiles: str) -> dict   # {"num_atoms", "atoms", "bonds"}

# ── Properties ─────────────────────────────────────────────
charges(smiles: str) -> ChargeResult          # .charges, .total_charge
sasa(elements, coords, probe_radius=1.4) -> SasaResult   # .total_sasa Å², .atom_sasa
population(elements, coords) -> PopulationResult         # .mulliken_charges, .lowdin_charges, .homo_energy, .lumo_energy
dipole(elements, coords) -> DipoleResult                 # .magnitude Debye, .vector [x,y,z]
dos(elements, coords, sigma=0.3, e_min=-30.0, e_max=5.0, n_points=500) -> DosResult
rmsd(coords, reference) -> AlignmentResult              # .rmsd Å, .aligned_coords
compute_esp(elements, coords, spacing=0.5, padding=3.0) -> EspGridPy

# ── Force Fields ───────────────────────────────────────────
uff_energy(smiles, coords) -> float                     # kcal/mol
uff_energy_with_aromatic_heuristics(smiles, coords) -> UffHeuristicEnergyPy
mmff94_energy(smiles, coords) -> float                  # kcal/mol

# ── PM3 — NDDO semi-empirical SCF (+ PM3(tm) transition metals) ─
pm3_calculate(elements: list[int], coords: list[float]) -> Pm3Result
# Pm3Result: .orbital_energies, .electronic_energy, .nuclear_repulsion, .total_energy,
#   .heat_of_formation, .n_basis, .n_electrons, .homo_energy, .lumo_energy,
#   .gap, .mulliken_charges, .scf_iterations, .converged

# ── GFN-xTB — GFN0 tight-binding (25 elements incl. transition metals)
xtb_calculate(elements: list[int], coords: list[float]) -> XtbResult
# XtbResult: .orbital_energies, .electronic_energy, .repulsive_energy, .total_energy,
#   .n_basis, .n_electrons, .homo_energy, .lumo_energy,
#   .gap, .mulliken_charges, .scc_iterations, .converged

# ── HF-3c + ANI ───────────────────────────────────────────
hf3c_calculate(elements, coords, max_scf_iter=300, n_cis_states=5, corrections=True) -> Hf3cResultPy
ani_calculate(elements, coords) -> AniResultPy

# ── EHT ────────────────────────────────────────────────────
eht_calculate(elements, coords, k=1.75) -> EhtResult
eht_orbital_mesh(elements, coords, mo_index, spacing=0.2, isovalue=0.02) -> dict

# ── ML Descriptors + Property Prediction ───────────────────
ml_descriptors(smiles: str) -> MolecularDescriptors
ml_predict(smiles: str) -> MlPropertyResult
# MlPropertyResult: .logp, .molar_refractivity, .log_solubility, .druglikeness,
#   .lipinski_violations, .lipinski_passes

# ── Spectroscopy ───────────────────────────────────────────
stda_uvvis(elements, coords, sigma=0.3, e_min=1.0, e_max=8.0, n_points=500, broadening="gaussian") -> StdaUvVisSpectrumPy
vibrational_analysis(elements, coords, method="pm3") -> VibrationalAnalysisPy
ir_spectrum(elements, coords, method="pm3", gamma=30.0, wn_min=500.0, wn_max=3500.0, n_points=1000, broadening="lorentzian") -> IrSpectrumPy

# ── NMR ────────────────────────────────────────────────────
nmr_shift(smiles, coords=[]) -> NmrShiftResultPy
j_coupling(smiles, coords=[]) -> list[JCouplingPy]
nmr_spectrum(smiles, coords=[], field_strength_mhz=400.0) -> NmrSpectrumPy
ensemble_j_couplings(smiles, conformer_coords, energies, temperature=298.15) -> list[JCouplingPy]

# ── Population & Bond Orders ──────────────────────────────
bond_orders(elements, coords) -> BondOrderResultPy
frontier_descriptors(elements, coords) -> FrontierDescriptorsPy

# ── Reactivity ────────────────────────────────────────────
fukui_descriptors(elements, coords) -> FukuiDescriptorsPy
reactivity_ranking(elements, coords) -> ReactivityRankingPy
empirical_pka(smiles) -> EmpiricalPkaResultPy

# ── Stereochemistry — CIP priorities, R/S, E/Z, helical, atropisomeric
stereo_analysis(smiles: str, coords: list[float] = []) -> StereoAnalysisPy
# StereoAnalysisPy: .stereocenters, .double_bonds, .atropisomeric_axes,
#   .prochiral_centers, .helical_chirality, .n_stereocenters, .n_double_bonds

# ── Solvation — Non-polar SASA + Generalized Born ────────
nonpolar_solvation(elements, coords, probe_radius=1.4) -> NonPolarSolvationPy
gb_solvation(elements, coords, charges, solvent_dielectric=78.5, solute_dielectric=1.0, probe_radius=1.4) -> GbSolvationPy

# ── Ring perception — SSSR ─────────────────────────────────
sssr(smiles: str) -> SssrResultPy

# ── ECFP fingerprints + Tanimoto ───────────────────────────
ecfp(smiles: str, radius=2, n_bits=2048) -> ECFPFingerprintPy
tanimoto(fp1, fp2) -> float

# ── Clustering — Butina RMSD ──────────────────────────────
butina_cluster(conformers: list[list[float]], rmsd_cutoff=1.0) -> ClusterResultPy
rmsd_matrix(conformers) -> list[list[float]]
filter_diverse(conformers, rmsd_cutoff=1.0) -> list[int]

# ── Materials ──────────────────────────────────────────────
unit_cell(a, b, c, alpha=90.0, beta=90.0, gamma=90.0) -> UnitCellPy
assemble_framework(topology="pcu", metal=30, geometry="octahedral", lattice_a=10.0, supercell=1) -> CrystalStructurePy

# ── Transport ──────────────────────────────────────────────
pack_conformers(results) -> RecordBatchPy
split_worker_tasks(smiles, n_workers=4, seed=42) -> list[dict]
estimate_workers(n_items, max_workers=8) -> int
```

### ConformerResult methods

```python
result.is_ok()          # True if no error
result.get_positions()  # list of (x, y, z) tuples
result.smiles           # input SMILES
result.num_atoms        # int
result.coords           # flat list[float]
result.elements         # list[int] of atomic numbers
result.time_ms          # float
```

### Typical pipeline

```python
from sci_form import embed, population, dipole, pm3_calculate, xtb_calculate, ml_predict

conf = embed("CCO", seed=42)
pop = population(conf.elements, conf.coords)
print(f"HOMO: {pop.homo_energy:.3f} eV, Gap: {pop.homo_lumo_gap:.3f} eV")

d = dipole(conf.elements, conf.coords)
print(f"|μ| = {d.magnitude:.3f} D")

# PM3 semi-empirical (supports TMs via PM3(tm))
pm3 = pm3_calculate(conf.elements, conf.coords)
print(f"PM3 HOF: {pm3.heat_of_formation:.2f} kcal/mol")

# GFN-xTB
xtb = xtb_calculate(conf.elements, conf.coords)
print(f"xTB gap: {xtb.gap:.3f} eV, converged: {xtb.converged}")

# ML properties (no 3D needed)
props = ml_predict("CCO")
print(f"LogP: {props.logp:.2f}, passes Lipinski: {props.lipinski_passes}")

# Stereochemistry (R/S, E/Z, helical, atropisomeric)
from sci_form import stereo_analysis
stereo = stereo_analysis("C(F)(Cl)(Br)I")
print(f"Stereocenters: {stereo.n_stereocenters}")

# Fingerprints + similarity
from sci_form import ecfp, tanimoto
fp1 = ecfp("c1ccccc1", radius=2, n_bits=2048)
fp2 = ecfp("Cc1ccccc1", radius=2, n_bits=2048)
print(f"Tanimoto: {tanimoto(fp1, fp2):.3f}")

# Solvation
from sci_form import nonpolar_solvation
solv = nonpolar_solvation(conf.elements, conf.coords)
print(f"Non-polar solvation: {solv.energy_kcal_mol:.2f} kcal/mol")

# Spectroscopy
from sci_form import ir_spectrum, stda_uvvis
ir = ir_spectrum(conf.elements, conf.coords, method="pm3")
uv = stda_uvvis(conf.elements, conf.coords)

# Parallel batch
from sci_form import embed_batch
results = embed_batch(["CCO", "c1ccccc1", "CC(=O)O"], seed=42, num_threads=4)
```

---

## TypeScript / JavaScript (WASM)

### Setup

```bash
npm install sci-form-wasm
```

**WASM cannot pass complex objects across the boundary.** All functions accept/return:
- `string` — JSON-encoded input/output
- `Float64Array` / `Float32Array` — typed arrays for performance-critical paths (avoid JSON overhead)

### Initialization (browser / Deno / Bun)

```typescript
import init, { embed, /* ... */ } from 'sci-form-wasm';
await init();
```

### Node.js CommonJS

```javascript
const sci = require('sci-form-wasm');
// no async init needed
```

### All exported functions

```typescript
// Conformer
embed(smiles: string, seed?: number): string                     // JSON ConformerResult
embed_coords(smiles: string, seed?: number): string              // JSON {coords, num_atoms}
embed_coords_typed(smiles: string, seed: number): Float64Array  // NO JSON — fastest
embed_batch(smiles_newline_sep: string, seed?: number): string   // JSON array
parse_smiles(smiles: string): string                             // JSON {num_atoms, num_bonds}

// EHT (elements/coords_flat are ALWAYS JSON strings)
eht_calculate(elements: string, coords_flat: string, k: number): string   // JSON EhtResult
eht_orbital_mesh(elements: string, coords_flat: string, mo_index: number, spacing: number, isovalue: number): string  // JSON OrbitalMesh
eht_orbital_grid_typed(elements: string, coords_flat: string, mo_index: number, spacing: number): Float32Array

// ESP
compute_esp_grid_typed(elements: string, coords_flat: string, spacing: number, padding: number): Float64Array
compute_esp_grid_info(elements: string, coords_flat: string, spacing: number, padding: number): string
// info JSON: {origin:[x,y,z], spacing:f, dims:[nx,ny,nz]}

// Properties
compute_charges(smiles: string): string          // JSON {charges, iterations, total_charge}
compute_sasa(elements: string, coords_flat: string, probe_radius: number): string
compute_population(elements: string, coords_flat: string): string
compute_bond_orders(elements: string, coords_flat: string): string
// JSON BondOrderResult: {atom_pairs, distances, wiberg, mayer, wiberg_valence, mayer_valence}
compute_dipole(elements: string, coords_flat: string): string
compute_dos(elements: string, coords_flat: string, sigma: number, e_min: number, e_max: number, n_points: number): string
compute_rmsd(coords: string, reference: string): string  // JSON {"rmsd": 0.034}
compute_uff_energy(smiles: string, coords: string): string
compute_mmff94_energy(smiles: string, coords: string): string  // JSON {"energy": 12.3}

// PM3 — NDDO semi-empirical SCF
compute_pm3(elements: string, coords_flat: string): string
// JSON Pm3Result: {orbital_energies, electronic_energy, nuclear_repulsion, total_energy,
//   heat_of_formation, n_basis, n_electrons, homo_energy, lumo_energy, gap,
//   mulliken_charges, scf_iterations, converged}

// GFN-xTB — GFN0 tight-binding (25 elements incl. transition metals)
compute_xtb(elements: string, coords_flat: string): string
// JSON XtbResult: {orbital_energies, electronic_energy, repulsive_energy, total_energy,
//   n_basis, n_electrons, homo_energy, lumo_energy, gap,
//   mulliken_charges, scc_iterations, converged}

// ML Properties + Descriptors (topology only, no 3D needed)
compute_ml_properties(smiles: string): string
// JSON: {logp, molar_refractivity, log_solubility, lipinski_violations, lipinski_passes, druglikeness}

compute_molecular_descriptors(smiles: string): string
// JSON MolecularDescriptors: {molecular_weight, n_heavy_atoms, n_hydrogens, n_bonds,
//   n_rotatable_bonds, n_hbd, n_hba, fsp3, wiener_index, n_rings, n_aromatic,
//   balaban_j, sum_electronegativity, sum_polarizability, ...}

// Spectroscopy — UV-Vis, IR, NMR
compute_stda_uvvis(elements: string, coords_flat: string, sigma: number, e_min: number, e_max: number, n_points: number, broadening: string): string
// JSON: {energies_ev, intensities, peaks, n_transitions}

compute_vibrational_analysis(elements: string, coords_flat: string, method: string, step_size: number): string
// JSON: {n_atoms, modes, zpve_hartree, zpve_kcal_mol, thermal_energy, entropy, ...}

compute_ir_spectrum(analysis_json: string, gamma: number, wn_min: number, wn_max: number, n_points: number): string
// JSON: {wavenumbers, intensities, peaks}

predict_nmr_shifts(smiles: string): string
// JSON: {h_shifts, c_shifts, f_shifts, p_shifts, n_shifts, ...}

predict_nmr_couplings(smiles: string, coords_flat: string): string
// JSON: [{h1_index, h2_index, j_hz, n_bonds, coupling_type}]

compute_nmr_spectrum(smiles: string, nucleus: string, gamma: number, ppm_min: number, ppm_max: number, n_points: number): string
// JSON: {ppm_axis, intensities, peaks, nucleus, gamma}

// Reactivity
compute_frontier_descriptors(elements: string, coords_flat: string): string
// JSON: {homo_atom_contributions, lumo_atom_contributions, dual_descriptor, homo_energy, lumo_energy, gap}

compute_fukui_descriptors(elements: string, coords_flat: string): string
// JSON: {condensed_atom_indices, condensed_f_plus, condensed_f_minus, condensed_f_radical, ...}

compute_reactivity_ranking(elements: string, coords_flat: string): string
// JSON: {nucleophilic_attack_sites, electrophilic_attack_sites, radical_attack_sites}

compute_empirical_pka(smiles: string): string
// JSON: {acidic_sites, basic_sites}

// Materials
create_unit_cell(a,b,c,alpha,beta,gamma: number): string
assemble_framework(topology: string, metal: number, geometry: string, lattice_a: number, supercell: number): string

// Transport
pack_batch_arrow(results_json: string): string
split_worker_tasks(smiles_json: string, n_workers: number, seed: number): string
estimate_workers(n_items: number, max_workers: number): number

// Stereochemistry — CIP priorities, R/S, E/Z
analyze_stereo(smiles: string, coords_flat: string): string
// JSON StereoAnalysis: {stereocenters, double_bonds, n_stereocenters, n_double_bonds}

// Solvation — Non-polar SASA + Generalized Born
compute_nonpolar_solvation(elements: string, coords_flat: string, probe_radius: number): string
// JSON NonPolarSolvation: {energy_kcal_mol, atom_contributions, atom_sasa, total_sasa}

compute_gb_solvation(elements: string, coords_flat: string, charges: string, solvent_dielectric: number, solute_dielectric: number, probe_radius: number): string
// JSON GbSolvation: {electrostatic_energy_kcal_mol, nonpolar_energy_kcal_mol, total_energy_kcal_mol, born_radii}

// Ring perception — SSSR
compute_sssr(smiles: string): string
// JSON SssrResult: {rings, atom_ring_count, atom_ring_sizes, ring_size_histogram}

// ECFP fingerprints + Tanimoto
compute_ecfp(smiles: string, radius: number, n_bits: number): string
// JSON ECFPFingerprint: {n_bits, on_bits, radius, raw_features}

compute_tanimoto(fp1_json: string, fp2_json: string): string
// JSON {"tanimoto": 0.85}

// Clustering — Butina RMSD
butina_cluster(conformers_json: string, rmsd_cutoff: number): string
// JSON ClusterResult: {n_clusters, assignments, centroid_indices, cluster_sizes}

compute_rmsd_matrix(conformers_json: string): string
// JSON 2D array
```

### JSON ↔ typed array pattern (critical for performance)

```typescript
// Convert embed result to JSON arrays for EHT/ESP:
const result = JSON.parse(embed('CCO', 42));
const elements = JSON.stringify(Array.from(result.elements));  // "[8,6,6,1,1,1,1,1,1]"
const coords = JSON.stringify(result.coords);                  // "[x0,y0,z0,...]"

// Or use typed array path (no JSON):
const coordsTyped: Float64Array = embed_coords_typed('CCO', 42);
const coordsJson = JSON.stringify(Array.from(coordsTyped));
```

### PM3 / xTB energy calculation

```typescript
import init, { embed, compute_pm3, compute_xtb, compute_ml_properties } from 'sci-form-wasm';
await init();

const result = JSON.parse(embed('CCO', 42));
const elements = JSON.stringify(result.elements);
const coords   = JSON.stringify(result.coords);

// PM3
const pm3 = JSON.parse(compute_pm3(elements, coords));
console.log(`PM3 HOF: ${pm3.heat_of_formation.toFixed(2)} kcal/mol, gap: ${pm3.gap.toFixed(3)} eV`);

// GFN-xTB
const xtb = JSON.parse(compute_xtb(elements, coords));
console.log(`xTB energy: ${xtb.total_energy.toFixed(4)} eV, converged: ${xtb.converged}`);

// ML properties (no 3D needed)
const props = JSON.parse(compute_ml_properties('CCO'));
console.log(`LogP: ${props.logp.toFixed(2)}, Lipinski: ${props.lipinski_passes}`);
```

### Three.js orbital visualization

```typescript
const eht = JSON.parse(eht_calculate(elements, coords, 0.0));
const mesh = JSON.parse(eht_orbital_mesh(elements, coords, eht.homo_index, 0.2, 0.02));

const geometry = new THREE.BufferGeometry();
geometry.setAttribute('position', new THREE.Float32BufferAttribute(mesh.vertices, 3));
geometry.setAttribute('normal',   new THREE.Float32BufferAttribute(mesh.normals, 3));
geometry.setIndex(mesh.indices);
```

### Web Worker pattern

```typescript
// Main thread
const nWorkers = estimate_workers(smiles.length, navigator.hardwareConcurrency ?? 4);
const tasks = JSON.parse(split_worker_tasks(JSON.stringify(smiles), nWorkers, 42));
// tasks[i] = { smiles: [...], seed: 42 }
const workers = tasks.map(task => {
  const w = new Worker(new URL('./worker.ts', import.meta.url), { type: 'module' });
  w.postMessage(task);
  return w;
});

// worker.ts
import init, { embed_batch } from 'sci-form-wasm';
await init();
self.onmessage = ({ data }) => {
  const results = JSON.parse(embed_batch(data.smiles.join('\n'), data.seed));
  self.postMessage(results);
};
```

---

## CLI

### Installation

```bash
cargo install sci-form-cli
# or download pre-built binary from GitHub Releases
```

### Commands

```bash
# Generate 3D conformer
sci-form embed "CCO"                        # JSON output
sci-form embed "c1ccccc1" -f xyz            # XYZ format
sci-form embed "CC(=O)O" -f sdf -s 123     # SDF with custom seed

# Batch from file (one SMILES per line)
sci-form batch molecules.smi
sci-form batch molecules.smi -f xyz -t 8   # 8 threads
cat molecules.smi | sci-form batch /dev/stdin

# Parse SMILES (no 3D)
sci-form parse "c1ccccc1"

# Gasteiger charges
sci-form charges "CCO"

# EHT calculation
sci-form eht "[6,1,1,1,1]" "[0,0,0,1,0,0,0,1,0,0,0,1]" --k 1.75

# Population analysis (Mulliken + Löwdin)
sci-form population "[6,8,1,1,1,1,1,1,1]" "[...]"

# Dipole moment
sci-form dipole "[6,8,1,1,1,1,1,1,1]" "[...]"

# Electrostatic potential grid
sci-form esp "[6,8,1,1,1,1,1,1,1]" "[...]" --spacing 0.5 --padding 3.0

# Density of states
sci-form dos "[6,8,1,1,1,1,1,1,1]" "[...]" --sigma 0.3 --e-min -30 --e-max 5 -n 500

# RMSD alignment
sci-form rmsd "[...]" "[...]"

# UFF force field energy
sci-form uff "CCO" "[...]"

# ANI ML potential
sci-form ani "[6,8,1,1,1,1,1,1,1]" "[...]"

# HF-3c quantum calculation
sci-form hf3c "[6,8,1,1,1,1,1,1,1]" "[...]"

# Stereochemistry analysis
sci-form stereo "C(F)(Cl)(Br)I"
sci-form stereo "C/C=C/C" --coords "[...]"

# Ring detection (SSSR)
sci-form sssr "c1ccccc1"

# ECFP fingerprint
sci-form ecfp "CCO" --radius 2 --n-bits 2048

# Tanimoto similarity
sci-form tanimoto "c1ccccc1" "Cc1ccccc1"

# Solvation energy (non-polar + optional GB)
sci-form solvation "[8,1,1]" "[0,0,0,0.757,0.586,0,-0.757,0.586,0]"
sci-form solvation "[8,1,1]" "[...]" --charges "[-0.8,0.4,0.4]"

# Unit cell
sci-form cell --a 5.43 --b 5.43 --c 5.43 --alpha 90 --beta 90 --gamma 90

# Framework assembly
sci-form assemble --topology pcu --metal 30 --geometry octahedral --a 10.0

# Build info
sci-form info
```

### Output formats: `json` (default) | `xyz` | `sdf`

### Exit codes: `0` success · `1` invalid SMILES / embedding failure · `2` I/O error

### Pipeline examples

```bash
# Filter only successful embeddings
sci-form batch molecules.smi | jq -c 'select(.error == null)'

# Batch statistics
sci-form batch molecules.smi | jq -s '{total:length, success:[.[]|select(.error==null)]|length}'

# Parallel with GNU Parallel
cat huge.smi | parallel -j8 --pipe -N100 sci-form batch /dev/stdin > results.json
```

---

## Common Pitfalls

- **`coords` is always flat** `[x0,y0,z0, x1,y1,z1, ...]` — never nested. Length = `num_atoms × 3`.
- **`elements` are atomic numbers** (u8/int): 1=H, 6=C, 7=N, 8=O, 9=F, 15=P, 16=S, 17=Cl, 35=Br, 53=I.
- **WASM: always JSON.stringify arrays** before passing to WASM functions (`elements` and `coords_flat` parameters).
- **`embed` never panics** — always check `result.error` (Rust) or `result.is_ok()` (Python).
- **EHT / PM3 / xTB all require 3D coordinates** — run `embed` first, then pass `elements` + `coords`.
- **`ml_predict` / `compute_ml_properties` do NOT need 3D** — they accept SMILES directly.
- **ESP/DOS are slow for large molecules** — prefer `spacing ≥ 0.3` Å for interactive use.
- **GFN-xTB transition metals**: supported elements are Ti(22), Cr(24), Mn(25), Fe(26), Co(27), Ni(28), Cu(29), Zn(30), Ru(44), Pd(46), Ag(47), Pt(78), Au(79). EHT also supports all of these.
- **PM3 `heat_of_formation`** is in kcal/mol; all energies (`electronic_energy`, `total_energy`) are in eV.
- **xTB field names**: `repulsive_energy` (not `repulsion_energy`), `scc_iterations` (not `scf_iterations`).
- **Parallel embed**: use `embed_batch` with `features = ["parallel"]` or `num_threads > 0` in Python; parallelize with `split_worker_tasks` in WASM/workers.
- **`topology` strings**: `"pcu"` (cubic), `"dia"` (diamond), `"sql"` (square lattice).
- **`geometry` strings**: `"linear"`, `"trigonal"`, `"tetrahedral"`, `"square_planar"`, `"octahedral"`.

---

## Units

| Quantity | Unit |
|----------|------|
| Coordinates | Å (ångströms) |
| Energies (EHT, PM3, xTB, DOS, HOMO/LUMO) | eV |
| PM3 heat of formation | kcal/mol |
| Force field energy (UFF, MMFF94) | kcal/mol |
| Solvation energy | kcal/mol |
| Dipole moment | Debye |
| SASA | Å² |
| Lattice parameters a,b,c | Å |
| Lattice angles α,β,γ | degrees |
| IR frequencies | cm⁻¹ |
| NMR shifts | ppm |
| J-coupling | Hz |
| UV-Vis energies | eV |

---

## SMILES Examples

```
CCO               ethanol
c1ccccc1          benzene
CC(=O)O           acetic acid
CC(=O)Nc1ccc(O)cc1  acetaminophen
c1ccc(cc1)C(=O)O  benzoic acid
[Zn]              zinc atom
C1CCCCC1          cyclohexane
C=C               ethylene
C#N               hydrogen cyanide
```
