# sci-form — Agent Instructions

**sci-form** is a Rust library (v0.7.0) for computational chemistry: 3D conformer generation, semi-empirical quantum chemistry (EHT, PM3, GFN-xTB), molecular properties, ML property prediction, stereochemistry, solvation, ring perception, fingerprints, clustering, and crystallographic materials. It is available as a Rust crate, Python package, WebAssembly/TypeScript module, and a CLI binary.

---

## Quick Reference — What the Library Does

| Module | Function |
|--------|----------|
| Conformer (ETKDG) | SMILES → 3D coordinates (ETKDGv2) |
| EHT | Semi-empirical Hückel MO theory, orbital grids, isosurfaces |
| PM3 | NDDO PM3 semi-empirical SCF (heat of formation, HOMO/LUMO) |
| GFN-xTB | GFN0 tight-binding SCC (25 elements incl. transition metals) |
| Charges | Gasteiger-Marsili partial charges |
| SASA | Solvent-accessible surface area (Shrake-Rupley) |
| Population | Mulliken & Löwdin population analysis |
| Dipole | Molecular dipole moment (Debye) |
| ESP | Electrostatic potential grid |
| DOS | Density of states with Gaussian smearing |
| Alignment | Kabsch SVD & quaternion RMSD alignment |
| Force Field | UFF and MMFF94 energy evaluation |
| ML Properties | LogP, MR, solubility, Lipinski Ro5, druglikeness |
| ML Descriptors | MW, Wiener index, Balaban J, FSP3, 17 descriptors |
| Materials | Unit cell, MOF framework assembly |
| Transport | Arrow columnar batch + Web Worker task splitting (parallelizable) |
| Stereochemistry | CIP priorities, R/S chirality, E/Z double bond configuration |
| Solvation | Non-polar SASA solvation, Generalized Born (HCT) electrostatic |
| SSSR | Smallest Set of Smallest Rings (Horton's algorithm) |
| ECFP | Extended-Connectivity Fingerprints (Morgan), Tanimoto similarity |
| Clustering | Butina (Taylor-Butina) RMSD clustering, diversity filtering |
| NMR | Chemical shifts (HOSE codes), J-coupling (Karplus), ensemble averaging |
| IR | Vibrational analysis, IR spectrum broadening, thermochemistry (RRHO) |

---

## Rust

### Setup

```toml
[dependencies]
sci-form = "0.7"
# parallel batch (rayon):
sci-form = { version = "0.7", features = ["parallel"] }
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
// Conformer
sci_form::embed(smiles: &str, seed: u64) -> ConformerResult
sci_form::embed_batch(smiles: &[&str], config: &ConformerConfig) -> Vec<ConformerResult>
sci_form::parse(smiles: &str) -> Result<graph::Molecule, String>

// Properties
sci_form::compute_charges(smiles: &str) -> Result<ChargeResult, String>
// ChargeResult { charges: Vec<f64>, total_charge: f64 }

sci_form::compute_sasa(elements: &[u8], coords_flat: &[f64], probe_radius: Option<f64>) -> Result<SasaResult, String>
// SasaResult { total_sasa: f64 /*Å²*/, per_atom: Vec<f64> }

sci_form::compute_population(elements: &[u8], positions: &[[f64;3]]) -> Result<PopulationResult, String>
// PopulationResult { mulliken_charges, lowdin_charges, homo_energy, lumo_energy, homo_lumo_gap }

sci_form::compute_dipole(elements: &[u8], positions: &[[f64;3]]) -> Result<DipoleResult, String>
// DipoleResult { vector: [f64;3] /*Debye*/, magnitude: f64 }

sci_form::compute_esp(elements: &[u8], positions: &[[f64;3]], spacing: f64, padding: f64) -> Result<EspGrid, String>
// EspGrid { values: Vec<f64>, origin: [f64;3], spacing: f64, nx, ny, nz: usize }

sci_form::compute_dos(elements: &[u8], positions: &[[f64;3]], sigma: f64, e_min: f64, e_max: f64, n_points: usize) -> Result<DosResult, String>
// DosResult { energies, dos, pdos, orbital_energies, homo_lumo_gap }

sci_form::compute_rmsd(coords: &[f64], reference: &[f64]) -> f64

sci_form::compute_uff_energy(smiles: &str, coords: &[f64]) -> Result<f64, String>

sci_form::compute_mmff94_energy(smiles: &str, coords: &[f64]) -> Result<f64, String>

// PM3 — NDDO semi-empirical SCF
sci_form::compute_pm3(elements: &[u8], positions: &[[f64;3]]) -> Result<pm3::Pm3Result, String>
// Pm3Result { orbital_energies, electronic_energy, nuclear_repulsion, total_energy,
//             heat_of_formation, n_basis, n_electrons, homo_energy, lumo_energy, gap,
//             mulliken_charges, scf_iterations, converged }

// GFN-xTB — GFN0 tight-binding SCC (25 elements)
sci_form::compute_xtb(elements: &[u8], positions: &[[f64;3]]) -> Result<xtb::XtbResult, String>
// XtbResult { orbital_energies, electronic_energy, repulsive_energy, total_energy,
//             n_basis, n_electrons, homo_energy, lumo_energy, gap,
//             mulliken_charges, scc_iterations, converged }
// Supported: H, B, C, N, O, F, Si, P, S, Cl, Br, I + Ti, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ru, Pd, Ag, Pt, Au

// ANI — Accurate Neural Network Potential (ANI-2x)
sci_form::compute_ani(elements: &[u8], positions: &[[f64;3]], config: &AniConfig, models: &HashMap<u8, FeedForwardNet>) -> Result<AniResult, String>
// AniResult { energy, forces: Vec<[f64;3]>, species, atomic_energies, aevs }

// HF-3c — Minimal basis Hartree-Fock with corrections (D3, gCP, SRB)
sci_form::solve_hf3c(elements: &[u8], positions: &[[f64;3]], config: &HfConfig) -> Result<Hf3cResult, String>
// Hf3cResult { total_energy, electronic_energy, nuclear_repulsion, scf_iterations, converged... }

// IR Spectroscopy
sci_form::compute_ir_spectrum(elements: &[u8], positions: &[[f64;3]], method: ScientificMethod) -> Result<IrSpectrum, String>
// IrSpectrum { peaks: Vec<IrPeak>, frequencies: Vec<f64>, intensities: Vec<f64> }
sci_form::compute_vibrational_analysis(elements: &[u8], positions: &[[f64;3]], method: ScientificMethod) -> Result<VibrationalAnalysis, String>
// VibrationalAnalysis { modes: Vec<VibrationalMode>, frequencies, intensities }

// ML Descriptors + Property Prediction
sci_form::compute_ml_descriptors(
    elements: &[u8],
    bonds: &[(usize, usize, String)],
    charges: &[f64],
    aromatic_atoms: &[bool],
) -> ml::MolecularDescriptors
// MolecularDescriptors { molecular_weight, n_heavy_atoms, n_hydrogens, n_bonds,
//   n_rotatable_bonds, n_hbd, n_hba, fsp3, total_abs_charge, max_charge, min_charge,
//   wiener_index, n_rings, n_aromatic, balaban_j, sum_electronegativity, sum_polarizability }

sci_form::predict_ml_properties(desc: &ml::MolecularDescriptors) -> ml::MlPropertyResult
// MlPropertyResult { logp, molar_refractivity, log_solubility, lipinski, druglikeness }
// LipinskiResult { mw_ok, logp_ok, hbd_ok, hba_ok, violations, passes }

// Materials
sci_form::create_unit_cell(a,b,c,alpha,beta,gamma: f64) -> materials::UnitCell
sci_form::assemble_framework(node: &Sbu, linker: &Sbu, topology: &Topology, cell: &UnitCell) -> CrystalStructure

// Stereochemistry — CIP priorities, R/S, E/Z
sci_form::analyze_stereo(smiles: &str, coords: &[f64]) -> Result<stereo::StereoAnalysis, String>
// StereoAnalysis { stereocenters: Vec<Stereocenter>, double_bonds: Vec<DoubleBondStereo>,
//                  n_stereocenters, n_double_bonds }
// Stereocenter { atom_index, element, substituent_indices, priorities, configuration: Option<"R"|"S"> }
// DoubleBondStereo { atom1, atom2, configuration: Option<"E"|"Z">, high_priority_sub1, high_priority_sub2 }

// Solvation — Non-polar SASA + Generalized Born
sci_form::compute_nonpolar_solvation(elements: &[u8], positions: &[[f64;3]], probe_radius: Option<f64>) -> NonPolarSolvation
// NonPolarSolvation { energy_kcal_mol, atom_contributions, atom_sasa, total_sasa }

sci_form::compute_gb_solvation(elements: &[u8], positions: &[[f64;3]], charges: &[f64],
    solvent_dielectric: Option<f64>, solute_dielectric: Option<f64>, probe_radius: Option<f64>) -> GbSolvation
// GbSolvation { electrostatic_energy_kcal_mol, nonpolar_energy_kcal_mol, total_energy_kcal_mol,
//               born_radii, charges, solvent_dielectric, solute_dielectric }

// Ring perception — SSSR (Horton's algorithm)
sci_form::compute_sssr(smiles: &str) -> Result<rings::sssr::SssrResult, String>
// SssrResult { rings: Vec<RingInfo>, atom_ring_count, atom_ring_sizes, ring_size_histogram }
// RingInfo { atoms: Vec<usize>, size, is_aromatic }

// Fingerprints — ECFP (Morgan) + Tanimoto similarity
sci_form::compute_ecfp(smiles: &str, radius: usize, n_bits: usize) -> Result<rings::ecfp::ECFPFingerprint, String>
// ECFPFingerprint { n_bits, on_bits: BTreeSet<usize>, radius, raw_features }

sci_form::compute_tanimoto(fp1: &ECFPFingerprint, fp2: &ECFPFingerprint) -> f64

// Clustering — Butina (Taylor-Butina) RMSD clustering
sci_form::butina_cluster(conformers: &[Vec<f64>], rmsd_cutoff: f64) -> clustering::ClusterResult
// ClusterResult { n_clusters, assignments, centroid_indices, cluster_sizes, rmsd_cutoff }

sci_form::compute_rmsd_matrix(conformers: &[Vec<f64>]) -> Vec<Vec<f64>>

// NMR — Ensemble-averaged J-couplings
sci_form::compute_ensemble_j_couplings(smiles: &str, conformer_coords: &[Vec<f64>],
    energies_kcal: &[f64], temperature_k: f64) -> Result<Vec<nmr::JCoupling>, String>

// IR — Broadened spectrum
sci_form::compute_ir_spectrum_broadened(analysis: &VibrationalAnalysis, gamma: f64,
    wn_min: f64, wn_max: f64, n_points: usize, broadening: &str) -> IrSpectrum

// Charges — Configurable Gasteiger-Marsili
sci_form::compute_charges_configured(smiles: &str, config: &GasteigerConfig) -> Result<ChargeResult, String>
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
let pop = sci_form::compute_population(&conf.elements, &pos).unwrap();
println!("Gap: {:.3} eV", pop.homo_lumo_gap);

let dipole = sci_form::compute_dipole(&conf.elements, &pos).unwrap();
println!("|μ| = {:.3} D", dipole.magnitude);

// PM3 semi-empirical SCF
let pm3 = sci_form::compute_pm3(&conf.elements, &pos).unwrap();
println!("PM3 HOF: {:.2} kcal/mol, gap: {:.3} eV", pm3.heat_of_formation, pm3.gap);

// GFN-xTB (also handles transition metals)
let xtb = sci_form::compute_xtb(&conf.elements, &pos).unwrap();
println!("xTB total energy: {:.4} eV, gap: {:.3} eV", xtb.total_energy, xtb.gap);

// ML properties (no 3D needed — only SMILES topology)
let desc = sci_form::compute_ml_descriptors(&conf.elements, &conf.bonds, &[], &[]);
let props = sci_form::predict_ml_properties(&desc);
println!("LogP: {:.2}, Druglikeness: {:.3}", props.logp, props.druglikeness);

// Parallel batch (requires features = ["parallel"])
let config = sci_form::ConformerConfig { seed: 42, num_threads: 0 };
let results = sci_form::embed_batch(&["CCO", "c1ccccc1", "CC(=O)O"], &config);

// Stereochemistry
let stereo = sci_form::analyze_stereo("C(F)(Cl)(Br)I", &conf.coords).unwrap();
println!("Stereocenters: {}", stereo.n_stereocenters);

// Ring perception + fingerprints
let rings = sci_form::compute_sssr("c1ccccc1").unwrap();
println!("Rings: {:?}", rings.ring_size_histogram);
let fp1 = sci_form::compute_ecfp("c1ccccc1", 2, 2048).unwrap();
let fp2 = sci_form::compute_ecfp("Cc1ccccc1", 2, 2048).unwrap();
println!("Tanimoto: {:.3}", sci_form::compute_tanimoto(&fp1, &fp2));

// Solvation
let solv = sci_form::compute_nonpolar_solvation(&conf.elements, &pos, None);
println!("Non-polar solvation: {:.2} kcal/mol", solv.energy_kcal_mol);
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
)

# Conformer
embed(smiles: str, seed: int = 42) -> ConformerResult
embed_batch(smiles_list: list[str], seed: int = 42, num_threads: int = 0) -> list[ConformerResult]
parse(smiles: str) -> dict   # {"num_atoms", "atoms", "bonds"}

# Properties (coords is flat list[float])
charges(smiles: str) -> ChargeResult          # .charges, .total_charge
sasa(elements, coords, probe_radius=1.4) -> SasaResult   # .total_sasa Å², .atom_sasa
population(elements, coords) -> PopulationResult         # .mulliken_charges, .lowdin_charges, .homo_energy, .lumo_energy
dipole(elements, coords) -> DipoleResult                 # .magnitude Debye, .vector [x,y,z]
dos(elements, coords, sigma=0.3, e_min=-30.0, e_max=5.0, n_points=500) -> DosResult  # .energies, .total_dos
rmsd(coords, reference) -> AlignmentResult              # .rmsd Å, .aligned_coords
uff_energy(smiles, coords) -> float                     # kcal/mol
mmff94_energy(smiles, coords) -> float                  # kcal/mol

# PM3 — NDDO semi-empirical SCF
pm3_calculate(elements: list[int], coords: list[float]) -> Pm3Result
# Pm3Result: .orbital_energies, .electronic_energy, .nuclear_repulsion, .total_energy
#            .heat_of_formation, .n_basis, .n_electrons, .homo_energy, .lumo_energy
#            .gap, .mulliken_charges, .scf_iterations, .converged

# GFN-xTB — GFN0 tight-binding (25 elements incl. transition metals)
extb_calculate(elements: list[int], coords: list[float]) -> XtbResult
# XtbResult: .orbital_energies, .electronic_energy, .repulsive_energy, .total_energy
#            .n_basis, .n_electrons, .homo_energy, .lumo_energy
#            .gap, .mulliken_charges, .scc_iterations, .converged

# ML Descriptors + Property Prediction (topology only, no 3D needed)
ml_descriptors(smiles: str) -> MolecularDescriptors
# MolecularDescriptors: .molecular_weight, .n_heavy_atoms, .n_hbd, .n_hba, .fsp3,
#   .wiener_index, .n_rings, .n_aromatic, .balaban_j, .sum_electronegativity,
#   .sum_polarizability, .n_hydrogens, .n_bonds, .n_rotatable_bonds

ml_predict(smiles: str) -> MlPropertyResult
# MlPropertyResult: .logp, .molar_refractivity, .log_solubility, .druglikeness
#   .lipinski_violations (int), .lipinski_passes (bool)

# EHT
eht_calculate(elements, coords, k=1.75) -> EhtResult     # .homo_energy, .gap, .energies, .homo_index
eht_orbital_mesh(elements, coords, mo_index, spacing=0.2, isovalue=0.02) -> dict
# dict keys: "vertices", "normals", "indices", "num_triangles", "isovalue"

# Materials
unit_cell(a,b,c, alpha=90,beta=90,gamma=90) -> UnitCell  # .volume, .lattice
assemble_framework(topology="pcu", metal=30, geometry="octahedral", lattice_a=10.0, supercell=1) -> CrystalStructure
# CrystalStructure: .num_atoms, .elements, .cart_coords, .frac_coords

# Transport
pack_conformers(results: list[ConformerResult]) -> RecordBatch
split_worker_tasks(smiles, n_workers=4, seed=42) -> list[dict]
estimate_workers(n_items, max_workers=8) -> int

# Stereochemistry — CIP priorities, R/S, E/Z
stereo_analysis(smiles: str, coords: list[float] = []) -> StereoAnalysisPy
# StereoAnalysisPy: .stereocenters, .double_bonds, .n_stereocenters, .n_double_bonds
# StereocenterPy: .atom_index, .element, .substituent_indices, .priorities, .configuration
# DoubleBondStereoPy: .atom1, .atom2, .configuration, .high_priority_sub1, .high_priority_sub2

# Solvation — Non-polar SASA + Generalized Born
nonpolar_solvation(elements, coords, probe_radius=1.4) -> NonPolarSolvationPy
# NonPolarSolvationPy: .energy_kcal_mol, .atom_contributions, .atom_sasa, .total_sasa

gb_solvation(elements, coords, charges, solvent_dielectric=78.5, solute_dielectric=1.0, probe_radius=1.4) -> GbSolvationPy
# GbSolvationPy: .electrostatic_energy_kcal_mol, .nonpolar_energy_kcal_mol, .total_energy_kcal_mol,
#                .born_radii, .charges, .solvent_dielectric, .solute_dielectric

# Ring perception — SSSR
sssr(smiles: str) -> SssrResultPy
# SssrResultPy: .rings, .atom_ring_count, .atom_ring_sizes, .ring_size_histogram
# RingInfoPy: .atoms, .size, .is_aromatic

# ECFP fingerprints + Tanimoto
ecfp(smiles: str, radius=2, n_bits=2048) -> ECFPFingerprintPy
# ECFPFingerprintPy: .n_bits, .on_bits, .radius, .density

tanimoto(fp1: ECFPFingerprintPy, fp2: ECFPFingerprintPy) -> float

# Clustering — Butina RMSD
butina_cluster(conformers: list[list[float]], rmsd_cutoff=1.0) -> ClusterResultPy
# ClusterResultPy: .n_clusters, .assignments, .centroid_indices, .cluster_sizes, .rmsd_cutoff

rmsd_matrix(conformers: list[list[float]]) -> list[list[float]]

filter_diverse(conformers: list[list[float]], rmsd_cutoff=1.0) -> list[int]
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
print(f"HOMO: {pop.homo_energy:.3f} eV, Gap: {pop.homo_energy - pop.lumo_energy:.3f} eV")

d = dipole(conf.elements, conf.coords)
print(f"|μ| = {d.magnitude:.3f} D")

# PM3 semi-empirical
pm3 = pm3_calculate(conf.elements, conf.coords)
print(f"PM3 HOF: {pm3.heat_of_formation:.2f} kcal/mol")

# GFN-xTB (also works on [Zn], [Fe], etc.)
extb = xtb_calculate(conf.elements, conf.coords)
print(f"xTB gap: {xtb.gap:.3f} eV, converged: {xtb.converged}")

# ML properties (no 3D needed)
props = ml_predict("CCO")
print(f"LogP: {props.logp:.2f}, passes Lipinski: {props.lipinski_passes}")

# Stereochemistry
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

# UFF force field energy
sci-form energy "CCO" --from-smiles
sci-form energy "CCO" --coords coords.json

# Build info
sci-form info

# Stereochemistry analysis
sci-form stereo "C(F)(Cl)(Br)I"
sci-form stereo "C/C=C/C" --coords "[...]"

# Ring detection (SSSR)
sci-form sssr "c1ccccc1"

# ECFP fingerprint
sci-form ecfp "CCO" --radius 2 --n-bits 2048

# Tanimoto similarity
sci-form tanimoto "c1ccccc1" "Cc1ccccc1"

# Solvation energy
sci-form solvation "[8,1,1]" "[0,0,0,0.757,0.586,0,-0.757,0.586,0]"
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
| Dipole moment | Debye |
| SASA | Å² |
| Lattice parameters a,b,c | Å |
| Lattice angles α,β,γ | degrees |

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
