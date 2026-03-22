# Rust API Reference

## Crate: `sci-form`

```toml
[dependencies]
sci-form = "0.9"
# For parallel batch processing:
sci-form = { version = "0.9", features = ["parallel"] }
# GPU acceleration:
sci-form = { version = "0.9", features = ["experimental-gpu"] }
```

---

## Top-Level Types

### `ConformerResult`

```rust
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConformerResult {
    pub smiles: String,
    pub num_atoms: usize,
    pub coords: Vec<f64>,                    // [x₀, y₀, z₀, x₁, y₁, z₁, ...]
    pub elements: Vec<u8>,                   // atomic numbers
    pub bonds: Vec<(usize, usize, String)>,  // (atom_a, atom_b, order)
    pub error: Option<String>,
    pub time_ms: f64,
}
```

| Field | Description |
|-------|-------------|
| `smiles` | Input SMILES string |
| `num_atoms` | Total atoms including implicit H |
| `coords` | Flat 3D coordinates, length = `num_atoms × 3` |
| `elements` | Atomic numbers: 1=H, 6=C, 7=N, 8=O, 9=F, … |
| `bonds` | `(start, end, order)` where order ∈ {"SINGLE","DOUBLE","TRIPLE","AROMATIC"} |
| `error` | `None` on success, `Some(msg)` on failure |
| `time_ms` | Generation time in milliseconds |

### `ConformerConfig`

```rust
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConformerConfig {
    pub seed: u64,         // default: 42
    pub num_threads: usize, // default: 0 = auto-detect
}
```

---

## Conformer Generation

### `embed`

```rust
pub fn embed(smiles: &str, seed: u64) -> ConformerResult
```

Generate a 3D conformer from a SMILES string using ETKDGv2. Always returns a `ConformerResult` — check `error` field for failures.

**Example:**

```rust
let result = sci_form::embed("CCO", 42);
assert!(result.error.is_none());
println!("{} atoms at ({:.2}, {:.2}, {:.2})",
    result.num_atoms, result.coords[0], result.coords[1], result.coords[2]);
```

### `embed_batch`

```rust
pub fn embed_batch(smiles_list: &[&str], config: &ConformerConfig) -> Vec<ConformerResult>
```

Batch-generate conformers. With `parallel` feature uses rayon; `config.num_threads = 0` auto-detects CPU count.

**Example:**

```rust
use sci_form::{embed_batch, ConformerConfig};

let smiles = vec!["CCO", "c1ccccc1", "CC(=O)O"];
let cfg = ConformerConfig { seed: 42, num_threads: 0 };
let results = embed_batch(&smiles, &cfg);
for r in &results {
    println!("{}: {} atoms", r.smiles, r.num_atoms);
}
```

### `parse`

```rust
pub fn parse(smiles: &str) -> Result<graph::Molecule, String>
```

Parse SMILES into a molecular graph without generating 3D coordinates. Returns `Ok(Molecule)` or `Err(msg)`.

---

## Gasteiger Charges

### `compute_charges`

```rust
pub fn compute_charges(smiles: &str) -> Result<charges::gasteiger::ChargeResult, String>
```

Compute Gasteiger-Marsili partial charges using default settings (6 iterations).

**Returns:**

```rust
pub struct ChargeResult {
    pub charges: Vec<f64>,   // per-atom partial charges
    pub total_charge: f64,   // sum (should be ~0 for neutral)
}
```

**Example:**

```rust
let result = sci_form::compute_charges("CCO").unwrap();
println!("O charge: {:.4}", result.charges[2]); // ~-0.38
```

### `compute_charges_configured`

```rust
pub fn compute_charges_configured(
    smiles: &str,
    config: &charges::gasteiger::GasteigerConfig,
) -> Result<charges::gasteiger::ChargeResult, String>
```

Compute Gasteiger-Marsili charges with fine-grained control.

```rust
pub struct GasteigerConfig {
    pub max_iter: usize,               // default: 6
    pub initial_damping: f64,          // default: 1.0
    pub convergence_threshold: f64,    // default: 1e-8
}
```

Supports all main-group elements up to period 4 (Li, Be, Na, Mg, Al, Si, …) and formal-charge-aware initialization for charged species.

**Example:**

```rust
use sci_form::charges::gasteiger::GasteigerConfig;

let config = GasteigerConfig {
    max_iter: 12,
    initial_damping: 0.9,
    convergence_threshold: 1e-10,
};
let result = sci_form::compute_charges_configured("[NH4+]", &config).unwrap();
println!("Total charge: {:.2}", result.total_charge); // ~+1.0
```

---

## Solvent-Accessible Surface Area

### `compute_sasa`

```rust
pub fn compute_sasa(
    elements: &[u8],
    coords_flat: &[f64],
    probe_radius: Option<f64>,
) -> Result<surface::sasa::SasaResult, String>
```

Shrake-Rupley SASA computation using a Fibonacci sphere (default 100 points), Bondi radii, and optional probe radius (default 1.4 Å for water).

**Returns:**

```rust
pub struct SasaResult {
    pub total_sasa: f64,      // Å²
    pub per_atom: Vec<f64>,   // per-atom SASA in Å²
}
```

**Example:**

```rust
let conf = sci_form::embed("CCO", 42);
let sasa = sci_form::compute_sasa(&conf.elements, &conf.coords, Some(1.4)).unwrap();
println!("SASA: {:.2} Å²", sasa.total_sasa);
```

---

## EHT & Population Analysis

### `compute_population`

```rust
pub fn compute_population(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<population::PopulationResult, String>
```

Runs EHT then computes Mulliken and Löwdin population analysis.

**Returns:**

```rust
pub struct PopulationResult {
    pub mulliken_charges: Vec<f64>,  // per-atom Mulliken charge
    pub lowdin_charges: Vec<f64>,    // per-atom Löwdin charge
    pub homo_energy: f64,            // eV
    pub lumo_energy: f64,            // eV
    pub homo_lumo_gap: f64,          // eV
}
```

**Example:**

```rust
let conf = sci_form::embed("CCO", 42);
let pos: Vec<[f64; 3]> = conf.coords.chunks(3)
    .map(|c| [c[0], c[1], c[2]])
    .collect();
let pop = sci_form::compute_population(&conf.elements, &pos).unwrap();
println!("HOMO: {:.3} eV, Gap: {:.3} eV", pop.homo_energy, pop.homo_lumo_gap);
```

### `compute_dipole`

```rust
pub fn compute_dipole(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<dipole::DipoleResult, String>
```

Runs EHT then computes the molecular dipole moment.

**Returns:**

```rust
pub struct DipoleResult {
    pub vector: [f64; 3],   // Debye, [x, y, z]
    pub magnitude: f64,     // Debye
}
```

---

## Electrostatic Potential

### `compute_esp`

```rust
pub fn compute_esp(
    elements: &[u8],
    positions: &[[f64; 3]],
    spacing: f64,
    padding: f64,
) -> Result<esp::EspGrid, String>
```

Compute an ESP grid using Mulliken charges from EHT. Grid extends `padding` Å beyond the molecular bounding box.

**Returns:**

```rust
pub struct EspGrid {
    pub values: Vec<f64>,       // grid values, row-major (z-innermost)
    pub origin: [f64; 3],       // Å
    pub spacing: f64,           // Å
    pub nx: usize, pub ny: usize, pub nz: usize,
}
```

**Advanced (parallel, color mapping):**

```rust
use sci_form::esp::{compute_esp_grid_parallel, esp_grid_to_colors, esp_color_map};

// Parallel grid (requires "parallel" feature)
let grid = compute_esp_grid_parallel(&elements, &pos, &mulliken_charges, 0.5, 3.0);

// Color-map the whole grid to RGB
let colors: Vec<u8> = esp_grid_to_colors(&grid, 0.05);  // range = ±0.05 a.u.

// Color-map a single value
let rgb: [u8; 3] = esp_color_map(0.03, 0.05);  // returns [r, g, b]

// Export to .cube file
use sci_form::esp::export_cube;
export_cube(&grid, &elements, &pos, "molecule.cube").unwrap();
```

---

## Density of States

### `compute_dos`

```rust
pub fn compute_dos(
    elements: &[u8],
    positions: &[[f64; 3]],
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> Result<dos::DosResult, String>
```

Runs EHT to get orbital energies, then computes DOS and per-atom PDOS with Gaussian smearing $\sigma$.

**Returns:**

```rust
pub struct DosResult {
    pub energies: Vec<f64>,          // eV, n_points values from e_min to e_max
    pub dos: Vec<f64>,               // total DOS
    pub pdos: Vec<Vec<f64>>,         // pdos[atom][energy_idx]
    pub orbital_energies: Vec<f64>,  // raw EHT orbital energies (eV)
    pub homo_lumo_gap: Option<f64>,  // eV
}
```

**Advanced:**

```rust
use sci_form::dos::{compute_pdos, dos_mse, export_dos_json};

let dos = compute_pdos(&elements, &coords_flat, &orbital_energies, &coefficients,
                       n_electrons, 0.2, -20.0, 5.0, 200);

// Mean squared error between two DOS curves
let mse = dos_mse(&dos.dos, &reference_dos);

// Export to JSON string
let json = export_dos_json(&dos);
```

---

## Molecular Alignment

### `compute_rmsd`

```rust
pub fn compute_rmsd(coords: &[f64], reference: &[f64]) -> f64
```

RMSD after optimal Kabsch alignment (SVD-based). Both arrays are flat `[x0,y0,z0,…]`.

**Advanced:**

```rust
use sci_form::alignment::{align_coordinates, align_quaternion, compute_rmsd};

// Kabsch SVD alignment
let result = align_coordinates(&mobile, &reference).unwrap();
// result.rotation: [[f64;3];3], result.translation: [f64;3], result.rmsd: f64

// Quaternion alignment (Coutsias 2004, faster for large N)
let result = align_quaternion(&mobile, &reference).unwrap();
```

---

## Force Field Energy

### `compute_uff_energy`

```rust
pub fn compute_uff_energy(smiles: &str, coords: &[f64]) -> Result<f64, String>
```

Evaluate the UFF force field energy (kcal/mol) for a molecule at the given coordinates.

**Advanced — MMFF94:**

```rust
use sci_form::forcefield::mmff94::Mmff94Builder;

let mol = sci_form::parse("CCO").unwrap();
let mmff = Mmff94Builder::build(&mol);
let energy = mmff.total_energy(&coords);
```

---

## Materials

### `create_unit_cell`

```rust
pub fn create_unit_cell(
    a: f64, b: f64, c: f64,
    alpha: f64, beta: f64, gamma: f64,
) -> materials::UnitCell
```

Create a periodic unit cell from lattice parameters (a, b, c in Å; α, β, γ in degrees).

### `assemble_framework`

```rust
pub fn assemble_framework(
    node: &materials::Sbu,
    linker: &materials::Sbu,
    topology: &materials::Topology,
    cell: &materials::UnitCell,
) -> materials::CrystalStructure
```

Assemble a MOF-type framework crystal structure from node/linker SBUs on a topology.

**Example:**

```rust
use sci_form::materials::{Sbu, Topology, CellParameters};

let cell = sci_form::create_unit_cell(26.3, 26.3, 26.3, 90.0, 90.0, 90.0);
let node = Sbu { /* ... */ };
let linker = Sbu { /* ... */ };
let topology = Topology::Pcu;  // primitive cubic
let crystal = sci_form::assemble_framework(&node, &linker, &topology, &cell);
```

---

## Spectroscopy (Track D)

### `compute_stda_uvvis`

```rust
pub fn compute_stda_uvvis(
    elements: &[u8],
    positions: &[[f64; 3]],
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
    broadening: reactivity::BroadeningType,
) -> Result<reactivity::StdaUvVisSpectrum, String>
```

Compute a UV-Vis absorption spectrum via the simplified Tamm-Dancoff approximation (sTDA) on EHT MO transitions.

| Parameter | Description |
|-----------|-------------|
| `elements` | Atomic numbers |
| `positions` | 3D Cartesian coordinates (Å) |
| `sigma` | Broadening width in eV (Gaussian σ or Lorentzian γ) |
| `e_min` / `e_max` | Energy window in eV |
| `n_points` | Wavelength grid size |
| `broadening` | `BroadeningType::Gaussian` or `BroadeningType::Lorentzian` |

**Returns** `StdaUvVisSpectrum { wavelengths_nm, intensities, excitations, homo_energy, lumo_energy, gap, n_transitions, broadening_type, notes }`.

```rust
use sci_form::{embed, compute_stda_uvvis};
use sci_form::reactivity::BroadeningType;

let conf = embed("c1ccccc1", 42);
let pos: Vec<[f64;3]> = conf.coords.chunks(3).map(|c| [c[0],c[1],c[2]]).collect();
let spec = compute_stda_uvvis(&conf.elements, &pos, 0.3, 1.0, 8.0, 500, BroadeningType::Gaussian).unwrap();
println!("Gap: {:.2} eV, {} transitions", spec.gap, spec.n_transitions);
```

---

### `compute_vibrational_analysis`

```rust
pub fn compute_vibrational_analysis(
    elements: &[u8],
    positions: &[[f64; 3]],
    method: ir::ScientificMethod,
) -> Result<ir::VibrationalAnalysis, String>
```

Compute a numerical Hessian ($6N$ energy evaluations with mass-aware step size) and diagonalize to obtain vibrational frequencies, normal modes, IR intensities, ZPVE, and thermochemistry. Translation/rotation modes are projected out automatically.

**`ScientificMethod` enum:** `ScientificMethod::Eht`, `ScientificMethod::Pm3`, `ScientificMethod::Xtb`

**Returns** `VibrationalAnalysis { n_atoms, modes, n_real_modes, zpve_ev, method, notes, thermochemistry }`.

Each `VibrationalMode` has `frequency_cm1`, `ir_intensity` (km/mol), `displacement` (3N), `is_real`, `label` (functional group annotation).

`Thermochemistry` (RRHO, 298.15 K): `zpve_kcal`, `thermal_energy_kcal`, `entropy_vib_cal`, `gibbs_correction_kcal`.

```rust
use sci_form::ir::ScientificMethod;

let vib = sci_form::compute_vibrational_analysis(&conf.elements, &pos, ScientificMethod::Xtb).unwrap();
println!("ZPVE: {:.4} eV", vib.zpve_ev);
println!("Gibbs correction: {:.3} kcal/mol", vib.thermochemistry.gibbs_correction_kcal);
let strongest = vib.modes.iter().filter(|m| m.is_real)
    .max_by(|a,b| a.ir_intensity.partial_cmp(&b.ir_intensity).unwrap());
println!("Strongest: {:.1} cm⁻¹  ({:.1} km/mol)  [{}]",
    strongest.map(|m| m.frequency_cm1).unwrap_or(0.0),
    strongest.map(|m| m.ir_intensity).unwrap_or(0.0),
    strongest.and_then(|m| m.label.as_deref()).unwrap_or("?"));
```

---

### `compute_ir_spectrum`

```rust
pub fn compute_ir_spectrum(
    elements: &[u8],
    positions: &[[f64; 3]],
    method: ir::ScientificMethod,
) -> Result<ir::IrSpectrum, String>
```

Convenience wrapper: runs `compute_vibrational_analysis` and applies default Lorentzian broadening (γ = 15 cm⁻¹, 600–4000 cm⁻¹, 1000 points).

**Returns** `IrSpectrum { wavenumbers, intensities, transmittance, peaks, gamma, notes }`.

```rust
use sci_form::ir::ScientificMethod;

let spec = sci_form::compute_ir_spectrum(&conf.elements, &pos, ScientificMethod::Xtb).unwrap();
let max_i = spec.intensities.iter().cloned().fold(0.0_f64, f64::max);
let idx   = spec.intensities.iter().position(|&v| v == max_i).unwrap();
println!("Dominant band: {:.1} cm⁻¹", spec.wavenumbers[idx]);
```

---

### `compute_ir_spectrum_broadened`

```rust
pub fn compute_ir_spectrum_broadened(
    analysis: &ir::VibrationalAnalysis,
    gamma: f64,
    wn_min: f64,
    wn_max: f64,
    n_points: usize,
    broadening: &str,
) -> ir::IrSpectrum
```

Generate a broadened IR spectrum from an existing `VibrationalAnalysis` with full control over parameters.

| Parameter | Description |
|-----------|-------------|
| `gamma` | Half-width (Lorentzian) or σ (Gaussian) in cm⁻¹ |
| `broadening` | `"lorentzian"` (default) or `"gaussian"` |

**Returns** `IrSpectrum` with both `intensities` (absorbance) and `transmittance` axes.

```rust
// Gaussian broadening for publication-quality spectra
let spec = sci_form::compute_ir_spectrum_broadened(&vib, 10.0, 400.0, 4000.0, 2000, "gaussian");
println!("{} peaks labelled", spec.peaks.len());
```

---

### `predict_nmr_shifts`

```rust
pub fn predict_nmr_shifts(smiles: &str) -> Result<nmr::NmrShiftResult, String>
```

Predict ¹H and ¹³C chemical shifts from SMILES using HOSE code environment matching.

**Returns** `NmrShiftResult { h_shifts, c_shifts, notes }`. Each `ChemicalShift` has `atom_index`, `element`, `shift_ppm`, `environment`, `confidence`.

```rust
let shifts = sci_form::predict_nmr_shifts("CCO").unwrap();
for h in &shifts.h_shifts {
    println!("H#{}: {:.2} ppm ({})", h.atom_index, h.shift_ppm, h.environment);
}
```

---

### `predict_nmr_couplings`

```rust
pub fn predict_nmr_couplings(
    smiles: &str,
    positions: &[[f64; 3]],
) -> Result<Vec<nmr::JCoupling>, String>
```

Predict $^2$J, $^3$J, and $^4$J H–H coupling constants. Pass 3D positions for Karplus dihedral-angle evaluation; pass `&[]` for topology-based free-rotation averages.

Each `JCoupling` has `h1_index`, `h2_index`, `j_hz`, `n_bonds`, `coupling_type`.

Parameterized pathways: H-C-C-H (Altona-Sundaralingam), H-C-N-H (Bystrov), H-C-O-H, H-C-S-H.

---

### `compute_ensemble_j_couplings`

```rust
pub fn compute_ensemble_j_couplings(
    smiles: &str,
    conformer_coords: &[Vec<f64>],
    energies_kcal: &[f64],
    temperature_k: f64,
) -> Result<Vec<nmr::JCoupling>, String>
```

Boltzmann-average ³J couplings over a conformer ensemble.

| Parameter | Description |
|-----------|-------------|
| `conformer_coords` | Vec of flat coordinate arrays (one per conformer) |
| `energies_kcal` | Relative energies in kcal/mol (UFF or xTB) |
| `temperature_k` | Temperature for Boltzmann weighting (default: 298.15 K) |

```rust
let confs: Vec<Vec<f64>> = conf_results.iter().map(|r| r.coords.clone()).collect();
let energies: Vec<f64> = conf_results.iter().map(|r| r.uff_energy).collect();
let couplings = sci_form::compute_ensemble_j_couplings("CCCC", &confs, &energies, 298.15).unwrap();
for jc in &couplings {
    println!("H{}–H{}: {:.2} Hz", jc.h1_index, jc.h2_index, jc.j_hz);
}
```

---

### `compute_nmr_spectrum`

```rust
pub fn compute_nmr_spectrum(
    smiles: &str,
    nucleus: &str,     // "1H" or "13C"
    gamma: f64,
    ppm_min: f64,
    ppm_max: f64,
    n_points: usize,
) -> Result<nmr::NmrSpectrum, String>
```

Full NMR spectrum pipeline: SMILES → shifts → couplings → multiplet splitting → Lorentzian broadening.

**Returns** `NmrSpectrum { ppm_axis, intensities, peaks, nucleus, gamma }`.

```rust
let spec = sci_form::compute_nmr_spectrum("CCO", "1H", 0.01, -2.0, 12.0, 2000).unwrap();
println!("{} ¹H multiplet lines in spectrum", spec.peaks.len());
```

---

### `compute_hose_codes`

```rust
pub fn compute_hose_codes(smiles: &str, max_radius: usize) -> Result<Vec<nmr::HoseCode>, String>
```

Generate HOSE code environment strings for all atoms in the molecule (breadth-first spherical atom environments).

```rust
let codes = sci_form::compute_hose_codes("c1ccccc1", 4).unwrap();
for hose in &codes {
    println!("Atom {}: {}", hose.atom_index, hose.code_string);
}
```

---

## Stereochemistry

### `analyze_stereo`

```rust
pub fn analyze_stereo(
    smiles: &str,
    coords: &[f64],
) -> Result<stereo::StereoAnalysis, String>
```

Assign CIP priorities and determine R/S configuration at stereocenters and E/Z geometry at double bonds using 3D coordinates.

```rust
pub struct StereoAnalysis {
    pub stereocenters: Vec<Stereocenter>,
    pub double_bonds:  Vec<DoubleBondStereo>,
    pub n_stereocenters: usize,
    pub n_double_bonds:  usize,
}

pub struct Stereocenter {
    pub atom_index: usize,
    pub element: u8,
    pub substituent_indices: Vec<usize>,
    pub priorities: Vec<usize>,           // CIP priority rank (1 = highest)
    pub configuration: Option<String>,    // "R" or "S"
}

pub struct DoubleBondStereo {
    pub atom1: usize,
    pub atom2: usize,
    pub configuration: Option<String>,    // "E" or "Z"
    pub high_priority_sub1: Option<usize>,
    pub high_priority_sub2: Option<usize>,
}
```

CIP logic: DFS priority tree, duplicate-atom expansion for multiple bonds, recursive tie-breaking (up to depth 10), triple-product sign test for R/S, dihedral angle for E/Z.

**Example:**

```rust
let conf = sci_form::embed("C(F)(Cl)(Br)I", 42);
let stereo = sci_form::analyze_stereo("C(F)(Cl)(Br)I", &conf.coords).unwrap();
println!("{} stereocenters", stereo.n_stereocenters);
for sc in &stereo.stereocenters {
    println!("Atom {}: {:?}", sc.atom_index, sc.configuration);
}
```

---

## Implicit Solvation

### `compute_nonpolar_solvation`

```rust
pub fn compute_nonpolar_solvation(
    elements: &[u8],
    positions: &[[f64; 3]],
    probe_radius: Option<f64>,
) -> solvation::NonPolarSolvation
```

Non-polar solvation free energy: $\Delta G_{\text{np}} = \sum_i \sigma_i \cdot A_i$ using per-element atomic solvation parameters (ASP) and Shrake-Rupley SASA.

```rust
pub struct NonPolarSolvation {
    pub energy_kcal_mol: f64,
    pub atom_contributions: Vec<f64>,  // kcal/mol per atom
    pub atom_sasa: Vec<f64>,           // Å² per atom
    pub total_sasa: f64,               // Å²
}
```

### `compute_gb_solvation`

```rust
pub fn compute_gb_solvation(
    elements: &[u8],
    positions: &[[f64; 3]],
    charges: &[f64],
    solvent_dielectric: Option<f64>,   // default: 78.5 (water)
    solute_dielectric: Option<f64>,    // default: 1.0
    probe_radius: Option<f64>,         // default: 1.4 Å
) -> solvation::GbSolvation
```

Generalized Born + Surface Area (GB/SA) electrostatic solvation using the HCT pairwise descreening model for effective Born radii and the Still GB equation.

```rust
pub struct GbSolvation {
    pub electrostatic_energy_kcal_mol: f64,
    pub nonpolar_energy_kcal_mol: f64,
    pub total_energy_kcal_mol: f64,
    pub born_radii: Vec<f64>,
    pub charges: Vec<f64>,
    pub solvent_dielectric: f64,
    pub solute_dielectric: f64,
}
```

**Example:**

```rust
let conf = sci_form::embed("CCO", 42);
let pos: Vec<[f64;3]> = conf.coords.chunks(3).map(|c| [c[0],c[1],c[2]]).collect();
let charges = sci_form::compute_charges("CCO").unwrap().charges;

let np = sci_form::compute_nonpolar_solvation(&conf.elements, &pos, None);
println!("Non-polar ΔG: {:.2} kcal/mol", np.energy_kcal_mol);

let gb = sci_form::compute_gb_solvation(&conf.elements, &pos, &charges, None, None, None);
println!("Total solvation: {:.2} kcal/mol", gb.total_energy_kcal_mol);
```

---

## Ring Perception (SSSR)

### `compute_sssr`

```rust
pub fn compute_sssr(smiles: &str) -> Result<rings::sssr::SssrResult, String>
```

Smallest Set of Smallest Rings using Horton's shortest-path algorithm with Hückel aromaticity perception.

```rust
pub struct SssrResult {
    pub rings: Vec<RingInfo>,
    pub atom_ring_count: Vec<usize>,              // rings containing each atom
    pub atom_ring_sizes: Vec<Vec<usize>>,         // ring sizes per atom
    pub ring_size_histogram: HashMap<usize, usize>, // size → count
}

pub struct RingInfo {
    pub atoms: Vec<usize>,   // atom indices (ring traversal order)
    pub size: usize,
    pub is_aromatic: bool,   // Hückel 4n+2 rule
}
```

**Example:**

```rust
let r = sci_form::compute_sssr("c1ccc2ccccc2c1").unwrap();  // naphthalene
println!("{} rings, sizes: {:?}", r.rings.len(),
    r.rings.iter().map(|r| r.size).collect::<Vec<_>>());
```

---

## Fingerprints (ECFP)

### `compute_ecfp`

```rust
pub fn compute_ecfp(
    smiles: &str,
    radius: usize,
    n_bits: usize,
) -> Result<rings::ecfp::ECFPFingerprint, String>
```

Extended-Connectivity Fingerprint (Morgan algorithm) with configurable radius and bit-vector size.

```rust
pub struct ECFPFingerprint {
    pub n_bits: usize,
    pub on_bits: BTreeSet<usize>,
    pub radius: usize,
    pub raw_features: Vec<u64>,
}
impl ECFPFingerprint {
    pub fn density(&self) -> f64  // fraction of set bits
}
```

Atom invariants: element, degree, valence, ring membership, aromaticity, formal charge.

### `compute_tanimoto`

```rust
pub fn compute_tanimoto(
    fp1: &rings::ecfp::ECFPFingerprint,
    fp2: &rings::ecfp::ECFPFingerprint,
) -> f64
```

Jaccard-Tanimoto similarity: $T = |A \cap B| / |A \cup B|$. Returns 1.0 for identical fingerprints, 0.0 for disjoint.

**Example:**

```rust
let fp_benz = sci_form::compute_ecfp("c1ccccc1", 2, 2048).unwrap();
let fp_tol  = sci_form::compute_ecfp("Cc1ccccc1", 2, 2048).unwrap();
let t = sci_form::compute_tanimoto(&fp_benz, &fp_tol);
println!("Tanimoto benzene/toluene: {:.3}", t);  // ~0.5–0.7
```

---

## Conformer Clustering

### `butina_cluster`

```rust
pub fn butina_cluster(
    conformers: &[Vec<f64>],
    rmsd_cutoff: f64,
) -> clustering::ClusterResult
```

Taylor-Butina single-linkage RMSD clustering. Builds an O(N²) RMSD matrix with Kabsch alignment, then assigns conformers to cluster seeds greedily (largest first).

```rust
pub struct ClusterResult {
    pub n_clusters: usize,
    pub assignments: Vec<usize>,       // cluster index per conformer
    pub centroid_indices: Vec<usize>,  // representative conformer per cluster
    pub cluster_sizes: Vec<usize>,
    pub rmsd_cutoff: f64,
}
```

### `compute_rmsd_matrix`

```rust
pub fn compute_rmsd_matrix(conformers: &[Vec<f64>]) -> Vec<Vec<f64>>
```

O(N²) pairwise RMSD matrix with Kabsch alignment.

### `filter_diverse_conformers`

```rust
pub fn filter_diverse_conformers(
    conformers: &[Vec<f64>],
    rmsd_cutoff: f64,
) -> Vec<usize>
```

Return indices of diverse representative conformers (one per cluster centroid).

**Example:**

```rust
use sci_form::{embed_batch, ConformerConfig, butina_cluster};

let smiles = vec!["CCCC"; 20];
let cfg = ConformerConfig { seed: 42, num_threads: 0 };
let results = embed_batch(&smiles, &cfg);
let coords: Vec<Vec<f64>> = results.iter().map(|r| r.coords.clone()).collect();

let clusters = sci_form::butina_cluster(&coords, 1.0);
println!("{} clusters from {} conformers", clusters.n_clusters, results.len());
```

---

## Internal Module Reference

| Module | Key Types / Functions |
|--------|-----------------------|
| `sci_form::graph` | `Molecule`, `Atom`, `Bond`, `BondOrder`, `Hybridization` |
| `sci_form::smiles` | SMILES parser (60+ elements, implicit H inferral) |
| `sci_form::smarts` | SMARTS pattern engine (846 CSD torsion patterns) |
| `sci_form::conformer` | `generate_3d_conformer()`, adaptive large-mol params |
| `sci_form::distgeom` | `BoundsMatrix`, `smooth_bounds()`, `compute_initial_coords_rdkit()` |
| `sci_form::distgeom::parallel` | `compute_bounds_parallel()`, `floyd_warshall_parallel()` (requires `parallel` feature) |
| `sci_form::etkdg` | Torsion FF, `TorsionPattern`, macrocycle detection |
| `sci_form::forcefield` | `build_uff_force_field()`, `Mmff94Builder`, UFF param tables |
| `sci_form::optimization` | L-BFGS minimizer |
| `sci_form::eht` | `solve_eht()`, `EhtResult`, `evaluate_orbital_on_grid_chunked()`, `marching_cubes()` |
| `sci_form::eht::band_structure` | `compute_band_structure()`, `BandStructure`, `BandStructureConfig` |
| `sci_form::eht::gradients` | `compute_eht_gradient()`, `optimize_geometry_eht()`, `EhtGradient`, `EhtOptResult` |
| `sci_form::esp` | `compute_esp_grid()`, `compute_esp_grid_parallel()`, `esp_color_map()`, `export_cube()`, `read_cube()` |
| `sci_form::dos` | `compute_dos()`, `compute_pdos()`, `dos_mse()`, `export_dos_json()` |
| `sci_form::dos::multi_method` | `compute_dos_multimethod()`, `DosMethod`, `MultiMethodDosResult` |
| `sci_form::alignment` | `align_coordinates()`, `align_quaternion()`, `compute_rmsd()` |
| `sci_form::charges::gasteiger` | `gasteiger_marsili_charges()`, `gasteiger_marsili_charges_configured()`, `ChargeResult`, `GasteigerConfig` |
| `sci_form::surface::sasa` | `compute_sasa()`, `SasaResult` |
| `sci_form::population` | `compute_population()`, `PopulationResult` |
| `sci_form::population::npa` | `compute_npa()`, `compute_nbo()`, `NpaResult`, `NboResult`, `NboBond`, `NboLonePair` |
| `sci_form::dipole` | `compute_dipole_from_eht()`, `DipoleResult` |
| `sci_form::pm3` | `compute_pm3()`, `Pm3Result` (PM3 + PM3(tm) transition metals) |
| `sci_form::xtb` | `compute_xtb()`, `XtbResult` (GFN0) |
| `sci_form::xtb::gfn1` | `solve_gfn1()`, `Gfn1Result` (shell charges, D3 dispersion) |
| `sci_form::xtb::gfn2` | `solve_gfn2()`, `Gfn2Result` (multipole, D4 dispersion, halogen bonding) |
| `sci_form::hf` | `solve_hf3c()`, `HfConfig`, `Hf3cResult` (D3 + gCP + SRB corrections) |
| `sci_form::hf::cisd` | `compute_cisd()`, `CisdResult`, `CisdExcitation` (excited states) |
| `sci_form::scf::extended_basis` | `build_basis_set()`, `BasisSetType` (STO-3G, 3-21G, 6-31G) |
| `sci_form::ani` | `compute_ani()`, `AniConfig`, `AniResult`, `FeedForwardNet` |
| `sci_form::ani::ani_tm` | `compute_aevs_tm()`, `species_index_tm()` (24 elements incl. transition metals) |
| `sci_form::materials` | `UnitCell`, `Sbu`, `Topology`, `CrystalStructure`, `assemble_framework()` |
| `sci_form::materials::space_groups` | `space_group_by_number()`, `space_group_by_symbol()`, `SpaceGroup` (230 ITC groups) |
| `sci_form::materials::geometry_opt` | `optimize_framework()`, `FrameworkOptConfig`, `FrameworkOptResult`, `OptMethod` |
| `sci_form::transport` | `pack_batch_arrow()`, `ChunkedIterator`, `WorkerTask`, `split_worker_tasks()` |
| `sci_form::reactivity` | `StdaUvVisSpectrum`, `StdaExcitation`, `BroadeningType`, `compute_stda_uvvis_spectrum()` |
| `sci_form::ir` | `VibrationalAnalysis`, `VibrationalMode`, `IrSpectrum`, `IrPeak`, `compute_vibrational_analysis()`, `compute_ir_spectrum()` |
| `sci_form::ir::peak_assignment` | `assign_peaks()`, `AssignmentResult`, `PeakAssignment` (functional group recognition) |
| `sci_form::nmr` | `NmrShiftResult`, `ChemicalShift`, `JCoupling`, `NmrSpectrum`, `NmrPeak`, `HoseCode`, `NmrNucleus` |
| `sci_form::stereo` | `analyze_stereo()`, `StereoAnalysis`, `Stereocenter`, `DoubleBondStereo`, `HelicalChirality` |
| `sci_form::solvation` | `compute_nonpolar_solvation()`, `compute_gb_solvation()`, `compute_born_radii()`, `NonPolarSolvation`, `GbSolvation` |
| `sci_form::rings::sssr` | `compute_sssr()`, `SssrResult`, `RingInfo` |
| `sci_form::rings::ecfp` | `compute_ecfp()`, `compute_tanimoto()`, `ECFPFingerprint` |
| `sci_form::clustering` | `butina_cluster()`, `compute_rmsd_matrix()`, `filter_diverse_conformers()`, `ClusterResult` |
| `sci_form::ml` | `compute_ml_descriptors()`, `predict_ml_properties()`, `MolecularDescriptors`, `MlPropertyResult` |
| `sci_form::ml::whim` | `compute_whim()`, `WhimDescriptors`, `WhimWeighting` (5 weighting schemes) |
| `sci_form::ml::rdf_descriptors` | `compute_rdf()`, `RdfDescriptors` (4 weighting types) |
| `sci_form::ml::getaway` | `compute_getaway()`, `GetawayDescriptors` |
| `sci_form::ml::advanced_models` | `train_random_forest()`, `train_gradient_boosting()`, `RandomForest`, `GradientBoosting` |
| `sci_form::periodic` | `build_periodic_molecule()`, `detect_hapticity()`, `PeriodicMolecule`, `HapticityAnalysis` |
| `sci_form::smirks` | `parse_smirks()`, `apply_smirks()`, `SmirksTransform`, `SmirksResult` |
| `sci_form::gpu::mmff94_gpu` | `compute_mmff94_nonbonded_gpu()` (requires `experimental-gpu` feature) |
