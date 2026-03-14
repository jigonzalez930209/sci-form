# Rust API Reference

## Crate: `sci-form`

```toml
[dependencies]
sci-form = "0.2"
# For parallel batch processing and ESP parallel:
sci-form = { version = "0.2", features = ["parallel"] }
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

Compute Gasteiger-Marsili partial charges (6 iterations of electronegativity equalization).

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

## Internal Module Reference

| Module | Key Types / Functions |
|--------|-----------------------|
| `sci_form::graph` | `Molecule`, `Atom`, `Bond`, `BondOrder`, `Hybridization` |
| `sci_form::smiles` | SMILES parser (60+ elements, implicit H inferral) |
| `sci_form::smarts` | SMARTS pattern engine (846 CSD torsion patterns) |
| `sci_form::conformer` | `generate_3d_conformer()`, adaptive large-mol params |
| `sci_form::distgeom` | `BoundsMatrix`, `smooth_bounds()`, `compute_initial_coords_rdkit()` |
| `sci_form::etkdg` | Torsion FF, `TorsionPattern`, macrocycle detection |
| `sci_form::forcefield` | `build_uff_force_field()`, `Mmff94Builder`, UFF param tables |
| `sci_form::optimization` | L-BFGS minimizer |
| `sci_form::eht` | `solve_eht()`, `EhtResult`, `evaluate_orbital_on_grid_chunked()`, `marching_cubes()` |
| `sci_form::esp` | `compute_esp_grid()`, `compute_esp_grid_parallel()`, `esp_color_map()`, `export_cube()`, `read_cube()` |
| `sci_form::dos` | `compute_dos()`, `compute_pdos()`, `dos_mse()`, `export_dos_json()` |
| `sci_form::alignment` | `align_coordinates()`, `align_quaternion()`, `compute_rmsd()` |
| `sci_form::charges::gasteiger` | `gasteiger_marsili_charges()`, `ChargeResult` |
| `sci_form::surface::sasa` | `compute_sasa()`, `SasaResult` |
| `sci_form::population` | `compute_population()`, `PopulationResult` |
| `sci_form::dipole` | `compute_dipole_from_eht()`, `DipoleResult` |
| `sci_form::materials` | `UnitCell`, `Sbu`, `Topology`, `CrystalStructure`, `assemble_framework()` |
| `sci_form::transport` | `pack_batch_arrow()`, `ChunkedIterator`, `WorkerTask`, `split_worker_tasks()` |
