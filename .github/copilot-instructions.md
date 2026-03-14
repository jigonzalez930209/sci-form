# sci-form — Agent Instructions

**sci-form** is a Rust library (v0.3.1) for computational chemistry: 3D conformer generation, semi-empirical quantum chemistry (EHT), molecular properties, and crystallographic materials. It is available as a Rust crate, Python package, WebAssembly/TypeScript module, and a CLI binary.

---

## Quick Reference — What the Library Does

| Module | Function |
|--------|----------|
| Conformer (ETKDG) | SMILES → 3D coordinates (ETKDGv2) |
| EHT | Semi-empirical Hückel MO theory, orbital grids, isosurfaces |
| Charges | Gasteiger-Marsili partial charges |
| SASA | Solvent-accessible surface area (Shrake-Rupley) |
| Population | Mulliken & Löwdin population analysis |
| Dipole | Molecular dipole moment (Debye) |
| ESP | Electrostatic potential grid |
| DOS | Density of states with Gaussian smearing |
| Alignment | Kabsch SVD & quaternion RMSD alignment |
| Force Field | UFF and MMFF94 energy evaluation |
| Materials | Unit cell, MOF framework assembly |
| Transport | Arrow columnar batch + Web Worker task splitting |

---

## Rust

### Setup

```toml
[dependencies]
sci-form = "0.3"
# parallel batch:
sci-form = { version = "0.3", features = ["parallel"] }
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

// Materials
sci_form::create_unit_cell(a,b,c,alpha,beta,gamma: f64) -> materials::UnitCell
sci_form::assemble_framework(node: &Sbu, linker: &Sbu, topology: &Topology, cell: &UnitCell) -> CrystalStructure
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
    population, dipole, dos, rmsd, uff_energy,
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
from sci_form import embed, population, dipole

conf = embed("CCO", seed=42)
pop = population(conf.elements, conf.coords)
print(f"HOMO: {pop.homo_energy:.3f} eV, Gap: {pop.homo_energy - pop.lumo_energy:.3f} eV")

d = dipole(conf.elements, conf.coords)
print(f"|μ| = {d.magnitude:.3f} D")
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

// Materials
create_unit_cell(a,b,c,alpha,beta,gamma: number): string
assemble_framework(topology: string, metal: number, geometry: string, lattice_a: number, supercell: number): string

// Transport
pack_batch_arrow(results_json: string): string
split_worker_tasks(smiles_json: string, n_workers: number, seed: number): string
estimate_workers(n_items: number, max_workers: number): number
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
- **EHT requires 3D coordinates** — run `embed` first, then pass `elements` + `coords`.
- **ESP/DOS are slow for large molecules** — prefer `spacing ≥ 0.3` Å for interactive use.
- **`topology` strings**: `"pcu"` (cubic), `"dia"` (diamond), `"sql"` (square lattice).
- **`geometry` strings**: `"linear"`, `"trigonal"`, `"tetrahedral"`, `"square_planar"`, `"octahedral"`.

---

## Units

| Quantity | Unit |
|----------|------|
| Coordinates | Å (ångströms) |
| Energies (EHT, DOS, HOMO/LUMO) | eV |
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
