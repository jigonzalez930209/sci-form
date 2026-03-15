# sci-form

**High-performance 3D molecular conformer generation and quantum-chemistry-inspired property computation**, written in pure Rust.

Generates chemically valid 3D coordinates from SMILES strings with RDKit-quality accuracy (ETKDGv2), and provides a full suite of computational chemistry tools: Extended Hückel Theory, electrostatic potentials, density of states, population analysis, molecular alignment, force field energy evaluation (UFF + MMFF94), partial charges, SASA, and materials framework assembly.

Native bindings for **Rust**, **Python**, **TypeScript/JavaScript (WASM)**, and a **cross-platform CLI**.

[![crates.io](https://img.shields.io/crates/v/sci-form)](https://crates.io/crates/sci-form)
[![PyPI](https://img.shields.io/pypi/v/sciforma)](https://pypi.org/project/sciforma)
[![npm](https://img.shields.io/npm/v/sci-form-wasm)](https://www.npmjs.com/package/sci-form-wasm)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

## Features

### Conformer Generation
- **ETKDGv2 Distance Geometry** — CSD torsion preferences (846 SMARTS patterns), matches RDKit accuracy
- **High Accuracy** — 0.00% heavy-atom RMSD > 0.5 Å vs RDKit on GDB-20 (2000 molecules, ensemble comparison)
- **Fast** — 60+ molecules/second in Rust, parallel batch processing via rayon
- **Full chemical coverage** — organics, stereocenters, macrocycles, fused rings, metals, halogens (He→Bi)

### Quantum Chemistry (EHT)
- **Extended Hückel Theory** — Wolfsberg-Helmholtz Hamiltonian, Löwdin orthogonalization, HOMO/LUMO gaps
- **Mulliken & Löwdin population analysis** — atomic orbital contributions per atom
- **Molecular dipole moment** — bond dipoles + lone-pair contributions in Debye
- **Volumetric orbital grids** — STO-3G basis, chunked evaluation for large molecules, marching cubes isosurfaces
- **Density of States** — total DOS + per-atom PDOS with Gaussian smearing, JSON export, MSE metric

### Electrostatics & Surface
- **Electrostatic Potential (ESP)** — Coulomb grid from Mulliken charges, color mapping (red/white/blue), `.cube` export
- **Parallel ESP** — rayon-accelerated grid computation (`parallel` feature)
- **Solvent-Accessible Surface Area (SASA)** — Shrake-Rupley algorithm, Fibonacci sphere, Bondi radii
- **Gasteiger-Marsili partial charges** — 6-iteration electronegativity equalization

### Force Fields
- **UFF** — Universal Force Field for 50+ element types (including transition metals Ti–Zn, Pd, Pt, Au)
- **MMFF94** — Merck force field (stretch, bend, torsion, van der Waals, electrostatics via 14-7 potential)
- **BFGS / L-BFGS minimizers** — dense BFGS for small molecules, L-BFGS for large systems

### Molecular Alignment
- **Kabsch alignment** — SVD-based optimal rotation, minimizes RMSD
- **Quaternion alignment** — Coutsias 2004 4×4 eigenproblem (faster for large molecules)
- **RMSD computation** — after optimal superposition

### Materials
- **Periodic unit cells** — lattice parameters (a, b, c, α, β, γ) to Cartesian tensor
- **Secondary Building Units (SBUs)** — node/linker topology with coordination sites
- **Framework assembly** — MOF-type crystal structure generation from SBUs + topology

### Platform
- **Multi-platform** — Rust lib, Python (PyO3), TypeScript/JS (WASM), CLI (Linux/macOS/Windows)
- **Zero runtime dependencies** — pure Rust, no C++ toolchain needed
- **SMILES + SMARTS** — full parser and substructure matcher; 60+ bracket elements

---

## Quick Start

### Rust

```toml
[dependencies]
sci-form = "0.2"
# For parallel batch + ESP:
sci-form = { version = "0.2", features = ["parallel"] }
```

```rust
use sci_form::{embed, compute_charges, compute_esp, compute_dos, compute_population};

// 3D conformer
let result = sci_form::embed("CCO", 42);
println!("Atoms: {}, Time: {:.1}ms", result.num_atoms, result.time_ms);

// Gasteiger–Marsili charges
let charges = sci_form::compute_charges("CCO").unwrap();
println!("Charges: {:?}", charges.charges);

// EHT → population analysis
let mol = sci_form::parse("CCO").unwrap();
let elements: Vec<u8> = /* ... */ vec![8, 6, 6, 1, 1, 1, 1, 1, 1];
let positions: Vec<[f64; 3]> = /* ... from result.coords */ vec![];
let pop = sci_form::compute_population(&elements, &positions).unwrap();
println!("HOMO: {:.3} eV", pop.homo_energy);

// ESP grid
let esp = sci_form::compute_esp(&elements, &positions, 0.5, 3.0).unwrap();

// DOS / PDOS
let dos = sci_form::compute_dos(&elements, &positions, 0.2, -20.0, 5.0, 200).unwrap();
println!("HOMO–LUMO gap: {:.3} eV", dos.homo_lumo_gap.unwrap_or(0.0));
```

→ [Full Rust API reference](https://jigonzalez930209.github.io/sci-form/api/rust) · [Guide](https://jigonzalez930209.github.io/sci-form/guide/rust)

---

### Python

```bash
pip install sciforma
```

```python
import sci_form

# 3D conformer
result = sci_form.embed("CCO")
print(f"Atoms: {result.num_atoms}, Time: {result.time_ms:.1f}ms")

# Batch
results = sci_form.embed_batch(["CCO", "c1ccccc1", "CC(=O)O"])

# Gasteiger charges
charges = sci_form.compute_charges("CCO")
print(charges["charges"])  # per-atom charges

# EHT + population analysis
elements = [8, 6, 6, 1, 1, 1, 1, 1, 1]
positions = [[0.0, 0.0, 0.0], ...]   # from result.coords
pop = sci_form.compute_population(elements, positions)
print(f"HOMO: {pop['homo_energy']:.3f} eV")

# DOS / PDOS
dos = sci_form.compute_dos(elements, positions, sigma=0.2, e_min=-20.0, e_max=5.0, n_points=200)
print(f"Gap: {dos['homo_lumo_gap']:.3f} eV")
```

→ [Full Python API reference](https://jigonzalez930209.github.io/sci-form/api/python) · [Guide](https://jigonzalez930209.github.io/sci-form/guide/python)

---

### TypeScript / JavaScript (WASM)

```bash
npm install sci-form-wasm
```

```typescript
import init, {
  embed, embed_coords_typed,
  compute_esp_grid_typed, compute_esp_grid_info,
  eht_calculate, eht_orbital_mesh, eht_orbital_grid_typed,
  compute_charges, compute_dos
} from 'sci-form-wasm';

await init();

// Conformer as JSON
const result = JSON.parse(embed('CCO', 42));
console.log(result.num_atoms);

// Typed-array conformer (faster, no JSON overhead)
const coords: Float64Array = embed_coords_typed('CCO', 42);

// EHT calculation
const eht = JSON.parse(eht_calculate('[6,8,6,1,1,1,1,1,1]', coords.toString(), 1.75));
console.log(`HOMO: ${eht.homo_energy} eV, LUMO: ${eht.lumo_energy} eV`);

// Orbital isosurface mesh
const mesh = JSON.parse(eht_orbital_mesh('[6,8,6,1,1,1,1,1,1]', coords.toString(), 1.75, 0, 0.02));
// mesh.vertices, mesh.normals, mesh.indices

// ESP grid (typed array)
const espData: Float64Array = compute_esp_grid_typed('CCO', 42, 0.5, 3.0);
const espInfo = JSON.parse(compute_esp_grid_info('CCO', 42, 0.5, 3.0));
console.log(`Grid: ${espInfo.nx}×${espInfo.ny}×${espInfo.nz}`);

// DOS
const dos = JSON.parse(compute_dos('[6,8,6,1,1,1,1,1,1]', coords.toString(), 0.2, -20.0, 5.0, 200));
```

→ [Full TypeScript API reference](https://jigonzalez930209.github.io/sci-form/api/typescript) · [Guide](https://jigonzalez930209.github.io/sci-form/guide/typescript)

---

### CLI

```bash
cargo install sci-form-cli
```

```bash
# Single molecule
sci-form embed "CCO" --format xyz

# Batch processing
sci-form batch -i molecules.smi -o output.sdf --format sdf --threads 8

# Parse only (no 3D)
sci-form parse "c1ccccc1"

# Gasteiger charges
sci-form charges "CCO"

# UFF energy
sci-form energy "CCO" --coords coords.json

# Version / features
sci-form info
```

Prebuilt binaries available at [GitHub Releases](https://github.com/jigonzalez930209/sci-form/releases):

| Platform | File |
|----------|------|
| Linux x86_64 | `sci-form-linux-x86_64` |
| Linux aarch64 | `sci-form-linux-aarch64` |
| macOS x86_64 | `sci-form-macos-x86_64` |
| macOS Apple Silicon | `sci-form-macos-aarch64` |
| Windows x86_64 | `sci-form-windows-x86_64.exe` |

→ [Full CLI reference](https://jigonzalez930209.github.io/sci-form/api/cli) · [Guide](https://jigonzalez930209.github.io/sci-form/guide/cli)

---

## Benchmark Results

### Conformer Generation — Diverse Molecules (131 molecules, all functional groups)

| Metric | Value |
|--------|-------|
| Parse success | 100% |
| Embed success | 97.7% |
| Geometry quality | 97.7% |
| Throughput | 60 mol/s |

### RDKit Comparison — Heavy-atom pairwise-distance RMSD

| Metric | Value |
|--------|-------|
| Average RMSD | 0.064 Å |
| Median RMSD | 0.011 Å |
| < 0.5 Å | 98.4% |
| < 0.3 Å | 94.4% |

### GDB-20 Ensemble (2000 molecules × 10 seeds vs 21 RDKit seeds)

| Metric | All-atom | Heavy-atom |
|--------|----------|------------|
| Avg RMSD | 0.035 Å | 0.018 Å |
| > 0.5 Å | 0.95% | **0.00%** |

---

## Module Overview

| Module | Description |
|--------|-------------|
| `sci_form::embed` | ETKDGv2 3D conformer generation from SMILES |
| `sci_form::embed_batch` | Parallel batch conformer generation (rayon) |
| `sci_form::parse` | SMILES → molecular graph |
| `sci_form::compute_charges` | Gasteiger-Marsili partial charges |
| `sci_form::compute_sasa` | Solvent-accessible surface area (SASA) |
| `sci_form::compute_population` | Mulliken & Löwdin population analysis (EHT) |
| `sci_form::compute_dipole` | Molecular dipole moment in Debye (EHT) |
| `sci_form::compute_esp` | Electrostatic potential grid (Mulliken charges) |
| `sci_form::compute_dos` | Density of states + PDOS (EHT orbital energies) |
| `sci_form::compute_rmsd` | RMSD after Kabsch alignment |
| `sci_form::compute_uff_energy` | UFF force field energy evaluation |
| `sci_form::create_unit_cell` | Periodic unit cell from lattice parameters |
| `sci_form::assemble_framework` | MOF-type framework assembly from SBUs |

### Sub-modules

| Module Path | Key API |
|-------------|---------|
| `sci_form::eht` | `solve_eht()`, `EhtResult`, `evaluate_orbital_on_grid_chunked()`, `marching_cubes()` |
| `sci_form::esp` | `compute_esp_grid_parallel()`, `esp_color_map()`, `esp_grid_to_colors()`, `export_cube()`, `read_cube()` |
| `sci_form::dos` | `compute_dos()`, `compute_pdos()`, `dos_mse()`, `export_dos_json()` |
| `sci_form::alignment` | `align_coordinates()`, `align_quaternion()`, `compute_rmsd()` |
| `sci_form::forcefield` | `build_uff_force_field()`, `Mmff94Builder::build()` |
| `sci_form::charges::gasteiger` | `gasteiger_marsili_charges()` |
| `sci_form::surface::sasa` | `compute_sasa()` |
| `sci_form::materials` | `UnitCell`, `Sbu`, `Topology`, `assemble_framework()` |
| `sci_form::transport` | `pack_batch_arrow()`, `ChunkedIterator`, `WorkerTask` |

---

## The Conformer Pipeline

sci-form implements ETKDGv2 (Experimental Torsion Knowledge Distance Geometry):

1. **SMILES Parsing** → Molecular graph with atoms, bonds, hybridization  
2. **Bounds Matrix** → 1-2, 1-3, 1-4, and VdW distance bounds from topology  
3. **Triangle Smoothing** → Floyd-Warshall triangle inequality enforcement  
4. **Distance Picking** → Random distances from smoothed bounds (MinstdRand)  
5. **Metric Matrix Embedding** → Eigendecomposition → 4D coordinates  
6. **Bounds Force Field** → BFGS minimization in 4D to satisfy distance constraints  
7. **Projection to 3D** → Drop lowest-variance dimension  
8. **ETKDG 3D Refinement** → Force field with CSD torsion preferences (846 patterns)  
9. **Validation** → Tetrahedral centers, planarity, double-bond geometry

See the [algorithm documentation](https://jigonzalez930209.github.io/sci-form/algorithm/overview) for mathematical derivations and diagrams.

---

## Building from Source

```bash
# Library + CLI
cargo build --release

# Python bindings
cd crates/python && pip install maturin && maturin develop --release

# WASM bindings
cd crates/wasm && wasm-pack build --target web --release

# With parallel feature
cargo build --release --features parallel
```

---

## Testing

```bash
# All unit tests
cargo test --lib

# Integration — diverse molecules (97.7% pass rate)
cargo test --release --test test_diverse_molecules -- --nocapture

# Integration — gradient correctness
cargo test --release --test test_gradient_check -- --nocapture

# Lint & format
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
```

---

## Releasing a New Version

Use the provided bump script. It updates all version strings, commits, tags, and pushes:

```bash
# Auto-increment patch (0.1.7 → 0.1.8)
./scripts/bump_version.sh

# Set a specific version
./scripts/bump_version.sh 0.2.0
```

This updates versions in:
- `Cargo.toml` (root lib)
- `crates/cli/Cargo.toml`
- `crates/python/Cargo.toml`
- `crates/wasm/Cargo.toml`
- `crates/python/pyproject.toml`
- `crates/wasm/pkg/package.json`
- `pkg/package.json` & `pkg-node/package.json`

Then creates a `vX.Y.Z` git tag, which triggers the [release workflow](.github/workflows/release.yml) to publish to **crates.io**, **PyPI**, and **npm** automatically.

**Required repository secrets:**

| Secret | Used for |
|--------|----------|
| `CARGO_REGISTRY_TOKEN` | Publishing to crates.io |
| `PYPI_API_TOKEN` | Publishing to PyPI (`sciforma`) |
| `NPM_TOKEN` | Publishing to npm (`sci-form-wasm`) — must be a **Granular Automation token** |

---

## License

MIT
