# sci-form

**High-performance computational chemistry in pure Rust**, spanning 3D conformer generation, semi-empirical quantum methods, spectroscopy, descriptors, and materials workflows.

Starting from SMILES, sci-form generates chemically valid 3D coordinates with ETKDGv2-quality distance geometry and exposes a broad scientific toolkit around them: EHT, PM3/PM3(tm), GFN0/GFN1/GFN2-xTB, HF-3c, ANI potentials, electrostatic potential grids, density of states, dipoles, population analysis, solvation, stereochemistry, fingerprints, clustering, spectroscopy, and periodic/materials utilities.

The same project ships native surfaces for **Rust**, **Python**, **TypeScript/JavaScript (WASM)**, and a **cross-platform CLI**.

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

### Electronic Structure & Spectroscopy
- **Extended Hückel Theory (EHT)** — orbital energies, population analysis, dipoles, ESP, DOS/PDOS, orbital grids, and meshes
- **PM3 / PM3(tm)** — NDDO SCF, heat of formation, HOMO/LUMO gaps, Mulliken charges, transition-metal support
- **GFN0 / GFN1 / GFN2-xTB** — SCC tight-binding with shell charges, multipoles, dispersion, and halogen-bond terms
- **HF-3c + CISD** — minimal-basis Hartree-Fock with D3/gCP/SRB corrections and excited-state workflows
- **UV-Vis, IR, and NMR** — sTDA spectra, vibrational analysis, broadened IR, HOSE-code shifts, and J-couplings

### Molecular Modeling & Screening
- **Charges, SASA, solvation, and reactivity** — Gasteiger, EEQ/ALPB/D4 experimental modules, GB/non-polar solvation, Fukui descriptors, pKa heuristics
- **Force fields** — UFF and MMFF94 energies, gradients, and geometry optimization helpers
- **Descriptors & ML** — WHIM, RDF, GETAWAY, topological descriptors, LogP/ESOL/druglikeness, Random Forest, Gradient Boosting
- **Fingerprints & clustering** — SSSR, ECFP/Morgan, Tanimoto similarity, RMSD matrices, and Butina clustering
- **Stereochemistry** — CIP priorities, R/S, E/Z, helical chirality, atropisomeric axes, and prochiral analysis

### Materials, Periodic Systems, and Transport
- **Crystallography & frameworks** — unit cells, 230 space groups, framework assembly, and periodic geometry optimization
- **Periodic graph analysis** — periodic bonding, hapticity detection, and metallocene-aware workflows
- **Transport helpers** — Arrow-style conformer packing, worker-task splitting, and batch-oriented browser workflows

### Platform & Performance
- **Multi-platform** — Rust lib, Python (PyO3), TypeScript/JS (WASM), CLI (Linux/macOS/Windows)
- **Parallel by default** — rayon-backed batch and property kernels with CPU-aware work scheduling
- **WebGPU validation tooling** — browser validation harness and benchmark flows for accelerated WASM volumetric workloads
- **Pure Rust** — no C++ runtime dependency for the library, bindings, or CLI

---

## Quick Start

### Rust

```toml
[dependencies]
sci-form = "0.10"
# For parallel batch + ESP:
sci-form = { version = "0.10", features = ["parallel"] }
```

```rust
use sci_form::{compute_dipole, compute_pm3, compute_population, embed};

// 3D conformer
let result = embed("CCO", 42);
assert!(result.error.is_none(), "embedding failed: {:?}", result.error);

let positions: Vec<[f64; 3]> = result
    .coords
    .chunks_exact(3)
    .map(|chunk| [chunk[0], chunk[1], chunk[2]])
    .collect();

let pop = compute_population(&result.elements, &positions).unwrap();
println!("HOMO: {:.3} eV", pop.homo_energy);

let dipole = compute_dipole(&result.elements, &positions).unwrap();
println!("Dipole: {:.3} D", dipole.magnitude);

let pm3 = compute_pm3(&result.elements, &positions).unwrap();
println!("PM3 gap: {:.3} eV (converged: {})", pm3.gap, pm3.converged);
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
result = sci_form.embed("CCO", seed=42)
if not result.is_ok():
    raise RuntimeError(result.error)

print(f"Atoms: {result.num_atoms}, Time: {result.time_ms:.1f}ms")

# Batch
results = sci_form.embed_batch(["CCO", "c1ccccc1", "CC(=O)O"])

# Population + dipole from the embedded geometry
population = sci_form.population(result.elements, result.coords)
print(f"HOMO: {population.homo_energy:.3f} eV")

dipole = sci_form.dipole(result.elements, result.coords)
print(f"Dipole: {dipole.magnitude:.3f} {dipole.unit}")

# Semi-empirical screening
pm3 = sci_form.pm3_calculate(result.elements, result.coords)
xtb = sci_form.xtb_calculate(result.elements, result.coords)
print(f"PM3 gap: {pm3.gap:.3f} eV | xTB gap: {xtb.gap:.3f} eV")
```

→ [Full Python API reference](https://jigonzalez930209.github.io/sci-form/api/python) · [Guide](https://jigonzalez930209.github.io/sci-form/guide/python)

---

### TypeScript / JavaScript (WASM)

```bash
npm install sci-form-wasm
```

```typescript
import init, {
  analyze_stereo,
  compute_pm3,
  compute_population,
  compute_xtb,
  embed,
} from 'sci-form-wasm';

await init();

const result = JSON.parse(embed('CCO', 42));
if (result.error) throw new Error(result.error);

const elements = JSON.stringify(result.elements);
const coords = JSON.stringify(result.coords);

const population = JSON.parse(compute_population(elements, coords));
console.log(`HOMO: ${population.homo_energy.toFixed(3)} eV`);

const pm3 = JSON.parse(compute_pm3(elements, coords));
const xtb = JSON.parse(compute_xtb(elements, coords));
console.log(`PM3 gap: ${pm3.gap.toFixed(3)} eV | xTB gap: ${xtb.gap.toFixed(3)} eV`);

const chiral = JSON.parse(embed('C(F)(Cl)(Br)I', 42));
const stereo = JSON.parse(analyze_stereo('C(F)(Cl)(Br)I', JSON.stringify(chiral.coords)));
console.log(`Stereocenters: ${stereo.n_stereocenters}`);
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

# Run property workflows on known coordinates
sci-form population "[8,1,1]" "[0,0,0,0.96,0,0,-0.24,0.93,0]"
sci-form pm3 "[8,1,1]" "[0,0,0,0.96,0,0,-0.24,0.93,0]"
sci-form xtb "[8,1,1]" "[0,0,0,0.96,0,0,-0.24,0.93,0]"

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

## Experimental Engine

The repository now includes an isolated experimental quantum-chemistry stack under `sci_form::experimental_2`.

### Phase Status

| Phase | Status | Scope |
|------|--------|-------|
| Phase 1 | Complete | GPU infrastructure scaffolding, aligned types, CPU fallback, WGSL interface stubs |
| Phase 2 | Complete | Gaussian basis, overlap/kinetic/nuclear/core matrices, ERIs, validation helpers |
| Phase 3 | Complete | RHF SCF loop, Löwdin orthogonalization, DIIS, Mulliken analysis, gradients, optimizer, parallel ERI |
| Phase 4 | Prototype complete | Experimental sTDA UV-Vis, GIAO-style NMR shielding, Hessian/IR workflows |
| Phase 5 | Prototype complete | Orbital grid evaluation, marching cubes, isosurface generation, GPU-ready rendering path |

### Current Practical Status

- Stable production APIs remain unchanged; the experimental work is isolated in `experimental_2`
- Real acceleration today is CPU parallelism via rayon in the ERI build path
- GPU execution is not enabled yet; `phase1_gpu_infrastructure` currently exposes a CPU fallback plus WGSL-ready interfaces
- The main known scientific limitation is absolute RHF/STO-3G total energy scaling in the experimental engine; comparative gaps and regression behavior are the primary validation target today

### Validation Snapshot

The experimental stack is covered by dedicated regression suites:

```bash
# Build all library + test targets
cargo check --tests

# Base experimental regression suite
cargo test --test regression test_experimental_comparison -- --nocapture

# Extended complex-molecule battery (fast active tests)
cargo test --test regression test_extended_molecules -- --nocapture

# Heavy experimental benchmarks and long-running tests
cargo test --release --test regression test_extended_molecules -- --include-ignored
```

Current verified results on this repository state:

| Command | Result |
|--------|--------|
| `cargo check --tests` | passes |
| `cargo test --test regression test_experimental_comparison` | **54 passed, 0 failed** |
| `cargo test --test regression test_extended_molecules` | **14 passed, 0 failed, 7 ignored** |

More detailed coverage notes live in [TESTING.md](TESTING.md), and the broader project plan remains in [ROADMAP.md](ROADMAP.md).

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
cd crates/wasm && ./build.sh --web-only --web-features "parallel experimental-gpu"

# With parallel feature
cargo build --release --features parallel
```

---

## Testing

```bash
# All unit tests
cargo test --lib

# Smoke battery (CI gate)
cargo test --release --test ci -- --nocapture

# Full integration suites
cargo test --release --test regression -- --nocapture
cargo test --release --test analysis -- --nocapture
cargo test --release --test debug -- --nocapture
cargo test --release --test benchmarks -- --nocapture

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
