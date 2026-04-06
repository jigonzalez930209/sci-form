# sci-form

**High-performance 3D molecular conformer generation and quantum-chemistry-inspired property computation**, written in pure Rust.

Generates chemically valid 3D coordinates from SMILES strings with RDKit-quality accuracy (ETKDGv2), and provides a full suite of computational chemistry tools: Extended Hückel Theory, semi-empirical quantum chemistry (PM3, PM3(tm), GFN0/GFN1/GFN2-xTB), ab-initio methods (HF-3c, UHF/ROHF, CISD), neural network potentials (ANI-2x, ANI-TM), electrostatic potentials, density of states, population analysis, molecular alignment, force field energy evaluation (UFF + MMFF94), partial charges, SASA, stereochemistry, NMR/IR/UV-Vis spectroscopy, CIF crystallography, and materials framework assembly.

Native bindings for **Rust**, **Python**, **TypeScript/JavaScript (WASM)**, and a **cross-platform CLI**.

[![crates.io](https://img.shields.io/crates/v/sci-form)](https://crates.io/crates/sci-form)
[![PyPI](https://img.shields.io/pypi/v/sciforma)](https://pypi.org/project/sciforma)
[![npm](https://img.shields.io/npm/v/sci-form-wasm)](https://www.npmjs.com/package/sci-form-wasm)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

See also: [CHANGELOG.md](CHANGELOG.md) · [TESTING.md](TESTING.md) · [ROADMAP_REACTION_DYNAMICS_3D.md](ROADMAP_REACTION_DYNAMICS_3D.md)

## Features

### Conformer Generation
- **ETKDGv2 Distance Geometry** — CSD torsion preferences (846 SMARTS patterns), matches RDKit accuracy
- **High Accuracy** — 0.00% heavy-atom RMSD > 0.5 Å vs RDKit on GDB-20 (2000 molecules, ensemble comparison)
- **Fast** — 60+ molecules/second in Rust, parallel batch processing via rayon
- **Full chemical coverage** — organics, stereocenters, macrocycles, fused rings, metals, halogens (He→Bi)

### Quantum Chemistry (EHT)
- **Extended Hückel Theory** — Wolfsberg-Helmholtz Hamiltonian, Löwdin orthogonalization, HOMO/LUMO gaps
- **Mulliken & Löwdin population analysis** — atomic orbital contributions per atom (parallel, Z=1–86 including lanthanides)
- **Molecular dipole moment** — bond dipoles + lone-pair contributions in Debye
- **Volumetric orbital grids** — STO-3G basis, chunked evaluation for large molecules, marching cubes isosurfaces
- **Density of States** — total DOS + per-atom PDOS with Gaussian smearing, JSON export, MSE metric

### Semi-Empirical & Tight-Binding Methods
- **PM3 / PM3(tm)** — NDDO semi-empirical SCF with Gaussian core-core corrections, heat of formation, HOMO/LUMO; transition metal support (Ti–Au)
- **GFN0/GFN1/GFN2-xTB** — Tight-binding DFT family with Broyden SCC mixing, shell-resolved charges, D3/D4 dispersion, halogen bonding
- **EEQ Charges** — Extended Electronegativity Equalization with improved Gaussian damping

### Ab-Initio Methods
- **HF-3c** — Minimal-basis Hartree-Fock with D3-BJ dispersion, gCP, and SRB corrections
- **UHF / ROHF** — Unrestricted and restricted open-shell Hartree-Fock SCF with configurable level shift and damping
- **CISD** — Configuration Interaction Singles+Doubles for excited states
- **AO→MO Integral Transform** — 4-index integral transform with 4-fold symmetry for post-HF methods

### Spectroscopy
- **UV-Vis** — sTDA-xTB vertical excitations, Gaussian/Lorentzian broadening
- **IR** — Numerical Hessian vibrational analysis, dipole intensities, thermochemistry (RRHO), peak assignment
- **NMR** — Chemical shifts via HOSE codes, J-coupling (Karplus 2J–5J including long-range), ensemble averaging

### Machine Learning
- **ANI-2x / ANI-TM** — Neural network potentials with analytical gradients (24 elements including transition metals)
- **ML Properties** — LogP, molar refractivity, solubility, Lipinski Ro5, druglikeness
- **3D Descriptors** — WHIM, RDF, GETAWAY molecular descriptors
- **ML Models** — Decision Trees, Random Forest, Gradient Boosting with cross-validation

### Experimental RHF Engine
- **Isolated `experimental_2` namespace** — next-generation Roothaan-Hall RHF engine without affecting stable APIs
- **Analytical Gaussian integrals** — overlap, kinetic, nuclear attraction, core Hamiltonian, and two-electron ERIs
- **SCF with DIIS + parallel ERI** — deterministic RHF/STO-3G workflow with rayon acceleration for the expensive $O(N^4)$ step
- **Experimental spectroscopy stack** — prototype sTDA UV-Vis, GIAO-like NMR shielding, and semi-numerical IR/Hessian workflows
- **GPU-oriented architecture** — WGSL shader stubs, orbital grid evaluation, and marching-cubes rendering pipeline prepared for future `wgpu` enablement

### Electrostatics & Surface
- **Electrostatic Potential (ESP)** — Coulomb grid from Mulliken charges, color mapping (red/white/blue), `.cube` export
- **Parallel ESP** — rayon-accelerated grid computation (`parallel` feature)
- **Solvent-Accessible Surface Area (SASA)** — Shrake-Rupley algorithm, Fibonacci sphere, Bondi radii, parallelized per-atom evaluation
- **Gasteiger-Marsili partial charges** — 6-iteration electronegativity equalization

### Parallel Computation
- **Automatic rayon thread pool** — all compute functions (DOS, ESP SASA, population, dipole, etc.) parallelized with work-stealing queue
- **Zero-configuration** — no API changes needed; parallelization enabled by default via `parallel` feature
- **Intra-library parallelism** — grid point loops, per-atom evaluation, force field term accumulation all use `par_iter()`
- **CPU-aware workload balancing** — handles small molecules (sequential) and large molecules (parallel) automatically

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
- **230 space groups** — full ITC space group library with equivalent position generation
- **CIF import/export** — parse and write CIF 1.1 crystallographic files with uncertainty notation support

### Platform
- **Multi-platform** — Rust lib, Python (PyO3), TypeScript/JS (WASM), CLI (Linux/macOS/Windows)
- **Zero runtime dependencies** — pure Rust, no C++ toolchain needed
- **SMILES, SMARTS, SMIRKS** — full parser, substructure matcher, and reaction transforms (multi-component); 60+ bracket elements

---

## Quick Start

### Rust

```toml
[dependencies]
sci-form = "0.13"
# For parallel batch + ESP:
sci-form = { version = "0.13", features = ["parallel"] }
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
| `sci_form::compute_pm3` | PM3/PM3(tm) semi-empirical SCF |
| `sci_form::compute_xtb` | GFN0-xTB tight-binding with Broyden SCC |
| `sci_form::compute_uhf` | Unrestricted Hartree-Fock SCF |
| `sci_form::compute_rohf` | Restricted open-shell Hartree-Fock SCF |
| `sci_form::compute_uff_energy` | UFF force field energy evaluation |
| `sci_form::compute_mmff94_energy` | MMFF94 force field energy evaluation |
| `sci_form::parse_cif` | CIF crystallographic file import |
| `sci_form::write_cif` | CIF crystallographic file export |
| `sci_form::create_unit_cell` | Periodic unit cell from lattice parameters |
| `sci_form::assemble_framework` | MOF-type framework assembly from SBUs |

### Sub-modules

| Module Path | Key API |
|-------------|---------|
| `sci_form::eht` | `solve_eht()`, `EhtResult`, `evaluate_orbital_on_grid_chunked()`, `marching_cubes()` |
| `sci_form::pm3` | `compute_pm3()`, PM3/PM3(tm) SCF, Gaussian core-core corrections |
| `sci_form::xtb` | `compute_xtb()`, GFN0 SCC with Broyden mixing; `xtb::gfn1`, `xtb::gfn2` |
| `sci_form::scf::uhf` | `run_uhf()`, `run_rohf()`, open-shell SCF with spin contamination analysis |
| `sci_form::hf` | `solve_hf3c()`, CISD, `mo_transform::ao_to_mo_transform()` |
| `sci_form::esp` | `compute_esp_grid_parallel()`, `esp_color_map()`, `esp_grid_to_colors()`, `export_cube()`, `read_cube()` |
| `sci_form::dos` | `compute_dos()`, `compute_pdos()`, `dos_mse()`, `export_dos_json()` |
| `sci_form::alignment` | `align_coordinates()`, `align_quaternion()`, `compute_rmsd()` |
| `sci_form::forcefield` | `build_uff_force_field()`, `Mmff94Builder::build()` |
| `sci_form::charges::gasteiger` | `gasteiger_marsili_charges()` |
| `sci_form::charges_eeq` | `compute_eeq_charges()` — improved Gaussian damping |
| `sci_form::surface::sasa` | `compute_sasa()` |
| `sci_form::materials` | `UnitCell`, `Sbu`, `Topology`, `assemble_framework()` |
| `sci_form::materials::cif` | `parse_cif()`, `write_cif()`, `CifStructure`, `CifAtomSite` |
| `sci_form::nmr` | `predict_nmr_shifts()`, `predict_j_couplings()` (2J–5J) |
| `sci_form::smirks` | `apply_smirks()`, `apply_smirks_multi()`, `parse_smirks()` |
| `sci_form::stereo` | `analyze_stereo()` — R/S, E/Z, helical, atropisomeric |
| `sci_form::population` | `compute_population()`, `compute_bond_orders()`, parallel variants |
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
