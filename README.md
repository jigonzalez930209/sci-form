# sci-form

**Comprehensive computational chemistry library in pure Rust** — from 3D conformer generation and semi-empirical quantum chemistry to neural network potentials, spectroscopy, cheminformatics, and crystallographic materials.

Native bindings for **Rust**, **Python**, **TypeScript/JavaScript (WASM)**, and a **cross-platform CLI**.

[![crates.io](https://img.shields.io/crates/v/sci-form)](https://crates.io/crates/sci-form)
[![PyPI](https://img.shields.io/pypi/v/sciforma)](https://pypi.org/project/sciforma)
[![npm](https://img.shields.io/npm/v/sci-form-wasm)](https://www.npmjs.com/package/sci-form-wasm)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

---

## Features

### Conformer Generation
- **ETKDGv2 Distance Geometry** — CSD torsion preferences (846 SMARTS patterns), matches RDKit accuracy
- **High Accuracy** — 0.00% heavy-atom RMSD > 0.5 Å vs RDKit on GDB-20 (2000 molecules)
- **Fast** — 60+ molecules/second, parallel batch via rayon
- **Full coverage** — organics, stereocenters, macrocycles, fused rings, metals, halogens (He→Bi)

### Quantum Chemistry
- **EHT** — Extended Hückel Theory, orbital grids, marching cubes isosurfaces, analytical gradients, geometry optimization, k-path band structure for periodic systems
- **PM3** — NDDO semi-empirical SCF: heat of formation, HOMO/LUMO, thermochemistry
- **PM3(tm)** — PM3 with transition-metal parameters (Ti–Au)
- **GFN0-xTB** — tight-binding SCC (25 elements including transition metals)
- **GFN1-xTB** — shell-resolved charges + D3-BJ dispersion
- **GFN2-xTB** — multipole electrostatics + D4 dispersion + halogen-bond correction
- **HF-3c** — minimal-basis Hartree-Fock with D3, gCP, and SRB corrections
- **CISD** — Configuration Interaction Singles+Doubles for excited states (via HF-3c)

### Neural Network Potentials
- **ANI-2x** — neural network potential for H, C, N, O, F, S, Cl (energies + forces)
- **ANI-TM** — ANI extended to 24 elements including Ti, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ru, Pd, Ag, Pt, Au

### Spectroscopy
- **UV-Vis (sTDA)** — oscillator strengths, Gaussian/Lorentzian broadening, multi-peak output
- **IR** — numerical Hessian, 3N normal modes, ZPVE, km/mol intensities, Lorentzian/Gaussian broadening, functional-group peak assignment
- **NMR** — ¹H/¹³C/¹⁹F/³¹P chemical shifts (HOSE codes), J-coupling (Karplus equation), spectrum simulation, ensemble-averaged J-couplings

### Molecular Properties
- **Charges** — Gasteiger-Marsili, EEQ (geometry-dependent)
- **SASA** — Shrake-Rupley solvent-accessible surface area
- **Population analysis** — Mulliken, Löwdin, NPA, NBO, bond orders
- **Dipole moment** — EHT-based, in Debye
- **ESP** — Coulomb grid, parallel evaluation, `.cube` export
- **DOS** — total + PDOS, multi-method (EHT/PM3/xTB/GFN1/GFN2/HF-3c)
- **Dispersion** — D4 energy and gradients
- **Reactivity** — Fukui descriptors, frontier orbital descriptors, empirical pKa

### Solvation
- **Non-polar SASA solvation** — probe-radius-based cavity energy
- **Generalized Born (HCT)** — electrostatic solvation (configurable dielectric)
- **ALPB** — Analytical Linearised Poisson-Boltzmann solvation

### Machine Learning
- **2D descriptors** — MW, Wiener index, Balaban J, FSP3, 17+ topology-based descriptors
- **3D descriptors** — WHIM (5 weightings), RDF, GETAWAY
- **Property prediction** — LogP, MR, solubility, Lipinski Ro5, druglikeness
- **ML models** — Random Forest (with variance estimation), Gradient Boosting, cross-validation

### Cheminformatics
- **Stereochemistry** — CIP priorities, R/S tetrahedral, E/Z double bonds, atropisomeric M/P axes, helical chirality
- **SSSR** — Smallest Set of Smallest Rings (Horton's algorithm)
- **ECFP fingerprints** — Extended-Connectivity Fingerprints (Morgan algorithm), Tanimoto similarity
- **Clustering** — Butina (Taylor-Butina) RMSD clustering, diversity filtering
- **SMIRKS** — reaction transform engine (atom-mapped SMIRKS patterns)
- **Canonical SMILES** — deterministic canonical output

### Materials & Crystallography
- **Unit cells** — lattice parameters (a, b, c, α, β, γ) to Cartesian tensor
- **230 ITC space groups** — full symmetry operations and equivalent positions
- **Framework assembly** — MOF-type crystal structures (pcu, dia, sql topologies)
- **Framework geometry optimization** — BFGS + steepest descent with PBC
- **Periodic systems** — PBC molecular graphs with bond tolerance, hapticity/metallocene detection

### Force Fields & Alignment
- **UFF** — Universal Force Field (50+ element types including transition metals)
- **MMFF94** — Merck force field (Halgren 14-7 vdW, quartic/cubic bends, Fourier torsions)
- **Kabsch / quaternion alignment** — SVD-based (Coutsias 2004) RMSD minimization
- **Molecular dynamics** — trajectory computation

### Transport & Parallel
- **Arrow columnar batch** — Apache Arrow transport for batch results
- **Web Worker task splitting** — parallelisable SMILES batches for WASM workers
- **Rayon parallelism** — grid loops, per-atom evaluation, force field accumulation via `parallel` feature

### Experimental / Research (feature-gated)
- **GPU** (`experimental-gpu`) — MMFF94 non-bonded, EEQ GPU kernels, WGSL shaders
- **Alpha** — CGA (conformal geometric algebra), GSM (growing string method), SDR (spectral dimensionality reduction)
- **Beta** — KPM, MBH, CPM, RandNLA, Riemannian optimization

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
let conf = sci_form::embed("CCO", 42);
assert!(conf.error.is_none());
let pos: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

// PM3 semi-empirical SCF
let pm3 = sci_form::compute_pm3(&conf.elements, &pos).unwrap();
println!("PM3 HOF: {:.2} kcal/mol, gap: {:.3} eV", pm3.heat_of_formation, pm3.gap);

// GFN2-xTB
let gfn2 = sci_form::xtb::gfn2::solve_gfn2(&conf.elements, &pos).unwrap();
println!("GFN2 gap: {:.3} eV, D4: {:.4} eV", gfn2.gap, gfn2.dispersion_energy);

// Population analysis
let pop = sci_form::compute_population(&conf.elements, &pos).unwrap();
println!("HOMO: {:.3} eV, gap: {:.3} eV", pop.homo_energy, pop.homo_lumo_gap);

// Stereochemistry
let stereo = sci_form::analyze_stereo("C(F)(Cl)(Br)I", &conf.coords).unwrap();
println!("Stereocenters: {}", stereo.n_stereocenters);

// ECFP fingerprints + Tanimoto
let fp1 = sci_form::compute_ecfp("c1ccccc1", 2, 2048).unwrap();
let fp2 = sci_form::compute_ecfp("Cc1ccccc1", 2, 2048).unwrap();
println!("Tanimoto: {:.3}", sci_form::compute_tanimoto(&fp1, &fp2));

// ML property prediction (no 3D needed)
let desc = sci_form::compute_ml_descriptors(&conf.elements, &conf.bonds, &[], &[]);
let props = sci_form::predict_ml_properties(&desc);
println!("LogP: {:.2}, passes Lipinski: {}", props.logp, props.lipinski_passes);

// 3D descriptors
use sci_form::ml::whim::{compute_whim, WhimWeighting};
let whim = compute_whim(&conf.elements, &pos, WhimWeighting::Mass);
println!("WHIM anisotropy: {:.3}", whim.anisotropy);

// Space groups
use sci_form::materials::space_groups::space_group_by_number;
let sg = space_group_by_number(225).unwrap(); // Fm-3m
println!("{}: {} ({})", sg.number, sg.symbol, sg.crystal_system);

// Parallel batch (requires features = ["parallel"])
let config = sci_form::ConformerConfig { seed: 42, num_threads: 0 };
let results = sci_form::embed_batch(&["CCO", "c1ccccc1", "CC(=O)O"], &config);
```

→ [Full Rust API reference](https://jigonzalez930209.github.io/sci-form/api/rust) · [Guide](https://jigonzalez930209.github.io/sci-form/guide/rust)

---

### Python

```bash
pip install sciforma
```

```python
from sci_form import embed, pm3_calculate, xtb_calculate, stereo_analysis, ecfp, tanimoto, ml_predict

conf = embed("CCO", seed=42)
pos = conf.get_positions()  # list of (x, y, z) tuples

# PM3 semi-empirical
pm3 = pm3_calculate(conf.elements, conf.coords)
print(f"PM3 HOF: {pm3.heat_of_formation:.2f} kcal/mol, converged: {pm3.converged}")

# GFN-xTB
xtb = xtb_calculate(conf.elements, conf.coords)
print(f"xTB gap: {xtb.gap:.3f} eV")

# Stereochemistry
stereo = stereo_analysis("C(F)(Cl)(Br)I")
print(f"Stereocenters: {stereo.n_stereocenters}")

# Fingerprints + similarity
fp1 = ecfp("c1ccccc1", radius=2, n_bits=2048)
fp2 = ecfp("Cc1ccccc1", radius=2, n_bits=2048)
print(f"Tanimoto: {tanimoto(fp1, fp2):.3f}")

# ML property prediction (no 3D needed)
props = ml_predict("CCO")
print(f"LogP: {props.logp:.2f}, passes Lipinski: {props.lipinski_passes}")

# IR spectrum
from sci_form import ir_spectrum
ir = ir_spectrum(conf.elements, conf.coords, method="pm3")

# NMR shifts
from sci_form import nmr_shift
nmr = nmr_shift("CCO")

# Solvation
from sci_form import nonpolar_solvation, gb_solvation
from sci_form import charges as compute_charges
solv = nonpolar_solvation(conf.elements, conf.coords)
print(f"Non-polar solvation: {solv.energy_kcal_mol:.2f} kcal/mol")

# Parallel batch
from sci_form import embed_batch
results = embed_batch(["CCO", "c1ccccc1", "CC(=O)O"], seed=42, num_threads=4)
```

→ [Full Python API reference](https://jigonzalez930209.github.io/sci-form/api/python) · [Guide](https://jigonzalez930209.github.io/sci-form/guide/python)

---

### TypeScript / JavaScript (WASM)

```bash
npm install sci-form-wasm
```

```typescript
import init, {
  embed, compute_pm3, compute_xtb, compute_ml_properties,
  analyze_stereo, compute_ecfp, compute_tanimoto,
  compute_ir_spectrum, compute_vibrational_analysis,
  eht_calculate, eht_orbital_mesh,
} from 'sci-form-wasm';

await init();

const conf = JSON.parse(embed('CCO', 42));
const elements = JSON.stringify(conf.elements);
const coords   = JSON.stringify(conf.coords);

// PM3 semi-empirical
const pm3 = JSON.parse(compute_pm3(elements, coords));
console.log(`PM3 HOF: ${pm3.heat_of_formation.toFixed(2)} kcal/mol, gap: ${pm3.gap.toFixed(3)} eV`);

// GFN-xTB
const xtb = JSON.parse(compute_xtb(elements, coords));
console.log(`xTB energy: ${xtb.total_energy.toFixed(4)} eV, converged: ${xtb.converged}`);

// ML properties (no 3D needed)
const props = JSON.parse(compute_ml_properties('CCO'));
console.log(`LogP: ${props.logp.toFixed(2)}, Lipinski: ${props.lipinski_passes}`);

// Stereochemistry
const stereo = JSON.parse(analyze_stereo('C(F)(Cl)(Br)I', '[]'));
console.log(`Stereocenters: ${stereo.n_stereocenters}`);

// EHT orbital isosurface (Three.js-ready)
const eht = JSON.parse(eht_calculate(elements, coords, 1.75));
const mesh = JSON.parse(eht_orbital_mesh(elements, coords, eht.homo_index, 0.2, 0.02));
// mesh.vertices, mesh.normals, mesh.indices → BufferGeometry

// ECFP fingerprints + Tanimoto
const fp1 = JSON.parse(compute_ecfp('c1ccccc1', 2, 2048));
const fp2 = JSON.parse(compute_ecfp('Cc1ccccc1', 2, 2048));
const sim  = JSON.parse(compute_tanimoto(JSON.stringify(fp1), JSON.stringify(fp2)));
console.log(`Tanimoto: ${sim.tanimoto.toFixed(3)}`);
```

→ [Full TypeScript API reference](https://jigonzalez930209.github.io/sci-form/api/typescript) · [Guide](https://jigonzalez930209.github.io/sci-form/guide/typescript)

---

### CLI

```bash
cargo install sci-form-cli
# or download a prebuilt binary from GitHub Releases
```

```bash
# 3D conformer
sci-form embed "CCO" --format xyz

# Batch from file (one SMILES per line)
sci-form batch molecules.smi --threads 8

# PM3 semi-empirical energy
sci-form pm3 "[6,8,1,1,1,1,1,1,1]" "[0,0,0,...]"

# GFN-xTB energy
sci-form xtb "[6,8,1,1,1,1,1,1,1]" "[0,0,0,...]"

# Stereochemistry analysis
sci-form stereo "C(F)(Cl)(Br)I"

# ECFP fingerprint + Tanimoto
sci-form ecfp "CCO" --radius 2 --n-bits 2048
sci-form tanimoto "c1ccccc1" "Cc1ccccc1"

# Solvation energy
sci-form solvation "[8,1,1]" "[0,0,0,0.757,0.586,0,-0.757,0.586,0]"

# Space group lookup
sci-form cell --a 5.43 --b 5.43 --c 5.43

# Build info
sci-form info
```

Prebuilt binaries at [GitHub Releases](https://github.com/jigonzalez930209/sci-form/releases):

| Platform | File |
|----------|------|
| Linux x86_64 | `sci-form-linux-x86_64` |
| Linux aarch64 | `sci-form-linux-aarch64` |
| macOS x86_64 | `sci-form-macos-x86_64` |
| macOS Apple Silicon | `sci-form-macos-aarch64` |
| Windows x86_64 | `sci-form-windows-x86_64.exe` |

→ [Full CLI reference](https://jigonzalez930209.github.io/sci-form/api/cli) · [Guide](https://jigonzalez930209.github.io/sci-form/guide/cli)

---

## Module Overview

### Top-level API

| Function | Description |
|----------|-------------|
| `sci_form::embed` | ETKDGv2 3D conformer generation from SMILES |
| `sci_form::embed_batch` | Parallel batch conformer generation (rayon) |
| `sci_form::parse` | SMILES → molecular graph |
| `sci_form::compute_charges` | Gasteiger-Marsili partial charges |
| `sci_form::compute_sasa` | Solvent-accessible surface area |
| `sci_form::compute_population` | Mulliken & Löwdin population analysis |
| `sci_form::compute_dipole` | Molecular dipole moment in Debye |
| `sci_form::compute_esp` | Electrostatic potential grid |
| `sci_form::compute_dos` | Density of states + PDOS |
| `sci_form::compute_pm3` | PM3 semi-empirical SCF |
| `sci_form::compute_xtb` | GFN0-xTB tight-binding |
| `sci_form::solve_hf3c` | HF-3c with D3/gCP/SRB corrections |
| `sci_form::compute_ani` | ANI-2x neural network potential |
| `sci_form::compute_ir_spectrum` | IR spectrum (numerical Hessian) |
| `sci_form::analyze_stereo` | CIP stereochemistry (R/S, E/Z, helical, atropisomeric) |
| `sci_form::compute_ecfp` | ECFP fingerprints (Morgan algorithm) |
| `sci_form::compute_tanimoto` | Tanimoto similarity |
| `sci_form::butina_cluster` | Butina RMSD clustering |
| `sci_form::compute_ml_descriptors` | 2D molecular descriptors |
| `sci_form::predict_ml_properties` | LogP, solubility, druglikeness, Lipinski Ro5 |
| `sci_form::compute_nonpolar_solvation` | Non-polar SASA solvation |
| `sci_form::compute_gb_solvation` | Generalized Born (HCT) solvation |
| `sci_form::compute_sssr` | Smallest Set of Smallest Rings |
| `sci_form::create_unit_cell` | Periodic unit cell from lattice parameters |
| `sci_form::assemble_framework` | MOF-type framework assembly |

### Sub-modules

| Path | Key API |
|------|---------|
| `sci_form::eht` | `solve_eht()`, `compute_band_structure()`, `optimize_geometry_eht()` |
| `sci_form::xtb::gfn1` | `solve_gfn1()` — D3-BJ dispersion, shell charges |
| `sci_form::xtb::gfn2` | `solve_gfn2()` — D4, multipole, halogen bonds |
| `sci_form::pm3` | `Pm3Result` — HOF, HOMO/LUMO, thermochemistry |
| `sci_form::hf` | `cisd::compute_cisd()` — excited states |
| `sci_form::ani` | `ani_tm::compute_aevs_tm()` — 24-element AEVs |
| `sci_form::population::npa` | `compute_npa()`, `compute_nbo()` |
| `sci_form::dos::multi_method` | `compute_dos_multimethod()`, `DosMethod` |
| `sci_form::ml::whim` | `compute_whim()`, `WhimWeighting` |
| `sci_form::ml::rdf_descriptors` | `compute_rdf()` |
| `sci_form::ml::getaway` | `compute_getaway()` |
| `sci_form::ml::advanced_models` | `train_random_forest()`, `train_gradient_boosting()` |
| `sci_form::ir::peak_assignment` | `assign_peaks()`, functional group recognition |
| `sci_form::nmr` | `compute_ensemble_j_couplings()` |
| `sci_form::smirks` | `parse_smirks()`, `apply_smirks()` |
| `sci_form::periodic` | `build_periodic_molecule()`, `detect_hapticity()` |
| `sci_form::materials::space_groups` | `space_group_by_number()`, `space_group_by_symbol()` |
| `sci_form::materials::geometry_opt` | `optimize_framework()`, `FrameworkOptConfig` |
| `sci_form::esp` | `compute_esp_grid_parallel()`, `export_cube()` |
| `sci_form::alignment` | `align_coordinates()`, `align_quaternion()` |
| `sci_form::forcefield` | `build_uff_force_field()`, `Mmff94Builder::build()` |
| `sci_form::transport` | `pack_batch_arrow()`, `split_worker_tasks()` |

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

## Experimental Engine

The repository includes an isolated experimental quantum-chemistry stack under `sci_form::experimental_2`.

### Phase Status

| Phase | Status | Scope |
|-------|--------|-------|
| Phase 1 | ✅ Complete | GPU infrastructure scaffolding, aligned types, CPU fallback, WGSL interface stubs |
| Phase 2 | ✅ Complete | Gaussian basis, overlap/kinetic/nuclear/core matrices, ERIs, validation helpers |
| Phase 3 | ✅ Complete | RHF SCF loop, Löwdin orthogonalization, DIIS, Mulliken analysis, gradients, parallel ERI |
| Phase 4 | 🔶 Prototype | Experimental sTDA UV-Vis, GIAO-style NMR shielding, Hessian/IR workflows |
| Phase 5 | 🔶 Prototype | Orbital grid evaluation, marching cubes, isosurface generation, GPU-ready rendering path |

- Stable production APIs remain unchanged; experimental work is isolated in `experimental_2`
- Real acceleration today is CPU parallelism via rayon in the ERI build path
- GPU execution is not yet enabled; `phase1_gpu_infrastructure` exposes a CPU fallback plus WGSL-ready interfaces

### Validation Snapshot

```bash
cargo check --tests
cargo test --test regression test_experimental_comparison -- --nocapture   # 54 passed
cargo test --test regression test_extended_molecules -- --nocapture         # 14 passed
```

More details in [TESTING.md](TESTING.md) and [ROADMAP.md](ROADMAP.md).

---

## Benchmark Results

### Conformer Generation — Diverse Molecules (131 molecules)

| Metric | Value |
|--------|-------|
| Parse success | 100% |
| Embed success | 97.7% |
| Throughput | 60+ mol/s |

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

## Building from Source

```bash
# Library + CLI
cargo build --release

# Python bindings
cd crates/python && pip install maturin && maturin develop --release

# WASM bindings
cd crates/wasm && ./build.sh --web-only --web-features "parallel experimental-gpu"

# With GPU acceleration
cargo build --release --features experimental-gpu
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
cargo test --release --test benchmarks -- --nocapture

# Lint & format
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
```

---

## Releasing a New Version

Use the provided bump script — it updates all version strings, commits, tags, and pushes:

```bash
# Auto-increment patch
./scripts/bump_version.sh

# Set a specific version
./scripts/bump_version.sh 0.11.0
```

This updates `Cargo.toml`, `crates/cli/Cargo.toml`, `crates/python/Cargo.toml`, `crates/wasm/Cargo.toml`, `crates/python/pyproject.toml`, and the npm package manifests, then creates a `vX.Y.Z` git tag. The [release workflow](.github/workflows/release.yml) automatically publishes to **crates.io**, **PyPI**, and **npm**.

**Required repository secrets:**

| Secret | Used for |
|--------|----------|
| `CARGO_REGISTRY_TOKEN` | crates.io |
| `PYPI_API_TOKEN` | PyPI (`sciforma`) |
| `NPM_TOKEN` | npm (`sci-form-wasm`) — Granular Automation token |

---

## License

MIT
