# sci-form

**High-performance 3D molecular conformer generation** using ETKDG distance geometry, written in Rust.

Generates chemically valid 3D coordinates from SMILES strings, matching RDKit's ETKDGv2 quality while offering native bindings for **Rust**, **Python**, **TypeScript/JavaScript (WASM)**, and a **cross-platform CLI**.

[![crates.io](https://img.shields.io/crates/v/sci-form)](https://crates.io/crates/sci-form)
[![PyPI](https://img.shields.io/pypi/v/sciforma)](https://pypi.org/project/sciforma)
[![npm](https://img.shields.io/npm/v/sci-form-wasm)](https://www.npmjs.com/package/sci-form-wasm)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

## Features

- **ETKDG Distance Geometry** — Cambridge Structural Database torsion preferences (837 SMARTS patterns)
- **High Accuracy** — 0.00% heavy-atom RMSD > 0.5 Å vs RDKit on GDB-20 (2000 molecules, ensemble comparison)
- **Fast** — 60+ molecules/second in Rust, parallel batch processing via rayon
- **Multi-platform** — Rust lib, Python (PyO3), TypeScript/JS (WASM), CLI (Linux/macOS/Windows)
- **Zero runtime dependencies** — pure Rust, no C++ toolchain needed
- **SMILES + SMARTS** — full SMILES parser and SMARTS pattern matching engine

---

## Installation

### Rust

```toml
[dependencies]
sci-form = "0.1"
```

```rust
let result = sci_form::embed("CCO", 42);
println!("Atoms: {}, Coords: {:?}", result.num_atoms, result.coords);
```

→ [Rust API docs](https://docs.rs/sci-form) · [Full guide](https://jigonzalez930209.github.io/sci-form/guide/rust)

---

### Python

The Python package is published as **`sciforma`** on PyPI.
The import name is `sci_form`.

```bash
pip install sciforma
```

```python
import sci_form

result = sci_form.embed("CCO")            # seed=42 is default
print(f"Atoms: {result.num_atoms}, Time: {result.time_ms:.1f}ms")
positions = result.get_positions()        # list of (x, y, z) tuples

# Batch
results = sci_form.embed_batch(["CCO", "c1ccccc1", "CC(=O)O"])
for r in results:
    if r.is_ok():
        print(f"{r.smiles}: {r.num_atoms} atoms")
```

→ [Full Python guide](https://jigonzalez930209.github.io/sci-form/guide/python)

---

### TypeScript / JavaScript (Node.js + Browser)

The npm package is **`sci-form-wasm`**.

```bash
npm install sci-form-wasm
```

**Node.js (CommonJS)**
```javascript
const sci = require('sci-form-wasm');
const result = JSON.parse(sci.embed('CCO', 42));
console.log(`Atoms: ${result.num_atoms}`);
```

**ES Module / TypeScript**
```typescript
import init, { embed } from 'sci-form-wasm';
await init();
const result = JSON.parse(embed('CC(=O)O', 42));
console.log(`Atoms: ${result.num_atoms}`);
```

→ [Full TypeScript guide](https://jigonzalez930209.github.io/sci-form/guide/typescript)

---

### CLI

**Download prebuilt binary** from [GitHub Releases](https://github.com/jigonzalez930209/sci-form/releases):

| Platform | File |
|----------|------|
| Linux x86_64 | `sci-form-linux-x86_64` |
| Linux aarch64 | `sci-form-linux-aarch64` |
| macOS x86_64 | `sci-form-macos-x86_64` |
| macOS Apple Silicon | `sci-form-macos-aarch64` |
| Windows x86_64 | `sci-form-windows-x86_64.exe` |

**Or install via cargo:**
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

# Show version / features
sci-form info
```

→ [Full CLI guide](https://jigonzalez930209.github.io/sci-form/guide/cli)

---

## Benchmark Results

### Diverse Molecules (131 molecules, all chemical functional groups)

| Metric | Value |
|--------|-------|
| Parse success | 100% |
| Embed success | 97.7% |
| Geometry quality | 97.7% |
| Throughput | 60 mol/s |

### RDKit Comparison (heavy-atom pairwise-distance RMSD)

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

## Algorithm

sci-form implements ETKDGv2 (Experimental Torsion Knowledge Distance Geometry):

1. **SMILES Parsing** → Molecular graph with atoms, bonds, hybridization
2. **Bounds Matrix** → 1-2, 1-3, 1-4, and VdW distance bounds from topology
3. **Triangle Smoothing** → Floyd-Warshall triangle inequality enforcement
4. **Distance Picking** → Random distances from smoothed bounds (MinstdRand)
5. **Metric Matrix Embedding** → Eigendecomposition → 4D coordinates
6. **Bounds Force Field** → BFGS minimization in 4D to satisfy distance constraints
7. **Projection to 3D** → Drop lowest-variance dimension
8. **ETKDG 3D Refinement** — Force field with CSD torsion preferences (837 patterns)
9. **Validation** — Tetrahedral centers, planarity, double-bond geometry

See the [algorithm documentation](https://jigonzalez930209.github.io/sci-form/algorithm/overview) for mathematical derivations and step-by-step diagrams.

---

## Building from Source

```bash
# Library + CLI
cargo build --release

# Python bindings
cd crates/python && pip install maturin && maturin develop --release

# WASM bindings
cd crates/wasm && wasm-pack build --target bundler --release
```

---

## Testing

```bash
# Unit tests
cargo test --lib

# Integration — diverse molecules
cargo test --release --test test_diverse_molecules -- --nocapture

# Integration — geometry quality (requires gdb20_reference.json, see scripts/)
cargo test --release --test test_geometry_quality -- --nocapture

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
