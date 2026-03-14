# sci-form

**High-performance 3D molecular conformer generation** using ETKDG distance geometry, written in Rust.

Generates chemically valid 3D coordinates from SMILES strings, matching RDKit's ETKDGv2 quality while offering native bindings for **Rust**, **Python**, **TypeScript/JavaScript (WASM)**, and a **cross-platform CLI**.

## Features

- **ETKDG Distance Geometry** — Cambridge Structural Database torsion preferences (837 SMARTS patterns)
- **High Accuracy** — 0.00% heavy-atom RMSD > 0.5 Å vs RDKit on GDB-20 (2000 molecules, ensemble comparison)
- **Fast** — 60+ molecules/second in Rust, parallel batch processing via rayon
- **Multi-platform** — Rust lib, Python (PyO3), TypeScript/JS (WASM), CLI (Linux/macOS/Windows)
- **Zero dependencies at runtime** — pure Rust, no C++ toolchain needed
- **SMILES + SMARTS** — full SMILES parser and SMARTS pattern matching engine

## Quick Start

### Rust

```toml
[dependencies]
sci-form = "0.1"
```

```rust
let result = sci_form::embed("CCO", 42);
println!("Atoms: {}, Coords: {:?}", result.num_atoms, result.coords);
```

### Python

```bash
pip install sci-form
```

```python
import sci_form

result = sci_form.embed("CCO")
print(f"Atoms: {result.num_atoms}, Time: {result.time_ms:.1f}ms")
positions = result.get_positions()  # [(x, y, z), ...]
```

### TypeScript / JavaScript

```bash
npm install sci-form
```

```typescript
import { embed } from 'sci-form';

const result = JSON.parse(embed("CCO", 42));
console.log(`Atoms: ${result.num_atoms}`);
```

### CLI

```bash
# Single molecule
sci-form embed "CCO" --format xyz

# Batch processing
sci-form batch -i molecules.smi -o output.sdf --format sdf --threads 8

# Parse only (no 3D)
sci-form parse "c1ccccc1"
```

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

## Algorithm

sci-form implements the ETKDGv2 (Experimental Torsion Knowledge Distance Geometry) algorithm:

1. **SMILES Parsing** → Molecular graph with atoms, bonds, hybridization
2. **Bounds Matrix** → 1-2, 1-3, 1-4, and VdW distance bounds from topology
3. **Triangle Smoothing** → Floyd-Warshall triangle inequality enforcement
4. **Distance Picking** → Random distances from smoothed bounds (MinstdRand)
5. **Metric Matrix Embedding** → Eigendecomposition → 4D coordinates
6. **Bounds Force Field** → BFGS minimization in 4D to satisfy distance constraints
7. **Projection to 3D** → Drop lowest-variance dimension
8. **ETKDG 3D Refinement** — Force field with CSD torsion preferences (837 patterns)
9. **Validation** — Tetrahedral centers, planarity, double-bond geometry

See [documentation](https://sci-form.dev) for detailed algorithm descriptions with mathematical derivations.

## Building from Source

```bash
# Library + CLI
cargo build --release

# Python bindings
cd crates/python && maturin develop --release

# WASM bindings
cd crates/wasm && wasm-pack build --target bundler --release
```

## Testing

```bash
# Unit tests
cargo test --lib

# Diverse molecule benchmark
cargo test --release --test test_diverse_molecules -- --nocapture

# Geometry quality (requires GDB20.50000.smi)
cargo test --release --test test_geometry_quality -- --nocapture

# Gradient correctness
cargo test --release --test test_gradient_check -- --nocapture
```

## License

MIT
