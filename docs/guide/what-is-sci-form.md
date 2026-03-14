# What is sci-form?

**sci-form** is a high-performance 3D molecular conformer generator that converts SMILES strings into chemically valid 3D coordinates.

It implements the **ETKDGv2** (Experimental Torsion Knowledge Distance Geometry) algorithm — the same algorithm used by RDKit — entirely in Rust, with no C++ dependencies.

## Why sci-form?

| Feature | sci-form | RDKit | OpenBabel |
|---------|----------|-------|-----------|
| Language | Pure Rust | C++ / Python | C++ |
| Install | `cargo add` / `pip install` / `npm install` | Complex build | Complex build |
| WASM support | ✅ Native | ❌ | ❌ |
| Accuracy vs RDKit | 0.064 Å avg | — | ~0.5 Å avg |
| Speed | 60 mol/s | ~50 mol/s | ~30 mol/s |
| Binary size | ~2 MB | ~100 MB | ~50 MB |

## How It Works

sci-form takes a SMILES string and produces 3D atomic coordinates through a 9-step pipeline:

1. **Parse** the SMILES into a molecular graph
2. **Build** a distance bounds matrix from topology
3. **Smooth** bounds via Floyd-Warshall triangle inequality
4. **Pick** random distances within the smoothed bounds
5. **Embed** via eigendecomposition of the metric matrix
6. **Minimize** in 4D using a bounds violation force field
7. **Project** to 3D by dropping the lowest-variance dimension
8. **Refine** using the ETKDG 3D force field with CSD torsion preferences
9. **Validate** tetrahedral centers, planarity, and double-bond geometry

Each step is described in detail in the [Algorithm section](/algorithm/overview).

## Chemical Coverage

sci-form handles a wide range of chemical structures:

- **Organic**: All functional groups (alcohols, amines, carbonyls, acids, esters, amides, ethers)
- **Stereochemistry**: Tetrahedral (R/S) and cis/trans (E/Z) with chiral volume constraints
- **Ring systems**: Fused, bridged, spiro, and macrocyclic structures
- **Heterocycles**: Nitrogen, oxygen, and sulfur heterocycles
- **Halogens**: F, Cl, Br, I
- **Extended elements**: Silicon, boron, selenium, phosphorus
- **Aromatics**: Benzene, naphthalene, pyridine, and complex polyaromatics

## Platforms

sci-form provides native bindings for 4 platforms:

- **[Rust](/guide/rust)** — Core library via crates.io
- **[Python](/guide/python)** — PyO3 bindings via PyPI
- **[TypeScript / JS](/guide/typescript)** — WASM bindings via npm
- **[CLI](/guide/cli)** — Prebuilt binaries for Linux, macOS, Windows
