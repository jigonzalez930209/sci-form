# What is sci-form?

**sci-form** is a high-performance computational chemistry library written in pure Rust that provides:

1. **3D Molecular Conformer Generation** — converts SMILES strings to chemically valid 3D coordinates using the ETKDGv2 algorithm
2. **Quantum Chemistry Properties** — Extended Hückel Theory (EHT), orbital analysis, population analysis, dipole moments
3. **Electrostatics** — ESP grids, Gasteiger partial charges, SASA
4. **Spectroscopy** — Density of States (DOS/PDOS) from EHT orbital energies
5. **Force Fields** — UFF and MMFF94 energy evaluation
6. **Molecular Alignment** — Kabsch + quaternion optimal superposition
7. **Materials** — Periodic unit cells, SBU topology, framework assembly

Everything is available from four entry points: **Rust**, **Python**, **TypeScript/WASM**, and **CLI** — with no C++ dependencies and native performance.

## Why sci-form?

| Feature | sci-form | RDKit | OpenBabel |
|---------|----------|-------|-----------|
| Language | Pure Rust | C++ / Python | C++ |
| Install | `cargo add` / `pip install` / `npm install` | Complex build | Complex build |
| WASM support | ✅ Native | ❌ | ❌ |
| Conformer accuracy vs RDKit | 0.064 Å avg | — | ~0.5 Å avg |
| Conformer speed | 60 mol/s | ~50 mol/s | ~30 mol/s |
| EHT / DOS / ESP | ✅ Built-in | Partial (custom) | ❌ |
| MMFF94 | ✅ Pure Rust | ✅ C++ only | ✅ C++ only |
| Materials (MOF) | ✅ Built-in | ❌ | Partial |
| Binary size | ~2 MB | ~100 MB | ~50 MB |

## The Conformer Pipeline

sci-form generates 3D conformers through a 9-step ETKDGv2 pipeline:

1. **Parse** the SMILES into a molecular graph (atoms, bonds, hybridization, stereo)
2. **Build** a distance bounds matrix from topology (1-2, 1-3, 1-4, VdW)
3. **Smooth** bounds via Floyd-Warshall triangle inequality
4. **Pick** random distances within the smoothed bounds (MinstdRand)
5. **Embed** via eigendecomposition of the metric matrix into 4D coordinates
6. **Minimize** in 4D using a bounds violation force field (BFGS)
7. **Project** to 3D by dropping the lowest-variance dimension
8. **Refine** using the ETKDG 3D force field with 846 CSD torsion patterns
9. **Validate** tetrahedral centers, planarity, and double-bond geometry

Each step is detailed in the [Algorithm section](/algorithm/overview).

## Quantum Chemistry (EHT)

The Extended Hückel Theory module provides semi-empirical molecular orbital calculations based on the Wolfsberg-Helmholtz approximation:

- **Overlap integrals** — STO-3G contracted Gaussians for H, C, N, O, F, S, Cl, and transition metals
- **Hamiltonian** — diagonal Hückel parameters, off-diagonal Wolfsberg-Helmholtz $H_{ij} = \frac{K}{2} S_{ij}(H_{ii} + H_{jj})$
- **Löwdin orthogonalization** — symmetric $S^{-1/2}$ transform for orthogonal MOs
- **Results** — orbital energies (eV), coefficients, electron count, HOMO/LUMO levels

From EHT results, sci-form derives:
- **Mulliken + Löwdin population analysis** — charge per atom and per orbital
- **Dipole moment** — bond + lone-pair contributions in Debye
- **Volumetric orbital grids** — STO-3G probability density on 3D grid, chunked for large molecules
- **Isosurface meshes** — marching cubes from orbital grids

## Electrostatic Potential

The ESP module computes the Coulomb electrostatic potential on a 3D grid:

$$V(\mathbf{r}) = \sum_i \frac{q_i}{|\mathbf{r} - \mathbf{r}_i|}$$

where $q_i$ are Mulliken charges from EHT. The grid can be:
- Exported as a Gaussian Cube file (`.cube`)
- Color-mapped to RGB (red = negative, white = zero, blue = positive)
- Computed in parallel via rayon (`parallel` feature)

## Chemical Coverage (SMILES)

sci-form handles the full periodic table in SMILES bracket notation:
- **Organic**: All functional groups (alcohols, amines, carbonyls, acids, esters, amides, ethers, sulfides, phosphates, …)
- **Stereochemistry**: Tetrahedral (R/S) and cis/trans (E/Z) with chiral volume constraints
- **Ring systems**: Fused, bridged, spiro, and macrocyclic structures (ETKDGv3 macrocycle support)
- **Heterocycles**: N, O, S-heterocycles
- **Halogens**: F, Cl, Br, I
- **Extended elements**: He→Bi (60+ bracket elements including Si, Ge, As, Se, Sn, Sb, Te)
- **Transition metals**: Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Mo, Pd, Ag, Pt, Au

## Platforms

sci-form provides native bindings for 4 platforms:

- **[Rust](/guide/rust)** — Core library via crates.io, all modules available
- **[Python](/guide/python)** — PyO3 bindings via PyPI (`sciforma`), all top-level functions
- **[TypeScript / JS](/guide/typescript)** — WASM bindings, including typed-array APIs for high-throughput use
- **[CLI](/guide/cli)** — Prebuilt binaries for Linux, macOS, Windows
