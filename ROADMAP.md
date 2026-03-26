# Roadmap: sci-form Core Conformer Engine and Orbital Platform

**Primary goal:** Build a native Rust molecular engine that reproduces RDKit-grade 3D conformer generation without heavy C++ dependencies, ships clean Rust, Python, and WASM APIs, and grows into an electronic-structure and volumetric-visualization platform for the web.

**Execution principles:**
- Single-conformer, RDKit-style embedding comes first. Accuracy should come from the algorithm, not brute-force seed farming.
- Port mathematical logic, not C++ object hierarchies.
- Keep the public API usable from CLI, Python, and WebAssembly.
- Optimize for browser deployment once correctness is stable.

**Technology stack:**
- **Language:** Rust, targeting `wasm32-unknown-unknown`
- **Linear algebra:** `nalgebra` today, with room to evaluate `ndarray-linalg` where generalized eigensolvers become a better fit
- **Graph layer:** `petgraph`
- **Reference oracle:** Python + native RDKit
- **Parallel execution:** `rayon` for batch and volumetric workloads

---

## Current Implementation Snapshot (March 2026)

- [x] Rust core crate for SMILES parsing, molecular graph construction, distance geometry, ETKDG-style refinement, and force-field minimization
- [x] Differential testing and benchmark coverage against RDKit across GDB-20 and ChemBL-style datasets
- [x] Public Rust API plus Python and WASM bindings
- [x] Release optimization with `lto = true` and `strip = true`
- [x] Reach RDKit-class throughput on the practical single-conformer path (optimized: pre-computed pair lists, stack-allocated diffs, adaptive restart limits)

**Current measured status:**
- Single-conformer internal validation: max RMSD 0.477 A, avg RMSD 0.069 A on the current 100-molecule pipeline test set
- Large-batch heavy-atom validation: 0.024 A average pairwise-distance RMSD on GDB-20, with 0.00% above 0.5 A in the multi-seed cross-check
- ChemBL 10K practical benchmark: 97.54% embed success, 97.18% perfect-geometry rate, about 2.1 mol/s
- Main remaining conformer-engine gap: throughput. The current single-conformer path is still around 106 ms/mol versus the 11 ms/mol target.

**Experimental engine snapshot:**
- `src/experimental_2/` is now wired into the crate as an isolated namespace
- Phase 1 GPU infrastructure, Phase 2 quantum engine, and Phase 3 SCF engine are implemented and compile cleanly
- Phase 4 spectroscopy and Phase 5 GPU rendering exist as prototype tracks behind the experimental namespace
- Verified regression status for the experimental stack: `54 passed` in `test_experimental_comparison` and `14 passed, 7 ignored` in `test_extended_molecules`

---

## Base Data Model Direction

The core data structures need to support both the current conformer engine and the future orbital, electrostatics, and materials stack.

```rust
use nalgebra::Vector3;

pub struct Atom {
    pub element: u8,
    pub position: Vector3<f64>,
    pub formal_charge: i8,
    pub partial_charge: Option<f64>,
    pub hybridization: Hybridization,
}

pub struct Bond {
    pub start: usize,
    pub end: usize,
    pub order: BondOrder,
}

pub struct AtomicOrbital {
    pub atom_index: usize,
    pub principal_n: u8,
    pub angular_l: u8,
    pub magnetic_m: i8,
    pub label: String,
    pub vsip: f64,
    pub zeta: f64,
}

pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub metadata: MoleculeMeta,
    pub unit_cell: Option<UnitCell>,
    pub basis: Option<Vec<AtomicOrbital>>,
}

pub struct VolumetricGrid {
    pub origin: Vector3<f64>,
    pub spacing: f64,
    pub dims: [usize; 3],
    pub values: Vec<f32>,
}
```

**Development note:** do not translate RDKit file-by-file. Extract the numerical logic and rebuild it in Rust with composition-friendly modules and data structures.

---

## Track A: 3D Conformer Engine

### Phase A1: Infrastructure and Molecular Data Model

- [x] Initialize the Rust workspace and split it into core, CLI, Python, and WASM crates
- [x] Configure `wasm-bindgen`, PyO3, packaging, and CI build coverage
- [x] Build the differential-testing pipeline against RDKit reference outputs
- [x] Implement the molecular graph, atom and bond enums, radii tables, and SMILES parsing
- [x] Expand parser coverage: 60+ bracket elements, expanded implicit H valence table, error on unknown elements

### Phase A2: Distance Geometry and Bounds Matrix

- [x] Implement 1-2, 1-3, and 1-4 bounds-matrix generation
- [x] Implement bounds smoothing with shortest-path propagation
- [x] Implement RDKit-compatible RNG plumbing and distance sampling variants
- [x] Implement metric-matrix construction, centering, and `SymmetricEigen` projection to 3D and 4D coordinates
- [x] Reduce random-coordinate fallback cost for very large molecules (adaptive thresholds, scaled box size)

### Phase A3: Chemical Refinement and Stereochemistry

- [x] Detect and enforce tetrahedral chirality during embedding
- [x] Validate chiral volume, planarity, and double-bond geometry after embedding
- [x] Reproduce the RDKit-style retry-on-failure embedding loop
- [x] Add ETKDG-lite and ETKDGv3-style torsion knowledge, including macrocycle pattern support
- [x] Close remaining fidelity gaps: added amide, sulfonamide, sulfoxide, sulfone, phosphonate, allylic torsion patterns

### Phase A4: Energy Minimization and Force Fields

- [x] Translate UFF parameter tables and atom typing for the supported organic subset
- [x] Implement bonded and non-bonded energy terms plus gradients
- [x] Implement bounds-space BFGS and L-BFGS minimizers
- [x] Implement the ETKDG 3D force field and BFGS refinement path
- [x] Move the numerically sensitive pipeline to `f64` for stability
- [x] Expand parameter coverage: 18 new UFF atom types (transition metals, main-group), MMFF94 with bond/angle/torsion/vdW terms and builder

### Phase A5: Public API, WASM, and Performance Delivery

- [x] Expose stable Rust entry points for single and batch embedding
- [x] Ship WASM bindings for full results, compact coordinates, batch embedding, and parsing
- [x] Ship Python bindings for single and batch embedding
- [x] Publish benchmark documentation and CI build steps for `wasm-pack`
- [x] Enable release optimizations in Cargo profiles
- [x] Benchmark browser-side performance against RDKit MinimalLib (browser benchmark HTML page)
- [x] Replace JSON-heavy transfer paths with typed-array oriented APIs (Float64Array/Float32Array WASM functions)
- [x] Tune final WASM size with `wasm-opt` and dedicated size profiling (wasm_opt_tune.sh script)

---

## Track B: Electronic Structure and Orbital Visualization (EHT)

This track adds semiempirical electronic-structure support and 3D orbital rendering on top of the conformer engine. The intended first method is Extended Huckel Theory (EHT), using Slater-type orbitals conceptually and Gaussian expansions for practical integral evaluation.

**Core orbital form:**

$$
\phi(r, \theta, \varphi) = N r^{n-1} e^{-\zeta r} Y_l^m(\theta, \varphi)
$$

### Phase B1: Parameterization and Orbital Basis

- [x] Build an EHT parameter table for each supported element, including valence-shell orbital definitions
- [x] Store per-orbital VSIP values to populate the diagonal Hamiltonian terms $H_{ii}$
- [x] Store per-orbital Slater exponents $\zeta$ for orbital extent and contraction
- [x] Define a Rust representation for Slater-type orbitals that is stable enough for evaluation and rendering
- [x] Add contracted Gaussian expansions such as STO-3G or STO-6G for overlap-integral evaluation

### Phase B2: Matrix Construction ($S$ and $H$)

- [x] Assemble the symmetric overlap matrix $S$ for all valence orbitals in the molecule
- [x] Enforce the diagonal condition $S_{ii} = 1$ and validate symmetry and conditioning
- [x] Fill the diagonal Hamiltonian with VSIP values for each atomic orbital
- [x] Compute off-diagonal Hamiltonian elements with the Wolfsberg-Helmholtz approximation
- [x] Make the empirical constant $K$ configurable, with 1.75 as the default starting point

**Wolfsberg-Helmholtz approximation:**

$$
H_{ij} = \frac{1}{2} K S_{ij} \left(H_{ii} + H_{jj}\right)
$$

### Phase B3: Generalized Eigenproblem and Orbital Coefficients

- [x] Solve the generalized Schrodinger equation $HC = SCE$
- [x] Implement Lowdin orthogonalization by diagonalizing $S$ and building $S^{-1/2}$
- [x] Build the orthogonalized Hamiltonian $H' = S^{-1/2} H S^{-1/2}$
- [x] Run the standard symmetric eigensolve on $H'$ to obtain orbital energies and transformed coefficients
- [x] Back-transform the coefficients with $C = S^{-1/2} C'$ and expose orbital occupations, HOMO, and LUMO metadata

### Phase B4: 3D Volumetric Mapping Engine

- [x] Build a padded bounding box around the optimized molecular geometry
- [x] Generate a configurable 3D voxel grid for orbital and density sampling
- [x] Parallelize grid evaluation with `rayon` so millions of points are tractable
- [x] Evaluate molecular-orbital amplitudes pointwise with the linear combination of atomic orbitals
- [x] Add chunked or streamed evaluation modes: x-slab chunked evaluation with VolumeSlab iterator

**Orbital evaluation in space:**

$$
\Psi_i(x,y,z) = \sum_{\mu=1}^{N} C_{\mu i} \phi_{\mu}(x,y,z)
$$

### Phase B5: Output Structures and Rendering Paths

- [x] Export raw volumetric data as `Vec<f32>` or `Float32Array` plus grid metadata for GPU volume rendering
- [x] Define a WASM-friendly transfer format for dimensions, bounds, spacing, and orbital-energy metadata
- [x] Implement a Marching Cubes path for isosurface extraction in Rust
- [x] Return mesh vertices, normals, and indices for fast WebGL or WebGPU rendering when volume textures are too heavy
- [x] Benchmark transfer cost, render cost, and memory usage for raw volumes versus isosurfaces (integration test)

---

## Track C: Post-Core Platform Expansion

### Phase C1: Electrostatics and Surface Analysis

- [x] Implement a native partial-charge solver, either Gasteiger-Marsili or an EHT-derived charge workflow
- [x] Add solvent-accessible surface area and cavity analysis with a Shrake-Rupley style algorithm
- [x] Connect electrostatic descriptors and surface outputs back into the conformer and orbital pipelines

### Phase C2: Population Analysis (Mulliken & Löwdin Charges)

- [x] Compute population matrices from density matrix $P$ and overlap matrix $S$
- [x] Implement Mulliken and Löwdin charge partitioning to extract per-atom partial charges from MO coefficients
- [x] Validate against PySCF/Multiwfn with charge-sum invariants ($\sum_i q_i = Q_{total}$)
- [x] Expose charge arrays via Rust lib, CLI, WASM, and Python bindings

### Phase C3: Dipole Moments

- [x] Calculate molecular dipole vector $\vec{\mu} = (\mu_x, \mu_y, \mu_z)$ from nuclear and electron-density contributions
- [x] Validate magnitude (in Debyes) and spatial orientation against PySCF/ORCA/Gaussian
- [x] Support conformer-dependent dipole calculations

### Phase C4: Electrostatic Potential Maps (ESP)

- [x] Generate volumetric ESP on 3D grids: $\Phi(\vec{r}) = \sum_A \frac{Z_A}{|\vec{r} - \vec{R}_A|} - \int \frac{\rho(\vec{r}')}{|\vec{r} - \vec{r}'|} d\vec{r}'$
- [x] Export/import `.cube` format for validation against PySCF cubegen or Gaussian
- [x] Color-code ESP visualization: `esp_color_map()` and `esp_grid_to_colors()` (red=negative, blue=positive)
- [x] Parallelize grid evaluation with rayon: `compute_esp_grid_parallel()` for million-point grids

### Phase C5: Density of States (DOS) and Projected DOS (PDOS)

- [x] Apply Gaussian smearing to orbital energies to generate smooth DOS curves
- [x] Implement PDOS: separate DOS by atomic contribution
- [x] Validate DOS curves: `dos_mse()` mean-squared error metric for comparison
- [x] Export as graph data for web visualization: `export_dos_json()` with energies, total_dos, pdos

### Phase C6: Molecular Alignment and RMSD (Kabsch/Quaternion)

- [x] Implement Kabsch algorithm for optimal 3D alignment of two molecular geometries
- [x] Support quaternion-based rotation fitting: `align_quaternion()` (Coutsias et al. 2004)
- [x] Compute heavy-atom, all-atom, and atom-type-specific RMSD variants
- [x] Validate: quaternion matches Kabsch RMSD to <1e-8, aligned coords to <1e-6
- [x] Use RMSD as metric for conformer diversity and benchmark comparisons

### Phase C7: Force Fields (UFF / MMFF94) and Strain Energy

- [x] Fully parameterize UFF with bonded (stretch, bend, torsion) and non-bonded (vdW, Coulomb) terms
- [x] Implement MMFF94: atom typing (28 types), bond stretch, angle bend, torsion, vdW, builder for full molecules
- [x] Compute strain energy: sum of all individual energy terms
- [x] Validate gradients: numerical vs analytical validation for all MMFF94 terms (<1e-3 max error)
- [x] Enable geometry optimization with selected force field using L-BFGS

### Phase C8: Complex Materials Assembly

- [x] Add SBU and coordination-aware assembly logic for metal centers and framework nodes
- [x] Extend the data model with periodic unit-cell vectors and periodic boundary conditions
- [x] Support browser-side generation and visualization of repeating crystalline networks

### Phase C9: Web Data Transport and Parallelism

- [x] Add zero-copy Arrow-compatible memory layouts for large coordinate, orbital, and grid datasets
- [x] Move heavy browser workloads to Web Workers with a WASM-safe threading strategy
- [x] Add chunked streaming for large volumetric and ensemble outputs to avoid JSON bottlenecks

---

## Track E: Machine Learning Potentials (ANI Family)

This track provides near-DFT accuracy at force-field speeds using physical neural networks.

- [x] **Step 1.1: Spatial Partitioning:** Implemented `CellList` for O(N) neighbor searching within a $5.2 \text{ \AA}$ cutoff radius (`src/ani/neighbor.rs`).
- [x] **Step 1.2: Atomic Environment Vectors (AEVs):** Encoded Behler-Parrinello radial and angular symmetry functions for translational/rotational invariance (`src/ani/aev.rs`, `src/ani/aev_params.rs`, `src/ani/cutoff.rs`).
- [x] **Step 1.3: Neural Network Engine:** Created a pure-Rust feed-forward inference engine supporting GELU/CELU activations with a binary weight loader (`src/ani/nn.rs`, `src/ani/weights.rs`).
- [x] **Step 1.4: Analytical Gradients:** Implemented exact backpropagation through the AEVs to yield analytical forces respecting Newton's third law (`src/ani/gradients.rs`, `src/ani/api.rs`).

---

## Track F: HF-3c Quantum Engine

This track introduces a specialized minimal-basis Hartree-Fock engine with semi-empirical corrections for rapid geometric and spectroscopic screening.

- [x] **Step 2.1: Gaussian Integrals:** Implemented Obara-Saika recurrence for analytical evaluation of overlap, kinetic, nuclear attraction, and electron repulsion integrals (`src/hf/overlap_kin.rs`, `src/hf/nuclear.rs`, `src/hf/integrals.rs`, `src/hf/basis.rs`).
- [x] **Step 2.2: SCF Cycle:** Created a Roothaan-Hall SCF solver with Löwdin orthogonalization and DIIS convergence acceleration (`src/hf/fock.rs`, `src/hf/scf.rs`).
- [x] **Step 2.3: Empirical Corrections:** Added Grimme D3-BJ dispersion, geometric counterpoise (gCP), and short-range basis (SRB) corrections (`src/hf/d3.rs`, `src/hf/gcp.rs`, `src/hf/srb.rs`).
- [x] **Step 2.4: CIS UV-Vis:** Implemented Configuration Interaction Singles (CIS) for vertical excitation energies and oscillator strengths (`src/hf/cis.rs`, `src/hf/api.rs`).

---

## Track D: Advanced Spectroscopy and Response Properties

For a detailed breakdown of the spectroscopy development plan, see [docs/spectroscopy_roadmap.md](file:///home/lestad/github/sci-form/docs/spectroscopy_roadmap.md).
Specific deep-dive for NMR: [docs/nmr_detailed_roadmap.md](file:///home/lestad/github/sci-form/docs/nmr_detailed_roadmap.md).

- [x] **UV-Vis Spectroscopy:** Vertical excitations via sTDA-xTB and spectral broadening (Phase D1)
- [x] **IR Spectroscopy:** Numerical Hessian, vibrational frequencies, and dipole intensities (Phase D2)
- [x] **NMR Spectroscopy:** Topological/ML-based 1H/13C prediction plus fast relative heteronuclear inference and expanded nucleus registry (Phase D3)

---

## Priority Order From Here

- [x] Reduce conformer-engine runtime: pre-computed active pair lists, stack-allocated diffs, adaptive restart limits
- [x] Harden first-try correctness: adaptive thresholds for large molecules, expanded torsion patterns
- [x] Finish the production-quality browser benchmarking and typed-array transport path
- [x] Extend EHT pipeline with population analysis (Phase C2–C3) for industrial chemical descriptors
- [x] Implement Kabsch/RMSD (Phase C6) as the production alignment and diversity metric
- [x] Complete force-field validation (Phase C7) before enabling conformer geometry optimization
- [x] Evaluate Arrow/WASM-Workers architecture (Phase C9) for browser-side ensemble rendering
- [x] **[NEW]** Implement sTDA-xTB module for ultrafast UV-Vis screening (Phase D1)
- [x] **[NEW]** Build numerical Hessian infrastructure for IR vibrational analysis (Phase D2)
- [x] **[NEW]** Integrate ML models for NMR chemical shift prediction (Phase D3)

---

## Track G: Large-Scale Algorithmic Validation

With the conformer engine, spectroscopic models, ML potentials (ANI), and quantum engines (HF-3c) essentially built, the system needs to be stress-tested against established standard libraries using comprehensive molecule sets.

- [ ] **Step 3.1: Assembly of Representative Molecule Bank:** Curate a diverse set of at least 100 representative molecules spanning different functional groups, stereocenters, and sizes.
- [ ] **Step 3.2: Automated Pipeline Benchmarking:** Construct automated Python-based test pipelines linking `sci-form` binaries against reference implementations (e.g., PySCF for quantum targets, RDKit for geometries/descriptors).
- [ ] **Step 3.3: Metric Evaluation:** Enforce a strict passing threshold (e.g., <1% deviation) across energy evaluations, structure gradients, and physical property predictions.
- [ ] **Step 3.4: Documentation of Scaling Limits:** Profile execution times and memory consumption on the 100-molecule test bank to pinpoint throughput limiters for optimization.
