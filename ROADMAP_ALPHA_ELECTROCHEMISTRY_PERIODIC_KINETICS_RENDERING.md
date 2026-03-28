# Alpha Roadmap: Electrochemical Interfaces, Periodic Linear-Scaling Solvers, HTST Kinetics, and Rendering Data Bridges

> Status: alpha-only planning document
>
> Scope: all work described here remains isolated from stable production APIs until numerical validation, benchmark coverage, and cross-language bindings are complete.

## Mission

Build the next alpha layer of sci-form around four connected capabilities:

1. Explicit electrical double layer models for solid-liquid interfaces.
2. Periodic boundary conditions and Brillouin-zone sampling for linear-scaling electronic structure.
3. Harmonic transition state theory and microkinetic rate evaluation on top of reaction-path tooling.
4. Structured data bridges for browser and GPU-ready scientific rendering.

This roadmap is intentionally grounded in what already exists in the library today:

- Dynamic geometry-dependent charges in `charges_eeq`.
- ALPB solvation in `solvation_alpb`.
- Constant-potential infrastructure in `beta::cpm` behind `experimental-cpm`.
- Periodic cell and periodic graph logic in `materials` and `periodic`.
- EHT band structure support with explicit k-point paths in `eht::band_structure`.
- Linear-scaling beta solvers in `beta::kpm` and `beta::rand_nla`.
- Reaction-path infrastructure in `alpha::gsm` and vibrational reduction in `beta::mbh`.
- GPU kernels and context management in `gpu::*`.
- Columnar transport and batch splitting in `transport::*`.
- Python, WASM, and CLI experimental binding patterns already used for current feature-gated modules.

The goal is not to duplicate those subsystems. The goal is to connect, extend, and harden them into a coherent alpha research platform for electrochemistry, porous materials, kinetics, and high-performance visualization.

## Alpha Execution Rules

- Every new module in this roadmap lives under `src/alpha/` first, even if it wraps existing beta or core infrastructure.
- Stable public APIs are not changed during alpha stage unless there is a direct bug fix.
- New feature flags should be additive and narrow.
- CPU code paths must exist first for every capability.
- GPU paths accelerate the same kernels, not alternate physics.
- Rust implementations remain the source of truth; Python, WASM, and CLI follow after the Rust surface is stable.
- Result types must stay serialization-friendly from the start.
- Large outputs must support both JSON and columnar/typed-array transport.

## Proposed Alpha Feature Flags and Namespaces

| Alpha Track | Namespace | Feature Flag | Depends On |
|-------------|-----------|--------------|------------|
| Electrical double layer | `sci_form::alpha::edl` | `alpha-edl` | `experimental-cpm`, `experimental-eeq`, `experimental-alpb` |
| Periodic linear-scaling solvers | `sci_form::alpha::periodic_linear` | `alpha-periodic-linear` | `experimental-kpm`, `experimental-randnla` |
| HTST and microkinetics | `sci_form::alpha::kinetics` | `alpha-kinetics` | `experimental-gsm`, `experimental-mbh` |
| Rendering data bridges | `sci_form::alpha::render_bridge` | `alpha-render-bridge` | `experimental-kpm`, `experimental-gsm`, `experimental-gpu-rendering` |

## Cross-Cutting Reuse From the Current Library

### Physics and chemistry kernels already available

- `src/charges_eeq/` for dynamic charge equilibration and charge-dependent energy pieces.
- `src/solvation_alpb/` for implicit solvent response and Born-radius logic.
- `src/beta/cpm/` for constant-potential charge equilibration and grand-potential style workflows.
- `src/periodic.rs` and `src/materials/cell.rs` for unit cells, image handling, and periodic graph construction.
- `src/eht/band_structure.rs` for Bloch-style Hamiltonian construction at k-points.
- `src/beta/kpm/` for Chebyshev DOS and density-matrix infrastructure.
- `src/beta/rand_nla/` for randomized low-rank eigensolvers.
- `src/alpha/gsm/` and `src/beta/mbh/` for reaction paths, saddle localization, and reduced vibrational treatments.

### Performance and transport infrastructure already available

- `rayon`-based CPU parallelism patterns already used in DOS, ESP, surfaces, and batch work.
- `src/gpu/context.rs` and `src/gpu/shader_registry.rs` for optional GPU execution and fallback behavior.
- `src/gpu/eeq_gpu.rs`, `src/gpu/cpm_gpu.rs`, and `src/gpu/alpb_born_gpu.rs` as direct starting points for electrochemical kernels.
- `src/transport/arrow.rs`, `src/transport/chunked.rs`, and `src/transport/worker.rs` for zero-copy-style and worker-ready transport.
- Existing experimental bindings in `crates/python/src/experimental.rs`, `crates/wasm/src/experimental.rs`, and `crates/cli/src/experimental_cmds.rs`.

## Program-Level Parallelization Strategy

The four tracks should not be executed as one serial chain. They should be split into parallel implementation lanes with a strict dependency graph.

### Lane A: Core physics kernels

- EDL field terms, Stern-layer parameterization, ionic profiles, capacitance extraction.
- Periodic KPM and periodic RandNLA operator kernels.
- HTST partition-function assembly and microkinetic rate matrices.

### Lane B: GPU acceleration

- GPU kernels only after CPU reference kernels pass deterministic regression tests.
- Shared memory layout for pairwise electrostatics, Chebyshev recurrences, and spectral accumulators.
- Columnar buffer layouts aligned with WASM typed arrays and WebGPU upload buffers.

### Lane C: Data model and serialization

- Shared result schemas for EDL profiles, k-mesh DOS, HTST rate tables, and trajectory animation payloads.
- Arrow-style column batches and typed-array exports.
- Browser-first payload shaping to avoid late-stage data marshaling redesign.

### Lane D: Validation and benchmarking

- Gold-standard regression molecules and periodic fixtures.
- CPU vs GPU consistency checks.
- Accuracy-vs-speed sweeps over system size, polynomial order, sketch rank, mesh density, and temperature.

### Lane E: Bindings and examples

- Python notebooks and simple scripts for numerical verification.
- WASM examples for interactive DOS/path/profile plots.
- CLI commands for reproducible batch validation.

## Sprint Prioritization

The roadmap is intentionally large. These sprint bands define the recommended execution order and make the alpha program parallelizable without losing dependency discipline.

### Sprint Alpha 1: Contracts, adapters, and CPU references

Primary goal: establish stable alpha-shaped interfaces before optimizing kernels.

- Finalize all result schemas for EDL, periodic linear solvers, HTST, and rendering bridges.
- Add adapter layers around existing CPM, ALPB, EEQ, band structure, KPM, RandNLA, GSM, and MBH code.
- Implement CPU reference paths for the first compact-layer EDL solver, the first periodic mesh utilities, the first HTST rate calculator, and the first rendering payload builders.
- Freeze transport format decisions for JSON, typed arrays, and Arrow-style buffers.
- Add baseline regression fixtures for slabs, periodic cells, GSM paths, and DOS outputs.

### Sprint Alpha 2: Core physics and scaling kernels

Primary goal: make each track scientifically useful on CPU.

- Complete Helmholtz and Gouy-Chapman-Stern EDL profiles with CPM coupling.
- Implement Monkhorst-Pack k-mesh generation and periodic operator assembly.
- Extend KPM and RandNLA to periodic operators with k-point loops.
- Build the GSM-to-MBH-to-HTST calculation chain and the small-network microkinetic solver.
- Add chart-ready and animation-ready payloads for KPM DOS, GSM trajectories, and EDL profiles.

### Sprint Alpha 3: Parallel CPU scaling and validation hardening

Primary goal: remove bottlenecks and prove the alpha outputs are numerically trustworthy.

- Parallelize bias sweeps, k-point loops, temperature sweeps, and data chunking with `rayon`.
- Add error estimation and convergence diagnostics for all four tracks.
- Validate analytic limits, exact small-system references, and conservation laws.
- Benchmark CPU scaling and establish workload-size thresholds for GPU promotion.
- Propagate the stable alpha surfaces to Python, WASM, and CLI examples.

### Sprint Alpha 4: Selective GPU acceleration

Primary goal: accelerate only the kernels that clearly benefit from it.

- GPU-accelerate EDL profile scans with large interface grids.
- GPU-batch periodic KPM and select RandNLA stages for many k-points or large basis sets.
- Keep HTST and microkinetic scalar work CPU-first unless large network sweeps justify GPU batching.
- Add GPU-backed rendering buffers only after the CPU data model is stable.

### Sprint Alpha 5: Binding rollout and beta-readiness gate

Primary goal: turn the alpha stack into a usable cross-language platform.

- Complete Python, WASM, and CLI adapters for the finalized alpha surfaces.
- Add end-to-end examples for EDL scans, periodic DOS, HTST rates, and rendering payload export.
- Lock down the final benchmark suite and promotion criteria.
- Decide which alpha modules have enough evidence to move toward beta staging.

## Track A: Explicit Electrical Double Layer Models for Solid-Liquid Interfaces

### Goal

Extend the existing constant-potential, ALPB, and EEQ building blocks into explicit interfacial electrochemistry models that can predict charge-potential response, potential-drop profiles, compact/diffuse-layer decomposition, and differential capacitance.

### Why this track matters

Without an explicit double-layer model, the current constant-potential and implicit-solvent pieces can produce electrode charge states, but they cannot yet translate those charges into interfacial structure-aware dielectric response suitable for comparison against electrochemical impedance workflows.

### Existing code to reuse directly

- `src/beta/cpm/` for electrode charge equilibration.
- `src/charges_eeq/` for geometry-sensitive charges and charge-energy derivatives.
- `src/solvation_alpb/` for implicit electrolyte screening baseline.
- `src/gpu/eeq_gpu.rs`, `src/gpu/cpm_gpu.rs`, `src/gpu/alpb_born_gpu.rs` for GPU-adjacent kernels.

### New alpha module structure

```text
src/alpha/edl/
  mod.rs
  models.rs
  stern.rs
  gouy_chapman.rs
  helmholtz.rs
  profiles.rs
  capacitance.rs
  coupling.rs
  gpu.rs
  serialize.rs
```

### A1. Interfacial data model and parameter surfaces

#### A1.1 Geometry and electrode partitioning

- Define `ElectrodeRegion`, `ElectrolyteRegion`, and `InterfacePlane` abstractions.
- Support slab-like interfaces first, porous interfaces second.
- Add atom-selection helpers that can map current molecular/materials structures into electrode and electrolyte partitions.
- Support both Cartesian interfaces and unit-cell-aware interface normals.

#### A1.2 EDL model enum and configuration layer

- Add `EdlModel::Helmholtz`, `EdlModel::GouyChapman`, `EdlModel::GouyChapmanStern`, and `EdlModel::CompactDiffuseHybrid`.
- Define configuration structs with explicit ionic strength, temperature, dielectric profile, Stern thickness, ion valence set, and reference potential.
- Separate physically fitted parameters from numerical control parameters.

#### A1.3 Result schema designed for downstream rendering

- Standardize `EdlProfileResult` with fields for distance axis, electrostatic potential, field strength, charge density, ion density, solvent dielectric profile, compact-layer drop, diffuse-layer drop, total interfacial drop, and differential capacitance.
- Include metadata fields: `backend`, `used_gpu`, `converged`, `n_iterations`, `residual`, `temperature_k`, `ionic_strength_m`, and `model_name`.

### A2. Baseline Helmholtz and Stern-layer models

#### A2.1 Helmholtz compact-layer solver

- Implement a first alpha compact-layer model that maps electrode surface charge to linear potential drop across a user-defined compact layer.
- Couple to CPM charge output so the compact-layer drop is computed from actual constant-potential charge states, not synthetic charge input.
- Support both scalar dielectric constants and distance-dependent compact-layer dielectric profiles.

#### A2.2 Stern-layer correction on top of ALPB response

- Build a split response model: compact Stern layer plus diffuse solvent/electrolyte response.
- Reuse ALPB-derived dielectric estimates as the far-field solvent baseline.
- Expose the decomposition explicitly in results so later EIS fitting can isolate compact and diffuse contributions.

#### A2.3 CPU optimization plan

- Parallelize profile assembly over interface grid points with `rayon`.
- Precompute reusable distance bins and atom-to-interface signed distances.
- Cache compact-layer Green-function pieces between repeated potential sweeps.

#### A2.4 GPU optimization plan

- Offload per-point potential and field assembly to a compute kernel when profile grids exceed threshold size.
- Reuse the existing GPU context and fallback semantics from the current GPU subsystem.
- Keep reduction of global observables on CPU first for determinism, then move reductions to GPU only after bitwise or tolerance-based validation.

### A3. Gouy-Chapman and Gouy-Chapman-Stern diffuse-layer models

#### A3.1 Diffuse ion profile solver

- Implement a one-dimensional Poisson-Boltzmann solver for symmetric electrolytes first.
- Extend to asymmetric electrolytes and mixed valence later in the same track.
- Allow either analytic Gouy-Chapman evaluation or a discretized finite-difference solve depending on the parameter regime.

#### A3.2 CPM-to-diffuse-layer coupling

- Use CPM charge output as the boundary condition at the interface plane.
- Allow potential-controlled and charge-controlled modes for cross-checking.
- Add sweep utilities to compute charge-potential curves and capacitance minima/maxima across electrode bias.

#### A3.3 Dynamic EEQ coupling

- Recompute EEQ charges for adsorbates and near-interface species during potential sweeps.
- Add a fixed-point coupling loop: CPM surface charge -> field profile -> local polarization response -> EEQ update -> recompute CPM.
- Provide convergence acceleration with damping or Anderson mixing.

#### A3.4 CPU optimization plan

- Parallelize over independent bias points in capacitance scans.
- Parallelize over electrolyte compositions in parameter sweeps.
- Use vector-friendly dense recurrence loops for one-dimensional Poisson-Boltzmann solves.

#### A3.5 GPU optimization plan

- GPU-accelerate profile solves when evaluating many potentials or many pore/interface slices in parallel.
- Use batched kernels where each workgroup handles one bias point or one interface slice.
- Keep a CPU reference solve for every kernel family and enforce regression parity.

### A4. Differential capacitance, impedance-ready observables, and validation

#### A4.1 Differential capacitance extraction

- Implement `C_diff = dQ/dV` with robust finite-difference and spline-smoothed modes.
- Add compact-layer, diffuse-layer, and total capacitance decomposition.
- Support temperature-dependent capacitance surfaces.

#### A4.2 EIS-oriented observable preparation

- Export quasi-static observables that can later seed frequency-domain impedance models.
- Add profiles of potential, charge accumulation, screening length, and effective dielectric response vs potential.
- Package outputs so downstream analysis can compare atomistic response curves against equivalent-circuit fits without re-running expensive kernels.

#### A4.3 Validation matrix

- Validate Gouy-Chapman limits against known analytic solutions.
- Validate Helmholtz capacitance against closed-form slab cases.
- Check monotonicity and symmetry where expected for symmetric electrolytes.
- Benchmark CPU and GPU consistency across bias scans.

### A5. Alpha surface, bindings, and CLI

#### A5.1 Rust alpha API

- Add alpha-only entry points such as `compute_edl_profile`, `scan_edl_capacitance`, and `couple_cpm_edl`.
- Keep result types fully serializable and detached from internal solver structs.

#### A5.2 Python and WASM adapters

- Python: expose dense arrays directly for plotting and fitting.
- WASM: expose both JSON and typed-array profile outputs for browser plots.
- CLI: add reproducible parameter-sweep commands for batch generation.

#### A5.3 Exit criteria for this track

- Deterministic CPU solver for Helmholtz and Gouy-Chapman-Stern.
- At least one GPU-accelerated profile path for dense scan workloads.
- Capacitance and potential-drop outputs available from Rust, Python, WASM, and CLI behind alpha flags.

## Track B: Periodic Boundary Conditions and k-Point Sampling for Linear-Scaling Solvers

### Goal

Extend KPM and RandNLA beyond finite clusters so they operate on periodic materials, surfaces, and porous frameworks with explicit unit cells and Brillouin-zone integration.

### Why this track matters

The existing linear-scaling and randomized solvers only reach their practical value for large crystalline and porous systems if they can consume periodic Hamiltonians, exploit translational symmetry, and integrate over k-space rather than working on finite molecular fragments alone.

### Existing code to reuse directly

- `src/periodic.rs` and `src/materials/cell.rs`.
- `src/eht/band_structure.rs` for Bloch matrix construction and k-path infrastructure.
- `src/beta/kpm/` for Chebyshev expansions and DOS assembly.
- `src/beta/rand_nla/` for low-rank sketches, approximations, and error estimation.
- GPU infrastructure in `src/gpu/` for future batched linear algebra kernels.

### New alpha module structure

```text
src/alpha/periodic_linear/
  mod.rs
  cell_ops.rs
  kmesh.rs
  bloch.rs
  kpm_periodic.rs
  randnla_periodic.rs
  dos.rs
  density.rs
  gpu.rs
  serialize.rs
```

### B1. Periodic operator layer

#### B1.1 Unit-cell-aware operator assembly

- Add a shared periodic operator layer that constructs Hamiltonian and overlap contributions by image shell and lattice translation.
- Normalize the interface so both periodic KPM and periodic RandNLA consume the same operator abstraction.
- Support both explicit neighbor images and cutoff-pruned sparse image lists.

#### B1.2 Bloch phase application

- Factor out Bloch phase application from current band-structure logic into reusable operator kernels.
- Allow both path-based k-point evaluation and full uniform mesh integration.
- Support real and complex matrix modes where needed, but keep minimal complex allocations.

#### B1.3 Sparse representation first

- Introduce sparse or block-sparse operator storage for periodic expansions.
- Keep dense fallback mode for small-cell validation.
- Ensure the operator API can expose matrix-vector products without materializing the full dense matrix for large systems.

### B2. k-point mesh generation and integration

#### B2.1 Monkhorst-Pack mesh engine

- Implement uniform k-point mesh generation with optional gamma-centered mode.
- Add weights and symmetry labels in the result object.
- Provide deterministic ordering so DOS and density outputs remain reproducible.

#### B2.2 Symmetry-reduction scaffolding

- Start with no symmetry reduction in alpha, but design the mesh structure to store reducible and irreducible sets.
- Add placeholders for future space-group reduction using the existing materials infrastructure.

#### B2.3 Brillouin-zone integration utilities

- Provide weighted integration helpers for DOS, electron count, band energy, and k-averaged observables.
- Make the utilities generic so KPM and RandNLA reuse the same accumulation logic.

#### B2.4 CPU optimization plan

- Parallelize over k-points as the top-level coarse-grained strategy.
- Within each k-point, parallelize sparse matrix-vector work or sketch products when profitable.
- Reuse operator factorizations and image lists across all k-points in the same mesh.

#### B2.5 GPU optimization plan

- Batch many k-points together for GPU operator application.
- Upload lattice/image/operator data once, then stream k-phase vectors or Chebyshev coefficients.
- Use GPU only when mesh size and basis size exceed threshold to avoid transfer-dominated runs.

### B3. Periodic KPM

#### B3.1 Bloch-aware Chebyshev recursion

- Extend the Chebyshev recurrence to periodic Hamiltonians evaluated at each k-point.
- Keep the recurrence in operator form so large systems use repeated matrix-vector products instead of dense diagonalization.
- Support DOS and electron-count reconstruction first; density matrix second.

#### B3.2 k-averaged DOS and projected DOS

- Build `PeriodicKpmDosResult` with axes for energy grid, total DOS, optional orbital/atom projections, and mesh metadata.
- Allow accumulation during recursion to avoid storing all per-k intermediate vectors.

#### B3.3 Error control and adaptive order

- Add polynomial-order heuristics based on target spectral resolution and estimated bandwidth.
- Track convergence of DOS moments and stop early when possible for coarse scans.
- Expose quality diagnostics in results.

#### B3.4 CPU optimization plan

- Parallelize over stochastic vectors and over k-points depending on the workload shape.
- Use cache-friendly alternating buffers for Chebyshev steps.
- Fuse accumulation operations to minimize passes over large vectors.

#### B3.5 GPU optimization plan

- Map batched Chebyshev matrix-vector products to GPU compute shaders.
- Reuse workgroup-local memory for small block tiles and coefficient accumulation.
- Keep stochastic trace estimation as a batch dimension for GPU occupancy.

### B4. Periodic RandNLA

#### B4.1 k-dependent sketching

- Extend randomized sketches so they operate on Bloch operators at each k-point.
- Reuse the same random seed control across k-space for reproducibility.
- Add per-k accuracy diagnostics and a global mesh-averaged error summary.

#### B4.2 Low-rank eigenspectrum extraction for periodic systems

- Extract only the occupied and near-gap subspace required for band edges, DOS refinement, and electronic observables.
- Use fallback to exact small-cell diagonalization when the basis size is too small for randomized methods to pay off.

#### B4.3 Hybrid KPM + RandNLA workflow

- Use periodic KPM for coarse global DOS.
- Use periodic RandNLA for band-edge refinement, low-energy subspace extraction, and validation near the Fermi level.
- Standardize cross-validation utilities between both solvers.

#### B4.4 CPU optimization plan

- Parallelize sketch generation, operator application, and per-k Rayleigh-Ritz refinement.
- Reuse orthogonalization buffers and scratch spaces.
- Add basis-size and k-mesh heuristics to choose between exact, RandNLA, and KPM paths.

#### B4.5 GPU optimization plan

- Prototype GPU random projection and block orthogonalization for large batched k-mesh workloads.
- Keep QR/SVD fallbacks on CPU initially if GPU numerical stability is not yet acceptable.
- Profile transfer costs before moving refinement stages to GPU.

### B5. Periodic outputs, validation, and promotion gates

#### B5.1 Structured outputs

- Add `KMesh`, `KWeightedDos`, `PeriodicBandEdgeSummary`, and `PeriodicSpectralDiagnostics` result types.
- Ensure compatibility with rendering bridge buffers from Track D.

#### B5.2 Validation matrix

- Small-cell exact-vs-approximate comparisons for cubic and porous reference systems.
- Electron-count conservation across k-mesh integration.
- DOS convergence with respect to mesh density, polynomial order, and sketch rank.
- CPU/GPU parity on periodic DOS for representative systems.

#### B5.3 Exit criteria for this track

- Periodic KPM available for k-mesh DOS on at least one periodic EHT workflow.
- Periodic RandNLA available for band-edge and subspace refinement on the same class of systems.
- Shared k-mesh data model reused by Rust, Python, WASM, and rendering bridges.

## Track C: HTST and Microkinetic Rate Calculator

### Goal

Close the reaction-mechanism loop by converting GSM transition-state searches and MBH vibrational information into rate constants, temperature-dependent kinetic observables, and microkinetic network simulations.

### Why this track matters

The current reaction-path tooling can identify pathways and reduced vibrational structure, but it stops before the rate law. HTST is the lowest-friction next step because the major numerical ingredients already exist or are nearly available.

### Existing code to reuse directly

- `src/alpha/gsm/` for reaction paths and saddle-point discovery.
- `src/beta/mbh/` for reduced vibrational frequency workflows.
- Existing gradient infrastructures in PM3, xTB, and EHT for follow-up refinement and validation.
- Existing dynamics module for future comparison against explicit dynamical simulations.

### New alpha module structure

```text
src/alpha/kinetics/
  mod.rs
  htst.rs
  partition.rs
  prefactor.rs
  tunneling.rs
  network.rs
  solver.rs
  sensitivity.rs
  serialize.rs
```

### C1. Unified reaction-step schema

#### C1.1 GSM-to-HTST adapter layer

- Define a stable adapter that extracts reactant minimum, transition-state geometry, reaction coordinate, forward barrier, reverse barrier, and path metadata from GSM outputs.
- Support both direct GSM outputs and externally supplied TS/minimum geometries.

#### C1.2 MBH-to-partition-function adapter

- Convert MBH reduced vibrational outputs into partition-function-ready vibrational mode sets.
- Mark imaginary modes explicitly and isolate the transition-state unstable mode.
- Provide frequency sanitation rules for low-frequency cutoffs and numerical artifacts.

#### C1.3 Standard result objects

- Define `HtstRateResult`, `ReactionStepKinetics`, `MicrokineticNetwork`, and `MicrokineticTrajectory`.
- Include thermodynamic metadata, chosen approximation level, solver tolerances, and diagnostics.

### C2. Harmonic transition state theory core

#### C2.1 Partition functions

- Implement translational, rotational, vibrational, and electronic partition-function components.
- Provide both molecule-like and surface-adsorbate style modes where the translational dimensionality differs.
- Keep unit handling explicit for temperature, pressure, standard state, and concentration conventions.

#### C2.2 HTST rate expression engine

- Compute forward and reverse `k(T)` using barrier energies from GSM and partition-function ratios from MBH-derived frequencies.
- Support batch evaluation over temperature grids.
- Track contributions from entropy, zero-point energy, and thermal corrections separately.

#### C2.3 Low-friction correction terms

- Add optional Wigner tunneling first.
- Reserve Eckart or more advanced tunneling for later alpha stages.
- Add symmetry-number hooks even if only simple user-provided values are supported initially.

#### C2.4 CPU optimization plan

- Parallelize rate evaluation over temperature points and elementary steps.
- Cache partition sub-terms that are independent of temperature or shared between steps.
- Use vectorized logarithmic formulations to avoid overflow at low and high temperatures.

#### C2.5 GPU optimization plan

- GPU acceleration is not a first priority for scalar HTST formulas.
- GPU work only becomes useful for large kinetic networks or massive temperature/parameter sweeps.
- If implemented, GPU kernels should target batched partition/rate evaluations rather than the GSM or MBH solvers themselves.

### C3. Microkinetic network layer

#### C3.1 Elementary-step graph

- Build a graph structure for intermediates, elementary reactions, barriers, and stoichiometric coefficients.
- Support reversible and irreversible steps.
- Allow simple catalyst-site occupancy bookkeeping in alpha.

#### C3.2 ODE system assembly

- Convert elementary-step rates into concentration/time evolution equations.
- Support steady-state solve and explicit transient integration.
- Reuse the project’s existing numerical style for deterministic time stepping where possible.

#### C3.3 Sensitivity and degree-of-rate-control analysis

- Add first-order sensitivity to barriers and prefactors.
- Surface which elementary steps dominate overall rate or selectivity.
- Make results easy to compare against future parameter fitting loops.

#### C3.4 CPU optimization plan

- Parallelize sensitivity sweeps over parameters and temperatures.
- Use sparse stoichiometric structures for larger networks.
- Separate compile-once network assembly from many-run parameter evaluation.

#### C3.5 GPU optimization plan

- GPU path only for large ensembles of independent microkinetic solves.
- Keep a structure-of-arrays layout if kernels are added later.
- Do not move ODE kernels to GPU until profiling shows real gain over parallel CPU.

### C4. Validation, interfaces, and alpha gates

#### C4.1 Validation matrix

- Reproduce known textbook HTST rates for simple benchmark reactions.
- Check monotonic Arrhenius behavior where expected.
- Verify forward/reverse consistency against barrier differences and equilibrium constants.
- Compare MBH-based reduced-frequency rates against small exact Hessian references where possible.

#### C4.2 Rust, Python, WASM, and CLI surfaces

- Rust: expose direct step-level HTST and small-network microkinetic APIs.
- Python: prioritize ergonomic temperature scans and pandas-friendly tables.
- WASM: focus on small educational/interactive rate scans rather than large networks first.
- CLI: enable CSV/JSON export for kinetic curves and sensitivity tables.

#### C4.3 Exit criteria for this track

- One complete GSM -> MBH -> HTST pipeline running under alpha flags.
- Temperature-dependent `k(T)` available for single-step reactions.
- Small microkinetic networks supported with documented assumptions and diagnostics.

## Track D: Structured Data Bridges for Rendering and Scientific Visualization

### Goal

Turn raw scientific outputs from KPM, GSM, and the new alpha modules into durable, browser-ready, GPU-friendly payloads so the visualization pipeline is tested before full browser dynamics and WebGPU applications arrive.

### Why this track matters

The repo already has strong numerical kernels and some transport primitives. What is still missing is a visualization-facing contract: stable, typed, chunkable payloads that can move efficiently into plotting engines and real-time viewers without last-minute reshaping or bespoke adapters.

### Existing code to reuse directly

- `src/transport/arrow.rs`, `src/transport/chunked.rs`, and `src/transport/worker.rs`.
- Existing DOS and spectroscopy result structures.
- Existing orbital volumetric and mesh outputs.
- Existing experimental Python and WASM binding patterns.

### New alpha module structure

```text
src/alpha/render_bridge/
  mod.rs
  schema.rs
  dos_buffers.rs
  path_buffers.rs
  profile_buffers.rs
  mesh_buffers.rs
  chunking.rs
  wasm.rs
```

### D1. Shared rendering schema

#### D1.1 Chart-ready series contracts

- Define canonical XY and multi-series schemas for DOS, capacitance curves, potential scans, and kinetic traces.
- Support labels, units, axis metadata, and uncertainty bands.
- Keep contiguous numeric arrays for direct upload into 2D charting engines.

#### D1.2 Path and animation contracts

- Define a stable buffer layout for GSM trajectories, transition-state neighborhoods, and future MD-like sequences.
- Include per-frame coordinates, per-frame scalar energies, per-frame annotations, and optional per-atom properties.
- Support both frame-major and atom-major layouts where appropriate.

#### D1.3 Grid and field contracts

- Add buffer schemas for one-dimensional EDL profiles, two-dimensional contour maps, and three-dimensional scalar fields.
- Make metadata explicit: origin, spacing, dimensions, units, and normalization.

### D2. KPM DOS bridge

#### D2.1 Chebyshev-aware DOS payloads

- Build a direct adapter from KPM results into chart-ready columnar buffers.
- Preserve diagnostic metadata such as polynomial order, stochastic vector count, kernel type, spectral bounds, and k-mesh size when periodic KPM is used.
- Provide downsampled and full-resolution modes.

#### D2.2 Incremental streaming mode

- Add chunked DOS emission so long-running KPM computations can stream partial spectra to a UI.
- Make chunk ordering deterministic and self-describing.
- Reuse existing worker and chunking patterns in the repo.

#### D2.3 CPU and GPU delivery plan

- CPU path assembles chunked series data in parallel when many projections are requested.
- GPU path is only for generating the scientific data, not for marshaling JSON.
- Use typed arrays as the default browser transport when payload size grows beyond comfortable JSON limits.

### D3. GSM and kinetics trajectory bridge

#### D3.1 Reaction path payloads

- Convert GSM path coordinates into animation-ready buffers with synchronized energy traces.
- Include a compact representation for browser playback and a richer debug representation for offline analysis.
- Add step labels for reactant basin, climbing region, TS window, and product basin.

#### D3.2 HTST and microkinetic plot payloads

- Emit Arrhenius plots, partition-function decomposition tables, and sensitivity curves using the same chart schema as DOS.
- Enable direct overlay of multiple temperatures, catalysts, or parameter sets.

#### D3.3 Worker-ready splitting

- Split large trajectories and multi-condition kinetic plots into browser worker chunks.
- Preserve metadata enough to reassemble in order on the UI side without extra lookup tables.

### D4. EDL and periodic rendering bridges

#### D4.1 EDL profile payloads

- Emit compact one-dimensional and optional two-dimensional profile buffers from Track A.
- Support potential, field, ion-density, dielectric, and capacitance series using one shared schema.
- Allow bias-indexed multi-profile bundles.

#### D4.2 Periodic spectral payloads

- Emit k-path and k-mesh spectral outputs from Track B in a rendering-oriented layout.
- Support heatmaps, line plots, and contour plots.
- Include explicit k-coordinate and weight buffers rather than forcing a front-end reconstruction step.

#### D4.3 Browser and WebGPU readiness

- Align buffer layouts with `Float32Array` and `Float64Array` upload paths.
- Add conversion helpers for tightly packed vertex and scalar buffers.
- Avoid nested object graphs in high-volume rendering payloads.

### D5. Validation and alpha gates

#### D5.1 Transport validation

- Verify lossless round-tripping for all buffer schemas.
- Benchmark JSON vs typed-array vs columnar transport for DOS, trajectories, and profiles.
- Measure browser reconstruction cost separately from compute cost.

#### D5.2 Binding surfaces

- Python: direct NumPy-friendly arrays.
- WASM: typed arrays first, JSON fallback second.
- CLI: optional binary or JSON export depending on output mode.

#### D5.3 Exit criteria for this track

- KPM DOS, GSM trajectories, and EDL profiles all export through one coherent rendering schema family.
- At least one browser example consumes typed-array payloads without extra reshaping.
- Transport cost is benchmarked and documented for large alpha outputs.

## Recommended Delivery Order

### Phase 1: Data model and adapters

- Build alpha namespaces, config structs, and result schemas for all four tracks.
- Add adapters around current CPM, ALPB, EEQ, band structure, KPM, RandNLA, GSM, and MBH modules.
- Finish serialization and typed-array design before expanding kernels.

### Phase 2: First working CPU references

- Helmholtz and Gouy-Chapman-Stern EDL on CPU.
- Monkhorst-Pack meshes and periodic operator abstraction on CPU.
- HTST single-step calculator on CPU.
- DOS/path/profile rendering payload generators on CPU.

### Phase 3: Parallel CPU scaling

- Multi-bias EDL scans.
- Multi-k periodic KPM and RandNLA workflows.
- Multi-temperature and multi-step kinetic sweeps.
- Chunked rendering exports for large result sets.

### Phase 4: Selective GPU acceleration

- EDL profile and scan kernels with strong arithmetic intensity.
- Batched periodic KPM kernels.
- Optional batched RandNLA sketch/application steps.
- No premature GPU porting of scalar or low-arithmetic stages.

### Phase 5: Binding propagation and examples

- Python validation scripts and notebooks.
- WASM typed-array examples.
- CLI batch commands and export tooling.

## Critical Dependencies and Blocking Risks

### Scientific risks

- EDL model fidelity may depend on assumptions that are too simple for specific porous or strongly adsorbing interfaces.
- Periodic KPM and periodic RandNLA may expose conditioning or sparse-operator issues not visible in molecular mode.
- HTST quality depends on the reliability of TS localization and reduced vibrational treatments.

### Engineering risks

- GPU transfer overhead may erase the benefit of acceleration for small and medium workloads.
- Serialization can become the bottleneck if result schemas are not designed around contiguous arrays from day one.
- Cross-language binding drift can appear quickly if alpha result types change too often.

### Mitigation plan

- Keep every alpha solver paired with a deterministic CPU reference path.
- Standardize result schemas before broad binding rollout.
- Add benchmark fixtures early for slabs, periodic cells, and reaction-path cases.
- Gate promotion on measured accuracy and transport cost, not only on feature completeness.

## Validation Matrix Summary

| Track | Numerical validation | Performance validation | Binding validation |
|-------|----------------------|------------------------|-------------------|
| EDL | analytic limits, capacitance sanity, charge-potential consistency | CPU scan scaling, GPU profile thresholds | Rust + Python + WASM + CLI |
| Periodic linear | exact-vs-approx small cells, electron conservation, DOS convergence | k-point scaling, order/rank scaling, GPU batch gain | Rust + Python + WASM |
| HTST/kinetics | textbook rate checks, forward/reverse consistency, low-frequency sanitation | temperature grid sweeps, multi-step network scaling | Rust + Python + CLI, then WASM |
| Render bridge | round-trip fidelity, schema consistency, unit metadata correctness | transport size, chunking throughput, UI reconstruction cost | Python + WASM + CLI |

## Definition of Alpha Completion

This roadmap should be considered complete at alpha stage when the following are true:

- All four tracks compile behind independent alpha feature flags.
- Every track has a deterministic CPU reference implementation.
- At least the highest-payoff kernels have optional GPU acceleration with documented thresholds.
- Rust result types are stable enough to support Python, WASM, and CLI adapters.
- Large outputs can be exported without mandatory nested JSON marshaling.
- Each track has at least one benchmarked end-to-end example and one validation case.

## Final Product Vision

If executed in this order, sci-form will gain an alpha research stack that can:

- Predict electrode charge response together with explicit compact and diffuse interfacial structure.
- Run linear-scaling electronic-structure approximations on periodic materials and porous frameworks.
- Convert reaction paths into temperature-dependent rate constants and small microkinetic models.
- Feed those results directly into browser-grade scientific rendering pipelines with minimal reshaping.

That combination turns the library from a collection of strong scientific kernels into a connected electrochemical, materials, and kinetics experimentation platform that is already architected for CPU and GPU execution.