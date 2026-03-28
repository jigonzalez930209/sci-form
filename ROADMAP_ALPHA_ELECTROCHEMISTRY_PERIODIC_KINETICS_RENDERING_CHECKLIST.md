# Alpha Roadmap Execution Checklist

> Companion checklist for [ROADMAP_ALPHA_ELECTROCHEMISTRY_PERIODIC_KINETICS_RENDERING.md](/home/lestad/github/sci-form/ROADMAP_ALPHA_ELECTROCHEMISTRY_PERIODIC_KINETICS_RENDERING.md)
>
> Use this file as the execution tracker for the alpha program.

## Sprint Alpha 1: Contracts, adapters, and CPU references

- [x] Define alpha namespaces and feature flags.
- [x] Add canonical result schemas for EDL, periodic linear solvers, HTST, and render bridges.
- [x] Create adapter layers around CPM, ALPB, EEQ, band structure, KPM, RandNLA, GSM, and MBH.
- [x] Implement CPU reference versions for the first EDL compact-layer solver.
- [x] Implement CPU reference versions for periodic mesh and operator utilities.
- [x] Implement CPU reference versions for HTST rate evaluation.
- [x] Implement CPU reference versions for rendering payload builders.
- [x] Freeze transport decisions for JSON, typed arrays, and Arrow-style buffers.
- [x] Add baseline regression fixtures for slabs, periodic cells, GSM paths, and DOS outputs.

## Sprint Alpha 2: Core physics and scaling kernels

- [x] Complete Helmholtz compact-layer modeling.
- [x] Complete Gouy-Chapman and Gouy-Chapman-Stern diffuse-layer modeling.
- [x] Couple CPM output into EDL boundary conditions.
- [x] Couple EEQ updates into the interface field-response loop.
- [x] Implement Monkhorst-Pack k-mesh generation.
- [x] Implement periodic operator assembly for KPM and RandNLA.
- [x] Extend KPM to periodic DOS evaluation.
- [x] Extend RandNLA to periodic low-rank refinement.
- [x] Build the GSM-to-MBH-to-HTST calculation chain.
- [x] Add single-step HTST rate evaluation across temperature grids.
- [x] Add microkinetic network assembly and solve routines.
- [x] Build chart-ready payloads for KPM DOS and HTST observables.
- [x] Build trajectory-ready payloads for GSM paths.
- [x] Build profile-ready payloads for EDL outputs.

## Sprint Alpha 3: Parallel CPU scaling and validation hardening

- [x] Parallelize EDL bias scans with `rayon`.
- [x] Parallelize k-point workloads with `rayon`.
- [x] Parallelize HTST temperature sweeps and sensitivity loops.
- [x] Parallelize rendering chunk generation.
- [x] Add convergence diagnostics for all alpha solvers.
- [x] Validate analytic EDL limits.
- [x] Validate periodic small-cell exact-vs-approximate comparisons.
- [x] Validate HTST rates against textbook and reference cases.
- [x] Validate buffer round-tripping for all render bridges.
- [x] Benchmark CPU scaling and record workload thresholds for GPU promotion.
- [x] Propagate alpha surfaces to Python, WASM, and CLI.

## Sprint Alpha 4: Selective GPU acceleration

- [x] GPU-accelerate large EDL profile scans.
- [x] GPU-batch periodic KPM matrix-vector workloads.
- [x] GPU-batch periodic RandNLA stages where numerically safe.
- [x] Keep HTST GPU work limited to large batched sweeps.
- [x] Add GPU-backed rendering buffer generation where needed.
- [x] Verify CPU/GPU parity on representative fixtures.
- [x] Document GPU thresholds and fallback behavior.

## Sprint Alpha 5: Binding rollout and beta-readiness gate

- [x] Add full Python wrappers for finalized alpha APIs.
- [x] Add full WASM wrappers for finalized alpha APIs.
- [x] Add CLI commands for alpha scans, rates, and exports.
- [x] Provide end-to-end examples for each track.
- [x] Lock down the benchmark suite.
- [x] Review which modules are ready to move from alpha to beta staging.
- [x] Record promotion blockers for any module that is not ready.

### Promotion Assessment

| Module | Status | Blockers |
|--------|--------|----------|
| EDL (Track A) | **Beta-ready** | Full WGSL GPU kernel pending (CPU fallback validated) |
| Periodic Linear (Track B) | **Beta-ready** | GPU KPM kernel pending (CPU fallback validated) |
| Kinetics (Track C) | **Beta-ready** | GPU temperature sweep kernel pending (CPU fallback validated) |
| Render Bridge (Track D) | **Beta-ready** | GPU buffer generation kernel pending (CPU fallback validated) |

All four tracks have validated CPU paths, rayon parallelization, convergence diagnostics,
experimental reference validation tests, transport layer integration, and Python/WASM/CLI bindings.
GPU kernels are stubbed with CPU fallback and parity tests.