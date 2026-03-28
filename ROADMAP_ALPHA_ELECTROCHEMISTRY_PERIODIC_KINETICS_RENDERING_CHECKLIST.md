# Alpha Roadmap Execution Checklist

> Companion checklist for [ROADMAP_ALPHA_ELECTROCHEMISTRY_PERIODIC_KINETICS_RENDERING.md](/home/lestad/github/sci-form/ROADMAP_ALPHA_ELECTROCHEMISTRY_PERIODIC_KINETICS_RENDERING.md)
>
> Use this file as the execution tracker for the alpha program.

## Sprint Alpha 1: Contracts, adapters, and CPU references

- [ ] Define alpha namespaces and feature flags.
- [ ] Add canonical result schemas for EDL, periodic linear solvers, HTST, and render bridges.
- [ ] Create adapter layers around CPM, ALPB, EEQ, band structure, KPM, RandNLA, GSM, and MBH.
- [ ] Implement CPU reference versions for the first EDL compact-layer solver.
- [ ] Implement CPU reference versions for periodic mesh and operator utilities.
- [ ] Implement CPU reference versions for HTST rate evaluation.
- [ ] Implement CPU reference versions for rendering payload builders.
- [ ] Freeze transport decisions for JSON, typed arrays, and Arrow-style buffers.
- [ ] Add baseline regression fixtures for slabs, periodic cells, GSM paths, and DOS outputs.

## Sprint Alpha 2: Core physics and scaling kernels

- [ ] Complete Helmholtz compact-layer modeling.
- [ ] Complete Gouy-Chapman and Gouy-Chapman-Stern diffuse-layer modeling.
- [ ] Couple CPM output into EDL boundary conditions.
- [ ] Couple EEQ updates into the interface field-response loop.
- [ ] Implement Monkhorst-Pack k-mesh generation.
- [ ] Implement periodic operator assembly for KPM and RandNLA.
- [ ] Extend KPM to periodic DOS evaluation.
- [ ] Extend RandNLA to periodic low-rank refinement.
- [ ] Build the GSM-to-MBH-to-HTST calculation chain.
- [ ] Add single-step HTST rate evaluation across temperature grids.
- [ ] Add microkinetic network assembly and solve routines.
- [ ] Build chart-ready payloads for KPM DOS and HTST observables.
- [ ] Build trajectory-ready payloads for GSM paths.
- [ ] Build profile-ready payloads for EDL outputs.

## Sprint Alpha 3: Parallel CPU scaling and validation hardening

- [ ] Parallelize EDL bias scans with `rayon`.
- [ ] Parallelize k-point workloads with `rayon`.
- [ ] Parallelize HTST temperature sweeps and sensitivity loops.
- [ ] Parallelize rendering chunk generation.
- [ ] Add convergence diagnostics for all alpha solvers.
- [ ] Validate analytic EDL limits.
- [ ] Validate periodic small-cell exact-vs-approximate comparisons.
- [ ] Validate HTST rates against textbook and reference cases.
- [ ] Validate buffer round-tripping for all render bridges.
- [ ] Benchmark CPU scaling and record workload thresholds for GPU promotion.
- [ ] Propagate alpha surfaces to Python, WASM, and CLI.

## Sprint Alpha 4: Selective GPU acceleration

- [ ] GPU-accelerate large EDL profile scans.
- [ ] GPU-batch periodic KPM matrix-vector workloads.
- [ ] GPU-batch periodic RandNLA stages where numerically safe.
- [ ] Keep HTST GPU work limited to large batched sweeps.
- [ ] Add GPU-backed rendering buffer generation where needed.
- [ ] Verify CPU/GPU parity on representative fixtures.
- [ ] Document GPU thresholds and fallback behavior.

## Sprint Alpha 5: Binding rollout and beta-readiness gate

- [ ] Add full Python wrappers for finalized alpha APIs.
- [ ] Add full WASM wrappers for finalized alpha APIs.
- [ ] Add CLI commands for alpha scans, rates, and exports.
- [ ] Provide end-to-end examples for each track.
- [ ] Lock down the benchmark suite.
- [ ] Review which modules are ready to move from alpha to beta staging.
- [ ] Record promotion blockers for any module that is not ready.