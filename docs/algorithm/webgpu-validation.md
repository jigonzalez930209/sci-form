# WebGPU Validation

::: warning Status
The **browser harness that earlier versions of this page described does not exist in this repository**. The page has been rewritten to state plainly what is runnable today, what is a documented gap, and where the WebGPU path *is* actually validated.
:::

## What exists today

| Artifact | Path | State |
|----------|------|-------|
| WASM WebGPU runtime | `crates/wasm/src/webgpu.rs` | Present |
| Browser benchmark harness | `crates/wasm/examples/vite-webgpu-hybrid/` | **Missing from the repository** |
| Playwright runner | `tests/integration/test_webgpu_harness.js` | Present, but depends on the missing harness |
| Report comparison tool | `tests/integration/compare_webgpu_reports.js` | Present and runnable |
| Regression policy | `tests/fixtures/webgpu_regression_policy.json` | Present |

## The missing harness is a documented gap

`crates/wasm/examples/vite-webgpu-hybrid/` is referenced by three npm scripts in the root `package.json`:

```json
"webgpu:benchmark:dev":     "npm --prefix crates/wasm/examples/vite-webgpu-hybrid run benchmark:dev",
"webgpu:benchmark:build":   "npm --prefix crates/wasm/examples/vite-webgpu-hybrid run build",
"webgpu:benchmark:preview": "npm --prefix crates/wasm/examples/vite-webgpu-hybrid run benchmark:preview"
```

None of these directory exists in git, so all three commands fail immediately. There is no `benchmark.html`, no `src/benchmark.js` and no `src/benchmarkCases.js` — the files that earlier revisions of this page pointed at.

For the same reason `npm run webgpu:test` cannot work: it runs `tests/integration/test_webgpu_harness.js`, whose first action is to run `npm --prefix <harnessDir> run build` on the missing directory. The runner itself is fine — it launches Chromium with WebGPU flags, waits for `window.__SCI_FORM_WEBGPU_REPORT__` and writes the JSON artifact to disk — but it has nothing to serve.

**No browser benchmark numbers are claimed anywhere in this repository, and none are reproducible from it.** This page deliberately does not publish GPU-vs-CPU timings from a browser run.

## What *is* runnable

### Comparing two reports

`npm run webgpu:compare` works today. It takes two report JSON files and checks the candidate against the baseline using `tests/fixtures/webgpu_regression_policy.json`:

```bash
npm run webgpu:compare -- \
  --baseline build/webgpu/baseline.json \
  --candidate build/webgpu/latest-report.json
```

What it checks:

- a comparison that previously passed must not regress to `fail` or `skip`
- numerical drift (`maxAbsDiff`, `rmsDiff`) must stay within the policy multipliers and floors
- GPU speedup retention against the previous baseline
- optional absolute mode-time regression checks (disabled by default in the policy)

The policy file:

```json
{
  "status": { "requirePassIfBaselinePassed": true, "failOnMissingComparison": true },
  "accuracy": { "maxAbsMultiplier": 2.0, "rmsMultiplier": 2.0, "maxAbsFloor": 1e-8, "rmsFloor": 1e-8 },
  "performance": { "minSpeedupRetention": 0.75, "enforceModeTiming": false, "maxModeSlowdownRatio": 1.35 }
}
```

The two input reports still have to come from somewhere else, because the harness that would produce them is missing.

## Where the WebGPU path is actually covered

The WebGPU compute path is covered on the **Rust side**, by the GPU-vs-CPU gates in `src/gpu/` and by the timing benchmarks in `tests/regression/test_gpu_candidates.rs`. Those run in CI (the `cargo test --lib` unit battery and the `regression` gate) and they are the gates that keep the GPU and CPU implementations consistent.

Concretely, the checks that exist today are:

| Check | Where | What it asserts |
|-------|-------|-----------------|
| `test_fixed_point_accumulation_matches_cpu_energy` | `src/gpu/mmff94_gpu.rs` | The fixed-point GPU accumulation reproduces the CPU MMFF94 energy |
| `test_tiled_cpu_matches_untiled` | `src/gpu/mmff94_gpu.rs` | The tiled CPU path matches the untiled one |
| `test_eri_packing_round_trip_matches_cpu` | `src/gpu/two_electron_gpu.rs` | The GPU ERI packing round-trips to the CPU integrals |
| `test_cpu_fallback_creation`, `test_best_available_never_panics` | `src/gpu/context.rs` | A missing or broken adapter falls back to CPU instead of panicking |
| `test_with_report_falls_back_to_cpu` | `src/gpu/mmff94_gpu.rs` | The activation report records the CPU fallback |
| `bench_*` | `tests/regression/test_gpu_candidates.rs` | Timing of the GPU candidate paths against their CPU baselines (ESP grid, orbital grid, density grid, Fock build, isosurface, sTDA, …) |

Important caveat, stated so it is not over-read: **CI machines have no GPU**. Those gates therefore execute the CPU-fallback branches and check the fallback contract, the kernel arithmetic reproduced on the host, and the shader registry — not a real GPU↔CPU numerical comparison on real hardware.

## The WASM WebGPU layer that the harness used to exercise

The WASM package does expose a WebGPU runtime (`crates/wasm/src/webgpu.rs`):

- `init_webgpu()` — async initialisation, returns a JSON status string
- `webgpu_status()` — current runtime status as a JSON string
- `WasmExecutionMode` — `auto` | `cpu` | `gpu` | `hybrid`, parsed by `parse_execution_mode`

The accelerated surfaces that accept an `execution_mode` argument are:

| Export | Module |
|--------|--------|
| `compute_esp_grid_accelerated`, `compute_esp_grid_accelerated_typed` | `crates/wasm/src/esp.rs` |
| `eht_orbital_grid_accelerated`, `eht_orbital_grid_accelerated_typed` | `crates/wasm/src/eht.rs` |
| `eht_density_grid_accelerated`, `eht_density_grid_accelerated_typed` | `crates/wasm/src/eht.rs` |
| `compute_orbital_mesh_accelerated` | `crates/wasm/src/mesh.rs` |

To build a WASM package with the WebGPU runtime compiled in:

```bash
cd crates/wasm
./build.sh --web-only --web-features "parallel experimental-gpu"
```

The CI `wasm` job builds exactly this configuration (`WASM_WEB_FEATURES: parallel experimental-gpu …`) and verifies the resulting Node.js package, but it does not execute the accelerated paths in a browser.

## If the harness is restored

The runner and the comparison tool are still in place, so restoring the harness is a matter of adding back the directory with:

- `benchmark.html` — the harness page that fills `window.__SCI_FORM_WEBGPU_REPORT__`
- `src/benchmark.js` — the harness logic
- `src/benchmarkCases.js` — the shared molecular cases
- a `package.json` exposing `build`, `preview` and `benchmark:dev` scripts

The report contract the runner expects is: for each case, for each operation, a list of comparisons with `status` (`pass` / `fail` / `skip`), `backend`, `cpuMeanMs`, `modeMeanMs`, `maxAbsDiff` and `rmsDiff`, plus a `summary` block with `passed`, `completed` and `skipped` counters and a `webgpu` block with `initialized` and `error`.

Until then, treat the browser harness as a **documented gap** and the Rust-side GPU-vs-CPU gates as the coverage that actually runs.
