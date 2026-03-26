# WebGPU Validation Harness

This harness is the browser-side validation layer for the WebGPU-accelerated WASM paths. It is designed to answer two questions with the same report format:

1. Does `gpu` or `hybrid` stay numerically aligned with the CPU baseline?
2. Where does WebGPU produce real runtime wins over CPU-only execution?

## Scope

The initial harness covers the browser-exposed paths that already implement `cpu`, `gpu`, and `hybrid` execution modes behind the WASM bindings built with `parallel experimental-gpu`:

- ESP grid
- EHT orbital grid
- EHT density grid
- EHT orbital mesh

The page is structured so additional WebGPU surfaces can be added later without changing the JSON report shape.

## Runner layout

- Browser harness page: `crates/wasm/examples/vite-webgpu-hybrid/benchmark.html`
- Harness logic: `crates/wasm/examples/vite-webgpu-hybrid/src/benchmark.js`
- Shared cases: `crates/wasm/examples/vite-webgpu-hybrid/src/benchmarkCases.js`

## How to run

Build the WASM package with WebGPU support and start the Vite example:

```bash
npm install
npm run webgpu:benchmark:dev
```

Equivalent example-local commands remain available inside `crates/wasm/examples/vite-webgpu-hybrid`.

## Automated browser run

The Playwright runner launches Chromium with WebGPU flags, opens the harness page, waits for `window.__SCI_FORM_WEBGPU_REPORT__`, and writes the JSON artifact to disk.

```bash
npm run webgpu:test -- --report-out build/webgpu/baseline.json
```

By default the runner uses the lighter `smoke` profile. Use `--profile full` when you want the full browser workload.

Useful flags:

- `--iterations 5` for more stable timings
- `--profile full` to include the full browser workload instead of the smoke profile
- `--case benzene` to isolate one fixture
- `--allow-skips` to permit CPU-only environments to keep reporting without failing on skipped GPU rows
- `--allow-no-webgpu` to preserve the artifact even when WebGPU initialization fails

Default output path:

- `build/webgpu/latest-report.json`

Then open the benchmark page. Useful query parameters:

- `?iterations=5` to increase timing stability
- `?case=benzene` to isolate one fixture

## Report contract

The page stores the full report in `window.__SCI_FORM_WEBGPU_REPORT__` and can download it as JSON. Each comparison row records:

- molecule case
- operation name
- requested mode (`gpu` or `hybrid`)
- resolved backend string from WASM
- CPU average time
- mode average time
- speedup versus CPU
- maximum absolute deviation
- RMS deviation
- runtime notes and fallback details

## Validation policy

Each operation uses CPU as the baseline and enforces operation-specific tolerances over the returned numeric payload:

- ESP grid: max absolute difference
- EHT orbital grid: max absolute difference
- EHT density grid: max absolute difference
- EHT orbital mesh: max absolute difference over the backing grid plus triangle-count parity

If WebGPU is unavailable, the harness marks GPU comparisons as skipped instead of treating the browser as a failure.

## Why this design

The harness is intentionally browser-first and report-driven:

- it exercises the real WASM + browser WebGPU stack rather than only native `wgpu`
- it produces a stable JSON artifact that a future Playwright or CI runner can scrape without parsing DOM text
- it keeps correctness and timing in the same run so regressions in speed or accuracy show up together

## Next integration step

The intended automation path is to add a browser runner that launches Chromium with WebGPU enabled, opens `benchmark.html`, waits for `window.__SCI_FORM_WEBGPU_REPORT__`, and persists the JSON artifact for comparison across CPU/GPU machines.

That runner now lives in `tests/integration/test_webgpu_harness.js`.

## Historical comparison

The second layer compares a new report against a previously accepted baseline:

```bash
npm run webgpu:compare -- \
	--baseline build/webgpu/baseline.json \
	--candidate build/webgpu/latest-report.json
```

Policy file:

- `tests/fixtures/webgpu_regression_policy.json`

What it checks:

- a comparison that previously passed should not regress to `fail` or `skip`
- numerical drift (`maxAbsDiff`, `rmsDiff`) must stay within policy multipliers and floors
- GPU speedup retention is checked against the previous baseline
- optional absolute mode-time regression checks can be enabled in the policy for same-machine benchmarking

Recommended workflow:

1. Generate a trusted baseline report on the target machine.
2. Store that baseline outside ignored artifact directories if you want repo-tracked history.
3. Compare each new run against that baseline with `npm run webgpu:compare`.