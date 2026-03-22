# sci-form Implementation Guide

## What this skill is for

Use this repo when you need to implement or fix chemistry functionality across the Rust core and the three public surfaces that mirror it: Python, WASM/TypeScript, and the CLI.

## Sources of truth

- Code and tests are the primary source of truth.
- `README.md`, `docs/`, and `TESTING.md` describe intended behavior, supported surfaces, and validation.
- The bundled examples in `.agents/skills/sci-form/examples/` show the intended calling conventions.
- Repo memory notes track release and testing conventions, especially the compact CI suite.

## Data conventions

- `elements` means atomic numbers, stored as `u8` in Rust and integer arrays in other languages.
- `coords` always means a flat coordinate array `[x0, y0, z0, x1, y1, z1, ...]` in Å.
- Any 3D method expects both topology and coordinates; topology-only APIs are the exception, not the rule.
- WASM crossings should use JSON strings for structured data and typed arrays only where the API already exposes them.

## Feature flags and maturity levels

- `parallel` is the default Rust feature and should be preserved unless a change explicitly targets single-threaded behavior.
- `experimental-gpu` implies the GPU rendering and compute stack; keep its code paths isolated from stable CPU logic.
- Alpha and beta features are intentionally non-stable and may move faster than core APIs, but they still need consistent tests and docs.

## How to implement a change

### New public function

1. Add the core Rust implementation first.
2. Expose it through the relevant binding layer(s).
3. Add or update docs and examples.
4. Add a targeted test that exercises the new behavior on the narrowest possible input.

### Bug fix

1. Reproduce the failure with the smallest failing test or smoke case.
2. Fix the root cause in the earliest layer where the bug appears.
3. Verify that the same defect does not exist in the mirrored bindings.

### Serialization or schema change

1. Update Rust return types and any JSON conversion logic.
2. Update Python dataclass-like wrappers.
3. Update WASM JSON schemas and typed-array bridges.
4. Update CLI output and examples if the user-visible shape changes.

### Performance change

1. Keep the behavior identical unless the task explicitly changes semantics.
2. Prefer allocation reductions, parallel iteration, and stable numeric tolerances over new algorithms.
3. Measure on the specific workload that motivated the change before broadening the scope.

## Public surface synchronization rules

- If a Rust function is public and user-facing, assume Python, WASM, CLI, and docs may need matching updates.
- If a function is experimental and only surfaced in one language, keep the scope confined to that language unless the task says otherwise.
- If a change alters naming, units, default parameters, or return ordering, update every dependent example and smoke test.

## What to check before editing

- Which module owns the behavior.
- Whether the code path is core, alpha, beta, or GPU-specific.
- Whether the output is consumed by more than one surface.
- Whether existing regression tests already cover the edge case.
- Whether the change touches release-critical paths like `tests/ci.rs`.

## Common implementation shape

- Core math and chemistry live in Rust modules under `src/`.
- Bindings should stay thin and preserve native semantics.
- Tests should assert scientific invariants, not just snapshot text.
- Docs should explain units, expected inputs, and any feature-gate requirements.

## Practical decision rules

- Prefer the smallest compatible API change.
- Preserve backwards compatibility unless the task explicitly requests a breaking change.
- If the repo already has a clear pattern for a subsystem, follow it instead of introducing a new style.
- Do not guess at chemical or numerical behavior when the source or tests can answer it.
