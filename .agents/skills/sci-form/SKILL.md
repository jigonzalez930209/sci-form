---
name: sci-form
description: "Use for implementing, debugging, reviewing, or extending sci-form across Rust, Python, WASM/TypeScript, and CLI. Covers conformer generation, quantum chemistry, spectroscopy, materials, ML, and the test/release workflow."
argument-hint: "implement or fix sci-form feature across Rust/Python/WASM/CLI"
user-invocable: true
disable-model-invocation: false
---

# sci-form

Use this skill for any task in the sci-form repo that touches core chemistry code, bindings, packaging, tests, docs, or performance-sensitive paths.

## Operating rules

- Start from the existing source and tests; do not invent signatures or behavior.
- Keep changes minimal and root-cause focused.
- Preserve the shared data model: `elements` are atomic numbers (`u8`), and `coords` are flat `[x0, y0, z0, ...]` arrays in Å.
- If a public API changes, update every exposed surface that mirrors it: Rust, Python, WASM/JS, CLI, docs, and examples.
- Treat default `parallel` behavior, `experimental-gpu`, and alpha/beta feature flags as part of the implementation contract.
- Prefer deterministic behavior, low-copy data flow, and existing numerical patterns over new abstractions.

## Workflow

1. Identify the subsystem, public entry point, and existing tests.
2. Read the relevant source, binding layer, and user-facing docs.
3. Make the smallest change that fixes the issue or adds the feature.
4. Propagate any public surface change across language bindings and docs.
5. Validate with targeted tests first, then broader smoke checks if the surface is public.
6. Stop and ask only if the task requires a breaking API decision or the repo has conflicting expectations.

## Surface map

- Rust core: `src/`
- Python bindings: `crates/python/`
- WASM/JS bindings: `crates/wasm/`
- CLI: `crates/cli/`
- Regression, analysis, debug, benchmark, and integration tests: `tests/`
- User docs and examples: `README.md`, `docs/`, `.agents/skills/sci-form/examples/`

## Major domains covered by this repo

- Conformer generation and topology parsing
- Charges, SASA, population analysis, dipole, ESP, and DOS
- UFF/MMFF94 force fields and RMSD/alignment
- PM3, GFN0/GFN1/GFN2-xTB, HF-3c, CISD
- EHT gradients, orbital grids, band structure, and multi-method DOS
- ANI-2x and ANI-TM atomic environment vectors
- ML descriptors, property prediction, random forest, and gradient boosting
- Stereochemistry, solvation, rings, fingerprints, and clustering
- IR, UV-Vis, and NMR spectroscopy
- Unit cells, space groups, periodic systems, hapticity, and framework optimization
- SMIRKS transforms and batch/transport helpers
- Experimental alpha/beta/gpu modules and feature-gated code paths

## Validation checklist

- Rust formatting and linting: `cargo fmt --check`, `cargo clippy --all-targets -- -D warnings`
- Unit tests: `cargo test --lib`
- Release smoke suite: `cargo test --release --test ci -- --nocapture`
- Targeted integration tests under `tests/regression/`, `tests/analysis/`, `tests/debug/`, `tests/benchmarks/`, and `tests/integration/`
- Python wheel smoke and integration tests when Python bindings change
- WASM/Node smoke tests when JS bindings change
- CLI smoke tests when command behavior or serialization changes

## References

- [Implementation guide](./references/implementation-guide.md)
- [API surface map](./references/api-surface.md)
- [Testing matrix](./references/testing.md)
- [Examples](./examples/)
