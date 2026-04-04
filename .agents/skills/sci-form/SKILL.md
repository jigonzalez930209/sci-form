---
name: sci-form
description: "Use for implementing, debugging, reviewing, or extending sci-form across Rust, Python, WASM/TypeScript, and CLI. Covers conformer generation, quantum chemistry, spectroscopy, materials, ML, and the test/release workflow."
argument-hint: "implement or fix sci-form feature across Rust/Python/WASM/CLI"
user-invocable: true
disable-model-invocation: false
---

# sci-form

Use this skill for any task in the sci-form repo that touches core chemistry code, bindings, packaging, tests, docs, or performance-sensitive paths. Covers stable production APIs, experimental (alpha/beta) modules with feature gates, GPU acceleration, and parallel execution via rayon.

## Operating rules

- Start from the existing source and tests; do not invent signatures or behavior.
- Keep changes minimal and root-cause focused.
- Preserve the shared data model: `elements` are atomic numbers (`u8`), and `coords` are flat `[x0, y0, z0, ...]` arrays in Å.
- If a public API changes, update every exposed surface that mirrors it: Rust, Python, WASM/JS, CLI, docs, and examples.
- Treat feature flags as part of the implementation contract: `parallel` (rayon), `experimental-gpu` (wgpu), `alpha-*` (proof-of-concept), `beta-*` (validated experimental).
- Alpha modules may change API between minor versions; beta modules are API-stable within minor versions.
- Subpath imports for WASM/TS (`sci-form-wasm/alpha`, `sci-form-wasm/beta`) and Python submodules (`sci_form.alpha`, `sci_form.beta`) must be kept synchronized with Rust feature gates.
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

### Production (always available)

- Conformer generation and topology parsing (ETKDG)
- Charges (Gasteiger), SASA, population analysis, dipole, ESP, and DOS
- UFF/MMFF94 force fields and RMSD/alignment
- PM3 (+ PM3(tm) transition metals), GFN0/GFN1/GFN2-xTB, HF-3c, CISD
- UHF/ROHF open-shell Hartree-Fock (spin contamination, level shift, separate α/β orbitals)
- EHT gradients, orbital grids, band structure, and multi-method DOS
- ANI-2x and ANI-TM atomic environment vectors
- ML descriptors, property prediction, random forest, and gradient boosting
- Stereochemistry (R/S, E/Z, helical, atropisomeric), solvation, rings, fingerprints, and clustering
- IR, UV-Vis, and NMR spectroscopy with peak assignment
- Unit cells, space groups, periodic systems, hapticity, and framework optimization
- CIF import/export with uncertainty notation and space group recognition
- AO→MO 4-index integral transform for post-HF methods
- SMIRKS transforms (single and multi-component) and batch/transport helpers

### Alpha modules (proof-of-concept, may break API)

- **DFT** — Kohn-Sham (SVWN/PBE) with minimal Gaussian basis
- **ReaxFF** — Reactive force field with bond-order-dependent potentials, EEM charges, parallelized central-difference gradients
- **MLFF** — Neural network FF with AEV (Behler-Parrinello) descriptors, parallelized outer atom loop
- **Obara-Saika** — Two-electron repulsion integrals, Boys function, Schwarz screening (ERI engine for DFT/HF)
- **CGA** — Conformal Geometric Algebra motors for exact rigid dihedral rotations
- **GSM** — Growing String Method for reaction path finding and transition state search
- **SDR** — Semidefinite relaxation embedding from distance constraints (convex bounds optimization)
- **Dynamics-Live** — Real-time interactive MD with per-frame callbacks and Berendsen thermostat
- **IMD** — IMD wire-protocol bridge for NAMD/GROMACS steering

### Beta modules (validated, API stable within minor version)

- **KPM** — Kernel Polynomial Method for O(N) DOS via Chebyshev expansion, parallelized random probe vectors
- **MBH** — Mobile Block Hessian for reduced-DOF vibrational analysis on large molecules
- **RandNLA** — Randomized Nyström EHT diagonalization, O(N k²) instead of O(N³)
- **Riemannian** — Geometry on PSD cone: affine-invariant distance, projection, gradient maps
- **CPM** — Constant Potential Method for electrochemical charge equilibration at fixed electrode potential

### GPU acceleration (experimental)

- MMFF94 nonbonded force field via wgpu backend with CPU fallback
- Optional AEV computation offload for ANI descriptors

### Parallelization (via rayon, activated with `features = ["parallel"]`)

- `compute_aevs` — outer atom loop (O(1) per atom with thread count)
- `compute_reaxff_gradient` — central differences (each coordinate displacement independent)
- `ChebyshevExpansion::from_matrix` — random probe vectors (KPM stochastic trace)

## Feature gate patterns

**Rust features** (in `Cargo.toml` or `cargo build --features "..."`):
```
alpha-dft, alpha-reaxff, alpha-mlff, alpha-obara-saika, alpha-cga, alpha-gsm, alpha-sdr, alpha-dynamics-live, alpha-imd
beta-kpm, beta-mbh, beta-randnla, beta-riemannian, beta-cpm
parallel, experimental-gpu
```

**Python invocation** (after `maturin develop`):
```python
from sci_form.alpha import dft_calculate, reaxff_gradient, alpha_compute_aevs
from sci_form.beta import kpm_dos, eht_randnla, cpm_charges
```

**WASM subpath imports** (after `wasm-pack build`):
```typescript
import { alpha_compute_dft, alpha_compute_reaxff_gradient } from 'sci-form-wasm/alpha';
import { beta_compute_kpm_dos, beta_solve_eht_randnla } from 'sci-form-wasm/beta';
```

**CLI**: Alpha and beta commands available when library built with corresponding features.

## Validation checklist

- Rust formatting and linting: `cargo fmt --check`, `cargo clippy --all-targets -- -D warnings`
- Unit tests: `cargo test --lib` (633+ tests)
- Alpha features: `cargo test --lib --features alpha-dft,alpha-reaxff,alpha-mlff`
- Beta features: `cargo test --lib --features beta-kpm,beta-mbh,beta-randnla`
- Parallel paths: `cargo test --lib --features alpha-mlff,beta-kpm,parallel`
- Release smoke suite: `cargo test --release --test ci -- --nocapture`
- Targeted integration tests under `tests/regression/`, `tests/analysis/`, `tests/debug/`, `tests/benchmarks/`, and `tests/integration/`
- Python wheel smoke when bindings change: `maturin develop && python -c "from sci_form.alpha import *; from sci_form.beta import *"`
- WASM/Node tests when JS bindings change: `wasm-pack build && npm test`
- CLI smoke tests: `sci-form embed CCO` and `sci-form compute-* ...`

## Multi-language surface uniformity

When adding or modifying public surfaces, ensure consistency across:
- **Rust**: Direct crate API under feature gates
- **Python**: PyO3 wrappers in `crates/python/src/{alpha,beta}.rs` → `sci_form.{alpha,beta}` submodules
- **WASM/TS**: wasm-bindgen functions in `crates/wasm/src/{alpha,beta}.rs` → JS re-exports in `crates/wasm/js/{alpha,beta}/index.{js,d.ts}` → `sci-form-wasm/{alpha,beta}` subpath imports
- **CLI**: Subcommands in `crates/cli/` under feature gates

JSON serialization boundaries (WASM) and PyO3 class wrapping maintain type safety and zero-copy where possible.

## References

- [Implementation guide](./references/implementation-guide.md)
- [API surface map](./references/api-surface.md)
- [Testing matrix](./references/testing.md)
- [Alpha module docs](../../../docs/experimental/alpha.md)
- [Beta module docs](../../../docs/experimental/beta.md)
- [Experimental overview](../../../docs/experimental/overview.md)
- [Examples](./examples/)
