# sci-form Testing Matrix

## Fast checks

- `cargo fmt --check`
- `cargo clippy --all-targets -- -D warnings`
- `cargo test --lib`

## Experimental feature smoke matrix

- Alpha core smoke: `cargo test --lib --features alpha-dft,alpha-reaxff,alpha-mlff`
- Alpha geometry/path smoke: `cargo test --lib --features alpha-obara-saika,alpha-cga,alpha-gsm,alpha-sdr`
- Alpha dynamics smoke: `cargo test --lib --features alpha-dynamics-live,alpha-imd`
- Beta core smoke: `cargo test --lib --features beta-kpm,beta-mbh,beta-randnla`
- Beta math/electrochem smoke: `cargo test --lib --features beta-riemannian,beta-cpm`
- Parallel smoke: `cargo test --lib --features parallel,alpha-mlff,beta-kpm`

Run narrower combinations first when touching only one subsystem; do not enable every experimental gate unless the task crosses modules.

## Release smoke gate

- `cargo test --release --test ci -- --nocapture`

## Integration areas

- `tests/regression/` for functional chemistry regressions
- `tests/analysis/` for comparison and validation helpers
- `tests/debug/` for investigation harnesses
- `tests/benchmarks/` for large-data timing and scaling checks
- `tests/integration/` for Python and Node.js interop scripts

## Language-specific smoke tests

### Python

- Build a wheel with `maturin build --release` from `crates/python/`
- Install the wheel and run a tiny `embed` smoke case
- When alpha/beta bindings change, also verify submodule import:
	- `python -c "from sci_form.alpha import *; from sci_form.beta import *"`
- Run `tests/integration/test_python_integration.py` when bindings change

### WASM / Node.js

- Build with `wasm-pack build --target nodejs --release --out-dir pkg-node`
- Run a `require(...)` smoke case against `pkg-node/sci_form_wasm.js`
- When alpha/beta exports change, verify subpath imports and generated declarations for `sci-form-wasm/alpha` and `sci-form-wasm/beta`
- Run `tests/integration/test_wasm_integration.js` when JS bindings change

### CLI

- Build the binary with `cargo build --release --package sci-form-cli`
- Run version, `embed`, and `parse` smoke commands
- If CLI exposes experimental flags or subcommands, smoke one command per touched module

## Validation order by change type

### Core chemistry or numerics

- Targeted unit tests in the owning module
- `cargo test --lib`
- Release smoke if the change is user-facing

### Experimental algorithm change

- Narrow feature-gated test run for the touched module
- Cross-surface smoke if bindings are exposed
- Parallel smoke if the touched path has rayon support

### Binding-only change

- Rust compile check for the crate owning the binding
- Python or WASM import smoke
- One end-to-end call through the changed binding

### Documentation-only change for skills

- Confirm relative links resolve from `.agents/skills/sci-form/`
- Keep examples executable and aligned with current public names

## Validation order for chemistry work

- Start with the smallest reproducible input.
- Validate geometry and topology first.
- Then validate the intermediate numeric outputs that feed the final observable.
- Only after the intermediates match should you check the end-user quantity.

## Release note

- Release flow is tag-driven: pushing `v*` tags triggers release.
- Keep the compact smoke suite in `tests/ci.rs` green before considering a change release-ready.
- For AI-agent-facing docs, stale examples are treated as release blockers because they cause bad downstream automation.
