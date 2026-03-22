# sci-form Testing Matrix

## Fast checks

- `cargo fmt --check`
- `cargo clippy --all-targets -- -D warnings`
- `cargo test --lib`

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
- Run `tests/integration/test_python_integration.py` when bindings change

### WASM / Node.js

- Build with `wasm-pack build --target nodejs --release --out-dir pkg-node`
- Run a `require(...)` smoke case against `pkg-node/sci_form_wasm.js`
- Run `tests/integration/test_wasm_integration.js` when JS bindings change

### CLI

- Build the binary with `cargo build --release --package sci-form-cli`
- Run version, `embed`, and `parse` smoke commands

## Validation order for chemistry work

- Start with the smallest reproducible input.
- Validate geometry and topology first.
- Then validate the intermediate numeric outputs that feed the final observable.
- Only after the intermediates match should you check the end-user quantity.

## Release note

- Release flow is tag-driven: pushing `v*` tags triggers release.
- Keep the compact smoke suite in `tests/ci.rs` green before considering a change release-ready.
